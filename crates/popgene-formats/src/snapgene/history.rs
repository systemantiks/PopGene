//! Native history trees and recoverable input snapshots. Decompression and action
//! execution are bounded. Unknown action types remain preserved, never guessed.
use super::{
    Document, Molecule, Node, Result, Source, bad, compression, packet, serialize_node, xml,
};
use crate::snapgene_preserve::packets;
use std::collections::BTreeMap;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NativeHistoryEntry {
    pub id: u32,
    pub parent: Option<u32>,
    pub name: String,
    pub operation: String,
    pub molecule: Molecule,
    /// Native duplex span (may include residues present on only one strand).
    pub sequence_length: usize,
    pub restorable: bool,
    pub attributes: Vec<(String, String)>,
    /// Original non-Node children, flattened into source paths and values.
    pub details: Vec<(String, String)>,
    /// More detail remains in the retained source than the bounded display holds.
    pub details_truncated: bool,
}

struct DetailBudget {
    bytes: usize,
    fields: usize,
}

fn history_details(node: &Node, budget: &mut DetailBudget) -> (Vec<(String, String)>, bool) {
    fn bounded(value: &str, limit: usize, truncated: &mut bool) -> String {
        if value.len() <= limit {
            return value.to_owned();
        }
        *truncated = true;
        let mut end = limit - " …".len();
        while !value.is_char_boundary(end) {
            end -= 1;
        }
        format!("{} …", &value[..end])
    }
    struct Display<'a> {
        rows: Vec<(String, String)>,
        bytes: usize,
        truncated: bool,
        budget: &'a mut DetailBudget,
    }
    impl Display<'_> {
        fn push(&mut self, path: String, value: &str) -> bool {
            let path = bounded(&path, 1024, &mut self.truncated);
            let value = bounded(value, 4096, &mut self.truncated);
            let size = path.len() + value.len();
            if self.rows.len() >= 512
                || self.bytes + size > 64 * 1024
                || self.budget.fields == 0
                || self.budget.bytes < size
            {
                self.truncated = true;
                return false;
            }
            self.bytes += size;
            self.budget.bytes -= size;
            self.budget.fields -= 1;
            self.rows.push((path, value));
            true
        }
        fn walk(&mut self, node: &Node, path: String) -> bool {
            for (key, value) in &node.attrs {
                let key = bounded(key, 512, &mut self.truncated);
                if !self.push(format!("{path}/@{key}"), value) {
                    return false;
                }
            }
            if !node.text.trim().is_empty() && !self.push(format!("{path}/text()"), &node.text) {
                return false;
            }
            if node.attrs.is_empty()
                && node.text.trim().is_empty()
                && node.children.is_empty()
                && !self.push(path.clone(), "(empty)")
            {
                return false;
            }
            for (index, child) in node.children.iter().enumerate() {
                let name = bounded(&child.name, 256, &mut self.truncated);
                let child_path = bounded(
                    &format!("{path}/{name}[{}]", index + 1),
                    1024,
                    &mut self.truncated,
                );
                if !self.walk(child, child_path) {
                    return false;
                }
            }
            true
        }
    }
    let mut display = Display {
        rows: Vec::new(),
        bytes: 0,
        truncated: false,
        budget,
    };
    for (index, child) in node
        .children
        .iter()
        .filter(|child| child.name != "Node")
        .enumerate()
    {
        let name = bounded(&child.name, 256, &mut display.truncated);
        if !display.walk(child, format!("{name}[{}]", index + 1)) {
            break;
        }
    }
    (display.rows, display.truncated)
}
struct History<'a> {
    root: Option<Node>,
    snapshots: BTreeMap<u32, &'a [u8]>,
}
fn parse(source: &Source) -> Result<History<'_>> {
    if source.format != "SnapGene" {
        return Ok(History {
            root: None,
            snapshots: BTreeMap::new(),
        });
    }
    if source.bytes.len() as u64 > crate::MAX_FILE {
        return Err(bad("History source exceeds file limit"));
    }
    let mut root = None;
    let mut snapshots = BTreeMap::new();
    for (kind, body) in packets(&source.bytes)? {
        match kind {
            7 => {
                if root.is_some() {
                    return Err(bad("Duplicate history trees"));
                }
                let decoded = if body.starts_with(b"<") {
                    body.to_vec()
                } else {
                    compression::xdecode(body, 16 * 1024 * 1024)?
                };
                let tree = xml(&decoded)?;
                if tree.name != "HistoryTree"
                    || tree.children.iter().filter(|c| c.name == "Node").count() != 1
                {
                    return Err(bad("Invalid native history tree"));
                }
                root = tree.children.into_iter().find(|c| c.name == "Node");
            }
            11 => {
                let id = u32::from_be_bytes(
                    body.get(..4)
                        .ok_or_else(|| bad("Truncated history snapshot ID"))?
                        .try_into()
                        .unwrap(),
                );
                if snapshots.insert(id, &body[4..]).is_some() {
                    return Err(bad("Duplicate native history snapshot ID"));
                }
            }
            _ => {}
        }
    }
    Ok(History { root, snapshots })
}
fn id(node: &Node) -> Result<u32> {
    node.attr("ID")
        .parse()
        .map_err(|_| bad("Invalid history node ID"))
}
fn molecule(node: &Node) -> Result<Molecule> {
    match node.attr("type") {
        "DNA" => Ok(Molecule::Dna),
        "RNA" => Ok(Molecule::Rna),
        "Protein" | "protein" => Ok(Molecule::Protein),
        _ => Err(bad("Unsupported history molecule type")),
    }
}
struct Snapshot {
    manipulation: Option<Node>,
    sequence: Option<String>,
    content: Vec<u8>,
}
fn snapshot(body: &[u8], budget: &mut usize, include_content: bool) -> Result<Snapshot> {
    let nested = packets(body)?;
    let mut manipulation = None;
    let mut content = None;
    let mut sequence = None;
    for (kind, body) in nested {
        match kind {
            0 | 1 | 21 | 32 => {
                if sequence.is_some() {
                    return Err(bad("Duplicate history snapshot sequence"));
                }
                sequence = Some(if kind == 1 {
                    super::packed::decode(body)?
                } else {
                    std::str::from_utf8(body)
                        .map_err(|_| bad("Invalid history sequence encoding"))?
                        .to_owned()
                });
                *budget = budget
                    .checked_sub(sequence.as_ref().unwrap().len())
                    .ok_or_else(|| bad("Aggregate history expansion exceeds limit"))?;
            }
            29 => {
                if manipulation.is_some() {
                    return Err(bad("Duplicate history manipulation"));
                }
                let decoded = compression::xdecode(body, (*budget).min(16 * 1024 * 1024))?;
                *budget -= decoded.len();
                manipulation = Some(xml(&decoded)?);
            }
            30 => {
                if content.is_some() {
                    return Err(bad("Duplicate history content"));
                }
                content = Some(if include_content {
                    let decoded = compression::xdecode(body, *budget)?;
                    *budget -= decoded.len();
                    decoded
                } else {
                    Vec::new()
                });
            }
            _ => {}
        }
    }
    if manipulation
        .as_ref()
        .is_some_and(|m| m.name != "Manipulation")
        || (manipulation.is_some() == sequence.is_some())
    {
        return Err(bad(
            "History snapshot requires exactly one sequence or manipulation",
        ));
    }
    Ok(Snapshot {
        manipulation,
        sequence,
        content: content.ok_or_else(|| bad("History snapshot lacks content"))?,
    })
}
fn actions_supported(manipulation: &Node) -> bool {
    manipulation.child("Undo").is_some_and(|undo| {
        !undo.children.is_empty()
            && undo.children.iter().all(|action| {
                action.name == "Action"
                    && match action.attr("type") {
                        "REMOVE" => true,
                        "INSERT" => {
                            action.children.len() == 1
                                && action.children[0].name == "Residues"
                                && action.children[0].attr("type") == "GENERIC"
                        }
                        _ => false,
                    }
            })
    })
}
fn apply_undo(sequence: &str, manipulation: &Node, molecule: Molecule) -> Result<String> {
    if !actions_supported(manipulation) {
        return Err(bad(
            "This history step uses an unsupported native undo action; original records are retained",
        ));
    }
    let mut result = sequence.to_owned();
    for action in &manipulation.child("Undo").unwrap().children {
        match action.attr("type") {
            "REMOVE" => {
                let (start, end) = action
                    .attr("range")
                    .split_once('-')
                    .ok_or_else(|| bad("Invalid history removal range"))?;
                let start: usize = start
                    .parse()
                    .map_err(|_| bad("Invalid history coordinate"))?;
                let end: usize = end.parse().map_err(|_| bad("Invalid history coordinate"))?;
                let end = end
                    .checked_add(1)
                    .ok_or_else(|| bad("History coordinate overflow"))?;
                if start >= end || end > result.len() {
                    return Err(bad("History removal is out of bounds"));
                }
                result.replace_range(start..end, "");
            }
            "INSERT" => {
                let position: usize = action
                    .attr("position")
                    .parse()
                    .map_err(|_| bad("Invalid history insertion position"))?;
                if position > result.len() {
                    return Err(bad("History insertion is out of bounds"));
                }
                let residues = molecule.normalize(action.children[0].attr("residues"))?;
                if result
                    .len()
                    .checked_add(residues.len())
                    .is_none_or(|n| n > popgene_core::MAX_BASES)
                {
                    return Err(bad("History sequence exceeds limits"));
                }
                result.insert_str(position, &residues);
            }
            _ => unreachable!(),
        }
    }
    Ok(result)
}
pub fn native_history(source: &Source) -> Result<Vec<NativeHistoryEntry>> {
    let history = parse(source)?;
    let Some(root) = &history.root else {
        return Ok(vec![]);
    };
    fn walk(
        node: &Node,
        parent: Option<u32>,
        path_supported: bool,
        history: &History<'_>,
        result: &mut Vec<NativeHistoryEntry>,
        budget: &mut usize,
        details_budget: &mut DetailBudget,
    ) -> Result<()> {
        let id = id(node)?;
        if result.len() >= 10_000 || result.iter().any(|entry| entry.id == id) {
            return Err(bad("Duplicate or excessive history node IDs"));
        }
        let restored = history
            .snapshots
            .get(&id)
            .and_then(|body| snapshot(body, budget, false).ok());
        let restorable = parent.is_none()
            || restored.as_ref().is_some_and(|snapshot| {
                snapshot.sequence.is_some()
                    || path_supported
                        && snapshot
                            .manipulation
                            .as_ref()
                            .is_some_and(actions_supported)
            });
        let length: usize = node
            .attr("seqLen")
            .parse()
            .map_err(|_| bad("Invalid history sequence length"))?;
        if length > popgene_core::MAX_BASES {
            return Err(bad("History sequence length exceeds limits"));
        }
        let (details, details_truncated) = history_details(node, details_budget);
        result.push(NativeHistoryEntry {
            id,
            parent,
            name: node.attr("name").into(),
            operation: node.attr("operation").into(),
            molecule: molecule(node)?,
            sequence_length: length,
            restorable,
            attributes: node
                .attrs
                .iter()
                .map(|(k, v)| (k.clone(), v.clone()))
                .collect(),
            details,
            details_truncated,
        });
        for child in node.children.iter().filter(|child| child.name == "Node") {
            walk(
                child,
                Some(id),
                restorable,
                history,
                result,
                budget,
                details_budget,
            )?;
        }
        Ok(())
    }
    let mut entries = vec![];
    walk(
        root,
        None,
        true,
        &history,
        &mut entries,
        &mut (crate::MAX_FILE as usize),
        &mut DetailBudget {
            bytes: 1024 * 1024,
            fields: 10_000,
        },
    )?;
    Ok(entries)
}
pub fn native_history_document(source: &Source, target: u32) -> Result<Document> {
    let history = parse(source)?;
    let root = history
        .root
        .as_ref()
        .ok_or_else(|| bad("No native history tree"))?;
    // Validate IDs/lengths once before resolving an unambiguous path.
    native_history(source)?;
    fn path<'a>(node: &'a Node, target: u32, result: &mut Vec<&'a Node>) -> Result<bool> {
        result.push(node);
        if id(node)? == target {
            return Ok(true);
        }
        for child in node.children.iter().filter(|child| child.name == "Node") {
            if path(child, target, result)? {
                return Ok(true);
            }
        }
        result.pop();
        Ok(false)
    }
    let mut nodes = vec![];
    if !path(root, target, &mut nodes)? {
        return Err(bad("Unknown history node"));
    }
    if nodes.len() == 1 {
        return super::read(&source.file_name, &source.bytes);
    }
    let source_packets = packets(&source.bytes)?;
    let sequence_packet = source_packets
        .iter()
        .find(|(kind, _)| matches!(kind, 0 | 21 | 32))
        .ok_or_else(|| bad("Missing source sequence"))?;
    let mut sequence = std::str::from_utf8(&sequence_packet.1[1..])
        .map_err(|_| bad("Invalid history source sequence"))?
        .to_ascii_uppercase();
    // A full snapshot is independently recoverable even when an earlier
    // ancestor uses an unsupported undo operator.
    let start = nodes
        .iter()
        .enumerate()
        .skip(1)
        .filter_map(|(i, node)| {
            history
                .snapshots
                .get(&id(node).ok()?)
                .and_then(|body| packets(body).ok())
                .and_then(|packets| {
                    packets
                        .iter()
                        .any(|(kind, _)| matches!(kind, 0 | 1 | 21 | 32))
                        .then_some(i)
                })
        })
        .next_back()
        .unwrap_or(1);
    let mut budget = crate::MAX_FILE as usize;
    let mut restored = None;
    for node in nodes.into_iter().skip(start) {
        let body = history
            .snapshots
            .get(&id(node)?)
            .ok_or_else(|| bad("History input has no recoverable snapshot"))?;
        let snapshot = snapshot(body, &mut budget, true)?;
        let mol = molecule(node)?;
        sequence = if let Some(sequence) = snapshot.sequence {
            mol.normalize(&sequence)?
        } else {
            apply_undo(&sequence, snapshot.manipulation.as_ref().unwrap(), mol)?
        };
        let content = snapshot.content;
        let mut encoded = vec![];
        let (kind, type_id) = match mol {
            Molecule::Dna => (0, 1_u16),
            Molecule::Rna => (32, 7),
            Molecule::Protein => (21, 2),
        };
        let mut cookie = b"SnapGene".to_vec();
        cookie.extend(type_id.to_be_bytes());
        cookie.extend([0, 15, 0, 20]);
        packet(&mut encoded, 9, &cookie)?;
        let flags = if mol == Molecule::Protein {
            0
        } else {
            u8::from(node.attr("circular") == "1")
                | if node.attr("strandedness") == "single" {
                    0
                } else {
                    2
                }
        };
        let mut seq = vec![flags];
        seq.extend(sequence.as_bytes());
        packet(&mut encoded, kind, &seq)?;
        let content_packets = packets(&content)?;
        let inherited_attachments = !content_packets.iter().any(|(kind, _)| *kind == 23)
            && source_packets.iter().any(|(kind, _)| *kind == 23);
        for (k, b) in content_packets {
            if !matches!(k, 0 | 1 | 7 | 9 | 11 | 21 | 32) {
                packet(&mut encoded, k, b)?;
            }
        }
        // Attachment packets form a file-wide pool rather than sequence undo
        // operations. Preserve access to that pool, without claiming when each
        // attachment was added to a historical sequence.
        if inherited_attachments {
            for (kind, body) in &source_packets {
                if *kind == 23 {
                    packet(&mut encoded, *kind, body)?;
                }
            }
        }
        let name = if node.attr("name").is_empty() {
            "Historical input"
        } else {
            node.attr("name")
        };
        let mut doc = super::read(name, &encoded)?;
        // Snapshot records establish topology/strand but do not document the
        // source sequence's methylation flag byte. Do not turn absent history
        // evidence into a recorded clear state through the reconstructed header.
        doc.methylation = None;
        if sequence.len()
            != node
                .attr("seqLen")
                .parse::<usize>()
                .map_err(|_| bad("Invalid history length"))?
        {
            return Err(bad("Recovered history length disagrees with the tree"));
        }
        doc.warnings.push("Recovered original history input. Open as a separate document and choose a new save path.".into());
        if mol == Molecule::Dna {
            doc.warnings.push("The history snapshot does not establish current Dam/Dcm/EcoKI flags; methylation assumptions are unknown. Any original flag values in Notes are historical provenance.".into());
        }
        if inherited_attachments {
            doc.warnings.push("Embedded files come from the source document's shared attachment pool; the history does not establish which files were present at this historical step.".into());
        }
        restored = Some(doc);
    }
    restored.ok_or_else(|| bad("No recoverable history input"))
}

pub(super) fn append(
    doc: &Document,
    original: &Document,
    source: &Source,
    new_native: &Document,
) -> Result<Vec<(u8, Vec<u8>)>> {
    if doc.molecule != original.molecule {
        return Err(bad("Native history cannot change molecule type"));
    }
    let parsed = parse(source)?;
    let source_packets = packets(&source.bytes)?;
    let (old_sequence_kind, old_sequence) = source_packets
        .iter()
        .find(|(kind, _)| matches!(kind, 0 | 21 | 32))
        .ok_or_else(|| bad("Missing original sequence"))?;
    let flags_note = (*old_sequence_kind == 0 || old_sequence[0] & !3 != 0).then(|| {
        format!(
            "Original imported sequence flags (historical provenance only; not reconstructed as the chemistry of this restored or edited molecule):\n{}",
            super::sequence_flags_text(*old_sequence_kind, old_sequence[0])
        )
    });
    let old_sequence =
        std::str::from_utf8(&old_sequence[1..]).map_err(|_| bad("Invalid original sequence"))?;
    let mut old = parsed
        .root
        .unwrap_or_else(|| node_for(original, 0, old_sequence.len()));
    let entries = native_history(source)?;
    let new_id = entries
        .iter()
        .map(|entry| entry.id)
        .max()
        .unwrap_or(0)
        .checked_add(1)
        .ok_or_else(|| bad("History ID overflow"))?;
    let old_id = id(&old)?;
    old.attrs.insert("resurrectable".into(), "2".into());
    if old.child("Features").is_none()
        && let Some((_, features)) = source_packets.iter().find(|(kind, _)| *kind == 10)
    {
        old.children.push(xml(features)?);
    }
    let mut root = node_for(doc, new_id, new_native.sequence.len());
    root.attrs.insert("operation".into(), "replace".into());
    let mut summary = Node {
        name: "InputSummary".into(),
        ..Default::default()
    };
    summary.attrs.extend([
        ("manipulation".into(), "replace".into()),
        ("val1".into(), "0".into()),
        ("val2".into(), (old_sequence.len() - 1).to_string()),
    ]);
    root.children.push(summary);
    root.children.push(old);
    let tree = Node {
        name: "HistoryTree".into(),
        children: vec![root],
        ..Default::default()
    };
    let mut tree_xml = String::new();
    serialize_node(&tree, &mut tree_xml);
    let mut manipulation = Node {
        name: "Manipulation".into(),
        ..Default::default()
    };
    if let Some((_, props)) = source_packets.iter().find(|(kind, _)| *kind == 8) {
        manipulation.children.push(xml(props)?);
    }
    fn actions(name: &str, remove: usize, insert: &str) -> Node {
        let mut remove_node = Node {
            name: "Action".into(),
            ..Default::default()
        };
        remove_node.attrs.extend([
            ("type".into(), "REMOVE".into()),
            ("range".into(), format!("0-{}", remove - 1)),
        ]);
        let mut insert_node = Node {
            name: "Action".into(),
            ..Default::default()
        };
        insert_node.attrs.extend([
            ("type".into(), "INSERT".into()),
            ("position".into(), "0".into()),
        ]);
        let mut residues = Node {
            name: "Residues".into(),
            ..Default::default()
        };
        residues.attrs.extend([
            ("type".into(), "GENERIC".into()),
            ("residues".into(), insert.into()),
        ]);
        insert_node.children.push(residues);
        Node {
            name: name.into(),
            children: vec![remove_node, insert_node],
            ..Default::default()
        }
    }
    manipulation
        .children
        .push(actions("Redo", old_sequence.len(), &new_native.sequence));
    manipulation
        .children
        .push(actions("Undo", new_native.sequence.len(), old_sequence));
    let mut manipulation_xml = String::new();
    serialize_node(&manipulation, &mut manipulation_xml);
    let mut content = vec![];
    let mut copied_notes = false;
    for (kind, body) in source_packets {
        if kind == 6
            && let Some(note) = &flags_note
        {
            let mut notes = xml(body)?;
            if let Some(comments) = notes.children.iter_mut().find(|n| n.name == "Comments") {
                comments.text.push_str("\n\n");
                comments.text.push_str(note);
            } else {
                notes.children.push(Node {
                    name: "Comments".into(),
                    text: note.clone(),
                    ..Default::default()
                });
            }
            let mut text = String::new();
            serialize_node(&notes, &mut text);
            packet(&mut content, kind, text.as_bytes())?;
            copied_notes = true;
            continue;
        }
        if !matches!(kind, 0 | 1 | 7 | 9 | 11 | 21 | 23 | 32) {
            packet(&mut content, kind, body)?;
        }
    }
    if !copied_notes && let Some(note) = flags_note {
        let notes = Node {
            name: "Notes".into(),
            children: vec![Node {
                name: "Comments".into(),
                text: note,
                ..Default::default()
            }],
            ..Default::default()
        };
        let mut text = String::new();
        serialize_node(&notes, &mut text);
        packet(&mut content, 6, text.as_bytes())?;
    }
    let mut snapshot = old_id.to_be_bytes().to_vec();
    packet(
        &mut snapshot,
        29,
        &compression::xencode(manipulation_xml.as_bytes())?,
    )?;
    packet(&mut snapshot, 30, &compression::xencode(&content)?)?;
    let mut result = vec![
        (7, compression::xencode(tree_xml.as_bytes())?),
        (11, snapshot),
    ];
    for (id, body) in parsed.snapshots {
        if id == old_id {
            return Err(bad(
                "Current native history root already has an input snapshot",
            ));
        }
        let mut value = id.to_be_bytes().to_vec();
        value.extend(body);
        result.push((11, value));
    }
    Ok(result)
}
fn node_for(doc: &Document, id: u32, native_length: usize) -> Node {
    let mut node = Node {
        name: "Node".into(),
        ..Default::default()
    };
    node.attrs.extend([
        ("ID".into(), id.to_string()),
        ("name".into(), doc.name.clone()),
        ("seqLen".into(), native_length.to_string()),
        (
            "type".into(),
            match doc.molecule {
                Molecule::Dna => "DNA",
                Molecule::Rna => "RNA",
                Molecule::Protein => "Protein",
            }
            .into(),
        ),
        (
            "strandedness".into(),
            if doc.single_stranded {
                "single"
            } else {
                "double"
            }
            .into(),
        ),
        ("circular".into(), u8::from(doc.circular).to_string()),
        ("operation".into(), "invalid".into()),
    ]);
    if let Some(ends) = &doc.duplex_ends {
        use popgene_core::restriction::EndKind;
        for (key, end) in [
            ("upstreamStickiness", &ends.left),
            ("downstreamStickiness", &ends.right),
        ] {
            let value = match end.kind {
                EndKind::FivePrime => end.overhang.len() as isize,
                EndKind::ThreePrime => -(end.overhang.len() as isize),
                _ => 0,
            };
            if value != 0 {
                node.attrs.insert(key.into(), value.to_string());
            }
        }
        for (key, value) in [
            ("upstreamModification", &ends.left_modification),
            ("downstreamModification", &ends.right_modification),
        ] {
            node.attrs.insert(
                key.into(),
                value.clone().unwrap_or_else(|| "Unmodified".into()),
            );
        }
    }
    node
}

#[cfg(test)]
mod tests {
    use super::*;
    use popgene_core::{Feature, Segment};
    #[test]
    fn history_operation_inputs_are_visible_without_reinterpreting_them() {
        let doc = Document::new("history source", "ACGT").unwrap();
        let mut bytes = super::super::write(&doc).unwrap();
        packet(&mut bytes, 7, br#"<HistoryTree><Node ID="0" name="PCR product" type="DNA" seqLen="4" operation="amplifyFragment"><InputSummary manipulation="amplify" val1="23"/><Oligo name="forward" sequence="gcatACGT" phosphorylated="1"/><Parameter name="annealingTemperature" val="58"/><Features><Feature name="old annotation"><Segment range="1-4"/></Feature></Features></Node></HistoryTree>"#).unwrap();
        let source = super::super::read("history.dna", &bytes)
            .unwrap()
            .source
            .unwrap();
        let entries = native_history(&source).unwrap();
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].operation, "amplifyFragment");
        assert!(!entries[0].details_truncated);
        for (path, value) in [
            ("InputSummary[1]/@manipulation", "amplify"),
            ("Oligo[2]/@sequence", "gcatACGT"),
            ("Parameter[3]/@val", "58"),
            ("Features[4]/Feature[1]/Segment[1]/@range", "1-4"),
        ] {
            assert!(
                entries[0]
                    .details
                    .iter()
                    .any(|(k, v)| k == path && v == value)
            );
        }
        let restored = native_history_document(&source, 0).unwrap();
        assert_eq!(restored.sequence, doc.sequence);
        assert!(restored.features.is_empty());
    }
    #[test]
    fn history_detail_limits_are_aggregate_and_utf8_safe() {
        let mut root = Node::default();
        root.children.push(Node {
            name: "Parameter".into(),
            text: "α".repeat(5000),
            ..Default::default()
        });
        root.children.extend((0..1000).map(|i| Node {
            name: "Parameter".into(),
            text: i.to_string(),
            ..Default::default()
        }));
        let mut budget = DetailBudget {
            bytes: 1024 * 1024,
            fields: 10_000,
        };
        let (details, truncated) = history_details(&root, &mut budget);
        assert!(truncated);
        assert_eq!(details.len(), 512);
        assert!(details[0].1.ends_with(" …"));
        assert!(details[0].1.len() <= 4096);
        let used: usize = details.iter().map(|(k, v)| k.len() + v.len()).sum();
        assert!(used <= 64 * 1024);
        assert_eq!(budget.bytes, 1024 * 1024 - used);
        assert_eq!(budget.fields, 10_000 - details.len());
        let mut budget = DetailBudget {
            bytes: 100,
            fields: 1,
        };
        let small = xml(b"<Node><Parameter name=\"temperature\" val=\"58\"/></Node>").unwrap();
        let (first, truncated) = history_details(&small, &mut budget);
        assert!(truncated);
        assert_eq!(first.len(), 1);
        let (second, truncated) = history_details(&small, &mut budget);
        assert!(truncated);
        assert!(second.is_empty());
    }
    #[test]
    fn recovered_input_keeps_file_wide_attachments_with_explicit_provenance() {
        let doc = Document::new("embedded history", "AAACCCGGGTTT").unwrap();
        let mut bytes = super::super::write(&doc).unwrap();
        packet(
            &mut bytes,
            23,
            b"<Files><File id=\"7\" name=\"evidence.bin\" size=\"3\" compressible=\"0\"/></Files>",
        )
        .unwrap();
        let mut payload = 7_u32.to_be_bytes().to_vec();
        payload.extend([0, 255, 1]);
        packet(&mut bytes, 23, &payload).unwrap();
        let mut edited = super::super::read("embedded.dna", &bytes).unwrap();
        edited.replace(0, 1, "T").unwrap();
        let bytes = super::super::write(&edited).unwrap();
        let edited = super::super::read("edited.dna", &bytes).unwrap();
        let recovered = native_history_document(edited.source.as_ref().unwrap(), 0).unwrap();
        assert_eq!(recovered.sequence, doc.sequence);
        assert_eq!(
            super::super::attachments::native_attachment_bytes(
                recovered.source.as_ref().unwrap(),
                7
            )
            .unwrap(),
            vec![0, 255, 1]
        );
        assert!(
            recovered
                .warnings
                .iter()
                .any(|warning| warning.contains("shared attachment pool"))
        );
    }
    #[test]
    fn successive_edits_restore_each_input_with_original_annotations_and_citations() {
        let mut original = Document::new("history test", "AAACCCGGGTTT").unwrap();
        original.references.push(popgene_core::metadata::Reference {
            title: "Synthetic citation".into(),
            ..Default::default()
        });
        original.features.push(Feature {
            id: 1,
            label: "original domain".into(),
            kind: "misc_feature".into(),
            parts: vec![Segment { start: 3, end: 9 }],
            location: None,
            reverse: false,
            qualifiers: vec![],
            edited: false,
        });
        let bytes = super::super::write(&original).unwrap();
        let mut first = super::super::read("first.dna", &bytes).unwrap();
        first.replace(0, 1, "T").unwrap();
        let bytes = super::super::write(&first).unwrap();
        let mut second = super::super::read("second.dna", &bytes).unwrap();
        second.replace(10, 12, "ACGT").unwrap();
        let bytes = super::super::write(&second).unwrap();
        let latest = super::super::read("latest.dna", &bytes).unwrap();
        let source = latest.source.as_ref().unwrap();
        let entries = native_history(source).unwrap();
        assert_eq!(
            entries.iter().map(|e| e.id).collect::<Vec<_>>(),
            vec![2, 1, 0]
        );
        assert!(entries.iter().all(|e| e.restorable));
        let restored = native_history_document(source, 0).unwrap();
        assert_eq!(restored.sequence, original.sequence);
        assert_eq!(restored.features[0].parts, original.features[0].parts);
        assert_eq!(restored.references[0].title, original.references[0].title);
        assert_eq!(
            native_history_document(source, 1).unwrap().sequence,
            first.sequence
        );
        assert_eq!(
            native_history_document(source, 2).unwrap().sequence,
            latest.sequence
        );
    }
    #[test]
    fn undo_bounds_and_unknown_operations_never_rewrite_sequence() {
        let m = xml(
            b"<Manipulation><Undo><Action type=\"REMOVE\" range=\"0-99\"/></Undo></Manipulation>",
        )
        .unwrap();
        assert!(apply_undo("ACGT", &m, Molecule::Dna).is_err());
        let m =
            xml(b"<Manipulation><Undo><Action type=\"MYSTERY\"/></Undo></Manipulation>").unwrap();
        assert!(apply_undo("ACGT", &m, Molecule::Dna).is_err());
    }
}
