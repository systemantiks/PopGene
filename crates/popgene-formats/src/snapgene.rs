//! Native .dna/.rna/.prot interchange using observed TLV/XML records.
//! RNA/protein packet identifiers and flags were checked with SnapGene 8.2.3.
//! Unchanged files retain their bytes. Edited files rebuild supported records; opaque
//! history, alignments and other source blocks must never describe the edited molecule.
use super::{Error, Result};
use popgene_core::{
    Document, Feature, FeatureDirection, FeatureLocation, Molecule, Segment, Source,
};
use quick_xml::{Reader, events::Event};
use std::{collections::BTreeMap, path::Path};

mod annotation_order;
mod attachments;
mod cleavage;
mod colors;
mod compression;
mod ends;
mod enzyme_panels;
#[cfg(test)]
mod exact_locations;
#[cfg(test)]
mod feature_colors;
mod history;
mod packed;
mod primer_envelope;
mod record_preservation;
mod references;
mod segment_types;
mod stored_records;
pub use attachments::{NativeAttachment, native_attachment_bytes, native_attachments};
pub use enzyme_panels::{NativeEnzymePanel, NativeEnzymePanels, native_enzyme_panels};
pub use history::{NativeHistoryEntry, native_history, native_history_document};
pub use stored_records::{NativeEmbeddedTrace, native_embedded_trace, native_embedded_traces};

const REVIEW_NOTE: &str = "Sequence changed in PopGene; annotation needs review.";

fn bad(message: impl Into<String>) -> Error {
    Error::Message(format!("SnapGene import: {}", message.into()))
}
#[derive(Clone, Default)]
struct Node {
    name: String,
    attrs: BTreeMap<String, String>,
    text: String,
    children: Vec<Node>,
}
impl Node {
    fn attr(&self, k: &str) -> &str {
        self.attrs.get(k).map(String::as_str).unwrap_or("")
    }
    fn child(&self, k: &str) -> Option<&Node> {
        self.children.iter().find(|n| n.name == k)
    }
}
fn xml(bytes: &[u8]) -> Result<Node> {
    if bytes.len() > 16 * 1024 * 1024 {
        return Err(bad("XML block exceeds 16 MiB"));
    }
    let mut r = Reader::from_reader(bytes);
    r.config_mut().expand_empty_elements = true;
    let mut stack = vec![Node::default()];
    let mut count = 0;
    loop {
        match r.read_event().map_err(|e| bad(e.to_string()))? {
            Event::Start(e) => {
                count += 1;
                if count > 100_000 || stack.len() > 64 {
                    return Err(bad("XML structure exceeds limits"));
                }
                let mut n = Node {
                    name: String::from_utf8_lossy(e.name().as_ref()).into_owned(),
                    ..Default::default()
                };
                for a in e.attributes() {
                    let a = a.map_err(|e| bad(e.to_string()))?;
                    n.attrs.insert(
                        String::from_utf8_lossy(a.key.as_ref()).into_owned(),
                        a.decoded_and_normalized_value(
                            quick_xml::XmlVersion::Implicit1_0,
                            r.decoder(),
                        )
                        .map_err(|e| bad(e.to_string()))?
                        .into_owned(),
                    );
                }
                stack.push(n);
            }
            Event::End(_) => {
                if stack.len() < 2 {
                    return Err(bad("Unbalanced XML"));
                }
                let n = stack.pop().unwrap();
                stack.last_mut().unwrap().children.push(n);
            }
            Event::Text(t) => stack
                .last_mut()
                .unwrap()
                .text
                .push_str(&t.decode().map_err(|e| bad(e.to_string()))?),
            Event::CData(t) => stack
                .last_mut()
                .unwrap()
                .text
                .push_str(&t.decode().map_err(|e| bad(e.to_string()))?),
            Event::GeneralRef(t) => {
                if let Some(c) = t.resolve_char_ref().map_err(|e| bad(e.to_string()))? {
                    stack.last_mut().unwrap().text.push(c);
                } else {
                    let name = t.decode().map_err(|e| bad(e.to_string()))?;
                    let c = match name.as_ref() {
                        "amp" => '&',
                        "lt" => '<',
                        "gt" => '>',
                        "quot" => '"',
                        "apos" => '\'',
                        _ => return Err(bad("Unknown XML entity")),
                    };
                    stack.last_mut().unwrap().text.push(c);
                }
            }
            Event::DocType(_) => return Err(bad("DTD declarations are not supported")),
            Event::Eof => break,
            _ => {}
        }
    }
    if stack.len() != 1 || stack[0].children.len() != 1 {
        return Err(bad("Invalid XML root"));
    }
    Ok(stack.pop().unwrap().children.remove(0))
}
fn plain(s: &str) -> String {
    // SnapGene rich text is an HTML document inside an XML string. Plain names and
    // qualifiers can legitimately contain '<', newlines or repeated spaces.
    if !s.trim_start().starts_with("<html") && !s.trim_start().starts_with("<!DOCTYPE HTML") {
        return s.to_string();
    }
    let mut in_tag = false;
    let mut out = String::new();
    for c in s.chars() {
        match c {
            '<' => in_tag = true,
            '>' => {
                in_tag = false;
                out.push(' ');
            }
            _ if !in_tag => out.push(c),
            _ => {}
        }
    }
    let out = out.trim().to_string();
    quick_xml::escape::unescape(&out)
        .map(|s| s.into_owned())
        .unwrap_or(out)
}
fn ranges(s: &str, n: usize, primer: bool, circular: bool) -> Result<Vec<Segment>> {
    let (a, b) = s
        .split_once('-')
        .ok_or_else(|| bad("Missing feature range"))?;
    let mut a = a.parse::<usize>().map_err(|_| bad("Invalid range start"))?;
    let mut b = b.parse::<usize>().map_err(|_| bad("Invalid range end"))?;
    if primer {
        b = b.checked_add(1).ok_or_else(|| bad("Range overflow"))?;
    } else {
        a = a
            .checked_sub(1)
            .ok_or_else(|| bad("Feature range starts at zero"))?;
    }
    if a >= n || b == 0 || b > n {
        return Err(bad("Feature range is outside the sequence"));
    }
    if a < b {
        Ok(vec![Segment { start: a, end: b }])
    } else if circular {
        Ok(vec![
            Segment { start: a, end: n },
            Segment { start: 0, end: b },
        ])
    } else {
        Err(bad("Wrapping range on a linear sequence"))
    }
}

/// Original Notes plus annotation/primer/end/color settings. These are provenance,
/// not editable current state or a substitute for the dedicated history parser.
fn sequence_flags_text(kind: u8, flags: u8) -> String {
    let mut value = format!("Raw sequence flags: 0x{flags:02x}");
    if kind == 0 {
        for (label, mask) in [
            ("Dam whole-template flag", 4),
            ("Dcm whole-template flag", 8),
            ("EcoKI whole-template flag", 16),
        ] {
            value.push_str(&format!(
                "\n{label}: {}",
                if flags & mask != 0 { "set" } else { "clear" }
            ));
        }
    }
    value
}

pub(super) fn source_metadata(bytes: &[u8]) -> Result<Vec<(String, String)>> {
    fn describe(node: &Node, output: &mut String, depth: usize) {
        if output.len() > 1024 * 1024 {
            return;
        }
        for (key, value) in &node.attrs {
            output.push_str(&format!("{key}: {}\n", plain(value)));
        }
        if !node.text.trim().is_empty() {
            output.push_str(&plain(&node.text));
            output.push('\n');
        }
        for child in &node.children {
            if output.len() > 1024 * 1024 {
                break;
            }
            output.push_str(&format!("{}{}:\n", "  ".repeat(depth), child.name));
            describe(child, output, depth + 1);
        }
    }
    let mut offset = 0usize;
    let mut count = 0;
    let mut fields = vec![];
    let mut total = 0;
    while offset < bytes.len() {
        if bytes.len() - offset < 5 {
            return Err(bad("Truncated block header"));
        }
        count += 1;
        if count > 10_000 {
            return Err(bad("Too many blocks"));
        }
        let kind = bytes[offset];
        let size = u32::from_be_bytes(bytes[offset + 1..offset + 5].try_into().unwrap()) as usize;
        offset += 5;
        let end = offset
            .checked_add(size)
            .filter(|end| *end <= bytes.len())
            .ok_or_else(|| bad("Truncated block payload"))?;
        if matches!(kind, 0 | 21 | 32)
            && let Some(flags) = bytes[offset..end].first()
        {
            let value = sequence_flags_text(kind, *flags);
            total += value.len() + "Original sequence flags".len();
            fields.push(("Original sequence flags".into(), value));
        }
        if kind == 6 {
            let root = xml(&bytes[offset..end])?;
            if root.name != "Notes" {
                return Err(bad("Invalid notes root"));
            }
            for node in &root.children {
                let mut value = String::new();
                describe(node, &mut value, 0);
                total += value.len() + node.name.len();
                if total > 1024 * 1024 || fields.len() >= 10_000 {
                    fields.push(("Display limit".into(), "Additional notes remain in the retained original file (metadata display limited to 1 MiB / 10,000 fields).".into()));
                    return Ok(fields);
                }
                fields.push((node.name.clone(), value.trim_end().to_owned()));
            }
        }
        offset = end;
    }
    // Notes/references take priority over potentially much larger annotation
    // settings. Surface unmodeled source attributes instead of hiding them.
    for (kind, body) in super::snapgene_preserve::packets(bytes)? {
        if !matches!(kind, 5 | 8 | 10 | 20) {
            continue;
        }
        let root = xml(body)?;
        let mut nodes: Vec<(&Node, String)> = if matches!(kind, 5 | 10) {
            root.children
                .iter()
                .enumerate()
                .map(|(index, node)| {
                    (
                        node,
                        format!(
                            "Original {} {} ({})",
                            node.name,
                            index + 1,
                            plain(node.attr("name"))
                        ),
                    )
                })
                .collect()
        } else {
            vec![(&root, format!("Original {}", root.name))]
        };
        let header = Node {
            attrs: root.attrs.clone(),
            ..Default::default()
        };
        if matches!(kind, 5 | 10) && !root.attrs.is_empty() {
            nodes.insert(0, (&header, format!("Original {} settings", root.name)));
        }
        for (node, label) in nodes {
            let mut value = String::new();
            describe(node, &mut value, 0);
            total += value.len() + label.len();
            if total > 1024 * 1024 || fields.len() >= 10_000 {
                fields.push(("Display limit".into(),"Additional native settings remain in the retained source file (metadata display limited to 1 MiB / 10,000 fields).".into()));
                return Ok(fields);
            }
            fields.push((label, value.trim_end().to_owned()));
        }
    }
    stored_records::append_metadata(bytes, &mut fields)?;
    Ok(fields)
}

pub fn read(name: &str, bytes: &[u8]) -> Result<Document> {
    let mut offset = 0usize;
    let mut blocks = vec![];
    while offset < bytes.len() {
        if bytes.len() - offset < 5 {
            return Err(bad("Truncated block header"));
        }
        let kind = bytes[offset];
        let size = u32::from_be_bytes(bytes[offset + 1..offset + 5].try_into().unwrap()) as usize;
        offset += 5;
        let end = offset
            .checked_add(size)
            .filter(|&e| e <= bytes.len())
            .ok_or_else(|| bad("Truncated block payload"))?;
        blocks.push((kind, &bytes[offset..end]));
        offset = end;
        if blocks.len() > 10_000 {
            return Err(bad("Too many blocks"));
        }
    }
    let cookie = blocks
        .first()
        .filter(|(k, b)| *k == 9 && b.len() == 14 && &b[..8] == b"SnapGene")
        .ok_or_else(|| bad("Invalid SnapGene header"))?
        .1;
    let (molecule, sequence_kind) = match u16::from_be_bytes(cookie[8..10].try_into().unwrap()) {
        1 => (Molecule::Dna, 0),
        2 => (Molecule::Protein, 21),
        7 => (Molecule::Rna, 32),
        value => return Err(bad(format!("Unsupported native document type {value}"))),
    };
    let sequence_blocks: Vec<_> = blocks
        .iter()
        .filter(|(kind, _)| matches!(kind, 0 | 1 | 21 | 32))
        .collect();
    if sequence_blocks.len() != 1
        || sequence_blocks[0].0 != sequence_kind
        || sequence_blocks[0].1.len() < 2
    {
        return Err(bad(
            "Expected exactly one sequence block matching the document type",
        ));
    }
    let data = sequence_blocks[0].1;
    let seq = std::str::from_utf8(&data[1..]).map_err(|_| bad("Invalid sequence encoding"))?;
    let title = Path::new(name)
        .file_stem()
        .unwrap_or_default()
        .to_string_lossy();
    // The observed RNA packet contains U; accepting T silently would change
    // the native sequence. Inconsistent alphabets are rejected by the model.
    let mut doc = Document::new_typed(title, seq, molecule)?;
    if molecule == Molecule::Protein && data[0] & 3 != 0 {
        return Err(bad(
            "Protein sequence has nucleotide topology/strandedness flags",
        ));
    }
    doc.circular = data[0] & 1 != 0;
    doc.single_stranded = molecule != Molecule::Protein && data[0] & 2 == 0;
    doc.methylation =
        (molecule == Molecule::Dna).then_some(popgene_core::metadata::DnaMethylation {
            dam: data[0] & 4 != 0,
            dcm: data[0] & 8 != 0,
            eco_ki: data[0] & 16 != 0,
        });
    let n = doc.sequence.len();
    let mut unsupported_segments = 0usize;
    for (kind, body) in &blocks {
        match kind {
            10 => {
                let root = xml(body)?;
                if root.name != "Features" {
                    return Err(bad("Invalid features root"));
                }
                for f in root.children.iter().filter(|f| f.name == "Feature") {
                    let cleavage_sites = cleavage::read_attribute(f, n)?;
                    let mut parts = vec![];
                    let mut qualifiers = vec![];
                    let mut unsupported_types = segment_types::SourceTypes::default();
                    let mut display_segments = vec![];
                    let native_segments: Vec<_> =
                        f.children.iter().filter(|s| s.name == "Segment").collect();
                    for seg in f.children.iter().filter(|s| s.name == "Segment") {
                        let ranges = ranges(seg.attr("range"), n, false, doc.circular)?;
                        if unsupported_types.retain(seg.attr("type"), &mut qualifiers) {
                            unsupported_segments += 1;
                        }
                        if seg.attr("type") == "gap" {
                            for range in ranges {
                                let color = if seg.attr("color").is_empty() {
                                    "#------"
                                } else {
                                    seg.attr("color")
                                };
                                let mut value =
                                    format!("{} .. {} / {color}", range.start + 1, range.end);
                                if !seg.attr("name").is_empty() {
                                    value.push_str(&format!(" / {}", plain(seg.attr("name"))));
                                }
                                qualifiers.push(("snapgene_gap".into(), Some(value)));
                            }
                            continue;
                        }
                        for range in &ranges {
                            display_segments.push((
                                range.clone(),
                                seg.attr("color"),
                                plain(seg.attr("name")),
                            ));
                        }
                        parts.extend(ranges);
                        if !seg.attr("color").is_empty() {
                            qualifiers
                                .push(("snapgene_color".into(), Some(seg.attr("color").into())));
                        }
                        if !seg.attr("name").is_empty() {
                            qualifiers.push(("segment_name".into(), Some(seg.attr("name").into())));
                        }
                    }
                    let location = if f.attr("classification") == "1" {
                        if native_segments.len() != 1 {
                            return Err(bad("Between-base feature needs one segment"));
                        }
                        let (left, right) = native_segments[0]
                            .attr("range")
                            .split_once('-')
                            .ok_or_else(|| bad("Invalid between-base range"))?;
                        let coordinate = |value: &str| {
                            value
                                .parse::<usize>()
                                .ok()
                                .and_then(|n| n.checked_sub(1))
                                .ok_or_else(|| bad("Invalid between-base coordinate"))
                        };
                        let mut location = FeatureLocation::Between {
                            left: coordinate(left)?,
                            right: coordinate(right)?,
                        };
                        if f.attr("directionality") == "2" {
                            location = FeatureLocation::Complement(Box::new(location));
                        }
                        location
                            .validate(n, doc.circular)
                            .map_err(|e| bad(e.to_string()))?;
                        parts.clear();
                        Some(location)
                    } else {
                        if parts.is_empty() {
                            return Err(bad("Feature has no location"));
                        }
                        None
                    };
                    for q in f.children.iter().filter(|q| q.name == "Q") {
                        if q.attr("name").is_empty() {
                            return Err(bad("Qualifier has no name"));
                        }
                        if q.children.is_empty() {
                            qualifiers.push((q.attr("name").into(), None));
                        }
                        for v in &q.children {
                            let value = v
                                .attrs
                                .get("text")
                                .or_else(|| v.attrs.get("predef"))
                                .or_else(|| v.attrs.get("int"));
                            qualifiers.push((q.attr("name").into(), value.map(|s| plain(s))));
                        }
                    }
                    let label = plain(f.attr("name"));
                    if !qualifiers.iter().any(|(k, _)| k == "label") {
                        qualifiers.push(("label".into(), Some(label.clone())));
                    }
                    qualifiers.push((
                        "snapgene_directionality".into(),
                        Some(f.attr("directionality").into()),
                    ));
                    // Signed feature-relative frame semantics were checked against
                    // whole-sequence and shifted SnapGene fixtures. Preserve raw values and
                    // conflicts, deriving an INSDC offset only when unambiguous.
                    if let Some(frame) = f.attrs.get("readingFrame") {
                        qualifiers.push(("snapgene_readingFrame".into(), Some(frame.clone())));
                        if f.attr("type") == "CDS"
                            && !qualifiers.iter().any(|(key, _)| key == "codon_start")
                            && qualifiers
                                .iter()
                                .filter(|(key, _)| key == "snapgene_readingFrame")
                                .all(|(_, value)| value.as_deref() == Some(frame.as_str()))
                            && let Ok(frame) = frame.parse::<i8>()
                            && matches!(frame, -3..=-1 | 1..=3)
                            && (frame < 0) == (f.attr("directionality") == "2")
                        {
                            qualifiers.push((
                                "codon_start".into(),
                                Some(frame.unsigned_abs().to_string()),
                            ));
                        }
                    }
                    // Encode explicit native segment styling in the published
                    // conventional note grammar. Coordinate associations survive
                    // missing names and wrapped ranges; no private blob is added.
                    if (native_segments.len() > 1
                        || display_segments.iter().any(|(_, _, name)| !name.is_empty()))
                        && display_segments.iter().all(|(_, color, _)| {
                            color.len() == 7
                                && color.starts_with('#')
                                && u32::from_str_radix(&color[1..], 16).is_ok()
                        })
                    {
                        let direction = match f.attr("directionality") {
                            "1" => "forward directional ",
                            "2" => "reverse directional ",
                            "3" => "bidirectional ",
                            _ => "",
                        };
                        let mut note = format!(
                            "This {direction}feature has {} segments:",
                            display_segments.len()
                        );
                        for (index, (range, color, name)) in display_segments.iter().enumerate() {
                            note.push_str(&format!(
                                " {}: {} .. {} / {}",
                                index + 1,
                                range.start + 1,
                                range.end,
                                color.to_ascii_lowercase()
                            ));
                            if !name.is_empty() {
                                note.push_str(&format!(" / {name}"));
                            }
                        }
                        cleavage::append_note(&mut note, &cleavage_sites);
                        if qualifiers
                            .iter()
                            .rev()
                            .find(|(key, _)| key == "note")
                            .and_then(|(_, value)| value.as_deref())
                            != Some(&note)
                        {
                            qualifiers.retain(|(key, value)| {
                                key != "note" || value.as_deref() != Some(&note)
                            });
                            qualifiers.push(("note".into(), Some(note)));
                        }
                    }
                    let reverse = location.is_none() && f.attr("directionality") == "2";
                    let mut feature = Feature {
                        location,
                        id: doc.features.len() as u64 + 1,
                        label,
                        kind: if f.attr("type").is_empty() {
                            "misc_feature".into()
                        } else {
                            f.attr("type").into()
                        },
                        parts,
                        reverse,
                        qualifiers,
                        edited: false,
                    };
                    cleavage::apply(&mut feature, &cleavage_sites)?;
                    doc.features.push(feature);
                }
            }
            5 => {
                let root = xml(body)?;
                if root.name != "Primers" {
                    return Err(bad("Invalid primers root"));
                }
                for primer in root.children.iter().filter(|p| p.name == "Primer") {
                    let mut seen = vec![];
                    let mut count = 0;
                    for site in primer.children.iter().filter(|s| s.name == "BindingSite") {
                        let parts = ranges(site.attr("location"), n, true, doc.circular)?;
                        let reverse = site.attr("boundStrand") == "1";
                        if seen.contains(&(parts.clone(), reverse)) {
                            continue;
                        }
                        seen.push((parts.clone(), reverse));
                        count += 1;
                        let mut qualifiers = vec![
                            ("label".into(), Some(plain(primer.attr("name")))),
                            (
                                "primer_sequence".into(),
                                Some(primer.attr("sequence").into()),
                            ),
                        ];
                        for (key, value) in &primer.attrs {
                            if key != "name" && key != "sequence" {
                                qualifiers.push((format!("snapgene_{key}"), Some(value.clone())));
                            }
                        }
                        for (key, value) in &site.attrs {
                            qualifiers.push((format!("binding_{key}"), Some(value.clone())));
                        }
                        doc.features.push(Feature {
                            location: None,
                            id: doc.features.len() as u64 + 1,
                            label: plain(primer.attr("name")),
                            kind: "primer_bind".into(),
                            parts,
                            reverse,
                            qualifiers,
                            edited: false,
                        });
                    }
                    if count == 0 {
                        let mut qualifiers = vec![
                            ("label".into(), Some(plain(primer.attr("name")))),
                            (
                                "primer_sequence".into(),
                                Some(primer.attr("sequence").into()),
                            ),
                        ];
                        for (key, value) in &primer.attrs {
                            if key != "name" && key != "sequence" {
                                qualifiers.push((format!("snapgene_{key}"), Some(value.clone())));
                            }
                        }
                        doc.features.push(Feature {
                            location: None,
                            id: doc.features.len() as u64 + 1,
                            label: plain(primer.attr("name")),
                            kind: "primer_bind".into(),
                            parts: vec![],
                            reverse: false,
                            qualifiers,
                            edited: false,
                        });
                    }
                }
            }
            6 => {
                let root = xml(body)?;
                if root.name != "Notes" {
                    return Err(bad("Invalid notes root"));
                }
                doc.references = references::read(&root);
                if let Some(v) = root.child("Description") {
                    doc.description = plain(&v.text);
                }
                if let Some(v) = root.child("CustomMapLabel") {
                    let label = plain(&v.text);
                    if !label.is_empty() {
                        doc.name = label;
                    }
                }
                if let Some(v) = root.child("AccessionNumber") {
                    doc.accession = Some(plain(&v.text));
                }
                if let Some(v) = root.child("Comments") {
                    doc.comments.push(plain(&v.text));
                }
            }
            _ => {}
        }
        if doc.features.len() > 20_000 {
            return Err(bad("More than 20,000 annotations"));
        }
    }
    colors::read(&mut doc, &blocks)?;
    ends::read(&mut doc, &blocks)?;
    if unsupported_segments != 0 {
        doc.warnings
            .push(segment_types::import_warning(unsupported_segments));
    }
    doc.warnings.push("SnapGene import: DNA/RNA/protein, topology, strandedness, exact feature ranges, strand, qualifiers, primer binding annotations and selected notes are interpreted. An unchanged native save retains the original bytes. Saving edits rebuilds sequence and annotations, preserves reference/Notes fields and embedded files, and archives sequence-dependent source records in recoverable native history when the molecule changes.".into());
    doc.warnings.push("Primer binding sites are imported without recalculating SnapGene hybridization filters; duplicate simplified sites are collapsed. Imported melting temperatures are source values, not PopGene calculations.".into());
    doc.source = Some(Source {
        file_name: name.into(),
        format: "SnapGene".into(),
        bytes: bytes.into(),
        retained_blocks: blocks.iter().map(|(k, b)| (*k, b.len())).collect(),
    });
    super::lineage::read(&mut doc, bytes, super::SaveFormat::SnapGene);
    doc.validate()?;
    Ok(doc)
}

/// Compare the supported model, not the presence of retained bytes or UI warnings.
/// Re-parsing the source prevents a source attachment from authorizing stale saves.
fn unchanged_source(doc: &Document) -> Option<&[u8]> {
    let source = doc.source.as_ref().filter(|s| s.format == "SnapGene")?;
    let mut original = read(&source.file_name, &source.bytes).ok()?;
    let mut current = doc.clone();
    for d in [&mut original, &mut current] {
        d.source = None;
        d.warnings.clear();
        for (index, feature) in d
            .features
            .iter_mut()
            .chain(&mut d.duplex_annotations.features)
            .enumerate()
        {
            feature.id = index as u64;
        }
    }
    (original == current).then_some(source.bytes.as_ref())
}

fn qualifier<'a>(feature: &'a Feature, name: &str) -> Option<&'a str> {
    feature
        .qualifiers
        .iter()
        .find_map(|(key, value)| (key == name).then_some(value.as_deref()).flatten())
}

// SnapGene 8.2.3 saved INSDC 6^7 as classification=1, range=6-7.
fn native_between(feature: &Feature) -> Option<(usize, usize, bool)> {
    match feature.location.as_ref()? {
        FeatureLocation::Between { left, right } => Some((*left, *right, false)),
        FeatureLocation::Complement(child) => {
            if let FeatureLocation::Between { left, right } = child.as_ref() {
                Some((*left, *right, true))
            } else {
                None
            }
        }
        _ => None,
    }
}

// Native segments represent an ordered concatenation on one strand. An explicit
// AST can carry that same meaning after workflow projection or nested joins;
// its presence alone does not require a different native field. Fuzzy bounds
// must be checked separately: uniform_projection deliberately permits them for
// display, whereas native standard segments would turn them into exact bounds.
fn native_exact_location(feature: &Feature, molecule: Molecule) -> bool {
    fn exact(location: &FeatureLocation, molecule: Molecule) -> bool {
        match location {
            FeatureLocation::Range { before, after, .. } => !before && !after,
            FeatureLocation::Complement(child) => {
                molecule != Molecule::Protein && exact(child, molecule)
            }
            FeatureLocation::Join(children) => children.iter().all(|child| exact(child, molecule)),
            _ => false,
        }
    }
    feature.location.as_ref().is_some_and(|location| {
        exact(location, molecule)
            && location
                .uniform_projection()
                .is_some_and(|(parts, reverse)| {
                    parts == feature.parts && reverse == feature.reverse
                })
    })
}

/// An oligo is written as a native Primer only when its exact full-length binding
/// can be established independently. Mismatches, tails and discontinuous sites
/// remain ordinary primer_bind annotations with their original qualifiers.
fn native_primer<'a>(feature: &'a Feature, doc: &Document) -> Option<&'a str> {
    if doc.molecule != Molecule::Dna || feature.kind != "primer_bind" || feature.edited {
        return None;
    }
    let sequence = qualifier(feature, "primer_sequence")?;
    if sequence.is_empty()
        || !sequence.bytes().all(|c| b"ACGTRYSWKMBDHVN".contains(&c))
        || feature
            .parts
            .iter()
            .any(|p| p.start >= p.end || p.end > doc.sequence.len())
        || (!feature.parts.is_empty() && sequence != feature.extract(&doc.sequence))
    {
        return None;
    }
    match feature.parts.as_slice() {
        [] | [_] => {}
        [a, b] if doc.circular && a.end == doc.sequence.len() && b.start == 0 => {}
        _ => return None,
    }
    // Native primer XML has no general qualifier bag. Keep any additional user
    // annotations as features instead of silently dropping them.
    feature
        .qualifiers
        .iter()
        .all(|(key, _)| {
            matches!(key.as_str(), "label" | "primer_sequence")
                || key.starts_with("snapgene_")
                || key.starts_with("binding_")
        })
        .then_some(sequence)
}

pub(crate) fn write_warnings(doc: &Document) -> Vec<String> {
    if unchanged_source(doc).is_some() {
        return vec![];
    }
    let original = doc
        .source
        .as_ref()
        .filter(|source| source.format == "SnapGene")
        .and_then(|source| read(&source.file_name, &source.bytes).ok());
    let same_molecule = original
        .as_ref()
        .is_some_and(|original| super::snapgene_preserve::same_molecule(doc, original));
    let source_packets = original
        .as_ref()
        .and(doc.source.as_ref())
        .and_then(|source| super::snapgene_preserve::packets(&source.bytes).ok())
        .unwrap_or_default();
    let preserved = original
        .as_ref()
        .and_then(|original| {
            record_preservation::unchanged_records(doc, original, &source_packets).ok()
        })
        .unwrap_or_default();
    let mut warnings = vec![];
    let preserve_annotations = same_molecule
        && original
            .as_ref()
            .is_some_and(|original| record_preservation::same_annotations(doc, original));
    if !preserve_annotations
        && let Ok((native, _)) = ends::encode(doc)
        && primer_envelope::plan(&native, &preserved, &source_packets, same_molecule).is_err()
    {
        warnings.push(primer_envelope::REFUSAL.into());
    }
    if segment_types::requires_rebuild(doc, original.as_ref(), &preserved) {
        warnings.push(segment_types::REBUILD_REFUSAL.into());
    }
    if doc.molecule == Molecule::Dna && doc.methylation.is_none() {
        warnings.push("Current whole-template methylation flags are unknown. Native DNA has only binary Dam/Dcm/EcoKI flags, so the exported flags are clear; unknown or partial methylation cannot be represented. Reopening those flags does not establish measured absence of methylation.".into());
    }
    if doc.sequence_qualities.is_some() || doc.opposite_sequence_qualities.is_some() {
        warnings.push("Native sequence export does not encode per-residue top- or opposite-strand quality scores; save a separate quality-bearing read for those values.".into());
    }
    if doc
        .features
        .iter()
        .chain(&doc.duplex_annotations.features)
        .any(|feature| {
            !preserved.contains(feature.id)
                && feature
                    .display_segments()
                    .iter()
                    .any(|segment| segment.style == popgene_core::FeatureSegmentStyle::Line)
        })
    {
        warnings.push("Native conversion of line appearance is refused until its native encoding is verified; GenBank/EMBL retains this style. Unchanged native annotation packets can be preserved during unrelated metadata edits.".into());
    }
    if doc.features.iter().any(|f| f.edited) {
        warnings.push("Annotations changed by sequence edits carry a review note; the in-memory review flag is not stored. Primers requiring review are saved as annotated primer_bind features.".into());
    }
    if doc.features.iter().any(|feature| {
        feature.kind == "CDS" && popgene_core::translation::cds_phase(feature).is_err()
    }) {
        warnings.push("Invalid or conflicting CDS readingFrame/codon_start values remain qualifiers; the native readingFrame attribute is omitted. Review these CDS offsets before using their translations.".into());
    }
    if doc.source.as_ref().is_some_and(|s| s.format == "SnapGene") {
        warnings.push(if same_molecule {
            "Unchanged native feature and primer records retain their full source attributes, binding details and chemistry. Changed annotations are rebuilt from supported fields. Original references, Notes fields, attachments, history, alignments and other source packets remain because the molecule is unchanged. Sequence letter case is normalized.".into()
        } else {
            "Original references, Notes fields and embedded files are retained. A consolidated replace step preserves the original molecule and its records in native history, and physical ends are rebuilt from the current strands. Source alignments and display records remain historical. Original methylation/other sequence flags are recorded as historical Notes, not reconstructed as current chemistry.".into()
        });
    }
    if doc
        .features
        .iter()
        .chain(&doc.duplex_annotations.features)
        .any(|f| !preserved.contains(f.id) && native_primer(f, doc).is_some())
    {
        warnings.push("Exact full-length primer matches are saved as native primers. For rebuilt primers, source melting temperatures, hybridization settings and detailed oligo chemistry are not copied; SnapGene may recalculate binding sites. Unchanged native primers retain those source fields.".into());
    }
    if doc.features.iter().any(|f| {
        !preserved.contains(f.id) && f.kind == "primer_bind" && native_primer(f, doc).is_none()
    }) {
        warnings.push("Primer sites without a verified full-length exact oligo match, or with additional qualifiers, are saved as primer_bind features. Their sequences and qualifiers remain annotations rather than native SnapGene primers.".into());
    }
    if doc.trace.is_some() {
        warnings.push("The .dna file includes called bases but not chromatogram channels, peaks or quality scores. Keep the original trace file.".into());
    }
    if references::supplemental(&doc.references).is_some() {
        warnings.push("Reference consortium, scope, remarks and unsupported bibliographic attributes are preserved as readable comments; native reference fields cannot retain their separate structure.".into());
    }
    if doc.comments.len() > 1 {
        warnings.push("Comments are combined into SnapGene's single comments field; their text is preserved, but separate comment boundaries are not.".into());
    }
    warnings
}

fn escaped(value: &str) -> String {
    quick_xml::escape::escape(value)
        .replace('\n', "&#10;")
        .replace('\r', "&#13;")
        .replace('\t', "&#9;")
}

fn serialize_node(node: &Node, out: &mut String) {
    out.push('<');
    out.push_str(&node.name);
    // Native history readers historically use this attribute order. Preserve
    // the observed convention rather than relying on XML attribute order freedom.
    const HISTORY_ATTRS: &[&str] = &[
        "name",
        "type",
        "seqLen",
        "strandedness",
        "ID",
        "circular",
        "upstreamStickiness",
        "downstreamStickiness",
        "upstreamModification",
        "downstreamModification",
        "resurrectable",
        "operation",
    ];
    if node.name == "Node" {
        for name in HISTORY_ATTRS {
            if let Some(value) = node.attrs.get(*name) {
                out.push_str(&format!(" {name}=\"{}\"", escaped(value)));
            }
        }
    }
    for (name, value) in &node.attrs {
        if node.name != "Node" || !HISTORY_ATTRS.contains(&name.as_str()) {
            out.push_str(&format!(" {name}=\"{}\"", escaped(value)));
        }
    }
    if node.text.is_empty() && node.children.is_empty() {
        out.push_str("/>");
        return;
    }
    out.push('>');
    out.push_str(&escaped(&node.text));
    for child in &node.children {
        serialize_node(child, out);
    }
    out.push_str("</");
    out.push_str(&node.name);
    out.push('>');
}

fn packet(out: &mut Vec<u8>, kind: u8, body: &[u8]) -> Result<()> {
    let size = u32::try_from(body.len())
        .map_err(|_| Error::Message("SnapGene export: packet exceeds 4 GiB".into()))?;
    if matches!(kind, 5 | 6 | 10) {
        // Apply the same size and structure limits to generated and imported XML.
        xml(body)?;
    }
    if out
        .len()
        .checked_add(5)
        .and_then(|n| n.checked_add(body.len()))
        .is_none_or(|n| n as u64 > crate::MAX_FILE)
    {
        return Err(bad("Generated native file exceeds 128 MiB"));
    }
    out.push(kind);
    out.extend_from_slice(&size.to_be_bytes());
    out.extend_from_slice(body);
    Ok(())
}

/// Write conventional SnapGene TLV records, never private packet types.
/// The format flags are documented by SnapGeneReader; feature/primer coordinates
/// and reverse-strand compound order are independently covered by Biopython.
pub(crate) fn write(doc: &Document) -> Result<Vec<u8>> {
    if let Some(original) = unchanged_source(doc) {
        return Ok(original.to_vec());
    }
    super::lineage::write(doc, super::SaveFormat::SnapGene, write_plain)
}

/// Current model IDs in the exact order in which native records decode. This
/// uses the same emission plan as the writer, including retained multi-site
/// Primer expansion. The reader partitions axes without changing those IDs.
pub(super) fn encoded_feature_ids(doc: &Document) -> Result<Vec<u64>> {
    let original = doc
        .source
        .as_ref()
        .filter(|s| s.format == "SnapGene")
        .and_then(|s| read(&s.file_name, &s.bytes).ok());
    let source_packets = original
        .as_ref()
        .and(doc.source.as_ref())
        .map(|s| super::snapgene_preserve::packets(&s.bytes))
        .transpose()?
        .unwrap_or_default();
    if let Some(original) = &original
        && super::snapgene_preserve::same_molecule(doc, original)
        && record_preservation::same_annotations(doc, original)
    {
        let mut pairs = original
            .features
            .iter()
            .chain(&original.duplex_annotations.features)
            .zip(doc.features.iter().chain(&doc.duplex_annotations.features))
            .map(|(old, current)| (old.id, current.id))
            .collect::<Vec<_>>();
        pairs.sort_by_key(|(old_id, _)| *old_id);
        return Ok(pairs.into_iter().map(|(_, current)| current).collect());
    }
    let records = original
        .as_ref()
        .map(|original| record_preservation::unchanged_records(doc, original, &source_packets))
        .transpose()?
        .unwrap_or_default();
    let (native, _) = ends::encode(doc)?;
    let plan = annotation_order::plan(&native, &records);
    Ok(plan
        .iter()
        .filter(|item| !item.is_primer(&records))
        .chain(plan.iter().filter(|item| item.is_primer(&records)))
        .flat_map(|item| item.model_ids(&native, &records).iter().copied())
        .collect())
}

pub(super) fn write_plain(doc: &Document) -> Result<Vec<u8>> {
    doc.validate()?;
    if doc.sequence.is_empty() {
        return Err(Error::Message(
            "SnapGene export requires a non-empty sequence.".into(),
        ));
    }
    if let Some(original) = unchanged_source(doc) {
        return Ok(original.to_vec());
    }
    if doc
        .features
        .iter()
        .chain(&doc.duplex_annotations.features)
        .any(|feature| {
            feature.location.is_some()
                && native_between(feature).is_none()
                && !native_exact_location(feature, doc.molecule)
        })
    {
        return Err(Error::Message("Native SnapGene export cannot currently represent this fuzzy, remote or nonuniform compound location without changing its meaning. Use GenBank or EMBL to preserve it.".into()));
    }
    if doc
        .features
        .iter()
        .chain(&doc.duplex_annotations.features)
        .any(|feature| {
            feature.parts.is_empty()
                && native_primer(feature, doc).is_none()
                && native_between(feature).is_none()
        })
    {
        return Err(Error::Message("Unlocated native records require a valid independent primer sequence with native-compatible qualifiers.".into()));
    }
    if doc.features.len() + doc.duplex_annotations.features.len() > 20_000 {
        return Err(Error::Message(
            "SnapGene export: more than 20,000 annotations".into(),
        ));
    }
    let texts = std::iter::once(doc.name.as_str())
        .chain(std::iter::once(doc.description.as_str()))
        .chain(doc.accession.as_deref())
        .chain(doc.comments.iter().map(String::as_str))
        .chain(
            doc.features
                .iter()
                .chain(&doc.duplex_annotations.features)
                .flat_map(|f| {
                    std::iter::once(f.label.as_str())
                        .chain(std::iter::once(f.kind.as_str()))
                        .chain(f.qualifiers.iter().flat_map(|(key, value)| {
                            std::iter::once(key.as_str()).chain(value.as_deref())
                        }))
                }),
        );
    if texts.into_iter().any(|text| text.chars().any(|c| {
        !matches!(c, '\t' | '\n' | '\r' | '\u{20}'..='\u{d7ff}' | '\u{e000}'..='\u{fffd}' | '\u{10000}'..='\u{10ffff}')
    })) {
        return Err(Error::Message("SnapGene export: metadata contains characters that XML cannot represent".into()));
    }
    if doc
        .features
        .iter()
        .chain(&doc.duplex_annotations.features)
        .any(|f| f.qualifiers.iter().any(|(key, _)| key.is_empty()))
    {
        return Err(Error::Message(
            "SnapGene export: a feature qualifier has an empty name".into(),
        ));
    }
    let original = doc
        .source
        .as_ref()
        .filter(|source| source.format == "SnapGene")
        .and_then(|source| read(&source.file_name, &source.bytes).ok());
    let source_packets = original
        .as_ref()
        .and(doc.source.as_ref())
        .map(|source| super::snapgene_preserve::packets(&source.bytes))
        .transpose()?
        .unwrap_or_default();
    let same_molecule = original
        .as_ref()
        .is_some_and(|original| super::snapgene_preserve::same_molecule(doc, original));
    let preserve_annotations = original.as_ref().is_some_and(|original| {
        super::snapgene_preserve::same_molecule(doc, original)
            && record_preservation::same_annotations(doc, original)
    });
    let mut preserved_records = original
        .as_ref()
        .map(|original| record_preservation::unchanged_records(doc, original, &source_packets))
        .transpose()?
        .unwrap_or_default();
    if segment_types::requires_rebuild(doc, original.as_ref(), &preserved_records) {
        return Err(Error::Message(segment_types::REBUILD_REFUSAL.into()));
    }
    if !preserve_annotations
        && doc
            .features
            .iter()
            .chain(&doc.duplex_annotations.features)
            .any(|feature| {
                !preserved_records.contains(feature.id)
                    && feature
                        .display_segments()
                        .iter()
                        .any(|segment| segment.style == popgene_core::FeatureSegmentStyle::Line)
            })
    {
        return Err(bad(
            "Native line appearance encoding has not been verified. Save GenBank or EMBL to preserve this display style; unrelated metadata-only native edits retain the original annotation packets.",
        ));
    }
    let (native, end_properties) = ends::encode(doc)?;
    let primer_envelope = if preserve_annotations {
        None
    } else {
        Some(primer_envelope::plan(
            &native,
            &preserved_records,
            &source_packets,
            same_molecule,
        )?)
    };
    if let Some(envelope) = &primer_envelope {
        preserved_records.next_primer_id = preserved_records
            .next_primer_id
            .max(envelope.minimum_next_id());
    }
    let mut out = Vec::with_capacity(native.sequence.len() + 1024);
    // Export version 15 / import version 20 are used by the installed 8.2.3
    // synthetic fixtures. This writer emits only uncompressed, established blocks.
    let (native_type, sequence_kind) = match doc.molecule {
        Molecule::Dna => (1_u16, 0),
        Molecule::Rna => (7, 32),
        Molecule::Protein => (2, 21),
    };
    let mut cookie = b"SnapGene".to_vec();
    cookie.extend_from_slice(&native_type.to_be_bytes());
    cookie.extend_from_slice(&[0, 15, 0, 20]);
    packet(&mut out, 9, &cookie)?;
    let mut sequence = Vec::with_capacity(doc.sequence.len() + 1);
    let preserved_flags = if same_molecule {
        source_packets
            .iter()
            .find(|(kind, _)| *kind == sequence_kind)
            .and_then(|(_, body)| body.first())
            .copied()
            .unwrap_or(0)
            & if doc.molecule == Molecule::Dna {
                !31
            } else {
                !3
            }
    } else {
        0
    };
    let methylation_flags = doc.methylation.map_or(0, |m| {
        (u8::from(m.dam) << 2) | (u8::from(m.dcm) << 3) | (u8::from(m.eco_ki) << 4)
    });
    sequence.push(
        preserved_flags
            | methylation_flags
            | if doc.molecule == Molecule::Protein {
                0
            } else {
                u8::from(doc.circular) | if doc.single_stranded { 0 } else { 2 }
            },
    );
    sequence.extend_from_slice(native.sequence.as_bytes());
    packet(&mut out, sequence_kind, &sequence)?;

    let mut features = String::new();
    let mut primers = String::new();
    let mut primer_id = 0;
    for emission in annotation_order::plan(&native, &preserved_records) {
        let feature_index = emission.feature_index;
        let feature = &native.features[feature_index];
        let effective_color = feature.effective_color_rgb(feature_index);
        if let annotation_order::Payload::Preserved(index) = emission.payload {
            let record = preserved_records.record(index);
            let (kind, mut node) = (record.kind, record.node.clone());
            if kind == 10 {
                node.attrs.insert(
                    "recentID".into(),
                    preserved_records.feature_id(feature.id)?.to_string(),
                );
                serialize_node(&node, &mut features);
            } else {
                node.attrs.insert(
                    "recentID".into(),
                    preserved_records.primer_id(feature.id)?.to_string(),
                );
                serialize_node(&node, &mut primers);
                primer_id += 1;
            }
            continue;
        }

        if let annotation_order::Payload::Primer(sequence) = emission.payload {
            let native_primer_id = preserved_records.primer_id(feature.id)?;
            primers.push_str(&format!(
                "<Primer recentID=\"{native_primer_id}\" name=\"{}\" sequence=\"{}\"",
                escaped(&feature.label),
                sequence
            ));
            for (key, target) in [
                ("snapgene_color", "color"),
                ("snapgene_dateAdded", "dateAdded"),
            ] {
                if let Some(value) = qualifier(feature, key) {
                    primers.push_str(&format!(" {target}=\"{}\"", escaped(value)));
                }
            }
            if let (Some(first), Some(last)) = (feature.parts.first(), feature.parts.last()) {
                let range = format!("{}-{}", first.start, last.end - 1);
                primers.push_str(&format!("><BindingSite location=\"{range}\" boundStrand=\"{}\" annealedBases=\"{sequence}\"><Component hybridizedRange=\"{range}\" bases=\"{sequence}\"/></BindingSite></Primer>", u8::from(feature.reverse)));
            } else {
                primers.push_str("/>");
            }
            primer_id += 1;
            continue;
        }
        let direction = if native_between(feature).is_some_and(|(_, _, reverse)| reverse) {
            "2"
        } else if doc.molecule == Molecule::Protein {
            "0"
        } else if feature.reverse {
            "2"
        } else {
            match feature.display_direction() {
                Some(FeatureDirection::Nondirectional) => "0",
                Some(FeatureDirection::Bidirectional) => "3",
                _ => "1",
            }
        };
        let reading_frame = (doc.molecule != Molecule::Protein && feature.kind == "CDS")
            .then(|| popgene_core::translation::cds_phase(feature).ok())
            .flatten()
            .map(|phase| (phase as i8 + 1) * if feature.reverse { -1 } else { 1 });
        let native_feature_id = preserved_records.feature_id(feature.id)?;
        features.push_str(&format!("<Feature recentID=\"{native_feature_id}\" name=\"{}\" type=\"{}\" directionality=\"{direction}\" allowSegmentOverlaps=\"1\"", escaped(&feature.label), escaped(&feature.kind)));
        if native_between(feature).is_some() {
            features.push_str(" classification=\"1\"");
        }
        if let Some(frame) = reading_frame {
            features.push_str(&format!(" readingFrame=\"{frame}\""));
        }
        let cleavage_sites = feature.display_cleavage_sites();
        if !cleavage_sites.is_empty() {
            features.push_str(&format!(
                " cleavageArrows=\"{}\"",
                cleavage_sites
                    .iter()
                    .map(usize::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
            ));
        }
        features.push('>');
        let display_segments = feature.display_segments();
        if let Some((left, right, _)) = native_between(feature) {
            features.push_str(&format!(
                "<Segment range=\"{}-{}\" type=\"standard\" color=\"#{:06x}\"/>",
                left + 1,
                right + 1,
                effective_color
            ));
        }
        // Keep forward coordinate order even on reverse features: the strand
        // applies to the whole concatenation, not to each segment independently.
        for display_segment in &display_segments {
            let part = &display_segment.range;
            features.push_str(&format!(
                "<Segment range=\"{}-{}\" type=\"standard\"",
                part.start + 1,
                part.end
            ));
            if let Some(color) = display_segment.color_rgb {
                features.push_str(&format!(" color=\"#{color:06x}\""));
            } else {
                features.push_str(&format!(" color=\"#{effective_color:06x}\""));
            }
            if let Some(name) = &display_segment.name {
                features.push_str(&format!(" name=\"{}\"", escaped(name)));
            }
            if feature.kind == "CDS" {
                features.push_str(" translated=\"1\"");
            }
            features.push_str("/>");
        }
        for gap in feature.display_gaps() {
            features.push_str(&format!(
                "<Segment range=\"{}-{}\" type=\"gap\"",
                gap.range.start + 1,
                gap.range.end
            ));
            if let Some(color) = gap.color_rgb {
                features.push_str(&format!(" color=\"#{color:06x}\""));
            }
            if let Some(name) = &gap.name {
                features.push_str(&format!(" name=\"{}\"", escaped(name)));
            }
            features.push_str("/>");
        }
        // Group equal keys so readers that use a dictionary do not overwrite
        // earlier values. Display metadata remains in standard Q records too.
        let mut qualifiers: BTreeMap<&str, Vec<Option<&str>>> = BTreeMap::new();
        for (key, value) in &feature.qualifiers {
            if key != "snapgene_color"
                && key != "snapgene_gap"
                && key != "snapgene_directionality"
                && !(key == "snapgene_readingFrame" && reading_frame.is_some())
                && !(key == "segment_name"
                    && display_segments
                        .iter()
                        .any(|segment| segment.name.is_some()))
            {
                qualifiers.entry(key).or_default().push(value.as_deref());
            }
        }
        if feature.edited {
            let notes = qualifiers.entry("note").or_default();
            if !notes.contains(&Some(REVIEW_NOTE)) {
                notes.push(Some(REVIEW_NOTE));
            }
        }
        for (key, values) in qualifiers {
            features.push_str(&format!("<Q name=\"{}\">", escaped(key)));
            for value in values {
                if let Some(value) = value {
                    let attr = if matches!(key, "codon_start" | "transl_table")
                        && value.parse::<u32>().is_ok()
                    {
                        "int"
                    } else {
                        "text"
                    };
                    features.push_str(&format!("<V {attr}=\"{}\"/>", escaped(value)));
                } else {
                    features.push_str("<V/>");
                }
            }
            features.push_str("</Q>");
        }
        features.push_str("</Feature>");
    }
    features = format!(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?><Features nextValidID=\"{}\">{features}</Features>",
        preserved_records.next_feature_id
    );
    if preserve_annotations {
        for (kind, body) in &source_packets {
            if matches!(kind, 5 | 10) {
                packet(&mut out, *kind, body)?;
            }
        }
    } else {
        packet(&mut out, 10, features.as_bytes())?;
        if let Some(envelope) = primer_envelope {
            envelope.emit(
                &mut out,
                &primers,
                primer_id,
                preserved_records.next_primer_id,
            )?;
        }
    }
    // Merge modeled values into original Notes, preserving reference subtrees,
    // sequence class, authorship, dates and other independent source metadata.
    let notes_source = source_packets
        .iter()
        .find(|(kind, _)| *kind == 6)
        .map(|(_, body)| xml(body))
        .transpose()?;
    let mut notes = notes_source.unwrap_or_else(|| Node {
        name: "Notes".into(),
        ..Default::default()
    });
    if original
        .as_ref()
        .is_none_or(|source| source.references != doc.references)
    {
        notes.children.retain(|child| child.name != "References");
        if !doc.references.is_empty() {
            notes.children.push(references::write(&doc.references));
        }
    }
    let mut comments = doc.comments.clone();
    if let Some(supplemental) = references::supplemental(&doc.references) {
        comments.push(supplemental);
    }
    for (key, value) in [
        ("CustomMapLabel", Some(doc.name.clone())),
        ("Description", Some(doc.description.clone())),
        ("AccessionNumber", doc.accession.clone()),
        (
            "Comments",
            (!comments.is_empty()).then(|| comments.join("\n")),
        ),
    ] {
        notes.children.retain(|child| child.name != key);
        if let Some(text) = value {
            notes.children.push(Node {
                name: key.into(),
                text,
                ..Default::default()
            });
        }
    }
    let mut serialized = String::from("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
    serialize_node(&notes, &mut serialized);
    packet(&mut out, 6, serialized.as_bytes())?;
    if let Some(properties) = end_properties {
        let mut serialized = String::new();
        serialize_node(&properties, &mut serialized);
        packet(&mut out, 8, serialized.as_bytes())?;
    }
    if !same_molecule && let (Some(original), Some(source)) = (&original, &doc.source) {
        for (kind, body) in history::append(doc, original, source, &native)? {
            packet(&mut out, kind, &body)?;
        }
    }
    let preserve_colors = same_molecule
        && original.as_ref().is_some_and(|original| {
            original.sequence_colors == doc.sequence_colors
                && original.duplex_annotations.sequence_colors
                    == doc.duplex_annotations.sequence_colors
        });
    if !preserve_colors && !native.sequence_colors.is_empty() {
        let mut serialized = String::new();
        serialize_node(&colors::write(&native), &mut serialized);
        packet(&mut out, 20, serialized.as_bytes())?;
    }
    for (kind, body) in source_packets {
        // Attachments have no sequence coordinates. Every other unmodeled block
        // remains exact only if sequence, molecule type and topology agree.
        if kind == 23
            || (same_molecule
                && !matches!(kind, 0 | 1 | 5 | 6 | 9 | 10 | 21 | 32)
                && !(kind == 8 && doc.duplex_ends.is_some())
                && (kind != 20 || preserve_colors))
        {
            packet(&mut out, kind, body)?;
        }
    }

    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    const DNA: &[u8] = include_bytes!("../../../validation/fixtures/annotated.dna");
    #[test]
    fn recorded_methylation_flags_import_write_and_change_independently() {
        use popgene_core::metadata::DnaMethylation;
        for flags in 0_u8..8 {
            let recorded = DnaMethylation {
                dam: flags & 1 != 0,
                dcm: flags & 2 != 0,
                eco_ki: flags & 4 != 0,
            };
            let mut doc = Document::new("recorded methylation", "GATCCAGCTG").unwrap();
            doc.methylation = Some(recorded);
            let bytes = write(&doc).unwrap();
            let blocks = crate::snapgene_preserve::packets(&bytes).unwrap();
            assert_eq!(
                blocks.iter().find(|(kind, _)| *kind == 0).unwrap().1[0],
                2 | flags << 2
            );
            let mut restored = read("methylation.dna", &bytes).unwrap();
            assert_eq!(restored.methylation, Some(recorded));
            assert_eq!(write(&restored).unwrap(), bytes);
            restored.methylation.as_mut().unwrap().dam = !recorded.dam;
            let changed = write(&restored).unwrap();
            assert_ne!(changed, bytes);
            let changed = read("changed-methylation.dna", &changed).unwrap();
            assert_eq!(changed.methylation.unwrap().dam, !recorded.dam);
            let historical = native_history_document(changed.source.as_ref().unwrap(), 0).unwrap();
            assert_eq!(historical.methylation, None);
            assert!(
                historical
                    .warnings
                    .iter()
                    .any(|w| w.contains("methylation assumptions are unknown"))
            );
        }
    }
    #[test]
    fn metadata_edits_preserve_sequence_flags_and_changed_molecules_archive_their_values() {
        let original = Document::new("flag provenance", "AAACCCGGGTTT").unwrap();
        let mut bytes = vec![];
        for (kind, body) in crate::snapgene_preserve::packets(&write(&original).unwrap()).unwrap() {
            let mut body = body.to_vec();
            if kind == 0 {
                body[0] |= 0xfc;
            }
            packet(&mut bytes, kind, &body).unwrap();
        }
        let mut doc = read("flags.dna", &bytes).unwrap();
        assert_eq!(
            doc.methylation,
            Some(popgene_core::metadata::DnaMethylation {
                dam: true,
                dcm: true,
                eco_ki: true
            })
        );
        assert!(source_metadata(&bytes).unwrap().iter().any(|(_, value)| {
            value.contains("Raw sequence flags: 0xfe")
                && value.contains("Dam whole-template flag: set")
                && value.contains("Dcm whole-template flag: set")
                && value.contains("EcoKI whole-template flag: set")
        }));
        doc.description = "Edited description".into();
        let written = write(&doc).unwrap();
        let blocks = crate::snapgene_preserve::packets(&written).unwrap();
        assert_eq!(
            blocks.iter().find(|(kind, _)| *kind == 0).unwrap().1[0],
            0xfe
        );
        doc.replace(0, 1, "T").unwrap();
        assert_eq!(doc.methylation, None);
        assert!(
            write_warnings(&doc)
                .iter()
                .any(|w| w.contains("methylation flags are unknown"))
        );
        let written = write(&doc).unwrap();
        let blocks = crate::snapgene_preserve::packets(&written).unwrap();
        assert_eq!(blocks.iter().find(|(kind, _)| *kind == 0).unwrap().1[0], 2);
        let edited = read("edited-flags.dna", &written).unwrap();
        let restored = native_history_document(edited.source.as_ref().unwrap(), 0).unwrap();
        assert_eq!(restored.sequence, original.sequence);
        assert_eq!(restored.methylation, None);
        assert!(restored.comments.iter().any(|comment| {
            comment.contains("historical provenance only")
                && comment.contains("Raw sequence flags: 0xfe")
                && comment.contains("Dam whole-template flag: set")
        }));
    }
    #[test]
    fn observed_native_rna_and_protein_layouts_roundtrip() {
        // Minimal synthetic packets using the independently observed 8.2.3
        // .rna/.prot cookie, sequence type and flags (not vendor-produced bytes).
        for (molecule, native_type, kind, flags, sequence, extension) in [
            (Molecule::Rna, 7_u16, 32, 3, "acguaugcuuagnnry", "rna"),
            (Molecule::Rna, 7, 32, 0, "augcuu", "rna"),
            (
                Molecule::Protein,
                2,
                21,
                0,
                "mkwvtfisllflfssaysrg*x",
                "prot",
            ),
        ] {
            let mut bytes = vec![];
            let mut cookie = b"SnapGene".to_vec();
            cookie.extend_from_slice(&native_type.to_be_bytes());
            cookie.extend_from_slice(&[0, 15, 0, 20]);
            packet(&mut bytes, 9, &cookie).unwrap();
            let mut payload = vec![flags];
            payload.extend_from_slice(sequence.as_bytes());
            packet(&mut bytes, kind, &payload).unwrap();
            packet(&mut bytes, 10, br##"<Features><Feature name="Synthetic domain" type="misc_feature"><Segment range="2-5" color="#a6acb3" type="standard"/><Q name="note"><V text="Synthetic metadata"/></Q></Feature></Features>"##).unwrap();
            let mut document = read(&format!("native.{extension}"), &bytes).unwrap();
            assert_eq!(document.molecule, molecule);
            assert_eq!(document.sequence.as_ref(), sequence.to_ascii_uppercase());
            assert_eq!(document.circular, flags & 1 != 0);
            assert_eq!(
                document.single_stranded,
                molecule != Molecule::Protein && flags & 2 == 0
            );
            assert_eq!(
                document.features[0].parts,
                vec![Segment { start: 1, end: 5 }]
            );
            assert_eq!(write(&document).unwrap(), bytes);
            document.description = "Updated description".into();
            let encoded = write(&document).unwrap();
            assert_eq!(encoded[14], native_type as u8);
            assert_eq!(encoded[19], kind);
            let rebuilt = read(&format!("edited.{extension}"), &encoded).unwrap();
            assert_eq!(rebuilt.molecule, molecule);
            assert_eq!(rebuilt.sequence, document.sequence);
            assert_eq!(rebuilt.circular, document.circular);
            assert_eq!(rebuilt.single_stranded, document.single_stranded);
            assert_eq!(
                rebuilt.features[0].qualifier("note"),
                Some("Synthetic metadata")
            );
            assert_eq!(rebuilt.description, "Updated description");
            assert_eq!(
                crate::SaveFormat::SnapGene.extension_for(&document),
                extension
            );
            let directory = tempfile::tempdir().unwrap();
            assert!(
                crate::save(
                    &directory.path().join("wrong.dna"),
                    &document,
                    crate::SaveFormat::SnapGene,
                    None
                )
                .is_err()
            );
            crate::save(
                &directory.path().join(format!("correct.{extension}")),
                &document,
                crate::SaveFormat::SnapGene,
                None,
            )
            .unwrap();
        }
    }
    #[test]
    fn native_cookie_and_sequence_blocks_must_agree() {
        for (native_type, kind, flags, sequence) in [
            (7_u16, 0, 0, "ACGT"),
            (2, 32, 0, "ACGU"),
            (2, 21, 3, "MKW"),
            (7, 32, 0, "ACGT"),
            (99, 0, 2, "ACGT"),
        ] {
            let mut bytes = vec![];
            let mut cookie = b"SnapGene".to_vec();
            cookie.extend_from_slice(&native_type.to_be_bytes());
            cookie.extend_from_slice(&[0, 15, 0, 20]);
            packet(&mut bytes, 9, &cookie).unwrap();
            let mut data = vec![flags];
            data.extend_from_slice(sequence.as_bytes());
            packet(&mut bytes, kind, &data).unwrap();
            assert!(read("native", &bytes).is_err());
        }
    }
    #[test]
    fn native_offset_fixture_matches_genbank_calculations_and_rebuilt_frames() {
        let gb = crate::import(
            "offsets.gb",
            include_bytes!("../../../validation/fixtures/regression/reading-frames.gb"),
        )
        .unwrap();
        let native_bytes =
            include_bytes!("../../../validation/fixtures/regression/reading-frames-native.dna");
        let mut native = read("offsets.dna", native_bytes).unwrap();
        let expected = [
            b"MK*PLFHT".as_slice(),
            b"*NNRYFI".as_slice(),
            b"EITVISY".as_slice(),
            b"V*NNGYF".as_slice(),
        ];
        for (index, expected) in expected.into_iter().enumerate() {
            assert_eq!(
                popgene_core::translation::translate_cds(&gb, &gb.features[index])
                    .unwrap()
                    .amino_acids,
                expected
            );
            assert_eq!(
                popgene_core::translation::translate_cds(&native, &native.features[index])
                    .unwrap()
                    .amino_acids,
                expected
            );
        }
        assert_eq!(
            native.features[1].qualifier("snapgene_readingFrame"),
            Some("2")
        );
        assert_eq!(
            native.features[2].qualifier("snapgene_readingFrame"),
            Some("3")
        );
        assert_eq!(
            native.features[3].qualifier("snapgene_readingFrame"),
            Some("-2")
        );
        assert_eq!(write(&native).unwrap(), native_bytes);
        let rebuilt = read("rebuilt.dna", &write(&gb).unwrap()).unwrap();
        for (index, expected) in ["1", "2", "3", "-2"].iter().enumerate() {
            assert_eq!(
                rebuilt.features[index].qualifier("snapgene_readingFrame"),
                Some(*expected)
            );
        }
        native.revcomp().unwrap();
        assert_eq!(
            native.features[1].qualifier("snapgene_readingFrame"),
            Some("-2")
        );
        assert_eq!(
            native.features[3].qualifier("snapgene_readingFrame"),
            Some("2")
        );
        for (index, expected) in expected.into_iter().enumerate() {
            assert_eq!(
                popgene_core::translation::translate_cds(&native, &native.features[index])
                    .unwrap()
                    .amino_acids,
                expected
            );
        }
        // Reading-frame-only input derives a conventional codon_start qualifier.
        let mut only_native = rebuilt;
        only_native.features[1]
            .qualifiers
            .retain(|(key, _)| key != "codon_start");
        let roundtrip = read("frame-only.dna", &write(&only_native).unwrap()).unwrap();
        assert_eq!(roundtrip.features[1].qualifier("codon_start"), Some("2"));
    }
    #[test]
    fn shifted_native_frames_are_feature_relative() {
        let gb = crate::import(
            "shifted.gb",
            include_bytes!("../../../validation/fixtures/regression/reading-frames-shifted.gb"),
        )
        .unwrap();
        let bytes = include_bytes!(
            "../../../validation/fixtures/regression/reading-frames-shifted-native.dna"
        );
        let mut native = read("shifted.dna", bytes).unwrap();
        let expected = ["EITVISY", "NNRYFI", "YEITVIS", "*NNGYF"];
        let frames = ["2", "3", "-2", "-3"];
        let assert_translations = |document: &Document| {
            for (feature, expected) in document.features.iter().zip(expected) {
                assert_eq!(
                    popgene_core::translation::translate_cds(document, feature)
                        .unwrap()
                        .amino_acids,
                    expected.as_bytes()
                );
            }
        };
        assert_translations(&gb);
        assert_translations(&native);
        for (feature, frame) in native.features.iter().zip(frames) {
            assert_eq!(feature.qualifier("snapgene_readingFrame"), Some(frame));
        }
        assert_eq!(write(&native).unwrap(), bytes);
        // A non-triplet insertion before every CDS changes all absolute frames
        // and document length, but must preserve the CDS-relative phase.
        native.replace(0, 0, "C").unwrap();
        assert_eq!(native.sequence.len(), 25);
        assert_translations(&native);
        let shifted = read("edited.dna", &write(&native).unwrap()).unwrap();
        assert_translations(&shifted);
        for (feature, frame) in shifted.features.iter().zip(frames) {
            assert_eq!(feature.qualifier("snapgene_readingFrame"), Some(frame));
        }
        native.revcomp().unwrap();
        assert_translations(&native);
        let reversed = read("reversed.dna", &write(&native).unwrap()).unwrap();
        assert_translations(&reversed);
        for (feature, frame) in reversed.features.iter().zip(["-2", "-3", "2", "3"]) {
            assert_eq!(feature.qualifier("snapgene_readingFrame"), Some(frame));
        }
    }
    #[test]
    fn native_reading_frame_attribute_is_preserved_without_inventing_codon_start() {
        let mut bytes = vec![];
        packet(&mut bytes, 9, b"SnapGene\x00\x01\x00\x0f\x00\x14").unwrap();
        packet(&mut bytes, 0, b"\x02AATGAAATAA").unwrap();
        packet(&mut bytes, 10, br##"<Features><Feature name="frame" type="CDS" directionality="1" readingFrame="1"><Segment range="1-10" type="standard" color="#993366"/><Q name="snapgene_readingFrame"><V text="conflicting source value"/></Q></Feature></Features>"##).unwrap();
        let mut doc = read("frame.dna", &bytes).unwrap();
        assert_eq!(doc.features[0].qualifier("codon_start"), None);
        let values: Vec<_> = doc.features[0]
            .qualifiers
            .iter()
            .filter(|(key, _)| key == "snapgene_readingFrame")
            .map(|(_, value)| value.as_deref())
            .collect();
        assert_eq!(values, vec![Some("conflicting source value"), Some("1")]);
        assert_eq!(write(&doc).unwrap(), bytes);
        doc.description = "force supported-record reconstruction".into();
        assert!(
            write_warnings(&doc)
                .iter()
                .any(|warning| warning.contains("readingFrame"))
        );
        let rebuilt = read("rebuilt.dna", &write(&doc).unwrap()).unwrap();
        let values: Vec<_> = rebuilt.features[0]
            .qualifiers
            .iter()
            .filter(|(key, _)| key == "snapgene_readingFrame")
            .map(|(_, value)| value.as_deref())
            .collect();
        assert_eq!(values, vec![Some("conflicting source value"), Some("1")]);
    }
    #[test]
    fn segmented_genbank_and_native_fixture_display_and_edits_agree() {
        let gb = crate::import(
            "segments.gb",
            include_bytes!("../../../validation/fixtures/regression/segment-colours.gb"),
        )
        .unwrap();
        let mut native = read(
            "segments.dna",
            include_bytes!("../../../validation/fixtures/regression/segment-colours-native.dna"),
        )
        .unwrap();
        assert_eq!(
            gb.features[0].parts,
            vec![Segment {
                start: 100,
                end: 400
            }]
        );
        assert_eq!(native.features[0].parts.len(), 2);
        assert_eq!(
            gb.features[0].display_segments(),
            native.features[0].display_segments()
        );
        assert_eq!(
            gb.features[0].extract(&gb.sequence),
            native.features[0].extract(&native.sequence)
        );
        assert_eq!(
            write(&native).unwrap(),
            include_bytes!("../../../validation/fixtures/regression/segment-colours-native.dna")
        );
        let saved = read("rebuilt.dna", &write(&gb).unwrap()).unwrap();
        assert_eq!(
            saved.features[0].display_segments(),
            gb.features[0].display_segments()
        );
        native.replace(120, 121, "T").unwrap();
        for bytes in [
            write(&native).unwrap(),
            crate::genbank_bytes(&native).unwrap(),
        ] {
            let reread = crate::import("roundtrip", &bytes).unwrap();
            assert_eq!(
                native.features[0].display_segments(),
                reread.features[0].display_segments()
            );
            assert!(reread.features[0].qualifiers.iter().any(|(_, value)| {
                value
                    .as_deref()
                    .is_some_and(|value| value.contains("Original explanatory note"))
            }));
        }
        native.revcomp().unwrap();
        let saved = read("reverse.dna", &write(&native).unwrap()).unwrap();
        assert_eq!(
            saved.features[0].display_segments(),
            native.features[0].display_segments()
        );
        assert_eq!(
            saved.features[0].display_segments()[0].name.as_deref(),
            Some("Blue domain")
        );
        assert!(saved.features[0].reverse);
    }
    #[test]
    fn native_parent_matches_independent_genbank() {
        let d = read("parent.dna", DNA).unwrap();
        let g = super::super::import(
            "fixture.gb",
            include_bytes!("../../../validation/fixtures/annotated.gb"),
        )
        .unwrap();
        assert_eq!(d.sequence, g.sequence);
        assert!(d.circular);
        assert_eq!(d.features.len(), 6);
        for (a, b) in d.features.iter().zip(g.features.iter()) {
            assert_eq!(a.label, b.label);
            assert_eq!(a.parts, b.parts);
            assert_eq!(a.reverse, b.reverse);
        }
        assert_eq!(d.accession.as_deref(), Some("POPGENE823"));
        assert_eq!(d.source.unwrap().bytes.as_ref(), DNA);
    }
    #[test]
    fn pcr_primers_use_zero_based_inclusive_source_ranges() {
        let d = read(
            "pcr.dna",
            include_bytes!("../../../validation/fixtures/synthetic-pcr.dna"),
        )
        .unwrap();
        assert_eq!(d.sequence.len(), 1024);
        assert!(!d.circular);
        assert_eq!(d.features.len(), 4);
        assert_eq!(d.features[0].parts, vec![Segment { start: 0, end: 24 }]);
        assert_eq!(
            d.features[1].parts,
            vec![Segment {
                start: 1000,
                end: 1024
            }]
        );
        assert!(d.features[1].reverse);
    }
    #[test]
    fn rejects_truncation_and_dtd() {
        for end in [0, 4, 10, DNA.len() - 1] {
            assert!(read("bad.dna", &DNA[..end]).is_err());
        }
        assert!(xml(b"<!DOCTYPE x [<!ENTITY a 'b'>]><x>&a;</x>").is_err());
    }
    #[test]
    fn unchanged_native_save_retains_all_source_bytes() {
        let mut d = read("parent.dna", DNA).unwrap();
        d.warnings.clear();
        d.features[0].id = 999;
        assert_eq!(write(&d).unwrap(), DNA);
        assert!(write_warnings(&d).is_empty());
    }

    #[test]
    fn feature_and_primer_review_status_survives_as_conventional_notes() {
        let mut d = read("parent.dna", DNA).unwrap();
        d.features[0].edited = true;
        let primer_index = d
            .features
            .iter()
            .position(|f| f.kind == "primer_bind")
            .unwrap();
        d.features[primer_index].edited = true;
        let labels = [
            d.features[0].label.clone(),
            d.features[primer_index].label.clone(),
        ];
        assert!(unchanged_source(&d).is_none());
        assert!(native_primer(&d.features[primer_index], &d).is_none());
        assert!(write_warnings(&d).iter().any(|w| w.contains("review note")));
        let bytes = write(&d).unwrap();
        let readback = read("review.dna", &bytes).unwrap();
        for label in &labels {
            let feature = readback
                .features
                .iter()
                .find(|f| &f.label == label)
                .unwrap();
            assert!(!feature.edited);
            assert!(
                feature
                    .qualifiers
                    .contains(&("note".into(), Some(REVIEW_NOTE.into())))
            );
        }
        // An ordinary external qualifier carries the review information; no
        // hidden extension or edited flag is needed, and later saves retain it.
        assert_eq!(write(&readback).unwrap(), bytes);
        if let Ok(directory) = std::env::var("POPGENE_INTEROP_OUT") {
            std::fs::create_dir_all(&directory).unwrap();
            std::fs::write(Path::new(&directory).join("popgene-review.dna"), bytes).unwrap();
        }
    }

    #[test]
    fn edited_native_save_rebuilds_supported_blocks_without_stale_history() {
        let mut source = DNA.to_vec();
        packet(&mut source, 250, b"synthetic opaque record").unwrap();
        let mut d = read("parent.dna", &source).unwrap();
        assert_eq!(write(&d).unwrap(), source);
        d.description = "Edited synthetic fixture".into();
        let bytes = write(&d).unwrap();
        assert_ne!(bytes, DNA);
        let readback = read("edited.dna", &bytes).unwrap();
        assert_eq!(readback.sequence, d.sequence);
        assert_eq!(readback.name, d.name);
        assert_eq!(readback.description, d.description);
        assert_eq!(readback.accession, d.accession);
        assert_eq!(readback.features.len(), d.features.len());
        for (before, after) in d.features.iter().zip(&readback.features) {
            assert_eq!(before.label, after.label);
            assert_eq!(before.parts, after.parts);
            assert_eq!(before.reverse, after.reverse);
            assert_eq!(
                before.extract(&d.sequence),
                after.extract(&readback.sequence)
            );
        }
        assert!(
            readback
                .source
                .unwrap()
                .retained_blocks
                .iter()
                .any(|(kind, _)| *kind == 250)
        );
        assert!(write_warnings(&d).iter().any(|w| w.contains("history")));
        // Optional artifact for independent readers and the installed SnapGene.
        if let Ok(directory) = std::env::var("POPGENE_INTEROP_OUT") {
            std::fs::create_dir_all(&directory).unwrap();
            std::fs::write(Path::new(&directory).join("popgene-edited.dna"), bytes).unwrap();
        }
    }

    #[test]
    fn synthetic_origin_reverse_ambiguity_and_metadata_roundtrip() {
        let mut d = Document::new(
            "Synthetic <map> & α",
            "ACGTRYSWKMBDHVNACGTACGTACGTACGTACGTA",
        )
        .unwrap();
        d.circular = true;
        d.description = "Literal <angle> & \"quotes\"\nSecond line".into();
        d.comments = vec!["Plain text with  two spaces".into()];
        d.accession = Some("SYNTHETIC_1".into());
        d.features = vec![
            Feature {
                location: None,
                id: 10,
                label: "reverse compound".into(),
                kind: "misc_feature".into(),
                parts: vec![
                    Segment { start: 30, end: 36 },
                    Segment { start: 0, end: 4 },
                    Segment { start: 8, end: 12 },
                ],
                reverse: true,
                edited: false,
                qualifiers: vec![
                    ("note".into(), Some("a < b & b > c".into())),
                    ("note".into(), Some("second value".into())),
                    ("pseudo".into(), None),
                ],
            },
            Feature {
                location: None,
                id: 11,
                label: "origin primer".into(),
                kind: "primer_bind".into(),
                parts: vec![Segment { start: 30, end: 36 }, Segment { start: 0, end: 4 }],
                reverse: false,
                edited: false,
                qualifiers: vec![],
            },
        ];
        let primer = d.features[1].extract(&d.sequence);
        d.features[1]
            .qualifiers
            .push(("primer_sequence".into(), Some(primer.clone())));
        let bytes = write(&d).unwrap();
        let readback = read("different-name.dna", &bytes).unwrap();
        assert_eq!(readback.sequence, d.sequence);
        assert!(readback.circular);
        assert_eq!(readback.name, d.name);
        assert_eq!(readback.description, d.description);
        assert_eq!(readback.comments, d.comments);
        assert_eq!(readback.features[0].parts, d.features[0].parts);
        assert!(readback.features[0].reverse);
        assert_eq!(
            readback.features[0].extract(&readback.sequence),
            d.features[0].extract(&d.sequence)
        );
        for q in &d.features[0].qualifiers {
            assert!(readback.features[0].qualifiers.contains(q));
        }
        assert_eq!(
            qualifier(&readback.features[1], "primer_sequence"),
            Some(primer.as_str())
        );
        assert_eq!(readback.features[1].parts, d.features[1].parts);
        if let Ok(directory) = std::env::var("POPGENE_INTEROP_OUT") {
            std::fs::create_dir_all(&directory).unwrap();
            std::fs::write(Path::new(&directory).join("popgene-synthetic.dna"), bytes).unwrap();
        }
    }

    #[test]
    fn strandedness_and_topology_are_independent_flags() {
        for single_stranded in [false, true] {
            for circular in [false, true] {
                let mut d = Document::new("flags", "ACGTRYSWKMBDHVN").unwrap();
                d.single_stranded = single_stranded;
                d.circular = circular;
                let readback = read("flags.dna", &write(&d).unwrap()).unwrap();
                assert_eq!(readback.single_stranded, single_stranded);
                assert_eq!(readback.circular, circular);
            }
        }
    }

    #[test]
    fn mismatched_primers_keep_annotations_without_claiming_native_binding() {
        let mut d = read("parent.dna", DNA).unwrap();
        let primer = d
            .features
            .iter_mut()
            .find(|f| f.kind == "primer_bind")
            .unwrap();
        let label = primer.label.clone();
        primer
            .qualifiers
            .iter_mut()
            .find(|(k, _)| k == "primer_sequence")
            .unwrap()
            .1 = Some("AAAA".into());
        let readback = read("edited.dna", &write(&d).unwrap()).unwrap();
        let primer = readback.features.iter().find(|f| f.label == label).unwrap();
        assert_eq!(qualifier(primer, "primer_sequence"), Some("AAAA"));
        assert!(
            write_warnings(&d)
                .iter()
                .any(|w| w.contains("remain annotations"))
        );
    }

    #[test]
    fn rejects_invalid_document_metadata() {
        let mut d = Document::new("bad metadata", "ACGT").unwrap();
        d.description = "invalid\0XML".into();
        assert!(write(&d).is_err());
        d.description.clear();
        d.sequence = "!".into();
        assert!(write(&d).is_err());
        let mut d = read("parent.dna", DNA).unwrap();
        d.source.as_mut().unwrap().bytes = b"invalid source".as_slice().into();
        assert_ne!(write(&d).unwrap(), b"invalid source");
    }
}
