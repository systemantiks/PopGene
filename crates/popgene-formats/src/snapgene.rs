//! Native .dna interchange using the public TLV/XML descriptions and Biopython.
//! Unchanged files retain their bytes. Edited files rebuild supported records; opaque
//! history, alignments and other source blocks must never describe the edited molecule.
use super::{Error, Result};
use popgene_core::{Document, Feature, FeatureDirection, Molecule, Segment, Source};
use quick_xml::{Reader, events::Event};
use std::{collections::BTreeMap, path::Path};

const REVIEW_NOTE: &str = "Sequence changed in PopGene; annotation needs review.";

fn bad(message: impl Into<String>) -> Error {
    Error::Message(format!("SnapGene import: {}", message.into()))
}
#[derive(Default)]
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

/// Notes/reference fields from the retained original container. This path does
/// not populate the editable molecule or claim to reconstruct operation history.
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
        if kind == 6 {
            let root = xml(&bytes[offset..end])?;
            if root.name != "Notes" {
                return Err(bad("Invalid notes root"));
            }
            for node in &root.children {
                let mut value = String::new();
                describe(node, &mut value, 0);
                total += value.len();
                if total > 1024 * 1024 || fields.len() >= 10_000 {
                    fields.push(("Display limit".into(), "Additional notes remain in the retained original file (metadata display limited to 1 MiB / 10,000 fields).".into()));
                    return Ok(fields);
                }
                fields.push((node.name.clone(), value.trim_end().to_owned()));
            }
        }
        offset = end;
    }
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
    if u16::from_be_bytes(cookie[8..10].try_into().unwrap()) != 1 {
        return Err(bad("Only DNA documents are currently supported"));
    }
    let dna: Vec<_> = blocks.iter().filter(|(k, _)| *k == 0).collect();
    if dna.len() != 1 || dna[0].1.len() < 2 {
        return Err(bad("Expected exactly one DNA block"));
    }
    let seq = std::str::from_utf8(&dna[0].1[1..]).map_err(|_| bad("Invalid DNA encoding"))?;
    let title = Path::new(name)
        .file_stem()
        .unwrap_or_default()
        .to_string_lossy();
    let mut doc = Document::new(title, seq)?;
    doc.circular = dna[0].1[0] & 1 != 0;
    doc.single_stranded = dna[0].1[0] & 2 == 0;
    let n = doc.sequence.len();
    let mut unbound = 0;
    for (kind, body) in &blocks {
        match kind {
            10 => {
                let root = xml(body)?;
                if root.name != "Features" {
                    return Err(bad("Invalid features root"));
                }
                for f in root.children.iter().filter(|f| f.name == "Feature") {
                    let mut parts = vec![];
                    let mut qualifiers = vec![];
                    let mut display_segments = vec![];
                    let native_segments: Vec<_> =
                        f.children.iter().filter(|s| s.name == "Segment").collect();
                    for seg in f.children.iter().filter(|s| s.name == "Segment") {
                        if seg.attr("type") == "gap" {
                            return Err(bad(
                                "Gapped feature segments need a richer location model",
                            ));
                        }
                        let ranges = ranges(seg.attr("range"), n, false, doc.circular)?;
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
                    if parts.is_empty() {
                        return Err(bad("Feature has no location"));
                    }
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
                    doc.features.push(Feature {
                        id: doc.features.len() as u64 + 1,
                        label,
                        kind: if f.attr("type").is_empty() {
                            "misc_feature".into()
                        } else {
                            f.attr("type").into()
                        },
                        parts,
                        reverse: f.attr("directionality") == "2",
                        qualifiers,
                        edited: false,
                    });
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
                        unbound += 1;
                    }
                }
            }
            6 => {
                let root = xml(body)?;
                if root.name != "Notes" {
                    return Err(bad("Invalid notes root"));
                }
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
                for key in ["Created", "LastModified", "Type"] {
                    if let Some(v) = root.child(key) {
                        doc.comments.push(format!(
                            "Imported {key}: {}{}",
                            v.text,
                            v.attrs
                                .get("UTC")
                                .map_or(String::new(), |t| format!(" {t} UTC"))
                        ));
                    }
                }
            }
            _ => {}
        }
        if doc.features.len() > 20_000 {
            return Err(bad("More than 20,000 annotations"));
        }
    }
    doc.warnings.push("SnapGene import: DNA, topology, strandedness, exact feature ranges, strand, qualifiers, primer binding annotations and selected notes are interpreted. An unchanged .dna save retains the original bytes. Saving edits rebuilds supported records; display settings, full oligo models, references, end chemistry, alignments and procedure history are not reconstructed.".into());
    doc.warnings.push("Primer binding sites are imported without recalculating SnapGene hybridization filters; duplicate simplified sites are collapsed. Imported melting temperatures are source values, not PopGene calculations.".into());
    if unbound > 0 {
        doc.warnings.push(format!(
            "{unbound} unbound primer(s) remain only in the retained original file."
        ));
    }
    doc.source = Some(Source {
        file_name: name.into(),
        format: "SnapGene".into(),
        bytes: bytes.into(),
        retained_blocks: blocks.iter().map(|(k, b)| (*k, b.len())).collect(),
    });
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
        for (index, feature) in d.features.iter_mut().enumerate() {
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

/// An oligo is written as a native Primer only when its exact full-length binding
/// can be established independently. Mismatches, tails and discontinuous sites
/// remain ordinary primer_bind annotations with their original qualifiers.
fn native_primer<'a>(feature: &'a Feature, doc: &Document) -> Option<&'a str> {
    if feature.kind != "primer_bind" || feature.edited {
        return None;
    }
    let sequence = qualifier(feature, "primer_sequence")?;
    if sequence.is_empty()
        || !sequence.bytes().all(|c| b"ACGTRYSWKMBDHVN".contains(&c))
        || feature.parts.is_empty()
        || feature
            .parts
            .iter()
            .any(|p| p.start >= p.end || p.end > doc.sequence.len())
        || sequence != feature.extract(&doc.sequence)
    {
        return None;
    }
    match feature.parts.as_slice() {
        [_] => {}
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
    let mut warnings = vec![];
    if doc.features.iter().any(|f| f.edited) {
        warnings.push("Annotations changed by sequence edits carry a review note; the in-memory review flag is not stored. Primers requiring review are saved as annotated primer_bind features.".into());
    }
    if doc.features.iter().any(|feature| {
        feature.kind == "CDS" && popgene_core::translation::cds_phase(feature).is_err()
    }) {
        warnings.push("Invalid or conflicting CDS readingFrame/codon_start values remain qualifiers; the native readingFrame attribute is omitted. Review these CDS offsets before using their translations.".into());
    }
    if doc.source.as_ref().is_some_and(|s| s.format == "SnapGene") {
        warnings.push("Edited SnapGene files are rebuilt from the supported sequence, annotations, primers and notes. Sequence letter case is normalized to uppercase. Original procedure history, alignments, references, unbound primers, end chemistry, methylation and display settings are omitted. Keep the original file if these are needed.".into());
    }
    if doc.features.iter().any(|f| native_primer(f, doc).is_some()) {
        warnings.push("Exact full-length primer matches are saved as native primers. Source melting temperatures, hybridization settings and detailed oligo chemistry are not copied; SnapGene may recalculate binding sites.".into());
    }
    if doc
        .features
        .iter()
        .any(|f| f.kind == "primer_bind" && native_primer(f, doc).is_none())
    {
        warnings.push("Primer sites without a verified full-length exact oligo match, or with additional qualifiers, are saved as primer_bind features. Their sequences and qualifiers remain annotations rather than native SnapGene primers.".into());
    }
    if doc.trace.is_some() {
        warnings.push("The .dna file includes called bases but not chromatogram channels, peaks or quality scores. Keep the original trace file.".into());
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

fn packet(out: &mut Vec<u8>, kind: u8, body: &[u8]) -> Result<()> {
    let size = u32::try_from(body.len())
        .map_err(|_| Error::Message("SnapGene export: packet exceeds 4 GiB".into()))?;
    if matches!(kind, 5 | 6 | 10) {
        // Apply the same size and structure limits to generated and imported XML.
        xml(body)?;
    }
    out.push(kind);
    out.extend_from_slice(&size.to_be_bytes());
    out.extend_from_slice(body);
    Ok(())
}

/// Write conventional SnapGene TLV records, never application-specific extensions.
/// The format flags are documented by SnapGeneReader; feature/primer coordinates
/// and reverse-strand compound order are independently covered by Biopython.
pub(crate) fn write(doc: &Document) -> Result<Vec<u8>> {
    doc.validate()?;
    if doc.molecule != Molecule::Dna || doc.sequence.is_empty() {
        return Err(Error::Message("SnapGene .dna export requires a non-empty DNA sequence; use GenBank or FASTA for RNA and protein.".into()));
    }
    if let Some(original) = unchanged_source(doc) {
        return Ok(original.to_vec());
    }
    if doc.features.len() > 20_000 {
        return Err(Error::Message(
            "SnapGene export: more than 20,000 annotations".into(),
        ));
    }
    let texts =
        std::iter::once(doc.name.as_str())
            .chain(std::iter::once(doc.description.as_str()))
            .chain(doc.accession.as_deref())
            .chain(doc.comments.iter().map(String::as_str))
            .chain(doc.features.iter().flat_map(|f| {
                std::iter::once(f.label.as_str())
                    .chain(std::iter::once(f.kind.as_str()))
                    .chain(f.qualifiers.iter().flat_map(|(key, value)| {
                        std::iter::once(key.as_str()).chain(value.as_deref())
                    }))
            }));
    if texts.into_iter().any(|text| text.chars().any(|c| {
        !matches!(c, '\t' | '\n' | '\r' | '\u{20}'..='\u{d7ff}' | '\u{e000}'..='\u{fffd}' | '\u{10000}'..='\u{10ffff}')
    })) {
        return Err(Error::Message("SnapGene export: metadata contains characters that XML cannot represent".into()));
    }
    if doc
        .features
        .iter()
        .any(|f| f.qualifiers.iter().any(|(key, _)| key.is_empty()))
    {
        return Err(Error::Message(
            "SnapGene export: a feature qualifier has an empty name".into(),
        ));
    }
    let mut out = Vec::with_capacity(doc.sequence.len() + 1024);
    // Export version 15 / import version 20 are used by the installed 8.2.3
    // synthetic fixtures. This writer emits only uncompressed, established blocks.
    packet(&mut out, 9, b"SnapGene\x00\x01\x00\x0f\x00\x14")?;
    let mut dna = Vec::with_capacity(doc.sequence.len() + 1);
    dna.push(u8::from(doc.circular) | if doc.single_stranded { 0 } else { 2 });
    dna.extend_from_slice(doc.sequence.as_bytes());
    packet(&mut out, 0, &dna)?;

    let mut features = format!(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?><Features nextValidID=\"{}\">",
        doc.features.len()
    );
    let mut primers = format!(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?><Primers nextValidID=\"{}\">",
        doc.features.len()
    );
    let mut primer_id = 0;
    for (index, feature) in doc.features.iter().enumerate() {
        if let Some(sequence) = native_primer(feature, doc) {
            let start = feature.parts[0].start;
            let end = feature.parts.last().unwrap().end - 1;
            let range = format!("{start}-{end}");
            primers.push_str(&format!(
                "<Primer recentID=\"{primer_id}\" name=\"{}\" sequence=\"{}\"",
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
            primers.push_str(&format!("><BindingSite location=\"{range}\" boundStrand=\"{}\" annealedBases=\"{sequence}\"><Component hybridizedRange=\"{range}\" bases=\"{sequence}\"/></BindingSite></Primer>", u8::from(feature.reverse)));
            primer_id += 1;
            continue;
        }
        let direction = if feature.reverse {
            "2"
        } else {
            match feature.display_direction() {
                Some(FeatureDirection::Nondirectional) => "0",
                Some(FeatureDirection::Bidirectional) => "3",
                _ => "1",
            }
        };
        let reading_frame = (feature.kind == "CDS")
            .then(|| popgene_core::translation::cds_phase(feature).ok())
            .flatten()
            .map(|phase| (phase as i8 + 1) * if feature.reverse { -1 } else { 1 });
        features.push_str(&format!("<Feature recentID=\"{index}\" name=\"{}\" type=\"{}\" directionality=\"{direction}\" allowSegmentOverlaps=\"1\"", escaped(&feature.label), escaped(&feature.kind)));
        if let Some(frame) = reading_frame {
            features.push_str(&format!(" readingFrame=\"{frame}\""));
        }
        features.push('>');
        let display_segments = feature.display_segments();
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
                features.push_str(" color=\"#5b9bd5\"");
            }
            if let Some(name) = &display_segment.name {
                features.push_str(&format!(" name=\"{}\"", escaped(name)));
            }
            if feature.kind == "CDS" {
                features.push_str(" translated=\"1\"");
            }
            features.push_str("/>");
        }
        // Group equal keys so readers that use a dictionary do not overwrite
        // earlier values. Display metadata remains in standard Q records too.
        let mut qualifiers: BTreeMap<&str, Vec<Option<&str>>> = BTreeMap::new();
        for (key, value) in &feature.qualifiers {
            if key != "snapgene_color"
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
    features.push_str("</Features>");
    primers.push_str("</Primers>");
    packet(&mut out, 10, features.as_bytes())?;
    if primer_id > 0 {
        packet(&mut out, 5, primers.as_bytes())?;
    }
    let mut notes = format!(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?><Notes><CustomMapLabel>{}</CustomMapLabel><Description>{}</Description>",
        escaped(&doc.name),
        escaped(&doc.description)
    );
    if let Some(accession) = &doc.accession {
        notes.push_str(&format!(
            "<AccessionNumber>{}</AccessionNumber>",
            escaped(accession)
        ));
    }
    if !doc.comments.is_empty() {
        notes.push_str(&format!(
            "<Comments>{}</Comments>",
            escaped(&doc.comments.join("\n"))
        ));
    }
    notes.push_str("</Notes>");
    packet(&mut out, 6, notes.as_bytes())?;
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    const DNA: &[u8] = include_bytes!("../../../validation/fixtures/annotated.dna");
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
                .all(|(kind, _)| matches!(kind, 0 | 5 | 6 | 9 | 10))
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
    fn rejects_invalid_document_metadata_and_unsupported_molecules() {
        for molecule in [Molecule::Rna, Molecule::Protein] {
            let d = Document::new_typed("not DNA", "ACG", molecule).unwrap();
            assert!(write(&d).is_err());
        }
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
