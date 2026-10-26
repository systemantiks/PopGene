//! Readable current-fragment supplement for GenBank/EMBL comments.
//! Ordinary sequence records remain the actual top strand. This convention is
//! independent of annotation ancestry and contains no ancestor sequence.
use crate::{Error, Result, SaveFormat};
use popgene_core::{
    Document, Molecule,
    metadata::{DnaMethylation, DuplexAnnotations, Reference, SequenceColor},
    workflow::DuplexEnds,
};
use serde_json::Value;
use sha2::{Digest, Sha256};

const START: &str = "[PopGene physical fragment v";
const SUFFIX: &str = " begin]";
const END: &str = "[PopGene physical fragment end]";
const VERSION: &str = "1";
const MAX_PAYLOAD: usize = 4 * 1024 * 1024;
const MAX_DEPTH: usize = 96;
const MAX_FEATURES: usize = 8192;

fn bad(message: &str) -> Error {
    Error::Message(format!("Physical fragment supplement: {message}"))
}
fn hash(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}
fn format_name(format: SaveFormat) -> &'static str {
    match format {
        SaveFormat::GenBank => "genbank",
        SaveFormat::Embl => "embl",
        _ => "unsupported",
    }
}

struct Fragment {
    ends: DuplexEnds,
    annotations: DuplexAnnotations,
    colors: Vec<SequenceColor>,
    top_qualities: Option<Vec<Option<u8>>>,
    opposite_qualities: Option<Vec<Option<u8>>>,
    methylation: Option<DnaMethylation>,
    references: Vec<Reference>,
}
/// Count borrowed fields before cloning or expanding them into a JSON tree.
/// Optional transport must not allocate a giant tree just to reject its size.
struct PayloadBudget(usize);
impl std::io::Write for PayloadBudget {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        if bytes.len() > MAX_PAYLOAD.saturating_sub(self.0) {
            return Err(std::io::Error::other(
                "physical fragment payload exceeds 4 MiB",
            ));
        }
        self.0 += bytes.len();
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn preflight(doc: &Document) -> Result<()> {
    if doc.duplex_annotations.features.len() > MAX_FEATURES {
        return Err(bad("full-extent annotation count exceeds 8192"));
    }
    let mut budget = PayloadBudget(0);
    serde_json::to_writer(&mut budget, &doc.duplex_ends)?;
    serde_json::to_writer(&mut budget, &doc.duplex_annotations)?;
    serde_json::to_writer(&mut budget, &doc.sequence_colors)?;
    serde_json::to_writer(&mut budget, &doc.sequence_qualities)?;
    serde_json::to_writer(&mut budget, &doc.opposite_sequence_qualities)?;
    serde_json::to_writer(&mut budget, &doc.methylation)?;
    serde_json::to_writer(&mut budget, &doc.references)?;
    Ok(())
}
impl Fragment {
    fn from_document(doc: &Document) -> Result<Self> {
        preflight(doc)?;
        let mut annotations = doc.duplex_annotations.clone();
        if annotations.features.len() > MAX_FEATURES {
            return Err(bad("full-extent annotation count exceeds 8192"));
        }
        for (index, feature) in annotations.features.iter_mut().enumerate() {
            // The ordinary parser assigns 1..N in actual emission order. The
            // supplemental occurrences follow those records, never matching by
            // label, sequence or otherwise identical annotation content.
            feature.id = (doc.features.len() + index + 1) as u64;
        }
        Ok(Self {
            ends: doc.duplex_ends.clone().ok_or_else(|| bad("missing ends"))?,
            annotations,
            colors: doc.sequence_colors.clone(),
            top_qualities: doc.sequence_qualities.clone(),
            opposite_qualities: doc.opposite_sequence_qualities.clone(),
            methylation: doc.methylation,
            references: doc.references.clone(),
        })
    }
    fn value(&self) -> Value {
        serde_json::json!({
            "ends": self.ends,
            "duplex_annotations": self.annotations,
            "sequence_colors": self.colors,
            "sequence_qualities": self.top_qualities,
            "opposite_sequence_qualities": self.opposite_qualities,
            "methylation": self.methylation,
            "references": self.references,
        })
    }
    fn parse(value: &Value) -> Result<Self> {
        let field = |key| {
            value
                .get(key)
                .cloned()
                .ok_or_else(|| bad("missing fragment field"))
        };
        let result = Self {
            ends: serde_json::from_value(field("ends")?)?,
            annotations: serde_json::from_value(field("duplex_annotations")?)?,
            colors: serde_json::from_value(field("sequence_colors")?)?,
            top_qualities: serde_json::from_value(field("sequence_qualities")?)?,
            opposite_qualities: serde_json::from_value(field("opposite_sequence_qualities")?)?,
            methylation: serde_json::from_value(field("methylation")?)?,
            references: serde_json::from_value(field("references")?)?,
        };
        if result.value() != *value {
            return Err(bad("unknown or noncanonical fragment fields"));
        }
        Ok(result)
    }
    fn apply(self, doc: &mut Document) -> Result<()> {
        if doc.molecule != Molecule::Dna
            || doc.single_stranded
            || doc.circular
            || doc.trace.is_some()
            || doc.duplex_ends.is_some()
            || !doc.duplex_annotations.is_empty()
        {
            return Err(bad(
                "requires an ordinary linear double-stranded DNA record",
            ));
        }
        if self.annotations.features.len() > MAX_FEATURES
            || self
                .annotations
                .features
                .iter()
                .enumerate()
                .any(|(index, f)| f.id != (doc.features.len() + index + 1) as u64)
        {
            return Err(bad(
                "invalid supplemental annotation occurrence IDs or count",
            ));
        }
        self.ends.validate(&doc.sequence)?;
        if self.references.len() != doc.references.len() {
            return Err(bad(
                "supplemental references do not preserve ordinary numeric slots",
            ));
        }
        doc.references = self.references;
        for feature in doc.features.iter().chain(&self.annotations.features) {
            for (_, value) in feature
                .qualifiers
                .iter()
                .filter(|(key, _)| key == "citation")
            {
                citation(doc, value.as_deref())?;
            }
        }
        doc.duplex_ends = Some(self.ends);
        doc.duplex_annotations = self.annotations;
        doc.sequence_colors = self.colors;
        doc.sequence_qualities = self.top_qualities;
        doc.opposite_sequence_qualities = self.opposite_qualities;
        doc.methylation = self.methylation;
        doc.validate()?;
        Ok(())
    }
}

fn citation<'a>(
    doc: &'a Document,
    value: Option<&str>,
) -> Result<&'a popgene_core::metadata::Reference> {
    value
        .and_then(|v| {
            v.trim()
                .trim_start_matches('[')
                .trim_end_matches(']')
                .parse::<usize>()
                .ok()
        })
        .and_then(|i| i.checked_sub(1))
        .and_then(|i| doc.references.get(i))
        .ok_or_else(|| bad("a full-extent annotation citation has no current reference"))
}

fn references_for_codec(
    source: &Document,
    canonical: &Document,
    format: SaveFormat,
) -> Result<Vec<Reference>> {
    // Validate source citations before a codec may append a nonbiological
    // software marker. A formerly invalid number must not become its citation.
    for feature in source
        .features
        .iter()
        .chain(&source.duplex_annotations.features)
    {
        for (_, value) in feature
            .qualifiers
            .iter()
            .filter(|(key, _)| key == "citation")
        {
            citation(source, value.as_deref())?;
        }
    }
    let mut references = source.references.clone();
    if canonical.references.len() == references.len() {
        return Ok(references);
    }
    if format == SaveFormat::GenBank
        && canonical.references.len() == references.len() + 1
        && crate::genbank_display::is_generated_reference(
            canonical.references.last().unwrap(),
            references.len(),
        )
    {
        references.push(canonical.references.last().unwrap().clone());
        Ok(references)
    } else {
        Err(bad(
            "ordinary codec changes reference slots beyond its verified display marker",
        ))
    }
}

fn context(doc: &Document) -> Result<crate::canonical_json::ContextHash> {
    // Field order also reconstructs the one pre-canonical preserve_order hash.
    #[derive(serde::Serialize)]
    struct Context<'a> {
        name: &'a str,
        description: &'a str,
        accession: &'a Option<String>,
        sequence_sha256: String,
        molecule: popgene_core::Molecule,
        circular: bool,
        single_stranded: bool,
        features: &'a [popgene_core::Feature],
        references: &'a [Reference],
    }
    crate::canonical_json::context(&Context {
        name: &doc.name,
        description: &doc.description,
        accession: &doc.accession,
        sequence_sha256: hash(doc.sequence.as_bytes()),
        molecule: doc.molecule,
        circular: doc.circular,
        single_stranded: doc.single_stranded,
        features: &doc.features,
        references: &doc.references,
    })
}

fn depth_check(text: &str) -> Result<()> {
    let (mut depth, mut quoted, mut escaped) = (0_usize, false, false);
    for byte in text.bytes() {
        if quoted {
            if escaped {
                escaped = false;
            } else if byte == b'\\' {
                escaped = true;
            } else if byte == b'"' {
                quoted = false;
            }
        } else {
            match byte {
                b'"' => quoted = true,
                b'{' | b'[' => {
                    depth += 1;
                    if depth > MAX_DEPTH {
                        return Err(bad("JSON nesting exceeds 96 levels"));
                    }
                }
                b'}' | b']' => depth = depth.saturating_sub(1),
                _ => {}
            }
        }
    }
    Ok(())
}
fn block(value: &Value) -> Result<String> {
    let text = serde_json::to_string(value)?;
    let mut compact = String::with_capacity(text.len().min(MAX_PAYLOAD));
    for ch in text.chars() {
        if ch.is_whitespace() {
            for unit in ch.encode_utf16(&mut [0; 2]) {
                use std::fmt::Write as _;
                write!(compact, "\\u{unit:04x}").unwrap();
            }
        } else {
            compact.push(ch);
        }
        if compact.len() > MAX_PAYLOAD {
            return Err(bad("readable payload exceeds 4 MiB"));
        }
    }
    depth_check(&compact)?;
    let header = format!(
        "bytes={};sha256={};json={compact}",
        compact.len(),
        hash(compact.as_bytes())
    );
    let mut result = format!("{START}{VERSION}{SUFFIX}\n");
    let mut width = 0;
    for ch in header.chars() {
        if width + ch.len_utf8() > 48 {
            result.push('\n');
            width = 0;
        }
        result.push(ch);
        width += ch.len_utf8();
    }
    result.push('\n');
    result.push_str(END);
    Ok(result)
}
struct Parsed {
    comment: usize,
    start: usize,
    end: usize,
    value: Value,
}
fn parse(comments: &[String]) -> Result<Option<Parsed>> {
    let mut found = None;
    let mut ends = 0;
    for (index, comment) in comments.iter().enumerate() {
        for (start, _) in comment.match_indices(START) {
            if found.replace((index, start)).is_some() {
                return Err(bad("duplicate or colliding markers"));
            }
        }
        ends += comment.matches(END).count();
    }
    let Some((comment, start)) = found else {
        return if ends == 0 {
            Ok(None)
        } else {
            Err(bad("unmatched end marker"))
        };
    };
    if ends != 1 {
        return Err(bad("missing or duplicate end marker"));
    }
    let text = &comments[comment];
    let version = start + START.len();
    let suffix = text[version..]
        .find(SUFFIX)
        .filter(|n| *n <= 10)
        .ok_or_else(|| bad("malformed version marker"))?;
    if &text[version..version + suffix] != VERSION {
        return Err(bad("unknown version"));
    }
    let payload_start = version + suffix + SUFFIX.len();
    let payload_end = text[payload_start..]
        .find(END)
        .map(|i| i + payload_start)
        .ok_or_else(|| bad("markers span separate comments"))?;
    let mut compact = String::new();
    for ch in text[payload_start..payload_end].chars() {
        if ch.is_ascii_whitespace() {
            continue;
        }
        if ch.is_whitespace() {
            return Err(bad("non-ASCII wrapping whitespace"));
        }
        compact.push(ch);
        if compact.len() > MAX_PAYLOAD + 128 {
            return Err(bad("payload exceeds 4 MiB"));
        }
    }
    let (length, rest) = compact
        .strip_prefix("bytes=")
        .and_then(|s| s.split_once(";sha256="))
        .ok_or_else(|| bad("invalid length header"))?;
    let length = length
        .parse::<usize>()
        .ok()
        .filter(|n| *n <= MAX_PAYLOAD)
        .ok_or_else(|| bad("invalid payload length"))?;
    let (checksum, payload) = rest
        .split_once(";json=")
        .ok_or_else(|| bad("missing JSON header"))?;
    if length != payload.len() || checksum != hash(payload.as_bytes()) {
        return Err(bad("length or checksum mismatch"));
    }
    depth_check(payload)?;
    Ok(Some(Parsed {
        comment,
        start,
        end: payload_end + END.len(),
        value: serde_json::from_str(payload)?,
    }))
}

/// Apply only after ordinary records are parsed, and before ancestry validation.
/// All checks occur on a clone: failed optional metadata cannot partly change a
/// sequence, annotation set, physical geometry or confidence vector.
pub(super) fn read(doc: &mut Document, bytes: &[u8], format: SaveFormat) {
    let restored = (|| -> Result<Option<(Parsed, Document)>> {
        let Some(parsed) = parse(&doc.comments)? else {
            return Ok(None);
        };
        let obj = parsed
            .value
            .as_object()
            .ok_or_else(|| bad("payload is not an object"))?;
        if obj.len() != 4
            || obj.get("format") != Some(&Value::from(format_name(format)))
            || obj.get("records_sha256")
                != Some(&Value::from(crate::lineage::records_hash(bytes, format)?))
            || !context(doc)?.matches(obj.get("context_sha256"))
        {
            return Err(bad(
                "current sequence, metadata or ordered annotation records changed",
            ));
        }
        let fragment =
            Fragment::parse(obj.get("fragment").ok_or_else(|| bad("missing fragment"))?)?;
        let mut restored = doc.clone();
        fragment.apply(&mut restored)?;
        Ok(Some((parsed, restored)))
    })();
    match restored {
        Ok(Some((parsed, mut restored))) => {
            let comment = &mut restored.comments[parsed.comment];
            comment.replace_range(parsed.start..parsed.end, "");
            if comment.ends_with('\n') {
                comment.pop();
            }
            if comment.is_empty() {
                restored.comments.remove(parsed.comment);
            }
            restored.warnings.push("Physical strands, end chemistry and full-extent annotations were restored from a readable PopGene comment convention. Other readers may show only the ordinary top-strand record; checksums establish integrity, not authenticity.".into());
            *doc = restored;
        }
        Ok(None) => {}
        Err(error) => doc.warnings.push(format!(
            "{error}. Original supplement text is retained; no physical metadata was activated."
        )),
    }
}

struct Prepared {
    bytes: Vec<u8>,
    omitted: Option<String>,
}
fn prepare(
    doc: &Document,
    format: SaveFormat,
    plain: fn(&Document) -> Result<Vec<u8>>,
) -> Result<Prepared> {
    if doc.duplex_ends.is_some() && crate::lineage::marker_candidates(&doc.comments, START, END) > 0
    {
        // Returning ordinary bytes containing a stale *valid* carrier could
        // resurrect old end chemistry on reopen. Active-state collisions need
        // explicit resolution rather than optional-fallback serialization.
        return Err(bad(
            "existing supplement text cannot be replaced safely; preserve it separately or remove it explicitly before saving current physical ends",
        ));
    }
    let first = plain(doc)?;
    if doc.duplex_ends.is_none() {
        // This fixed prefix also catches split START/END wording that ordinary
        // comment wrapping may normalize back into a recognized marker.
        if crate::lineage::marker_candidates(&doc.comments, START, END) > 0
            && crate::import(&format!("physical.{}", format.extension_for(doc)), &first)?
                .duplex_ends
                .is_some()
        {
            return Err(bad(
                "previously inactive supplement text would become active after serialization; preserve it separately or remove it explicitly before saving",
            ));
        }
        return Ok(Prepared {
            bytes: first,
            omitted: None,
        });
    }
    let attempt = (|| -> Result<Vec<u8>> {
        preflight(doc)?;
        let canonical = crate::import(&format!("physical.{}", format.extension_for(doc)), &first)?;
        if doc.sequence != canonical.sequence
            || doc.molecule != canonical.molecule
            || doc.circular != canonical.circular
            || doc.single_stranded != canonical.single_stranded
            || doc.features.len() != canonical.features.len()
        {
            return Err(bad(
                "ordinary codec changes sequence, molecule or annotation occurrence count",
            ));
        }
        let mut fragment = Fragment::from_document(doc)?;
        fragment.references = references_for_codec(doc, &canonical, format)?;
        let payload = serde_json::json!({
            "format": format_name(format),
            "records_sha256": crate::lineage::records_hash(&first, format)?,
            "context_sha256": context(&canonical)?.canonical,
            "fragment": fragment.value(),
        });
        // Validate the complete staged physical state before serializing it.
        let mut expected = canonical.clone();
        Fragment::parse(&payload["fragment"])?.apply(&mut expected)?;
        let mut with_carrier = doc.clone();
        with_carrier.comments.push(block(&payload)?);
        let bytes = plain(&with_carrier)?;
        let reopened = crate::import(&format!("physical.{}", format.extension_for(doc)), &bytes)?;
        if reopened.duplex_ends != expected.duplex_ends
            || reopened.duplex_annotations != expected.duplex_annotations
            || reopened.sequence_colors != expected.sequence_colors
            || reopened.sequence_qualities != expected.sequence_qualities
            || reopened.opposite_sequence_qualities != expected.opposite_sequence_qualities
            || reopened.methylation != expected.methylation
            || reopened.references != expected.references
        {
            return Err(bad(
                "written supplement did not pass independent import validation",
            ));
        }
        Ok(bytes)
    })();
    match attempt {
        Ok(bytes) => Ok(Prepared {
            bytes,
            omitted: None,
        }),
        Err(reason) => Ok(Prepared {
            bytes: first,
            omitted: Some(reason.to_string()),
        }),
    }
}
pub(super) fn write(
    doc: &Document,
    format: SaveFormat,
    plain: fn(&Document) -> Result<Vec<u8>>,
) -> Result<Vec<u8>> {
    Ok(prepare(doc, format, plain)?.bytes)
}
/// Whether this exact conversion carries the supplement, for accurate loss
/// warnings. This uses the same preparation and read-back path as the writer.
pub(super) fn preservation(doc: &Document, format: SaveFormat) -> (bool, Vec<String>) {
    if doc.duplex_ends.is_none() || !matches!(format, SaveFormat::GenBank | SaveFormat::Embl) {
        return (false, vec![]);
    }
    let plain = if format == SaveFormat::GenBank {
        crate::genbank_bytes_plain
    } else {
        crate::embl::write_plain
    };
    match prepare(doc, format, plain) {
        Ok(Prepared { omitted: None, .. }) => (true, vec!["Physical strands, end chemistry, qualities and full-extent annotations use a readable PopGene comment extension. Standard readers see the actual top strand and may ignore this extension; external rewriting makes inconsistent metadata inactive.".into()]),
        Ok(Prepared { omitted: Some(reason), .. }) => (false, vec![format!("Physical fragment supplement is omitted: {reason}. Only ordinary supported top-strand records survive; cohesive geometry and full-extent metadata cannot be recovered from names or ancestry.")]),
        Err(error) => (false, vec![format!("Physical fragment supplement cannot be written: {error}")]),
    }
}

#[cfg(test)]
mod tests;
