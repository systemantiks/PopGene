//! Public, readable annotation-lineage metadata carried in ordinary comments.
//!
//! This is a PopGene text convention, not a native ancestry packet. Checksums
//! detect corruption, not authenticity. A checksum alone never activates it.
use crate::{Error, Result};
use popgene_core::{
    Document, Feature, FeatureDirection, FeatureLocation, Segment,
    lineage::{Axis, Binding, Registry},
};
use serde_json::Value;
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, BTreeSet};

use crate::SaveFormat;

const START_PREFIX: &str = "[PopGene lineage v";
const START_SUFFIX: &str = " begin]";
const END: &str = "[PopGene lineage end]";
const VERSION: u32 = 1;
const MAX_PAYLOAD: usize = 4 * 1024 * 1024;
const MAX_DEPTH: usize = 96;
const CHUNK_BYTES: usize = 48;

/// Conservative marker candidates before ordinary comment whitespace folding.
/// Work per candidate stops after the fixed marker length; arbitrary comments
/// are never copied into a normalized string. The literal word prefix cannot
/// be introduced by whitespace folding.
pub(super) fn marker_candidates(comments: &[String], start: &str, end: &str) -> usize {
    comments
        .iter()
        .map(|comment| {
            comment
                .match_indices("[PopGene")
                .filter(|(index, _)| {
                    [start, end].iter().any(|target| {
                        let mut input = comment[*index..].chars().filter(|ch| !ch.is_whitespace());
                        target
                            .chars()
                            .filter(|ch| !ch.is_whitespace())
                            .all(|expected| input.next() == Some(expected))
                    })
                })
                .count()
        })
        .sum()
}

fn bad(message: &str) -> Error {
    Error::Message(format!("Annotation lineage: {message}"))
}
fn hash(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

/// Every whitespace character belonging to a JSON string is escaped, including
/// Unicode separators. Consequently ASCII wrapping can be removed at any byte
/// boundary without changing a value. Non-ASCII whitespace in the physical
/// carrier is refused instead of being silently reinterpreted as wrapping.
fn compact_payload(value: &Value) -> Result<String> {
    let text = serde_json::to_string(value)?;
    let mut output = String::with_capacity(text.len().min(MAX_PAYLOAD));
    for ch in text.chars() {
        if ch.is_whitespace() {
            for unit in ch.encode_utf16(&mut [0; 2]) {
                use std::fmt::Write as _;
                write!(output, "\\u{unit:04x}").unwrap();
            }
        } else {
            output.push(ch);
        }
        if output.len() > MAX_PAYLOAD {
            return Err(bad("readable provenance exceeds the 4 MiB payload limit"));
        }
    }
    depth_check(&output)?;
    Ok(output)
}
fn block(value: &Value) -> Result<String> {
    let payload = compact_payload(value)?;
    let header = format!(
        "bytes={};sha256={};json={payload}",
        payload.len(),
        hash(payload.as_bytes())
    );
    let mut output = format!("{START_PREFIX}{VERSION}{START_SUFFIX}\n");
    let mut width = 0;
    for ch in header.chars() {
        if width + ch.len_utf8() > CHUNK_BYTES {
            output.push('\n');
            width = 0;
        }
        output.push(ch);
        width += ch.len_utf8();
    }
    output.push('\n');
    output.push_str(END);
    Ok(output)
}
fn depth_check(text: &str) -> Result<()> {
    let mut depth = 0_usize;
    let mut quoted = false;
    let mut escaped = false;
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

struct ParsedBlock {
    comment: usize,
    start: usize,
    end: usize,
    value: Value,
}

/// Parse exactly one marker pair, without deleting or mutating source text.
/// Invalid, colliding and unsupported carriers must stay ordinary readable
/// comments. The caller decides whether to activate a semantically valid value.
fn parse(comments: &[String]) -> Result<Option<ParsedBlock>> {
    let mut found = None;
    let mut starts = 0_usize;
    let mut ends = 0_usize;
    for (index, comment) in comments.iter().enumerate() {
        for (start, _) in comment.match_indices(START_PREFIX) {
            starts += 1;
            if starts > 1 {
                return Err(bad(
                    "multiple provenance blocks or colliding markers; text retained without active ancestry",
                ));
            }
            found = Some((index, start));
        }
        ends = ends.saturating_add(comment.matches(END).count());
    }
    let Some((index, start)) = found else {
        return if ends == 0 {
            Ok(None)
        } else {
            Err(bad("unmatched provenance end marker; text retained"))
        };
    };
    if ends != 1 {
        return Err(bad(
            "missing or repeated provenance end marker; text retained",
        ));
    }
    let comment = &comments[index];
    let version_start = start + START_PREFIX.len();
    let suffix = comment[version_start..]
        .find(START_SUFFIX)
        .filter(|offset| *offset <= 10)
        .ok_or_else(|| bad("malformed provenance version marker; text retained"))?;
    if comment[version_start..version_start + suffix] != VERSION.to_string() {
        return Err(bad(
            "unsupported provenance version; text retained without active ancestry",
        ));
    }
    let payload_start = version_start + suffix + START_SUFFIX.len();
    let relative_end = comment[payload_start..]
        .find(END)
        .ok_or_else(|| bad("provenance markers span separate comments; text retained"))?;
    let payload_end = payload_start + relative_end;
    // No allocation proportional to arbitrary source comments: only the one
    // bounded carrier is compacted, and it is not decompressed or decoded.
    let mut compact = String::new();
    for ch in comment[payload_start..payload_end].chars() {
        if ch.is_ascii_whitespace() {
            continue;
        }
        if ch.is_whitespace() {
            return Err(bad(
                "non-ASCII wrapping whitespace is unsupported; text retained",
            ));
        }
        compact.push(ch);
        if compact.len() > MAX_PAYLOAD + 128 {
            return Err(bad("provenance payload exceeds 4 MiB; text retained"));
        }
    }
    let (length, rest) = compact
        .strip_prefix("bytes=")
        .and_then(|text| text.split_once(";sha256="))
        .ok_or_else(|| bad("malformed provenance length/checksum header; text retained"))?;
    let length = length
        .parse::<usize>()
        .ok()
        .filter(|length| *length <= MAX_PAYLOAD)
        .ok_or_else(|| bad("invalid provenance payload length; text retained"))?;
    let (checksum, payload) = rest
        .split_once(";json=")
        .ok_or_else(|| bad("missing provenance JSON; text retained"))?;
    if length != payload.len()
        || checksum.len() != 64
        || !checksum
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
        || checksum != hash(payload.as_bytes())
    {
        return Err(bad(
            "provenance length or checksum mismatch; text retained without active ancestry",
        ));
    }
    depth_check(payload)?;
    let value = serde_json::from_str(payload)?;
    Ok(Some(ParsedBlock {
        comment: index,
        start,
        end: payload_end + END.len(),
        value,
    }))
}

fn format_name(format: SaveFormat) -> &'static str {
    match format {
        SaveFormat::SnapGene => "native",
        SaveFormat::GenBank => "genbank",
        SaveFormat::Embl => "embl",
        _ => "unsupported",
    }
}
fn geometry(doc: &Document) -> (usize, Segment, Segment) {
    if let Some(ends) = &doc.duplex_ends {
        let g = ends.geometry(doc.sequence.len());
        (g.extent, g.top, g.bottom)
    } else {
        let n = doc.sequence.len();
        (
            n,
            Segment { start: 0, end: n },
            Segment {
                start: 0,
                end: if doc.single_stranded { 0 } else { n },
            },
        )
    }
}
fn physical_sequence(doc: &Document) -> Result<String> {
    Ok(match &doc.duplex_ends {
        Some(ends) => ends.reference_sequence(&doc.sequence)?,
        None => doc.sequence.to_string(),
    })
}
fn context(doc: &Document) -> Result<crate::canonical_json::ContextHash> {
    // Field order also reconstructs the one pre-canonical preserve_order hash.
    #[derive(serde::Serialize)]
    struct Context<'a> {
        sequence_sha256: String,
        physical_reference_sha256: String,
        geometry: (usize, Segment, Segment),
        molecule: popgene_core::Molecule,
        single_stranded: bool,
        circular: bool,
        features: &'a [Feature],
        duplex_features: &'a [Feature],
        sequence_colors: &'a [popgene_core::metadata::SequenceColor],
        duplex_colors: &'a [popgene_core::metadata::SequenceColor],
        references: &'a [popgene_core::metadata::Reference],
    }
    crate::canonical_json::context(&Context {
        sequence_sha256: hash(doc.sequence.as_bytes()),
        physical_reference_sha256: hash(physical_sequence(doc)?.as_bytes()),
        geometry: geometry(doc),
        molecule: doc.molecule,
        single_stranded: doc.single_stranded,
        circular: doc.circular,
        features: &doc.features,
        duplex_features: &doc.duplex_annotations.features,
        sequence_colors: &doc.sequence_colors,
        duplex_colors: &doc.duplex_annotations.sequence_colors,
        references: &doc.references,
    })
}
/// Bind to the actual serialized annotation records as well as their parsed
/// meaning. Reordering or rewriting records externally invalidates the carrier;
/// it is never repaired by matching otherwise identical feature fingerprints.
pub(super) fn records_hash(bytes: &[u8], format: SaveFormat) -> Result<String> {
    let mut digest = Sha256::new();
    if format == SaveFormat::SnapGene {
        for (kind, body) in crate::snapgene_preserve::packets(bytes)? {
            if matches!(kind, 5 | 10) {
                digest.update([kind]);
                digest.update((body.len() as u64).to_be_bytes());
                digest.update(body);
            }
        }
    } else {
        let text =
            std::str::from_utf8(bytes).map_err(|_| bad("annotation records are not UTF-8"))?;
        let mut in_features = false;
        for line in text.split_inclusive('\n') {
            let include = if format == SaveFormat::GenBank {
                if line.starts_with("FEATURES") {
                    in_features = true;
                } else if line.starts_with("ORIGIN") || line.starts_with("CONTIG") {
                    in_features = false;
                }
                in_features
            } else {
                line.starts_with("FT")
            };
            if include {
                digest.update(line.as_bytes());
            }
        }
    }
    Ok(format!("{:x}", digest.finalize()))
}

fn features(doc: &Document) -> impl Iterator<Item = (Axis, &Feature)> {
    doc.features.iter().map(|f| (Axis::Top, f)).chain(
        doc.duplex_annotations
            .features
            .iter()
            .map(|f| (Axis::Duplex, f)),
    )
}
struct ProjectionView<'a> {
    doc: &'a Document,
    sequence: String,
    features: BTreeMap<(Axis, u64), (usize, &'a Feature)>,
}
impl<'a> ProjectionView<'a> {
    fn new(doc: &'a Document) -> Result<Self> {
        Ok(Self {
            doc,
            sequence: physical_sequence(doc)?,
            features: features(doc)
                .enumerate()
                .map(|(index, (axis, f))| ((axis, f.id), (index, f)))
                .collect(),
        })
    }
}
fn full_feature(doc: &Document, axis: Axis, f: &Feature) -> Result<Feature> {
    let mut full = f.clone();
    if axis == Axis::Top {
        popgene_core::workflow::shift_feature(&mut full, geometry(doc).1.start as isize)?;
    }
    Ok(full)
}
fn effective_location(f: &Feature) -> FeatureLocation {
    f.location.clone().unwrap_or_else(|| {
        let parts = f
            .parts
            .iter()
            .map(|p| FeatureLocation::Range {
                start: p.start,
                end: p.end,
                before: false,
                after: false,
            })
            .collect::<Vec<_>>();
        let location = if parts.len() == 1 {
            parts[0].clone()
        } else {
            FeatureLocation::Join(parts)
        };
        if f.reverse {
            FeatureLocation::Complement(Box::new(location))
        } else {
            location
        }
    })
}
fn exact(location: &FeatureLocation) -> bool {
    match location {
        FeatureLocation::Range { before, after, .. } => !before && !after,
        FeatureLocation::Complement(child) => exact(child),
        FeatureLocation::Join(children) => children.iter().all(exact),
        _ => false,
    }
}
fn same_location(a: &Feature, b: &Feature) -> bool {
    let a = effective_location(a);
    let b = effective_location(b);
    a == b
        || (exact(&a)
            && exact(&b)
            && a.uniform_projection().is_some()
            && a.uniform_projection() == b.uniform_projection())
}
fn direction(feature: &Feature) -> FeatureDirection {
    feature.display_direction().unwrap_or(if feature.reverse {
        FeatureDirection::Reverse
    } else {
        FeatureDirection::Forward
    })
}
fn referenced<'a>(
    doc: &'a Document,
    feature: &Feature,
) -> Option<Vec<&'a popgene_core::metadata::Reference>> {
    feature
        .qualifiers
        .iter()
        .filter(|(key, _)| key == "citation")
        .map(|(_, value)| {
            value
                .as_ref()?
                .trim()
                .trim_start_matches('[')
                .trim_end_matches(']')
                .parse::<usize>()
                .ok()?
                .checked_sub(1)
                .and_then(|index| doc.references.get(index))
        })
        .collect()
}
const REVIEW_NOTE: &str = "Sequence changed in PopGene; annotation needs review.";
fn conventional_binding(feature: &Feature, key: &str, value: Option<&str>) -> bool {
    match key {
        "binding_boundStrand" => value == Some(if feature.reverse { "1" } else { "0" }),
        "binding_annealedBases" => value.is_some() && value == feature.primer_sequence().as_deref(),
        "binding_location" => match (feature.parts.first(), feature.parts.last()) {
            (Some(first), Some(last)) if last.end > 0 => {
                value == Some(format!("{}-{}", first.start, last.end - 1).as_str())
            }
            _ => false,
        },
        _ => false,
    }
}
fn equivalent_display_value(feature: &Feature, key: &str, value: Option<&str>) -> bool {
    if !matches!(key, "snapgene_color" | "snapgene_directionality") {
        return false;
    }
    let Some(value) = value.filter(|value| value.len() <= 7) else {
        return false;
    };
    // Interpret this one value through the shared vocabulary. A valid effective
    // display value elsewhere must not authorize discarding an unknown or
    // conflicting duplicate. Known color names and #RRGGBB fit this small bound.
    let isolated = Feature {
        id: 0,
        label: String::new(),
        kind: String::new(),
        parts: vec![],
        reverse: false,
        location: None,
        qualifiers: vec![(key.into(), Some(value.into()))],
        edited: false,
    };
    match key {
        "snapgene_color" => isolated
            .display_color_rgb()
            .is_some_and(|color| feature.display_color_rgb() == Some(color)),
        "snapgene_directionality" => isolated
            .display_direction()
            .is_some_and(|direction| feature.display_direction() == Some(direction)),
        _ => false,
    }
}
fn qualifiers(feature: &Feature, references_valid: bool) -> Vec<(String, Option<String>)> {
    let terminal = feature.terminal_display_note();
    let mut result = feature
        .qualifiers
        .iter()
        .filter(|(key, value)| {
            !((key == "label" && value.as_deref() == Some(feature.label.as_str()))
                || equivalent_display_value(feature, key, value.as_deref())
                || key == "snapgene_recentID"
                || conventional_binding(feature, key, value.as_deref())
                || (key == "citation" && references_valid)
                || (feature.kind == "CDS"
                    && matches!(key.as_str(), "codon_start" | "snapgene_readingFrame")
                    && popgene_core::translation::cds_phase(feature).is_ok())
                || (key == "note" && value.as_deref().is_some_and(|v| Some(v) == terminal)))
        })
        .cloned()
        .collect::<Vec<_>>();
    result.sort();
    result
}
fn same_feature(
    before: &ProjectionView<'_>,
    after: &ProjectionView<'_>,
    binding: &Binding,
) -> Result<bool> {
    let Some(&(before_index, original)) =
        before.features.get(&(binding.from_axis, binding.from_id))
    else {
        return Ok(false);
    };
    let Some(&(after_index, current)) = after.features.get(&(binding.to_axis, binding.to_id))
    else {
        return Ok(false);
    };
    let a = full_feature(before.doc, binding.from_axis, original)?;
    let b = full_feature(after.doc, binding.to_axis, current)?;
    if a.label != b.label
        || a.kind != b.kind
        || a.parts != b.parts
        || a.reverse != b.reverse
        || !same_location(&a, &b)
        || direction(&a) != direction(&b)
        || original.effective_color_rgb(before_index) != current.effective_color_rgb(after_index)
        || a.display_cleavage_sites() != b.display_cleavage_sites()
    {
        return Ok(false);
    }
    let segments = |f: &Feature, index| {
        f.display_segments()
            .into_iter()
            .map(|mut s| {
                s.color_rgb = Some(s.color_rgb.unwrap_or_else(|| f.effective_color_rgb(index)));
                s
            })
            .collect::<Vec<_>>()
    };
    if segments(&a, before_index) != segments(&b, after_index) {
        return Ok(false);
    }
    if a.kind == "CDS"
        && popgene_core::translation::cds_phase(&a).ok()
            != popgene_core::translation::cds_phase(&b).ok()
    {
        return Ok(false);
    }
    let refs_a = referenced(before.doc, &a);
    let refs_b = referenced(after.doc, &b);
    if refs_a != refs_b {
        return Ok(false);
    }
    let qa = qualifiers(&a, refs_a.is_some());
    let mut qb = qualifiers(&b, refs_b.is_some());
    if a.edited && !b.edited {
        // Exactly one conventional review note is an established transport of
        // the in-memory review flag. Existing identical notes remain counted.
        let note = ("note".to_owned(), Some(REVIEW_NOTE.to_owned()));
        if let Some(index) = qb.iter().position(|q| q == &note) {
            qb.remove(index);
        } else {
            return Ok(false);
        }
    } else if a.edited != b.edited {
        return Ok(false);
    }
    let extracted_a = effective_location(&a)
        .extract(&before.sequence, before.doc.molecule)
        .map_err(|_| bad("tracked source extraction is unresolved"))?;
    let extracted_b = effective_location(&b)
        .extract(&after.sequence, after.doc.molecule)
        .map_err(|_| bad("tracked target extraction is unresolved"))?;
    Ok(qa == qb && extracted_a == extracted_b)
}

fn bindings(before: &Document, after: &Document, format: SaveFormat) -> Result<Vec<Binding>> {
    let ids = if format == SaveFormat::SnapGene {
        crate::snapgene::encoded_feature_ids(before)?
    } else {
        before
            .features
            .iter()
            .chain(
                before
                    .duplex_annotations
                    .features
                    .iter()
                    .filter(|_| after.duplex_ends.is_some()),
            )
            .map(|f| f.id)
            .collect()
    };
    if ids.len() != features(after).count() {
        return Err(bad("codec annotation occurrence count changed"));
    }
    let mut seen = BTreeSet::new();
    let from_axes = features(before)
        .map(|(axis, f)| (f.id, axis))
        .collect::<BTreeMap<_, _>>();
    let to_axes = features(after)
        .map(|(axis, f)| (f.id, axis))
        .collect::<BTreeMap<_, _>>();
    ids.into_iter()
        .enumerate()
        .map(|(index, from_id)| {
            if !seen.insert(from_id) {
                return Err(bad("codec annotation occurrence association is ambiguous"));
            }
            let from_axis = from_axes
                .get(&from_id)
                .copied()
                .ok_or_else(|| bad("codec source annotation occurrence is missing"))?;
            let to_id = index as u64 + 1;
            let to_axis = to_axes
                .get(&to_id)
                .copied()
                .ok_or_else(|| bad("codec target annotation occurrence is missing"))?;
            Ok(Binding {
                from_axis,
                from_id,
                to_axis,
                to_id,
            })
        })
        .collect()
}

/// Decode after ordinary biological records are parsed. Invalid text is retained
/// verbatim and never changes the active registry or ordinary record import.
pub(super) fn read(doc: &mut Document, bytes: &[u8], format: SaveFormat) {
    let restored = (|| -> Result<Option<(ParsedBlock, Registry)>> {
        let Some(parsed) = parse(&doc.comments)? else {
            return Ok(None);
        };
        let object = parsed
            .value
            .as_object()
            .ok_or_else(|| bad("payload is not an object"))?;
        if object.len() != 4
            || object.get("format") != Some(&Value::from(format_name(format)))
            || object.get("records_sha256") != Some(&Value::from(records_hash(bytes, format)?))
            || !context(doc)?.matches(object.get("context_sha256"))
        {
            return Err(bad(
                "record association or canonical context changed; readable provenance is inactive",
            ));
        }
        let wire = object
            .get("registry")
            .ok_or_else(|| bad("registry is missing"))?;
        let registry: Registry = serde_json::from_value(wire.clone())?;
        if serde_json::to_value(&registry)? != *wire {
            return Err(bad(
                "unknown or noncanonical registry fields; readable provenance is inactive",
            ));
        }
        registry.validate(doc)?;
        let (extent, top, bottom) = geometry(doc);
        for piece in &registry.pieces {
            if !registry.actionable_piece(doc, piece)
                || piece.runs.iter().any(|run| {
                    let Some(end) = run.target_start.checked_add(run.length) else {
                        return true;
                    };
                    end > extent
                        || (run.strands & 1 != 0 && (run.target_start < top.start || end > top.end))
                        || (run.strands & 2 != 0
                            && (run.target_start < bottom.start || end > bottom.end))
                })
            {
                return Err(bad(
                    "piece baseline, sequence, references or physical mapping changed; readable provenance is inactive",
                ));
            }
        }
        Ok(Some((parsed, registry)))
    })();
    match restored {
        Ok(Some((parsed, registry))) => {
            let comment = &mut doc.comments[parsed.comment];
            comment.replace_range(parsed.start..parsed.end, "");
            if comment.ends_with('\n') {
                comment.pop();
            }
            if comment.is_empty() {
                doc.comments.remove(parsed.comment);
            }
            if registry.tracking_limited {
                doc.warnings.push("Annotation ancestry was limited or omitted by an earlier operation or save. Untracked annotations are not inferred to be complete ancestors.".into());
            }
            doc.annotation_lineage = registry;
        }
        Ok(None) => {}
        Err(error) => doc.warnings.push(format!(
            "{error}. Original provenance text is retained; no ancestry was activated."
        )),
    }
}

pub(super) fn write(
    doc: &Document,
    format: SaveFormat,
    plain: fn(&Document) -> Result<Vec<u8>>,
) -> Result<Vec<u8>> {
    Ok(prepare(doc, format, plain)?.bytes)
}

struct Prepared {
    bytes: Vec<u8>,
    omitted: Option<Omitted>,
}

enum Omitted {
    LimitedMarker(String),
    NoMarker {
        reason: String,
        marker_reason: String,
    },
}

fn prepare(
    doc: &Document,
    format: SaveFormat,
    plain: fn(&Document) -> Result<Vec<u8>>,
) -> Result<Prepared> {
    let retained_seed = doc.annotation_lineage.is_empty()
        && doc.source.as_ref().is_some_and(|source| {
            source
                .bytes
                .windows(START_PREFIX.len())
                .any(|bytes| bytes == START_PREFIX.as_bytes())
                && crate::import(&source.file_name, &source.bytes).is_ok_and(|original| {
                    original.annotation_lineage.instance_id == doc.annotation_lineage.instance_id
                        && !original.comments.iter().any(|c| c.contains(START_PREFIX))
                })
        });
    if doc.annotation_lineage.is_empty() && !retained_seed {
        let bytes = plain(doc)?;
        let candidate_count = |doc: &Document| marker_candidates(&doc.comments, START_PREFIX, END);
        let candidates = candidate_count(doc);
        if candidates > 0 {
            let reopened =
                crate::import(&format!("lineage.{}", format.extension_for(doc)), &bytes)?;
            if !reopened.annotation_lineage.is_empty() || candidate_count(&reopened) < candidates {
                return Err(bad(
                    "previously inactive provenance text would become active after serialization; preserve it separately or remove it explicitly before saving",
                ));
            }
        }
        return Ok(Prepared {
            bytes,
            omitted: None,
        });
    }
    if marker_candidates(&doc.comments, START_PREFIX, END) > 0 {
        return Err(bad(
            "existing provenance text cannot be replaced safely; preserve it separately or remove it explicitly before writing active lineage",
        ));
    }
    doc.annotation_lineage.validate(doc)?;
    let mut base = doc.clone();
    base.annotation_lineage = Registry::default();
    let first = plain(&base)?;
    let canonical = crate::import(&format!("lineage.{}", format.extension_for(doc)), &first)?;
    let attempted = projected_registry(doc, &canonical, format)
        .and_then(|registry| write_registry(&base, &canonical, &first, format, plain, &registry));
    match attempted {
        Ok(bytes) => Ok(Prepared {
            bytes,
            omitted: None,
        }),
        Err(reason) => {
            // Optional evidence must not disable an otherwise valid conversion.
            // Keep a small, explicit limitation marker where the target can fit
            // it; never rebase failed mappings or resurrect dropped descriptors.
            let limited = Registry {
                instance_id: doc.annotation_lineage.instance_id.clone(),
                tracking_limited: true,
                ..Registry::default()
            };
            match write_registry(&base, &canonical, &first, format, plain, &limited) {
                Ok(bytes) => Ok(Prepared {
                    bytes,
                    omitted: Some(Omitted::LimitedMarker(reason.to_string())),
                }),
                Err(marker_reason) => Ok(Prepared {
                    bytes: first,
                    omitted: Some(Omitted::NoMarker {
                        reason: reason.to_string(),
                        marker_reason: marker_reason.to_string(),
                    }),
                }),
            }
        }
    }
}

fn projected_registry(
    doc: &Document,
    canonical: &Document,
    format: SaveFormat,
) -> Result<Registry> {
    if doc.molecule != canonical.molecule
        || doc.circular != canonical.circular
        || doc.single_stranded != canonical.single_stranded
        || doc.sequence != canonical.sequence
        || geometry(doc) != geometry(canonical)
    {
        return Err(bad(
            "target loses sequence or physical duplex geometry needed by annotation ancestry; this text carrier cannot replace missing end geometry",
        ));
    }
    let before_view = ProjectionView::new(doc)?;
    let after_view = ProjectionView::new(canonical)?;
    if before_view.sequence != after_view.sequence {
        return Err(bad(
            "target changes the physical residues needed by annotation ancestry",
        ));
    }
    let bindings = bindings(doc, canonical, format)?;
    for piece in &doc.annotation_lineage.pieces {
        if !doc.annotation_lineage.actionable_piece(doc, piece) {
            continue;
        }
        let binding = bindings
            .iter()
            .find(|b| b.from_axis == piece.axis && b.from_id == piece.feature_id)
            .ok_or_else(|| bad("an active piece has no emitted annotation occurrence"))?;
        if !same_feature(&before_view, &after_view, binding)? {
            return Err(bad(
                "target changes tracked annotation location, phase, references, display or qualifiers; provenance cannot be rebound to altered data",
            ));
        }
    }
    Ok(Registry::transport_for_codec(doc, canonical, &bindings)?)
}

fn write_registry(
    base: &Document,
    canonical: &Document,
    first: &[u8],
    format: SaveFormat,
    plain: fn(&Document) -> Result<Vec<u8>>,
    registry: &Registry,
) -> Result<Vec<u8>> {
    let value = serde_json::json!({
        "format": format_name(format),
        "records_sha256": records_hash(first, format)?,
        "context_sha256": context(canonical)?.canonical,
        "registry": registry,
    });
    let mut with_carrier = base.clone();
    with_carrier.comments.push(block(&value)?);
    let bytes = plain(&with_carrier)?;
    let reopened = crate::import(&format!("lineage.{}", format.extension_for(base)), &bytes)?;
    if &reopened.annotation_lineage != registry
        || reopened.annotation_lineage.instance_id != registry.instance_id
    {
        return Err(bad(
            "written provenance did not pass independent import validation",
        ));
    }
    Ok(bytes)
}

pub(super) fn warnings(doc: &Document, format: SaveFormat) -> Vec<String> {
    if doc.annotation_lineage.is_empty() {
        return vec![];
    }
    if crate::save::original_matches(doc, crate::save::source_format(format)) {
        return vec![];
    }
    if !matches!(
        format,
        SaveFormat::SnapGene | SaveFormat::GenBank | SaveFormat::Embl
    ) {
        return vec!["Annotation ancestry is omitted by this format. Reopening or the canonical saved document cannot reconstruct ancestry from matching sequence or names.".into()];
    }
    let mut result = vec!["Annotation ancestry uses an explicit, readable PopGene provenance convention in comments. It is not native ancestry support; third-party rewriting or removing annotation records or provenance makes it inactive. Checksums establish integrity, not authenticity.".into()];
    let plain = match format {
        SaveFormat::SnapGene => crate::snapgene::write_plain,
        SaveFormat::GenBank => crate::genbank_bytes_physical,
        SaveFormat::Embl => crate::embl::write_physical,
        _ => unreachable!(),
    };
    match prepare(doc, format, plain) {
        Ok(prepared) => {
            match prepared.omitted {
                Some(Omitted::LimitedMarker(reason)) => result.push(format!("Active annotation ancestry is omitted from this saved representation: {reason}. Supported biological records are still saved, with an explicit ancestry-limited marker; no failed baseline is refreshed.")),
                Some(Omitted::NoMarker { reason, marker_reason }) => result.push(format!("Active annotation ancestry is omitted from this saved representation: {reason}. Supported biological records are still saved, but no ancestry marker could be stored: {marker_reason}. Reopening cannot retain this tracking limitation or the source identity; no failed baseline is refreshed.")),
                None => {}
            }
        }
        Err(reason) => result.push(format!(
            "Provenance cannot be written while the target conversion is refused: {reason}"
        )),
    }
    let inactive = doc
        .annotation_lineage
        .pieces
        .iter()
        .filter(|p| !doc.annotation_lineage.actionable_piece(doc, p))
        .count();
    if inactive > 0 {
        result.push(format!("{inactive} edited or inconsistent ancestry piece(s) are not carried as active provenance; saving will not revalidate their baselines."));
    }
    result
}

#[cfg(test)]
mod text_tests {
    use super::*;

    #[test]
    fn readable_unicode_payload_survives_ascii_line_wrapping() {
        let original = serde_json::json!({
            "label": "α promoter\u{a0}one\tline\nnext  two spaces",
            "note": "literal \\u0020 and [PopGene lineage end] and \"quotes\"",
            "large_label": "𐐀DNA".repeat(50),
        });
        let encoded = block(&original).unwrap();
        assert!(encoded.contains("label"));
        assert!(encoded.contains("\\u00a0"));
        assert!(
            encoded
                .lines()
                .skip(1)
                .all(|line| line.len() <= CHUNK_BYTES)
        );
        // Ordinary flat-file readers normalize comments to single spaces.
        let wrapped = encoded.split_whitespace().collect::<Vec<_>>().join(" ");
        let comments = vec!["ordinary source note".into(), wrapped];
        let read = parse(&comments).unwrap().unwrap();
        assert_eq!(read.value, original);
        assert_eq!(read.comment, 1);
        assert_eq!(read.start, 0);
        assert_eq!(read.end, comments[1].len());
    }

    #[test]
    fn corrupt_unknown_colliding_and_deep_text_is_not_accepted() {
        let encoded = block(&serde_json::json!({"origin":"distinct"})).unwrap();
        for bad_text in [
            encoded.replacen("v1 begin", "v2 begin", 1),
            encoded.replacen("distinct", "modified", 1),
            encoded.replace(END, ""),
            format!("{encoded}\n{encoded}"),
            encoded.replacen("\n", "\u{a0}", 1),
        ] {
            let comments = vec![bad_text.clone()];
            assert!(parse(&comments).is_err());
            assert_eq!(
                comments[0], bad_text,
                "failure must not mutate readable text"
            );
        }
        assert!(depth_check(&format!("{}0{}", "[".repeat(97), "]".repeat(97))).is_err());
        assert!(parse(&[END.into()]).is_err());
    }
}

#[cfg(test)]
mod tests;
