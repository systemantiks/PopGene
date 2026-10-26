//! Bounded annotation ancestry for verified sequence constructions.
//!
//! The registry retains annotation descriptors and coordinate evidence, never
//! ancestor sequence or a recursive Document. A matching label is not evidence.
use crate::{
    Document, Error, Feature, FeatureLocation as L, Molecule, Segment, metadata::Reference,
};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, BTreeSet};

pub const VERSION: u32 = 1;
pub const MAX_ORIGINS: usize = 2048;
pub const MAX_PIECES: usize = 8192;
pub const MAX_RUNS: usize = 65_536;
pub const MAX_BYTES: usize = 8 * 1024 * 1024;
pub const MAX_MAPPED_BASES: usize = 2_000_000;
pub const MAX_VALIDATION_WORK: usize = 16_000_000;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum Axis {
    Top,
    Duplex,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Origin {
    pub id: String,
    pub source_id: String,
    pub occurrence_id: String,
    pub top: Segment,
    pub bottom: Segment,
    pub feature: Feature,
    pub references: Vec<Reference>,
    pub extent: usize,
    pub molecule: Molecule,
    pub circular: bool,
    pub source_digest: String,
    pub sequence_digest: String,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Run {
    /// DFS Range occurrence in the original AST, including repeated leaves.
    pub leaf: u32,
    /// Offset from this original Range's low coordinate, irrespective of strand.
    pub offset: usize,
    pub length: usize,
    /// Low coordinate on the current full physical extent.
    pub target_start: usize,
    pub reverse: bool,
    /// Actual physical residue coverage: top=1, opposite=2.
    pub strands: u8,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Piece {
    pub mapping_digest: String,
    pub origin_id: String,
    pub copy_id: String,
    pub feature_id: u64,
    pub axis: Axis,
    pub runs: Vec<Run>,
    pub baseline: Feature,
    pub references_digest: String,
    pub sequence_digest: String,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Registry {
    pub version: u32,
    #[serde(default)]
    pub tracking_limited: bool,
    /// Persist even when no pieces exist: later digests of a recovered source
    /// must share identity with earlier fragments. This is not a content hash.
    pub instance_id: String,
    #[serde(default)]
    pub feature_epochs: Vec<(u64, String)>,
    pub origins: Vec<Origin>,
    pub pieces: Vec<Piece>,
}
impl Default for Registry {
    fn default() -> Self {
        Self {
            version: VERSION,
            tracking_limited: false,
            instance_id: uuid::Uuid::new_v4().to_string(),
            feature_epochs: vec![],
            origins: vec![],
            pieces: vec![],
        }
    }
}
impl PartialEq for Registry {
    fn eq(&self, other: &Self) -> bool {
        self.version == other.version
            && self.tracking_limited == other.tracking_limited
            && self.origins == other.origins
            && self.pieces == other.pieces
            && (self.is_empty() && other.is_empty() || self.feature_epochs == other.feature_epochs)
            && (self.is_empty() && other.is_empty() || self.instance_id == other.instance_id)
    }
}
impl Eq for Registry {}
#[derive(Clone, Copy, Debug)]
pub struct Binding {
    pub from_axis: Axis,
    pub from_id: u64,
    pub to_axis: Axis,
    pub to_id: u64,
}
#[derive(Clone, Debug)]
pub(crate) struct Map {
    pub source: Segment,
    pub target_start: usize,
    pub reverse: bool,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Merge {
    Physical,
    VerifiedOverlap,
}
fn error(message: &str) -> Error {
    Error::Operation(format!("Annotation lineage: {message}"))
}
fn digest(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}
fn encoded<T: Serialize>(value: &T) -> Vec<u8> {
    serde_json::to_vec(value).expect("lineage has only serializable scalar fields")
}
fn hash<T: Serialize>(value: &T) -> String {
    digest(&encoded(value))
}
fn origin_id(origin: &Origin) -> String {
    hash(&(
        &origin.source_id,
        &origin.occurrence_id,
        &origin.feature,
        &origin.references,
        origin.extent,
        origin.molecule,
        origin.circular,
        &origin.source_digest,
        &origin.sequence_digest,
        &origin.top,
        &origin.bottom,
    ))
}
fn mapping_digest(piece: &Piece, doc: &Document) -> String {
    hash(&(
        &piece.runs,
        piece.axis,
        geometry(doc),
        doc.molecule,
        doc.circular,
        doc.single_stranded,
    ))
}
fn is_digest(value: &str) -> bool {
    value.len() == 64 && value.bytes().all(|b| b.is_ascii_hexdigit())
}
fn location(feature: &Feature) -> L {
    feature.location.clone().unwrap_or_else(|| {
        let mut leaves = feature
            .parts
            .iter()
            .map(|p| L::Range {
                start: p.start,
                end: p.end,
                before: false,
                after: false,
            })
            .collect::<Vec<_>>();
        let value = if leaves.len() == 1 {
            leaves.remove(0)
        } else {
            L::Join(leaves)
        };
        if feature.reverse {
            L::Complement(Box::new(value))
        } else {
            value
        }
    })
}
fn leaves(value: &L) -> Option<Vec<Segment>> {
    fn walk(value: &L, out: &mut Vec<Segment>) -> Option<()> {
        match value {
            L::Range { start, end, .. } => out.push(Segment {
                start: *start,
                end: *end,
            }),
            L::Complement(child) => walk(child, out)?,
            L::Join(children) => {
                for child in children {
                    walk(child, out)?;
                }
            }
            _ => return None,
        }
        Some(())
    }
    let mut out = vec![];
    walk(value, &mut out)?;
    (!out.is_empty()).then_some(out)
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
                end: if doc.single_stranded || doc.molecule == Molecule::Protein {
                    0
                } else {
                    n
                },
            },
        )
    }
}
fn feature(doc: &Document, axis: Axis, id: u64) -> Option<&Feature> {
    match axis {
        Axis::Top => &doc.features,
        Axis::Duplex => &doc.duplex_annotations.features,
    }
    .iter()
    .find(|f| f.id == id)
}
fn full_feature(doc: &Document, axis: Axis, f: &Feature) -> Result<Feature, Error> {
    let mut f = f.clone();
    if axis == Axis::Top {
        crate::workflow::shift_feature(&mut f, geometry(doc).1.start as isize)?;
    }
    Ok(f)
}
fn extract_current(doc: &Document, feature: &Feature) -> Option<String> {
    fn range(doc: &Document, start: usize, end: usize) -> Option<String> {
        let (extent, top, bottom) = geometry(doc);
        if start > end || end > extent {
            return None;
        }
        let mut result = String::with_capacity(end - start);
        for (a, b, is_top) in [
            (start, end.min(top.start), false),
            (start.max(top.start), end.min(top.end), true),
            (start.max(top.end), end, false),
        ] {
            if a >= b {
                continue;
            }
            if is_top {
                result.push_str(doc.sequence.get(a - top.start..b - top.start)?);
            } else {
                if a < bottom.start || b > bottom.end {
                    return None;
                }
                let opposite = &doc.duplex_ends.as_ref()?.opposite_sequence;
                result.push_str(&doc.complement(opposite.get(bottom.end - b..bottom.end - a)?));
            }
        }
        Some(result)
    }
    fn walk(doc: &Document, value: &L) -> Option<String> {
        match value {
            L::Range { start, end, .. } => range(doc, *start, *end),
            L::Complement(child) => Some(doc.complement(&walk(doc, child)?)),
            L::Join(children) => {
                let mut out = String::new();
                for child in children {
                    out.push_str(&walk(doc, child)?);
                }
                Some(out)
            }
            _ => None,
        }
    }
    walk(doc, &location(feature))
}
fn serialized_size<T: Serialize>(value: &T) -> usize {
    struct Count(usize);
    impl std::io::Write for Count {
        fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
            self.0 = self.0.saturating_add(b.len());
            Ok(b.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut count = Count(0);
    let _ = serde_json::to_writer(&mut count, value);
    count.0
}
fn cited_references(f: &Feature, refs: &[Reference]) -> Vec<Reference> {
    let indices = f
        .qualifiers
        .iter()
        .filter(|(k, _)| k == "citation")
        .filter_map(|(_, v)| {
            v.as_ref()?
                .trim()
                .trim_start_matches('[')
                .trim_end_matches(']')
                .parse::<usize>()
                .ok()?
                .checked_sub(1)
        })
        .filter(|i| *i < refs.len())
        .collect::<BTreeSet<_>>();
    let mut selected = vec![Reference::default(); indices.last().map_or(0, |i| i + 1)];
    for index in indices {
        selected[index] = refs[index].clone();
    }
    selected
}
fn references_digest(f: &Feature, refs: &[Reference]) -> String {
    let targets = f
        .qualifiers
        .iter()
        .filter(|(k, _)| k == "citation")
        .map(|(_, v)| {
            let reference = v
                .as_ref()
                .and_then(|v| {
                    v.trim()
                        .trim_start_matches('[')
                        .trim_end_matches(']')
                        .parse::<usize>()
                        .ok()
                })
                .and_then(|n| n.checked_sub(1))
                .and_then(|n| refs.get(n));
            (v, reference)
        })
        .collect::<Vec<_>>();
    hash(&targets)
}
impl Registry {
    pub(crate) fn record_feature(&mut self, id: u64) {
        self.feature_epochs.retain(|(old, _)| *old != id);
        if self.feature_epochs.len() >= MAX_PIECES {
            self.tracking_limited = true;
            return;
        }
        self.feature_epochs
            .push((id, uuid::Uuid::new_v4().to_string()));
    }
    pub fn is_empty(&self) -> bool {
        self.origins.is_empty() && self.pieces.is_empty() && !self.tracking_limited
    }
    pub fn snapshot_bytes(&self) -> usize {
        serialized_size(self)
    }
    /// Structural validation is separate from actionability: an edited annotation
    /// can remain a valid document without permitting ancestral restoration.
    pub fn validate(&self, doc: &Document) -> Result<(), Error> {
        if self.version != VERSION || uuid::Uuid::parse_str(&self.instance_id).is_err() {
            return Err(error("unsupported version or invalid source identity"));
        }
        if self.origins.len() > MAX_ORIGINS
            || self.pieces.len() > MAX_PIECES
            || self.pieces.iter().map(|p| p.runs.len()).sum::<usize>() > MAX_RUNS
            || self.snapshot_bytes() > MAX_BYTES
        {
            return Err(error("registry exceeds its bounded storage limit"));
        }
        if self.feature_epochs.len() > MAX_PIECES {
            return Err(error("feature occurrence identity limit exceeded"));
        }
        let mut epochs = BTreeSet::new();
        for (id, epoch) in &self.feature_epochs {
            if !epochs.insert(*id) || uuid::Uuid::parse_str(epoch).is_err() {
                return Err(error("invalid feature occurrence identity"));
            }
        }
        let mut ids = BTreeSet::new();
        let mut total = 0usize;
        let mut original_bases = 0usize;
        for origin in &self.origins {
            if !is_digest(&origin.id)
                || !is_digest(&origin.source_digest)
                || !is_digest(&origin.sequence_digest)
                || !ids.insert(&origin.id)
                || origin.extent == 0
                || origin.extent > crate::MAX_BASES
                || uuid::Uuid::parse_str(&origin.source_id).is_err()
                || !is_digest(&origin.occurrence_id)
                || origin.id != origin_id(origin)
                || origin.top.start > origin.top.end
                || origin.top.end > origin.extent
                || origin.bottom.start > origin.bottom.end
                || origin.bottom.end > origin.extent
            {
                return Err(error("invalid or duplicate origin"));
            }
            location(&origin.feature)
                .validate(origin.extent, origin.circular)
                .map_err(|_| error("invalid original location"))?;
            let ast = location(&origin.feature);
            let (parts, reverse) = ast
                .uniform_projection()
                .unwrap_or_else(|| (ast.local_parts(), false));
            if parts != origin.feature.parts || reverse != origin.feature.reverse {
                return Err(error("inconsistent original location projection"));
            }
            let ranges = leaves(&location(&origin.feature))
                .ok_or_else(|| error("unsupported original location operator"))?;
            for range in ranges {
                original_bases = original_bases
                    .checked_add(range.end - range.start)
                    .ok_or_else(|| error("original occurrence budget overflow"))?;
            }
            if original_bases > MAX_MAPPED_BASES {
                return Err(error("original occurrence work limit exceeded"));
            }
        }
        let mut current = BTreeSet::new();
        let mut work = 0usize;
        for piece in &self.pieces {
            let origin = self
                .origins
                .iter()
                .find(|o| o.id == piece.origin_id)
                .ok_or_else(|| error("piece has no origin"))?;
            let original = leaves(&location(&origin.feature))
                .ok_or_else(|| error("invalid original location"))?;
            work = work
                .checked_add(original.iter().map(|r| r.end - r.start).sum::<usize>())
                .ok_or_else(|| error("validation work overflow"))?;
            if work > MAX_VALIDATION_WORK {
                return Err(error("validation work limit exceeded"));
            }
            if !is_digest(&piece.mapping_digest)
                || uuid::Uuid::parse_str(&piece.copy_id).is_err()
                || !is_digest(&piece.references_digest)
                || !is_digest(&piece.sequence_digest)
                || !current.insert((piece.axis, piece.feature_id))
            {
                return Err(error("invalid or duplicate piece association"));
            }
            for run in &piece.runs {
                total = total
                    .checked_add(run.length)
                    .ok_or_else(|| error("mapping budget overflow"))?;
                let leaf = original
                    .get(run.leaf as usize)
                    .ok_or_else(|| error("missing original leaf"))?;
                if run.length == 0
                    || run.strands == 0
                    || run.strands > 3
                    || run
                        .offset
                        .checked_add(run.length)
                        .is_none_or(|v| v > leaf.end - leaf.start)
                    || run
                        .target_start
                        .checked_add(run.length)
                        .is_none_or(|v| v > crate::MAX_BASES * 2)
                {
                    return Err(error("invalid occurrence mapping"));
                }
            }
        }
        if total > MAX_MAPPED_BASES {
            return Err(error("mapping work limit exceeded"));
        }
        let _ = doc; // Current edits may make a structurally valid piece inactive.
        Ok(())
    }
    pub fn actionable_piece(&self, doc: &Document, piece: &Piece) -> bool {
        let Some(f) = feature(doc, piece.axis, piece.feature_id) else {
            return false;
        };
        if piece.mapping_digest != mapping_digest(piece, doc)
            || f != &piece.baseline
            || references_digest(f, &doc.references) != piece.references_digest
        {
            return false;
        }
        let Ok(f) = full_feature(doc, piece.axis, f) else {
            return false;
        };

        let (extent, top, bottom) = geometry(doc);
        let mut mapped = Vec::new();
        let mut source_seen = BTreeSet::new();
        for run in &piece.runs {
            if run.length > MAX_MAPPED_BASES
                || mapped.len().saturating_add(run.length) > MAX_MAPPED_BASES
                || run
                    .target_start
                    .checked_add(run.length)
                    .is_none_or(|end| end > extent)
            {
                return false;
            }
            for i in 0..run.length {
                let target = run.target_start + i;
                if run.strands & 1 != 0 && !(top.start..top.end).contains(&target)
                    || run.strands & 2 != 0 && !(bottom.start..bottom.end).contains(&target)
                    || !source_seen.insert((run.leaf, run.offset + i))
                {
                    return false;
                }
                mapped.push(target);
            }
        }
        let expected_len = f.parts.iter().map(|p| p.end - p.start).sum::<usize>();
        if expected_len != mapped.len() {
            return false;
        }
        let Some(origin) = self.origins.iter().find(|o| o.id == piece.origin_id) else {
            return false;
        };
        let Some(original_path) = biological_path(&location(&origin.feature)) else {
            return false;
        };
        let mut source_map = BTreeMap::new();
        for run in &piece.runs {
            for i in 0..run.length {
                source_map.insert(
                    (run.leaf as usize, run.offset + i),
                    (
                        run.target_start + if run.reverse { run.length - 1 - i } else { i },
                        run.reverse,
                    ),
                );
            }
        }
        let observed = original_path
            .into_iter()
            .filter_map(|(leaf, offset, strand)| {
                source_map
                    .get(&(leaf, offset))
                    .map(|&(target, reversed)| (target, strand ^ reversed))
            })
            .collect::<Vec<_>>();
        let Some(current_ranges) = leaves(&location(&f)) else {
            return false;
        };
        let Some(current_path) = biological_path(&location(&f)) else {
            return false;
        };
        let expected = current_path
            .into_iter()
            .map(|(leaf, offset, strand)| (current_ranges[leaf].start + offset, strand))
            .collect::<Vec<_>>();
        expected == observed
            && extract_current(doc, &f)
                .is_some_and(|s| digest(s.as_bytes()) == piece.sequence_digest)
    }
    /// Caller proves codec semantics; this bridge only accepts explicit unique
    /// record associations after validating the pre-export current baselines.
    pub fn transport_for_codec(
        before: &Document,
        canonical: &Document,
        bindings: &[Binding],
    ) -> Result<Self, Error> {
        before.annotation_lineage.validate(before)?;
        let mut result = before.annotation_lineage.clone();
        let mut seen = BTreeSet::new();
        result.pieces.clear();
        result.feature_epochs = before
            .annotation_lineage
            .feature_epochs
            .iter()
            .filter_map(|(id, epoch)| {
                bindings
                    .iter()
                    .find(|b| b.from_id == *id)
                    .map(|b| (b.to_id, epoch.clone()))
            })
            .collect();
        for piece in &before.annotation_lineage.pieces {
            if !before.annotation_lineage.actionable_piece(before, piece) {
                result.tracking_limited = true;
                continue;
            }
            let matches = bindings
                .iter()
                .filter(|b| b.from_axis == piece.axis && b.from_id == piece.feature_id)
                .collect::<Vec<_>>();
            if matches.len() != 1 {
                return Err(error(
                    "codec needs exactly one binding for each active piece",
                ));
            }
            let b = matches[0];
            if !seen.insert((b.to_axis, b.to_id))
                || feature(canonical, b.to_axis, b.to_id).is_none()
            {
                return Err(error("codec binding is missing or duplicated"));
            }
            let mut p = piece.clone();
            p.axis = b.to_axis;
            p.feature_id = b.to_id;
            result.pieces.push(p);
        }
        if result.pieces.iter().any(|p| {
            feature(canonical, p.axis, p.feature_id)
                .is_none_or(|f| f.len() != p.runs.iter().map(|r| r.length).sum::<usize>())
        }) {
            return Err(error("codec changed mapped occurrence coverage"));
        }
        refresh(&mut result, canonical)?;
        if result
            .pieces
            .iter()
            .any(|p| !result.actionable_piece(canonical, p))
        {
            return Err(error("codec changed mapped occurrence order or geometry"));
        }
        result.prune();
        result.validate(canonical)?;
        Ok(result)
    }
    fn prune(&mut self) {
        self.origins
            .retain(|o| self.pieces.iter().any(|p| p.origin_id == o.id));
    }
}
fn refresh(registry: &mut Registry, doc: &Document) -> Result<(), Error> {
    registry.pieces.retain_mut(|piece| {
        let axis = if doc.features.iter().any(|f| f.id == piece.feature_id) {
            Axis::Top
        } else {
            Axis::Duplex
        };
        let Some(f) = feature(doc, axis, piece.feature_id) else {
            return false;
        };
        let Ok(full) = full_feature(doc, axis, f) else {
            return false;
        };
        let Some(extracted) = extract_current(doc, &full) else {
            return false;
        };
        piece.axis = axis;
        piece.mapping_digest = mapping_digest(piece, doc);
        piece.baseline = f.clone();
        piece.references_digest = references_digest(f, &doc.references);
        piece.sequence_digest = digest(extracted.as_bytes());
        true
    });
    registry.prune();
    Ok(())
}
/// Obtain actionable evidence, seeding original occurrences only where no prior
/// association exists. An invalidated piece is never silently re-seeded.
pub(crate) fn tracked(doc: &Document) -> Result<Registry, Error> {
    doc.annotation_lineage.validate(doc)?;
    let mut registry = doc.annotation_lineage.clone();
    registry
        .pieces
        .retain(|p| doc.annotation_lineage.actionable_piece(doc, p));
    if registry.pieces.len() != doc.annotation_lineage.pieces.len() {
        registry.tracking_limited = true;
    }
    registry.prune();
    let (extent, top, bottom) = geometry(doc);
    let mut source_digest = None::<String>;
    let mut original_bases = registry
        .origins
        .iter()
        .map(|o| o.feature.len())
        .sum::<usize>();
    let mut mapped_bases = registry
        .pieces
        .iter()
        .flat_map(|p| &p.runs)
        .map(|r| r.length)
        .sum::<usize>();
    let mut run_count = registry.pieces.iter().map(|p| p.runs.len()).sum::<usize>();
    let mut bytes = registry.snapshot_bytes();
    let mut work = registry
        .pieces
        .iter()
        .map(|p| {
            registry
                .origins
                .iter()
                .find(|o| o.id == p.origin_id)
                .map_or(0, |o| o.feature.len())
        })
        .sum::<usize>();
    for (axis, f) in doc.features.iter().map(|f| (Axis::Top, f)).chain(
        doc.duplex_annotations
            .features
            .iter()
            .map(|f| (Axis::Duplex, f)),
    ) {
        if doc
            .annotation_lineage
            .pieces
            .iter()
            .any(|p| p.feature_id == f.id)
        {
            continue;
        }
        if registry.tracking_limited {
            break;
        }
        let Some(ranges) = leaves(&location(f)) else {
            continue;
        };
        let bases = ranges.iter().map(|r| r.end - r.start).sum::<usize>();
        if original_bases.saturating_add(bases) > MAX_MAPPED_BASES
            || mapped_bases.saturating_add(bases) > MAX_MAPPED_BASES
            || work.saturating_add(bases) > MAX_VALIDATION_WORK
            || registry.origins.len() >= MAX_ORIGINS
            || registry.pieces.len() >= MAX_PIECES
        {
            registry.tracking_limited = true;
            break;
        }
        let full = full_feature(doc, axis, f)?;
        let ranges = leaves(&location(&full)).expect("same supported operators");
        let Some(extracted) = extract_current(doc, &full) else {
            continue;
        };
        let references = cited_references(f, &doc.references);
        let occurrence_id = hash(&(
            &registry.instance_id,
            f.id,
            registry
                .feature_epochs
                .iter()
                .find(|(id, _)| *id == f.id)
                .map(|(_, epoch)| epoch),
        ));
        let mut origin = Origin {
            id: String::new(),
            source_id: registry.instance_id.clone(),
            occurrence_id,
            top: top.clone(),
            bottom: bottom.clone(),
            feature: full,
            references,
            extent,
            molecule: doc.molecule,
            circular: doc.circular,
            source_digest: source_digest
                .get_or_insert_with(|| {
                    hash(&(
                        doc.sequence.as_ref(),
                        doc.duplex_ends.as_ref().map(|e| &e.opposite_sequence),
                        doc.molecule,
                        doc.circular,
                    ))
                })
                .clone(),
            sequence_digest: digest(extracted.as_bytes()),
        };
        if bytes.saturating_add(serialized_size(&origin)) > MAX_BYTES {
            registry.tracking_limited = true;
            break;
        }
        origin.id = origin_id(&origin);
        let id = origin.id.clone();
        let runs: Vec<_> = ranges
            .iter()
            .enumerate()
            .flat_map(|(leaf, r)| {
                split_strands(
                    Run {
                        leaf: leaf as u32,
                        offset: 0,
                        length: r.end - r.start,
                        target_start: r.start,
                        reverse: false,
                        strands: 3,
                    },
                    &top,
                    &bottom,
                )
            })
            .collect();
        let mut piece = Piece {
            mapping_digest: String::new(),
            origin_id: id,
            copy_id: registry.instance_id.clone(),
            feature_id: f.id,
            axis,
            runs,
            baseline: f.clone(),
            references_digest: references_digest(f, &doc.references),
            sequence_digest: digest(extracted.as_bytes()),
        };
        piece.mapping_digest = mapping_digest(&piece, doc);
        let additional = serialized_size(&origin)
            .saturating_add(serialized_size(&piece))
            .saturating_add(2);
        if bytes.saturating_add(additional) > MAX_BYTES
            || run_count.saturating_add(piece.runs.len()) > MAX_RUNS
        {
            registry.tracking_limited = true;
            break;
        }
        bytes += additional;
        original_bases += bases;
        work += bases;
        mapped_bases += bases;
        run_count += piece.runs.len();
        registry.origins.push(origin);
        registry.pieces.push(piece);
    }
    refresh(&mut registry, doc)?;
    registry.prune();
    registry.validate(doc)?;
    Ok(registry)
}
fn split_strands(run: Run, top: &Segment, bottom: &Segment) -> Vec<Run> {
    let end = run.target_start + run.length;
    let mut breaks = vec![run.target_start, end];
    for p in [top.start, top.end, bottom.start, bottom.end] {
        if p > run.target_start && p < end {
            breaks.push(p);
        }
    }
    breaks.sort_unstable();
    breaks.dedup();
    breaks
        .windows(2)
        .filter_map(|p| {
            let bits = u8::from((top.start..top.end).contains(&p[0]))
                | (u8::from((bottom.start..bottom.end).contains(&p[0])) << 1);
            let strands = run.strands & bits;
            if strands == 0 {
                return None;
            }
            Some(Run {
                offset: run.offset
                    + if run.reverse {
                        end - p[1]
                    } else {
                        p[0] - run.target_start
                    },
                length: p[1] - p[0],
                target_start: p[0],
                strands,
                ..run.clone()
            })
        })
        .collect()
}
fn map_runs(runs: &[Run], maps: &[Map]) -> Vec<Run> {
    let mut result = vec![];
    for run in runs {
        for map in maps {
            let a = run.target_start.max(map.source.start);
            let b = (run.target_start + run.length).min(map.source.end);
            if a >= b {
                continue;
            }
            result.push(Run {
                leaf: run.leaf,
                offset: run.offset
                    + if run.reverse {
                        run.target_start + run.length - b
                    } else {
                        a - run.target_start
                    },
                length: b - a,
                target_start: map.target_start
                    + if map.reverse {
                        map.source.end - b
                    } else {
                        a - map.source.start
                    },
                reverse: run.reverse ^ map.reverse,
                strands: if map.reverse {
                    ((run.strands & 1) << 1) | ((run.strands & 2) >> 1)
                } else {
                    run.strands
                },
            });
        }
    }
    result
}
pub(crate) fn limited_warning(doc: &mut Document) {
    const MESSAGE: &str = "Some annotation ancestry is unavailable because tracking evidence changed or exceeded bounded resources. Current annotations and physical sequence were retained; missing ancestry will not be guessed.";
    if doc.annotation_lineage.tracking_limited && !doc.warnings.iter().any(|w| w == MESSAGE) {
        doc.warnings.push(MESSAGE.into());
    }
}
fn bounded_finish(doc: &mut Document) -> Result<(), Error> {
    if let Err(error) = doc.annotation_lineage.validate(doc) {
        let text = error.to_string();
        if text.contains("limit") || text.contains("bounded storage") {
            doc.annotation_lineage.origins.clear();
            doc.annotation_lineage.pieces.clear();
            doc.annotation_lineage.tracking_limited = true;
        } else {
            return Err(error);
        }
    }
    limited_warning(doc);
    Ok(())
}
pub(crate) fn project(source: &Document, target: &mut Document, maps: &[Map]) -> Result<(), Error> {
    let mut registry = tracked(source)?;
    registry.instance_id = target.annotation_lineage.instance_id.clone();
    registry.feature_epochs = target.annotation_lineage.feature_epochs.clone();
    for piece in &mut registry.pieces {
        piece.runs = map_runs(&piece.runs, maps);
    }
    registry.pieces.retain(|p| !p.runs.is_empty());
    refresh(&mut registry, target)?;
    target.annotation_lineage = registry;
    bounded_finish(target)
}
pub(crate) fn to_reference(source: &Document, target: &mut Document) -> Result<(), Error> {
    let mut registry = tracked(source)?;
    refresh(&mut registry, target)?;
    target.annotation_lineage = registry;
    bounded_finish(target)
}
pub(crate) fn finish(doc: &mut Document, clip_strands: bool) -> Result<(), Error> {
    let mut registry = std::mem::take(&mut doc.annotation_lineage);
    if clip_strands {
        let (_, top, bottom) = geometry(doc);
        for p in &mut registry.pieces {
            p.runs = p
                .runs
                .iter()
                .flat_map(|r| split_strands(r.clone(), &top, &bottom))
                .collect();
        }
    }
    refresh(&mut registry, doc)?;
    doc.annotation_lineage = registry;
    bounded_finish(doc)
}
pub(crate) fn append(
    source: &Document,
    target: &mut Document,
    offset: usize,
    bindings: &[(u64, u64)],
) -> Result<(), Error> {
    if target.annotation_lineage.tracking_limited && target.annotation_lineage.pieces.is_empty() {
        limited_warning(target);
        return Ok(());
    }
    let registry = tracked(source)?;
    target.annotation_lineage.tracking_limited |= registry.tracking_limited;
    for origin in registry.origins {
        if let Some(old) = target
            .annotation_lineage
            .origins
            .iter()
            .find(|o| o.id == origin.id)
        {
            if old != &origin {
                return Err(error("conflicting original descriptors"));
            }
        } else {
            target.annotation_lineage.origins.push(origin);
        }
    }
    for mut piece in registry.pieces {
        if let Some((_, id)) = bindings.iter().find(|(id, _)| *id == piece.feature_id) {
            piece.feature_id = *id;
            piece.axis = Axis::Top;
            for r in &mut piece.runs {
                r.target_start = r
                    .target_start
                    .checked_add(offset)
                    .ok_or(Error::Coordinates)?;
            }
            target.annotation_lineage.pieces.push(piece);
        }
    }
    bounded_finish(target)
}
pub(crate) fn wrap(doc: &mut Document, origin: usize) {
    let n = doc.sequence.len();
    for piece in &mut doc.annotation_lineage.pieces {
        piece.runs = piece
            .runs
            .iter()
            .flat_map(|run| {
                let mut out = vec![];
                let mut used = 0;
                while used < run.length {
                    let start = (run.target_start + used + n - origin) % n;
                    let length = (run.length - used).min(n - start);
                    out.push(Run {
                        offset: run.offset
                            + if run.reverse {
                                run.length - used - length
                            } else {
                                used
                            },
                        target_start: start,
                        length,
                        ..run.clone()
                    });
                    used += length;
                }
                out
            })
            .collect();
    }
}
pub(crate) fn reverse(before: &Document, after: &mut Document) -> Result<(), Error> {
    let n = geometry(before).0;
    project(
        before,
        after,
        &[Map {
            source: Segment { start: 0, end: n },
            target_start: 0,
            reverse: true,
        }],
    )?;
    // Returning an original source to its exact original state needs only its
    // persisted seed; avoid introducing a dirty delta after two reversals.
    let current_digest = hash(&(
        after.sequence.as_ref(),
        after.duplex_ends.as_ref().map(|e| &e.opposite_sequence),
        after.molecule,
        after.circular,
    ));
    let mut registry = std::mem::take(&mut after.annotation_lineage);
    registry.pieces.retain(|p| {
        let Some(o) = registry.origins.iter().find(|o| o.id == p.origin_id) else {
            return false;
        };
        !(p.copy_id == registry.instance_id
            && o.source_digest == current_digest
            && feature(after, p.axis, p.feature_id)
                .is_some_and(|f| full_feature(after, p.axis, f).is_ok_and(|f| f == o.feature))
            && p.runs.iter().all(|r| {
                !r.reverse
                    && leaves(&location(&o.feature)).is_some_and(|ranges| {
                        r.target_start == ranges[r.leaf as usize].start + r.offset
                    })
            }))
    });
    registry.prune();
    after.annotation_lineage = registry;
    Ok(())
}
pub(crate) fn before_edit(
    doc: &Document,
    start: usize,
    end: usize,
    replacement: usize,
) -> Result<Registry, Error> {
    let mut registry = doc.annotation_lineage.clone();
    let (extent, top, _) = geometry(doc);
    let a = start + top.start;
    let b = end + top.start;
    registry.pieces.retain(|p| {
        registry_actionable(doc, p)
            && !p.runs.iter().any(|r| {
                if a == b {
                    r.target_start < a && a < r.target_start + r.length
                } else {
                    r.target_start < b && a < r.target_start + r.length
                }
            })
    });
    if registry.pieces.len() != doc.annotation_lineage.pieces.len() {
        registry.tracking_limited = true;
    }
    let maps = [
        Map {
            source: Segment { start: 0, end: a },
            target_start: 0,
            reverse: false,
        },
        Map {
            source: Segment {
                start: b,
                end: extent,
            },
            target_start: a + replacement,
            reverse: false,
        },
    ];
    for p in &mut registry.pieces {
        p.runs = map_runs(&p.runs, &maps);
    }
    registry.prune();
    Ok(registry)
}
fn registry_actionable(doc: &Document, p: &Piece) -> bool {
    doc.annotation_lineage.actionable_piece(doc, p)
}
pub(crate) fn new_copy(doc: &mut Document) {
    let id = uuid::Uuid::new_v4().to_string();
    for piece in &mut doc.annotation_lineage.pieces {
        piece.copy_id = id.clone();
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Cell {
    target: usize,
    reverse: bool,
    strands: u8,
}
type Coverage = Vec<Vec<Option<Cell>>>;
fn coverage(origin: &Origin, pieces: &[&Piece], mode: Merge) -> Option<Coverage> {
    let leaves = leaves(&location(&origin.feature))?;
    let mut cells = leaves
        .iter()
        .map(|p| vec![None::<Cell>; p.end - p.start])
        .collect::<Vec<_>>();
    for piece in pieces {
        for run in &piece.runs {
            for i in 0..run.length {
                let target = run.target_start + if run.reverse { run.length - 1 - i } else { i };
                let cell = cells.get_mut(run.leaf as usize)?.get_mut(run.offset + i)?;
                if let Some(old) = cell {
                    if old.target != target
                        || old.reverse != run.reverse
                        || mode == Merge::Physical && old.strands & run.strands != 0
                    {
                        return None;
                    }
                    old.strands |= run.strands;
                } else {
                    *cell = Some(Cell {
                        target,
                        reverse: run.reverse,
                        strands: run.strands,
                    });
                }
            }
        }
    }
    Some(cells)
}
fn toggle(value: L) -> L {
    if let L::Complement(child) = value {
        *child
    } else {
        L::Complement(Box::new(value))
    }
}
fn join(mut values: Vec<L>) -> L {
    if values.len() == 1 {
        values.remove(0)
    } else {
        L::Join(values)
    }
}
fn mark_ends(value: &mut L, first: bool, last: bool) {
    match value {
        L::Range { before, after, .. } => {
            *before |= first;
            *after |= last;
        }
        L::Complement(child) => mark_ends(child, last, first),
        L::Join(children) => {
            if let Some(c) = children.first_mut() {
                mark_ends(c, first, false);
            }
            if let Some(c) = children.last_mut() {
                mark_ends(c, false, last);
            }
        }
        _ => {}
    }
}
/// Each original Range must stay contiguous. Splitting is allowed only at an
/// actual circular product seam; a foreign insertion cannot become a new join.
fn mapped_location(original: &L, cells: &Coverage, n: usize, circular: bool) -> Option<L> {
    fn walk(
        original: &L,
        cells: &Coverage,
        leaf: &mut usize,
        n: usize,
        circular: bool,
    ) -> Option<Option<L>> {
        Some(match original {
            L::Range { before, after, .. } => {
                let row = cells.get(*leaf)?;
                *leaf += 1;
                let Some(first) = row.iter().position(Option::is_some) else {
                    return Some(None);
                };
                let last = row.iter().rposition(Option::is_some)?;
                if row[first..=last].iter().any(Option::is_none) {
                    return None;
                }
                let reverse = row[first]?.reverse;
                for pair in row[first..=last].windows(2) {
                    let (a, b) = (pair[0]?, pair[1]?);
                    if a.reverse != reverse || b.reverse != reverse {
                        return None;
                    }
                    let expected = if reverse {
                        a.target
                            .checked_sub(1)
                            .or_else(|| circular.then_some(n - 1))
                    } else {
                        a.target
                            .checked_add(1)
                            .filter(|v| *v < n)
                            .or_else(|| circular.then_some(0))
                    };
                    if expected != Some(b.target) {
                        return None;
                    }
                }
                let mut ranges = vec![];
                let mut a = first;
                while a <= last {
                    let mut b = a;
                    while b < last && row[b + 1]?.target.abs_diff(row[b]?.target) == 1 {
                        b += 1;
                    }
                    let low = row[a]?.target.min(row[b]?.target);
                    let high = row[a]?.target.max(row[b]?.target) + 1;
                    let begin = (*before || first > 0) && a == first;
                    let end = (*after || last + 1 < row.len()) && b == last;
                    let range = L::Range {
                        start: low,
                        end: high,
                        before: if reverse { end } else { begin },
                        after: if reverse { begin } else { end },
                    };
                    ranges.push(if reverse { toggle(range) } else { range });
                    a = b + 1;
                }
                Some(join(ranges))
            }
            L::Complement(child) => walk(child, cells, leaf, n, circular)?.map(toggle),
            L::Join(children) => {
                let mapped = children
                    .iter()
                    .map(|c| walk(c, cells, leaf, n, circular))
                    .collect::<Option<Vec<_>>>()?;
                let Some(first) = mapped.iter().position(Option::is_some) else {
                    return Some(None);
                };
                let last = mapped.iter().rposition(Option::is_some)?;
                // Dropped internal occurrences cannot be represented as an exact join.
                if mapped[first..=last].iter().any(Option::is_none) {
                    return None;
                }
                let mut kept = mapped.into_iter().flatten().collect::<Vec<_>>();
                if first > 0 {
                    mark_ends(kept.first_mut()?, true, false);
                }
                if last + 1 < children.len() {
                    mark_ends(kept.last_mut()?, false, true);
                }
                Some(L::Join(kept))
            }
            _ => return None,
        })
    }
    walk(original, cells, &mut 0, n, circular)?
}
fn biological_path(value: &L) -> Option<Vec<(usize, usize, bool)>> {
    let ranges = leaves(value)?;
    if ranges
        .iter()
        .try_fold(0usize, |n, r| n.checked_add(r.end - r.start))?
        > MAX_MAPPED_BASES
    {
        return None;
    }
    fn walk(value: &L, index: &mut usize) -> Option<Vec<(usize, usize, bool)>> {
        match value {
            L::Range { start, end, .. } => {
                let leaf = *index;
                *index += 1;
                Some((0..end - start).map(|i| (leaf, i, false)).collect())
            }
            L::Complement(child) => {
                let mut out = walk(child, index)?;
                out.reverse();
                for (_, _, strand) in &mut out {
                    *strand = !*strand;
                }
                Some(out)
            }
            L::Join(children) => {
                let mut out = vec![];
                for c in children {
                    out.extend(walk(c, index)?);
                }
                Some(out)
            }
            _ => None,
        }
    }
    walk(value, &mut 0)
}
fn biological_occurrences(value: &L) -> Option<Vec<(usize, usize)>> {
    Some(
        biological_path(value)?
            .into_iter()
            .map(|(l, o, _)| (l, o))
            .collect(),
    )
}

fn remap_citations(feature: &mut Feature, original: &[Reference], target: &mut Vec<Reference>) {
    for (key, value) in &mut feature.qualifiers {
        if key == "citation" {
            let reference = value
                .as_ref()
                .and_then(|v| {
                    v.trim()
                        .trim_start_matches('[')
                        .trim_end_matches(']')
                        .parse::<usize>()
                        .ok()
                })
                .and_then(|n| n.checked_sub(1))
                .and_then(|n| original.get(n));
            if let Some(reference) = reference {
                let index = target
                    .iter()
                    .position(|r| r == reference)
                    .unwrap_or_else(|| {
                        target.push(reference.clone());
                        target.len() - 1
                    });
                *value = Some(format!("[{}]", index + 1));
            } else {
                *key = "source_citation".into();
            }
        }
    }
}
fn reconstructed(origin: &Origin, cells: &Coverage, doc: &Document) -> Option<Feature> {
    let original = location(&origin.feature);
    let mapped = mapped_location(&original, cells, geometry(doc).0, doc.circular)?;
    let complete = cells.iter().all(|r| r.iter().all(Option::is_some));
    let mut feature = origin.feature.clone();
    let ranges = leaves(&original)?;
    let mut coords = BTreeMap::new();
    for (leaf, row) in cells.iter().enumerate() {
        for (offset, cell) in row.iter().enumerate() {
            if let Some(cell) = cell {
                let p = ranges[leaf].start + offset;
                if coords
                    .insert(p, cell.target)
                    .is_some_and(|old| old != cell.target)
                {
                    return None;
                }
            }
        }
    }
    let reversed = cells.iter().flatten().flatten().next()?.reverse;
    if cells
        .iter()
        .flatten()
        .flatten()
        .any(|c| c.reverse != reversed)
    {
        return None;
    }
    if reversed {
        feature.reverse_display_direction(origin.extent);
        feature.parts = feature
            .parts
            .iter()
            .rev()
            .map(|p| Segment {
                start: origin.extent - p.end,
                end: origin.extent - p.start,
            })
            .collect();
        feature.reverse = !feature.reverse;
    }
    let gaps = feature.display_gaps();
    let mut gap_maps = BTreeMap::new();
    if !gaps.is_empty() {
        if !complete {
            return None;
        }
        let n = geometry(doc).0;
        let mut delta = None;
        for (&source, &target) in &coords {
            let value = if reversed {
                target as isize + source as isize
            } else {
                target as isize - source as isize
            };
            let value = if doc.circular {
                value.rem_euclid(n as isize)
            } else {
                value
            };
            if delta.is_some_and(|old| old != value) {
                return None;
            }
            delta = Some(value);
        }
        let delta = delta?
            - if reversed {
                origin.extent as isize - 1
            } else {
                0
            };
        for gap in gaps {
            let start = gap.range.start as isize + delta;
            let length = gap.range.end - gap.range.start;
            let mapped = if doc.circular {
                let start = start.rem_euclid(n as isize) as usize;
                if length > n {
                    return None;
                }
                if start + length <= n {
                    vec![Segment {
                        start,
                        end: start + length,
                    }]
                } else {
                    vec![
                        Segment { start, end: n },
                        Segment {
                            start: 0,
                            end: start + length - n,
                        },
                    ]
                }
            } else {
                let start = usize::try_from(start).ok()?;
                if start.checked_add(length)? > n {
                    return None;
                }
                vec![Segment {
                    start,
                    end: start + length,
                }]
            };
            gap_maps.insert((gap.range.start, gap.range.end), mapped);
        }
    }
    feature.project_display_segments(|range| {
        if let Some(mapped) = gap_maps.get(&(range.start, range.end)) {
            return mapped.clone();
        }
        let mut out: Vec<Segment> = vec![];
        for p in range.start..range.end {
            let original_p = if reversed { origin.extent - 1 - p } else { p };
            if let Some(&target) = coords.get(&original_p) {
                if let Some(last) = out.last_mut() {
                    if last.end == target {
                        last.end += 1;
                        continue;
                    }
                    if target + 1 == last.start {
                        last.start = target;
                        continue;
                    }
                }
                out.push(Segment {
                    start: target,
                    end: target + 1,
                });
            }
        }
        out
    });
    (feature.parts, feature.reverse) = mapped
        .uniform_projection()
        .unwrap_or_else(|| (mapped.local_parts(), false));
    feature.location = Some(mapped);
    crate::translation::project_translation_exceptions(&mut feature, |p| coords.get(&p).copied());
    let extracted = extract_current(doc, &feature)?;
    if complete {
        if digest(extracted.as_bytes()) != origin.sequence_digest {
            return None;
        }
    } else {
        feature.edited = true;
        feature.qualifiers.retain(|(k, _)| k != "translation");
        if !feature.qualifiers.iter().any(|(k, _)| k == "partial") {
            feature.qualifiers.push(("partial".into(), None));
        }
        let order = biological_occurrences(&original)?;
        let first = order.iter().position(|&(l, o)| cells[l][o].is_some())?;
        let last = order.iter().rposition(|&(l, o)| cells[l][o].is_some())?;
        if order[first..=last]
            .iter()
            .any(|&(l, o)| cells[l][o].is_none())
        {
            return None;
        }
        if feature.kind.eq_ignore_ascii_case("CDS")
            && let Ok(phase) = crate::translation::cds_phase(&origin.feature)
        {
            let phase = (phase + 3 - first % 3) % 3 + 1;
            feature.qualifiers.retain(|(k, _)| k != "codon_start");
            feature
                .qualifiers
                .push(("codon_start".into(), Some(phase.to_string())));
            for (_, v) in feature
                .qualifiers
                .iter_mut()
                .filter(|(k, _)| k == "snapgene_readingFrame")
            {
                *v = Some(format!("{}{phase}", if feature.reverse { "-" } else { "" }));
            }
        }
    }
    // Native readingFrame sign follows the restored current strand, while phase
    // remains feature-relative. Do not repair malformed/conflicting source data.
    if complete && let Ok(phase) = crate::translation::cds_phase(&origin.feature) {
        for (_, v) in feature
            .qualifiers
            .iter_mut()
            .filter(|(k, _)| k == "snapgene_readingFrame")
        {
            *v = Some(format!(
                "{}{}",
                if feature.reverse { "-" } else { "" },
                phase + 1
            ));
        }
    }
    Some(feature)
}
/// Consolidate only proven occurrence paths after the physical construction has
/// succeeded. Input qualities, ends, methylation and sequence are never restored.
pub(crate) fn reconstruct(doc: &mut Document, mode: Merge) -> Result<(), Error> {
    finish(doc, false)?;
    doc.annotation_lineage.validate(doc)?;
    let registry = doc.annotation_lineage.clone();
    let mut groups: BTreeMap<(String, String), Vec<&Piece>> = BTreeMap::new();
    for p in &registry.pieces {
        if !registry.actionable_piece(doc, p) {
            continue;
        }
        groups
            .entry((
                p.origin_id.clone(),
                if mode == Merge::VerifiedOverlap {
                    String::new()
                } else {
                    p.copy_id.clone()
                },
            ))
            .or_default()
            .push(p);
    }
    let mut restored = 0;
    let mut partial = 0;
    let mut unresolved = 0;
    for ((origin_id, _), pieces) in groups {
        if pieces.len() < 2 {
            continue;
        }
        let origin = registry
            .origins
            .iter()
            .find(|o| o.id == origin_id)
            .expect("validated origin");
        let candidates = if let Some(cells) = coverage(origin, &pieces, mode) {
            vec![(pieces.clone(), cells)]
        } else {
            // Multiple complete copies can be separated only when each piece has
            // one consistent placement transform. Do not guess joins by labels.
            let ranges = leaves(&location(&origin.feature)).expect("validated leaves");
            let mut by_placement: BTreeMap<(bool, isize), Vec<&Piece>> = BTreeMap::new();
            let mut ambiguous = false;
            for p in &pieces {
                let mut key = None;
                for r in &p.runs {
                    let source = ranges[r.leaf as usize].start + r.offset;
                    let target = r.target_start + if r.reverse { r.length - 1 } else { 0 };
                    let delta = if r.reverse {
                        target as isize + source as isize
                    } else {
                        target as isize - source as isize
                    };
                    let delta = if doc.circular {
                        delta.rem_euclid(doc.sequence.len() as isize)
                    } else {
                        delta
                    };
                    if key.is_some_and(|old| old != (r.reverse, delta)) {
                        ambiguous = true;
                        break;
                    }
                    key = Some((r.reverse, delta));
                }
                if let Some(key) = key {
                    by_placement.entry(key).or_default().push(p);
                }
            }
            if ambiguous {
                vec![]
            } else {
                by_placement
                    .into_values()
                    .filter_map(|p| coverage(origin, &p, mode).map(|c| (p, c)))
                    .collect()
            }
        };
        let mut changed = false;
        for (component, cells) in candidates {
            if component.len() < 2 {
                continue;
            }
            let Some(mut feature) = reconstructed(origin, &cells, doc) else {
                continue;
            };
            let complete = cells.iter().all(|r| r.iter().all(Option::is_some));
            // Do not consolidate a partial component onto another physical copy.
            if !complete && component.len() != pieces.len() {
                continue;
            }
            let ids = component
                .iter()
                .map(|p| p.feature_id)
                .collect::<BTreeSet<_>>();
            feature.id = *ids.first().expect("nonempty component");
            remap_citations(&mut feature, &origin.references, &mut doc.references);
            let mut piece = component[0].clone();
            piece.feature_id = feature.id;
            piece.runs.clear();
            for (leaf, row) in cells.iter().enumerate() {
                for (offset, cell) in row.iter().enumerate() {
                    if let Some(c) = cell {
                        piece.runs.push(Run {
                            leaf: leaf as u32,
                            offset,
                            length: 1,
                            target_start: c.target,
                            reverse: c.reverse,
                            strands: c.strands,
                        });
                    }
                }
            }
            coalesce(&mut piece.runs);
            doc.features.retain(|f| !ids.contains(&f.id));
            doc.duplex_annotations
                .features
                .retain(|f| !ids.contains(&f.id));
            let top = geometry(doc).1;
            if feature
                .parts
                .iter()
                .all(|p| p.start >= top.start && p.end <= top.end)
            {
                crate::workflow::shift_feature(&mut feature, -(top.start as isize))?;
                piece.axis = Axis::Top;
                doc.features.push(feature);
            } else {
                piece.axis = Axis::Duplex;
                doc.duplex_annotations.features.push(feature);
            }
            doc.annotation_lineage
                .pieces
                .retain(|p| !ids.contains(&p.feature_id));
            doc.annotation_lineage.pieces.push(piece);
            if complete {
                restored += 1
            } else {
                partial += 1
            };
            changed = true;
        }
        if !changed {
            unresolved += 1;
        }
    }
    finish(doc, false)?;
    if mode == Merge::VerifiedOverlap {
        new_copy(doc);
    }
    if restored + partial > 0 {
        doc.warnings.push(format!("Annotation ancestry reconstructed {restored} complete feature(s) and consolidated {partial} partial feature(s); current product sequence and physical strands were retained."));
    }
    if unresolved > 0 {
        doc.warnings.push(format!("{unresolved} annotation ancestry group(s) remain separate: missing, conflicting, duplicated or noncontiguous occurrence paths prevent safe reconstruction."));
    }
    doc.annotation_lineage.validate(doc)
}
fn coalesce(runs: &mut Vec<Run>) {
    let mut out: Vec<Run> = vec![];
    for run in std::mem::take(runs) {
        if let Some(last) = out.last_mut()
            && last.leaf == run.leaf
            && last.offset + last.length == run.offset
            && last.reverse == run.reverse
            && last.strands == run.strands
            && (if run.reverse {
                run.target_start + run.length == last.target_start
            } else {
                last.target_start + last.length == run.target_start
            })
        {
            if run.reverse {
                last.target_start = run.target_start;
            }
            last.length += run.length;
        } else {
            out.push(run);
        }
    }
    *runs = out;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        cloning::{self, AssemblyInput, LigationOptions},
        restriction::{self, DigestConditions, ENZYMES},
        translation, workflow,
    };
    fn source(sequence: &str, phase: usize) -> Document {
        let mut d = Document::new("source", sequence).unwrap();
        d.annotate(
            "same label".into(),
            Segment {
                start: 0,
                end: sequence.len(),
            },
            false,
        )
        .unwrap();
        d.features[0].kind = "CDS".into();
        d.features[0]
            .qualifiers
            .push(("codon_start".into(), Some(phase.to_string())));
        d
    }
    fn digest_products(d: &Document) -> Vec<Document> {
        let i = ENZYMES.iter().position(|e| e.name == "EcoRI").unwrap();
        let c = DigestConditions::default();
        let result = restriction::digest_document(d, &[i], c).unwrap();
        workflow::digest_products(d, &result.fragments, &["EcoRI".into()], c).unwrap()
    }
    fn ligate(docs: &[Document], reverse: bool, circular: bool) -> Document {
        cloning::ligate(
            &docs
                .iter()
                .map(|d| AssemblyInput {
                    document: d.clone(),
                    reverse,
                })
                .collect::<Vec<_>>(),
            LigationOptions {
                circular,
                phosphorylate: true,
                allow_nicks: false,
            },
            "product",
        )
        .unwrap()
    }
    fn aa(d: &Document, f: &Feature) -> String {
        String::from_utf8(translation::translate_cds(d, f).unwrap().amino_acids).unwrap()
    }
    #[test]
    fn ecori_complete_reassembly_restores_original_cds_and_all_phases() {
        for (phase, expected) in [(1, "MEFK*"), (2, "WNSS"), (3, "GIQV")] {
            let d = source("ATGGAATTCAAGTAA", phase);
            let fragments = digest_products(&d);
            assert_eq!(fragments[0].sequence.as_ref(), "ATGG");
            assert_eq!(
                fragments[0].duplex_ends.as_ref().unwrap().opposite_sequence,
                "AATTCCAT"
            );
            let product = ligate(&fragments, false, false);
            assert_eq!(product.sequence, d.sequence);
            assert_eq!(product.features.len(), 1, "{:#?}", product.features);
            assert_eq!(product.features[0].location_label(), "1..15");
            assert_eq!(aa(&product, &product.features[0]), expected);
            assert!(!product.features[0].edited);
            assert_eq!(product.features[0].qualifier("partial"), None);
        }
    }
    #[test]
    fn reassembly_handles_reverse_and_circular_cohesive_occurrences() {
        let d = source("ATGGAATTCAAGTAA", 1);
        let mut fragments = digest_products(&d);
        fragments.reverse();
        let product = ligate(&fragments, true, false);
        assert_eq!(product.sequence.as_ref(), "TTACTTGAATTCCAT");
        assert_eq!(product.features.len(), 1);
        assert_eq!(
            product.feature_sequence(&product.features[0]).unwrap(),
            d.sequence.as_ref()
        );
        assert_eq!(aa(&product, &product.features[0]), "MEFK*");
        let mut d = d;
        d.circular = true;
        let fragments = digest_products(&d);
        let product = ligate(&fragments, false, true);
        assert_eq!(product.sequence.as_ref(), "AATTCAAGTAAATGG");
        assert_eq!(product.features.len(), 1, "{:#?}", product.features);
        assert_eq!(product.features[0].location_label(), "join(12..15,1..11)");
        assert_eq!(aa(&product, &product.features[0]), "MEFK*");
    }
    #[test]
    fn staged_partial_reassembly_retains_phase_then_completes() {
        let d = source("ATGGAATTCAAGGAATTCTAA", 1);
        let pieces = digest_products(&d);
        assert_eq!(pieces.len(), 3);
        let partial = ligate(&pieces[..2], false, false);
        let reference = workflow::duplex_reference(&partial).unwrap();
        assert_eq!(reference.features.len(), 1, "{:#?}", reference.features);
        assert_eq!(reference.features[0].location_label(), "1..>17");
        assert_eq!(aa(&reference, &reference.features[0]), "MEFKE");
        assert!(reference.features[0].edited);
        assert_eq!(reference.features[0].qualifier("translation"), None);
        let complete = ligate(&[partial, pieces[2].clone()], false, false);
        assert_eq!(complete.features.len(), 1);
        assert_eq!(aa(&complete, &complete.features[0]), "MEFKEF*");
    }
    #[test]
    fn mixed_source_feature_restores_without_joining_foreign_interruption() {
        let a = source("ATGGAATTCAAGTAA", 1);
        let a = digest_products(&a);
        let b = source("CCGTTT", 1);
        let left = workflow::extract_region(&b, 0, 3, "left").unwrap();
        let right = workflow::extract_region(&b, 3, 3, "right").unwrap();
        let product = ligate(&[left, a[0].clone(), a[1].clone(), right], false, false);
        assert_eq!(product.sequence.as_ref(), "CCGATGGAATTCAAGTAATTT");
        assert_eq!(product.features.len(), 3);
        let restored = product
            .features
            .iter()
            .find(|f| f.parts == vec![Segment { start: 3, end: 18 }])
            .unwrap();
        assert_eq!(aa(&product, restored), "MEFK*");
        assert_eq!(product.features.iter().filter(|f| f.edited).count(), 2);
    }
    #[test]
    fn reverse_restoration_maps_exceptions_arrows_and_interbase_sites() {
        let mut d = source("ATGGAATTCTGATAA", 1);
        d.features[0].qualifiers.extend([
            ("transl_except".into(), Some("(pos:10..12,aa:Sec)".into())),
            ("snapgene_directionality".into(), Some("1".into())),
            ("note".into(), Some("This forward directional feature has 2 segments: 1: 1 .. 8 / #ff0000 / first 2: 9 .. 15 / #0000ff / second Cleavage site after base 4".into())),
        ]);
        let mut fragments = digest_products(&d);
        fragments.reverse();
        let product = ligate(&fragments, true, false);
        assert_eq!(product.sequence.as_ref(), "TTATCAGAATTCCAT");
        assert_eq!(product.features.len(), 1);
        let f = &product.features[0];
        assert_eq!(f.qualifier("snapgene_directionality"), Some("2"));
        assert_eq!(f.display_cleavage_sites(), vec![11]);
        let segments = f.display_segments();
        assert_eq!(segments.len(), 2);
        assert_eq!(segments[0].range, Segment { start: 0, end: 7 });
        assert_eq!(segments[0].name.as_deref(), Some("second"));
        assert_eq!(segments[0].color_rgb, Some(0x0000ff));
        assert_eq!(segments[1].range, Segment { start: 7, end: 15 });
        assert_eq!(segments[1].name.as_deref(), Some("first"));
        assert_eq!(
            f.qualifier("transl_except"),
            Some("(pos:complement(4..6),aa:Sec)")
        );
        assert_eq!(aa(&product, f), "MEFU*");
    }
    #[test]
    fn duplicate_join_leaves_remain_distinct_and_genuine_partialness_survives() {
        let mut d = source("ATGGAATTCAAGTAA", 2);
        let ast = L::Join(vec![
            L::Range {
                start: 0,
                end: 8,
                before: true,
                after: false,
            },
            L::Range {
                start: 0,
                end: 8,
                before: false,
                after: false,
            },
            L::Range {
                start: 4,
                end: 12,
                before: false,
                after: true,
            },
        ]);
        d.features[0].parts = ast.local_parts();
        d.features[0].location = Some(ast.clone());
        let product = ligate(&digest_products(&d), false, false);
        assert_eq!(product.features.len(), 1, "{:#?}", product.features);
        assert_eq!(
            product.feature_sequence(&product.features[0]).unwrap(),
            "ATGGAATTATGGAATTAATTCAAG"
        );
        assert_eq!(product.features[0].location, Some(ast));
        assert_eq!(product.features[0].qualifier("codon_start"), Some("2"));
        assert!(!product.features[0].edited);
    }
    #[test]
    fn references_and_same_label_sources_do_not_conflate() {
        let mut a = source("ATGGAATTCAAGTAA", 1);
        a.references.push(Reference {
            title: "A citation".into(),
            ..Default::default()
        });
        a.features[0]
            .qualifiers
            .push(("citation".into(), Some("[1]".into())));
        let mut b = source("CCGTTT", 1);
        b.references.push(Reference {
            title: "B citation".into(),
            ..Default::default()
        });
        b.features[0]
            .qualifiers
            .push(("citation".into(), Some("[1]".into())));
        let fragments = digest_products(&a);
        let product = ligate(
            &[b, fragments[0].clone(), fragments[1].clone()],
            false,
            false,
        );
        assert_eq!(
            product
                .references
                .iter()
                .map(|r| r.title.as_str())
                .collect::<Vec<_>>(),
            vec!["B citation", "A citation"]
        );
        let f = product
            .features
            .iter()
            .find(|f| f.parts == vec![Segment { start: 6, end: 21 }])
            .unwrap();
        assert_eq!(f.qualifier("citation"), Some("[2]"));
        let mut changed = a.clone();
        changed.references[0].title = "Changed citation".into();
        let changed = digest_products(&changed);
        let product = ligate(&[fragments[0].clone(), changed[1].clone()], false, false);
        assert_eq!(
            product.features.len() + product.duplex_annotations.features.len(),
            2
        );
        let mut unused = a.clone();
        unused.references.push(Reference {
            title: "Unrelated".into(),
            ..Default::default()
        });
        let unused = digest_products(&unused);
        let product = ligate(&[fragments[0].clone(), unused[1].clone()], false, false);
        assert_eq!(product.features.len(), 1);
    }
    #[test]
    fn edited_corrupt_and_independent_pieces_fail_closed() {
        let a = source("ATGGAATTCAAGTAA", 1);
        let b = source("ATGGAATTCAAGTAA", 1);
        let af = digest_products(&a);
        let bf = digest_products(&b);
        assert_eq!(
            ligate(&[af[0].clone(), bf[1].clone()], false, false)
                .features
                .len(),
            2
        );
        for mutation in 0..4 {
            let mut altered = af[1].clone();
            match mutation {
                0 => altered.features[0].label = "edited".into(),
                1 => altered.sequence = "AATTCAAGCAA".into(),
                2 => altered.annotation_lineage.pieces[0].runs[0].target_start += 1,
                _ => altered.annotation_lineage.pieces[0].runs[0].strands = 2,
            };
            let p = &altered.annotation_lineage.pieces[0];
            assert!(
                !altered.annotation_lineage.actionable_piece(&altered, p),
                "mutation {mutation}"
            );
            if mutation != 1 {
                let product = ligate(&[af[0].clone(), altered], false, false);
                assert_eq!(product.features.len(), 2, "mutation {mutation}");
            }
        }
    }
    #[test]
    fn session_seed_redigestion_undo_and_occurrence_limits() {
        let d = source("ATGGAATTCAAGTAA", 1);
        let fragments = digest_products(&d);
        let restored: Document = serde_json::from_slice(&serde_json::to_vec(&d).unwrap()).unwrap();
        let again = digest_products(&restored);
        assert_eq!(
            ligate(&[fragments[0].clone(), again[1].clone()], false, false)
                .features
                .len(),
            1
        );
        let mut decoded: Document =
            serde_json::from_slice(&serde_json::to_vec(&fragments[1]).unwrap()).unwrap();
        let previous = decoded.annotation_lineage.clone();
        decoded.replace(5, 6, "T").unwrap();
        assert!(decoded.annotation_lineage.pieces.is_empty());
        decoded.annotation_lineage = previous;
        assert!(
            !decoded
                .annotation_lineage
                .actionable_piece(&decoded, &decoded.annotation_lineage.pieces[0])
        );
        let mut bad = fragments[1].clone();
        bad.annotation_lineage.origins[0].feature.parts = vec![Segment {
            start: 0,
            end: MAX_MAPPED_BASES + 1,
        }];
        bad.annotation_lineage.origins[0].feature.location = None;
        bad.annotation_lineage.origins[0].extent = MAX_MAPPED_BASES + 1;
        bad.annotation_lineage.origins[0].id = origin_id(&bad.annotation_lineage.origins[0]);
        bad.annotation_lineage.pieces[0].origin_id = bad.annotation_lineage.origins[0].id.clone();
        assert!(
            bad.validate()
                .unwrap_err()
                .to_string()
                .contains("occurrence work limit")
        );
        let mut legacy = serde_json::to_value(&d).unwrap();
        legacy.as_object_mut().unwrap().remove("annotation_lineage");
        let old: Document = serde_json::from_value(legacy).unwrap();
        assert!(old.annotation_lineage.is_empty());
        old.validate().unwrap();
    }

    #[test]
    fn repeated_dna_cannot_hide_swapped_original_occurrences() {
        let mut d = source("AAAAAA", 1);
        d.features[0].parts = vec![Segment { start: 0, end: 3 }, Segment { start: 3, end: 6 }];
        d.features[0].location = Some(L::Join(
            d.features[0]
                .parts
                .iter()
                .map(|p| L::Range {
                    start: p.start,
                    end: p.end,
                    before: false,
                    after: false,
                })
                .collect(),
        ));
        d.annotation_lineage = tracked(&d).unwrap();
        let p = &d.annotation_lineage.pieces[0];
        assert!(d.annotation_lineage.actionable_piece(&d, p));
        d.annotation_lineage.pieces[0].runs[0].target_start = 3;
        d.annotation_lineage.pieces[0].runs[1].target_start = 0;
        let checksum = mapping_digest(&d.annotation_lineage.pieces[0], &d);
        d.annotation_lineage.pieces[0].mapping_digest = checksum;
        assert!(
            !d.annotation_lineage
                .actionable_piece(&d, &d.annotation_lineage.pieces[0])
        );
    }

    #[test]
    fn two_copies_missing_internal_leaves_and_linear_reorder_are_distinct() {
        let d = source("ATGGAATTCAAGTAA", 1);
        let f = digest_products(&d);
        let product = ligate(
            &[f[0].clone(), f[1].clone(), f[0].clone(), f[1].clone()],
            false,
            false,
        );
        assert_eq!(product.features.len(), 2);
        assert_eq!(
            product
                .features
                .iter()
                .map(|f| f.parts.clone())
                .collect::<Vec<_>>(),
            vec![
                vec![Segment { start: 0, end: 15 }],
                vec![Segment { start: 15, end: 30 }]
            ]
        );
        let reordered = ligate(&[f[1].clone(), f[0].clone()], false, false);
        assert_eq!(
            reordered.features.len() + reordered.duplex_annotations.features.len(),
            2
        );
        let mut joined = d.clone();
        let ast = L::Join(vec![
            L::Range {
                start: 0,
                end: 3,
                before: false,
                after: false,
            },
            L::Range {
                start: 6,
                end: 9,
                before: false,
                after: false,
            },
            L::Range {
                start: 12,
                end: 15,
                before: false,
                after: false,
            },
        ]);
        joined.features[0].parts = ast.local_parts();
        joined.features[0].location = Some(ast);
        let first = workflow::extract_region(&joined, 0, 3, "first").unwrap();
        let last = workflow::extract_region(&joined, 12, 3, "last").unwrap();
        let product = ligate(&[first, last], false, false);
        assert_eq!(product.sequence.as_ref(), "ATGTAA");
        assert_eq!(product.features.len(), 2);
        assert!(product.features.iter().all(|f| f.edited));
    }
    #[test]
    fn verified_overlap_reconverges_synthesized_copies_but_not_unrelated_origins() {
        let mut d = source("ATGCCGTAACGTTGACCTAGGCTAACGT", 1);
        d.features[0].kind = "misc_feature".into();
        d.sequence_qualities = Some(vec![Some(30); 28]);
        d.methylation = Some(Default::default());
        let mut a = workflow::extract_region(&d, 0, 20, "a").unwrap();
        let mut b = workflow::extract_region(&d, 10, 18, "b").unwrap();
        new_copy(&mut a);
        new_copy(&mut b);
        assert_ne!(
            a.annotation_lineage.pieces[0].copy_id,
            b.annotation_lineage.pieces[0].copy_id
        );
        let product = cloning::assemble_overlaps(
            &[
                AssemblyInput {
                    document: a.clone(),
                    reverse: false,
                },
                AssemblyInput {
                    document: b,
                    reverse: false,
                },
            ],
            10,
            10,
            false,
            "assembled",
        )
        .unwrap();
        assert_eq!(product.sequence, d.sequence);
        assert_eq!(product.features.len(), 1);
        assert!(
            product.sequence_qualities.is_none()
                && product.opposite_sequence_qualities.is_none()
                && product.methylation.is_none()
        );
        assert!(
            product
                .annotation_lineage
                .actionable_piece(&product, &product.annotation_lineage.pieces[0])
        );
        let mut independent = source(&d.sequence, 1);
        independent.features[0].kind = "misc_feature".into();
        let b = workflow::extract_region(&independent, 10, 18, "unrelated").unwrap();
        let product = cloning::assemble_overlaps(
            &[
                AssemblyInput {
                    document: a,
                    reverse: false,
                },
                AssemblyInput {
                    document: b,
                    reverse: false,
                },
            ],
            10,
            10,
            false,
            "separate",
        )
        .unwrap();
        assert_eq!(product.features.len(), 2);
    }
    #[test]
    fn asymmetric_physical_strands_rejoin_full_extent_annotation() {
        use crate::restriction::{EndKind, FragmentEnd};
        let axis = "TTAACCGGAATTCTGGAACCGGTT";
        for left in [EndKind::FivePrime, EndKind::ThreePrime] {
            for right in [EndKind::FivePrime, EndKind::ThreePrime] {
                let top_start = usize::from(left == EndKind::ThreePrime) * 2;
                let top_end = axis.len() - usize::from(right == EndKind::FivePrime) * 3;
                let bottom_start = usize::from(left == EndKind::FivePrime) * 2;
                let bottom_end = axis.len() - usize::from(right == EndKind::ThreePrime) * 3;
                let mut d = source(axis, 1);
                d.features[0].kind = "misc_feature".into();
                d.features[0].reverse = true;
                d.sequence = axis[top_start..top_end].into();
                d.duplex_ends = Some(workflow::DuplexEnds {
                    left: FragmentEnd {
                        kind: left,
                        overhang: if left == EndKind::FivePrime {
                            axis[..2].into()
                        } else {
                            crate::reverse_complement(&axis[..2])
                        },
                        enzymes: vec![],
                    },
                    right: FragmentEnd {
                        kind: right,
                        overhang: if right == EndKind::ThreePrime {
                            axis[axis.len() - 3..].into()
                        } else {
                            crate::reverse_complement(&axis[axis.len() - 3..])
                        },
                        enzymes: vec![],
                    },
                    opposite_sequence: crate::reverse_complement(&axis[bottom_start..bottom_end]),
                    left_modification: Some("5′ phosphate".into()),
                    right_modification: Some("5′ phosphate".into()),
                });
                workflow::partition_duplex_annotations(&mut d).unwrap();
                let product = ligate(&digest_products(&d), false, false);
                assert_eq!(product.sequence, d.sequence);
                assert_eq!(product.duplex_ends, d.duplex_ends);
                let reference = workflow::duplex_reference(&product).unwrap();
                assert_eq!(reference.features.len(), 1);
                assert_eq!(
                    reference.feature_sequence(&reference.features[0]).unwrap(),
                    crate::reverse_complement(axis)
                );
                assert!(
                    product
                        .annotation_lineage
                        .pieces
                        .iter()
                        .all(|p| product.annotation_lineage.actionable_piece(&product, p))
                );
            }
        }
    }
    #[test]
    fn pcr_tail_shifts_and_redigestion_retain_valid_copy_evidence() {
        let sequence =
            include_str!("../../../validation/fixtures/primer-design-template.txt").trim();
        let mut d = source(sequence, 1);
        d.features[0].kind = "misc_feature".into();
        fn primer(id: u64, sequence: String, start: usize, reverse: bool) -> Feature {
            Feature {
                id,
                label: format!("primer{id}"),
                kind: "primer_bind".into(),
                parts: vec![Segment {
                    start,
                    end: start + 20,
                }],
                location: None,
                reverse,
                qualifiers: vec![("primer_sequence".into(), Some(sequence))],
                edited: false,
            }
        }
        let f = primer(100, format!("GGG{}", &sequence[139..159]), 139, false);
        let r = primer(
            101,
            crate::reverse_complement(&sequence[394..414]),
            394,
            true,
        );
        let product = cloning::pcr_product(&d, &f, &r, Default::default(), "copy").unwrap();
        assert_eq!(
            product.sequence.as_ref(),
            format!("GGG{}", &sequence[139..414])
        );
        let inherited = product
            .annotation_lineage
            .pieces
            .iter()
            .find(|p| p.baseline.kind == "misc_feature")
            .unwrap();
        assert!(
            product
                .annotation_lineage
                .actionable_piece(&product, inherited)
        );
        assert_ne!(inherited.copy_id, d.annotation_lineage.instance_id);
        assert_eq!(inherited.runs[0].target_start, 3);
        let recovered: Document =
            serde_json::from_slice(&serde_json::to_vec(&product).unwrap()).unwrap();
        let redigested = digest_products(&recovered);
        assert!(!redigested.is_empty());
        assert!(redigested.iter().all(|d| {
            d.annotation_lineage
                .pieces
                .iter()
                .all(|p| d.annotation_lineage.actionable_piece(d, p))
        }));
    }
    #[test]
    fn large_optional_ancestry_falls_back_without_lowering_sequence_limits() {
        let mut d = source(&"A".repeat(3_000_000), 1);
        d.features[0].kind = "misc_feature".into();
        d.references.push(Reference {
            title: "large source".into(),
            ..Default::default()
        });
        d.features[0]
            .qualifiers
            .push(("citation".into(), Some("[1]".into())));
        let product = workflow::extract_region(&d, 100, 100, "large cut").unwrap();
        assert_eq!(product.sequence.as_ref(), "A".repeat(100));
        assert_eq!(product.features.len(), 1);
        assert_eq!(product.references, d.references);
        assert!(product.features[0].edited);
        assert!(product.annotation_lineage.tracking_limited);
        assert!(
            product.annotation_lineage.origins.is_empty()
                && product.annotation_lineage.pieces.is_empty()
        );
        assert!(
            product
                .warnings
                .iter()
                .any(|w| w.contains("bounded resources"))
        );
        let again = workflow::extract_region(&product, 0, 50, "second cut").unwrap();
        assert!(again.annotation_lineage.tracking_limited);
        assert!(again.annotation_lineage.origins.is_empty());
        assert_eq!(again.features.len(), 1);
        d.validate().unwrap();
    }
    #[test]
    fn integrity_digests_detect_origin_mask_and_geometry_changes() {
        let d = source("ATGGAATTCAAGTAA", 1);
        let pieces = digest_products(&d);
        let mut altered = pieces[1].clone();
        altered.annotation_lineage.pieces[0]
            .runs
            .last_mut()
            .unwrap()
            .strands = 1;
        assert!(
            !altered
                .annotation_lineage
                .actionable_piece(&altered, &altered.annotation_lineage.pieces[0])
        );
        let mut origin = pieces[1].clone();
        origin.annotation_lineage.origins[0].feature.label = "corrupt original label".into();
        assert!(origin.validate().is_err());
        let mut removed = d.clone();
        let first = digest_products(&removed);
        removed.features.clear();
        removed
            .annotate("same label".into(), Segment { start: 0, end: 15 }, false)
            .unwrap();
        removed.features[0].kind = "CDS".into();
        removed.features[0]
            .qualifiers
            .push(("codon_start".into(), Some("1".into())));
        let second = digest_products(&removed);
        assert_ne!(
            first[0].annotation_lineage.origins[0].id,
            second[0].annotation_lineage.origins[0].id
        );
    }
    #[test]
    fn exhausted_epoch_budget_keeps_new_annotation_available() {
        let mut d = source("ACGT", 1);
        d.features[0].id = MAX_PIECES as u64;
        let epoch = uuid::Uuid::new_v4().to_string();
        d.annotation_lineage.feature_epochs = (0..MAX_PIECES as u64)
            .map(|id| (id, epoch.clone()))
            .collect();
        d.annotate("extra".into(), Segment { start: 1, end: 3 }, false)
            .unwrap();
        d.validate().unwrap();
        assert_eq!(d.features.len(), 2);
        assert_eq!(d.features[1].label, "extra");
        assert!(d.annotation_lineage.tracking_limited);
        assert_eq!(d.annotation_lineage.feature_epochs.len(), MAX_PIECES);
        assert!(d.warnings.iter().any(|w| w.contains("bounded resources")));
    }
    #[test]
    fn exact_join_reassembly_preserves_nonbiological_gap_connectors() {
        let mut d = source("ATGGAATTCAAGTAA", 1);
        d.features[0].kind = "misc_feature".into();
        let ast = L::Join(vec![
            L::Range {
                start: 0,
                end: 6,
                before: false,
                after: false,
            },
            L::Range {
                start: 10,
                end: 15,
                before: false,
                after: false,
            },
        ]);
        d.features[0].parts = ast.local_parts();
        d.features[0].location = Some(ast);
        d.features[0].qualifiers.push((
            "snapgene_gap".into(),
            Some("7 .. 10 / #aabbcc / bridge".into()),
        ));
        let fragments = digest_products(&d);
        let product = ligate(&fragments, false, false);
        assert_eq!(product.features.len(), 1);
        let gaps = product.features[0].display_gaps();
        assert_eq!(gaps.len(), 1);
        assert_eq!(gaps[0].range, Segment { start: 6, end: 10 });
        assert_eq!(gaps[0].name.as_deref(), Some("bridge"));
        assert_eq!(gaps[0].color_rgb, Some(0xaabbcc));
        let mut fragments = fragments;
        fragments.reverse();
        let reversed = ligate(&fragments, true, false);
        assert_eq!(reversed.features.len(), 1);
        assert_eq!(
            reversed.features[0].display_gaps()[0].range,
            Segment { start: 5, end: 9 }
        );
        assert_eq!(
            reversed.feature_sequence(&reversed.features[0]).unwrap(),
            d.feature_sequence(&d.features[0]).unwrap()
        );
    }
}
