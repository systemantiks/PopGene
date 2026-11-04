//! Authoritative sequence state. Coordinates are zero-based, half-open.
use serde::{Deserialize, Serialize};
use std::{collections::VecDeque, sync::Arc};
use thiserror::Error;

pub const MAX_BASES: usize = 10_000_000;
pub mod alignment;
pub mod cloning;
pub mod creation;
pub mod gel;
pub mod lineage;
pub mod location;
pub mod metadata;
pub use location::FeatureLocation;
mod feature_metadata;
pub use feature_metadata::{FeatureDirection, FeatureDisplaySegment, FeatureSegmentStyle};
pub mod primer;
pub mod protein;
pub mod restriction;
pub mod trace;
pub mod translation;
pub mod workflow;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum Molecule {
    #[default]
    Dna,
    Rna,
    Protein,
}
impl Molecule {
    pub fn label(self) -> &'static str {
        match self {
            Self::Dna => "DNA",
            Self::Rna => "RNA",
            Self::Protein => "Protein",
        }
    }
    pub fn unit(self) -> &'static str {
        match self {
            Self::Protein => "aa",
            _ => "nt",
        }
    }
    pub fn alphabet(self) -> &'static str {
        match self {
            Self::Dna => "ACGTRYSWKMBDHVN",
            Self::Rna => "ACGURYSWKMBDHVN",
            Self::Protein => "ACDEFGHIKLMNPQRSTVWYBXZJUO*",
        }
    }
    pub fn normalize(self, text: &str) -> Result<String, Error> {
        let mut result = String::with_capacity(text.len().min(MAX_BASES));
        for c in text.chars().filter(|c| !c.is_ascii_whitespace()) {
            let c = c.to_ascii_uppercase();
            if !self.alphabet().contains(c) {
                return Err(Error::MoleculeAlphabet(self.label(), c));
            }
            result.push(c);
            if result.len() > MAX_BASES {
                return Err(Error::Length);
            }
        }
        Ok(result)
    }
}

#[derive(Debug, Error)]
pub enum Error {
    #[error("Invalid DNA character '{0}'. Use IUPAC DNA bases (A C G T R Y S W K M B D H V N).")]
    Alphabet(char),
    #[error("Invalid {0} character '{1}'. Check the selected molecule type and alphabet.")]
    MoleculeAlphabet(&'static str, char),
    #[error(
        "The sequence must contain 1–10,000,000 symbols (or be an explicit empty placeholder)."
    )]
    Length,
    #[error("{0}")]
    Operation(String),
    #[error("Coordinates are outside the sequence or reversed.")]
    Coordinates,
    #[error("Feature '{0}' has an invalid location.")]
    Feature(String),
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Segment {
    pub start: usize,
    pub end: usize,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Feature {
    pub id: u64,
    pub label: String,
    pub kind: String,
    /// Ordered forward-strand segments; `reverse` applies to their whole concatenation.
    pub parts: Vec<Segment>,
    /// Full INSDC location when the simple local segment projection is insufficient.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub location: Option<FeatureLocation>,
    pub reverse: bool,
    pub qualifiers: Vec<(String, Option<String>)>,
    pub edited: bool,
}

impl Feature {
    pub fn len(&self) -> usize {
        self.parts.iter().map(|p| p.end - p.start).sum()
    }
    pub fn is_empty(&self) -> bool {
        self.parts.is_empty()
    }
    pub fn extract(&self, sequence: &str) -> String {
        let joined: String = self
            .parts
            .iter()
            .map(|p| &sequence[p.start..p.end])
            .collect();
        if self.reverse {
            reverse_complement(&joined)
        } else {
            joined
        }
    }
    pub fn location_label(&self) -> String {
        if let Some(location) = &self.location {
            return location.expression();
        }
        if self.parts.is_empty() && self.kind == "primer_bind" {
            return "Unbound oligo".into();
        }
        let s = self
            .parts
            .iter()
            .map(|p| format!("{}..{}", p.start + 1, p.end))
            .collect::<Vec<_>>()
            .join(", ");
        if self.reverse {
            format!("rev {s}")
        } else {
            format!("fwd {s}")
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Source {
    pub file_name: String,
    pub format: String,
    #[serde(skip)]
    pub bytes: Arc<[u8]>,
    #[serde(default)]
    pub retained_blocks: Vec<(u8, usize)>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct TraceCallSet {
    pub label: String,
    pub version: u32,
    pub sequence: String,
    pub qualities: Vec<u8>,
    #[serde(default = "trace_qualities_known")]
    pub qualities_known: bool,
    pub peaks: Vec<usize>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Trace {
    pub calls_version: u32,
    pub qualities: Vec<u8>,
    #[serde(default = "trace_qualities_known")]
    pub qualities_known: bool,
    pub peaks: Vec<usize>,
    pub channel_order: String,
    pub channels: Vec<Vec<u16>>,
    /// Unmodified instrument call families, independent of the current edited calls.
    #[serde(default)]
    pub original_calls: Vec<TraceCallSet>,
    #[serde(default)]
    pub edited: bool,
    /// Optional visible sample interval; raw channel arrays are never cropped.
    #[serde(default)]
    pub sample_range: Option<Segment>,
    /// One flag per active call; empty means no manually introduced calls.
    #[serde(default)]
    pub manual_calls: Vec<bool>,
}
fn trace_qualities_known() -> bool {
    true
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Document {
    #[serde(default)]
    pub annotation_lineage: lineage::Registry,
    pub name: String,
    #[serde(default)]
    pub molecule: Molecule,
    #[serde(default)]
    pub single_stranded: bool,
    /// Recorded whole-template DNA methylation assumptions, not measured occupancy.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub methylation: Option<metadata::DnaMethylation>,
    #[serde(default)]
    pub placeholder: bool,
    pub sequence: Arc<str>,
    /// Per-symbol Phred scores from quality-bearing sequence files, independent of traces.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sequence_qualities: Option<Vec<Option<u8>>>,
    /// Optional measured Phred scores in the opposite strand's 5′→3′ order.
    /// None means unavailable, independently of top-strand score availability.
    #[serde(default)]
    pub opposite_sequence_qualities: Option<Vec<Option<u8>>>,
    #[serde(default)]
    pub sequence_colors: Vec<metadata::SequenceColor>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub duplex_ends: Option<workflow::DuplexEnds>,
    /// Annotations using the full physical duplex extent, including bottom-only residues.
    #[serde(default)]
    pub duplex_annotations: metadata::DuplexAnnotations,
    pub circular: bool,
    pub description: String,
    pub accession: Option<String>,
    pub comments: Vec<String>,
    #[serde(default)]
    pub references: Vec<metadata::Reference>,
    pub features: Vec<Feature>,
    pub warnings: Vec<String>,
    pub source: Option<Source>,
    #[serde(default)]
    pub trace: Option<Trace>,
}

pub fn normalize(text: &str) -> Result<String, Error> {
    let mut result = String::with_capacity(text.len().min(MAX_BASES));
    for c in text.chars().filter(|c| !c.is_ascii_whitespace()) {
        let c = c.to_ascii_uppercase();
        if !"ACGTRYSWKMBDHVN".contains(c) {
            return Err(Error::Alphabet(c));
        }
        result.push(c);
        if result.len() > MAX_BASES {
            return Err(Error::Length);
        }
    }
    Ok(result)
}

pub fn reverse_complement(text: &str) -> String {
    text.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'T' => 'A',
            b'C' => 'G',
            b'G' => 'C',
            b'R' => 'Y',
            b'Y' => 'R',
            b'K' => 'M',
            b'M' => 'K',
            b'B' => 'V',
            b'V' => 'B',
            b'D' => 'H',
            b'H' => 'D',
            b'S' => 'S',
            b'W' => 'W',
            _ => 'N',
        })
        .collect()
}

impl Document {
    pub fn new(name: impl Into<String>, sequence: &str) -> Result<Self, Error> {
        Self::new_typed(name, sequence, Molecule::Dna)
    }
    pub fn new_typed(
        name: impl Into<String>,
        sequence: &str,
        molecule: Molecule,
    ) -> Result<Self, Error> {
        let result = Self {
            annotation_lineage: Default::default(),
            name: name.into(),
            molecule,
            single_stranded: molecule == Molecule::Rna,
            methylation: None,
            placeholder: false,
            sequence: molecule.normalize(sequence)?.into(),
            sequence_qualities: None,
            opposite_sequence_qualities: None,
            sequence_colors: vec![],
            duplex_ends: None,
            duplex_annotations: Default::default(),
            circular: false,
            description: String::new(),
            accession: None,
            comments: vec![],
            references: vec![],
            features: vec![],
            warnings: vec![],
            source: None,
            trace: None,
        };
        result.validate()?;
        Ok(result)
    }
    pub fn validate(&self) -> Result<(), Error> {
        if (self.sequence.is_empty() && !self.placeholder)
            || self.sequence.len() > MAX_BASES
            || (self.placeholder && !self.sequence.is_empty())
        {
            return Err(Error::Length);
        }
        if self.molecule.normalize(&self.sequence)? != self.sequence.as_ref() {
            return Err(Error::Alphabet('?'));
        }
        if self.methylation.is_some() && self.molecule != Molecule::Dna {
            return Err(Error::Operation(
                "Recorded DNA methylation flags require a DNA document.".into(),
            ));
        }
        if self.molecule == Molecule::Protein
            && (self.circular || self.single_stranded || self.features.iter().any(|f| f.reverse))
        {
            return Err(Error::Operation(
                "Protein documents have no nucleotide strands or circular topology in PopGene."
                    .into(),
            ));
        }
        if let Some(qualities) = &self.sequence_qualities
            && (qualities.len() != self.sequence.len()
                || qualities.iter().flatten().any(|q| *q > 93))
        {
            return Err(Error::Operation(
                "Sequence quality scores must be Phred 0–93 or unknown, one per symbol.".into(),
            ));
        }
        if let Some(qualities) = &self.opposite_sequence_qualities {
            let length = self
                .duplex_ends
                .as_ref()
                .map_or(self.sequence.len(), |ends| ends.opposite_sequence.len());
            if self.molecule == Molecule::Protein
                || self.single_stranded
                || self.trace.is_some()
                || qualities.len() != length
                || qualities.iter().flatten().any(|q| *q > 93)
            {
                return Err(Error::Operation("Opposite qualities require paired DNA/RNA and Phred0–93 or unknown scores for every opposite-strand base".into()));
            }
        }
        if self.sequence_colors.len() > 100_000
            || self.sequence_colors.iter().any(|color| {
                self.molecule == Molecule::Protein
                    || color.color_rgb > 0xffffff
                    || color.range.start >= color.range.end
                    || color.range.end > self.sequence.len()
            })
        {
            return Err(Error::Operation(
                "Nucleotide color ranges must be within the sequence and use 24-bit RGB".into(),
            ));
        }
        trace::validate(self)?;
        self.annotation_lineage.validate(self)?;
        if let Some(ends) = &self.duplex_ends {
            if self.molecule != Molecule::Dna
                || self.single_stranded
                || self.circular
                || self.trace.is_some()
            {
                return Err(Error::Operation(
                    "Physical ends require linear double-stranded DNA without trace data".into(),
                ));
            }
            ends.validate(&self.sequence)?;
        }
        let extent = self
            .duplex_ends
            .as_ref()
            .map(|ends| ends.geometry(self.sequence.len()).extent);
        if !self.duplex_annotations.is_empty() && extent.is_none() {
            return Err(Error::Operation(
                "Duplex-extent annotations require physical duplex ends".into(),
            ));
        }
        if self.duplex_annotations.sequence_colors.len() > 100_000
            || self.duplex_annotations.sequence_colors.iter().any(|color| {
                color.color_rgb > 0xffffff
                    || color.range.start >= color.range.end
                    || color.range.end > extent.unwrap_or(0)
            })
        {
            return Err(Error::Operation(
                "Duplex color ranges lie outside the physical extent".into(),
            ));
        }
        let mut ids = std::collections::HashSet::new();
        for (f, length) in self
            .features
            .iter()
            .map(|f| (f, self.sequence.len()))
            .chain(
                self.duplex_annotations
                    .features
                    .iter()
                    .map(|f| (f, extent.unwrap_or(0))),
            )
        {
            if let Some(location) = &f.location {
                location
                    .validate(length, self.circular)
                    .map_err(|e| Error::Operation(e.to_string()))?;
                let (parts, reverse) = location
                    .uniform_projection()
                    .unwrap_or_else(|| (location.local_parts(), false));
                if parts != f.parts || reverse != f.reverse {
                    return Err(Error::Operation(format!(
                        "Feature '{}' has an inconsistent location projection",
                        f.label
                    )));
                }
            }
            if !ids.insert(f.id)
                || (f.parts.is_empty()
                    && f.location.is_none()
                    && !(f.kind == "primer_bind" && f.primer_sequence().is_some()))
                || f.parts.iter().any(|p| p.start >= p.end || p.end > length)
            {
                return Err(Error::Feature(f.label.clone()));
            }
        }
        Ok(())
    }
    pub fn gc(&self) -> f32 {
        let canonical = self
            .sequence
            .bytes()
            .filter(|b| b"ACGTU".contains(b))
            .count();
        let gc = self.sequence.bytes().filter(|b| b"GC".contains(b)).count();
        if canonical == 0 {
            0.0
        } else {
            100.0 * gc as f32 / canonical as f32
        }
    }
    /// Conservative retained size for bounded undo. Count signal/call arrays and
    /// annotation text, rather than treating a chromatogram as base calls alone.
    pub fn snapshot_bytes(&self) -> usize {
        let mut bytes = self.annotation_lineage.snapshot_bytes()
            + (self.sequence_colors.len() + self.duplex_annotations.sequence_colors.len())
                * std::mem::size_of::<metadata::SequenceColor>()
            + self.sequence.len()
            + self.name.len()
            + self.description.len()
            + self
                .sequence_qualities
                .as_ref()
                .map_or(0, |scores| std::mem::size_of_val(scores.as_slice()))
            + self
                .opposite_sequence_qualities
                .as_ref()
                .map_or(0, |scores| std::mem::size_of_val(scores.as_slice()))
            + self.accession.as_ref().map_or(0, String::len);
        bytes = bytes.saturating_add(
            self.comments
                .iter()
                .chain(&self.warnings)
                .map(String::len)
                .sum::<usize>(),
        );
        for feature in self
            .features
            .iter()
            .chain(&self.duplex_annotations.features)
        {
            bytes = bytes.saturating_add(
                feature.label.len()
                    + feature.kind.len()
                    + feature.parts.len() * std::mem::size_of::<Segment>(),
            );
            bytes = bytes.saturating_add(
                feature
                    .qualifiers
                    .iter()
                    .map(|(key, value)| key.len() + value.as_ref().map_or(0, String::len))
                    .sum::<usize>(),
            );
            if let Some(location) = &feature.location {
                bytes = bytes.saturating_add(location.expression().len());
            }
        }
        for reference in &self.references {
            bytes = bytes.saturating_add(reference.title.len());
            for field in [
                &reference.authors,
                &reference.consortium,
                &reference.journal,
                &reference.pubmed,
                &reference.remark,
                &reference.scope,
            ] {
                bytes = bytes.saturating_add(field.as_ref().map_or(0, String::len));
            }
            bytes = bytes.saturating_add(
                reference
                    .extra
                    .iter()
                    .map(|(k, v)| k.len() + v.len())
                    .sum::<usize>(),
            );
        }
        if let Some(ends) = &self.duplex_ends {
            bytes = bytes.saturating_add(
                ends.opposite_sequence.len() + ends.left.overhang.len() + ends.right.overhang.len(),
            );
        }
        if let Some(trace) = &self.trace {
            bytes = bytes.saturating_add(
                trace
                    .channels
                    .iter()
                    .map(|c| c.len() * std::mem::size_of::<u16>())
                    .sum::<usize>(),
            );
            bytes = bytes.saturating_add(
                trace.qualities.len()
                    + trace.peaks.len() * std::mem::size_of::<usize>()
                    + trace.manual_calls.len(),
            );
            for calls in &trace.original_calls {
                bytes = bytes.saturating_add(
                    calls.label.len()
                        + calls.sequence.len()
                        + calls.qualities.len()
                        + calls.peaks.len() * std::mem::size_of::<usize>(),
                );
            }
        }
        bytes
    }
    /// Literal normalized DNA search; circular results may wrap past the origin.
    pub fn find(&self, query: &str) -> Result<Vec<Segment>, Error> {
        let q = self.molecule.normalize(query)?;
        if q.is_empty() || q.len() > self.sequence.len() {
            return Ok(vec![]);
        }
        let mut haystack = self.sequence.to_string();
        if self.circular {
            haystack.push_str(&self.sequence[..q.len() - 1]);
        }
        Ok(haystack
            .as_bytes()
            .windows(q.len())
            .enumerate()
            .filter(|(_, w)| *w == q.as_bytes())
            .take(1000)
            .map(|(start, _)| Segment {
                start,
                end: start + q.len(),
            })
            .collect())
    }
    /// Edit policy: boundary insertions are outside features; internal insertions extend them.
    /// Overlapping annotations are marked changed and their cached translations invalidated.
    pub fn replace(&mut self, start: usize, end: usize, replacement: &str) -> Result<(), Error> {
        if self.trace.is_some() {
            return Err(Error::Operation("Trace sequences are read-only.".into()));
        }
        if start > end || end > self.sequence.len() {
            return Err(Error::Coordinates);
        }
        let bases = self.molecule.normalize(replacement)?;
        let new_len = self.sequence.len() - (end - start) + bases.len();
        if new_len == 0 || new_len > MAX_BASES {
            return Err(Error::Length);
        }
        let ends = self
            .duplex_ends
            .as_ref()
            .map(|ends| ends.replaced(&self.sequence, start, end, &bases))
            .transpose()?;
        if self.sequence[start..end] == bases {
            return Ok(());
        }
        let lineage = lineage::before_edit(self, start, end, bases.len())?;
        let opposite_qualities = self.opposite_sequence_qualities.as_ref().map(|scores| {
            let (a, b) = if let Some(ends) = &self.duplex_ends {
                let (top_prefix, top_suffix, bottom_prefix, _) = ends.overhang_lengths();
                let paired = self.sequence.len() - top_prefix - top_suffix;
                (
                    bottom_prefix + paired - (end - top_prefix),
                    bottom_prefix + paired - (start - top_prefix),
                )
            } else {
                (self.sequence.len() - end, self.sequence.len() - start)
            };
            let mut result = scores.clone();
            result.splice(a..b, std::iter::repeat_n(None, bases.len()));
            result
        });
        let features = metadata::replaced_features(&self.features, start, end, bases.len())?;
        let duplex_annotations = if let Some(duplex) = &self.duplex_ends {
            let offset = duplex.geometry(self.sequence.len()).top.start;
            self.duplex_annotations
                .replaced(start + offset, end + offset, bases.len())?
        } else {
            self.duplex_annotations.clone()
        };
        let mut seq = String::with_capacity(new_len);
        seq.push_str(&self.sequence[..start]);
        seq.push_str(&bases);
        seq.push_str(&self.sequence[end..]);
        if let Some(ends) = &ends {
            ends.validate(&seq)?;
        }
        self.duplex_ends = ends;
        self.opposite_sequence_qualities = opposite_qualities;
        self.duplex_annotations = duplex_annotations;
        self.sequence_colors = self
            .sequence_colors
            .iter()
            .filter_map(|color| color.replaced(start, end, bases.len()))
            .collect();
        if let Some(qualities) = &mut self.sequence_qualities {
            qualities.splice(start..end, std::iter::repeat_n(None, bases.len()));
        }
        self.invalidate_methylation();
        self.sequence = seq.into();
        self.placeholder = false;
        self.features = features;
        self.annotation_lineage = lineage;
        lineage::finish(self, false)?;
        Ok(())
    }
    pub(crate) fn invalidate_methylation(&mut self) {
        if let Some(recorded) = self.methylation.take() {
            self.comments.push(format!(
                "Historical whole-template methylation flags before sequence editing: Dam {}; Dcm {}; EcoKI {}. Current methylation is unknown; these flags do not establish measured occupancy.",
                recorded.dam, recorded.dcm, recorded.eco_ki
            ));
        }
    }
    pub fn revcomp(&mut self) -> Result<(), Error> {
        if self.molecule == Molecule::Protein || self.trace.is_some() {
            return Err(Error::Operation(
                "Reverse complement requires a DNA or RNA sequence without trace data.".into(),
            ));
        }
        if self.duplex_ends.is_some() {
            *self = workflow::reverse_duplex(self)?;
            return Ok(());
        }
        let before = self.clone();
        let n = self.sequence.len();
        let mut features = self.features.clone();
        for f in &mut features {
            translation::reverse_translation_exceptions(f, n);
            f.reverse_display_direction(n);
            if let Some(location) = &f.location {
                let location = location
                    .reverse_complemented(n)
                    .map_err(|e| Error::Operation(e.to_string()))?;
                (f.parts, f.reverse) = location
                    .uniform_projection()
                    .unwrap_or_else(|| (location.local_parts(), false));
                f.location = Some(location);
                continue;
            }
            f.parts = f
                .parts
                .iter()
                .rev()
                .map(|p| Segment {
                    start: n - p.end,
                    end: n - p.start,
                })
                .collect();
            f.reverse = !f.reverse;
        }
        self.sequence = self.complement(&self.sequence).into();
        for color in &mut self.sequence_colors {
            color.range = Segment {
                start: n - color.range.end,
                end: n - color.range.start,
            };
            color.reverse = !color.reverse;
        }
        if self.single_stranded {
            if let Some(qualities) = &mut self.sequence_qualities {
                qualities.reverse();
            }
        } else {
            std::mem::swap(
                &mut self.sequence_qualities,
                &mut self.opposite_sequence_qualities,
            );
        }
        self.features = features;
        lineage::reverse(&before, self)?;
        Ok(())
    }
    pub fn complement(&self, sequence: &str) -> String {
        if self.molecule == Molecule::Rna {
            reverse_complement(&sequence.replace('U', "T")).replace('T', "U")
        } else {
            reverse_complement(sequence)
        }
    }
    /// Exact extraction never interprets order/bond/alternative sites as joins.
    /// Remote accessions require FeatureLocation::extract_with and an explicit resolver.
    pub fn feature_sequence(&self, feature: &Feature) -> Result<String, Error> {
        if let Some(location) = &feature.location {
            return location
                .extract(&self.sequence, self.molecule)
                .map_err(|e| Error::Operation(e.to_string()));
        }
        if feature
            .parts
            .iter()
            .any(|p| p.start >= p.end || p.end > self.sequence.len())
        {
            return Err(Error::Coordinates);
        }
        let sequence: String = feature
            .parts
            .iter()
            .map(|p| &self.sequence[p.start..p.end])
            .collect();
        Ok(if feature.reverse {
            self.complement(&sequence)
        } else {
            sequence
        })
    }
    pub fn unit(&self) -> &'static str {
        if self.molecule == Molecule::Dna && !self.single_stranded {
            "bp"
        } else {
            self.molecule.unit()
        }
    }
    pub fn summary(&self) -> String {
        if self.placeholder {
            return format!("{} placeholder · no sequence", self.molecule.label());
        }
        if self.molecule == Molecule::Protein {
            return format!("{} aa · protein", self.sequence.len());
        }
        if let Some(ends) = &self.duplex_ends {
            return format!(
                "{} nt top · {} nt opposite · {}-column duplex · linear DNA",
                self.sequence.len(),
                ends.opposite_sequence.len(),
                ends.geometry(self.sequence.len()).extent
            );
        }
        format!(
            "{} {} · {} {} · {} · {:.1}% GC",
            self.sequence.len(),
            self.unit(),
            if self.single_stranded { "ss" } else { "ds" },
            self.molecule.label(),
            if self.circular { "circular" } else { "linear" },
            self.gc()
        )
    }
    pub fn annotate(&mut self, label: String, part: Segment, reverse: bool) -> Result<(), Error> {
        if self.molecule == Molecule::Protein && reverse {
            return Err(Error::Operation(
                "Protein annotations have no reverse strand.".into(),
            ));
        }
        if part.start >= part.end || part.end > self.sequence.len() {
            return Err(Error::Coordinates);
        }
        let id = self
            .features
            .iter()
            .chain(&self.duplex_annotations.features)
            .map(|f| f.id)
            .chain(self.annotation_lineage.pieces.iter().map(|p| p.feature_id))
            .max()
            .unwrap_or(0)
            .checked_add(1)
            .ok_or_else(|| Error::Operation("Feature identifiers are exhausted".into()))?;
        self.annotation_lineage.record_feature(id);
        self.features.push(Feature {
            location: None,
            id,
            label: label.clone(),
            kind: "misc_feature".into(),
            parts: vec![part],
            reverse,
            qualifiers: vec![("label".into(), Some(label))],
            edited: false,
        });
        lineage::limited_warning(self);
        Ok(())
    }
}

pub struct Editor {
    pub document: Document,
    undo: VecDeque<Document>,
    redo: Vec<Document>,
}
impl Editor {
    pub fn new(document: Document) -> Self {
        Self {
            document,
            undo: VecDeque::new(),
            redo: vec![],
        }
    }
    pub fn apply(
        &mut self,
        operation: impl FnOnce(&mut Document) -> Result<(), Error>,
    ) -> Result<(), Error> {
        let mut next = self.document.clone();
        operation(&mut next)?;
        next.validate()?;
        if next == self.document {
            return Ok(());
        }
        self.undo
            .push_back(std::mem::replace(&mut self.document, next));
        self.redo.clear();
        while self.undo.len() > 32
            || self
                .undo
                .iter()
                .map(Document::snapshot_bytes)
                .sum::<usize>()
                > 32 * 1024 * 1024
        {
            self.undo.pop_front();
        }
        Ok(())
    }
    pub fn can_undo(&self) -> bool {
        !self.undo.is_empty()
    }
    pub fn can_redo(&self) -> bool {
        !self.redo.is_empty()
    }
    pub fn undo(&mut self) {
        if let Some(doc) = self.undo.pop_back() {
            self.redo.push(std::mem::replace(&mut self.document, doc));
        }
    }
    pub fn redo(&mut self) {
        if let Some(doc) = self.redo.pop() {
            self.undo
                .push_back(std::mem::replace(&mut self.document, doc));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn sample() -> Document {
        let mut d = Document::new("test", "ACGTACGTACGT").unwrap();
        d.annotate("gene".into(), Segment { start: 2, end: 8 }, false)
            .unwrap();
        d
    }
    #[test]
    fn complement_iupac() {
        assert_eq!(reverse_complement("ACGTRYMKBDHVN"), "NBDHVMKRYACGT");
    }
    #[test]
    fn recorded_methylation_survives_rotation_and_noops_but_not_sequence_edits() {
        let mut document = sample();
        assert_eq!(document.methylation, None);
        document.methylation = Some(metadata::DnaMethylation {
            dam: true,
            dcm: false,
            eco_ki: true,
        });
        document.sequence_qualities = Some(vec![Some(30); document.sequence.len()]);
        let original = document.clone();
        document.replace(0, 4, "acgt").unwrap();
        assert_eq!(document, original);
        document.revcomp().unwrap();
        assert_eq!(document.methylation, original.methylation);
        document.revcomp().unwrap();
        assert_eq!(document, original);
        let mut editor = Editor::new(document);
        editor.apply(|d| d.replace(0, 1, "T")).unwrap();
        assert_eq!(editor.document.methylation, None);
        assert!(
            editor
                .document
                .comments
                .iter()
                .any(|c| c.contains("Dam true; Dcm false; EcoKI true"))
        );
        editor.undo();
        assert_eq!(editor.document, original);
        let mut legacy = serde_json::to_value(original).unwrap();
        legacy.as_object_mut().unwrap().remove("methylation");
        assert_eq!(
            serde_json::from_value::<Document>(legacy)
                .unwrap()
                .methylation,
            None
        );
        let mut rna = Document::new_typed("RNA", "ACGU", Molecule::Rna).unwrap();
        rna.methylation = Some(Default::default());
        assert!(rna.validate().is_err());
    }
    #[test]
    fn invalid_edit_is_atomic() {
        let mut e = Editor::new(sample());
        let old = e.document.clone();
        assert!(e.apply(|d| d.replace(0, 2, "XYZ")).is_err());
        assert_eq!(old, e.document);
        assert!(!e.can_undo());
    }
    #[test]
    fn internal_insert_extends_and_marks_feature() {
        let mut d = sample();
        d.replace(4, 4, "AA").unwrap();
        assert_eq!(d.features[0].parts, vec![Segment { start: 2, end: 10 }]);
        assert!(d.features[0].edited);
    }
    #[test]
    fn boundary_insert_stays_outside_feature() {
        let mut d = sample();
        d.replace(2, 2, "AA").unwrap();
        assert_eq!(d.features[0].parts, vec![Segment { start: 4, end: 10 }]);
        assert!(!d.features[0].edited);
    }
    #[test]
    fn deletion_removes_only_deleted_features() {
        let mut d = sample();
        d.replace(2, 8, "").unwrap();
        assert!(d.features.is_empty());
        assert_eq!(d.sequence.as_ref(), "ACACGT");
    }
    #[test]
    fn undo_restores_annotations() {
        let mut e = Editor::new(sample());
        let old = e.document.clone();
        e.apply(|d| d.replace(0, 5, "N")).unwrap();
        let new = e.document.clone();
        e.undo();
        assert_eq!(old, e.document);
        e.redo();
        assert_eq!(new, e.document);
    }
    #[test]
    fn reverse_compound_extraction_is_invariant() {
        let mut d = sample();
        d.features[0].parts = vec![Segment { start: 10, end: 12 }, Segment { start: 0, end: 3 }];
        let old = d.clone();
        let extracted = d.features[0].extract(&d.sequence);
        d.revcomp().unwrap();
        assert_eq!(d.features[0].extract(&d.sequence), extracted);
        d.revcomp().unwrap();
        assert_eq!(d, old);
    }
    #[test]
    fn circular_search_crosses_origin() {
        let mut d = Document::new("circle", "AAAACCCC").unwrap();
        d.circular = true;
        assert_eq!(d.find("CCAA").unwrap(), vec![Segment { start: 6, end: 10 }]);
        d.circular = false;
        assert!(d.find("CCAA").unwrap().is_empty());
    }
    #[test]
    fn trace_undo_budget_counts_signal_arrays_and_restores_latest_calls() {
        let mut document = Document::new("large trace", "AAAAA").unwrap();
        document.trace = Some(Trace {
            calls_version: 2,
            qualities: vec![40; 5],
            qualities_known: true,
            peaks: vec![10, 20, 30, 40, 50],
            channel_order: "ACGT".into(),
            channels: (0..4)
                .map(|channel| vec![100 + channel; 1_000_000])
                .collect(),
            original_calls: vec![],
            edited: false,
            sample_range: None,
            manual_calls: vec![],
        });
        let mut editor = Editor::new(document);
        for index in 0..5 {
            editor
                .apply(|doc| trace::replace_call(doc, index, "C"))
                .unwrap();
        }
        assert!(
            editor.undo.len() < 5,
            "Large raw signal snapshots must evict old states despite the five-base sequence"
        );
        assert!(
            editor
                .undo
                .iter()
                .map(Document::snapshot_bytes)
                .sum::<usize>()
                <= 32 * 1024 * 1024
        );
        assert_eq!(editor.document.sequence.as_ref(), "CCCCC");
        editor.undo();
        assert_eq!(editor.document.sequence.as_ref(), "CCCCA");
        let trace = editor.document.trace.as_ref().unwrap();
        assert_eq!(trace.qualities, vec![0, 0, 0, 0, 40]);
        assert_eq!(trace.peaks, vec![10, 20, 30, 40, 50]);
        for (channel, samples) in trace.channels.iter().enumerate() {
            assert_eq!(samples.len(), 1_000_000);
            assert!(samples.iter().all(|&value| value == 100 + channel as u16));
        }
        editor.redo();
        assert_eq!(editor.document.sequence.as_ref(), "CCCCC");
    }
}
