//! Authoritative sequence state. Coordinates are zero-based, half-open.
use serde::{Deserialize, Serialize};
use std::{collections::VecDeque, sync::Arc};
use thiserror::Error;

pub const MAX_BASES: usize = 10_000_000;
pub mod creation;
mod feature_metadata;
pub use feature_metadata::{FeatureDirection, FeatureDisplaySegment};
pub mod protein;
pub mod restriction;
pub mod translation;

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
pub struct Trace {
    pub calls_version: u32,
    pub qualities: Vec<u8>,
    pub peaks: Vec<usize>,
    pub channel_order: String,
    pub channels: Vec<Vec<u16>>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Document {
    pub name: String,
    #[serde(default)]
    pub molecule: Molecule,
    #[serde(default)]
    pub single_stranded: bool,
    #[serde(default)]
    pub placeholder: bool,
    pub sequence: Arc<str>,
    pub circular: bool,
    pub description: String,
    pub accession: Option<String>,
    pub comments: Vec<String>,
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
            name: name.into(),
            molecule,
            single_stranded: molecule == Molecule::Rna,
            placeholder: false,
            sequence: molecule.normalize(sequence)?.into(),
            circular: false,
            description: String::new(),
            accession: None,
            comments: vec![],
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
        if self.molecule == Molecule::Protein
            && (self.circular || self.single_stranded || self.features.iter().any(|f| f.reverse))
        {
            return Err(Error::Operation(
                "Protein documents have no nucleotide strands or circular topology in PopGene."
                    .into(),
            ));
        }
        if let Some(trace) = &self.trace {
            let samples = trace.channels.first().map_or(0, Vec::len);
            if self.molecule != Molecule::Dna
                || self.placeholder
                || trace.channels.len() != 4
                || samples == 0
                || trace.channels.iter().any(|c| c.len() != samples)
                || trace.peaks.len() != self.sequence.len()
                || trace.qualities.len() != self.sequence.len()
                || trace.peaks.iter().any(|&p| p >= samples)
                || trace.peaks.windows(2).any(|p| p[0] > p[1])
                || trace.channel_order.len() != 4
                || !b"ACGT"
                    .iter()
                    .all(|b| trace.channel_order.as_bytes().contains(b))
            {
                return Err(Error::Feature("Invalid trace arrays".into()));
            }
        }
        let mut ids = std::collections::HashSet::new();
        for f in &self.features {
            if !ids.insert(f.id)
                || f.parts.is_empty()
                || f.parts
                    .iter()
                    .any(|p| p.start >= p.end || p.end > self.sequence.len())
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
        let delta = bases.len() as isize - (end - start) as isize;
        let shift = |v: usize| (v as isize + delta) as usize;
        let mut features = self.features.clone();
        for feature in &mut features {
            feature.shift_display_segments(start, end, bases.len());
            let mut changed = false;
            for part in &mut feature.parts {
                if part.end <= start {
                    continue;
                }
                if part.start >= end {
                    part.start = shift(part.start);
                    part.end = shift(part.end);
                    continue;
                }
                changed = true;
                part.start = part.start.min(start);
                part.end = if part.end > end {
                    shift(part.end)
                } else {
                    start + bases.len()
                };
            }
            feature.parts.retain(|p| p.start < p.end);
            if changed {
                feature.edited = true;
                feature.qualifiers.retain(|(k, _)| k != "translation");
            }
        }
        features.retain(|f| !f.parts.is_empty());
        let mut seq = String::with_capacity(new_len);
        seq.push_str(&self.sequence[..start]);
        seq.push_str(&bases);
        seq.push_str(&self.sequence[end..]);
        self.sequence = seq.into();
        self.placeholder = false;
        self.features = features;
        Ok(())
    }
    pub fn revcomp(&mut self) -> Result<(), Error> {
        if self.molecule == Molecule::Protein || self.trace.is_some() {
            return Err(Error::Operation(
                "Reverse complement requires a DNA or RNA sequence without trace data.".into(),
            ));
        }
        let n = self.sequence.len();
        self.sequence = self.complement(&self.sequence).into();
        for f in &mut self.features {
            f.reverse_display_direction(n);
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
        Ok(())
    }
    pub fn complement(&self, sequence: &str) -> String {
        if self.molecule == Molecule::Rna {
            reverse_complement(&sequence.replace('U', "T")).replace('T', "U")
        } else {
            reverse_complement(sequence)
        }
    }
    pub fn feature_sequence(&self, feature: &Feature) -> String {
        let sequence: String = feature
            .parts
            .iter()
            .map(|p| &self.sequence[p.start..p.end])
            .collect();
        if feature.reverse {
            self.complement(&sequence)
        } else {
            sequence
        }
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
        let id = self.features.iter().map(|f| f.id).max().unwrap_or(0) + 1;
        self.features.push(Feature {
            id,
            label: label.clone(),
            kind: "misc_feature".into(),
            parts: vec![part],
            reverse,
            qualifiers: vec![("label".into(), Some(label))],
            edited: false,
        });
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
            || self.undo.iter().map(|d| d.sequence.len()).sum::<usize>() > 32 * 1024 * 1024
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
}
