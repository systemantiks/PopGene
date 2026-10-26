//! Bounded, deterministic nucleotide translation.
//!
//! All 27 active NCBI genetic codes, gc.prt v4.6 (verified 2026-10-25):
//! https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
//! Exact local transl_except and joined ribosomal-slippage coordinates follow
//! https://www.insdc.org/submitting-standards/feature-table/ (v11.4, April 2026).
use crate::{Document, Feature, Molecule};
use std::{fmt, ops::Range};

mod exceptions;
mod flanks;
mod genetic_codes;
pub use flanks::{FlankingCodon, cds_flanking_codons};
mod orfs;
#[cfg(test)]
mod regression;
pub use exceptions::{
    project_translation_exceptions, remap_translation_exceptions, reverse_translation_exceptions,
};
pub use orfs::{Orf, OrfOptions, OrfSearch, find_orfs};

/// One of the 27 active NCBI genetic codes. Retired IDs are not accepted.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GeneticCode(usize);
impl GeneticCode {
    pub const STANDARD: Self = Self(0);
    pub const VERTEBRATE_MITOCHONDRIAL: Self = Self(1);
    pub const BACTERIAL: Self = Self(8);
    pub fn from_id(id: u16) -> Option<Self> {
        genetic_codes::CODES
            .iter()
            .position(|code| code.id == id)
            .map(Self)
    }
    pub fn all() -> impl Iterator<Item = Self> {
        (0..genetic_codes::CODES.len()).map(Self)
    }
    pub fn id(self) -> u16 {
        genetic_codes::CODES[self.0].id
    }
    pub fn name(self) -> &'static str {
        genetic_codes::CODES[self.0].name
    }
    pub fn has_context_dependent_stops(self) -> bool {
        let data = &genetic_codes::CODES[self.0];
        data.starts
            .iter()
            .zip(data.amino_acids)
            .any(|(&flag, &aa)| flag == b'*' && aa != b'*')
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TranslationError(pub String);
impl fmt::Display for TranslationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}
impl std::error::Error for TranslationError {}

fn mask(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'T' | b'U' => 1,
        b'C' => 2,
        b'A' => 4,
        b'G' => 8,
        b'Y' => 3,
        b'W' => 5,
        b'K' => 9,
        b'M' => 6,
        b'S' => 10,
        b'R' => 12,
        b'H' => 7,
        b'B' => 11,
        b'D' => 13,
        b'V' => 14,
        b'N' => 15,
        _ => 0,
    }
}
fn concrete(index: usize, code: GeneticCode) -> u8 {
    genetic_codes::CODES[code.0].amino_acids[index]
}

/// A predicate must hold for every possible concrete codon.
fn all_codons(codon: [u8; 3], predicate: impl Fn(usize) -> bool) -> bool {
    let masks = codon.map(mask);
    if masks.contains(&0) {
        return false;
    }
    if masks.iter().all(|mask| mask.count_ones() == 1) {
        return predicate(
            (masks[0].trailing_zeros() * 16
                + masks[1].trailing_zeros() * 4
                + masks[2].trailing_zeros()) as usize,
        );
    }
    for first in 0..4 {
        for second in 0..4 {
            for third in 0..4 {
                if masks[0] & (1 << first) != 0
                    && masks[1] & (1 << second) != 0
                    && masks[2] & (1 << third) != 0
                    && !predicate(first * 16 + second * 4 + third)
                {
                    return false;
                }
            }
        }
    }
    true
}

pub fn is_start_codon(codon: [u8; 3], code: GeneticCode, alternative_starts: bool) -> bool {
    all_codons(codon, |index| {
        if alternative_starts {
            genetic_codes::CODES[code.0].starts[index] == b'M'
        } else {
            index == 35
        }
    })
}
/// Terminal flags include dual-coding codons in NCBI tables 27, 28 and 31.
/// Such codons translate as amino acids internally; their terminal use requires
/// an annotated CDS boundary or an explicitly labelled ORF prediction.
pub fn is_stop_codon(codon: [u8; 3], code: GeneticCode) -> bool {
    all_codons(codon, |index| {
        genetic_codes::CODES[code.0].starts[index] == b'*'
    })
}

/// Ambiguity resolves only when every concrete codon has the same assignment.
/// Otherwise X is returned, including stop/amino-acid mixtures and invalid input.
pub fn translate_codon(codon: [u8; 3], code: GeneticCode) -> u8 {
    let masks = codon.map(mask);
    if masks.contains(&0) {
        return b'X';
    }
    if masks.iter().all(|mask| mask.count_ones() == 1) {
        return concrete(
            (masks[0].trailing_zeros() * 16
                + masks[1].trailing_zeros() * 4
                + masks[2].trailing_zeros()) as usize,
            code,
        );
    }
    let mut result = None;
    for first in 0..4 {
        for second in 0..4 {
            for third in 0..4 {
                if masks[0] & (1 << first) == 0
                    || masks[1] & (1 << second) == 0
                    || masks[2] & (1 << third) == 0
                {
                    continue;
                }
                let amino_acid = concrete(first * 16 + second * 4 + third, code);
                if result.is_some_and(|previous| previous != amino_acid) {
                    return b'X';
                }
                result = Some(amino_acid);
            }
        }
    }
    result.unwrap_or(b'X')
}

pub fn complement(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' | b'U' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        b'S' => b'S',
        b'W' => b'W',
        _ => b'N',
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrameCodon {
    pub center: usize,
    /// A contiguous display interval; wrapped codons highlight their center only.
    pub start: usize,
    pub end: usize,
    pub amino_acid: u8,
    pub positions: [usize; 3],
}
impl FrameCodon {
    pub fn selection(&self) -> Vec<crate::Segment> {
        let mut positions = self.positions;
        positions.sort_unstable();
        let mut result: Vec<crate::Segment> = vec![];
        for p in positions {
            if let Some(last) = result.last_mut()
                && last.end == p
            {
                last.end = p + 1;
                continue;
            }
            result.push(crate::Segment {
                start: p,
                end: p + 1,
            });
        }
        result
    }
}

/// Backward-compatible linear standard-code frame calculation.
pub fn frame_codons(sequence: &str, frame: i8, range: Range<usize>) -> Vec<FrameCodon> {
    frame_codons_with_options(sequence, frame, range, GeneticCode::STANDARD, false)
}
/// Each frame is anchored at its strand's 5′ end. Circular codons may cross the
/// origin once; when length is not divisible by three, the phase changes across
/// the seam. Codons are returned only when their central base is in the row.
pub fn frame_codons_with_options(
    sequence: &str,
    frame: i8,
    range: Range<usize>,
    code: GeneticCode,
    circular: bool,
) -> Vec<FrameCodon> {
    let n = sequence.len();
    if !matches!(frame,-3..=-1|1..=3) || n < 3 {
        return vec![];
    }
    let phase = usize::from(frame.unsigned_abs() - 1);
    let mut result = vec![];
    for center in range.start.min(n)..range.end.min(n) {
        let coding_center = if frame < 0 { n - 1 - center } else { center };
        let offset = if circular {
            (coding_center + n - 1) % n
        } else {
            let Some(offset) = coding_center.checked_sub(1) else {
                continue;
            };
            offset
        };
        if offset < phase || (offset - phase) % 3 != 0 || (!circular && offset + 3 > n) {
            continue;
        }
        let positions = std::array::from_fn(|i| {
            let p = (offset + i) % n;
            if frame < 0 { n - 1 - p } else { p }
        });
        let codon = positions.map(|p| {
            if frame < 0 {
                complement(sequence.as_bytes()[p])
            } else {
                sequence.as_bytes()[p]
            }
        });
        let lo = *positions.iter().min().unwrap();
        let hi = *positions.iter().max().unwrap();
        let (start, end) = if hi - lo == 2 {
            (lo, hi + 1)
        } else {
            (center, center + 1)
        };
        result.push(FrameCodon {
            center,
            start,
            end,
            amino_acid: translate_codon(codon, code),
            positions,
        });
    }
    result
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CdsTranslation {
    pub amino_acids: Vec<u8>,
    pub phase: usize,
    pub genetic_code: u16,
    pub trailing_bases: usize,
    pub exception_codons: Vec<usize>,
    pub terminal_partial_stop: bool,
}

fn qualifier<'a>(feature: &'a Feature, name: &str) -> Result<Option<&'a str>, TranslationError> {
    let mut result = None;
    for (_, value) in feature.qualifiers.iter().filter(|(key, _)| key == name) {
        let value = value
            .as_deref()
            .ok_or_else(|| TranslationError(format!("/{name} has no value")))?
            .trim();
        if result.is_some_and(|previous| previous != value) {
            return Err(TranslationError(format!("conflicting /{name} values")));
        }
        result = Some(value);
    }
    Ok(result)
}

/// INSDC phase and observed SnapGene signed readingFrame must agree. Native
/// values +1/+2/+3 denote the forward strand and -1/-2/-3 the reverse strand.
/// Native 8.2.3 fixtures with both whole-sequence and shifted CDS bounds establish
/// that the magnitude is feature-relative codon_start, not a global frame.
pub fn cds_phase(feature: &Feature) -> Result<usize, TranslationError> {
    let recorded = match qualifier(feature, "codon_start")? {
        None => None,
        Some("1") => Some(0),
        Some("2") => Some(1),
        Some("3") => Some(2),
        _ => {
            return Err(TranslationError(
                "invalid /codon_start (expected 1, 2 or 3)".into(),
            ));
        }
    };
    let native = match qualifier(feature, "snapgene_readingFrame")? {
        None => None,
        Some(value) => {
            let frame = value
                .parse::<i8>()
                .ok()
                .filter(|frame| matches!(frame, -3..=-1 | 1..=3))
                .ok_or_else(|| {
                    TranslationError(
                        "invalid native readingFrame (expected signed 1, 2 or 3)".into(),
                    )
                })?;
            if (frame < 0) != feature.reverse {
                return Err(TranslationError(
                    "native readingFrame conflicts with the CDS strand".into(),
                ));
            }
            Some(usize::from(frame.unsigned_abs() - 1))
        }
    };
    if recorded
        .zip(native)
        .is_some_and(|(recorded, native)| recorded != native)
    {
        return Err(TranslationError(
            "native readingFrame conflicts with /codon_start".into(),
        ));
    }
    Ok(recorded.or(native).unwrap_or(0))
}
/// Calculate from the current sequence, never from a cached /translation value.
/// The result includes terminal/internal '*' and ignores incomplete trailing
/// codons. Whole ordered parts are concatenated then reverse-complemented.
pub fn translate_cds(
    document: &Document,
    feature: &Feature,
) -> Result<CdsTranslation, TranslationError> {
    if document.molecule == Molecule::Protein || !feature.kind.eq_ignore_ascii_case("CDS") {
        return Err(TranslationError("not a nucleotide CDS".into()));
    }
    for (key, value) in &feature.qualifiers {
        let needs_review = match key.as_str() {
            "exception" => {
                value.as_deref().map(str::trim) != Some("ribosomal slippage")
                    || feature.parts.len() < 2
            }
            "ribosomal_slippage" => feature.parts.len() < 2,
            "codon" | "pseudo" | "pseudogene" | exceptions::REVIEW => true,
            _ => false,
        };
        if needs_review {
            return Err(TranslationError(format!(
                "/{key} requires review; no calculated translation is asserted"
            )));
        }
    }
    let (partial_start, partial_end) = if let Some(location) = &feature.location {
        location
            .validate(document.sequence.len(), document.circular)
            .map_err(|e| TranslationError(e.to_string()))?;
        let (parts,reverse)=location.uniform_projection().ok_or_else(||TranslationError("CDS location is not a resolved local join on one strand; remote, mixed-strand, order, gap or alternative locations need review".into()))?;
        if parts != feature.parts || reverse != feature.reverse {
            return Err(TranslationError(
                "CDS location and displayed coordinate projection disagree".into(),
            ));
        }
        location.partial_ends()
    } else {
        (false, false)
    };
    let phase = cds_phase(feature)?;
    let table = qualifier(feature, "transl_table")?
        .unwrap_or("1")
        .parse::<u16>()
        .map_err(|_| TranslationError("invalid /transl_table".into()))?;
    let code = GeneticCode::from_id(table).ok_or_else(|| {
        TranslationError(format!(
            "genetic code table {table} is not an active NCBI code"
        ))
    })?;
    if feature.parts.is_empty()
        || feature
            .parts
            .iter()
            .any(|part| part.start >= part.end || part.end > document.sequence.len())
    {
        return Err(TranslationError("invalid CDS location".into()));
    }
    if feature.len() > crate::MAX_BASES {
        return Err(TranslationError(
            "CDS exceeds the 10-million-base translation limit".into(),
        ));
    }
    let mut coding = Vec::with_capacity(feature.len());
    for part in &feature.parts {
        coding.extend_from_slice(&document.sequence.as_bytes()[part.start..part.end]);
    }
    if feature.reverse {
        coding.reverse();
        coding.iter_mut().for_each(|base| *base = complement(*base));
    }
    let coding = coding.get(phase..).unwrap_or_default();
    let mut amino_acids = Vec::with_capacity(coding.len() / 3);
    let context_dependent_stops = code.has_context_dependent_stops();
    for (index, &codon) in coding.as_chunks::<3>().0.iter().enumerate() {
        let aa = if context_dependent_stops
            && !partial_end
            && coding.len() % 3 == 0
            && (index + 1) * 3 == coding.len()
            && is_stop_codon(codon, code)
        {
            b'*'
        } else if index == 0 && phase == 0 && !partial_start && is_start_codon(codon, code, true) {
            b'M'
        } else {
            translate_codon(codon, code)
        };
        amino_acids.push(aa);
    }
    let mut targets = std::collections::BTreeMap::new();
    for (_, value) in feature
        .qualifiers
        .iter()
        .filter(|(key, _)| key == "transl_except")
    {
        if targets.len() >= 1024 {
            return Err(TranslationError(
                "too many translation exceptions (limit 1024)".into(),
            ));
        }
        let exception = exceptions::parse(
            value
                .as_deref()
                .ok_or_else(|| TranslationError("/transl_except has no value".into()))?,
        )?;
        if exception.reverse != feature.reverse {
            return Err(TranslationError(
                "/transl_except strand does not match the CDS".into(),
            ));
        }
        if let Some(&(aa, _)) = targets.get(&exception.positions)
            && aa != exception.amino_acid
        {
            return Err(TranslationError(
                "conflicting /transl_except residues at the same codon".into(),
            ));
        }
        targets.insert(exception.positions, (exception.amino_acid, 0_usize));
    }
    let mut exception_codons = vec![];
    let mut terminal_partial_stop = false;
    if !targets.is_empty() {
        let mut positions: Vec<_> = feature
            .parts
            .iter()
            .flat_map(|part| part.start..part.end)
            .collect();
        if feature.reverse {
            positions.reverse();
        }
        for (index, codon_positions) in positions
            .get(phase..)
            .unwrap_or_default()
            .chunks(3)
            .enumerate()
        {
            if let Some((aa, hits)) = targets.get_mut(codon_positions) {
                *hits += 1;
                if codon_positions.len() == 3 {
                    amino_acids[index] = *aa;
                } else {
                    // INSDC allows a final one/two-base TERM completed by poly-A.
                    // Only trailing T or TA can become TAA after adding A bases.
                    let partial = &coding[index * 3..];
                    if *aa != b'*'
                        || !matches!(partial, [b'T'] | [b'T', b'A'] | [b'U'] | [b'U', b'A'])
                    {
                        return Err(TranslationError("partial /transl_except TERM must be trailing T or TA completed by polyadenylation".into()));
                    }
                    amino_acids.push(b'*');
                    terminal_partial_stop = true;
                }
                exception_codons.push(index);
            }
        }
        if targets.values().any(|(_, hits)| *hits != 1) {
            return Err(TranslationError(
                "/transl_except must identify exactly one codon in the CDS reading frame".into(),
            ));
        }
    }
    Ok(CdsTranslation {
        amino_acids,
        phase,
        genetic_code: code.id(),
        trailing_bases: if terminal_partial_stop {
            0
        } else {
            coding.len() % 3
        },
        exception_codons,
        terminal_partial_stop,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Segment;
    fn cds(sequence: &str, reverse: bool) -> (Document, Feature) {
        let doc = Document::new("translation", sequence).unwrap();
        let feature = Feature {
            location: None,
            id: 1,
            label: "CDS".into(),
            kind: "CDS".into(),
            parts: vec![Segment {
                start: 0,
                end: sequence.len(),
            }],
            reverse,
            qualifiers: vec![],
            edited: false,
        };
        (doc, feature)
    }
    #[test]
    fn standard_assignments_cover_all_codons_by_amino_acid_family() {
        // Independently enumerated families, not the implementation's TCAG order.
        for (amino, codons) in [
            (b'A', "GCT GCC GCA GCG"),
            (b'C', "TGT TGC"),
            (b'D', "GAT GAC"),
            (b'E', "GAA GAG"),
            (b'F', "TTT TTC"),
            (b'G', "GGT GGC GGA GGG"),
            (b'H', "CAT CAC"),
            (b'I', "ATT ATC ATA"),
            (b'K', "AAA AAG"),
            (b'L', "TTA TTG CTT CTC CTA CTG"),
            (b'M', "ATG"),
            (b'N', "AAT AAC"),
            (b'P', "CCT CCC CCA CCG"),
            (b'Q', "CAA CAG"),
            (b'R', "CGT CGC CGA CGG AGA AGG"),
            (b'S', "TCT TCC TCA TCG AGT AGC"),
            (b'T', "ACT ACC ACA ACG"),
            (b'V', "GTT GTC GTA GTG"),
            (b'W', "TGG"),
            (b'Y', "TAT TAC"),
            (b'*', "TAA TAG TGA"),
        ] {
            for codon in codons.split_whitespace() {
                assert_eq!(
                    translate_codon(codon.as_bytes().try_into().unwrap(), GeneticCode::STANDARD),
                    amino,
                    "{codon}"
                );
            }
        }
    }
    #[test]
    fn ambiguity_rna_and_alternative_codes_are_explicit() {
        for (codon, aa) in [
            ("GCN", b'A'),
            ("TAR", b'*'),
            ("ATH", b'I'),
            ("TAN", b'X'),
            ("NNN", b'X'),
            ("aug", b'M'),
        ] {
            assert_eq!(
                translate_codon(codon.as_bytes().try_into().unwrap(), GeneticCode::STANDARD),
                aa
            );
        }
        for (codon, aa) in [("TGA", b'W'), ("ATA", b'M'), ("AGA", b'*'), ("AGG", b'*')] {
            assert_eq!(
                translate_codon(
                    codon.as_bytes().try_into().unwrap(),
                    GeneticCode::VERTEBRATE_MITOCHONDRIAL
                ),
                aa
            );
        }
        assert_eq!(GeneticCode::from_id(7), None);
        assert_eq!(GeneticCode::all().count(), 27);
    }
    #[test]
    fn six_frames_use_complete_codons_and_reverse_complement_coordinates() {
        // ATG AA A TAA: +1 = MK*, +2 = *N, +3 = EI.
        let sequence = "ATGAAATAA";
        for (frame, expected) in [
            (1, "MK*"),
            (2, "*N"),
            (3, "EI"),
            (-1, "LFH"),
            (-2, "YF"),
            (-3, "IS"),
        ] {
            let mut glyphs = frame_codons(sequence, frame, 0..sequence.len());
            if frame < 0 {
                glyphs.reverse();
            }
            assert_eq!(
                glyphs
                    .iter()
                    .map(|c| c.amino_acid as char)
                    .collect::<String>(),
                expected,
                "frame {frame}"
            );
        }
        assert_eq!(
            frame_codons(sequence, 1, 4..5),
            vec![FrameCodon {
                center: 4,
                start: 3,
                end: 6,
                amino_acid: b'K',
                positions: [3, 4, 5]
            }]
        );
        assert!(frame_codons("AT", 1, 0..2).is_empty());
        assert!(frame_codons(sequence, 0, 0..99).is_empty());
    }
    #[test]
    fn cds_reads_current_sequence_frame_and_ordered_compounds() {
        let (mut doc, mut feature) = cds("CCATGAAATAAGGG", false);
        feature.parts = vec![Segment { start: 2, end: 5 }, Segment { start: 5, end: 11 }];
        feature
            .qualifiers
            .push(("translation".into(), Some("WRONG".into())));
        assert_eq!(translate_cds(&doc, &feature).unwrap().amino_acids, b"MK*");
        feature.edited = true;
        doc.sequence = "CCATGCAATAAGGG".into();
        assert_eq!(translate_cds(&doc, &feature).unwrap().amino_acids, b"MQ*");
        let (doc, mut feature) = cds("AATGAAATAA", false);
        feature
            .qualifiers
            .push(("codon_start".into(), Some("2".into())));
        let result = translate_cds(&doc, &feature).unwrap();
        assert_eq!((result.amino_acids, result.phase), (b"MK*".to_vec(), 1));
        let (doc, feature) = cds("TTATTTCAT", true);
        assert_eq!(translate_cds(&doc, &feature).unwrap().amino_acids, b"MK*");
        let (mut doc, mut feature) = cds("AAATAAATG", false);
        doc.circular = true;
        feature.parts = vec![Segment { start: 6, end: 9 }, Segment { start: 0, end: 6 }];
        assert_eq!(translate_cds(&doc, &feature).unwrap().amino_acids, b"MK*");
        let (doc, mut feature) = cds("TTAGGGTTTCAT", true);
        feature.parts = vec![Segment { start: 0, end: 3 }, Segment { start: 6, end: 12 }];
        assert_eq!(translate_cds(&doc, &feature).unwrap().amino_acids, b"MK*");
        let doc = Document::new_typed("RNA", "UUAUUUCAU", Molecule::Rna).unwrap();
        feature.parts = vec![Segment { start: 0, end: 9 }];
        assert_eq!(translate_cds(&doc, &feature).unwrap().amino_acids, b"MK*");
    }
    #[test]
    fn cds_initiation_trailing_bases_and_unsupported_metadata() {
        let (start_doc, start_feature) = cds("YTGCTGTAA", false);
        assert_eq!(
            translate_cds(&start_doc, &start_feature)
                .unwrap()
                .amino_acids,
            b"ML*"
        );
        let (doc, mut feature) = cds("GTGAAATAAC", false);
        feature
            .qualifiers
            .push(("transl_table".into(), Some("11".into())));
        let result = translate_cds(&doc, &feature).unwrap();
        assert_eq!(result.amino_acids, b"MK*");
        assert_eq!(result.trailing_bases, 1);
        for (key, value) in [
            ("transl_table", "7"),
            ("codon_start", "0"),
            ("transl_except", "(pos:2..4,aa:Sec)"),
            ("exception", "RNA editing"),
            ("pseudo", ""),
        ] {
            feature.qualifiers = vec![(key.into(), Some(value.into()))];
            assert!(translate_cds(&doc, &feature).is_err(), "{key}");
        }
    }

    #[test]
    fn native_reading_frame_uses_observed_phase_and_preserves_conflicts() {
        let (doc, mut feature) = cds("AATGAAATAA", false);
        feature
            .qualifiers
            .push(("snapgene_readingFrame".into(), Some("2".into())));
        let original = feature.clone();
        assert_eq!(translate_cds(&doc, &feature).unwrap().amino_acids, b"MK*");
        assert_eq!(feature, original);
        feature
            .qualifiers
            .push(("codon_start".into(), Some("1".into())));
        assert!(
            translate_cds(&doc, &feature)
                .unwrap_err()
                .to_string()
                .contains("conflicts with /codon_start")
        );
        feature.qualifiers.pop();
        feature.reverse = true;
        assert!(
            translate_cds(&doc, &feature)
                .unwrap_err()
                .to_string()
                .contains("conflicts with the CDS strand")
        );
        feature.reverse = false;
        feature
            .qualifiers
            .push(("snapgene_readingFrame".into(), Some("3".into())));
        assert!(
            translate_cds(&doc, &feature)
                .unwrap_err()
                .to_string()
                .contains("conflicting /snapgene_readingFrame")
        );
    }
}
