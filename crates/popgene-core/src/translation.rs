//! Bounded, deterministic nucleotide translation.
//!
//! Codon assignments and initiators: NCBI Genetic Codes, tables 1, 2 and 11,
//! https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi (2024-09-23 revision).
//! CDS offsets follow INSDC codon_start. Exceptions are refused explicitly.
//! Frame translation does not infer ORFs or translate across a circular seam.
use crate::{Document, Feature, Molecule};
use std::{fmt, ops::Range};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GeneticCode {
    Standard,
    VertebrateMitochondrial,
    Bacterial,
}
impl GeneticCode {
    pub fn from_id(id: u16) -> Option<Self> {
        match id {
            1 => Some(Self::Standard),
            2 => Some(Self::VertebrateMitochondrial),
            11 => Some(Self::Bacterial),
            _ => None,
        }
    }
    pub fn id(self) -> u16 {
        match self {
            Self::Standard => 1,
            Self::VertebrateMitochondrial => 2,
            Self::Bacterial => 11,
        }
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

/// NCBI's TCAG codon order. Table 11 has the same internal codons as table 1.
const STANDARD: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
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
    if code == GeneticCode::VertebrateMitochondrial {
        match index {
            14 => return b'W',      // TGA
            34 => return b'M',      // ATA
            46 | 47 => return b'*', // AGA / AGG
            _ => {}
        }
    }
    STANDARD[index]
}

/// Ambiguity resolves only when every concrete codon has the same assignment.
/// Otherwise X is returned, including stop/amino-acid mixtures and invalid input.
pub fn translate_codon(codon: [u8; 3], code: GeneticCode) -> u8 {
    let masks = codon.map(mask);
    if masks.contains(&0) {
        return b'X';
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
    pub start: usize,
    pub end: usize,
    pub amino_acid: u8,
}

/// Codons centered in a visible forward-coordinate range. Reverse frame -1
/// starts at the final nucleotide, -2 one base before it, -3 two before it.
/// Complete codons can straddle a visual row boundary, but never sequence ends.
pub fn frame_codons(sequence: &str, frame: i8, range: Range<usize>) -> Vec<FrameCodon> {
    if !matches!(frame, -3..=-1 | 1..=3) {
        return vec![];
    }
    let bases = sequence.as_bytes();
    let phase = usize::from(frame.unsigned_abs() - 1);
    let mut result = vec![];
    for center in range.start.min(bases.len())..range.end.min(bases.len()) {
        let coding_center = if frame < 0 {
            bases.len() - 1 - center
        } else {
            center
        };
        if coding_center < phase + 1
            || (coding_center - phase) % 3 != 1
            || center == 0
            || center + 1 >= bases.len()
        {
            continue;
        }
        let codon = if frame < 0 {
            [
                complement(bases[center + 1]),
                complement(bases[center]),
                complement(bases[center - 1]),
            ]
        } else {
            [bases[center - 1], bases[center], bases[center + 1]]
        };
        result.push(FrameCodon {
            center,
            start: center - 1,
            end: center + 2,
            amino_acid: translate_codon(codon, GeneticCode::Standard),
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
fn initiator(codon: [u8; 3], code: GeneticCode) -> bool {
    let masks = codon.map(mask);
    if masks.contains(&0) {
        return false;
    }
    let allowed: &[usize] = match code {
        GeneticCode::Standard => &[3, 19, 35],
        GeneticCode::VertebrateMitochondrial => &[32, 33, 34, 35, 51],
        GeneticCode::Bacterial => &[3, 19, 32, 33, 34, 35, 51],
    };
    for first in 0..4 {
        for second in 0..4 {
            for third in 0..4 {
                if masks[0] & (1 << first) != 0
                    && masks[1] & (1 << second) != 0
                    && masks[2] & (1 << third) != 0
                    && !allowed.contains(&(first * 16 + second * 4 + third))
                {
                    return false;
                }
            }
        }
    }
    true
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
    for (key, _) in &feature.qualifiers {
        if matches!(
            key.as_str(),
            "transl_except"
                | "exception"
                | "ribosomal_slippage"
                | "codon"
                | "pseudo"
                | "pseudogene"
        ) {
            return Err(TranslationError(format!(
                "/{key} requires interpretation; calculated CDS translation omitted"
            )));
        }
    }
    let phase = cds_phase(feature)?;
    let table = qualifier(feature, "transl_table")?
        .unwrap_or("1")
        .parse::<u16>()
        .map_err(|_| TranslationError("invalid /transl_table".into()))?;
    let code = GeneticCode::from_id(table).ok_or_else(|| {
        TranslationError(format!(
            "genetic code table {table} is unsupported (supported: 1, 2, 11)"
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
    for (index, &codon) in coding.as_chunks::<3>().0.iter().enumerate() {
        let aa = if index == 0 && phase == 0 && initiator(codon, code) {
            b'M'
        } else {
            translate_codon(codon, code)
        };
        amino_acids.push(aa);
    }
    Ok(CdsTranslation {
        amino_acids,
        phase,
        genetic_code: code.id(),
        trailing_bases: coding.len() % 3,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Segment;
    fn cds(sequence: &str, reverse: bool) -> (Document, Feature) {
        let doc = Document::new("translation", sequence).unwrap();
        let feature = Feature {
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
                    translate_codon(codon.as_bytes().try_into().unwrap(), GeneticCode::Standard),
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
                translate_codon(codon.as_bytes().try_into().unwrap(), GeneticCode::Standard),
                aa
            );
        }
        for (codon, aa) in [("TGA", b'W'), ("ATA", b'M'), ("AGA", b'*'), ("AGG", b'*')] {
            assert_eq!(
                translate_codon(
                    codon.as_bytes().try_into().unwrap(),
                    GeneticCode::VertebrateMitochondrial
                ),
                aa
            );
        }
        assert_eq!(GeneticCode::from_id(4), None);
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
                amino_acid: b'K'
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
            ("transl_table", "4"),
            ("codon_start", "0"),
            ("transl_except", "(pos:1..3,aa:Sec)"),
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
