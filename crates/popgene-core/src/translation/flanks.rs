//! Hypothetical upstream/downstream translations, never imported CDS claims.
use super::{
    Document, Feature, FrameCodon, GeneticCode, TranslationError, complement, is_stop_codon,
    translate_cds, translate_codon,
};
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FlankingCodon {
    pub codon: FrameCodon,
    pub upstream: bool,
}
/// Extend the annotated reading phase up to a stop or `max_codons` (capped at
/// 300) on each side. Circular extension stops before revisiting CDS bases.
/// Joined CDS parts remain joined; only the outer biological ends are extended.
pub fn cds_flanking_codons(
    document: &Document,
    feature: &Feature,
    max_codons: usize,
) -> Result<Vec<FlankingCodon>, TranslationError> {
    let translation = translate_cds(document, feature)?;
    let code = GeneticCode::from_id(translation.genetic_code).unwrap();
    let n = document.sequence.len();
    let length = feature.len();
    if n < 3 || length == 0 {
        return Ok(vec![]);
    }
    let mut offsets = Vec::with_capacity(feature.parts.len());
    let mut total = 0;
    for part in &feature.parts {
        offsets.push(total);
        total += part.end - part.start;
    }
    let inside = |coding_offset: usize| {
        let offset = if feature.reverse {
            length - 1 - coding_offset
        } else {
            coding_offset
        };
        let part_index = offsets.partition_point(|start| *start <= offset) - 1;
        feature.parts[part_index].start + offset - offsets[part_index]
    };
    let first = inside(0);
    let last = inside(length - 1);
    let coordinate = |offset: isize| -> Option<usize> {
        if (0..length as isize).contains(&offset) {
            return Some(inside(offset as usize));
        }
        let raw = if offset < 0 {
            first as isize + offset * if feature.reverse { -1 } else { 1 }
        } else {
            last as isize + (offset - length as isize + 1) * if feature.reverse { -1 } else { 1 }
        };
        if !document.circular && !(0..n as isize).contains(&raw) {
            return None;
        }
        let p = raw.rem_euclid(n as isize) as usize;
        if feature
            .parts
            .iter()
            .any(|part| (part.start..part.end).contains(&p))
        {
            return None;
        }
        Some(p)
    };
    let mut result = vec![];
    for upstream in [true, false] {
        let first_offset = if upstream {
            translation.phase as isize - 3
        } else {
            (translation.phase + (length.saturating_sub(translation.phase) / 3) * 3) as isize
        };
        for index in 0..max_codons.min(300) {
            let offset = first_offset + index as isize * if upstream { -3 } else { 3 };
            if offset.unsigned_abs() > length + n + 3 {
                break;
            }
            let Some(positions) = (0..3)
                .map(|i| coordinate(offset + i))
                .collect::<Option<Vec<_>>>()
            else {
                break;
            };
            let positions: [usize; 3] = positions.try_into().unwrap();
            if positions[0] == positions[1]
                || positions[0] == positions[2]
                || positions[1] == positions[2]
            {
                break;
            }
            let bases = positions.map(|p| {
                if feature.reverse {
                    complement(document.sequence.as_bytes()[p])
                } else {
                    document.sequence.as_bytes()[p]
                }
            });
            let amino_acid = translate_codon(bases, code);
            let center = positions[1];
            let low = *positions.iter().min().unwrap();
            let high = *positions.iter().max().unwrap();
            let (start, end) = if high - low == 2 {
                (low, high + 1)
            } else {
                (center, center + 1)
            };
            result.push(FlankingCodon {
                codon: FrameCodon {
                    center,
                    start,
                    end,
                    amino_acid,
                    positions,
                },
                upstream,
            });
            if is_stop_codon(bases, code) {
                break;
            }
        }
    }
    Ok(result)
}
