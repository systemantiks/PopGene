//! Bounded start-to-stop ORF predictions; these are not asserted annotations.
use super::{
    GeneticCode, all_codons, complement, genetic_codes, is_start_codon, is_stop_codon,
    translate_codon,
};
use crate::{Document, Molecule, Segment};
use std::collections::VecDeque;

#[derive(Clone, Copy, Debug)]
pub struct OrfOptions {
    pub genetic_code: GeneticCode,
    /// Residues before the stop codon; the stop is not counted.
    pub min_codons: usize,
    pub alternative_starts: bool,
    pub include_partial: bool,
    pub max_results: usize,
}
impl Default for OrfOptions {
    fn default() -> Self {
        Self {
            genetic_code: GeneticCode::STANDARD,
            min_codons: 30,
            alternative_starts: false,
            include_partial: false,
            max_results: 10_000,
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Orf {
    /// First base in coding direction, in forward genomic coordinates.
    pub start: usize,
    pub reverse: bool,
    /// Includes a terminal stop, when present. Never exceeds one sequence lap.
    pub nucleotide_length: usize,
    pub amino_acid_length: usize,
    pub complete: bool,
    pub genetic_code: u16,
    pub context_dependent_stop: bool,
}
impl Orf {
    pub fn segments(&self, sequence_length: usize) -> Vec<Segment> {
        if sequence_length == 0 || self.nucleotide_length == 0 {
            return vec![];
        }
        let start = if self.reverse {
            (self.start + sequence_length + 1 - self.nucleotide_length) % sequence_length
        } else {
            self.start
        };
        let end = start + self.nucleotide_length;
        if end <= sequence_length {
            vec![Segment { start, end }]
        } else {
            vec![
                Segment {
                    start,
                    end: sequence_length,
                },
                Segment {
                    start: 0,
                    end: end - sequence_length,
                },
            ]
        }
    }
    pub fn frame(&self, sequence_length: usize) -> i8 {
        let coding_start = if self.reverse {
            sequence_length - 1 - self.start
        } else {
            self.start
        };
        (coding_start % 3 + 1) as i8 * if self.reverse { -1 } else { 1 }
    }
    pub fn protein(&self, sequence: &str) -> String {
        let Some(code) = GeneticCode::from_id(self.genetic_code) else {
            return String::new();
        };
        let n = sequence.len();
        if n == 0 || self.start >= n || self.nucleotide_length > n {
            return String::new();
        }
        (0..self.nucleotide_length / 3)
            .map(|index| {
                if index == 0 {
                    return 'M';
                }
                if self.complete && index + 1 == self.nucleotide_length / 3 {
                    return '*';
                }
                let codon = std::array::from_fn(|i| {
                    let offset = index * 3 + i;
                    if self.reverse {
                        complement(sequence.as_bytes()[(self.start + n - offset % n) % n])
                    } else {
                        sequence.as_bytes()[(self.start + offset) % n]
                    }
                });
                translate_codon(codon, code) as char
            })
            .collect()
    }
}
#[derive(Clone, Debug, Default)]
pub struct OrfSearch {
    pub orfs: Vec<Orf>,
    pub truncated: bool,
}

/// Longest start-to-stop candidate in each uninterrupted interval on each strand.
/// All three phases are scanned. Circular candidates cross the origin, but never
/// traverse more than one circumference. Ambiguous possible stops break a search
/// interval instead of asserting a continuous open frame through unknown bases.
pub fn find_orfs(document: &Document, options: OrfOptions) -> OrfSearch {
    let n = document.sequence.len();
    let mut result = OrfSearch::default();
    if n < 3 || document.molecule == Molecule::Protein {
        return result;
    }
    let code = options.genetic_code;
    let max_results = options.max_results.min(100_000);
    let emit = |result: &mut OrfSearch, start: usize, end: usize, reverse: bool, complete: bool| {
        if start >= n || end <= start || end - start > n {
            return;
        }
        let residues = (end - start) / 3 - usize::from(complete);
        if residues < options.min_codons {
            return;
        }
        if result.orfs.len() >= max_results {
            result.truncated = true;
            return;
        }
        result.orfs.push(Orf {
            start: if reverse { n - 1 - start } else { start },
            reverse,
            nucleotide_length: end - start,
            amino_acid_length: residues,
            complete,
            genetic_code: code.id(),
            context_dependent_stop: complete && code.has_context_dependent_stops(),
        });
    };
    for reverse in [false, true] {
        for phase in 0..3 {
            let limit = if document.circular {
                n.saturating_mul(2)
            } else {
                n
            };
            let mut starts = VecDeque::new();
            let mut last_end = phase;
            let mut partial_emitted = false;
            for offset in (phase..limit.saturating_sub(2)).step_by(3) {
                while starts.front().is_some_and(|start| offset + 3 - *start > n) {
                    let start = starts.pop_front().unwrap();
                    if options.include_partial && !partial_emitted {
                        partial_emitted = true;
                        emit(&mut result, start, offset, reverse, false);
                    }
                }
                let codon = std::array::from_fn(|i| {
                    let p = (offset + i) % n;
                    if reverse {
                        complement(document.sequence.as_bytes()[n - 1 - p])
                    } else {
                        document.sequence.as_bytes()[p]
                    }
                });
                let stop = is_stop_codon(codon, code);
                let possible_stop = !all_codons(codon, |index| {
                    genetic_codes::CODES[code.0].starts[index] != b'*'
                });
                if stop {
                    if let Some(start) = starts.pop_front() {
                        emit(&mut result, start, offset + 3, reverse, true);
                    }
                    starts.clear();
                } else if possible_stop {
                    // Unknown sequence cannot establish a complete ORF.
                    starts.clear();
                } else if offset < n && is_start_codon(codon, code, options.alternative_starts) {
                    starts.push_back(offset);
                }
                last_end = offset + 3;
                if result.truncated {
                    break;
                }
            }
            if options.include_partial
                && !document.circular
                && let Some(start) = starts.front()
            {
                emit(&mut result, *start, last_end, reverse, false);
            }
            if result.truncated {
                break;
            }
        }
        if result.truncated {
            break;
        }
    }
    result
        .orfs
        .sort_by_key(|orf| (orf.start, orf.reverse, orf.nucleotide_length));
    result.orfs.dedup();
    result
}
