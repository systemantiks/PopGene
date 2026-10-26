//! Offline recognition-site analysis for a versioned, licensed enzyme catalogue.
//!
//! Coordinates are zero-based boundaries on the document's forward strand. Recognition
//! is not a prediction of digestion efficiency. The digest module adds explicit
//! methylation/condition assessment and geometric complete-digest products.
//! Catalogue provenance: Biopython1.88 / REBASE404 (2024), plus vendor EcoP15I.

mod digest;
pub use digest::*;

/// Cuts are signed boundary offsets from the start of the written 5′→3′ motif.
/// For example, EcoRI has offsets 1/5; SapI GCTCTTC(1/4) has offsets 8/11.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Enzyme {
    pub name: &'static str,
    pub recognition: &'static str,
    pub top_cut: i32,
    pub bottom_cut: i32,
    /// Unknown offsets use zero storage placeholders; callers must check this.
    pub cleavage_known: bool,
    /// Additional pair of cuts for enzymes that excise their recognition site.
    pub second_cut: Option<(i32, i32)>,
    pub rebase_id: Option<u32>,
    /// Reference catalogue value, not a guarantee for a particular supplier/buffer.
    pub catalogue_temperature_c: i16,
}
impl Enzyme {
    /// Shortest recognition alternative, for catalogue length filters. Individual
    /// matches carry their actual length in `RestrictionSite::recognition_length`.
    pub fn recognition_length(&self) -> usize {
        self.recognition.split('|').map(str::len).min().unwrap_or(0)
    }
}
// All entries in the pinned Biopython1.88 / REBASE404 snapshot plus EcoP15I.
// Original catalogue indices remain stable. Unknown cuts are explicitly marked.
include!("restriction_catalogue.rs");

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RestrictionSite {
    pub enzyme_index: usize,
    /// Leftmost motif position in forward document coordinates, including reverse sites.
    pub start: usize,
    pub recognition_length: usize,
    pub reverse: bool,
    /// Cleavage boundary on the forward document strand; None outside a linear molecule.
    pub top_cut: Option<usize>,
    /// Cleavage boundary on the opposite strand, in the same document coordinate system.
    pub bottom_cut: Option<usize>,
    pub second_top_cut: Option<usize>,
    pub second_bottom_cut: Option<usize>,
}

/// Keep output memory bounded even for repetitive chromosomes. Counts remain exact.
pub const MAX_SITES_PER_ENZYME: usize = 512;
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RestrictionAnalysis {
    /// Grouped by catalogue order, then recognition-site start. Capped per enzyme.
    pub sites: Vec<RestrictionSite>,
    /// Exact recognition-site counts, indexed in the same order as ENZYMES.
    pub counts: Vec<usize>,
    pub total_sites: usize,
    /// At least one enzyme has more positions than the retained per-enzyme limit.
    pub truncated: bool,
}

/// Search an ASCII DNA sequence. Lowercase is accepted. Any non-ACGT target symbol,
/// including N/R/Y, breaks a match rather than claiming a definite recognition site.
/// Circular motifs may cross the origin but must not exceed the molecule's length.
pub fn find_sites(sequence: &str, circular: bool) -> RestrictionAnalysis {
    find_sites_cancellable(sequence, circular, || false).expect("cancellation is disabled")
}

/// Same analysis with cooperative cancellation. Call from a worker for large documents.
/// The callback is checked before each enzyme and every 16,384 input symbols. Cancellation
/// returns None rather than exposing partial counts as if they were a complete analysis.
pub fn find_sites_cancellable(
    sequence: &str,
    circular: bool,
    should_cancel: impl FnMut() -> bool,
) -> Option<RestrictionAnalysis> {
    analyze_selected(
        sequence,
        circular,
        &(0..ENZYMES.len()).collect::<Vec<_>>(),
        MAX_SITES_PER_ENZYME,
        should_cancel,
    )
}

/// Analyze a selected enzyme set without scanning the rest of the catalogue.
pub fn find_sites_for(
    sequence: &str,
    circular: bool,
    enzyme_indices: &[usize],
) -> Result<RestrictionAnalysis, String> {
    if enzyme_indices.iter().any(|&i| i >= ENZYMES.len()) {
        return Err("Invalid enzyme index".into());
    }
    Ok(analyze_selected(
        sequence,
        circular,
        enzyme_indices,
        MAX_SITES_PER_ENZYME,
        || false,
    )
    .unwrap())
}

fn analyze_selected(
    sequence: &str,
    circular: bool,
    enzyme_indices: &[usize],
    site_limit: usize,
    mut should_cancel: impl FnMut() -> bool,
) -> Option<RestrictionAnalysis> {
    let sequence = sequence.as_bytes();
    let n = sequence.len();
    let mut result = RestrictionAnalysis {
        sites: Vec::new(),
        counts: vec![0; ENZYMES.len()],
        total_sites: 0,
        truncated: false,
    };
    let mut groups = std::collections::BTreeMap::<&str, Vec<usize>>::new();
    for &index in enzyme_indices {
        let indices = groups.entry(ENZYMES[index].recognition).or_default();
        if !indices.contains(&index) {
            indices.push(index);
        }
    }
    for (recognition, indices) in groups {
        if should_cancel() {
            return None;
        }
        for motif in recognition.split('|') {
            let length = motif.len();
            if length > n {
                continue;
            }
            // Shift-And: each bit encodes a matched prefix of the motif. The pattern
            // masks contain all canonical alternatives for degenerate IUPAC symbols.
            let forward = pattern_masks(motif.as_bytes(), false);
            let reverse = pattern_masks(motif.as_bytes(), true);
            let palindrome = forward == reverse;
            let done = 1u64 << (length - 1);
            let mut forward_state = 0;
            let mut reverse_state = 0;
            let scan_length = n + if circular { length - 1 } else { 0 };
            for i in 0..scan_length {
                if i & 0x3fff == 0 && should_cancel() {
                    return None;
                }
                let base = sequence[if i < n { i } else { i - n }];
                let Some(base) = canonical_index(base) else {
                    forward_state = 0;
                    reverse_state = 0;
                    continue;
                };
                forward_state = ((forward_state << 1) | 1) & forward[base];
                if !palindrome {
                    reverse_state = ((reverse_state << 1) | 1) & reverse[base];
                }
                let found_forward = forward_state & done != 0;
                let found_reverse = !palindrome && reverse_state & done != 0;
                if !found_forward && !found_reverse {
                    continue;
                }
                let start = i + 1 - length;
                for is_reverse in [false, true] {
                    if (!is_reverse && !found_forward) || (is_reverse && !found_reverse) {
                        continue;
                    }
                    for &enzyme_index in &indices {
                        let enzyme = &ENZYMES[enzyme_index];
                        let count = &mut result.counts[enzyme_index];
                        *count += 1;
                        result.total_sites += 1;
                        if *count > site_limit {
                            result.truncated = true;
                            continue;
                        }
                        let (top, bottom) = if is_reverse {
                            (
                                length as i32 - enzyme.bottom_cut,
                                length as i32 - enzyme.top_cut,
                            )
                        } else {
                            (enzyme.top_cut, enzyme.bottom_cut)
                        };
                        let secondary = enzyme.second_cut.map(|(top, bottom)| {
                            if is_reverse {
                                (length as i32 - bottom, length as i32 - top)
                            } else {
                                (top, bottom)
                            }
                        });
                        result.sites.push(RestrictionSite {
                            enzyme_index,
                            start,
                            recognition_length: length,
                            reverse: is_reverse,
                            top_cut: enzyme
                                .cleavage_known
                                .then(|| boundary(start, top, n, circular))
                                .flatten(),
                            bottom_cut: enzyme
                                .cleavage_known
                                .then(|| boundary(start, bottom, n, circular))
                                .flatten(),
                            second_top_cut: secondary
                                .and_then(|(top, _)| boundary(start, top, n, circular)),
                            second_bottom_cut: secondary
                                .and_then(|(_, bottom)| boundary(start, bottom, n, circular)),
                        });
                    }
                }
            }
        }
    }
    result
        .sites
        .sort_by_key(|site| (site.enzyme_index, site.start, site.reverse));
    Some(result)
}

fn boundary(start: usize, offset: i32, length: usize, circular: bool) -> Option<usize> {
    let cut = start as i128 + i128::from(offset);
    if circular {
        Some(cut.rem_euclid(length as i128) as usize)
    } else if (0..=length as i128).contains(&cut) {
        Some(cut as usize)
    } else {
        None
    }
}
fn canonical_index(base: u8) -> Option<usize> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}
fn iupac_mask(base: u8) -> u8 {
    match base {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' => 8,
        b'R' => 5,
        b'Y' => 10,
        b'S' => 6,
        b'W' => 9,
        b'K' => 12,
        b'M' => 3,
        b'B' => 14,
        b'D' => 13,
        b'H' => 11,
        b'V' => 7,
        b'N' => 15,
        _ => 0,
    }
}
fn complement_mask(mask: u8) -> u8 {
    ((mask & 1) << 3) | ((mask & 2) << 1) | ((mask & 4) >> 1) | ((mask & 8) >> 3)
}
fn pattern_masks(pattern: &[u8], reverse: bool) -> [u64; 4] {
    let mut masks = [0; 4];
    for position in 0..pattern.len() {
        let mask = if reverse {
            complement_mask(iupac_mask(pattern[pattern.len() - 1 - position]))
        } else {
            iupac_mask(pattern[position])
        };
        for (base, bits) in masks.iter_mut().enumerate() {
            if mask & (1 << base) != 0 {
                *bits |= 1 << position;
            }
        }
    }
    masks
}

#[cfg(test)]
mod tests {
    use super::*;
    fn sites(sequence: &str, circular: bool, name: &str) -> Vec<RestrictionSite> {
        find_sites(sequence, circular)
            .sites
            .into_iter()
            .filter(|s| ENZYMES[s.enzyme_index].name == name)
            .collect()
    }
    #[test]
    fn catalogue_has_valid_unique_names_and_known_palindromic_cut_offsets() {
        let mut names = std::collections::HashSet::new();
        for enzyme in ENZYMES {
            assert!(names.insert(enzyme.name));
            assert!((1..=64).contains(&enzyme.recognition.len()));
            assert!(
                enzyme
                    .recognition
                    .bytes()
                    .all(|b| b == b'|' || iupac_mask(b) != 0)
            );
            if enzyme.cleavage_known
                && enzyme.second_cut.is_none()
                && (0..=enzyme.recognition.len() as i32).contains(&enzyme.top_cut)
                && (0..=enzyme.recognition.len() as i32).contains(&enzyme.bottom_cut)
                && pattern_masks(enzyme.recognition.as_bytes(), false)
                    == pattern_masks(enzyme.recognition.as_bytes(), true)
            {
                assert_eq!(
                    enzyme.top_cut + enzyme.bottom_cut,
                    enzyme.recognition.len() as i32
                );
            }
        }
    }
    #[test]
    fn palindromic_ecori_is_reported_once_with_correct_boundaries() {
        let found = sites("AAAGAATTCTTT", false, "EcoRI");
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].start, 3);
        assert!(!found[0].reverse);
        assert_eq!((found[0].top_cut, found[0].bottom_cut), (Some(4), Some(8)));
        assert_eq!(sites("aaagaattcttt", false, "EcoRI"), found);
    }
    #[test]
    fn same_motif_can_have_distinct_cleavage_positions() {
        let a = sites("AAAGGTACCTTT", false, "Acc65I");
        let k = sites("AAAGGTACCTTT", false, "KpnI");
        assert_eq!(a[0].start, k[0].start);
        assert_eq!((a[0].top_cut, a[0].bottom_cut), (Some(4), Some(8)));
        assert_eq!((k[0].top_cut, k[0].bottom_cut), (Some(8), Some(4)));
    }
    #[test]
    fn four_base_cutter_preserves_terminal_boundaries() {
        let found = sites("GATC", false, "DpnII");
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].start, 0);
        assert_eq!((found[0].top_cut, found[0].bottom_cut), (Some(0), Some(4)));
        let circular = sites("GATC", true, "DpnII");
        assert_eq!(
            (circular[0].top_cut, circular[0].bottom_cut),
            (Some(0), Some(0))
        );
    }
    #[test]
    fn type_iis_forward_and_reverse_offsets_are_mirrored_and_swapped() {
        let forward = sites("AAAAGCTCTTCAAAAAAAA", false, "SapI");
        let reverse = sites("AAAAAAAAGAAGAGCAAAA", false, "SapI");
        assert_eq!((forward[0].start, forward[0].reverse), (4, false));
        assert_eq!(
            (forward[0].top_cut, forward[0].bottom_cut),
            (Some(12), Some(15))
        );
        assert_eq!((reverse[0].start, reverse[0].reverse), (8, true));
        assert_eq!(
            (reverse[0].top_cut, reverse[0].bottom_cut),
            (Some(4), Some(7))
        );
    }
    #[test]
    fn linear_edge_sites_remain_visible_but_impossible_cuts_are_none() {
        let forward = sites("GCTCTTC", false, "SapI");
        let reverse = sites("GAAGAGC", false, "SapI");
        assert_eq!((forward[0].top_cut, forward[0].bottom_cut), (None, None));
        assert_eq!((reverse[0].top_cut, reverse[0].bottom_cut), (None, None));
        let partial = sites("GCTCTTCA", false, "SapI");
        assert_eq!((partial[0].top_cut, partial[0].bottom_cut), (Some(8), None));
    }
    #[test]
    fn circular_origin_and_outside_motif_cuts_wrap_at_boundaries() {
        assert!(sites("AATTCAAAAG", false, "EcoRI").is_empty());
        let origin = sites("AATTCAAAAG", true, "EcoRI");
        assert_eq!(origin.len(), 1);
        assert_eq!(
            (origin[0].start, origin[0].top_cut, origin[0].bottom_cut),
            (9, Some(0), Some(4))
        );
        let reverse = sites("GAAGAGC", true, "SapI");
        assert_eq!(
            (reverse[0].top_cut, reverse[0].bottom_cut),
            (Some(3), Some(6))
        );
        assert!(sites("AAT", true, "EcoRI").is_empty());
    }
    #[test]
    fn motif_degeneracy_accepts_canonical_bases_but_never_unknown_targets() {
        for sequence in ["ATGCGCAC", "GTGCGCAT"] {
            assert_eq!(sites(sequence, false, "FspAI").len(), 1);
        }
        for sequence in ["CTGCGCAC", "ATGCGCAA", "RTGCGCAY", "NTGCGCAC"] {
            assert!(sites(sequence, false, "FspAI").is_empty());
        }
        assert_eq!(sites("GGTGACC", false, "BstEII").len(), 1);
        assert!(sites("GGTNACC", false, "BstEII").is_empty());
        assert!(sites("GAANTTC", false, "EcoRI").is_empty());
    }
    #[test]
    fn overlapping_sites_are_counted_and_output_cap_does_not_hide_unique_cutters() {
        let sequence = format!("{}GGATCC", "GCGC".repeat(2_000));
        let result = find_sites(&sequence, false);
        let bss = ENZYMES.iter().position(|e| e.name == "BssHII").unwrap();
        let bam = ENZYMES.iter().position(|e| e.name == "BamHI").unwrap();
        assert_eq!(result.counts[bss], 3_998);
        assert_eq!(
            result
                .sites
                .iter()
                .filter(|s| s.enzyme_index == bss)
                .count(),
            MAX_SITES_PER_ENZYME
        );
        assert_eq!(result.counts[bam], 1);
        assert!(result.sites.iter().any(|s| s.enzyme_index == bam));
        assert!(result.truncated);
        assert_eq!(result.total_sites, result.counts.iter().sum::<usize>());
    }
    #[test]
    fn taqii_uses_revised_specificity_and_ecop15i_has_documented_offsets() {
        assert!(sites("CACCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", false, "TaqII").is_empty());
        let taq = sites("GACCGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", false, "TaqII");
        assert_eq!((taq[0].top_cut, taq[0].bottom_cut), (Some(17), Some(15)));
        let eco = sites("CAGCAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", false, "EcoP15I");
        assert_eq!((eco[0].top_cut, eco[0].bottom_cut), (Some(31), Some(33)));
    }
    #[test]
    fn cancellation_never_returns_partial_counts() {
        let mut calls = 0;
        let sequence = "ACGT".repeat(20_000);
        assert!(
            find_sites_cancellable(&sequence, true, || {
                calls += 1;
                calls == 3
            })
            .is_none()
        );
        assert_eq!(calls, 3);
        assert_eq!(find_sites("", true).total_sites, 0);
    }
}
