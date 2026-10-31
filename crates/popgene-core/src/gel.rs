//! Agarose gel display model for double-stranded DNA products.
//! Efficient separation ranges for standard agarose: Ausubel et al. (1995)
//! Current Protocols in Molecular Biology, as tabulated (retrieved 2026-10-30) at
//! https://www.thermofisher.cn/cn/en/home/references/ambion-tech-support/rna-electrophoresis-markers/general-articles/gel-electrophoresis-tables.html
//! Positions place bands for comparison between lanes of the same gel; they are not
//! measured mobilities. Voltage, buffer, run time, stain, load and DNA topology are
//! not modeled. Equimolar lanes carry mass in proportion to fragment length.
use crate::restriction::DigestFragment;

/// (agarose %, smallest bp, largest bp) efficiently separated as linear dsDNA.
const SEPARATION: &[(f64, f64, f64)] = &[
    (0.5, 1_000., 30_000.),
    (0.7, 800., 12_000.),
    (1.0, 500., 10_000.),
    (1.2, 400., 7_000.),
    (1.5, 200., 3_000.),
];
/// Run fractions from the well at which the largest and smallest tabulated
/// lengths are drawn. Lengths beyond the range are compressed towards the ends.
const RANGE_TOP: f64 = 0.15;
const RANGE_BOTTOM: f64 = 0.85;
/// Fragments closer than this fraction of the run are drawn as one band.
pub const RESOLUTION: f64 = 0.012;
pub const MAX_LANES: usize = 12;
const MAX_COMIGRATION_WARNINGS: usize = 8;
pub const MODEL_NOTE: &str = "Band positions follow a logistic in log10(length) through the tabulated separation limits of standard agarose. They compare lanes on one gel; they are not measured mobilities. Voltage, buffer, run time, stain, load and DNA topology are not modeled.";

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Ladder {
    pub name: &'static str,
    /// Fragment lengths in bp, largest first.
    pub sizes: &'static [usize],
    /// Manufacturer-designated orientation bands.
    pub reference: &'static [usize],
    /// A complete digest has equal moles per fragment; purified mixtures do not.
    pub equimolar: bool,
    pub note: &'static str,
    pub source: &'static str,
}

/// Fragment lists are quoted from the manufacturer's product manuals.
pub const LADDERS: &[Ladder] = &[
    Ladder {
        name: "GeneRuler 1 kb",
        sizes: &[
            10_000, 8_000, 6_000, 5_000, 4_000, 3_500, 3_000, 2_500, 2_000, 1_500, 1_000, 750, 500,
            250,
        ],
        reference: &[6_000, 3_000, 1_000],
        equimolar: false,
        note: "Fourteen purified fragments. Band masses other than the three reference bands are not modeled.",
        source: "https://documents.thermofisher.com/TFS-Assets/LSG/manuals/MAN0013004_GeneRuler_1kb_DNALadder_250ug_UG.pdf",
    },
    Ladder {
        name: "λ DNA/HindIII",
        sizes: &[23_130, 9_416, 6_557, 4_361, 2_322, 2_027, 564, 125],
        reference: &[],
        equimolar: true,
        note: "Unless heated before loading, the cohesive 23,130 and 4,361 bp fragments can anneal into an additional 27,491 bp band.",
        source: "https://documents.thermofisher.com/TFS-Assets/LSG/manuals/MAN0012986_Lambda_DNA_HindIII_Marker2_250ug_UG.pdf",
    },
];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LaneFragment {
    pub length: usize,
    pub circular: bool,
    pub reference: bool,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Lane {
    pub label: String,
    pub fragments: Vec<LaneFragment>,
    pub equimolar: bool,
}
impl Lane {
    pub fn ladder(ladder: &Ladder) -> Self {
        Self {
            label: ladder.name.into(),
            fragments: ladder
                .sizes
                .iter()
                .map(|&length| LaneFragment {
                    length,
                    circular: false,
                    reference: ladder.reference.contains(&length),
                })
                .collect(),
            equimolar: ladder.equimolar,
        }
    }
    /// A user-entered ladder has unknown band masses, so bands are drawn uniformly.
    pub fn custom_ladder(label: impl Into<String>, sizes: &[usize]) -> Self {
        Self {
            label: label.into(),
            fragments: sizes
                .iter()
                .map(|&length| LaneFragment {
                    length,
                    circular: false,
                    reference: false,
                })
                .collect(),
            equimolar: false,
        }
    }
    /// Complete-digest products are equimolar. Sizes are top-strand lengths.
    pub fn digest(label: impl Into<String>, fragments: &[DigestFragment]) -> Self {
        Self {
            label: label.into(),
            fragments: fragments
                .iter()
                .map(|f| LaneFragment {
                    length: f.sequence.len(),
                    circular: f.circular,
                    reference: false,
                })
                .collect(),
            equimolar: true,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Band {
    /// Lengths of every fragment drawn in this band, largest first.
    pub lengths: Vec<usize>,
    /// Run fractions from the well of the slowest and fastest member.
    pub top: f64,
    pub bottom: f64,
    /// Mass relative to the lane's most intense band; uniform for unknown masses.
    pub intensity: f64,
    pub circular: bool,
    pub reference: bool,
}
#[derive(Clone, Debug, PartialEq)]
pub struct GelLane {
    pub label: String,
    pub bands: Vec<Band>,
}
#[derive(Clone, Debug, PartialEq)]
pub struct Gel {
    pub agarose_percent: f64,
    /// Efficient separation range for this gel, in bp.
    pub separation: (usize, usize),
    pub lanes: Vec<GelLane>,
    pub warnings: Vec<String>,
}

/// Bounds between tabulated concentrations are interpolated geometrically.
pub fn separation_range(percent: f64) -> Result<(f64, f64), String> {
    let (first, last) = (SEPARATION[0].0, SEPARATION[SEPARATION.len() - 1].0);
    if !percent.is_finite() || !(first..=last).contains(&percent) {
        return Err(format!(
            "Agarose must be between {first}% and {last}%; separation outside the tabulated standard-agarose range is not modeled"
        ));
    }
    let i = SEPARATION
        .windows(2)
        .position(|w| percent <= w[1].0)
        .unwrap_or(0);
    let (a, b) = (SEPARATION[i], SEPARATION[i + 1]);
    let t = (percent - a.0) / (b.0 - a.0);
    let between = |x: f64, y: f64| 10f64.powf(x.log10() + t * (y.log10() - x.log10()));
    Ok((between(a.1, b.1), between(a.2, b.2)))
}

/// Logistic in log10(length) that passes through the tabulated limits.
fn position(length: usize, small: f64, large: f64) -> f64 {
    let (lo, hi) = (small.log10(), large.log10());
    let (at_small, at_large) = ((1. / RANGE_BOTTOM - 1.).ln(), (1. / RANGE_TOP - 1.).ln());
    let slope = (at_large - at_small) / (hi - lo);
    1. / (1. + (at_small + slope * ((length as f64).log10() - lo)).exp())
}

/// Displayed run fraction from the well (0) towards the gel's end (1).
pub fn migration(length: usize, percent: f64) -> Result<f64, String> {
    if length == 0 {
        return Err("A fragment must contain at least one base pair".into());
    }
    let (small, large) = separation_range(percent)?;
    Ok(position(length, small, large))
}

pub fn format_bp(length: usize) -> String {
    let digits = length.to_string();
    let mut output = String::new();
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i).is_multiple_of(3) {
            output.push(',');
        }
        output.push(c);
    }
    output
}

fn comigration(label: &str, band: &Band) -> String {
    let (largest, smallest) = (band.lengths[0], band.lengths[band.lengths.len() - 1]);
    if largest == smallest {
        format!(
            "{label}: {} × {} bp fragments run as one band",
            band.lengths.len(),
            format_bp(largest)
        )
    } else {
        format!(
            "{label}: {} fragments of {}–{} bp may not resolve into separate bands",
            band.lengths.len(),
            format_bp(smallest),
            format_bp(largest)
        )
    }
}

/// Place every lane on one gel. Bands within `RESOLUTION` of each other merge
/// and are reported; circular molecules are drawn at their linear length.
pub fn simulate(lanes: &[Lane], percent: f64) -> Result<Gel, String> {
    if lanes.is_empty() || lanes.len() > MAX_LANES {
        return Err(format!("A gel needs between 1 and {MAX_LANES} lanes"));
    }
    let (small, large) = separation_range(percent)?;
    let mut warnings = vec![];
    let mut output = vec![];
    for lane in lanes {
        if lane.fragments.is_empty() {
            return Err(format!("{} has no fragments to load", lane.label));
        }
        if lane.fragments.iter().any(|f| f.length == 0) {
            return Err(format!("{} contains a zero-length fragment", lane.label));
        }
        let mut placed = lane
            .fragments
            .iter()
            .map(|f| (position(f.length, small, large), f))
            .collect::<Vec<_>>();
        placed.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut bands: Vec<Band> = vec![];
        for (p, f) in placed {
            let mass = if lane.equimolar { f.length as f64 } else { 1. };
            match bands.last_mut() {
                Some(band) if band.circular == f.circular && p - band.bottom < RESOLUTION => {
                    band.lengths.push(f.length);
                    band.bottom = p;
                    band.intensity += mass;
                    band.reference |= f.reference;
                }
                _ => bands.push(Band {
                    lengths: vec![f.length],
                    top: p,
                    bottom: p,
                    intensity: mass,
                    circular: f.circular,
                    reference: f.reference,
                }),
            }
        }
        let strongest = bands.iter().map(|b| b.intensity).fold(0., f64::max);
        for band in &mut bands {
            band.intensity /= strongest;
        }
        let outside = lane
            .fragments
            .iter()
            .filter(|f| !f.circular && !(small..=large).contains(&(f.length as f64)))
            .count();
        if outside > 0 {
            warnings.push(format!(
                "{}: {outside} fragment{} outside the {}–{} bp separation range of {percent}% agarose; positions there are compressed and poorly resolved",
                lane.label,
                if outside == 1 { " lies" } else { "s lie" },
                format_bp(small.round() as usize),
                format_bp(large.round() as usize)
            ));
        }
        if bands.iter().any(|b| b.circular) {
            warnings.push(format!(
                "{}: circular DNA migrates according to its topology (supercoiled, nicked or relaxed), not its length; it is drawn at its linear length for reference only",
                lane.label
            ));
        }
        let merged = bands
            .iter()
            .filter(|b| b.lengths.len() > 1)
            .collect::<Vec<_>>();
        warnings.extend(
            merged
                .iter()
                .take(MAX_COMIGRATION_WARNINGS)
                .map(|b| comigration(&lane.label, b)),
        );
        if merged.len() > MAX_COMIGRATION_WARNINGS {
            warnings.push(format!(
                "{}: {} further merged bands",
                lane.label,
                merged.len() - MAX_COMIGRATION_WARNINGS
            ));
        }
        output.push(GelLane {
            label: lane.label.clone(),
            bands,
        });
    }
    Ok(Gel {
        agarose_percent: percent,
        separation: (small.round() as usize, large.round() as usize),
        lanes: output,
        warnings,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::restriction::{DigestConditions, ENZYMES, digest};

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }
    fn fragments(lengths: &[usize]) -> Lane {
        Lane {
            label: "Sample".into(),
            fragments: lengths
                .iter()
                .map(|&length| LaneFragment {
                    length,
                    circular: false,
                    reference: false,
                })
                .collect(),
            equimolar: true,
        }
    }

    #[test]
    fn separation_follows_tabulated_ranges_and_refuses_extrapolation() {
        for &(percent, small, large) in SEPARATION {
            let (s, l) = separation_range(percent).unwrap();
            assert!(close(s, small) && close(l, large), "{percent}%");
        }
        let (s, l) = separation_range(1.1).unwrap();
        assert!(s < 500. && s > 400. && l < 10_000. && l > 7_000.);
        for percent in [0.4, 1.6, 2.0, f64::NAN, f64::INFINITY] {
            assert!(separation_range(percent).is_err(), "{percent}");
        }
    }

    #[test]
    fn tabulated_limits_are_drawn_at_fixed_run_fractions() {
        assert!((migration(10_000, 1.0).unwrap() - RANGE_TOP).abs() < 1e-9);
        assert!((migration(500, 1.0).unwrap() - RANGE_BOTTOM).abs() < 1e-9);
        let mut previous = 0.;
        for length in [50_000, 10_000, 3_000, 1_000, 500, 100, 10] {
            let p = migration(length, 1.0).unwrap();
            assert!(p > previous && p < 1., "{length}");
            previous = p;
        }
        // A denser gel retards the same fragment.
        assert!(migration(1_000, 0.7).unwrap() > migration(1_000, 1.5).unwrap());
        assert!(migration(0, 1.0).is_err());
    }

    #[test]
    fn equimolar_lanes_scale_intensity_with_length() {
        let gel = simulate(&[fragments(&[4_000, 1_000])], 1.0).unwrap();
        let bands = &gel.lanes[0].bands;
        assert_eq!(bands[0].lengths, [4_000]);
        assert!(close(bands[0].intensity, 1.) && close(bands[1].intensity, 0.25));
        assert!(gel.warnings.is_empty());
        let ladder = simulate(&[Lane::custom_ladder("Custom", &[4_000, 1_000])], 1.0).unwrap();
        assert!(ladder.lanes[0].bands.iter().all(|b| close(b.intensity, 1.)));
    }

    #[test]
    fn unresolved_and_identical_fragments_merge_with_explanations() {
        let gel = simulate(&[fragments(&[3_000, 1_010, 1_000, 500, 500])], 1.0).unwrap();
        let lengths = gel.lanes[0]
            .bands
            .iter()
            .map(|b| b.lengths.clone())
            .collect::<Vec<_>>();
        assert_eq!(lengths, [vec![3_000], vec![1_010, 1_000], vec![500, 500]]);
        let merged = &gel.lanes[0].bands[1];
        assert!(merged.bottom > merged.top);
        assert!(gel.warnings.iter().any(|w| w.contains("1,000–1,010 bp")));
        assert!(gel.warnings.iter().any(|w| w.contains("2 × 500 bp")));
    }

    #[test]
    fn out_of_range_and_circular_molecules_are_disclosed() {
        let gel = simulate(&[fragments(&[20_000, 1_000, 100])], 1.0).unwrap();
        assert!(gel.warnings[0].contains("2 fragments lie outside the 500–10,000 bp"));
        let lane = Lane {
            label: "Uncut".into(),
            fragments: vec![LaneFragment {
                length: 5_000,
                circular: true,
                reference: false,
            }],
            equimolar: true,
        };
        let gel = simulate(&[lane], 1.0).unwrap();
        assert!(gel.lanes[0].bands[0].circular);
        assert!(gel.warnings.iter().any(|w| w.contains("topology")));
        assert!(gel.warnings.iter().all(|w| !w.contains("outside")));
    }

    #[test]
    fn invalid_gels_are_refused() {
        assert!(simulate(&[], 1.0).is_err());
        assert!(simulate(&vec![fragments(&[1_000]); MAX_LANES + 1], 1.0).is_err());
        assert!(simulate(&[fragments(&[])], 1.0).is_err());
        assert!(simulate(&[fragments(&[0])], 1.0).is_err());
        assert!(simulate(&[fragments(&[1_000])], 3.0).is_err());
    }

    #[test]
    fn ladder_catalogue_is_ordered_and_marks_reference_bands() {
        for ladder in LADDERS {
            assert!(
                ladder.sizes.windows(2).all(|w| w[0] > w[1]),
                "{}",
                ladder.name
            );
            assert!(ladder.reference.iter().all(|r| ladder.sizes.contains(r)));
            let lane = Lane::ladder(ladder);
            assert_eq!(
                lane.fragments.iter().filter(|f| f.reference).count(),
                ladder.reference.len()
            );
        }
        let gel = simulate(&[Lane::ladder(&LADDERS[0])], 1.0).unwrap();
        assert_eq!(gel.lanes[0].bands.len(), 14);
        assert_eq!(gel.lanes[0].bands.iter().filter(|b| b.reference).count(), 3);
    }

    #[test]
    fn digest_lanes_use_product_lengths_and_topology() {
        let ecori = ENZYMES.iter().position(|e| e.name == "EcoRI").unwrap();
        let mut sequence = "A".repeat(3_000);
        sequence.replace_range(100..106, "GAATTC");
        sequence.replace_range(1_100..1_106, "GAATTC");
        let cut = digest(&sequence, true, &[ecori], DigestConditions::default()).unwrap();
        let lane = Lane::digest("EcoRI", &cut.fragments);
        let mut lengths = lane.fragments.iter().map(|f| f.length).collect::<Vec<_>>();
        lengths.sort();
        assert_eq!(lengths, [1_000, 2_000]);
        assert!(lane.equimolar && lane.fragments.iter().all(|f| !f.circular));
        let uncut = digest(
            &"A".repeat(3_000),
            true,
            &[ecori],
            DigestConditions::default(),
        )
        .unwrap();
        let lane = Lane::digest("EcoRI", &uncut.fragments);
        assert_eq!(lane.fragments.len(), 1);
        assert!(lane.fragments[0].circular && lane.fragments[0].length == 3_000);
    }
}
