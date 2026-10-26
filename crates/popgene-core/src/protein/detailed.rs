//! Explicit sequence scope, alternate charge models, and declared modifications.
//! Mass deltas are factual Unimod records 1, 2, 21, 35, 2020. Original EMBOSS
//! pKa values: https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/emboss/data/Epk.dat
//! Grimsley 2009 measured protein-group means (doi:10.1002/pro.19), completed
//! with Arg 12.04 exactly as published in Kozlowski 2016 Table 4
//! (doi:10.1186/s13062-016-0159-9). Arg was not measured in the former summary.
//! No sequence processing or modification is inferred from annotations.
use super::{ProteinMass, ProteinProperties, STANDARD};
use serde::{Deserialize, Serialize};
use std::{collections::HashSet, ops::Range};
mod secondary;
pub use secondary::{SecondaryPrediction, predict_secondary};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum PkaModel {
    #[default]
    Bjellqvist,
    Emboss,
    Grimsley2009,
}
impl PkaModel {
    pub fn label(self) -> &'static str {
        match self {
            Self::Bjellqvist => "Bjellqvist",
            Self::Emboss => "Original EMBOSS",
            Self::Grimsley2009 => "Grimsley 2009 (IPC Arg)",
        }
    }

    /// Identifies the parameter source separately from the shared charge solver.
    pub fn description(self) -> &'static str {
        match self {
            Self::Bjellqvist => "Bjellqvist pKa values with residue-specific N/C termini.",
            Self::Emboss => "Original EMBOSS Epk.dat pKa values; N-terminus 7.5, C-terminus 3.6.",
            Self::Grimsley2009 => {
                "Grimsley et al. 2009 protein-group mean pKa values, with Arg 12.04 supplied by Kozlowski 2016 Table 4: N-terminus 7.7, C-terminus 3.3, Cys 6.8, Asp 3.5, Glu 4.2, His 6.6, Lys 10.5, Tyr 10.3. One independent ionizable group per residue; unmodified cysteines are reduced."
            }
        }
    }

    pub fn references(self) -> &'static [&'static str] {
        match self {
            Self::Bjellqvist => &["https://web.expasy.org/compute_pi/pi_tool-doc.html"],
            Self::Emboss => &[
                "https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/emboss/data/Epk.dat",
            ],
            Self::Grimsley2009 => &[
                "https://doi.org/10.1002/pro.19",
                "https://doi.org/10.1186/s13062-016-0159-9",
            ],
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Modification {
    NAcetyl,
    CAmide,
    /// Zero-based position on the original sequence, not relative to a selection.
    LysAcetyl(usize),
    MetOxidation(usize),
    Phosphorylation(usize),
    /// Two distinct cysteines, in original-sequence coordinates.
    Disulfide(usize, usize),
}
impl Modification {
    pub fn label(&self) -> String {
        match self {
            Self::NAcetyl => "N-terminal acetylation (Unimod:1)".into(),
            Self::CAmide => "C-terminal amidation (Unimod:2)".into(),
            Self::LysAcetyl(p) => format!("K{} acetylation (Unimod:1)", p.saturating_add(1)),
            Self::MetOxidation(p) => format!("M{} oxidation (Unimod:35)", p.saturating_add(1)),
            Self::Phosphorylation(p) => format!(
                "Residue {} phosphorylation (Unimod:21)",
                p.saturating_add(1)
            ),
            Self::Disulfide(a, b) => format!(
                "C{}–C{} disulfide (Unimod:2020)",
                a.saturating_add(1),
                b.saturating_add(1)
            ),
        }
    }
    fn mass_delta(&self) -> ProteinMass {
        let (average_da, monoisotopic_da) = match self {
            Self::NAcetyl | Self::LysAcetyl(_) => (42.0367, 42.010565),
            Self::CAmide => (-0.9848, -0.984016),
            Self::MetOxidation(_) => (15.9994, 15.994915),
            Self::Phosphorylation(_) => (79.9799, 79.966331),
            Self::Disulfide(_, _) => (-2.0159, -2.015650),
        };
        ProteinMass {
            average_da,
            monoisotopic_da,
        }
    }
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ProteinOptions {
    /// A contiguous half-open interval; None means the full sequence. A selected
    /// interval is treated as an excised linear peptide with new free termini.
    pub range: Option<Range<usize>>,
    pub model: PkaModel,
    pub ph: f64,
    pub modifications: Vec<Modification>,
    pub secondary_structure: bool,
}
impl Default for ProteinOptions {
    fn default() -> Self {
        Self {
            range: None,
            model: PkaModel::Bjellqvist,
            ph: 7.,
            modifications: vec![],
            secondary_structure: false,
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct DetailedProteinProperties {
    /// Unmodified composition and reference properties, retained for comparison.
    pub unmodified: ProteinProperties,
    pub range: Range<usize>,
    pub model: PkaModel,
    pub ph: f64,
    pub modifications: Vec<Modification>,
    /// Includes explicitly requested mass deltas.
    pub mass: Option<ProteinMass>,
    pub theoretical_pi: Option<f64>,
    pub charge: Option<f64>,
    /// Charge sampled at pH 0, 0.1, ..., 14; empty if charge is unsupported.
    pub charge_curve: Vec<(f64, f64)>,
    /// Trp/Tyr + explicitly specified cystine bonds, rather than assuming all pair.
    pub extinction_280: Option<usize>,
    pub absorbance_1mg_ml: Option<f64>,
    pub secondary: Option<SecondaryPrediction>,
    pub warnings: Vec<String>,
}

pub fn analyze_with_options(
    sequence: &str,
    options: &ProteinOptions,
) -> Result<DetailedProteinProperties, String> {
    if !options.ph.is_finite() || !(0.0..=14.0).contains(&options.ph) {
        return Err("pH must be between 0 and 14".into());
    }
    let range = options.range.clone().unwrap_or(0..sequence.len());
    if range.is_empty() || range.end > sequence.len() || !sequence.is_ascii() {
        return Err("Choose a nonempty, valid protein interval".into());
    }
    if range.len() > 1_000_000 {
        return Err("Analyze at most 1,000,000 residues at once".into());
    }
    if options.modifications.len() > 10_000 {
        return Err("At most 10,000 explicit modifications are supported".into());
    }
    let selected = sequence
        .get(range.clone())
        .ok_or("Invalid protein interval")?
        .to_ascii_uppercase();
    let unmodified = super::analyze(&selected);
    let mut warnings=vec!["Sequence estimates assume one linear chain. Charge models omit solvent, temperature and residue-environment effects; they are not experimental measurements.".into()];
    if range.start != 0 || range.end != sequence.len() {
        warnings.push("The selected interval is modeled as an excised peptide with new N/C termini; this is not its charge while embedded in the intact protein.".into());
    }
    let mut counts = [0usize; 26];
    for &b in selected.as_bytes() {
        if b.is_ascii_uppercase() {
            counts[(b - b'A') as usize] += 1;
        }
    }
    let (mut n_terminal, mut c_terminal) = (true, true);
    let mut occupied = HashSet::new();
    let mut bonds = 0;
    let mut charge_supported = true;
    let mut optics_supported = true;
    let mut mass = unmodified.mass;
    for modification in &options.modifications {
        let positions: Vec<(usize, &[u8])> = match modification {
            Modification::NAcetyl => {
                if !n_terminal {
                    return Err("Duplicate N-terminal modification".into());
                }
                n_terminal = false;
                vec![]
            }
            Modification::CAmide => {
                if !c_terminal {
                    return Err("Duplicate C-terminal modification".into());
                }
                c_terminal = false;
                vec![]
            }
            Modification::LysAcetyl(p) => vec![(*p, b"K")],
            Modification::MetOxidation(p) => vec![(*p, b"M")],
            Modification::Phosphorylation(p) => {
                charge_supported = false;
                optics_supported = false;
                vec![(*p, b"STY")]
            }
            Modification::Disulfide(a, b) => {
                bonds += 1;
                vec![(*a, b"C"), (*b, b"C")]
            }
        };
        for (position, allowed) in positions {
            if !range.contains(&position) {
                return Err(format!(
                    "Modification at residue {} is outside the selected interval",
                    position.saturating_add(1)
                ));
            }
            if !allowed.contains(&sequence.as_bytes()[position].to_ascii_uppercase()) {
                return Err(format!(
                    "Modification is incompatible with residue {}",
                    position.saturating_add(1)
                ));
            }
            if !occupied.insert(position) {
                return Err(format!(
                    "Residue {} has conflicting or duplicate modifications",
                    position.saturating_add(1)
                ));
            }
        }
        if let Some(mass) = mass.as_mut() {
            let delta = modification.mass_delta();
            mass.average_da += delta.average_da;
            mass.monoisotopic_da += delta.monoisotopic_da;
        }
        match modification {
            Modification::LysAcetyl(_) => counts[(b'K' - b'A') as usize] -= 1,
            Modification::Disulfide(_, _) => counts[(b'C' - b'A') as usize] -= 2,
            _ => {}
        }
    }
    if !charge_supported {
        warnings.push("Phosphorylation mass is included; pI, charge and optical estimates are withheld because modified-residue pKa/extinction values are not modeled.".into());
    }
    if !options.modifications.is_empty() {
        warnings.push("Only the declared chemical changes are included. Ionizable groups removed by acetylation, amidation or disulfide bonding are excluded; all remaining pKa values stay in the selected reference model.".into());
    }
    let mut result = DetailedProteinProperties {
        unmodified,
        range,
        model: options.model,
        ph: options.ph,
        modifications: options.modifications.clone(),
        mass,
        theoretical_pi: None,
        charge: None,
        charge_curve: vec![],
        extinction_280: None,
        absorbance_1mg_ml: None,
        secondary: None,
        warnings,
    };
    if let Some(mass) = result.mass {
        if charge_supported {
            let charge = |ph| {
                modeled_charge(
                    &counts,
                    selected.as_bytes()[0],
                    selected.as_bytes()[selected.len() - 1],
                    ph,
                    options.model,
                    n_terminal,
                    c_terminal,
                )
            };
            result.charge = Some(charge(options.ph));
            result.charge_curve = (0..=140)
                .map(|i| {
                    let ph = i as f64 / 10.;
                    (ph, charge(ph))
                })
                .collect();
            let (mut lo, mut hi) = (0., 14.);
            if charge(lo) > 0. && charge(hi) < 0. {
                for _ in 0..48 {
                    let mid = (lo + hi) / 2.;
                    if charge(mid) > 0. {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                result.theoretical_pi = Some((lo + hi) / 2.);
            } else {
                result.warnings.push("No unique charge-balance root exists within pH 0–14 for this model; pI is unavailable.".into());
            }
        }
        if optics_supported {
            let extinction = counts[(b'W' - b'A') as usize] * 5500
                + counts[(b'Y' - b'A') as usize] * 1490
                + bonds * 125;
            result.extinction_280 = Some(extinction);
            result.absorbance_1mg_ml = Some(extinction as f64 / mass.average_da);
        }
        if options.secondary_structure {
            if options.modifications.is_empty() {
                result.secondary = Some(predict_secondary(&selected)?);
            } else {
                result.warnings.push("Secondary structure is withheld for modified chains; the sequence predictor has no modification model.".into());
            }
        }
    }
    Ok(result)
}

fn modeled_charge(
    counts: &[usize; 26],
    first: u8,
    last: u8,
    ph: f64,
    model: PkaModel,
    n_terminal: bool,
    c_terminal: bool,
) -> f64 {
    let positive = |pka: f64| 1. / (1. + 10_f64.powf(ph - pka));
    let negative = |pka: f64| 1. / (1. + 10_f64.powf(pka - ph));
    let (n, c) = match model {
        PkaModel::Bjellqvist => (
            match first {
                b'A' => 7.59,
                b'M' => 7.,
                b'S' => 6.93,
                b'P' => 8.36,
                b'T' => 6.82,
                b'V' => 7.44,
                b'E' => 7.7,
                _ => 7.5,
            },
            match last {
                b'D' => 4.55,
                b'E' => 4.75,
                _ => 3.55,
            },
        ),
        PkaModel::Emboss => (7.5, 3.6),
        PkaModel::Grimsley2009 => (7.7, 3.3),
    };
    let (basic, acidic) = match model {
        PkaModel::Bjellqvist => (
            [(b'K', 10.), (b'R', 12.), (b'H', 5.98)],
            [(b'D', 4.05), (b'E', 4.45), (b'C', 9.), (b'Y', 10.)],
        ),
        PkaModel::Emboss => (
            [(b'K', 10.8), (b'R', 12.5), (b'H', 6.5)],
            [(b'D', 3.9), (b'E', 4.1), (b'C', 8.5), (b'Y', 10.1)],
        ),
        PkaModel::Grimsley2009 => (
            [(b'K', 10.5), (b'R', 12.04), (b'H', 6.6)],
            [(b'D', 3.5), (b'E', 4.2), (b'C', 6.8), (b'Y', 10.3)],
        ),
    };
    basic
        .iter()
        .map(|&(b, p)| counts[(b - b'A') as usize] as f64 * positive(p))
        .sum::<f64>()
        - acidic
            .iter()
            .map(|&(b, p)| counts[(b - b'A') as usize] as f64 * negative(p))
            .sum::<f64>()
        + if n_terminal { positive(n) } else { 0. }
        - if c_terminal { negative(c) } else { 0. }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn grimsley_published_table_matches_independent_charge_oracle() {
        // Factual Grimsley/IPC table evaluated independently with Biopython's
        // Henderson–Hasselbalch implementation and an 80-step root search.
        let sequence = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYKK";
        let result = analyze_with_options(
            sequence,
            &ProteinOptions {
                model: PkaModel::Grimsley2009,
                ..ProteinOptions::default()
            },
        )
        .unwrap();
        assert!((result.theoretical_pi.unwrap() - 7.603081545360311).abs() < 1e-11);
        assert!((result.charge.unwrap() - 1.1785952157817565).abs() < 1e-12);
        assert_eq!(format!("{:.2}", result.theoretical_pi.unwrap()), "7.60");
        assert_eq!(format!("{:.2}", result.charge.unwrap()), "1.18");
        assert_eq!(result.extinction_280, Some(13980));
        assert_eq!(result.mass, super::super::analyze(sequence).mass);
        assert_eq!(ProteinOptions::default().model, PkaModel::Bjellqvist);
        assert!(result.charge_curve.windows(2).all(|p| p[0].1 >= p[1].1));
    }

    #[test]
    fn grimsley_selected_peptide_and_declared_modifications_use_same_group_model() {
        let mut options = ProteinOptions {
            range: Some(1..3),
            model: PkaModel::Grimsley2009,
            ..ProteinOptions::default()
        };
        let selected = analyze_with_options("KGGD", &options).unwrap();
        // Only the selected GG is present, with free termini: (7.7 + 3.3)/2.
        assert!((selected.theoretical_pi.unwrap() - 5.5).abs() < 1e-11);
        assert!((selected.charge.unwrap() + 0.1661380443878404).abs() < 1e-12);
        options.range = None;
        options.modifications = vec![
            Modification::NAcetyl,
            Modification::CAmide,
            Modification::LysAcetyl(8),
            Modification::Disulfide(1, 21),
        ];
        let modified =
            analyze_with_options("ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYKK", &options).unwrap();
        assert!((modified.theoretical_pi.unwrap() - 9.805616953021431).abs() < 1e-11);
        assert!((modified.charge.unwrap() - 1.5713230282735124).abs() < 1e-12);
        assert_eq!(modified.extinction_280, Some(14105));
        options.range = Some(0..20);
        assert!(
            analyze_with_options("ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYKK", &options).is_err()
        );
    }

    #[test]
    fn grimsley_does_not_invent_roots_or_modified_residue_parameters() {
        let mut options = ProteinOptions {
            model: PkaModel::Grimsley2009,
            modifications: vec![Modification::NAcetyl, Modification::CAmide],
            ..ProteinOptions::default()
        };
        for sequence in ["GG", "DG", "KG"] {
            let p = analyze_with_options(sequence, &options).unwrap();
            assert!(p.theoretical_pi.is_none());
            assert!(p.charge.is_some());
        }
        options.modifications = vec![Modification::Phosphorylation(0)];
        let p = analyze_with_options("SY", &options).unwrap();
        assert!(p.mass.is_some());
        assert!(p.charge.is_none() && p.theoretical_pi.is_none());
        options.modifications.clear();
        for sequence in ["AX", "AZ", "A*", "A*A"] {
            let p = analyze_with_options(sequence, &options).unwrap();
            assert!(p.mass.is_none() && p.charge.is_none() && p.theoretical_pi.is_none());
        }
    }

    #[test]
    fn models_and_selection_use_declared_termini() {
        let reference = super::super::analyze("ACDEFGHIKLMNPQRSTVWY");
        let result =
            analyze_with_options("ACDEFGHIKLMNPQRSTVWY", &ProteinOptions::default()).unwrap();
        assert!((result.theoretical_pi.unwrap() - reference.theoretical_pi.unwrap()).abs() < 1e-10);
        let selected = analyze_with_options(
            "KGGD",
            &ProteinOptions {
                range: Some(1..3),
                model: PkaModel::Emboss,
                ..ProteinOptions::default()
            },
        )
        .unwrap();
        assert!((selected.theoretical_pi.unwrap() - 5.55).abs() < 1e-10); // two termini: mean of 7.5 and 3.6
        assert_eq!(selected.unmodified.symbol_count, 2);
        assert!(selected.charge_curve.windows(2).all(|p| p[0].1 >= p[1].1));
    }
    #[test]
    fn chemical_deltas_and_invalid_modifications() {
        let opts = ProteinOptions {
            modifications: vec![
                Modification::NAcetyl,
                Modification::CAmide,
                Modification::MetOxidation(0),
                Modification::Disulfide(1, 2),
            ],
            ..ProteinOptions::default()
        };
        let p = analyze_with_options("MCCY", &opts).unwrap();
        let original = super::super::analyze("MCCY").mass.unwrap();
        assert!(
            (p.mass.unwrap().monoisotopic_da - original.monoisotopic_da - 55.005814).abs() < 1e-9
        );
        assert_eq!(p.extinction_280, Some(1615)); // one Tyr and one explicit cystine
        assert!(analyze_with_options("MAC", &opts).is_err());
        assert!(
            analyze_with_options(
                "MCC",
                &ProteinOptions {
                    modifications: vec![Modification::Disulfide(1, 1)],
                    ..ProteinOptions::default()
                }
            )
            .is_err()
        );
        let phosphorylated = analyze_with_options(
            "SY",
            &ProteinOptions {
                modifications: vec![Modification::Phosphorylation(0)],
                ..ProteinOptions::default()
            },
        )
        .unwrap();
        assert!(phosphorylated.mass.is_some());
        assert!(phosphorylated.charge.is_none());
        assert!(phosphorylated.extinction_280.is_none());
        let neutral = analyze_with_options(
            "GG",
            &ProteinOptions {
                modifications: vec![Modification::NAcetyl, Modification::CAmide],
                ..ProteinOptions::default()
            },
        )
        .unwrap();
        assert_eq!(neutral.charge, Some(0.));
        assert!(neutral.theoretical_pi.is_none());
    }
    #[test]
    fn unresolved_symbols_never_gain_invented_physical_values() {
        let p = analyze_with_options(
            "AX*",
            &ProteinOptions {
                secondary_structure: true,
                ..ProteinOptions::default()
            },
        )
        .unwrap();
        assert!(p.mass.is_none());
        assert!(p.charge.is_none());
        assert!(p.secondary.is_none());
        assert_eq!(p.unmodified.symbol_count, 3);
    }
}
