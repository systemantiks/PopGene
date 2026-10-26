//! Sequence-only properties of one unmodified, linear polypeptide.
//!
//! Neutral molecular mass is the sum of free amino-acid masses minus one water
//! per peptide bond. Factual average/monoisotopic masses match the public
//! Biopython 1.88 IUPACData tables. The theoretical pI uses the Bjellqvist pKa
//! set, including residue-dependent terminal values, and Henderson–Hasselbalch
//! charge balance. This is an independent implementation, not a structure or
//! post-translational modification model.
//!
//! Methods: https://web.expasy.org/compute_pi/pi_tool-doc.html
//! Tables: https://github.com/biopython/biopython/blob/biopython-188/Bio/Data/IUPACData.py
//! pKa reference: https://biopython.org/docs/latest/api/Bio.SeqUtils.IsoelectricPoint.html

mod detailed;
pub use detailed::*;

/// All 20 standard amino acids, in the same order as the mass tables.
const STANDARD: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";
const NAMES: [&str; 20] = [
    "Alanine",
    "Cysteine",
    "Aspartate",
    "Glutamate",
    "Phenylalanine",
    "Glycine",
    "Histidine",
    "Isoleucine",
    "Lysine",
    "Leucine",
    "Methionine",
    "Asparagine",
    "Proline",
    "Glutamine",
    "Arginine",
    "Serine",
    "Threonine",
    "Valine",
    "Tryptophan",
    "Tyrosine",
];
const AVERAGE_FREE: [f64; 20] = [
    89.0932, 121.1582, 133.1027, 147.1293, 165.1891, 75.0666, 155.1546, 131.1729, 146.1876,
    131.1729, 149.2113, 132.1179, 115.1305, 146.1445, 174.2010, 105.0926, 119.1192, 117.1463,
    204.2252, 181.1885,
];
const MONO_FREE: [f64; 20] = [
    89.047678, 121.019749, 133.037508, 147.053158, 165.078979, 75.032028, 155.069477, 131.094629,
    146.105528, 131.094629, 149.051049, 132.053492, 115.063329, 146.069142, 174.111676, 105.042593,
    119.058243, 117.078979, 204.089878, 181.073893,
];

#[derive(Clone, Debug, PartialEq)]
pub struct ResidueCount {
    pub symbol: char,
    pub name: &'static str,
    pub count: usize,
    /// Percentage of all stored sequence symbols, including stops/unresolved residues.
    pub percent: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ProteinMass {
    pub average_da: f64,
    pub monoisotopic_da: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ProteinAbsorbance {
    /// Estimated molar extinction at 280 nm, M^-1 cm^-1; all cysteines reduced.
    pub reduced_extinction: usize,
    /// Alternative extinction if every possible cysteine pair forms a cystine.
    pub paired_cysteine_extinction: usize,
    /// Estimated A280 for 1 mg/mL over 1 cm, using reduced extinction/average mass.
    pub reduced_absorbance: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ProteinProperties {
    pub symbol_count: usize,
    pub residue_count: usize,
    pub stop_count: usize,
    pub composition: Vec<ResidueCount>,
    pub mass: Option<ProteinMass>,
    pub theoretical_pi: Option<f64>,
    pub charge_at_ph7: Option<f64>,
    pub absorbance: Option<ProteinAbsorbance>,
    pub unavailable_reason: Option<String>,
}

/// Analyze stored sequence without discarding ambiguity, stops or unsupported
/// residues. Only the 20 standard residues support numeric mass and pI here.
/// Lowercase is accepted; no formatting characters are silently removed.
pub fn analyze(sequence: &str) -> ProteinProperties {
    let mut counts = [0usize; 26];
    let mut stop_count = 0;
    let mut invalid_count = 0;
    let mut symbol_count = 0;
    for c in sequence.chars() {
        symbol_count += 1;
        let c = c.to_ascii_uppercase();
        if c.is_ascii_uppercase() {
            counts[(c as u8 - b'A') as usize] += 1;
        } else if c == '*' {
            stop_count += 1;
        } else {
            invalid_count += 1;
        }
    }
    let residue_count = symbol_count - stop_count;
    let percent = |count| {
        if symbol_count == 0 {
            0.
        } else {
            100. * count as f64 / symbol_count as f64
        }
    };
    let mut composition: Vec<_> = STANDARD
        .iter()
        .zip(NAMES)
        .map(|(&symbol, name)| {
            let count = counts[(symbol - b'A') as usize];
            ResidueCount {
                symbol: symbol as char,
                name,
                count,
                percent: percent(count),
            }
        })
        .collect();
    let extras = [
        (b'B', "Aspartate or asparagine"),
        (b'J', "Isoleucine or leucine"),
        (b'O', "Pyrrolysine (not modeled)"),
        (b'U', "Selenocysteine (not modeled)"),
        (b'X', "Unknown residue"),
        (b'Z', "Glutamate or glutamine"),
    ];
    let mut unsupported = Vec::new();
    for (symbol, name) in extras {
        let count = counts[(symbol - b'A') as usize];
        if count > 0 {
            unsupported.push(symbol as char);
            composition.push(ResidueCount {
                symbol: symbol as char,
                name,
                count,
                percent: percent(count),
            });
        }
    }
    if invalid_count > 0 {
        composition.push(ResidueCount {
            symbol: '?',
            name: "Invalid sequence symbols",
            count: invalid_count,
            percent: percent(invalid_count),
        });
    }
    let mut unavailable_reason = if symbol_count == 0 {
        Some("Add a protein sequence to calculate its properties.".into())
    } else if !unsupported.is_empty() || stop_count > 0 || invalid_count > 0 {
        let mut reasons = Vec::new();
        if !unsupported.is_empty() {
            reasons.push(format!(
                "residues {} are ambiguous or outside this model",
                unsupported
                    .iter()
                    .map(char::to_string)
                    .collect::<Vec<_>>()
                    .join(", ")
            ));
        }
        if stop_count > 0 {
            reasons.push(format!(
                "{stop_count} stop symbol(s) do not define one continuous polypeptide"
            ));
        }
        if invalid_count > 0 {
            reasons.push(format!("{invalid_count} invalid sequence symbol(s)"));
        }
        Some(format!(
            "Physical properties unavailable: {}. No sequence symbols were discarded.",
            reasons.join("; ")
        ))
    } else {
        None
    };
    let (mass, theoretical_pi, charge_at_ph7, absorbance) = if unavailable_reason.is_none() {
        let sum = |table: &[f64; 20], water: f64| {
            STANDARD
                .iter()
                .enumerate()
                .map(|(i, symbol)| counts[(symbol - b'A') as usize] as f64 * table[i])
                .sum::<f64>()
                - (residue_count - 1) as f64 * water
        };
        let first = sequence.as_bytes()[0].to_ascii_uppercase();
        let last = sequence.as_bytes()[sequence.len() - 1].to_ascii_uppercase();
        let pi = isoelectric_point(&counts, first, last);
        if pi.is_none() {
            unavailable_reason = Some("Theoretical pI unavailable: this model has no charge-balance root between pH 0 and 14. Molecular mass is still available.".into());
        }
        let mass = ProteinMass {
            average_da: sum(&AVERAGE_FREE, 18.0153),
            monoisotopic_da: sum(&MONO_FREE, 18.010565),
        };
        // Pace/Edelhoch coefficients, as documented by ExPASy ProtParam.
        let reduced_extinction =
            counts[(b'W' - b'A') as usize] * 5500 + counts[(b'Y' - b'A') as usize] * 1490;
        let paired_cysteine_extinction =
            reduced_extinction + (counts[(b'C' - b'A') as usize] / 2) * 125;
        (
            Some(mass),
            pi,
            Some(charge(&counts, first, last, 7.)),
            Some(ProteinAbsorbance {
                reduced_extinction,
                paired_cysteine_extinction,
                reduced_absorbance: reduced_extinction as f64 / mass.average_da,
            }),
        )
    } else {
        (None, None, None, None)
    };
    ProteinProperties {
        symbol_count,
        residue_count,
        stop_count,
        composition,
        mass,
        theoretical_pi,
        charge_at_ph7,
        absorbance,
        unavailable_reason,
    }
}

fn charge(counts: &[usize; 26], first: u8, last: u8, ph: f64) -> f64 {
    let n_term = match first {
        b'A' => 7.59,
        b'M' => 7.0,
        b'S' => 6.93,
        b'P' => 8.36,
        b'T' => 6.82,
        b'V' => 7.44,
        b'E' => 7.7,
        _ => 7.5,
    };
    let c_term = match last {
        b'D' => 4.55,
        b'E' => 4.75,
        _ => 3.55,
    };
    let basic = |pka: f64| 1. / (1. + 10_f64.powf(ph - pka));
    let acidic = |pka: f64| 1. / (1. + 10_f64.powf(pka - ph));
    let positive: f64 = [(b'K', 10.0), (b'R', 12.0), (b'H', 5.98)]
        .iter()
        .map(|&(aa, pka)| counts[(aa - b'A') as usize] as f64 * basic(pka))
        .sum();
    let negative: f64 = [(b'D', 4.05), (b'E', 4.45), (b'C', 9.0), (b'Y', 10.0)]
        .iter()
        .map(|&(aa, pka)| counts[(aa - b'A') as usize] as f64 * acidic(pka))
        .sum();
    positive + basic(n_term) - negative - acidic(c_term)
}

fn isoelectric_point(counts: &[usize; 26], first: u8, last: u8) -> Option<f64> {
    // Wider than the usual 4.05–12 search: acidic short/long peptides can have
    // charge-balance roots outside that interval. The pKa model is unchanged.
    let (mut low, mut high) = (0., 14.);
    if charge(counts, first, last, low) < 0. || charge(counts, first, last, high) > 0. {
        return None;
    }
    for _ in 0..48 {
        let midpoint = (low + high) / 2.;
        if charge(counts, first, last, midpoint) > 0. {
            low = midpoint;
        } else {
            high = midpoint;
        }
    }
    Some((low + high) / 2.)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn peptide_mass_includes_exactly_one_pair_of_termini() {
        let a = analyze("A").mass.unwrap();
        let aa = analyze("AA").mass.unwrap();
        assert!((a.average_da - 89.0932).abs() < 1e-8);
        assert!((aa.average_da - 160.1711).abs() < 1e-8);
        assert!((aa.monoisotopic_da - 160.084791).abs() < 1e-8);
        assert!((analyze("GG").mass.unwrap().average_da - 132.1179).abs() < 1e-8);
    }

    #[test]
    fn independent_canonical_fixture_values_and_composition() {
        let p = analyze("ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYKK");
        assert_eq!(p.residue_count, 42);
        let k = p.composition.iter().find(|r| r.symbol == 'K').unwrap();
        assert_eq!(k.count, 4);
        assert!((k.percent - 100. * 4. / 42.).abs() < 1e-9);
        assert!((p.mass.unwrap().average_da - 5029.7561).abs() < 1e-7);
        assert!((p.mass.unwrap().monoisotopic_da - 5026.429165).abs() < 1e-7);
        assert!((p.theoretical_pi.unwrap() - 8.764493).abs() < 1e-4);
        assert!((p.charge_at_ph7.unwrap() - 1.9522694241973433).abs() < 1e-10);
        let absorbance = p.absorbance.unwrap();
        assert_eq!(absorbance.reduced_extinction, 13980);
        assert_eq!(absorbance.paired_cysteine_extinction, 14105);
        assert!((absorbance.reduced_absorbance - 13980. / 5029.7561).abs() < 1e-10);
        assert!((p.composition.iter().map(|r| r.percent).sum::<f64>() - 100.).abs() < 1e-9);
    }

    #[test]
    fn bjellqvist_terminal_values_and_extreme_acidic_sequences() {
        // A single alanine has no charged side chain: pI=(7.59+3.55)/2.
        assert!((analyze("A").theoretical_pi.unwrap() - 5.57).abs() < 1e-9);
        assert!((analyze("PETER").theoretical_pi.unwrap() - 4.532085).abs() < 1e-4);
        assert!((analyze("INGAR").theoretical_pi.unwrap() - 9.75).abs() < 1e-4);
        assert!((analyze(&"D".repeat(100)).theoretical_pi.unwrap() - 2.052985).abs() < 1e-4);
        let outside = analyze(&"D".repeat(100_000));
        assert!(outside.theoretical_pi.is_none());
        assert!(outside.mass.is_some());
    }

    #[test]
    fn ambiguity_stop_nonstandard_and_invalid_symbols_are_never_cleaned_away() {
        for sequence in [
            "", "AX", "AB", "AZ", "AJ", "AU", "AO", "A*", "A*A", "A A", "A🧬",
        ] {
            let p = analyze(sequence);
            assert!(p.mass.is_none(), "{sequence}");
            assert!(p.theoretical_pi.is_none(), "{sequence}");
            assert!(p.charge_at_ph7.is_none(), "{sequence}");
            assert!(p.absorbance.is_none(), "{sequence}");
            assert!(p.unavailable_reason.is_some(), "{sequence}");
        }
        let p = analyze("AXZ*");
        assert_eq!((p.symbol_count, p.residue_count, p.stop_count), (4, 3, 1));
        assert_eq!(p.composition.iter().map(|r| r.count).sum::<usize>(), 3);
        assert_eq!(
            p.composition
                .iter()
                .find(|r| r.symbol == 'X')
                .unwrap()
                .percent,
            25.
        );
        assert_eq!(analyze("acdef"), analyze("ACDEF"));
    }

    #[test]
    fn extinction_distinguishes_reduced_and_paired_cysteine_scenarios() {
        let a = analyze("WYC").absorbance.unwrap();
        assert_eq!(
            (a.reduced_extinction, a.paired_cysteine_extinction),
            (6990, 6990)
        );
        let b = analyze("WYCCC").absorbance.unwrap();
        assert_eq!(
            (b.reduced_extinction, b.paired_cysteine_extinction),
            (6990, 7115)
        );
        assert_eq!(analyze("A").absorbance.unwrap().reduced_absorbance, 0.);
    }
}
