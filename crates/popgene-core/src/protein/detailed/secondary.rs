//! Sequence-only Chou–Fasman H/E/C prediction using the 1978 29-protein
//! propensity tables reproduced by ExPASy ProtScale. Rules: nucleation 4/6
//! helix or 3/5 sheet; extend while four-residue mean >=1; retain regions with
//! mean >1.03/1.05 and exceeding the other propensity; resolve overlaps by
//! their mean propensity. Proline is allowed only within the first three helix
//! residues. Turns are not separately classified; C means unassigned/coil.
//! This is an explicitly specified implementation of heuristic rules, not the
//! manual Chou–Fasman–Prevelige interpretation or a validated modern predictor.
//! https://web.expasy.org/protscale/pscale/alpha-helixFasman.html
//! https://web.expasy.org/protscale/pscale/beta-sheetFasman.html
//! https://doi.org/10.1002/9780470122921.ch2
use super::STANDARD;
use std::ops::Range;
const HELIX: [f64; 20] = [
    1.42, 0.70, 1.01, 1.51, 1.13, 0.57, 1., 1.08, 1.16, 1.21, 1.45, 0.67, 0.57, 1.11, 0.98, 0.77,
    0.83, 1.06, 1.08, 0.69,
];
const SHEET: [f64; 20] = [
    0.83, 1.19, 0.54, 0.37, 1.38, 0.75, 0.87, 1.60, 0.74, 1.30, 1.05, 0.89, 0.55, 1.10, 0.93, 0.75,
    1.19, 1.70, 1.37, 1.47,
];
#[derive(Clone, Debug, PartialEq)]
pub struct SecondaryPrediction {
    /// H=predicted helix, E=predicted sheet, C=unassigned/coil; one symbol/residue.
    pub assignments: String,
    /// Half-open coordinates relative to the analyzed sequence.
    pub segments: Vec<(char, Range<usize>)>,
    pub helix_count: usize,
    pub sheet_count: usize,
    pub coil_count: usize,
    pub method: &'static str,
}
fn prefix(values: &[f64]) -> Vec<f64> {
    let mut result = Vec::with_capacity(values.len() + 1);
    result.push(0.);
    for v in values {
        result.push(result.last().unwrap() + v);
    }
    result
}
fn mean(prefix: &[f64], start: usize, end: usize) -> f64 {
    (prefix[end] - prefix[start]) / (end - start) as f64
}
fn candidates(sequence: &[u8], propensity: &[f64], other: &[f64], helix: bool) -> Vec<bool> {
    let n = sequence.len();
    let width = if helix { 6 } else { 5 };
    let required = if helix { 4 } else { 3 };
    if n < width {
        return vec![false; n];
    }
    let p = prefix(propensity);
    let q = prefix(other);
    let threshold = if helix { 1.03 } else { 1.05 };
    // Reach arrays avoid quadratic extension on long homogeneous proteins.
    let mut left: Vec<usize> = (0..=n).collect();
    for i in 1..=n {
        if i + 2 < n && mean(&p, i - 1, i + 3) >= 1. - 1e-12 {
            left[i] = left[i - 1];
        }
    }
    let mut right: Vec<usize> = (0..=n).collect();
    for i in (3..n).rev() {
        if mean(&p, i - 3, i + 1) >= 1. - 1e-12 {
            right[i] = right[i + 1];
        }
    }
    let mut next_proline = vec![n; n + 1];
    let mut previous_proline = vec![None; n + 1];
    for i in (0..n).rev() {
        next_proline[i] = if sequence[i] == b'P' {
            i
        } else {
            next_proline[i + 1]
        };
    }
    for i in 0..n {
        previous_proline[i + 1] = if sequence[i] == b'P' {
            Some(i)
        } else {
            previous_proline[i]
        };
    }
    let mut marks = vec![0i64; n + 1];
    for start in 0..=n - width {
        let end = start + width;
        if propensity[start..end]
            .iter()
            .filter(|&&v| v > threshold)
            .count()
            < required
        {
            continue;
        }
        if helix && next_proline[start + 3] < end {
            continue;
        }
        let (mut a, mut b) = (left[start], right[end]);
        if helix {
            b = b.min(next_proline[end]);
            if let Some(proline) = previous_proline[b] {
                a = a.max(proline.saturating_sub(2));
            }
        }
        if b - a < width {
            continue;
        }
        if mean(&p, a, b) > threshold && mean(&p, a, b) > mean(&q, a, b) {
            marks[a] += 1;
            marks[b] -= 1;
        }
    }
    let mut level = 0;
    marks
        .into_iter()
        .take(n)
        .map(|delta| {
            level += delta;
            level > 0
        })
        .collect()
}

/// Bounded linear-time prediction. No coordinates, contacts or confidence
/// probabilities are inferred. Modifications require a different model.
pub fn predict_secondary(sequence: &str) -> Result<SecondaryPrediction, String> {
    if sequence.is_empty() || sequence.len() > 100_000 {
        return Err("Secondary prediction requires 1–100,000 residues".into());
    }
    let sequence = sequence.to_ascii_uppercase();
    let mut alpha = Vec::with_capacity(sequence.len());
    let mut beta = Vec::with_capacity(sequence.len());
    for b in sequence.bytes() {
        let i = STANDARD
            .iter()
            .position(|&aa| aa == b)
            .ok_or("Secondary prediction requires the 20 standard amino acids without stops")?;
        alpha.push(HELIX[i]);
        beta.push(SHEET[i]);
    }
    let h = candidates(sequence.as_bytes(), &alpha, &beta, true);
    let e = candidates(sequence.as_bytes(), &beta, &alpha, false);
    let mut assignments = vec![b'C'; sequence.len()];
    let mut i = 0;
    while i < sequence.len() {
        if h[i] && e[i] {
            let start = i;
            while i < sequence.len() && h[i] && e[i] {
                i += 1;
            }
            let a: f64 = alpha[start..i].iter().sum();
            let b: f64 = beta[start..i].iter().sum();
            // Exact ties remain unassigned rather than choosing one arbitrarily.
            let code = if a > b {
                b'H'
            } else if b > a {
                b'E'
            } else {
                b'C'
            };
            assignments[start..i].fill(code);
        } else {
            assignments[i] = if h[i] {
                b'H'
            } else if e[i] {
                b'E'
            } else {
                b'C'
            };
            i += 1;
        }
    }
    // Unioning adjacent helix candidates must not turn an allowed N-cap
    // proline into an internal helix proline. Break the merged assignment there.
    let mut helix_start = 0;
    for (i, code) in assignments.iter_mut().enumerate() {
        if *code != b'H' {
            helix_start = i + 1;
        } else if sequence.as_bytes()[i] == b'P' && i >= helix_start + 3 {
            *code = b'C';
            helix_start = i + 1;
        }
    }
    let mut segments = Vec::new();
    let mut start = 0;
    for end in 1..=assignments.len() {
        if end == assignments.len() || assignments[end] != assignments[start] {
            segments.push((assignments[start] as char, start..end));
            start = end;
        }
    }
    Ok(SecondaryPrediction {
        helix_count: assignments.iter().filter(|&&c| c == b'H').count(),
        sheet_count: assignments.iter().filter(|&&c| c == b'E').count(),
        coil_count: assignments.iter().filter(|&&c| c == b'C').count(),
        assignments: String::from_utf8(assignments).unwrap(),
        segments,
        method: "Chou–Fasman 1978 H/E/C heuristic (ExPASy 29-protein propensity tables); turns unassigned; no 3D structure or calibrated confidence",
    })
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn published_extreme_propensities_and_proline_breaker() {
        assert_eq!(
            predict_secondary("AAAAAAAAAA").unwrap().assignments,
            "HHHHHHHHHH"
        );
        assert_eq!(
            predict_secondary("VVVVVVVVVV").unwrap().assignments,
            "EEEEEEEEEE"
        );
        assert_eq!(
            predict_secondary("GGGGGGGGGG").unwrap().assignments,
            "CCCCCCCCCC"
        );
        let p = predict_secondary("AAAAAAAAAAPAAAAAAAAAA").unwrap();
        assert!(
            p.segments
                .iter()
                .filter(|(c, _)| *c == 'H')
                .all(|(_, r)| !r.contains(&10) || 10 < r.start + 3)
        );
        assert!(predict_secondary("AAXAAA").is_err());
    }
    #[test]
    fn bounds_and_linear_scaling_on_homopolymer() {
        assert!(predict_secondary("").is_err());
        let p = predict_secondary(&"A".repeat(100_000)).unwrap();
        assert_eq!(p.helix_count, 100_000);
    }
}
