//! Primer3 2.6.1 thermodynamics and primer-pair design, compiled from the
//! GPL-2.0-or-later sources distributed with primer3-py 2.3.1. No Python/runtime
//! executable is required. See vendor/primer3/provenance.json and source notices.
//! Methods: https://primer3.org/manual.html; doi:10.1073/pnas.95.4.1460;
//! doi:10.1021/bi702363u. Parameters use nM DNA, mM ions/dNTP, % DMSO, M formamide.
//! Run design on a worker. Calls serialize access to Primer3's global caches.
use serde::{Deserialize, Serialize};
use std::{
    ffi::{CStr, CString, c_char, c_int},
    ops::Range,
    sync::Mutex,
};

pub const MAX_DESIGN_BASES: usize = 100_000;
pub const METHOD: &str = "Primer3 2.6.1 (primer3-py 2.3.1 source distribution)";
static ENGINE: Mutex<()> = Mutex::new(());

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
#[repr(i32)]
pub enum MeltingModel {
    Breslauer = 0,
    #[default]
    SantaLucia = 1,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
#[repr(i32)]
pub enum SaltCorrection {
    Schildkraut = 0,
    #[default]
    SantaLucia = 1,
    Owczarzy = 2,
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Reaction {
    pub dna_nm: f64,
    pub monovalent_mm: f64,
    pub magnesium_mm: f64,
    pub dntp_mm: f64,
    pub dmso_percent: f64,
    pub formamide_m: f64,
    /// Temperature used for structure ΔG, not the predicted duplex Tm.
    pub temperature_c: f64,
    pub melting_model: MeltingModel,
    pub salt_correction: SaltCorrection,
}
impl Default for Reaction {
    fn default() -> Self {
        Self {
            dna_nm: 50.,
            monovalent_mm: 50.,
            magnesium_mm: 1.5,
            dntp_mm: 0.6,
            dmso_percent: 0.,
            formamide_m: 0.,
            temperature_c: 37.,
            melting_model: MeltingModel::SantaLucia,
            salt_correction: SaltCorrection::SantaLucia,
        }
    }
}
impl Reaction {
    pub fn validate(self) -> Result<(), String> {
        for (name, value, min, max) in [
            ("DNA concentration (nM)", self.dna_nm, 0.001, 1_000_000.),
            ("Monovalent ions (mM)", self.monovalent_mm, 0., 1000.),
            ("Magnesium (mM)", self.magnesium_mm, 0., 100.),
            ("dNTPs (mM)", self.dntp_mm, 0., 100.),
            ("DMSO (%)", self.dmso_percent, 0., 20.),
            ("Formamide (M)", self.formamide_m, 0., 10.),
            ("Structure temperature (°C)", self.temperature_c, 0., 100.),
        ] {
            if !value.is_finite() || !(min..=max).contains(&value) {
                return Err(format!("{name} must be finite and between {min} and {max}"));
            }
        }
        if self.monovalent_mm == 0. && self.magnesium_mm <= self.dntp_mm {
            return Err("At least one positive effective salt concentration is required".into());
        }
        Ok(())
    }
    fn native(self) -> NativeReaction {
        NativeReaction {
            dna_nm: self.dna_nm,
            monovalent_mm: self.monovalent_mm,
            magnesium_mm: self.magnesium_mm,
            dntp_mm: self.dntp_mm,
            dmso_percent: self.dmso_percent,
            formamide_m: self.formamide_m,
            temperature_c: self.temperature_c,
            tm_method: self.melting_model as i32,
            salt_method: self.salt_correction as i32,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Structure {
    pub found: bool,
    pub tm_c: f64,
    pub enthalpy_cal_mol: f64,
    pub entropy_cal_mol_k: f64,
    pub free_energy_cal_mol: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct OligoAnalysis {
    pub tm_c: f64,
    pub gc_percent: f64,
    pub hairpin: Structure,
    pub homodimer: Structure,
    pub warnings: Vec<String>,
}

fn dna(sequence: &str, max: usize) -> Result<String, String> {
    if sequence.len() < 2 || sequence.len() > max {
        return Err(format!("Expected 2–{max} DNA bases"));
    }
    if !sequence
        .bytes()
        .all(|base| matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
    {
        return Err(
            "Thermodynamics requires unambiguous A/C/G/T DNA; no symbols were discarded".into(),
        );
    }
    Ok(sequence.to_ascii_uppercase())
}
fn message(buffer: &[c_char]) -> String {
    // The C bridge writes bounded snprintf strings into a zero-filled buffer.
    unsafe { CStr::from_ptr(buffer.as_ptr()) }
        .to_string_lossy()
        .into_owned()
}
fn structure_locked(
    a: &CString,
    b: &CString,
    r: &NativeReaction,
    hairpin: bool,
) -> Result<Structure, String> {
    let mut result = NativeStructure::default();
    let mut error = [0; 1024];
    let status = unsafe {
        pg_primer_structure(
            a.as_ptr(),
            b.as_ptr(),
            r,
            hairpin as i32,
            &mut result,
            error.as_mut_ptr(),
            error.len(),
        )
    };
    if status != 0 {
        return Err(message(&error));
    }
    Ok(Structure {
        found: result.found != 0,
        tm_c: result.tm,
        enthalpy_cal_mol: result.dh,
        entropy_cal_mol_k: result.ds,
        free_energy_cal_mol: result.dg,
    })
}

/// Perfect-match duplex Tm for an unmodified DNA annealing region (2–60 bases).
/// It is not an annealing-temperature recommendation or a mismatch model.
pub fn melting_temperature(sequence: &str, reaction: Reaction) -> Result<f64, String> {
    reaction.validate()?;
    let sequence = CString::new(dna(sequence, 60)?).unwrap();
    let _guard = ENGINE
        .lock()
        .map_err(|_| "Primer3 worker lock is poisoned")?;
    let tm = unsafe { pg_primer_tm(sequence.as_ptr(), &reaction.native()) };
    if !tm.is_finite() || tm <= -273.15 {
        return Err("Primer3 could not calculate this duplex Tm".into());
    }
    Ok(tm)
}

pub fn analyze(sequence: &str, reaction: Reaction) -> Result<OligoAnalysis, String> {
    let sequence = dna(sequence, 60)?;
    let tm_c = melting_temperature(&sequence, reaction)?;
    let gc_percent = 100.
        * sequence
            .bytes()
            .filter(|b| matches!(b, b'C' | b'G'))
            .count() as f64
        / sequence.len() as f64;
    let c_sequence = CString::new(sequence).unwrap();
    let _guard = ENGINE
        .lock()
        .map_err(|_| "Primer3 worker lock is poisoned")?;
    let native = reaction.native();
    let hairpin = structure_locked(&c_sequence, &c_sequence, &native, true)?;
    let homodimer = structure_locked(&c_sequence, &c_sequence, &native, false)?;
    let mut warnings = vec![];
    if reaction.dmso_percent != 0. || reaction.formamide_m != 0. {
        warnings.push("DMSO/formamide corrections apply to duplex Tm; Primer3 structure calculations do not model these additives".into());
    }
    Ok(OligoAnalysis {
        tm_c,
        gc_percent,
        hairpin,
        homodimer,
        warnings,
    })
}

pub fn heterodimer(first: &str, second: &str, reaction: Reaction) -> Result<Structure, String> {
    reaction.validate()?;
    let first = CString::new(dna(first, 60)?).unwrap();
    let second = CString::new(dna(second, 60)?).unwrap();
    let _guard = ENGINE
        .lock()
        .map_err(|_| "Primer3 worker lock is poisoned")?;
    structure_locked(&first, &second, &reaction.native(), false)
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DesignOptions {
    pub reaction: Reaction,
    pub min_length: usize,
    pub optimal_length: usize,
    pub max_length: usize,
    pub min_tm: f64,
    pub optimal_tm: f64,
    pub max_tm: f64,
    pub min_gc: f64,
    pub max_gc: f64,
    pub max_tm_difference: f64,
    pub product_size: Range<usize>,
    /// Required internal target, zero-based half-open in the supplied template.
    pub target: Option<Range<usize>>,
    pub num_return: usize,
}
impl Default for DesignOptions {
    fn default() -> Self {
        Self {
            reaction: Reaction::default(),
            min_length: 18,
            optimal_length: 20,
            max_length: 27,
            min_tm: 57.,
            optimal_tm: 60.,
            max_tm: 63.,
            min_gc: 20.,
            max_gc: 80.,
            max_tm_difference: 3.,
            product_size: 100..301,
            target: None,
            num_return: 5,
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct DesignedPrimer {
    pub sequence: String,
    /// Forward-axis half-open interval. Circular designs may have end > n;
    /// use binding_segments(n) before writing document features.
    pub binding: Range<usize>,
    pub reverse: bool,
    pub tm_c: f64,
    pub gc_percent: f64,
    pub hairpin_tm_c: f64,
}
impl DesignedPrimer {
    pub fn binding_segments(&self, length: usize) -> Result<Vec<crate::Segment>, String> {
        if length == 0
            || self.binding.is_empty()
            || self.binding.start >= length
            || self.binding.len() > length
        {
            return Err("Invalid designed-primer binding interval".into());
        }
        if self.binding.end <= length {
            Ok(vec![crate::Segment {
                start: self.binding.start,
                end: self.binding.end,
            }])
        } else {
            Ok(vec![
                crate::Segment {
                    start: self.binding.start,
                    end: length,
                },
                crate::Segment {
                    start: 0,
                    end: self.binding.end - length,
                },
            ])
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct PrimerPair {
    pub forward: DesignedPrimer,
    pub reverse: DesignedPrimer,
    /// Circular products use an unwrapped range so `.len()` remains exact.
    pub product: Range<usize>,
    pub penalty: f64,
    pub heterodimer_tm_c: f64,
    pub end_heterodimer_tm_c: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct DesignResult {
    pub pairs: Vec<PrimerPair>,
    pub warnings: Vec<String>,
}

/// Ranked circular PCR candidates from four linearizations. Target intervals
/// may cross the origin by using an end greater than template length. A product
/// spans at most one circumference. Each accepted oligo must have exactly one
/// perfect binding site on the whole circle (including the origin), in its
/// intended orientation. Near-match specificity remains Primer3's per-rotation
/// template screening; this is not genome-wide or an exhaustive circular search.
pub fn design_pairs_circular(
    template: &str,
    options: &DesignOptions,
) -> Result<DesignResult, String> {
    let template = dna(template, MAX_DESIGN_BASES)?;
    let n = template.len();
    if !(1..=20).contains(&options.num_return) {
        return Err("Result count must be 1–20".into());
    }
    if options.product_size.start > n {
        return Err("Minimum product length exceeds the circle circumference".into());
    }
    if let Some(target) = &options.target
        && (target.is_empty() || target.start >= n || target.len() >= n)
    {
        return Err(
            "Circular target must start within the molecule and span less than one circumference"
                .into(),
        );
    }
    let mut seams = if let Some(target) = &options.target {
        let remaining = n - target.len();
        (1..=4)
            .map(|i| (target.end + remaining * i / 5) % n)
            .collect::<Vec<_>>()
    } else {
        (0..4).map(|i| n * i / 4).collect::<Vec<_>>()
    };
    seams.sort_unstable();
    seams.dedup();
    let mut result=DesignResult {pairs:vec![],warnings:vec!["Circular candidates were searched at four template seams, then mapped back to the original origin; the candidate search is bounded, not exhaustive.".into(),"Accepted oligos have one perfect binding site on the entire circle. Near-match screening uses the linearized templates and does not establish genome-wide specificity or reaction yield.".into()]};
    let mut seen = std::collections::HashSet::new();
    let mut ambiguous = 0;
    for seam in seams {
        let rotated = format!("{}{}", &template[seam..], &template[..seam]);
        let mut selected = options.clone();
        selected.product_size.end = selected.product_size.end.min(n + 1);
        selected.num_return = 20;
        selected.target = options.target.as_ref().map(|r| {
            let start = (r.start + n - seam) % n;
            start..start + r.len()
        });
        let native = design_pairs(&rotated, &selected)?;
        for warning in native.warnings {
            if !warning.starts_with("No primer pairs")
                && !warning.starts_with("Specificity screening")
                && !result.warnings.contains(&warning)
            {
                result.warnings.push(warning);
            }
        }
        for mut pair in native.pairs {
            let map = |r: &Range<usize>| {
                let start = (r.start + seam) % n;
                start..start + r.len()
            };
            pair.forward.binding = map(&pair.forward.binding);
            pair.reverse.binding = map(&pair.reverse.binding);
            pair.product = map(&pair.product);
            if !unique_circular_binding(&template, &pair.forward)
                || !unique_circular_binding(&template, &pair.reverse)
            {
                ambiguous += 1;
                continue;
            }
            if seen.insert((
                pair.forward.binding.start,
                pair.forward.sequence.clone(),
                pair.reverse.binding.start,
                pair.reverse.sequence.clone(),
            )) {
                result.pairs.push(pair);
            }
        }
    }
    result.pairs.sort_by(|a, b| {
        a.penalty
            .total_cmp(&b.penalty)
            .then(a.product.start.cmp(&b.product.start))
            .then(a.product.end.cmp(&b.product.end))
    });
    result.pairs.truncate(options.num_return);
    if result.pairs.is_empty() {
        result.warnings.push(
            "No circular primer pairs satisfy the supplied constraints and exact-binding screen"
                .into(),
        );
    }
    if options.product_size.end > n + 1 {
        result
            .warnings
            .push("Maximum product length was limited to one circumference.".into());
    }
    if ambiguous > 0 {
        result.warnings.push(format!("Discarded {ambiguous} rotation candidates with additional perfect-match oligo binding sites on the circle."));
    }
    Ok(result)
}
fn unique_circular_binding(template: &str, primer: &DesignedPrimer) -> bool {
    let forward = primer.sequence.as_bytes();
    let reverse = forward
        .iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            _ => b'C',
        })
        .collect::<Vec<_>>();
    let mut count = 0;
    for start in 0..template.len() {
        for (is_reverse, query) in [(false, forward), (true, reverse.as_slice())] {
            if query
                .iter()
                .enumerate()
                .all(|(i, b)| template.as_bytes()[(start + i) % template.len()] == *b)
            {
                if start != primer.binding.start || is_reverse != primer.reverse {
                    return false;
                }
                count += 1;
            }
        }
    }
    count == 1
}

/// Primer3 thermodynamic oligo screening and ranked PCR pairs on a linear
/// unambiguous template. Specificity is against this template only, not a genome.
/// Product-size range is half-open, consistent with Rust's Range type.
pub fn design_pairs(template: &str, options: &DesignOptions) -> Result<DesignResult, String> {
    options.reaction.validate()?;
    let template = dna(template, MAX_DESIGN_BASES)?;
    if !(2..=36).contains(&options.min_length)
        || options.min_length > options.optimal_length
        || options.optimal_length > options.max_length
        || options.max_length > 36
    {
        return Err("Primer lengths must satisfy 2 ≤ minimum ≤ optimum ≤ maximum ≤ 36".into());
    }
    if ![
        options.min_tm,
        options.optimal_tm,
        options.max_tm,
        options.min_gc,
        options.max_gc,
        options.max_tm_difference,
    ]
    .iter()
    .all(|v| v.is_finite())
        || options.min_tm > options.optimal_tm
        || options.optimal_tm > options.max_tm
        || options.min_tm < 0.
        || options.max_tm > 100.
        || options.min_gc < 0.
        || options.min_gc > options.max_gc
        || options.max_gc > 100.
        || options.max_tm_difference < 0.
        || options.max_tm_difference > 100.
    {
        return Err("Invalid primer temperature or GC constraints".into());
    }
    if options.product_size.start < 2
        || options.product_size.end <= options.product_size.start
        || options.product_size.end > MAX_DESIGN_BASES + 1
        || !(1..=20).contains(&options.num_return)
    {
        return Err("Invalid product-size range or result count (1–20)".into());
    }
    if let Some(target) = &options.target
        && (target.start >= target.end || target.end > template.len())
    {
        return Err("Target interval is outside the template".into());
    }
    let native = NativeDesign {
        reaction: options.reaction.native(),
        min_tm: options.min_tm,
        optimal_tm: options.optimal_tm,
        max_tm: options.max_tm,
        min_gc: options.min_gc,
        max_gc: options.max_gc,
        max_tm_difference: options.max_tm_difference,
        min_length: options.min_length as i32,
        optimal_length: options.optimal_length as i32,
        max_length: options.max_length as i32,
        product_min: options.product_size.start as i32,
        product_max: (options.product_size.end - 1) as i32,
        num_return: options.num_return as i32,
        target_start: options.target.as_ref().map_or(0, |r| r.start as i32),
        target_length: options.target.as_ref().map_or(0, |r| r.len() as i32),
    };
    let c_template = CString::new(template.as_str()).unwrap();
    let mut raw = [NativePair::default(); 20];
    let (mut count, mut error, mut warning) = (0, [0; 1024], [0; 1024]);
    let _guard = ENGINE
        .lock()
        .map_err(|_| "Primer3 worker lock is poisoned")?;
    let status = unsafe {
        pg_design_primers(
            c_template.as_ptr(),
            &native,
            raw.as_mut_ptr(),
            raw.len() as i32,
            &mut count,
            error.as_mut_ptr(),
            error.len(),
            warning.as_mut_ptr(),
            warning.len(),
        )
    };
    if status != 0 {
        return Err(message(&error));
    }
    if !(0..=20).contains(&count) {
        return Err("Invalid Primer3 result count".into());
    }
    let mut pairs = vec![];
    for pair in &raw[..count as usize] {
        let f_start =
            usize::try_from(pair.forward_start).map_err(|_| "Invalid forward coordinate")?;
        let f_len = usize::try_from(pair.forward_length).map_err(|_| "Invalid forward length")?;
        let r_end =
            usize::try_from(pair.reverse_end).map_err(|_| "Invalid reverse coordinate")? + 1;
        let r_len = usize::try_from(pair.reverse_length).map_err(|_| "Invalid reverse length")?;
        if f_len == 0
            || r_len == 0
            || f_start + f_len > r_end
            || r_end > template.len()
            || r_len > r_end
        {
            return Err("Primer3 returned an invalid binding interval".into());
        }
        if r_end - f_start != pair.product_size as usize {
            return Err("Primer3 product length mismatch".into());
        }
        let reverse_sequence: String = template.as_bytes()[r_end - r_len..r_end]
            .iter()
            .rev()
            .map(|b| match b {
                b'A' => 'T',
                b'T' => 'A',
                b'C' => 'G',
                _ => 'C',
            })
            .collect();
        pairs.push(PrimerPair {
            forward: DesignedPrimer {
                sequence: template[f_start..f_start + f_len].into(),
                binding: f_start..f_start + f_len,
                reverse: false,
                tm_c: pair.forward_tm,
                gc_percent: pair.forward_gc,
                hairpin_tm_c: pair.forward_hairpin,
            },
            reverse: DesignedPrimer {
                sequence: reverse_sequence,
                binding: r_end - r_len..r_end,
                reverse: true,
                tm_c: pair.reverse_tm,
                gc_percent: pair.reverse_gc,
                hairpin_tm_c: pair.reverse_hairpin,
            },
            product: f_start..r_end,
            penalty: pair.penalty,
            heterodimer_tm_c: pair.pair_any,
            end_heterodimer_tm_c: pair.pair_end,
        });
    }
    let mut warnings = vec!["Specificity screening is limited to the supplied linear template; no genome or repeat library was supplied".into()];
    let native_warning = message(&warning);
    if !native_warning.is_empty() {
        warnings.push(native_warning);
    }
    if pairs.is_empty() {
        warnings.push("No primer pairs satisfy the supplied constraints".into());
    }
    if options.reaction.dmso_percent != 0. || options.reaction.formamide_m != 0. {
        warnings.push(
            "Additives correct duplex Tm only; hairpin/dimer screening does not model them".into(),
        );
    }
    Ok(DesignResult { pairs, warnings })
}

#[repr(C)]
struct NativeReaction {
    dna_nm: f64,
    monovalent_mm: f64,
    magnesium_mm: f64,
    dntp_mm: f64,
    dmso_percent: f64,
    formamide_m: f64,
    temperature_c: f64,
    tm_method: c_int,
    salt_method: c_int,
}
#[derive(Default)]
#[repr(C)]
struct NativeStructure {
    tm: f64,
    dh: f64,
    ds: f64,
    dg: f64,
    found: c_int,
}
#[repr(C)]
struct NativeDesign {
    reaction: NativeReaction,
    min_tm: f64,
    optimal_tm: f64,
    max_tm: f64,
    min_gc: f64,
    max_gc: f64,
    max_tm_difference: f64,
    min_length: c_int,
    optimal_length: c_int,
    max_length: c_int,
    product_min: c_int,
    product_max: c_int,
    num_return: c_int,
    target_start: c_int,
    target_length: c_int,
}
#[derive(Clone, Copy, Default)]
#[repr(C)]
struct NativePair {
    forward_start: c_int,
    forward_length: c_int,
    reverse_end: c_int,
    reverse_length: c_int,
    product_size: c_int,
    forward_tm: f64,
    reverse_tm: f64,
    forward_gc: f64,
    reverse_gc: f64,
    penalty: f64,
    forward_hairpin: f64,
    reverse_hairpin: f64,
    pair_any: f64,
    pair_end: f64,
}
unsafe extern "C" {
    fn pg_primer_tm(sequence: *const c_char, reaction: *const NativeReaction) -> f64;
    fn pg_primer_structure(
        first: *const c_char,
        second: *const c_char,
        reaction: *const NativeReaction,
        hairpin: c_int,
        output: *mut NativeStructure,
        error: *mut c_char,
        capacity: usize,
    ) -> c_int;
    fn pg_design_primers(
        template: *const c_char,
        settings: *const NativeDesign,
        pairs: *mut NativePair,
        capacity: c_int,
        count: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
        warning: *mut c_char,
        warning_capacity: usize,
    ) -> c_int;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pinned_primer3_thermodynamics_and_structure_reference() {
        let result = analyze("CGTTCCAAAGATGTGGGCATGAGCTTAC", Reaction::default()).unwrap();
        assert!((result.tm_c - 66.30992381310529).abs() < 1e-8);
        assert!((result.hairpin.tm_c - 41.69969027913402).abs() < 1e-8);
        assert!(result.hairpin.found);
        assert!(
            (result.homodimer.enthalpy_cal_mol + 38400.).abs() < 1e-6,
            "{:?}",
            result.homodimer
        );
        let complementary = heterodimer(
            "CGTTCCAAAGATGTGGGCATGAGCTTAC",
            "GTAAGCTCATGCCCACATCTTTGGAACG",
            Reaction::default(),
        )
        .unwrap();
        assert!(complementary.tm_c > result.homodimer.tm_c + 60.);
    }

    #[test]
    fn design_reproduces_pinned_reference_and_valid_strand_coordinates() {
        let sequence =
            include_str!("../../../validation/fixtures/primer-design-template.txt").trim();
        let options = DesignOptions {
            target: Some(200..250),
            ..DesignOptions::default()
        };
        let result = design_pairs(sequence, &options).unwrap();
        assert_eq!(result.pairs.len(), 5);
        let pair = &result.pairs[0];
        assert_eq!(pair.forward.sequence, "ATACACTGCCGCAAGTTGGA");
        assert_eq!(pair.reverse.sequence, "GGGCAAGCTGAGATGTCCTT");
        assert_eq!(pair.forward.binding, 139..159);
        assert_eq!(pair.reverse.binding, 394..414);
        assert_eq!(pair.product, 139..414);
        assert!((pair.penalty - 0.07083419805553604).abs() < 1e-8);
        for pair in &result.pairs {
            assert!(pair.forward.binding.end <= 200 && pair.reverse.binding.start >= 250);
            assert!((100..301).contains(&pair.product.len()));
            assert!(pair.forward.tm_c >= 57. && pair.reverse.tm_c <= 63.);
        }
    }

    #[test]
    fn circular_target_crosses_origin_without_losing_binding_or_product_coordinates() {
        let original =
            include_str!("../../../validation/fixtures/primer-design-template.txt").trim();
        let sequence = format!("{}{}", &original[225..], &original[..225]);
        let n = sequence.len();
        let target = n - 25..n + 25;
        let result = design_pairs_circular(
            &sequence,
            &DesignOptions {
                target: Some(target.clone()),
                ..DesignOptions::default()
            },
        )
        .unwrap();
        assert_eq!(result.pairs.len(), 5);
        for pair in &result.pairs {
            assert!(pair.product.end > n);
            let relative = (target.start + n - pair.product.start) % n;
            assert!(relative + target.len() <= pair.product.len());
            assert_eq!(pair.product.start, pair.forward.binding.start);
            assert_eq!(pair.product.end % n, pair.reverse.binding.end % n);
            for p in [&pair.forward, &pair.reverse] {
                let parts = p.binding_segments(n).unwrap();
                let bound = parts
                    .iter()
                    .map(|s| &sequence[s.start..s.end])
                    .collect::<String>();
                let oligo = if p.reverse {
                    bound
                        .chars()
                        .rev()
                        .map(|b| match b {
                            'A' => 'T',
                            'T' => 'A',
                            'C' => 'G',
                            _ => 'C',
                        })
                        .collect()
                } else {
                    bound
                };
                assert_eq!(oligo, p.sequence);
            }
        }
        assert!(
            design_pairs_circular(
                &sequence,
                &DesignOptions {
                    num_return: 0,
                    ..DesignOptions::default()
                }
            )
            .is_err()
        );
    }
    #[test]
    fn wrapped_binding_segments_preserve_molecule_order() {
        let primer = DesignedPrimer {
            sequence: "AACCGG".into(),
            binding: 8..14,
            reverse: true,
            tm_c: 0.,
            gc_percent: 0.,
            hairpin_tm_c: 0.,
        };
        assert_eq!(
            primer.binding_segments(10).unwrap(),
            vec![
                crate::Segment { start: 8, end: 10 },
                crate::Segment { start: 0, end: 4 }
            ]
        );
        assert!(primer.binding_segments(5).is_err());
    }
    #[test]
    fn invalid_inputs_do_not_reach_native_code() {
        for sequence in ["", "A", "ACGN", "AC GU", "AUCG"] {
            assert!(analyze(sequence, Reaction::default()).is_err());
        }
        assert!(melting_temperature(&"A".repeat(61), Reaction::default()).is_err());
        let bad = Reaction {
            dna_nm: f64::NAN,
            ..Reaction::default()
        };
        assert!(analyze("ACGT", bad).is_err());
        let zero_salt = Reaction {
            monovalent_mm: 0.,
            magnesium_mm: 0.,
            ..Reaction::default()
        };
        assert!(analyze("ACGT", zero_salt).is_err());
        let bad_target = DesignOptions {
            target: Some(10..20),
            ..DesignOptions::default()
        };
        assert!(design_pairs("ACGT", &bad_target).is_err());
        let impossible = design_pairs(&"A".repeat(400), &DesignOptions::default()).unwrap();
        assert!(impossible.pairs.is_empty());
        assert!(
            impossible
                .warnings
                .iter()
                .any(|w| w.contains("No primer pairs"))
        );
    }
}
