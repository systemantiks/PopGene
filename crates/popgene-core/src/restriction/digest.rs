//! Geometric complete-digest products and explicit methylation assumptions.
//! Reference: https://www.neb.com/en/tools-and-resources/selection-charts/dam-dcm-and-cpg-methylation
//! Sensitivity is a guide, not a kinetic model. Unknown/conditional cases remain
//! unresolved; no cleavage efficiency or partial-digest probability is invented.
use super::{ENZYMES, RestrictionSite, analyze_selected};
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Methylation {
    /// Fully methylated GATC motifs on both strands (not hemimethylation).
    pub dam: bool,
    /// Fully methylated CCWGG motifs on both strands.
    pub dcm: bool,
    /// Every CpG is methylated on both strands; not a tissue-specific methylome.
    pub cpg: bool,
}
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct DigestConditions {
    pub temperature_c: f64,
    pub ph: f64,
    pub magnesium_mm: f64,
    pub glycerol_percent: f64,
    pub methylation: Methylation,
}
impl Default for DigestConditions {
    fn default() -> Self {
        Self {
            temperature_c: 37.,
            ph: 7.9,
            magnesium_mm: 10.,
            glycerol_percent: 5.,
            methylation: Methylation::default(),
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CleavageAssessment {
    Permitted,
    Blocked(String),
    Unresolved(String),
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum EndKind {
    Native,
    Blunt,
    FivePrime,
    ThreePrime,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct FragmentEnd {
    pub kind: EndKind,
    /// Actual protruding strand, written 5′→3′. Complementary fragment ends
    /// therefore carry reverse-complementary strings.
    pub overhang: String,
    pub enzymes: Vec<String>,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct DigestFragment {
    /// Forward-strand coordinates; end may be below start when crossing origin.
    pub start: usize,
    pub end: usize,
    pub wraps_origin: bool,
    pub circular: bool,
    pub sequence: String,
    /// Complementary strand written 5′→3′, including its own staggered ends.
    pub opposite_sequence: String,
    pub left_end: Option<FragmentEnd>,
    pub right_end: Option<FragmentEnd>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub left_modification: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub right_modification: Option<String>,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AssessedSite {
    pub site: RestrictionSite,
    pub assessment: CleavageAssessment,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DigestResult {
    pub fragments: Vec<DigestFragment>,
    pub sites: Vec<AssessedSite>,
    /// False if any recognized site cannot be resolved by the chosen model.
    pub all_sites_resolved: bool,
    pub warnings: Vec<String>,
}

fn overlaps_modification(
    sequence: &[u8],
    circular: bool,
    site: &RestrictionSite,
    motif_length: usize,
    kind: usize,
) -> bool {
    let n = sequence.len() as isize;
    let (pattern, modified): (&[u8], &[usize]) = match kind {
        0 => (b"GATC", &[1, 2]),
        1 => (b"CCWGG", &[1, 3]),
        _ => (b"CG", &[0, 1]),
    };
    for offset in -(pattern.len() as isize - 1)..motif_length as isize {
        let mut matches = true;
        for (i, expected) in pattern.iter().enumerate() {
            let p = site.start as isize + offset + i as isize;
            if !circular && !(0..n).contains(&p) {
                matches = false;
                break;
            }
            let actual = sequence[p.rem_euclid(n) as usize].to_ascii_uppercase();
            if !(*expected == actual || *expected == b'W' && matches!(actual, b'A' | b'T')) {
                matches = false;
                break;
            }
        }
        if matches
            && modified
                .iter()
                .any(|p| (0..motif_length as isize).contains(&(offset + *p as isize)))
        {
            return true;
        }
    }
    false
}

/// Assessment for fully methylated Dam/Dcm/CpG profiles only. Supplier-specific
/// mutants, partial/hemimethylation and other modified bases are not inferred.
pub fn assess_site(
    sequence: &str,
    circular: bool,
    site: &RestrictionSite,
    methylation: Methylation,
) -> CleavageAssessment {
    let Some(enzyme) = ENZYMES.get(site.enzyme_index) else {
        return CleavageAssessment::Unresolved("Invalid enzyme index".into());
    };
    if !enzyme.cleavage_known {
        return CleavageAssessment::Unresolved(
            "Cleavage offsets are unspecified in the reference catalogue".into(),
        );
    }
    if sequence.is_empty() || site.start >= sequence.len() {
        return CleavageAssessment::Unresolved("Invalid recognition-site position".into());
    }
    if enzyme.name == "DpnI" && !methylation.dam {
        return CleavageAssessment::Blocked(
            "DpnI requires methylated adenine in GATC; the selected profile is Dam-negative".into(),
        );
    }
    if matches!(
        enzyme.name,
        "MspJI" | "FspEI" | "LpnPI" | "AspBHI" | "AbaSI" | "SgeI" | "BisI" | "GlaI"
    ) {
        return CleavageAssessment::Unresolved(
            "This modification-dependent enzyme requires a site-specific modified-base model"
                .into(),
        );
    }
    for (kind, active) in [methylation.dam, methylation.dcm, methylation.cpg]
        .into_iter()
        .enumerate()
    {
        if !active
            || !overlaps_modification(
                sequence.as_bytes(),
                circular,
                site,
                site.recognition_length,
                kind,
            )
        {
            continue;
        }
        let name = ["Dam", "Dcm", "CpG"][kind];
        let blocked = match kind {
            0 => matches!(enzyme.name, "AlwI" | "BclI" | "DpnII" | "MboI" | "HphI"),
            1 => matches!(enzyme.name, "HphI" | "PspGI" | "SexAI"),
            _ => matches!(
                enzyme.name,
                "AatII"
                    | "AciI"
                    | "AclI"
                    | "AfeI"
                    | "AgeI"
                    | "AscI"
                    | "AsiSI"
                    | "AvaI"
                    | "BceAI"
                    | "BmgBI"
                    | "BsaAI"
                    | "BsaHI"
                    | "BsiEI"
                    | "BsiWI"
                    | "BsmBI"
                    | "BsrFI"
                    | "BssHII"
                    | "BstBI"
                    | "BstUI"
                    | "ClaI"
                    | "EagI"
                    | "Esp3I"
                    | "FauI"
                    | "FseI"
                    | "FspI"
                    | "HaeII"
                    | "HgaI"
                    | "HhaI"
                    | "HinP1I"
                    | "HpaII"
                    | "HpyCH4IV"
                    | "Hpy99I"
                    | "KasI"
                    | "MluI"
                    | "NaeI"
                    | "NarI"
                    | "NgoMIV"
                    | "NotI"
                    | "NruI"
                    | "PaeR7I"
                    | "PmlI"
                    | "PvuI"
                    | "RsrII"
                    | "SacII"
                    | "SalI"
                    | "SmaI"
                    | "XhoI"
            ),
        };
        if blocked {
            return CleavageAssessment::Blocked(format!(
                "{name} methylation blocks this enzyme in the reference sensitivity table"
            ));
        }
        let insensitive = match kind {
            0 => matches!(
                enzyme.name,
                "DpnI"
                    | "Acc65I"
                    | "AcuI"
                    | "AflII"
                    | "AleI"
                    | "AluI"
                    | "ApaI"
                    | "AvrII"
                    | "BamHI"
                    | "BbsI"
                    | "BmrI"
                    | "BsaI"
                    | "BsgI"
                    | "BsmBI"
                    | "BspQI"
                    | "BssHII"
                    | "BstAPI"
                    | "BstEII"
                    | "BstZ17I"
                    | "Bsu36I"
                    | "DraI"
                    | "EagI"
                    | "EcoP15I"
                    | "EcoRI"
                    | "EcoRV"
                    | "HaeIII"
                    | "HindIII"
                    | "HpaII"
                    | "KpnI"
                    | "MspI"
                    | "NaeI"
                    | "NgoMIV"
                    | "NheI"
                    | "NmeAIII"
                    | "NotI"
                    | "PacI"
                    | "PaeR7I"
                    | "PaqCI"
                    | "PciI"
                    | "PmlI"
                    | "PshAI"
                    | "PspOMI"
                    | "PstI"
                    | "RsrII"
                    | "SalI"
                    | "SapI"
                    | "SmaI"
                    | "SphI"
                    | "SspI"
                    | "TaqI"
                    | "XhoI"
                    | "Sau3AI"
            ),
            1 => matches!(
                enzyme.name,
                "DpnI"
                    | "AcuI"
                    | "AflII"
                    | "AleI"
                    | "AluI"
                    | "AvrII"
                    | "BamHI"
                    | "BbsI"
                    | "BmrI"
                    | "BsgI"
                    | "BsmBI"
                    | "BspQI"
                    | "BssHII"
                    | "BstAPI"
                    | "BstEII"
                    | "BstZ17I"
                    | "Bsu36I"
                    | "DpnII"
                    | "DraI"
                    | "EagI"
                    | "EcoP15I"
                    | "EcoRI"
                    | "EcoRV"
                    | "HaeIII"
                    | "HindIII"
                    | "HpaII"
                    | "KpnI"
                    | "MspI"
                    | "NaeI"
                    | "NgoMIV"
                    | "NheI"
                    | "NmeAIII"
                    | "NotI"
                    | "PacI"
                    | "PaeR7I"
                    | "PaqCI"
                    | "PciI"
                    | "PmlI"
                    | "PshAI"
                    | "PstI"
                    | "RsrII"
                    | "SalI"
                    | "SapI"
                    | "SmaI"
                    | "SphI"
                    | "SspI"
                    | "TaqI"
                    | "XbaI"
                    | "XhoI"
                    | "MboI"
                    | "Sau3AI"
            ),
            _ => matches!(
                enzyme.name,
                "AcuI"
                    | "AflII"
                    | "AluI"
                    | "AvrII"
                    | "BamHI"
                    | "BbsI"
                    | "BmrI"
                    | "BsgI"
                    | "BspQI"
                    | "BstEII"
                    | "Bsu36I"
                    | "DpnII"
                    | "DraI"
                    | "HaeIII"
                    | "HindIII"
                    | "KpnI"
                    | "MspI"
                    | "NmeAIII"
                    | "PacI"
                    | "PciI"
                    | "PstI"
                    | "SapI"
                    | "SphI"
                    | "SspI"
                    | "XbaI"
                    | "MboI"
            ),
        };
        if !insensitive {
            return CleavageAssessment::Unresolved(format!(
                "{name} methylation overlaps this site; sensitivity is conditional or not curated here"
            ));
        }
    }
    CleavageAssessment::Permitted
}

#[derive(Clone)]
struct Cut {
    top: usize,
    bottom: usize,
    stagger: i32,
    enzymes: Vec<String>,
}
fn span(sequence: &str, start: usize, length: usize) -> String {
    let n = sequence.len();
    (0..length)
        .map(|i| sequence.as_bytes()[(start + i) % n] as char)
        .collect()
}
fn reverse_complement(sequence: &str) -> String {
    sequence
        .bytes()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => 'T',
            b'T' => 'A',
            b'C' => 'G',
            b'G' => 'C',
            b'R' => 'Y',
            b'Y' => 'R',
            b'K' => 'M',
            b'M' => 'K',
            b'B' => 'V',
            b'V' => 'B',
            b'D' => 'H',
            b'H' => 'D',
            b'S' => 'S',
            b'W' => 'W',
            _ => 'N',
        })
        .collect()
}
fn end(sequence: &str, cut: &Cut, left: bool) -> FragmentEnd {
    let (kind, start) = if cut.stagger > 0 {
        (EndKind::FivePrime, cut.top)
    } else if cut.stagger < 0 {
        (EndKind::ThreePrime, cut.bottom)
    } else {
        (EndKind::Blunt, cut.top)
    };
    let mut overhang = span(sequence, start, cut.stagger.unsigned_abs() as usize);
    if (cut.stagger > 0 && !left) || (cut.stagger < 0 && left) {
        overhang = reverse_complement(&overhang);
    }
    FragmentEnd {
        kind,
        overhang,
        enzymes: cut.enzymes.clone(),
    }
}

/// All permitted cleavage sites are applied simultaneously. Results include
/// both strand sequences and cohesive ends. Intersecting cuts that would form
/// isolated single strands are rejected rather than assigned fictitious ends.
pub fn digest(
    sequence: &str,
    circular: bool,
    enzyme_indices: &[usize],
    conditions: DigestConditions,
) -> Result<DigestResult, String> {
    digest_internal(sequence, circular, enzyme_indices, conditions, None)
}

/// Digest a duplex document, respecting its actual paired interval and existing
/// cohesive ends. Existing unpaired sequence cannot establish a duplex site.
/// Outer strand sequences, ends and terminal modifications are retained exactly.
pub fn digest_document(
    document: &crate::Document,
    enzyme_indices: &[usize],
    conditions: DigestConditions,
) -> Result<DigestResult, String> {
    if document.molecule != crate::Molecule::Dna
        || document.single_stranded
        || document.trace.is_some()
    {
        return Err("Digestion requires a double-stranded DNA document".into());
    }
    let Some(ends) = &document.duplex_ends else {
        return digest(
            &document.sequence,
            document.circular,
            enzyme_indices,
            conditions,
        );
    };
    if document.circular {
        return Err("A circular molecule cannot carry linear duplex ends".into());
    }
    ends.validate(&document.sequence)
        .map_err(|e| e.to_string())?;
    let (tp, ts, bp, bs) = ends.overhang_lengths();
    let n = document.sequence.len();
    let mut result = digest_internal(
        &document.sequence,
        false,
        enzyme_indices,
        conditions,
        Some(tp..n - ts),
    )?;
    let bottom_axis_start = tp as isize - bs as isize;
    let bottom_axis_end = (n - ts + bp) as isize;
    let bottom_axis_sequence = reverse_complement(&ends.opposite_sequence);
    let stagger = |end: &FragmentEnd| match end.kind {
        EndKind::FivePrime => end.overhang.len() as isize,
        EndKind::ThreePrime => -(end.overhang.len() as isize),
        _ => 0,
    };
    for fragment in &mut result.fragments {
        let left_native = fragment
            .left_end
            .as_ref()
            .is_some_and(|e| e.kind == EndKind::Native);
        let right_native = fragment
            .right_end
            .as_ref()
            .is_some_and(|e| e.kind == EndKind::Native);
        let bottom_start = if left_native {
            bottom_axis_start
        } else {
            fragment.start as isize + fragment.left_end.as_ref().map_or(0, stagger)
        };
        let bottom_end = if right_native {
            bottom_axis_end
        } else {
            fragment.end as isize + fragment.right_end.as_ref().map_or(0, stagger)
        };
        let (a, b) = (
            bottom_start - bottom_axis_start,
            bottom_end - bottom_axis_start,
        );
        if a < 0 || b <= a || b as usize > bottom_axis_sequence.len() {
            return Err("Cleavage intersects an existing single-stranded end; a duplex product cannot be established".into());
        }
        fragment.opposite_sequence =
            reverse_complement(&bottom_axis_sequence[a as usize..b as usize]);
        if left_native {
            fragment.left_end = Some(ends.left.clone());
            fragment.left_modification = ends.left_modification.clone();
        }
        if right_native {
            fragment.right_end = Some(ends.right.clone());
            fragment.right_modification = ends.right_modification.clone();
        }
        crate::workflow::DuplexEnds {
            left: fragment.left_end.clone().ok_or("Missing left end")?,
            right: fragment.right_end.clone().ok_or("Missing right end")?,
            opposite_sequence: fragment.opposite_sequence.clone(),
            left_modification: fragment.left_modification.clone(),
            right_modification: fragment.right_modification.clone(),
        }
        .validate(&fragment.sequence)
        .map_err(|e| e.to_string())?;
    }
    result.warnings.push("Existing single-stranded ends are excluded from duplex recognition. Original terminal geometry and modifications are retained on outer products.".into());
    Ok(result)
}

fn digest_internal(
    sequence: &str,
    circular: bool,
    enzyme_indices: &[usize],
    conditions: DigestConditions,
    paired: Option<std::ops::Range<usize>>,
) -> Result<DigestResult, String> {
    if sequence.is_empty()
        || sequence.len() > 10_000_000
        || !sequence
            .bytes()
            .all(|b| super::iupac_mask(b.to_ascii_uppercase()) != 0)
    {
        return Err("Digest requires 1–10,000,000 IUPAC DNA bases".into());
    }
    if enzyme_indices.len() > 64 || enzyme_indices.iter().any(|&i| i >= ENZYMES.len()) {
        return Err("Choose up to 64 valid enzymes".into());
    }
    for (v, lo, hi) in [
        (conditions.temperature_c, 0., 100.),
        (conditions.ph, 0., 14.),
        (conditions.magnesium_mm, 0., 100.),
        (conditions.glycerol_percent, 0., 100.),
    ] {
        if !v.is_finite() || !(lo..=hi).contains(&v) {
            return Err("Invalid digestion conditions".into());
        }
    }
    let analysis = analyze_selected(sequence, circular, enzyme_indices, 100_000, || false).unwrap();
    if analysis.truncated || analysis.total_sites > 100_000 {
        return Err("Digest exceeds 100,000 recognition sites; analyze a smaller region".into());
    }
    let mut warnings=vec!["Assumes complete cleavage at permitted sites. Buffer activity, enzyme units, incubation kinetics, flanking-base requirements and multi-site enzyme interactions do not establish experimental yield.".into()];
    if conditions.glycerol_percent > 5. {
        warnings.push(
            "Glycerol exceeds 5%; star activity may alter specificity and is not simulated".into(),
        );
    }
    if !(7.0..=8.5).contains(&conditions.ph) {
        warnings.push(
            "pH is outside common restriction-buffer conditions; activity is not established"
                .into(),
        );
    }
    for &i in enzyme_indices {
        if (conditions.temperature_c - ENZYMES[i].catalogue_temperature_c as f64).abs() > 0.5 {
            warnings.push(format!("{}: reference catalogue temperature is {} °C; activity at the requested temperature is not established",ENZYMES[i].name,ENZYMES[i].catalogue_temperature_c));
        }
    }
    let mut cuts = vec![];
    let mut sites = vec![];
    let mut all_sites_resolved = true;
    let paired_sequence = paired.as_ref().map(|r| {
        format!(
            "{}{}{}",
            "N".repeat(r.start),
            &sequence[r.clone()],
            "N".repeat(sequence.len() - r.end)
        )
    });
    for site in analysis.sites {
        let enzyme = &ENZYMES[site.enzyme_index];
        let mut assessment = assess_site(
            paired_sequence.as_deref().unwrap_or(sequence),
            circular,
            &site,
            conditions.methylation,
        );
        if let Some(paired) = &paired
            && (site.start < paired.start || site.start + site.recognition_length > paired.end)
        {
            assessment=CleavageAssessment::Blocked("Recognition overlaps an existing single-stranded end; this duplex model requires a fully paired recognition site".into());
        }
        if conditions.magnesium_mm == 0. && assessment == CleavageAssessment::Permitted {
            assessment = CleavageAssessment::Unresolved(
                "No Mg²⁺ supplied; cleavage activity is not established".into(),
            );
        }
        if assessment == CleavageAssessment::Permitted {
            let mut site_cuts = Vec::new();
            for (top, bottom, stagger) in [
                (
                    site.top_cut,
                    site.bottom_cut,
                    enzyme.bottom_cut - enzyme.top_cut,
                ),
                (
                    site.second_top_cut,
                    site.second_bottom_cut,
                    enzyme.second_cut.map_or(0, |(t, b)| b - t),
                ),
            ]
            .into_iter()
            .take(if enzyme.second_cut.is_some() { 2 } else { 1 })
            {
                if let (Some(top), Some(bottom)) = (top, bottom) {
                    if paired.as_ref().is_some_and(|r| {
                        top <= r.start || bottom <= r.start || top >= r.end || bottom >= r.end
                    }) {
                        assessment=CleavageAssessment::Unresolved("Cleavage reaches an existing unpaired end or paired-region boundary; a complete duplex cut is not established".into());
                        break;
                    }
                    if stagger.unsigned_abs() as usize >= sequence.len() {
                        assessment = CleavageAssessment::Unresolved(
                            "Cleavage stagger is as long as the molecule".into(),
                        );
                        break;
                    }
                    site_cuts.push(Cut {
                        top,
                        bottom,
                        stagger,
                        enzymes: vec![enzyme.name.into()],
                    });
                } else {
                    assessment = CleavageAssessment::Unresolved(
                        "Cleavage boundary lies outside the linear molecule".into(),
                    );
                }
            }
            // Do not infer partial activity at a site with an unresolved pair.
            if assessment == CleavageAssessment::Permitted {
                cuts.extend(site_cuts);
            }
        }
        if matches!(assessment, CleavageAssessment::Unresolved(_)) {
            all_sites_resolved = false;
        }
        sites.push(AssessedSite { site, assessment });
    }
    cuts.sort_by_key(|c| (c.top, c.bottom));
    let mut unique: Vec<Cut> = vec![];
    for cut in cuts {
        if let Some(previous) = unique.last_mut()
            && previous.top == cut.top
        {
            if previous.bottom != cut.bottom {
                return Err(
                    "Overlapping cuts at one top-strand boundary have incompatible staggers".into(),
                );
            }
            for name in cut.enzymes {
                if !previous.enzymes.contains(&name) {
                    previous.enzymes.push(name);
                }
            }
        } else {
            unique.push(cut);
        }
    }
    let native = || FragmentEnd {
        kind: EndKind::Native,
        overhang: String::new(),
        enzymes: vec![],
    };
    if unique.is_empty() {
        return Ok(DigestResult {
            fragments: vec![DigestFragment {
                start: 0,
                end: sequence.len(),
                wraps_origin: false,
                circular,
                sequence: sequence.into(),
                opposite_sequence: reverse_complement(sequence),
                left_end: (!circular).then(native),
                right_end: (!circular).then(native),
                left_modification: None,
                right_modification: None,
            }],
            sites,
            all_sites_resolved,
            warnings,
        });
    }
    let n = sequence.len();
    let mut fragments = vec![];
    let count = if circular {
        unique.len()
    } else {
        unique.len() + 1
    };
    for i in 0..count {
        let left = if circular {
            Some(&unique[i])
        } else {
            i.checked_sub(1).map(|j| &unique[j])
        };
        let right = if circular {
            Some(&unique[(i + 1) % unique.len()])
        } else {
            unique.get(i)
        };
        let start = left.map_or(0, |c| c.top);
        let stop = right.map_or(n, |c| c.top);
        let top_len = if circular {
            if stop > start {
                stop - start
            } else {
                n - start + stop
            }
        } else {
            stop - start
        };
        let bottom_len = top_len as i64 + right.map_or(0, |c| c.stagger as i64)
            - left.map_or(0, |c| c.stagger as i64);
        if top_len == 0 && bottom_len == 0 {
            continue;
        }
        if top_len == 0 {
            return Err("A terminal staggered cut would create a single-stranded fragment; use a longer sequence".into());
        }
        if bottom_len <= 0 || bottom_len > n as i64 {
            return Err("Overlapping staggered cuts would create single-stranded fragments; use separate digests".into());
        }
        let bottom = span(sequence, left.map_or(0, |c| c.bottom), bottom_len as usize);
        fragments.push(DigestFragment {
            start,
            end: stop,
            wraps_origin: circular && stop <= start,
            circular: false,
            sequence: span(sequence, start, top_len),
            opposite_sequence: reverse_complement(&bottom),
            left_end: Some(left.map_or_else(native, |cut| end(sequence, cut, true))),
            right_end: Some(right.map_or_else(native, |cut| end(sequence, cut, false))),
            left_modification: None,
            right_modification: None,
        });
    }
    Ok(DigestResult {
        fragments,
        sites,
        all_sites_resolved,
        warnings,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn enzyme(name: &str) -> usize {
        ENZYMES.iter().position(|e| e.name == name).unwrap()
    }
    #[test]
    fn digest_document_preserves_all_native_end_geometries_and_ignores_unpaired_sites() {
        use crate::{Document, workflow::DuplexEnds};
        let paired = "AAAGAATTCTTT";
        for left_kind in [EndKind::FivePrime, EndKind::ThreePrime] {
            for right_kind in [EndKind::FivePrime, EndKind::ThreePrime] {
                let left = FragmentEnd {
                    kind: left_kind,
                    overhang: "GAATTC".into(),
                    enzymes: vec!["original left".into()],
                };
                let right = FragmentEnd {
                    kind: right_kind,
                    overhang: "GAATTC".into(),
                    enzymes: vec!["original right".into()],
                };
                let top = format!(
                    "{}{}{}",
                    if left_kind == EndKind::FivePrime {
                        "GAATTC"
                    } else {
                        ""
                    },
                    paired,
                    if right_kind == EndKind::ThreePrime {
                        "GAATTC"
                    } else {
                        ""
                    }
                );
                let opposite = format!(
                    "{}{}{}",
                    if right_kind == EndKind::FivePrime {
                        "GAATTC"
                    } else {
                        ""
                    },
                    reverse_complement(paired),
                    if left_kind == EndKind::ThreePrime {
                        "GAATTC"
                    } else {
                        ""
                    }
                );
                let ends = DuplexEnds {
                    left: left.clone(),
                    right: right.clone(),
                    opposite_sequence: opposite.clone(),
                    left_modification: Some("5′ phosphate".into()),
                    right_modification: Some("explicit terminal label".into()),
                };
                ends.validate(&top).unwrap();
                let mut document = Document::new("cohesive substrate", &top).unwrap();
                document.duplex_ends = Some(ends);
                let result =
                    digest_document(&document, &[enzyme("EcoRI")], DigestConditions::default())
                        .unwrap();
                assert!(result.all_sites_resolved);
                assert_eq!(result.fragments.len(), 2);
                assert_eq!(result.fragments[0].left_end.as_ref(), Some(&left));
                assert_eq!(result.fragments[1].right_end.as_ref(), Some(&right));
                assert_eq!(
                    result.fragments[0].left_modification.as_deref(),
                    Some("5′ phosphate")
                );
                assert_eq!(
                    result.fragments[1].right_modification.as_deref(),
                    Some("explicit terminal label")
                );
                assert_eq!(
                    result
                        .fragments
                        .iter()
                        .map(|f| f.sequence.as_str())
                        .collect::<String>(),
                    top
                );
                assert_eq!(
                    result
                        .fragments
                        .iter()
                        .rev()
                        .map(|f| f.opposite_sequence.as_str())
                        .collect::<String>(),
                    opposite
                );
                let blocked = result
                    .sites
                    .iter()
                    .filter(|s| matches!(s.assessment, CleavageAssessment::Blocked(_)))
                    .count();
                assert_eq!(
                    blocked,
                    usize::from(left_kind == EndKind::FivePrime)
                        + usize::from(right_kind == EndKind::ThreePrime)
                );
                let intact = digest_document(&document, &[], DigestConditions::default()).unwrap();
                assert_eq!(intact.fragments[0].opposite_sequence, opposite);
            }
        }
    }
    #[test]
    fn linear_and_circular_staggered_products_preserve_both_strands() {
        let sequence = "AAAGAATTCTTT";
        let linear = digest(
            sequence,
            false,
            &[enzyme("EcoRI")],
            DigestConditions::default(),
        )
        .unwrap();
        assert_eq!(
            linear
                .fragments
                .iter()
                .map(|f| f.sequence.as_str())
                .collect::<Vec<_>>(),
            ["AAAG", "AATTCTTT"]
        );
        assert_eq!(
            linear
                .fragments
                .iter()
                .map(|f| f.opposite_sequence.as_str())
                .collect::<Vec<_>>(),
            ["AATTCTTT", "AAAG"]
        );
        assert_eq!(
            linear.fragments[0].right_end.as_ref().unwrap().overhang,
            "AATT"
        );
        let circle = digest(
            sequence,
            true,
            &[enzyme("EcoRI")],
            DigestConditions::default(),
        )
        .unwrap();
        assert_eq!(circle.fragments.len(), 1);
        assert_eq!(circle.fragments[0].sequence.len(), sequence.len());
        assert_eq!(
            circle.fragments[0].left_end.as_ref().unwrap().kind,
            EndKind::FivePrime
        );
    }
    #[test]
    fn methylation_changes_actual_digest_for_dpn_and_hpa_isoschizomers() {
        let seq = "AAAGATCCCCGGTTT";
        let dam = DigestConditions {
            methylation: Methylation {
                dam: true,
                ..Methylation::default()
            },
            ..DigestConditions::default()
        };
        assert_eq!(
            digest(seq, false, &[enzyme("DpnI")], dam)
                .unwrap()
                .fragments
                .len(),
            2
        );
        assert_eq!(
            digest(seq, false, &[enzyme("DpnII")], dam)
                .unwrap()
                .fragments
                .len(),
            1
        );
        assert_eq!(
            digest(seq, false, &[enzyme("DpnI")], DigestConditions::default())
                .unwrap()
                .fragments
                .len(),
            1
        );
        let cpg = DigestConditions {
            methylation: Methylation {
                cpg: true,
                ..Methylation::default()
            },
            ..DigestConditions::default()
        };
        assert_eq!(
            digest(seq, false, &[enzyme("HpaII")], cpg)
                .unwrap()
                .fragments
                .len(),
            1
        );
        assert_eq!(
            digest(seq, false, &[enzyme("MspI")], cpg)
                .unwrap()
                .fragments
                .len(),
            2
        );
    }
    #[test]
    fn unknown_cleavage_remains_unresolved_and_type_iib_cuts_twice() {
        let unknown = ENZYMES
            .iter()
            .position(|e| !e.cleavage_known && e.recognition == "GATC");
        if let Some(index) = unknown {
            assert!(
                !digest("AAAGATCTTT", false, &[index], DigestConditions::default())
                    .unwrap()
                    .all_sites_resolved
            );
        }
        let sequence = format!(
            "{}CGA{}TGC{}",
            "A".repeat(30),
            "A".repeat(6),
            "A".repeat(30)
        );
        let result = digest(
            &sequence,
            false,
            &[enzyme("BcgI")],
            DigestConditions::default(),
        )
        .unwrap();
        assert_eq!(result.fragments.len(), 3);
        assert_eq!(
            result
                .fragments
                .iter()
                .map(|f| f.sequence.len())
                .sum::<usize>(),
            sequence.len()
        );
    }
}
