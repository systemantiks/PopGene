//! Agarose gel display model for double-stranded DNA products.
//! Resolution ranges for standard agarose in 1X TAE and 1X TBE: Thermo Fisher
//! UltraPure Agarose protocol (retrieved 2026-10-31), also printed in the
//! Invitrogen 1 Kb Plus DNA Ladder manual (MAN0000898):
//! https://www.thermofisher.com/jp/ja/home/references/protocols/nucleic-acid-purification-and-analysis/dna-protocol/ultrapure-agarose.html
//! Positions place bands for comparison between lanes of the same gel; they are not
//! measured mobilities. Voltage, run time in minutes, stain, load and DNA topology
//! are not modeled. Equimolar lanes carry mass in proportion to fragment length.
use crate::restriction::{self, DigestConditions, ENZYMES, Methylation};
use crate::{Document, MAX_BASES, workflow::DuplexEnds};
use serde::{Deserialize, Serialize};

/// (agarose %, smallest bp, largest bp) resolved as linear dsDNA.
const TAE: &[(f64, f64, f64)] = &[
    (0.6, 1_000., 20_000.),
    (0.8, 800., 10_000.),
    (1.0, 400., 8_000.),
    (1.2, 300., 7_000.),
    (1.5, 200., 4_000.),
    (2.0, 100., 3_000.),
];
const TBE: &[(f64, f64, f64)] = &[
    (0.5, 1_000., 20_000.),
    (0.7, 800., 10_000.),
    (0.85, 400., 8_000.),
    (1.0, 300., 7_000.),
    (1.25, 200., 4_000.),
    (1.75, 100., 3_000.),
];
/// At run length 1 the largest and smallest tabulated lengths are drawn at these
/// fractions of the gel. Lengths beyond the range are compressed towards the ends.
const RANGE_TOP: f64 = 0.15;
const RANGE_BOTTOM: f64 = 0.85;
/// Fragments closer than this fraction of the gel are drawn as one band.
pub const RESOLUTION: f64 = 0.012;
pub const MAX_LANES: usize = 12;
pub const MAX_LADDER_SIZES: usize = 64;
/// Relative electrophoresis time. Distances scale with it; nothing else changes.
pub const RUN_LENGTH: std::ops::RangeInclusive<f64> = 0.5..=2.;
const MAX_COMIGRATION_WARNINGS: usize = 8;
pub const MODEL_NOTE: &str = "Band positions follow a logistic in log10(length) through the tabulated resolution range of standard agarose in the chosen buffer. They compare lanes on one gel; they are not measured mobilities. Voltage, stain, load and DNA topology are not modeled.";
pub const GEL_FORMAT: &str = "popgene-gel";
pub const GEL_VERSION: u32 = 1;
pub const GEL_EXTENSION: &str = "popgel";

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum Buffer {
    #[default]
    Tae,
    Tbe,
}
impl Buffer {
    pub fn label(self) -> &'static str {
        match self {
            Self::Tae => "1X TAE",
            Self::Tbe => "1X TBE",
        }
    }
    fn table(self) -> &'static [(f64, f64, f64)] {
        match self {
            Self::Tae => TAE,
            Self::Tbe => TBE,
        }
    }
    /// Lowest and highest tabulated agarose percentages.
    pub fn agarose_range(self) -> (f64, f64) {
        let table = self.table();
        (table[0].0, table[table.len() - 1].0)
    }
}

/// Which part of a drawn band sits at the calculated position.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum BandEdge {
    Top,
    #[default]
    Middle,
    Bottom,
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct GelSettings {
    pub agarose_percent: f64,
    #[serde(default)]
    pub buffer: Buffer,
    #[serde(default = "unit_run")]
    pub run_length: f64,
    #[serde(default)]
    pub band_edge: BandEdge,
}
fn unit_run() -> f64 {
    1.
}
impl Default for GelSettings {
    fn default() -> Self {
        Self {
            agarose_percent: 1.,
            buffer: Buffer::Tae,
            run_length: 1.,
            band_edge: BandEdge::Middle,
        }
    }
}
impl GelSettings {
    /// The resolution range in bp for these settings.
    pub fn validate(&self) -> Result<(f64, f64), String> {
        if !self.run_length.is_finite() || !RUN_LENGTH.contains(&self.run_length) {
            return Err(format!(
                "Run length must be between {}× and {}×",
                RUN_LENGTH.start(),
                RUN_LENGTH.end()
            ));
        }
        separation_range(self.agarose_percent, self.buffer)
    }
}

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

/// Fragment lists are quoted from the manufacturers' product manuals. Ladders
/// still to be sourced are listed in docs/gel-simulation.md.
pub const LADDERS: &[Ladder] = &[
    Ladder {
        name: "GeneRuler 1 kb",
        sizes: &[
            10_000, 8_000, 6_000, 5_000, 4_000, 3_500, 3_000, 2_500, 2_000, 1_500, 1_000, 750, 500,
            250,
        ],
        reference: &[6_000, 3_000, 1_000],
        equimolar: false,
        note: "Fourteen purified fragments. Band masses other than the reference bands are not modeled.",
        source: "https://documents.thermofisher.com/TFS-Assets/LSG/manuals/MAN0013004_GeneRuler_1kb_DNALadder_250ug_UG.pdf",
    },
    Ladder {
        name: "GeneRuler 1 kb Plus",
        sizes: &[
            20_000, 10_000, 7_000, 5_000, 4_000, 3_000, 2_000, 1_500, 1_000, 700, 500, 400, 300,
            200, 75,
        ],
        reference: &[5_000, 1_500, 500],
        equimolar: false,
        note: "Fifteen purified fragments. Band masses other than the reference bands are not modeled.",
        source: "https://documents.thermofisher.com/TFS-Assets/LSG/manuals/MAN0013047_GeneRuler_1kb_Plus_DNALadder_250ug_UG.pdf",
    },
    Ladder {
        name: "GeneRuler 50 bp",
        sizes: &[
            1_000, 900, 800, 700, 600, 500, 400, 300, 250, 200, 150, 100, 50,
        ],
        reference: &[500, 250],
        equimolar: false,
        note: "Thirteen purified fragments. Band masses other than the reference bands are not modeled.",
        source: "https://documents.thermofisher.com/TFS-Assets/LSG/manuals/MAN0013016_GeneRuler_50bp_DNALadder_50ug_UG.pdf",
    },
    Ladder {
        name: "Invitrogen 1 Kb Plus",
        sizes: &[
            15_000, 10_000, 8_000, 7_000, 6_000, 5_000, 4_000, 3_000, 2_000, 1_500, 1_000, 850,
            650, 500, 400, 300, 200, 100,
        ],
        reference: &[1_500],
        equimolar: false,
        note: "Eighteen purified fragments from manual MAN0000898 Rev. A.0; earlier versions of this product had a different band set.",
        source: "https://documents.thermofisher.com/TFS-Assets/LSG/manuals/1Kb_Plus_DNA_ladder_man.pdf",
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
pub fn ladder(name: &str) -> Option<&'static Ladder> {
    LADDERS.iter().find(|l| l.name == name)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LaneFragment {
    pub length: usize,
    pub circular: bool,
    pub reference: bool,
}
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Lane {
    pub label: String,
    pub fragments: Vec<LaneFragment>,
    pub equimolar: bool,
    pub ladder: bool,
    /// Requested enzymes without a recognition site in this sample.
    pub noncutters: Vec<String>,
    /// Some recognized sites could not be resolved by the digestion model.
    pub provisional: bool,
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
            ladder: true,
            ..Default::default()
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
            ladder: true,
            ..Default::default()
        }
    }
    /// Complete-digest products are equimolar. Sizes are top-strand lengths.
    pub fn digest(label: impl Into<String>, fragments: &[restriction::DigestFragment]) -> Self {
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
            ..Default::default()
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Band {
    /// Lengths of every fragment drawn in this band, largest first.
    pub lengths: Vec<usize>,
    /// Indices into the lane's fragments, in the same order as `lengths`.
    pub members: Vec<usize>,
    /// Gel fractions from the well of the slowest and fastest member.
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
    pub ladder: bool,
    pub noncutters: Vec<String>,
    /// Lengths that have run past the end of the gel.
    pub ran_off: Vec<usize>,
}
#[derive(Clone, Debug, PartialEq)]
pub struct Gel {
    pub settings: GelSettings,
    /// Resolution range for this gel, in bp.
    pub separation: (usize, usize),
    pub lanes: Vec<GelLane>,
    pub warnings: Vec<String>,
}

/// Bounds between tabulated concentrations are interpolated geometrically.
pub fn separation_range(percent: f64, buffer: Buffer) -> Result<(f64, f64), String> {
    let table = buffer.table();
    let (first, last) = buffer.agarose_range();
    if !percent.is_finite() || !(first..=last).contains(&percent) {
        return Err(format!(
            "Agarose in {} must be between {first}% and {last}%; resolution outside the tabulated standard-agarose range is not modeled",
            buffer.label()
        ));
    }
    let i = table
        .windows(2)
        .position(|w| percent <= w[1].0)
        .unwrap_or(0);
    let (a, b) = (table[i], table[i + 1]);
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

/// Displayed gel fraction from the well (0); beyond 1 the fragment has run off.
pub fn migration(length: usize, settings: &GelSettings) -> Result<f64, String> {
    if length == 0 {
        return Err("A fragment must contain at least one base pair".into());
    }
    let (small, large) = settings.validate()?;
    Ok(position(length, small, large) * settings.run_length)
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
pub fn simulate(lanes: &[Lane], settings: &GelSettings) -> Result<Gel, String> {
    if lanes.is_empty() || lanes.len() > MAX_LANES {
        return Err(format!("A gel needs between 1 and {MAX_LANES} lanes"));
    }
    let (small, large) = settings.validate()?;
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
            .enumerate()
            .map(|(i, f)| (position(f.length, small, large) * settings.run_length, i, f))
            .collect::<Vec<_>>();
        placed.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut bands: Vec<Band> = vec![];
        let mut ran_off = vec![];
        for (p, i, f) in placed {
            if p > 1. {
                ran_off.push(f.length);
                continue;
            }
            let mass = if lane.equimolar { f.length as f64 } else { 1. };
            match bands.last_mut() {
                Some(band) if band.circular == f.circular && p - band.bottom < RESOLUTION => {
                    band.lengths.push(f.length);
                    band.members.push(i);
                    band.bottom = p;
                    band.intensity += mass;
                    band.reference |= f.reference;
                }
                _ => bands.push(Band {
                    lengths: vec![f.length],
                    members: vec![i],
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
        if lane.provisional {
            warnings.push(format!(
                "{}: products are provisional; some recognized sites could not be resolved under the chosen digestion assumptions",
                lane.label
            ));
        }
        let outside = lane
            .fragments
            .iter()
            .filter(|f| !f.circular && !(small..=large).contains(&(f.length as f64)))
            .count();
        if outside > 0 {
            warnings.push(format!(
                "{}: {outside} fragment{} outside the {}–{} bp resolution range of {}% agarose in {}; positions there are compressed and poorly resolved",
                lane.label,
                if outside == 1 { " lies" } else { "s lie" },
                format_bp(small.round() as usize),
                format_bp(large.round() as usize),
                settings.agarose_percent,
                settings.buffer.label()
            ));
        }
        if !ran_off.is_empty() {
            warnings.push(format!(
                "{}: {} fragment{} ran off the end of the gel; shorten the run to keep {}",
                lane.label,
                ran_off.len(),
                if ran_off.len() == 1 { "" } else { "s" },
                if ran_off.len() == 1 { "it" } else { "them" }
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
            ladder: lane.ladder,
            noncutters: lane.noncutters.clone(),
            ran_off,
        });
    }
    Ok(Gel {
        settings: *settings,
        separation: (small.round() as usize, large.round() as usize),
        lanes: output,
        warnings,
    })
}

/// The molecule loaded in a sample lane, independent of any open document.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct SampleDocument {
    pub name: String,
    pub sequence: String,
    pub circular: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub duplex_ends: Option<DuplexEnds>,
}
impl SampleDocument {
    pub fn from_document(document: &Document) -> Result<Self, String> {
        if document.molecule != crate::Molecule::Dna
            || document.single_stranded
            || document.trace.is_some()
            || document.placeholder
        {
            return Err(format!(
                "{} is not double-stranded DNA; only dsDNA can be loaded on this gel",
                document.name
            ));
        }
        Ok(Self {
            name: document.name.clone(),
            sequence: document.sequence.to_string(),
            circular: document.circular,
            duplex_ends: document.duplex_ends.clone(),
        })
    }
    pub fn document(&self) -> Result<Document, String> {
        let mut document = Document::new(self.name.clone(), &self.sequence)
            .map_err(|e| format!("{}: {e}", self.name))?;
        document.circular = self.circular;
        document.duplex_ends = self.duplex_ends.clone();
        if let Some(ends) = &document.duplex_ends {
            if document.circular {
                return Err(format!(
                    "{}: a circular molecule cannot carry linear duplex ends",
                    self.name
                ));
            }
            ends.validate(&document.sequence)
                .map_err(|e| format!("{}: {e}", self.name))?;
        }
        Ok(document)
    }
}

/// A lane as the user configured it; computing it may require a digestion.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum LaneSpec {
    Ladder {
        name: String,
    },
    CustomLadder {
        label: String,
        sizes: Vec<usize>,
    },
    Sample {
        label: String,
        sample: Box<SampleDocument>,
        /// Empty for uncut DNA.
        #[serde(default)]
        enzymes: Vec<String>,
        #[serde(default)]
        methylation: Methylation,
    },
}
impl LaneSpec {
    pub fn label(&self) -> &str {
        match self {
            Self::Ladder { name } => name,
            Self::CustomLadder { label, .. } | Self::Sample { label, .. } => label,
        }
    }
    pub fn is_ladder(&self) -> bool {
        !matches!(self, Self::Sample { .. })
    }
    pub fn enzyme_indices(names: &[String]) -> Result<Vec<usize>, String> {
        let mut indices = vec![];
        for name in names {
            let i = ENZYMES
                .iter()
                .position(|e| e.name.eq_ignore_ascii_case(name))
                .ok_or_else(|| format!("Unknown enzyme: {name}"))?;
            if !indices.contains(&i) {
                indices.push(i);
            }
        }
        if indices.len() > 64 {
            return Err("Choose at most 64 enzymes per lane".into());
        }
        Ok(indices)
    }
    /// Validates the lane without running a digestion.
    pub fn validate(&self) -> Result<(), String> {
        match self {
            Self::Ladder { name } => ladder(name).map(|_| ()).ok_or_else(|| {
                format!("This version of PopGene does not include the {name} ladder")
            }),
            Self::CustomLadder { label, sizes } => {
                if sizes.is_empty() || sizes.len() > MAX_LADDER_SIZES {
                    return Err(format!(
                        "{label}: enter between 1 and {MAX_LADDER_SIZES} ladder sizes"
                    ));
                }
                if sizes.iter().any(|s| !(1..=MAX_BASES).contains(s)) {
                    return Err(format!(
                        "{label}: ladder sizes must be from 1 to {} bp",
                        format_bp(MAX_BASES)
                    ));
                }
                Ok(())
            }
            Self::Sample {
                sample, enzymes, ..
            } => {
                sample.document()?;
                Self::enzyme_indices(enzymes).map(|_| ())
            }
        }
    }
    /// Digestion uses the default reaction with this lane's methylation assumptions.
    pub fn lane(&self) -> Result<Lane, String> {
        self.validate()?;
        Ok(match self {
            Self::Ladder { name } => Lane::ladder(ladder(name).unwrap()),
            Self::CustomLadder { label, sizes } => Lane::custom_ladder(label.clone(), sizes),
            Self::Sample {
                label,
                sample,
                enzymes,
                methylation,
            } => {
                let document = sample.document()?;
                let indices = Self::enzyme_indices(enzymes)?;
                if indices.is_empty() {
                    return Ok(Lane {
                        label: label.clone(),
                        fragments: vec![LaneFragment {
                            length: document.sequence.len(),
                            circular: document.circular,
                            reference: false,
                        }],
                        equimolar: true,
                        ..Default::default()
                    });
                }
                let conditions = DigestConditions {
                    methylation: *methylation,
                    ..Default::default()
                };
                let result = restriction::digest_document(&document, &indices, conditions)
                    .map_err(|e| format!("{label}: {e}"))?;
                let mut lane = Lane::digest(label.clone(), &result.fragments);
                lane.noncutters = indices
                    .iter()
                    .filter(|&&i| !result.sites.iter().any(|s| s.site.enzyme_index == i))
                    .map(|&i| ENZYMES[i].name.to_string())
                    .collect();
                lane.provisional = !result.all_sites_resolved;
                lane
            }
        })
    }
}

/// A saved gel: settings and lane definitions, including every loaded sequence.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct GelSpec {
    pub format: String,
    pub version: u32,
    pub settings: GelSettings,
    pub lanes: Vec<LaneSpec>,
}
impl GelSpec {
    pub fn new(settings: GelSettings, lanes: Vec<LaneSpec>) -> Self {
        Self {
            format: GEL_FORMAT.into(),
            version: GEL_VERSION,
            settings,
            lanes,
        }
    }
    pub fn validate(&self) -> Result<(), String> {
        if self.format != GEL_FORMAT {
            return Err("This file is not a PopGene gel".into());
        }
        if self.version != GEL_VERSION {
            return Err(format!(
                "Gel file version {} is not supported by this version of PopGene",
                self.version
            ));
        }
        if self.lanes.is_empty() || self.lanes.len() > MAX_LANES {
            return Err(format!("A gel needs between 1 and {MAX_LANES} lanes"));
        }
        self.settings.validate()?;
        self.lanes.iter().try_for_each(LaneSpec::validate)
    }
    pub fn to_json(&self) -> Result<String, String> {
        self.validate()?;
        serde_json::to_string_pretty(self).map_err(|e| e.to_string())
    }
    pub fn from_json(text: &str) -> Result<Self, String> {
        let spec: Self =
            serde_json::from_str(text).map_err(|e| format!("The gel file is malformed: {e}"))?;
        spec.validate()?;
        Ok(spec)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::restriction::digest;

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
            ..Default::default()
        }
    }
    fn settings(agarose_percent: f64) -> GelSettings {
        GelSettings {
            agarose_percent,
            ..Default::default()
        }
    }
    fn sample(sequence: &str, circular: bool, enzymes: &[&str]) -> LaneSpec {
        LaneSpec::Sample {
            label: "Plasmid".into(),
            sample: Box::new(SampleDocument {
                name: "Plasmid".into(),
                sequence: sequence.into(),
                circular,
                duplex_ends: None,
            }),
            enzymes: enzymes.iter().map(|e| e.to_string()).collect(),
            methylation: Methylation::default(),
        }
    }
    fn plasmid() -> String {
        let mut sequence = "A".repeat(3_000);
        sequence.replace_range(100..106, "GAATTC");
        sequence.replace_range(1_100..1_106, "GAATTC");
        sequence
    }

    #[test]
    fn resolution_follows_buffer_tables_and_refuses_extrapolation() {
        for (buffer, table) in [(Buffer::Tae, TAE), (Buffer::Tbe, TBE)] {
            for &(percent, small, large) in table {
                let (s, l) = separation_range(percent, buffer).unwrap();
                assert!(close(s, small) && close(l, large), "{percent}%");
            }
        }
        let (s, l) = separation_range(1.1, Buffer::Tae).unwrap();
        assert!(s < 400. && s > 300. && l < 8_000. && l > 7_000.);
        // The same agarose percentage resolves smaller fragments in TBE.
        assert!(
            separation_range(1.0, Buffer::Tbe).unwrap().0
                < separation_range(1.0, Buffer::Tae).unwrap().0
        );
        assert!(separation_range(2.0, Buffer::Tae).is_ok());
        for (percent, buffer) in [
            (0.5, Buffer::Tae),
            (2.1, Buffer::Tae),
            (2.0, Buffer::Tbe),
            (0.4, Buffer::Tbe),
            (f64::NAN, Buffer::Tae),
        ] {
            assert!(separation_range(percent, buffer).is_err(), "{percent}");
        }
    }

    #[test]
    fn tabulated_limits_are_drawn_at_fixed_fractions_and_scale_with_run_length() {
        let gel = settings(1.0);
        assert!((migration(8_000, &gel).unwrap() - RANGE_TOP).abs() < 1e-9);
        assert!((migration(400, &gel).unwrap() - RANGE_BOTTOM).abs() < 1e-9);
        let mut previous = 0.;
        for length in [50_000, 10_000, 3_000, 1_000, 500, 100, 10] {
            let p = migration(length, &gel).unwrap();
            assert!(p > previous && p < 1., "{length}");
            previous = p;
        }
        assert!(
            migration(1_000, &settings(0.8)).unwrap() > migration(1_000, &settings(1.5)).unwrap()
        );
        let longer = GelSettings {
            run_length: 1.5,
            ..gel
        };
        assert!(close(
            migration(1_000, &longer).unwrap(),
            1.5 * migration(1_000, &gel).unwrap()
        ));
        for run_length in [0.4, 2.1, f64::NAN] {
            assert!(migration(1_000, &GelSettings { run_length, ..gel }).is_err());
        }
        assert!(migration(0, &gel).is_err());
    }

    #[test]
    fn long_runs_lose_small_fragments_off_the_gel() {
        let long = GelSettings {
            run_length: 2.,
            ..settings(1.0)
        };
        let gel = simulate(&[fragments(&[5_000, 1_000, 200])], &long).unwrap();
        assert_eq!(gel.lanes[0].ran_off, [1_000, 200]);
        assert_eq!(gel.lanes[0].bands.len(), 1);
        assert!(
            gel.warnings
                .iter()
                .any(|w| w.contains("2 fragments ran off"))
        );
    }

    #[test]
    fn equimolar_lanes_scale_intensity_with_length() {
        let gel = simulate(&[fragments(&[4_000, 1_000])], &settings(1.0)).unwrap();
        let bands = &gel.lanes[0].bands;
        assert_eq!(bands[0].lengths, [4_000]);
        assert_eq!(bands[1].members, [1]);
        assert!(close(bands[0].intensity, 1.) && close(bands[1].intensity, 0.25));
        assert!(gel.warnings.is_empty());
        let ladder = simulate(
            &[Lane::custom_ladder("Custom", &[4_000, 1_000])],
            &settings(1.0),
        )
        .unwrap();
        assert!(ladder.lanes[0].ladder);
        assert!(ladder.lanes[0].bands.iter().all(|b| close(b.intensity, 1.)));
    }

    #[test]
    fn unresolved_and_identical_fragments_merge_with_explanations() {
        let gel = simulate(
            &[fragments(&[3_000, 1_010, 1_000, 500, 500])],
            &settings(1.0),
        )
        .unwrap();
        let lengths = gel.lanes[0]
            .bands
            .iter()
            .map(|b| b.lengths.clone())
            .collect::<Vec<_>>();
        assert_eq!(lengths, [vec![3_000], vec![1_010, 1_000], vec![500, 500]]);
        let merged = &gel.lanes[0].bands[1];
        assert!(merged.bottom > merged.top);
        assert_eq!(merged.members, [1, 2]);
        assert!(gel.warnings.iter().any(|w| w.contains("1,000–1,010 bp")));
        assert!(gel.warnings.iter().any(|w| w.contains("2 × 500 bp")));
    }

    #[test]
    fn out_of_range_and_circular_molecules_are_disclosed() {
        let gel = simulate(&[fragments(&[20_000, 1_000, 100])], &settings(1.0)).unwrap();
        assert!(gel.warnings[0].contains("2 fragments lie outside the 400–8,000 bp"));
        assert!(gel.warnings[0].contains("1X TAE"));
        let lane = Lane {
            label: "Uncut".into(),
            fragments: vec![LaneFragment {
                length: 5_000,
                circular: true,
                reference: false,
            }],
            equimolar: true,
            ..Default::default()
        };
        let gel = simulate(&[lane], &settings(1.0)).unwrap();
        assert!(gel.lanes[0].bands[0].circular);
        assert!(gel.warnings.iter().any(|w| w.contains("topology")));
        assert!(gel.warnings.iter().all(|w| !w.contains("outside")));
    }

    #[test]
    fn invalid_gels_are_refused() {
        let gel = settings(1.0);
        assert!(simulate(&[], &gel).is_err());
        assert!(simulate(&vec![fragments(&[1_000]); MAX_LANES + 1], &gel).is_err());
        assert!(simulate(&[fragments(&[])], &gel).is_err());
        assert!(simulate(&[fragments(&[0])], &gel).is_err());
        assert!(simulate(&[fragments(&[1_000])], &settings(3.0)).is_err());
    }

    #[test]
    fn ladder_catalogue_is_ordered_unique_and_marks_reference_bands() {
        for (i, ladder) in LADDERS.iter().enumerate() {
            assert!(
                ladder.sizes.windows(2).all(|w| w[0] > w[1]),
                "{}",
                ladder.name
            );
            assert!(ladder.reference.iter().all(|r| ladder.sizes.contains(r)));
            assert!(ladder.source.starts_with("https://"));
            assert!(LADDERS[..i].iter().all(|l| l.name != ladder.name));
            let lane = Lane::ladder(ladder);
            assert_eq!(
                lane.fragments.iter().filter(|f| f.reference).count(),
                ladder.reference.len()
            );
        }
        let gel = simulate(&[Lane::ladder(&LADDERS[0])], &settings(1.0)).unwrap();
        assert_eq!(gel.lanes[0].bands.len(), 14);
        assert_eq!(gel.lanes[0].bands.iter().filter(|b| b.reference).count(), 3);
        assert!(ladder("λ DNA/HindIII").is_some_and(|l| l.equimolar));
        assert!(ladder("Unknown").is_none());
    }

    #[test]
    fn digest_lanes_use_product_lengths_and_topology() {
        let ecori = ENZYMES.iter().position(|e| e.name == "EcoRI").unwrap();
        let cut = digest(&plasmid(), true, &[ecori], DigestConditions::default()).unwrap();
        let lane = Lane::digest("EcoRI", &cut.fragments);
        let mut lengths = lane.fragments.iter().map(|f| f.length).collect::<Vec<_>>();
        lengths.sort();
        assert_eq!(lengths, [1_000, 2_000]);
        assert!(lane.equimolar && lane.fragments.iter().all(|f| !f.circular));
    }

    #[test]
    fn sample_lanes_digest_report_noncutters_or_load_uncut_dna() {
        let lane = sample(&plasmid(), true, &["EcoRI", "BamHI"])
            .lane()
            .unwrap();
        assert_eq!(lane.fragments.len(), 2);
        assert_eq!(lane.noncutters, ["BamHI"]);
        assert!(!lane.provisional && !lane.ladder);
        let uncut = sample(&plasmid(), true, &[]).lane().unwrap();
        assert_eq!(uncut.fragments.len(), 1);
        assert!(uncut.fragments[0].circular && uncut.fragments[0].length == 3_000);
        let linear = sample(&plasmid(), false, &[]).lane().unwrap();
        assert!(!linear.fragments[0].circular);
        assert!(sample(&plasmid(), true, &["NotAnEnzyme"]).lane().is_err());
        assert!(sample("ACGU", false, &[]).lane().is_err());
        assert!(
            LaneSpec::Ladder {
                name: "Unknown".into()
            }
            .lane()
            .is_err()
        );
        assert!(
            LaneSpec::CustomLadder {
                label: "Custom".into(),
                sizes: vec![]
            }
            .validate()
            .is_err()
        );
    }

    #[test]
    fn gel_files_round_trip_and_reject_foreign_or_invalid_content() {
        let spec = GelSpec::new(
            GelSettings {
                agarose_percent: 1.25,
                buffer: Buffer::Tbe,
                run_length: 1.2,
                band_edge: BandEdge::Top,
            },
            vec![
                LaneSpec::Ladder {
                    name: LADDERS[1].name.into(),
                },
                sample(&plasmid(), true, &["EcoRI"]),
                LaneSpec::CustomLadder {
                    label: "Custom".into(),
                    sizes: vec![5_000, 500],
                },
            ],
        );
        let text = spec.to_json().unwrap();
        assert!(text.contains("\"format\": \"popgene-gel\""));
        assert_eq!(GelSpec::from_json(&text).unwrap(), spec);
        let foreign = text.replace("popgene-gel", "other");
        assert!(GelSpec::from_json(&foreign).is_err());
        let future = text.replace("\"version\": 1", "\"version\": 2");
        assert!(
            GelSpec::from_json(&future)
                .unwrap_err()
                .contains("version 2")
        );
        let unknown = text.replace(LADDERS[1].name, "Unknown ladder");
        assert!(GelSpec::from_json(&unknown).is_err());
        assert!(GelSpec::from_json("{").is_err());
        assert!(
            GelSpec::new(GelSettings::default(), vec![])
                .to_json()
                .is_err()
        );
        let minimal = r#"{"format":"popgene-gel","version":1,"settings":{"agarose_percent":1.0},"lanes":[{"kind":"ladder","name":"GeneRuler 1 kb"}]}"#;
        assert_eq!(
            GelSpec::from_json(minimal).unwrap().settings,
            GelSettings::default()
        );
    }
}
