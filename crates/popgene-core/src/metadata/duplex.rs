//! Annotations that need the full physical duplex extent, including residues
//! that occur only on the opposite strand. They never invent top-strand bases.
use crate::{Error, Feature, Segment, metadata::SequenceColor};
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct DuplexAnnotations {
    pub features: Vec<Feature>,
    pub sequence_colors: Vec<SequenceColor>,
}
impl DuplexAnnotations {
    pub fn is_empty(&self) -> bool {
        self.features.is_empty() && self.sequence_colors.is_empty()
    }
    pub(crate) fn replaced(
        &self,
        start: usize,
        end: usize,
        replacement: usize,
    ) -> Result<Self, Error> {
        Ok(Self {
            features: replaced_features(&self.features, start, end, replacement)?,
            sequence_colors: self
                .sequence_colors
                .iter()
                .filter_map(|color| color.replaced(start, end, replacement))
                .collect(),
        })
    }
}

/// Both strand spans are measured left to right on one forward reference axis.
/// `opposite_sequence` runs from bottom.end toward bottom.start (5′→3′).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DuplexGeometry {
    pub extent: usize,
    pub top: Segment,
    pub bottom: Segment,
}
impl crate::workflow::DuplexEnds {
    /// Geometry for an already validated pair of physical strands.
    pub fn geometry(&self, top_length: usize) -> DuplexGeometry {
        let (top_prefix, top_suffix, bottom_prefix, bottom_suffix) = self.overhang_lengths();
        let extent = bottom_suffix + top_length + bottom_prefix;
        DuplexGeometry {
            extent,
            top: Segment {
                start: bottom_suffix,
                end: bottom_suffix + top_length,
            },
            bottom: Segment {
                start: top_prefix,
                end: extent.saturating_sub(top_suffix),
            },
        }
    }
    /// A coordinate reference for serialization/extraction. Outside the actual
    /// top span these are complements of existing bottom residues, not physical
    /// top-strand residues; use geometry() when displaying the strands.
    pub fn reference_sequence(&self, top: &str) -> Result<String, Error> {
        self.validate(top)?;
        let (_, _, bottom_prefix, bottom_suffix) = self.overhang_lengths();
        Ok(format!(
            "{}{}{}",
            crate::reverse_complement(
                &self.opposite_sequence[self.opposite_sequence.len() - bottom_suffix..]
            ),
            top,
            crate::reverse_complement(&self.opposite_sequence[..bottom_prefix])
        ))
    }
}

pub(crate) fn replaced_features(
    input: &[Feature],
    start: usize,
    end: usize,
    replacement: usize,
) -> Result<Vec<Feature>, Error> {
    let delta = replacement as isize - (end - start) as isize;
    let shift = |value: usize| (value as isize + delta) as usize;
    let mut features = input.to_vec();
    for feature in &mut features {
        crate::translation::remap_translation_exceptions(feature, start, end, replacement);
        feature.shift_display_segments(start, end, replacement);
        let mut changed = false;
        if let Some(location) = &feature.location {
            changed = feature.parts.iter().any(|part| {
                part.start < end && part.end > start
                    || (start == end && part.start < start && part.end > start)
            });
            feature.location = location
                .replaced(start, end, replacement)
                .map_err(|e| Error::Operation(e.to_string()))?;
            if let Some(location) = &feature.location {
                (feature.parts, feature.reverse) = location
                    .uniform_projection()
                    .unwrap_or_else(|| (location.local_parts(), false));
            } else {
                feature.parts.clear();
            }
            if changed {
                feature.edited = true;
                feature.qualifiers.retain(|(key, _)| key != "translation");
            }
            continue;
        }
        for part in &mut feature.parts {
            if part.end <= start {
                continue;
            }
            if part.start >= end {
                part.start = shift(part.start);
                part.end = shift(part.end);
                continue;
            }
            changed = true;
            part.start = part.start.min(start);
            part.end = if part.end > end {
                shift(part.end)
            } else {
                start + replacement
            };
        }
        feature.parts.retain(|part| part.start < part.end);
        if changed {
            feature.edited = true;
            feature.qualifiers.retain(|(key, _)| key != "translation");
        }
    }
    features.retain(|f| {
        !f.parts.is_empty()
            || f.location.is_some()
            || (f.kind == "primer_bind" && f.primer_sequence().is_some())
    });
    Ok(features)
}
