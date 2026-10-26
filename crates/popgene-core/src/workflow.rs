//! Sequence products retain physical ends and annotation coordinates.
use crate::restriction::{DigestConditions, DigestFragment, EndKind, FragmentEnd};
use crate::{Document, Error, Feature, FeatureLocation, Molecule, Segment, reverse_complement};
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct DuplexEnds {
    pub left: FragmentEnd,
    pub right: FragmentEnd,
    /// Actual opposite strand, 5′→3′, including its own protruding ends.
    pub opposite_sequence: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub left_modification: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub right_modification: Option<String>,
}
impl DuplexEnds {
    pub fn validate(&self, top: &str) -> Result<(), Error> {
        let invalid = || Error::Operation("Duplex ends do not match the strand sequences".into());
        for end in [&self.left, &self.right] {
            if !end
                .overhang
                .bytes()
                .all(|b| b"ACGTRYSWKMBDHVN".contains(&b))
                || matches!(end.kind, EndKind::Blunt | EndKind::Native) && !end.overhang.is_empty()
                || matches!(end.kind, EndKind::FivePrime | EndKind::ThreePrime)
                    && end.overhang.is_empty()
            {
                return Err(invalid());
            }
        }
        if !self
            .opposite_sequence
            .bytes()
            .all(|b| b"ACGTRYSWKMBDHVN".contains(&b))
        {
            return Err(invalid());
        }
        let (tp, ts, bp, bs) = self.overhang_lengths();
        let bottom = &self.opposite_sequence;
        if tp + ts >= top.len() || bp + bs >= bottom.len() {
            return Err(invalid());
        }
        if reverse_complement(&top[tp..top.len() - ts]) != bottom[bp..bottom.len() - bs] {
            return Err(invalid());
        }
        let left_actual = match self.left.kind {
            EndKind::FivePrime => &top[..tp],
            EndKind::ThreePrime => &bottom[bottom.len() - bs..],
            _ => "",
        };
        let right_actual = match self.right.kind {
            EndKind::FivePrime => &bottom[..bp],
            EndKind::ThreePrime => &top[top.len() - ts..],
            _ => "",
        };
        if left_actual != self.left.overhang || right_actual != self.right.overhang {
            return Err(invalid());
        }
        Ok(())
    }
    /// Top prefix/suffix and opposite strand prefix/suffix, all 5′→3′.
    pub fn overhang_lengths(&self) -> (usize, usize, usize, usize) {
        let len = |end: &FragmentEnd, kind| {
            if end.kind == kind {
                end.overhang.len()
            } else {
                0
            }
        };
        (
            len(&self.left, EndKind::FivePrime),
            len(&self.right, EndKind::ThreePrime),
            len(&self.right, EndKind::FivePrime),
            len(&self.left, EndKind::ThreePrime),
        )
    }
    pub(crate) fn replaced(
        &self,
        top: &str,
        start: usize,
        end: usize,
        replacement: &str,
    ) -> Result<Self, Error> {
        self.validate(top)?;
        let (tp, ts, bp, _) = self.overhang_lengths();
        if start < tp || end > top.len() - ts {
            return Err(Error::Operation("This edit overlaps a cohesive end. Prepare or remove the end before editing its bases.".into()));
        }
        let mut result = self.clone();
        let paired = top.len() - tp - ts;
        let bottom_start = bp + paired - (end - tp);
        let bottom_end = bp + paired - (start - tp);
        result
            .opposite_sequence
            .replace_range(bottom_start..bottom_end, &reverse_complement(replacement));
        Ok(result)
    }
}

/// Translate local annotation coordinates without changing remote accession
/// coordinates, biological order, or the interpretation of location operators.
pub fn shift_feature(feature: &mut Feature, delta: isize) -> Result<(), Error> {
    fn position(value: usize, delta: isize) -> Result<usize, Error> {
        value.checked_add_signed(delta).ok_or(Error::Coordinates)
    }
    fn location(node: &mut FeatureLocation, delta: isize) -> Result<(), Error> {
        use FeatureLocation::*;
        match node {
            Range { start, end, .. } => {
                *start = position(*start, delta)?;
                *end = position(*end, delta)?;
            }
            Between { left, right } => {
                *left = position(*left, delta)?;
                *right = position(*right, delta)?;
            }
            Complement(child) => location(child, delta)?,
            Join(children) | Order(children) | Bond(children) | OneOf(children) => {
                for child in children {
                    location(child, delta)?;
                }
            }
            Remote { .. } | Gap { .. } => {}
        }
        Ok(())
    }
    let mut result = feature.clone();
    let mut parts = result.parts.clone();
    for part in &mut parts {
        part.start = position(part.start, delta)?;
        part.end = position(part.end, delta)?;
    }
    let mut mapped_location = result.location.clone();
    if let Some(node) = &mut mapped_location {
        location(node, delta)?;
    }
    // Check display-only coordinates too before calling the infallible mapper.
    if result
        .display_segments()
        .iter()
        .chain(result.display_gaps().iter())
        .any(|p| {
            p.range.start.checked_add_signed(delta).is_none()
                || p.range.end.checked_add_signed(delta).is_none()
        })
        || result
            .display_cleavage_sites()
            .iter()
            .any(|p| p.checked_add_signed(delta).is_none())
    {
        return Err(Error::Coordinates);
    }
    // The display metadata is still on the original axis at this point.
    result.project_display_segments(|p| {
        vec![Segment {
            start: p
                .start
                .checked_add_signed(delta)
                .expect("checked display start"),
            end: p
                .end
                .checked_add_signed(delta)
                .expect("checked display end"),
        }]
    });
    result.parts = parts;
    result.location = mapped_location;
    crate::translation::project_translation_exceptions(&mut result, |p| {
        p.checked_add_signed(delta)
    });
    *feature = result;
    Ok(())
}

/// Temporary coordinate reference covering the union of both physical strands.
/// Outside the top span its bases are complements of bottom residues. This is
/// for annotation projection, never an additional physical strand or UI record.
pub fn duplex_reference(source: &Document) -> Result<Document, Error> {
    source.validate()?;
    let Some(ends) = &source.duplex_ends else {
        let mut reference = source.clone();
        crate::lineage::to_reference(source, &mut reference)?;
        return Ok(reference);
    };
    let geometry = ends.geometry(source.sequence.len());
    let mut reference = source.clone();
    for feature in &mut reference.features {
        shift_feature(feature, geometry.top.start as isize)?;
    }
    for color in &mut reference.sequence_colors {
        color.range.start += geometry.top.start;
        color.range.end += geometry.top.start;
    }
    reference
        .features
        .append(&mut reference.duplex_annotations.features);
    reference
        .sequence_colors
        .append(&mut reference.duplex_annotations.sequence_colors);
    reference.sequence = ends.reference_sequence(&source.sequence)?.into();
    reference.sequence_qualities = None;
    reference.opposite_sequence_qualities = None;
    reference.duplex_ends = None;
    crate::lineage::to_reference(source, &mut reference)?;
    reference.validate()?;
    Ok(reference)
}

/// Partition annotations already on the full extent axis onto the actual top
/// strand or the full-extent collection. Cross-boundary annotations stay whole.
pub fn partition_duplex_annotations(document: &mut Document) -> Result<(), Error> {
    let Some(ends) = &document.duplex_ends else {
        return Ok(());
    };
    let top = ends.geometry(document.sequence.len()).top;
    if !document.duplex_annotations.is_empty() {
        return Err(Error::Operation(
            "Partition requires annotations on one extent axis".into(),
        ));
    }
    let mut features = Vec::new();
    for mut feature in std::mem::take(&mut document.features) {
        let outside = feature
            .parts
            .iter()
            .any(|p| p.start < top.start || p.end > top.end)
            || feature
                .display_gaps()
                .iter()
                .any(|p| p.range.start < top.start || p.range.end > top.end)
            || feature
                .display_cleavage_sites()
                .iter()
                .any(|&p| p < top.start || p > top.end)
            || feature.location.as_ref().is_some_and(|l| {
                l.local_boundaries()
                    .iter()
                    .any(|&p| p <= top.start || p >= top.end)
            });
        if outside {
            document.duplex_annotations.features.push(feature);
        } else {
            shift_feature(&mut feature, -(top.start as isize))?;
            features.push(feature);
        }
    }
    document.features = features;
    let mut colors = Vec::new();
    for mut color in std::mem::take(&mut document.sequence_colors) {
        if color.range.start < top.start || color.range.end > top.end {
            document.duplex_annotations.sequence_colors.push(color);
        } else {
            color.range.start -= top.start;
            color.range.end -= top.start;
            colors.push(color);
        }
    }
    document.sequence_colors = colors;
    crate::lineage::finish(document, true)?;
    document.validate()
}

/// Rotate the complete physical molecule: the original opposite strand becomes
/// the new top strand, with both ends and all full-extent annotations retained.
pub fn reverse_duplex(source: &Document) -> Result<Document, Error> {
    let ends = source
        .duplex_ends
        .as_ref()
        .ok_or_else(|| Error::Operation("No physical duplex ends".into()))?;
    let mut result = duplex_reference(source)?;
    result.revcomp()?;
    result.sequence = ends.opposite_sequence.clone().into();
    result.sequence_qualities = source.opposite_sequence_qualities.clone();
    result.opposite_sequence_qualities = source.sequence_qualities.clone();
    result.duplex_ends = Some(DuplexEnds {
        left: ends.right.clone(),
        right: ends.left.clone(),
        opposite_sequence: source.sequence.to_string(),
        left_modification: ends.right_modification.clone(),
        right_modification: ends.left_modification.clone(),
    });
    partition_duplex_annotations(&mut result)?;
    Ok(result)
}

struct Projection {
    spans: Vec<Segment>,
}
impl Projection {
    fn mark_partial_ends(location: &mut FeatureLocation, first: bool, last: bool) {
        use FeatureLocation::*;
        match location {
            Range { before, after, .. } => {
                *before |= first;
                *after |= last;
            }
            Complement(child) => Self::mark_partial_ends(child, last, first),
            Join(children) => {
                if first && let Some(child) = children.first_mut() {
                    Self::mark_partial_ends(child, true, false);
                }
                if last && let Some(child) = children.last_mut() {
                    Self::mark_partial_ends(child, false, true);
                }
            }
            // Nonlocal/ambiguous operators already have unresolved endpoints;
            // never alter remote accession coordinates to express local clipping.
            _ => {}
        }
    }
    fn position(&self, position: usize) -> Option<usize> {
        let mut offset = 0;
        for span in &self.spans {
            if (span.start..span.end).contains(&position) {
                return Some(offset + position - span.start);
            }
            offset += span.end - span.start;
        }
        None
    }
    fn range(&self, range: &Segment) -> Vec<Segment> {
        let mut mapped = vec![];
        let mut offset = 0;
        for span in &self.spans {
            let start = range.start.max(span.start);
            let end = range.end.min(span.end);
            if start < end {
                mapped.push((
                    start,
                    Segment {
                        start: offset + start - span.start,
                        end: offset + end - span.start,
                    },
                ));
            }
            offset += span.end - span.start;
        }
        mapped.sort_by_key(|(start, _)| *start);
        mapped.into_iter().map(|(_, s)| s).collect()
    }
    fn location(&self, location: &FeatureLocation) -> Option<FeatureLocation> {
        use FeatureLocation::*;
        match location {
            Range {
                start,
                end,
                before,
                after,
            } => {
                let ranges = self.range(&Segment {
                    start: *start,
                    end: *end,
                });
                let count = ranges.len();
                let nodes = ranges
                    .into_iter()
                    .enumerate()
                    .map(|(i, s)| Range {
                        start: s.start,
                        end: s.end,
                        before: i == 0 && (*before || self.position(*start).is_none()),
                        after: i + 1 == count && (*after || self.position(end - 1).is_none()),
                    })
                    .collect::<Vec<_>>();
                match nodes.len() {
                    0 => None,
                    1 => nodes.into_iter().next(),
                    _ => Some(Join(nodes)),
                }
            }
            Between { left, right } => {
                let (left, right) = (self.position(*left)?, self.position(*right)?);
                (left.checked_add(1) == Some(right)).then_some(Between { left, right })
            }
            Complement(c) => self.location(c).map(|n| Complement(Box::new(n))),
            Join(children) | Order(children) | Bond(children) | OneOf(children) => {
                let mut mapped = children
                    .iter()
                    .enumerate()
                    .filter_map(|(index, c)| self.location(c).map(|location| (index, location)))
                    .collect::<Vec<_>>();
                if mapped.is_empty() {
                    return None;
                }
                if matches!(location, Join(_)) {
                    // A whole excluded terminal part leaves no Range on which
                    // to record clipping. Carry that partial boundary onto the
                    // retained end, including through complement/nested joins.
                    let (index, first) = mapped.first_mut().unwrap();
                    Self::mark_partial_ends(first, *index > 0, false);
                    let (index, last) = mapped.last_mut().unwrap();
                    Self::mark_partial_ends(last, false, *index + 1 < children.len());
                }
                let children = mapped.into_iter().map(|(_, location)| location).collect();
                Some(match location {
                    Join(_) => Join(children),
                    Order(_) => Order(children),
                    Bond(_) => Bond(children),
                    _ => OneOf(children),
                })
            }
            Remote { .. } | Gap { .. } => Some(location.clone()),
        }
    }
}

/// Extract a contiguous region, optionally across a circular origin. Features
/// remain in biological order; clipped boundaries are fuzzy and marked for review.
pub fn extract_region(
    source: &Document,
    start: usize,
    length: usize,
    name: impl Into<String>,
) -> Result<Document, Error> {
    source.validate()?;
    let n = source.sequence.len();
    if length == 0 || length > n || start >= n || (!source.circular && start + length > n) {
        return Err(Error::Coordinates);
    }
    let end = start + length;
    let spans = if end <= n {
        vec![Segment { start, end }]
    } else {
        vec![
            Segment { start, end: n },
            Segment {
                start: 0,
                end: end - n,
            },
        ]
    };
    let projection = Projection { spans };
    let sequence = projection
        .spans
        .iter()
        .map(|s| &source.sequence[s.start..s.end])
        .collect::<String>();
    let mut result = Document::new_typed(name, &sequence, source.molecule)?;
    result.single_stranded = source.single_stranded;
    result.description = source.description.clone();
    result.references = source.references.clone();
    result.sequence_colors = source
        .sequence_colors
        .iter()
        .flat_map(|color| {
            projection
                .range(&color.range)
                .into_iter()
                .map(|range| crate::metadata::SequenceColor {
                    range,
                    ..color.clone()
                })
        })
        .collect();
    if let Some(qualities) = &source.sequence_qualities {
        result.sequence_qualities = Some(
            projection
                .spans
                .iter()
                .flat_map(|s| qualities[s.start..s.end].iter().copied())
                .collect(),
        );
    }
    if let Some(scores) = &source.opposite_sequence_qualities {
        let geometry = source.duplex_ends.as_ref().map(|ends| ends.geometry(n));
        let top_origin = geometry.as_ref().map_or(0, |g| g.top.start);
        let bottom = geometry
            .as_ref()
            .map_or(Segment { start: 0, end: n }, |g| g.bottom.clone());
        let positions = projection
            .spans
            .iter()
            .flat_map(|s| s.start..s.end)
            .collect::<Vec<_>>();
        result.opposite_sequence_qualities = Some(
            positions
                .into_iter()
                .rev()
                .map(|p| {
                    let axis = p + top_origin;
                    if (bottom.start..bottom.end).contains(&axis) {
                        scores[bottom.end - 1 - axis]
                    } else {
                        None
                    }
                })
                .collect(),
        );
    }
    for feature in &source.features {
        // Unbound oligos belong to the source inventory, not a physical fragment.
        if feature.parts.is_empty() && feature.location.is_none() {
            continue;
        }
        let original = feature.location.clone().unwrap_or_else(|| {
            let mut parts = feature
                .parts
                .iter()
                .map(|s| FeatureLocation::Range {
                    start: s.start,
                    end: s.end,
                    before: false,
                    after: false,
                })
                .collect::<Vec<_>>();
            let location = if parts.len() == 1 {
                parts.remove(0)
            } else {
                FeatureLocation::Join(parts)
            };
            if feature.reverse {
                FeatureLocation::Complement(Box::new(location))
            } else {
                location
            }
        });
        let Some(location) = projection.location(&original) else {
            continue;
        };
        let mut mapped = feature.clone();
        mapped.project_display_segments(|range| projection.range(range));
        let (parts, reverse) = location
            .uniform_projection()
            .unwrap_or_else(|| (location.local_parts(), false));
        let retained: usize = parts.iter().map(|s| s.end - s.start).sum();
        if retained < feature.len() {
            mapped.edited = true;
            mapped.qualifiers.retain(|(key, _)| key != "translation");
            if !mapped.qualifiers.iter().any(|(key, _)| key == "partial") {
                mapped.qualifiers.push(("partial".into(), None));
            }
            if feature.kind.eq_ignore_ascii_case("CDS")
                && original.uniform_projection().is_some()
                && let Ok(offset) = crate::translation::cds_phase(feature)
            {
                let positions = feature
                    .parts
                    .iter()
                    .flat_map(|s| s.start..s.end)
                    .collect::<Vec<_>>();
                let removed = if feature.reverse {
                    positions
                        .iter()
                        .rev()
                        .take_while(|&&p| projection.position(p).is_none())
                        .count()
                } else {
                    positions
                        .iter()
                        .take_while(|&&p| projection.position(p).is_none())
                        .count()
                };
                let phase = (offset + 3 - removed % 3) % 3 + 1;
                mapped.qualifiers.retain(|(key, _)| key != "codon_start");
                mapped
                    .qualifiers
                    .push(("codon_start".into(), Some(phase.to_string())));
                // Native phase is feature-relative too. Update both recorded
                // forms only after validating the source values agree; malformed
                // or conflicting values remain intact and untranslatable.
                for (_, value) in mapped
                    .qualifiers
                    .iter_mut()
                    .filter(|(key, _)| key == "snapgene_readingFrame")
                {
                    *value = Some(format!("{}{phase}", if feature.reverse { "-" } else { "" }));
                }
            }
        }
        mapped.parts = parts;
        mapped.reverse = reverse;
        mapped.location = Some(location);
        crate::translation::project_translation_exceptions(&mut mapped, |p| projection.position(p));
        result.features.push(mapped);
    }
    let mut offset = 0;
    let source_top = source
        .duplex_ends
        .as_ref()
        .map_or(0, |e| e.geometry(source.sequence.len()).top.start);
    let maps = projection
        .spans
        .iter()
        .map(|span| {
            let map = crate::lineage::Map {
                source: Segment {
                    start: span.start + source_top,
                    end: span.end + source_top,
                },
                target_start: offset,
                reverse: false,
            };
            offset += span.end - span.start;
            map
        })
        .collect::<Vec<_>>();
    crate::lineage::project(source, &mut result, &maps)?;
    result.comments.push(format!(
        "Extracted from {}: {}..{}{} (original coordinates).",
        source.name,
        start + 1,
        if end <= n { end } else { end - n },
        if end > n {
            " across circular origin"
        } else {
            ""
        }
    ));
    if result.features.iter().any(|f| f.edited) {
        result.warnings.push("Annotations cut by product boundaries are partial; review their functions and translations.".into());
    }
    result.validate()?;
    Ok(result)
}

/// Project up to one circumference plus the short duplicated cohesive span.
/// The extra reference bases represent physical residues on opposite ends of
/// different strands, not a synthesized extension of either strand.
fn digest_reference_region(
    source: &Document,
    start: usize,
    length: usize,
    name: String,
) -> Result<Document, Error> {
    let n = source.sequence.len();
    if length <= n {
        return extract_region(source, start, length, name);
    }
    if !source.circular || length > crate::MAX_BASES || length >= n * 2 {
        return Err(Error::Coordinates);
    }
    let mut result = extract_region(source, start, n, name)?;
    let prefix = extract_region(&result, 0, length - n, "cohesive reference overlap")?;
    result.replace(n, n, &prefix.sequence)?;
    let mut lineage_bindings = Vec::new();
    let mut next_id = result.features.iter().map(|f| f.id).max().unwrap_or(0);
    for mut feature in prefix.features.clone() {
        shift_feature(&mut feature, n as isize)?;
        next_id = next_id.checked_add(1).ok_or(Error::Coordinates)?;
        lineage_bindings.push((feature.id, next_id));
        feature.id = next_id;
        result.features.push(feature);
    }
    result
        .sequence_colors
        .extend(prefix.sequence_colors.iter().cloned().map(|mut color| {
            color.range.start += n;
            color.range.end += n;
            color
        }));
    crate::lineage::append(&prefix, &mut result, n, &lineage_bindings)?;
    crate::lineage::finish(&mut result, false)?;
    result.validate()?;
    Ok(result)
}

pub fn digest_products(
    source: &Document,
    fragments: &[DigestFragment],
    enzymes: &[String],
    conditions: DigestConditions,
) -> Result<Vec<Document>, Error> {
    source.validate()?;
    if source.molecule != Molecule::Dna || source.single_stranded || source.trace.is_some() {
        return Err(Error::Operation(
            "Digestion products require double-stranded DNA".into(),
        ));
    }
    let reference = duplex_reference(source)?;
    let source_origin = source
        .duplex_ends
        .as_ref()
        .map_or(0, |ends| ends.geometry(source.sequence.len()).top.start);
    fragments.iter().enumerate().map(|(i, fragment)| {
        let name = format!("{} digest {}", source.name, i + 1);
        let mut product = if fragment.circular {
            let mut result = extract_region(&reference, fragment.start, fragment.sequence.len(), name)?;
            result.circular = true;
            result
        } else {
            let ends = DuplexEnds {
                left: fragment.left_end.clone().ok_or_else(|| Error::Operation("Missing left product end".into()))?,
                right: fragment.right_end.clone().ok_or_else(|| Error::Operation("Missing right product end".into()))?,
                opposite_sequence: fragment.opposite_sequence.clone(),
                left_modification: fragment.left_modification.clone(), right_modification: fragment.right_modification.clone(),
            };
            ends.validate(&fragment.sequence)?;
            let geometry = ends.geometry(fragment.sequence.len());
            let axis_start = source_origin as isize + fragment.start as isize - geometry.top.start as isize;
            let axis_start = if source.circular { axis_start.rem_euclid(reference.sequence.len() as isize) as usize } else { usize::try_from(axis_start).map_err(|_| Error::Coordinates)? };
            let mut result = digest_reference_region(&reference, axis_start, geometry.extent, name)?;
            if result.sequence.get(geometry.top.start..geometry.top.end) != Some(fragment.sequence.as_str()) {
                return Err(Error::Operation("Digest physical extent does not match its source sequence".into()));
            }
            if let Some(scores) = result.sequence_qualities.as_mut() {
                *scores = scores[geometry.top.start..geometry.top.end].to_vec();
            }
            // A color belongs to a real residue on its declared strand. Avoid
            // duplicating a top color onto a bottom-only terminal reference.
            result.sequence_colors = result.sequence_colors.into_iter().filter_map(|mut color| {
                let strand = if color.reverse { &geometry.bottom } else { &geometry.top };
                color.range.start = color.range.start.max(strand.start);
                color.range.end = color.range.end.min(strand.end);
                (color.range.start < color.range.end).then_some(color)
            }).collect();
            result.sequence = fragment.sequence.clone().into();
            if let Some(scores) = &source.opposite_sequence_qualities {
                let bottom = source.duplex_ends.as_ref().map_or(Segment{start:0,end:source.sequence.len()},|ends|ends.geometry(source.sequence.len()).bottom);
                result.opposite_sequence_qualities = Some((0..fragment.opposite_sequence.len()).map(|i| {
                    let axis = axis_start+geometry.bottom.end-1-i;
                    let axis = if source.circular {axis%reference.sequence.len()}else{axis};
                    if !(bottom.start..bottom.end).contains(&axis) { return Err(Error::Coordinates); }
                    scores.get(bottom.end-1-axis).copied().ok_or(Error::Coordinates)
                }).collect::<Result<Vec<_>,_>>()?);
            }
            result.duplex_ends = Some(ends);
            partition_duplex_annotations(&mut result)?;
            result
        };
        if product.sequence.as_ref() != fragment.sequence { return Err(Error::Operation("Digest result does not match its source sequence".into())); }
        // Cleavage changes physical ends, not the template's recorded chemistry.
        product.methylation = source.methylation;
        product.sequence_qualities = source.sequence_qualities.as_ref().map(|scores| {
            (0..fragment.sequence.len()).map(|i|scores[(fragment.start+i)%scores.len()]).collect()
        });
        product.comments.extend(source.comments.iter().cloned());
        for warning in &source.warnings { if !product.warnings.contains(warning) { product.warnings.push(warning.clone()); } }
        product.comments.push(format!("Restriction digestion: {}. Temperature {} °C; pH {}; magnesium {} mM; glycerol {}%; Dam {}; Dcm {}; CpG {}. Assumes complete cleavage at permitted sites.", enzymes.join(", "), conditions.temperature_c, conditions.ph, conditions.magnesium_mm, conditions.glycerol_percent, conditions.methylation.dam, conditions.methylation.dcm, conditions.methylation.cpg));
        product.validate()?;
        Ok(product)
    }).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn unknown_quality_positions_remain_distinct_from_phred_zero_through_edits() {
        let mut legacy =
            serde_json::to_value(Document::new("legacy scores", "ACGT").unwrap()).unwrap();
        legacy["sequence_qualities"] = serde_json::json!([0, 10, 20, 30]);
        let legacy: Document = serde_json::from_value(legacy).unwrap();
        assert_eq!(
            legacy.sequence_qualities,
            Some(vec![Some(0), Some(10), Some(20), Some(30)])
        );
        let mut source = Document::new("measured", "ACGTAC").unwrap();
        source.sequence_qualities =
            Some(vec![Some(0), None, Some(12), Some(13), Some(14), Some(15)]);
        source.opposite_sequence_qualities =
            Some(vec![Some(30), Some(31), Some(32), None, Some(34), Some(35)]);
        source.replace(2, 2, "TT").unwrap();
        let serialized = serde_json::to_value(&source).unwrap();
        assert_eq!(serialized["sequence_qualities"][0], serde_json::json!(0));
        assert!(serialized["sequence_qualities"][1].is_null());
        assert_eq!(
            source.sequence_qualities,
            Some(vec![
                Some(0),
                None,
                None,
                None,
                Some(12),
                Some(13),
                Some(14),
                Some(15)
            ])
        );
        assert_eq!(
            source.opposite_sequence_qualities,
            Some(vec![
                Some(30),
                Some(31),
                Some(32),
                None,
                None,
                None,
                Some(34),
                Some(35)
            ])
        );
        let original = source.clone();
        source.revcomp().unwrap();
        assert_eq!(
            source.sequence_qualities,
            original.opposite_sequence_qualities
        );
        assert_eq!(
            source.opposite_sequence_qualities,
            original.sequence_qualities
        );
        source.revcomp().unwrap();
        assert_eq!(source.sequence_qualities, original.sequence_qualities);
        let mut unknown = Document::new("one strand observed", "ACGT").unwrap();
        unknown.sequence_qualities = Some(vec![Some(10); 4]);
        unknown.revcomp().unwrap();
        assert!(unknown.sequence_qualities.is_none());
        assert_eq!(unknown.opposite_sequence_qualities, Some(vec![Some(10); 4]));
    }
    #[test]
    fn physical_reversal_preserves_both_strands_and_end_only_annotations() {
        use crate::metadata::SequenceColor;
        for left in [EndKind::FivePrime, EndKind::ThreePrime] {
            for right in [EndKind::FivePrime, EndKind::ThreePrime] {
                let axis = "TAACCCGGGTTT";
                let top_start = usize::from(left == EndKind::ThreePrime) * 2;
                let top_end = axis.len() - usize::from(right == EndKind::FivePrime) * 3;
                let bottom_start = usize::from(left == EndKind::FivePrime) * 2;
                let bottom_end = axis.len() - usize::from(right == EndKind::ThreePrime) * 3;
                let mut doc = Document::new("physical", &axis[top_start..top_end]).unwrap();
                doc.duplex_ends = Some(DuplexEnds {
                    left: FragmentEnd {
                        kind: left,
                        overhang: if left == EndKind::FivePrime {
                            axis[..2].into()
                        } else {
                            reverse_complement(&axis[..2])
                        },
                        enzymes: vec!["left enzyme".into()],
                    },
                    right: FragmentEnd {
                        kind: right,
                        overhang: if right == EndKind::ThreePrime {
                            axis[axis.len() - 3..].into()
                        } else {
                            reverse_complement(&axis[axis.len() - 3..])
                        },
                        enzymes: vec!["right enzyme".into()],
                    },
                    opposite_sequence: reverse_complement(&axis[bottom_start..bottom_end]),
                    left_modification: Some("Phosphorylated".into()),
                    right_modification: None,
                });
                // Start in full-extent coordinates, including both one-strand ends.
                for (id, (start, end)) in [(0, 1), (1, 5), (5, 7), (10, 12)].into_iter().enumerate()
                {
                    doc.features.push(Feature {
                        id: id as u64 + 1,
                        label: format!("f{id}"),
                        kind: "misc_feature".into(),
                        parts: vec![Segment { start, end }],
                        reverse: id % 2 == 0,
                        location: None,
                        qualifiers: vec![],
                        edited: false,
                    });
                }
                doc.sequence_colors = vec![SequenceColor {
                    range: Segment { start: 0, end: 4 },
                    reverse: true,
                    color_rgb: 0x123456,
                }];
                partition_duplex_annotations(&mut doc).unwrap();
                let before = duplex_reference(&doc).unwrap();
                let original_ends = doc.duplex_ends.clone().unwrap();
                doc.revcomp().unwrap();
                assert_eq!(doc.sequence.as_ref(), original_ends.opposite_sequence);
                assert_eq!(
                    doc.duplex_ends.as_ref().unwrap().opposite_sequence,
                    axis[top_start..top_end]
                );
                assert_eq!(doc.duplex_ends.as_ref().unwrap().left, original_ends.right);
                let mut expected = before.clone();
                expected.revcomp().unwrap();
                let actual = duplex_reference(&doc).unwrap();
                for feature in &expected.features {
                    let found = actual.features.iter().find(|f| f.id == feature.id).unwrap();
                    assert_eq!(found.parts, feature.parts);
                    assert_eq!(found.reverse, feature.reverse);
                }
                assert_eq!(actual.sequence_colors, expected.sequence_colors);
                doc.revcomp().unwrap();
                assert_eq!(doc.duplex_ends.as_ref().unwrap(), &original_ends);
                let restored = duplex_reference(&doc).unwrap();
                for feature in &before.features {
                    assert_eq!(
                        restored
                            .features
                            .iter()
                            .find(|f| f.id == feature.id)
                            .unwrap(),
                        feature
                    );
                }
                assert_eq!(restored.sequence_colors, before.sequence_colors);
            }
        }
    }
    #[test]
    fn circular_extraction_keeps_joined_feature_sequence_and_display() {
        let mut source = Document::new("ring", "AAAACCCCGGGGTTTT").unwrap();
        source.circular = true;
        source.sequence_colors.push(crate::metadata::SequenceColor {
            range: Segment { start: 0, end: 3 },
            reverse: true,
            color_rgb: 0x123456,
        });
        source
            .annotate("origin".into(), Segment { start: 14, end: 16 }, false)
            .unwrap();
        source.features[0].parts.push(Segment { start: 0, end: 4 });
        source.features[0].qualifiers.push(("note".into(), Some("This forward directional feature has 2 segments: 1: 15 .. 16 / #00ff00 / tail 2: 1 .. 4 / #ff0000 / head".into())));
        let product = extract_region(&source, 12, 8, "product").unwrap();
        assert_eq!(&*product.sequence, "TTTTAAAA");
        assert_eq!(
            product.sequence_colors[0].range,
            Segment { start: 4, end: 7 }
        );
        assert!(product.sequence_colors[0].reverse);
        assert_eq!(
            product.feature_sequence(&product.features[0]).unwrap(),
            source.feature_sequence(&source.features[0]).unwrap()
        );
        assert_eq!(
            product.features[0].parts,
            vec![Segment { start: 2, end: 4 }, Segment { start: 4, end: 8 }]
        );
        assert_eq!(
            product.features[0].display_segments()[0].name.as_deref(),
            Some("tail")
        );
        assert_eq!(
            product.features[0].display_segments()[1].color_rgb,
            Some(0xff0000)
        );
        assert!(!product.features[0].edited);
    }
    #[test]
    fn clipped_cds_keeps_frame_and_partial_boundaries() {
        let mut source = Document::new("cds", "ATGAAACCCTAA").unwrap();
        source
            .annotate("coding".into(), Segment { start: 0, end: 12 }, false)
            .unwrap();
        source.features[0].kind = "CDS".into();
        let product = extract_region(&source, 4, 8, "partial").unwrap();
        assert_eq!(product.features[0].qualifier("codon_start"), Some("3"));
        assert!(product.features[0].location_label().contains('<'));
        assert!(product.features[0].edited);
    }
    #[test]
    fn clipping_uses_validated_cds_phase_and_keeps_native_phase_consistent() {
        for reverse in [false, true] {
            for explicit_phase in [false, true] {
                let mut source = Document::new(
                    "phase 2",
                    if reverse {
                        "TTACTTGAATTCCAT"
                    } else {
                        "ATGGAATTCAAGTAA"
                    },
                )
                .unwrap();
                source
                    .annotate("coding".into(), Segment { start: 0, end: 15 }, reverse)
                    .unwrap();
                let feature = &mut source.features[0];
                feature.kind = "CDS".into();
                feature.qualifiers.push((
                    "snapgene_readingFrame".into(),
                    Some(if reverse { "-2" } else { "2" }.into()),
                ));
                if explicit_phase {
                    feature
                        .qualifiers
                        .push(("codon_start".into(), Some("2".into())));
                }
                // The first four coding-strand bases are removed. Source phase2
                // starts at its second base, so the remaining AAT begins a codon.
                let product =
                    extract_region(&source, if reverse { 0 } else { 4 }, 11, "partial").unwrap();
                let feature = &product.features[0];
                assert_eq!(feature.qualifier("codon_start"), Some("1"));
                assert_eq!(
                    feature.qualifier("snapgene_readingFrame"),
                    Some(if reverse { "-1" } else { "1" })
                );
                let translated = crate::translation::translate_cds(&product, feature).unwrap();
                assert_eq!(translated.amino_acids, b"NSS");
                assert_eq!(translated.trailing_bases, 2);
            }
        }
    }
    #[test]
    fn clipping_does_not_repair_invalid_or_conflicting_source_phase() {
        for qualifiers in [
            vec![("codon_start", Some("4"))],
            vec![("codon_start", None)],
            vec![("codon_start", Some("1")), ("codon_start", Some("3"))],
            vec![
                ("codon_start", Some("1")),
                ("snapgene_readingFrame", Some("2")),
            ],
        ] {
            let mut source = Document::new("invalid phase", "ATGGAATTCAAGTAA").unwrap();
            source
                .annotate("coding".into(), Segment { start: 0, end: 15 }, false)
                .unwrap();
            source.features[0].kind = "CDS".into();
            source.features[0].qualifiers = qualifiers
                .into_iter()
                .map(|(k, v)| (k.into(), v.map(str::to_string)))
                .collect();
            assert!(crate::translation::cds_phase(&source.features[0]).is_err());
            let product = extract_region(&source, 4, 11, "partial").unwrap();
            for qualifier in &source.features[0].qualifiers {
                assert!(product.features[0].qualifiers.contains(qualifier));
            }
            assert!(crate::translation::cds_phase(&product.features[0]).is_err());
            assert!(crate::translation::translate_cds(&product, &product.features[0]).is_err());
        }
    }
    #[test]
    fn rotation_keeps_fuzziness_only_at_the_original_range_endpoints() {
        let mut source = Document::new("fuzzy ring", "ACGTACGTAC").unwrap();
        source.circular = true;
        source
            .annotate("partial".into(), Segment { start: 1, end: 9 }, false)
            .unwrap();
        source.features[0].location = Some(FeatureLocation::Range {
            start: 1,
            end: 9,
            before: true,
            after: true,
        });
        let product = extract_region(&source, 5, 10, "rotated").unwrap();
        assert_eq!(
            product.features[0].location,
            Some(FeatureLocation::Join(vec![
                FeatureLocation::Range {
                    start: 6,
                    end: 10,
                    before: true,
                    after: false
                },
                FeatureLocation::Range {
                    start: 0,
                    end: 4,
                    before: false,
                    after: true
                },
            ]))
        );
        assert_eq!(
            product.feature_sequence(&product.features[0]).unwrap(),
            source.feature_sequence(&source.features[0]).unwrap()
        );
    }
    #[test]
    fn excluded_compound_end_parts_do_not_create_initiators_or_terminal_stops() {
        for reverse in [false, true] {
            for (sequence, parts, code, start, expected, partial) in [
                (
                    "ATGAAAAAAGTGAAA",
                    vec![Segment { start: 0, end: 3 }, Segment { start: 9, end: 15 }],
                    "11",
                    9,
                    b"VK".as_slice(),
                    (true, false),
                ),
                (
                    "ATGTGAAAAAAATGA",
                    vec![Segment { start: 0, end: 6 }, Segment { start: 12, end: 15 }],
                    "27",
                    0,
                    b"MW".as_slice(),
                    (false, true),
                ),
            ] {
                let mut source = Document::new("joined CDS", sequence).unwrap();
                source
                    .annotate("coding".into(), parts[0].clone(), false)
                    .unwrap();
                let feature = &mut source.features[0];
                feature.kind = "CDS".into();
                feature.parts = parts;
                feature.location = Some(FeatureLocation::Join(
                    feature
                        .parts
                        .iter()
                        .map(|p| FeatureLocation::Range {
                            start: p.start,
                            end: p.end,
                            before: false,
                            after: false,
                        })
                        .collect(),
                ));
                feature
                    .qualifiers
                    .push(("transl_table".into(), Some(code.into())));
                if reverse {
                    source.revcomp().unwrap();
                }
                let start = if reverse {
                    source.sequence.len() - start - 6
                } else {
                    start
                };
                let product = extract_region(&source, start, 6, "partial").unwrap();
                let feature = &product.features[0];
                assert_eq!(feature.location.as_ref().unwrap().partial_ends(), partial);
                assert_eq!(
                    crate::translation::translate_cds(&product, feature)
                        .unwrap()
                        .amino_acids,
                    expected
                );
            }
        }
    }
    #[test]
    fn digest_products_validate_both_strands_and_shared_base_edits() {
        let mut source = Document::new("template", "AAAAGAATTCCCCCGGTACCAAAA").unwrap();
        source.methylation = Some(crate::metadata::DnaMethylation {
            dam: true,
            dcm: false,
            eco_ki: true,
        });
        let enzymes = ["EcoRI", "KpnI"].map(|name| {
            crate::restriction::ENZYMES
                .iter()
                .position(|e| e.name == name)
                .unwrap()
        });
        let conditions = DigestConditions::default();
        let result =
            crate::restriction::digest(&source.sequence, false, &enzymes, conditions).unwrap();
        let mut products = digest_products(
            &source,
            &result.fragments,
            &["EcoRI".into(), "KpnI".into()],
            conditions,
        )
        .unwrap();
        assert_eq!(products.len(), 3);
        for product in &products {
            product.validate().unwrap();
            assert!(product.source.is_none());
            assert_eq!(product.methylation, source.methylation);
        }
        let middle = &mut products[1];
        let (prefix, suffix, _, _) = middle.duplex_ends.as_ref().unwrap().overhang_lengths();
        assert!(prefix + suffix < middle.sequence.len());
        middle.replace(prefix, prefix + 1, "G").unwrap();
        assert_eq!(middle.methylation, None);
        middle.validate().unwrap();
        let unchanged = middle.clone();
        assert!(middle.replace(0, 1, "A").is_err());
        assert_eq!(*middle, unchanged);
    }
}
