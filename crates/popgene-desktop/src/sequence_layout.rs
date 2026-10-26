//! Coordinate-only helpers for the virtualized, annotated sequence viewer.
//!
//! A sorted interval tree finds only annotation parts touching a visible row.
//! Calculated CDS translations are projected from their recorded reading frame.
//! Complete reading-frame codons are calculated only for visible rows.
use popgene_core::translation::translate_cds;
use popgene_core::{Document, FeatureDisplaySegment, FeatureSegmentStyle, Molecule, Segment};

/// Continuous symbol spacing, excluding the view's left/right gutters.
pub fn base_x(column: usize, cell_width: f32) -> f32 {
    column as f32 * cell_width
}

/// Sparse, disjoint color runs. The last imported range wins overlapping colors
/// on the same strand. Building once avoids scanning up to100k ranges per glyph.
#[derive(Default)]
pub struct SequenceColorIndex {
    strands: [Vec<(Segment, u32)>; 2],
}
impl SequenceColorIndex {
    pub fn new(document: &Document) -> Self {
        let mut result = Self::default();
        for strand in 0..2 {
            let mut events = vec![];
            for (index, color) in document.sequence_colors.iter().enumerate() {
                if color.reverse == (strand == 1) {
                    events.push((color.range.start, true, index));
                    events.push((color.range.end, false, index));
                }
            }
            events.sort_unstable();
            let mut active = std::collections::BTreeSet::new();
            let mut i = 0;
            while i < events.len() {
                let coordinate = events[i].0;
                while i < events.len() && events[i].0 == coordinate {
                    let (_, start, index) = events[i];
                    if start {
                        active.insert(index);
                    } else {
                        active.remove(&index);
                    }
                    i += 1;
                }
                if let Some(&index) = active.last()
                    && let Some(&(end, _, _)) = events.get(i)
                    && coordinate < end
                {
                    let color = document.sequence_colors[index].color_rgb;
                    if let Some((last, previous)) = result.strands[strand].last_mut()
                        && last.end == coordinate
                        && *previous == color
                    {
                        last.end = end;
                    } else {
                        result.strands[strand].push((
                            Segment {
                                start: coordinate,
                                end,
                            },
                            color,
                        ));
                    }
                }
            }
        }
        result
    }
    pub fn color(&self, position: usize, reverse: bool) -> Option<u32> {
        let runs = &self.strands[usize::from(reverse)];
        let index = runs.partition_point(|(range, _)| range.end <= position);
        runs.get(index)
            .filter(|(range, _)| range.start <= position)
            .map(|(_, color)| *color)
    }
}

/// Hit-test a symbol; gutters and the space after the final symbol are not bases.
pub fn base_at_x(x: f32, cell_width: f32, count: usize) -> Option<usize> {
    if !x.is_finite() || !cell_width.is_finite() || x < 0. || cell_width <= 0. {
        return None;
    }
    let index = (x / cell_width).floor() as usize;
    (index < count).then_some(index)
}

/// Keep ten-base ruler intervals while fitting a bounded number of symbols.
pub fn bases_per_row(width: f32, cell_width: f32, gutter: f32) -> usize {
    if !width.is_finite() || !cell_width.is_finite() || cell_width <= 0. {
        return 60;
    }
    let available = (width - gutter.max(0.)).max(0.);
    ((available / cell_width).floor() as usize / 10 * 10).clamp(20, 120)
}

pub fn complement(base: u8, molecule: Molecule) -> u8 {
    if molecule == Molecule::Protein {
        return base;
    }
    match base.to_ascii_uppercase() {
        b'A' if molecule == Molecule::Rna => b'U',
        b'A' => b'T',
        b'T' | b'U' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'D' => b'H',
        b'H' => b'D',
        b'V' => b'B',
        _ => b'N',
    }
}

/// Actual opposite-strand coverage projected onto the top sequence axis.
pub fn opposite_span(document: &Document) -> Segment {
    let n = document.sequence.len();
    if document.single_stranded || document.molecule == Molecule::Protein {
        return Segment { start: 0, end: 0 };
    }
    if let Some(ends) = &document.duplex_ends {
        let geometry = ends.geometry(n);
        Segment {
            start: geometry
                .bottom
                .start
                .saturating_sub(geometry.top.start)
                .min(n),
            end: geometry
                .bottom
                .end
                .saturating_sub(geometry.top.start)
                .min(n),
        }
    } else {
        Segment { start: 0, end: n }
    }
}
/// Only actual lower residues are returned for an explicitly stored duplex.
pub fn opposite_base(document: &Document, top_position: usize) -> Option<u8> {
    let span = opposite_span(document);
    if !(span.start..span.end).contains(&top_position) {
        return None;
    }
    if let Some(ends) = &document.duplex_ends {
        let geometry = ends.geometry(document.sequence.len());
        let axis = geometry.top.start + top_position;
        ends.opposite_sequence
            .as_bytes()
            .get(geometry.bottom.end - 1 - axis)
            .copied()
    } else {
        document
            .sequence
            .as_bytes()
            .get(top_position)
            .map(|base| complement(*base, document.molecule))
    }
}
pub fn opposite_quality(document: &Document, top_position: usize) -> Option<u8> {
    let span = opposite_span(document);
    if !(span.start..span.end).contains(&top_position) {
        return None;
    }
    let position = if let Some(ends) = &document.duplex_ends {
        let geometry = ends.geometry(document.sequence.len());
        geometry.bottom.end - 1 - (geometry.top.start + top_position)
    } else {
        document.sequence.len() - 1 - top_position
    };
    document
        .opposite_sequence_qualities
        .as_ref()?
        .get(position)
        .copied()
        .flatten()
}
pub fn duplex_display_notice(document: &Document) -> Option<String> {
    let ends = document.duplex_ends.as_ref()?;
    let geometry = ends.geometry(document.sequence.len());
    Some(format!(
        "Top-strand coordinates shown. Physical duplex extent: {} positions; {} bottom-only residues and {} full-extent annotations are outside this view. Blank lower cells have no opposite strand. Use Analyze → Full duplex strands & annotations to inspect both complete strands.",
        geometry.extent,
        geometry.extent - document.sequence.len(),
        document.duplex_annotations.features.len()
    ))
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FeatureSpan {
    pub feature_index: usize,
    pub feature_id: u64,
    /// Half-open coordinates relative to the requested row.
    pub start: usize,
    pub end: usize,
    pub reverse: bool,
    /// The part is clipped at this edge or continues to another ordered part.
    /// Terminal arrows should be drawn only when the relevant flag is false.
    pub continued_left: bool,
    pub continued_right: bool,
    pub lane: usize,
    pub color_rgb: Option<u32>,
    pub segment_name: Option<String>,
    pub style: FeatureSegmentStyle,
    /// Interbase boundaries relative to the row, not symbols to translate.
    pub cleavage_sites: Vec<usize>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TranslationGlyph {
    pub feature_index: usize,
    pub feature_id: u64,
    /// The codon's central base, relative to the row, where its letter is drawn.
    pub base: usize,
    /// The codon interval clipped to the row.
    pub start: usize,
    pub end: usize,
    pub amino_acid: char,
    pub lane: usize,
    pub genetic_code: u16,
    /// Absolute forward coordinates in coding order, including joined codons.
    pub codon_positions: [usize; 3],
    pub split_codon: bool,
    pub nucleotide_count: usize,
    pub exception: bool,
}

#[derive(Clone, Default, Debug)]
pub struct RowLayout {
    pub spans: Vec<FeatureSpan>,
    pub translations: Vec<TranslationGlyph>,
    /// Lanes containing an enabled calculated CDS translation track.
    pub translation_lanes: Vec<usize>,
    pub lanes: usize,
    /// Number of annotation parts omitted because all allowed lanes overlapped.
    pub omitted: usize,
}

#[derive(Clone, Default, Debug, PartialEq, Eq)]
pub struct RowSummary {
    pub lanes: usize,
    pub omitted: usize,
    pub translation_lanes: Vec<usize>,
}

#[derive(Clone, Debug)]
struct IndexedPart {
    feature_index: usize,
    part_index: usize,
    start: usize,
    end: usize,
    /// Forward concatenation offset, before applying the feature's strand.
    offset: usize,
    gap_index: Option<usize>,
    point_boundary: Option<usize>,
}

#[derive(Clone, Debug)]
struct Translation {
    amino_acids: Vec<u8>,
    phase: usize,
    genetic_code: u16,
    exception_codons: Vec<usize>,
    terminal_partial_stop: bool,
}

#[derive(Clone, Debug)]
struct IndexedFeature {
    id: u64,
    reverse: bool,
    primer: bool,
    parts: Vec<Segment>,
    /// Forward concatenation offset plus styling, preserving repeated locations.
    display_parts: Vec<(usize, FeatureDisplaySegment)>,
    display_gaps: Vec<FeatureDisplaySegment>,
    cleavage_sites: Vec<usize>,
    offsets: Vec<usize>,
    length: usize,
    translation: Option<Translation>,
}

/// An annotation index independent of current row wrapping and display filters.
#[derive(Default)]
pub struct SequenceLayout {
    features: Vec<IndexedFeature>,
    parts: Vec<IndexedPart>,
    /// Segment-tree maxima; the root is at index 1.
    max_ends: Vec<usize>,
    document_length: usize,
    pub translation_notices: Vec<String>,
    pub omitted_translations: usize,
}

impl SequenceLayout {
    pub fn new(document: &Document) -> Self {
        Self::build(document, true)
    }
    /// Build only annotation geometry, avoiding synthetic sequence translations
    /// when the caller displays separate physical strands on a duplex axis.
    pub fn annotations_only(document: &Document) -> Self {
        Self::build(document, false)
    }
    fn build(document: &Document, calculate_translations: bool) -> Self {
        let mut result = Self {
            document_length: document.sequence.len(),
            ..Self::default()
        };
        let mut translation_budget = popgene_core::MAX_BASES;
        for (feature_index, feature) in document.features.iter().enumerate() {
            let valid = feature
                .parts
                .iter()
                .all(|part| part.start < part.end && part.end <= result.document_length);
            let mut offsets = Vec::with_capacity(feature.parts.len());
            let mut length = 0;
            for (part_index, part) in feature.parts.iter().enumerate() {
                offsets.push(length);
                if valid {
                    result.parts.push(IndexedPart {
                        feature_index,
                        part_index,
                        start: part.start,
                        end: part.end,
                        offset: length,
                        gap_index: None,
                        point_boundary: None,
                    });
                    length += part.end - part.start;
                }
            }
            let display_gaps = feature
                .display_gaps()
                .into_iter()
                .filter(|span| {
                    span.range.start < span.range.end && span.range.end <= result.document_length
                })
                .collect::<Vec<_>>();
            if valid {
                for (gap_index, gap) in display_gaps.iter().enumerate() {
                    result.parts.push(IndexedPart {
                        feature_index,
                        part_index: feature.parts.len() + gap_index,
                        start: gap.range.start,
                        end: gap.range.end,
                        offset: 0,
                        gap_index: Some(gap_index),
                        point_boundary: None,
                    });
                }
            }
            if let Some(location) = &feature.location {
                for boundary in location.local_boundaries() {
                    if boundary > 0 && boundary <= result.document_length {
                        let at = boundary.min(result.document_length - 1);
                        result.parts.push(IndexedPart {
                            feature_index,
                            part_index: usize::MAX,
                            start: at,
                            end: at + 1,
                            offset: 0,
                            gap_index: None,
                            point_boundary: Some(boundary),
                        });
                    }
                }
            }
            let translation = if calculate_translations
                && valid
                && document.molecule != Molecule::Protein
                && feature.kind.eq_ignore_ascii_case("CDS")
            {
                let calculated = if length <= translation_budget {
                    translation_budget -= length;
                    translate_cds(document, feature)
                } else {
                    Err(popgene_core::translation::TranslationError(
                        "document CDS translation budget exceeded (10 million bases)".into(),
                    ))
                };
                match calculated {
                    Ok(translation) => Some(Translation {
                        amino_acids: translation.amino_acids,
                        phase: translation.phase,
                        genetic_code: translation.genetic_code,
                        exception_codons: translation.exception_codons,
                        terminal_partial_stop: translation.terminal_partial_stop,
                    }),
                    Err(error) => {
                        result.omitted_translations += 1;
                        if result.translation_notices.len() < 8 {
                            result
                                .translation_notices
                                .push(format!("{}: {error}", feature.label));
                        }
                        None
                    }
                }
            } else {
                None
            };
            result.features.push(IndexedFeature {
                id: feature.id,
                reverse: feature.reverse,
                primer: feature.kind == "primer_bind",
                parts: feature.parts.clone(),
                display_parts: feature
                    .display_segments()
                    .into_iter()
                    .scan(0, |offset, display| {
                        let start = *offset;
                        *offset += display.range.end - display.range.start;
                        Some((start, display))
                    })
                    .collect(),
                offsets,
                display_gaps,
                cleavage_sites: feature.display_cleavage_sites(),
                length,
                translation,
            });
        }
        result
            .parts
            .sort_by_key(|part| (part.start, part.end, part.feature_index));
        result
            .max_ends
            .resize(result.parts.len().saturating_mul(4).max(2), 0);
        if !result.parts.is_empty() {
            result.build_max_ends(1, 0, result.parts.len());
        }
        result
    }

    fn build_max_ends(&mut self, node: usize, lo: usize, hi: usize) -> usize {
        let max = if hi - lo == 1 {
            self.parts[lo].end
        } else {
            let middle = lo + (hi - lo) / 2;
            self.build_max_ends(node * 2, lo, middle)
                .max(self.build_max_ends(node * 2 + 1, middle, hi))
        };
        self.max_ends[node] = max;
        max
    }

    fn touching<'a>(
        &'a self,
        node: usize,
        interval: std::ops::Range<usize>,
        range: std::ops::Range<usize>,
        output: &mut Vec<&'a IndexedPart>,
    ) {
        if interval.is_empty()
            || self.max_ends[node] <= range.start
            || self.parts[interval.start].start >= range.end
        {
            return;
        }
        if interval.len() == 1 {
            output.push(&self.parts[interval.start]);
        } else {
            let middle = interval.start + interval.len() / 2;
            self.touching(node * 2, interval.start..middle, range.clone(), output);
            self.touching(node * 2 + 1, middle..interval.end, range, output);
        }
    }

    /// Calculate only the tracks needed for a visible row.
    pub fn row(
        &self,
        start: usize,
        end: usize,
        max_lanes: usize,
        show_features: bool,
        show_primers: bool,
        show_translations: bool,
    ) -> RowLayout {
        self.layout_row(
            start..end,
            max_lanes,
            show_features,
            show_primers,
            show_translations,
            true,
        )
    }

    /// Measure row tracks without constructing spans or expanding amino-acid glyphs.
    pub fn row_summary(
        &self,
        start: usize,
        end: usize,
        max_lanes: usize,
        show_features: bool,
        show_primers: bool,
        show_translations: bool,
    ) -> RowSummary {
        let row = self.layout_row(
            start..end,
            max_lanes,
            show_features,
            show_primers,
            show_translations,
            false,
        );
        RowSummary {
            lanes: row.lanes,
            omitted: row.omitted,
            translation_lanes: row.translation_lanes,
        }
    }

    fn layout_row(
        &self,
        range: std::ops::Range<usize>,
        max_lanes: usize,
        show_features: bool,
        show_primers: bool,
        show_translations: bool,
        content: bool,
    ) -> RowLayout {
        let start = range.start;
        let end = range.end.min(self.document_length);
        if start >= end || self.parts.is_empty() {
            return RowLayout::default();
        }
        let mut touching = Vec::new();
        self.touching(1, 0..self.parts.len(), start..end, &mut touching);
        // Longer spans first at identical starts keep encompassing features legible.
        touching.sort_by_key(|part| {
            (
                part.start.max(start),
                std::cmp::Reverse(part.end),
                part.feature_index,
                part.part_index,
            )
        });
        let mut result = RowLayout::default();
        let mut lane_ends = Vec::<usize>::new();
        for part in touching {
            let feature = &self.features[part.feature_index];
            let visible_span = if feature.primer {
                show_primers
            } else {
                show_features
            };
            let visible_translation = show_translations
                && feature.translation.is_some()
                && part.gap_index.is_none()
                && part.point_boundary.is_none();
            if !visible_span && !visible_translation {
                continue;
            }
            let left = part.start.max(start);
            let right = part.end.min(end);
            let lane = if let Some(lane) = lane_ends.iter().position(|end| *end <= left) {
                lane
            } else if lane_ends.len() < max_lanes {
                lane_ends.push(0);
                lane_ends.len() - 1
            } else {
                result.omitted += 1;
                continue;
            };
            lane_ends[lane] = right;
            if visible_translation && !result.translation_lanes.contains(&lane) {
                result.translation_lanes.push(lane);
            }
            if content
                && visible_span
                && let Some(boundary) = part.point_boundary
            {
                result.spans.push(FeatureSpan {
                    feature_index: part.feature_index,
                    feature_id: feature.id,
                    start: boundary - start,
                    end: boundary - start,
                    reverse: feature.reverse,
                    continued_left: true,
                    continued_right: true,
                    lane,
                    color_rgb: None,
                    segment_name: None,
                    style: FeatureSegmentStyle::Line,
                    cleavage_sites: vec![boundary - start],
                });
            }
            if content
                && visible_span
                && let Some(gap_index) = part.gap_index
            {
                let gap = &feature.display_gaps[gap_index];
                result.spans.push(FeatureSpan {
                    feature_index: part.feature_index,
                    feature_id: feature.id,
                    start: left - start,
                    end: right - start,
                    reverse: feature.reverse,
                    continued_left: true,
                    continued_right: true,
                    lane,
                    color_rgb: gap.color_rgb,
                    segment_name: gap.name.clone(),
                    style: FeatureSegmentStyle::Gap,
                    cleavage_sites: vec![],
                });
            }
            if content && visible_span && part.gap_index.is_none() && part.point_boundary.is_none()
            {
                // Match occurrences in the ordered concatenation, not just
                // genomic coordinates: joined parts may overlap or repeat.
                let visible_offset_start = part.offset + left - part.start;
                let visible_offset_end = part.offset + right - part.start;
                for (display_offset, display) in &feature.display_parts {
                    let overlap_start = visible_offset_start.max(*display_offset);
                    let overlap_end = visible_offset_end
                        .min(display_offset + display.range.end - display.range.start);
                    if overlap_start >= overlap_end {
                        continue;
                    }
                    let display_left = part.start + overlap_start - part.offset;
                    let display_right = part.start + overlap_end - part.offset;
                    result.spans.push(FeatureSpan {
                        feature_index: part.feature_index,
                        feature_id: feature.id,
                        start: display_left - start,
                        end: display_right - start,
                        reverse: feature.reverse,
                        continued_left: part.start < display_left || part.part_index > 0,
                        continued_right: part.end > display_right
                            || part.part_index + 1 < feature.parts.len(),
                        lane,
                        color_rgb: display.color_rgb,
                        segment_name: display.name.clone(),
                        style: display.style,
                        cleavage_sites: feature
                            .cleavage_sites
                            .iter()
                            .filter(|&&p| p >= display_left && p <= display_right)
                            .map(|p| p - start)
                            .collect(),
                    });
                }
            }
            if content
                && show_translations
                && part.gap_index.is_none()
                && part.point_boundary.is_none()
                && let Some(translation) = &feature.translation
            {
                for coordinate in left..right {
                    let forward_offset = part.offset + coordinate - part.start;
                    let coding_offset = if feature.reverse {
                        feature.length - 1 - forward_offset
                    } else {
                        forward_offset
                    };
                    if coding_offset < translation.phase {
                        continue;
                    }
                    let in_cds = coding_offset - translation.phase;
                    if in_cds % 3 != 1
                        && !(translation.terminal_partial_stop
                            && in_cds % 3 == 0
                            && coding_offset + 1 == feature.length)
                    {
                        continue;
                    }
                    let Some(&amino_acid) = translation.amino_acids.get(in_cds / 3) else {
                        continue;
                    };
                    let codon_offset = translation.phase + (in_cds / 3) * 3;
                    let nucleotide_count = (feature.length - codon_offset).min(3);
                    let codon_positions = std::array::from_fn(|index| {
                        feature
                            .coordinate(codon_offset + index.min(nucleotide_count - 1))
                            .unwrap()
                    });
                    let before = codon_positions[0];
                    let after = codon_positions[2];
                    // Never stretch a translation across an intron or circular origin.
                    let contiguous = codon_positions[..nucleotide_count].windows(2).all(|pair| {
                        if feature.reverse {
                            pair[1].checked_add(1) == Some(pair[0])
                        } else {
                            pair[0].checked_add(1) == Some(pair[1])
                        }
                    });
                    result.translations.push(TranslationGlyph {
                        feature_index: part.feature_index,
                        feature_id: feature.id,
                        base: coordinate - start,
                        start: if contiguous {
                            before.min(after).max(start)
                        } else {
                            coordinate
                        } - start,
                        end: if contiguous {
                            (before.max(after) + 1).min(end)
                        } else {
                            coordinate + 1
                        } - start,
                        amino_acid: amino_acid as char,
                        lane,
                        genetic_code: translation.genetic_code,
                        codon_positions,
                        nucleotide_count,
                        exception: translation
                            .exception_codons
                            .binary_search(&(in_cds / 3))
                            .is_ok(),
                        split_codon: !contiguous,
                    });
                }
            }
        }
        result.lanes = lane_ends.len();
        result
    }
}

impl IndexedFeature {
    fn coordinate(&self, coding_offset: usize) -> Option<usize> {
        if coding_offset >= self.length {
            return None;
        }
        let offset = if self.reverse {
            self.length - 1 - coding_offset
        } else {
            coding_offset
        };
        let index = self
            .offsets
            .partition_point(|start| *start <= offset)
            .checked_sub(1)?;
        self.parts
            .get(index)
            .map(|part| part.start + offset - self.offsets[index])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popgene_core::Feature;

    fn feature(id: u64, parts: &[(usize, usize)], reverse: bool) -> Feature {
        Feature {
            location: None,
            id,
            label: format!("feature-{id}"),
            kind: "CDS".into(),
            parts: parts
                .iter()
                .map(|&(start, end)| Segment { start, end })
                .collect(),
            reverse,
            qualifiers: vec![("translation".into(), Some("MKV*".into()))],
            edited: false,
        }
    }
    fn document(features: Vec<Feature>) -> Document {
        let mut document = Document::new("Layout", &"A".repeat(120)).unwrap();
        document.features = features;
        document
    }

    #[test]
    fn complement_preserves_ambiguity_and_molecule_direction() {
        let dna: Vec<_> = b"ACGTRYSWKMBDHVN"
            .iter()
            .map(|base| complement(*base, Molecule::Dna))
            .collect();
        assert_eq!(dna, b"TGCAYRSWMKVHDBN");
        assert_eq!(complement(b'A', Molecule::Rna), b'U');
        assert_eq!(complement(b'U', Molecule::Rna), b'A');
        assert_eq!(complement(b'M', Molecule::Protein), b'M');
    }

    #[test]
    fn width_and_hit_testing_exclude_gutters_and_last_boundary() {
        assert_eq!(bases_per_row(780., 10., 100.), 60);
        assert_eq!(bases_per_row(200., 10., 100.), 20);
        assert_eq!(bases_per_row(5000., 10., 100.), 120);
        assert_eq!(base_x(10, 10.5), 105.);
        assert_eq!(base_at_x(104.9, 10.5, 60), Some(9));
        assert_eq!(base_at_x(105., 10.5, 60), Some(10));
        assert_eq!(base_at_x(630., 10.5, 60), None);
        assert_eq!(base_at_x(-0.1, 10.5, 60), None);
        assert_eq!(base_at_x(f32::NAN, 10.5, 60), None);
    }

    #[test]
    fn display_subdivisions_keep_names_colors_and_coding_frame() {
        let mut feature = feature(1, &[(0, 12)], false);
        feature.qualifiers.push(("note".into(), Some("This forward directional feature has 2 segments: 1: 1 .. 5 / #ff0000 / Red domain 2: 6 .. 12 / #0000ff / Blue domain".into())));
        let layout = SequenceLayout::new(&document(vec![feature]));
        let row = layout.row(0, 12, 2, true, true, true);
        assert_eq!(row.spans.len(), 2);
        assert_eq!(
            (row.spans[0].start, row.spans[0].end, row.spans[0].color_rgb),
            (0, 5, Some(0xff0000))
        );
        assert_eq!(row.spans[1].segment_name.as_deref(), Some("Blue domain"));
        assert_eq!(
            row.translations
                .iter()
                .map(|glyph| glyph.base)
                .collect::<Vec<_>>(),
            vec![1, 4, 7, 10]
        );
    }

    #[test]
    fn repeated_locations_keep_segment_occurrences_and_terminal_edges() {
        let mut feature = feature(1, &[(0, 6), (0, 6)], false);
        feature.qualifiers.push(("note".into(), Some("This forward directional feature has 2 segments: 1: 1 .. 6 / #ff0000 / First occurrence 2: 1 .. 6 / #0000ff / Second occurrence".into())));
        let doc = document(vec![feature]);
        let original = doc.clone();
        let row = SequenceLayout::new(&doc).row(0, 6, 2, true, true, true);
        assert_eq!(row.spans.len(), 2);
        for (index, (color, name)) in [
            (0xff0000, "First occurrence"),
            (0x0000ff, "Second occurrence"),
        ]
        .iter()
        .enumerate()
        {
            let span = &row.spans[index];
            assert_eq!((span.start, span.end, span.lane), (0, 6, index));
            assert_eq!(span.color_rgb, Some(*color));
            assert_eq!(span.segment_name.as_deref(), Some(*name));
        }
        assert!(!row.spans[0].continued_left);
        assert!(row.spans[0].continued_right);
        assert!(row.spans[1].continued_left);
        assert!(!row.spans[1].continued_right);
        assert_eq!(
            row.translations
                .iter()
                .map(|glyph| (glyph.lane, glyph.base))
                .collect::<Vec<_>>(),
            [(0, 1), (0, 4), (1, 1), (1, 4)]
        );
        assert_eq!(doc, original);
    }

    #[test]
    fn uncached_genbank_cds_calculates_both_strands_and_additional_code() {
        let document = popgene_formats::import(
            "calculated-translations.gb",
            include_bytes!("../../../validation/fixtures/regression/calculated-translations.gb"),
        )
        .unwrap();
        assert!(
            document
                .features
                .iter()
                .all(|feature| feature.qualifier("translation").is_none())
        );
        let layout = SequenceLayout::new(&document);
        let row = layout.row(0, 21, 4, true, true, true);
        let glyphs = |index| {
            row.translations
                .iter()
                .filter(|glyph| glyph.feature_index == index)
                .map(|glyph| (glyph.base, glyph.amino_acid))
                .collect::<Vec<_>>()
        };
        assert_eq!(glyphs(0), vec![(1, 'M'), (4, 'K'), (7, '*')]);
        assert_eq!(glyphs(1), vec![(13, '*'), (16, 'K'), (19, 'M')]);
        assert_eq!(glyphs(2), vec![(10, 'P')]);
        assert_eq!(layout.omitted_translations, 0);
        assert!(layout.translation_notices.is_empty());
    }

    #[test]
    fn recoded_and_partial_stop_glyphs_preserve_actual_coordinates() {
        for (sequence, location, center, count) in [("ATGT", "4", 3, 1), ("ATGTA", "4..5", 4, 2)] {
            let mut document = Document::new("partial", sequence).unwrap();
            let mut cds = feature(1, &[(0, sequence.len())], false);
            cds.qualifiers.push((
                "transl_except".into(),
                Some(format!("(pos:{location},aa:TERM)")),
            ));
            document.features.push(cds);
            let layout = SequenceLayout::new(&document);
            let row = layout.row(0, sequence.len(), 4, true, false, true);
            let stop = row
                .translations
                .iter()
                .find(|glyph| glyph.amino_acid == '*')
                .unwrap();
            assert_eq!(
                (stop.base, stop.nucleotide_count, stop.start, stop.end),
                (center, count, 3, sequence.len())
            );
            assert!(stop.exception);
        }
    }

    #[test]
    fn compounds_are_separate_clipped_and_tracks_are_bounded() {
        let doc = document(vec![
            feature(1, &[(110, 120), (0, 12)], true),
            feature(2, &[(3, 20)], false),
            feature(3, &[(6, 30)], false),
        ]);
        let layout = SequenceLayout::new(&doc);
        let row = layout.row(5, 15, 2, true, true, false);
        assert_eq!(row.lanes, 2);
        assert_eq!(row.omitted, 1);
        assert!(
            row.spans
                .iter()
                .all(|span| span.lane < 2 && span.start < span.end && span.end <= 10)
        );
        let origin = row.spans.iter().find(|span| span.feature_id == 1).unwrap();
        assert!(origin.reverse && origin.continued_left && !origin.continued_right);
        assert_eq!((origin.start, origin.end), (0, 7));
        let all = layout.row(0, 120, 5, true, true, false);
        assert_eq!(
            all.spans.iter().filter(|span| span.feature_id == 1).count(),
            2
        );
        let first_part = all
            .spans
            .iter()
            .find(|span| span.feature_id == 1 && span.start == 110)
            .unwrap();
        let last_part = all
            .spans
            .iter()
            .find(|span| span.feature_id == 1 && span.start == 0)
            .unwrap();
        assert!(!first_part.continued_left && first_part.continued_right);
        assert!(last_part.continued_left && !last_part.continued_right);
        assert!(layout.row(30, 110, 4, true, true, false).spans.is_empty());
        let touching = SequenceLayout::new(&document(vec![
            feature(10, &[(0, 10)], false),
            feature(11, &[(10, 20)], false),
        ]));
        assert_eq!(touching.row(0, 20, 1, true, true, false).omitted, 0);
    }

    #[test]
    fn translation_uses_recorded_frame_strand_and_compound_order() {
        let mut forward = feature(1, &[(10, 22)], false);
        forward
            .qualifiers
            .push(("codon_start".into(), Some("2".into())));
        let reverse = feature(2, &[(30, 42)], true);
        let origin = feature(3, &[(117, 120), (0, 9)], false);
        let layout = SequenceLayout::new(&document(vec![forward, reverse, origin]));
        let row = layout.row(0, 120, 6, true, true, true);
        let letters = |id| {
            row.translations
                .iter()
                .filter(|glyph| glyph.feature_id == id)
                .map(|glyph| (glyph.base, glyph.amino_acid))
                .collect::<Vec<_>>()
        };
        assert_eq!(letters(1), vec![(12, 'K'), (15, 'K'), (18, 'K')]);
        assert_eq!(letters(2), vec![(31, 'F'), (34, 'F'), (37, 'F'), (40, 'F')]);
        assert_eq!(letters(3), vec![(1, 'K'), (4, 'K'), (7, 'K'), (118, 'K')]);
    }

    #[test]
    fn calculated_translations_recompute_edits_and_reject_invalid_metadata() {
        let split = feature(1, &[(10, 12), (20, 30)], false);
        let adjacent = feature(2, &[(40, 42), (42, 52)], false);
        let mut edited = feature(3, &[(60, 72)], false);
        edited.edited = true;
        let mut invalid = feature(4, &[(80, 92)], false);
        invalid
            .qualifiers
            .push(("codon_start".into(), Some("4".into())));
        let mut primer = feature(5, &[(95, 107)], false);
        primer.kind = "primer_bind".into();
        let layout = SequenceLayout::new(&document(vec![split, adjacent, edited, invalid, primer]));
        let row = layout.row(0, 120, 8, true, false, true);
        let split_codon = row
            .translations
            .iter()
            .find(|glyph| glyph.feature_id == 1 && glyph.base == 11)
            .unwrap();
        assert!(split_codon.split_codon);
        assert_eq!(split_codon.codon_positions, [10, 11, 20]);
        assert_eq!((split_codon.start, split_codon.end), (11, 12));
        assert!(
            row.translations
                .iter()
                .any(|glyph| glyph.feature_id == 2 && glyph.base == 41 && glyph.amino_acid == 'K')
        );
        assert!(row.translations.iter().any(|glyph| glyph.feature_id == 3));
        assert!(row.translations.iter().all(|glyph| glyph.feature_id < 4));
        assert_eq!(layout.omitted_translations, 1);
        assert!(layout.translation_notices[0].contains("codon_start"));
        assert!(!row.spans.iter().any(|span| span.feature_id == 5));
        let primers = layout.row(0, 120, 8, false, true, false);
        assert_eq!(primers.spans.len(), 1);
        assert!(primers.translations.is_empty());
    }

    #[test]
    fn translation_only_layers_measure_identically_and_omit_only_visible_tracks() {
        let first = feature(1, &[(0, 12)], false);
        let second = feature(2, &[(0, 12)], false);
        let mut no_translation = feature(3, &[(0, 20)], false);
        no_translation.kind = "misc_feature".into();
        let mut primer = feature(4, &[(0, 20)], false);
        primer.kind = "primer_bind".into();
        let layout = SequenceLayout::new(&document(vec![first, second, no_translation, primer]));
        let translations = layout.row(0, 30, 1, false, false, true);
        assert!(translations.spans.is_empty());
        assert_eq!(translations.lanes, 1);
        assert_eq!(translations.omitted, 1);
        assert_eq!(translations.translations.len(), 4);
        assert!(
            translations
                .translations
                .iter()
                .all(|glyph| glyph.feature_id == 1 && glyph.lane == 0)
        );
        assert_eq!(
            layout.row_summary(0, 30, 1, false, false, true),
            RowSummary {
                lanes: 1,
                omitted: 1,
                translation_lanes: vec![0],
            }
        );
        assert_eq!(
            layout.row_summary(0, 30, 1, false, false, false),
            RowSummary::default()
        );
        assert_eq!(
            layout.row_summary(0, 30, 1, false, true, false),
            RowSummary {
                lanes: 1,
                omitted: 0,
                translation_lanes: vec![],
            }
        );
    }

    #[test]
    fn translated_codons_at_row_edges_are_clipped_without_repositioning() {
        let layout = SequenceLayout::new(&document(vec![feature(1, &[(0, 12)], false)]));
        let left = layout.row(1, 5, 2, true, true, true);
        assert_eq!(
            left.translations
                .iter()
                .map(|glyph| (glyph.base, glyph.start, glyph.end))
                .collect::<Vec<_>>(),
            vec![(0, 0, 2), (3, 2, 4)]
        );
        let right = layout.row(5, 12, 2, true, true, true);
        assert_eq!(
            right
                .translations
                .iter()
                .map(|glyph| glyph.base)
                .collect::<Vec<_>>(),
            vec![2, 5]
        );
    }
    #[test]
    fn annotation_line_gap_and_cleavage_styling_does_not_add_translated_bases() {
        let mut f = feature(1, &[(0, 6), (9, 15)], false);
        f.qualifiers.push(("note".into(),Some("This forward directional feature has 2 segments: 1: 1..6 / #------ / first 2: 10..15 / #ff0000 / second Cleavage site after base 3".into())));
        f.qualifiers
            .push(("snapgene_gap".into(), Some("7..9 / #777777".into())));
        let d = document(vec![f]);
        let layout = SequenceLayout::new(&d);
        let row = layout.row(0, 18, 4, true, true, true);
        assert!(
            row.spans
                .iter()
                .any(|s| s.style == FeatureSegmentStyle::Gap && s.start == 6 && s.end == 9)
        );
        assert!(
            row.spans
                .iter()
                .any(|s| s.style == FeatureSegmentStyle::Line)
        );
        assert!(row.spans.iter().any(|s| s.cleavage_sites == vec![3]));
        assert!(
            row.translations
                .iter()
                .all(|g| g.codon_positions.iter().all(|p| !(6..9).contains(p)))
        );
    }
    #[test]
    fn between_and_alternative_points_are_visible_without_invented_bases() {
        use popgene_core::FeatureLocation;
        let mut f = feature(1, &[], false);
        f.kind = "misc_feature".into();
        f.location = Some(FeatureLocation::OneOf(vec![
            FeatureLocation::Between { left: 2, right: 3 },
            FeatureLocation::Between { left: 5, right: 6 },
            FeatureLocation::Remote {
                accession: "OTHER.1".into(),
                location: Some(Box::new(FeatureLocation::Between { left: 8, right: 9 })),
            },
        ]));
        let d = document(vec![f]);
        let row = SequenceLayout::new(&d).row(0, 12, 4, true, true, true);
        assert_eq!(row.spans.len(), 2);
        assert!(row.spans.iter().all(|s| s.start == s.end));
        assert_eq!(
            row.spans
                .iter()
                .flat_map(|s| s.cleavage_sites.clone())
                .collect::<Vec<_>>(),
            vec![3, 6]
        );
        assert!(row.translations.is_empty());
    }
    #[test]
    fn custom_sequence_colors_are_sparse_strand_specific_and_last_range_wins() {
        use popgene_core::metadata::SequenceColor;
        let mut d = document(vec![]);
        d.sequence_colors = vec![
            SequenceColor {
                range: Segment { start: 1, end: 8 },
                reverse: false,
                color_rgb: 0xff0000,
            },
            SequenceColor {
                range: Segment { start: 3, end: 5 },
                reverse: false,
                color_rgb: 0x00ff00,
            },
            SequenceColor {
                range: Segment { start: 4, end: 7 },
                reverse: true,
                color_rgb: 0x0000ff,
            },
        ];
        let colors = SequenceColorIndex::new(&d);
        assert_eq!(colors.color(0, false), None);
        assert_eq!(colors.color(1, false), Some(0xff0000));
        assert_eq!(colors.color(4, false), Some(0x00ff00));
        assert_eq!(colors.color(5, false), Some(0xff0000));
        assert_eq!(colors.color(8, false), None);
        assert_eq!(colors.color(4, true), Some(0x0000ff));
        assert_eq!(colors.color(3, true), None);
        assert_eq!(colors.color(7, true), None);
        assert_eq!(colors.strands[0].len(), 3);
    }
}
