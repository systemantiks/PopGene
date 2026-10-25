//! Display subdivisions never replace the feature's biological location.
use super::{Feature, FeatureDirection, rgb_color};
use crate::Segment;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FeatureDisplaySegment {
    pub range: Segment,
    pub color_rgb: Option<u32>,
    pub name: Option<String>,
}

pub(super) struct SegmentMetadata {
    pub direction: Option<FeatureDirection>,
    pub spans: Vec<FeatureDisplaySegment>,
}

fn number(input: &str) -> Option<(usize, &str)> {
    let input = input.trim_start();
    let length = input.bytes().take_while(u8::is_ascii_digit).count();
    (length > 0).then_some(())?;
    Some((input[..length].parse().ok()?, &input[length..]))
}

fn entry(input: &str) -> Option<(Segment, u32, &str)> {
    let (start, rest) = number(input)?;
    let (end, rest) = number(rest.trim_start().strip_prefix("..")?)?;
    if start == 0 || start > end || end > crate::MAX_BASES {
        return None;
    }
    let rest = rest.trim_start().strip_prefix('/')?.trim_start();
    let color = rest.get(..7)?;
    if !color.starts_with('#')
        || !rest
            .get(7..)?
            .starts_with(|c: char| c.is_whitespace() || c == '/')
            && rest.len() != 7
    {
        return None;
    }
    Some((
        Segment {
            start: start - 1,
            end,
        },
        rgb_color(color)?,
        &rest[7..],
    ))
}

fn ordered_runs(parts: impl IntoIterator<Item = Segment>) -> Vec<Segment> {
    let mut result: Vec<Segment> = vec![];
    for part in parts {
        if let Some(last) = result.last_mut()
            && part.start == last.end
        {
            last.end = part.end;
        } else {
            result.push(part);
        }
    }
    result
}

fn parse(note: &str) -> Option<SegmentMetadata> {
    // A malformed or excessively large display note stays ordinary qualifier text.
    if note.len() > 1024 * 1024 || note.contains("Cleavage site") {
        return None;
    }
    let (direction, rest) = [
        (
            "This forward directional feature has ",
            Some(FeatureDirection::Forward),
        ),
        (
            "This reverse directional feature has ",
            Some(FeatureDirection::Reverse),
        ),
        (
            "This bidirectional feature has ",
            Some(FeatureDirection::Bidirectional),
        ),
        ("This feature has ", None),
    ]
    .into_iter()
    .find_map(|(prefix, direction)| note.strip_prefix(prefix).map(|rest| (direction, rest)))?;
    let (count, rest) = rest.split_once(" segments:")?;
    let count: usize = count.parse().ok()?;
    if count == 0 || count > 20_000 {
        return None;
    }
    // Find entry boundaries once. Rescanning the remaining note for every
    // segment would be quadratic for large imported feature tables.
    let rest = rest.trim_start();
    let mut rows = vec![];
    for (colon, _) in rest.match_indices(':') {
        let start = rest[..colon]
            .bytes()
            .rposition(|byte| !byte.is_ascii_digit())
            .map_or(0, |index| index + 1);
        if start == colon || (start > 0 && !rest[..start].chars().next_back()?.is_whitespace()) {
            continue;
        }
        if entry(&rest[colon + 1..]).is_some() {
            rows.push((start, colon, rest[start..colon].parse::<usize>().ok()?));
            if rows.len() > count {
                return None;
            }
        }
    }
    if rows.len() != count || rows.first()?.0 != 0 {
        return None;
    }
    let mut spans = Vec::with_capacity(count);
    for (index, (_, colon, ordinal)) in rows.iter().enumerate() {
        if *ordinal != index + 1 {
            return None;
        }
        let end = rows.get(index + 1).map_or(rest.len(), |row| row.0);
        let (range, color, name) = entry(&rest[colon + 1..end])?;
        let name = name.trim();
        let name = if name.is_empty() {
            None
        } else {
            let name = name.strip_prefix('/')?.trim();
            if name.is_empty() {
                None
            } else {
                Some(name.to_owned())
            }
        };
        spans.push(FeatureDisplaySegment {
            range,
            color_rgb: Some(color),
            name,
        });
    }
    Some(SegmentMetadata { direction, spans })
}

fn note(metadata: &SegmentMetadata) -> String {
    let direction = match metadata.direction {
        Some(FeatureDirection::Forward) => "forward directional ",
        Some(FeatureDirection::Reverse) => "reverse directional ",
        Some(FeatureDirection::Bidirectional) => "bidirectional ",
        _ => "",
    };
    let mut result = format!(
        "This {direction}feature has {} segments:",
        metadata.spans.len()
    );
    for (index, span) in metadata.spans.iter().enumerate() {
        result.push_str(&format!(
            " {}: {} .. {} / #{:06x}",
            index + 1,
            span.range.start + 1,
            span.range.end,
            span.color_rgb.unwrap()
        ));
        if let Some(name) = &span.name {
            result.push_str(&format!(" / {name}"));
        }
    }
    result
}

impl Feature {
    /// The validated terminal conventional note, for writers that need to keep
    /// it last when inserting a separate annotation-review note.
    pub fn segment_display_note(&self) -> Option<&str> {
        self.segment_metadata()?;
        self.qualifiers
            .iter()
            .rev()
            .find(|(key, _)| key == "note")?
            .1
            .as_deref()
    }

    pub(super) fn segment_metadata(&self) -> Option<SegmentMetadata> {
        let raw = self
            .qualifiers
            .iter()
            .rev()
            .find(|(key, _)| key == "note")?
            .1
            .as_deref()?;
        let parsed = parse(raw)?;
        // Display metadata cannot hide, add, reorder, or change biological bases.
        (ordered_runs(parsed.spans.iter().map(|p| p.range.clone()))
            == ordered_runs(self.parts.clone()))
        .then_some(parsed)
    }

    /// Display subdivisions, or one span per biological part if no validated
    /// metadata is present. The feature's extraction/translation is unchanged.
    pub fn display_segments(&self) -> Vec<FeatureDisplaySegment> {
        if let Some(metadata) = self.segment_metadata() {
            return metadata.spans;
        }
        let color = self.display_color_rgb();
        self.parts
            .iter()
            .map(|range| FeatureDisplaySegment {
                range: range.clone(),
                color_rgb: color,
                name: if self.parts.len() == 1 {
                    self.qualifier("segment_name").map(str::to_owned)
                } else {
                    None
                },
            })
            .collect()
    }

    pub(super) fn reverse_display_segments(&mut self, length: usize) {
        let Some(mut metadata) = self.segment_metadata() else {
            return;
        };
        for span in &mut metadata.spans {
            span.range = Segment {
                start: length - span.range.end,
                end: length - span.range.start,
            };
        }
        metadata.spans.reverse();
        metadata.direction = match metadata.direction {
            Some(FeatureDirection::Forward) => Some(FeatureDirection::Reverse),
            Some(FeatureDirection::Reverse) => Some(FeatureDirection::Forward),
            other => other,
        };
        self.replace_segment_note(&metadata);
    }

    pub(crate) fn shift_display_segments(&mut self, start: usize, end: usize, replacement: usize) {
        let Some(mut metadata) = self.segment_metadata() else {
            return;
        };
        let delta = replacement as isize - (end - start) as isize;
        let shift = |v: usize| (v as isize + delta) as usize;
        for span in &mut metadata.spans {
            let part = &mut span.range;
            if part.end <= start {
                continue;
            }
            if part.start >= end {
                part.start = shift(part.start);
                part.end = shift(part.end);
                continue;
            }
            part.start = part.start.min(start);
            part.end = if part.end > end {
                shift(part.end)
            } else {
                start + replacement
            };
        }
        metadata
            .spans
            .retain(|span| span.range.start < span.range.end);
        if !metadata.spans.is_empty() {
            self.replace_segment_note(&metadata);
        }
    }

    fn replace_segment_note(&mut self, metadata: &SegmentMetadata) {
        if let Some((_, value)) = self
            .qualifiers
            .iter_mut()
            .rev()
            .find(|(key, _)| key == "note")
        {
            *value = Some(note(metadata));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn sample() -> Feature {
        Feature { id: 1, label: "split".into(), kind: "misc_feature".into(), parts: vec![Segment { start: 100, end: 400 }], reverse: false, edited: false,
            qualifiers: vec![("note".into(), Some("This forward directional feature has 2 segments: 1: 101 .. 250 / #ff0000 / Red domain 2: 251 .. 400 / #0000ff / Blue domain".into()))] }
    }
    #[test]
    fn subdivisions_preserve_biological_location_and_read_names_colors() {
        let f = sample();
        let original = f.clone();
        let spans = f.display_segments();
        assert_eq!(spans.len(), 2);
        assert_eq!(
            spans[0].range,
            Segment {
                start: 100,
                end: 250
            }
        );
        assert_eq!(spans[1].name.as_deref(), Some("Blue domain"));
        assert_eq!(spans[1].color_rgb, Some(0x0000ff));
        assert_eq!(f.display_direction(), Some(FeatureDirection::Forward));
        assert_eq!(f, original);
    }
    #[test]
    fn invalid_or_partial_metadata_never_hides_bases() {
        for invalid in [
            "1: 101 .. 250 / #ff0000 2: 251 .. 399 / #0000ff",
            "1: 101 .. 250 / red 2: 251 .. 400 / #0000ff",
            "1: 101 .. 250 / #ff0000 3: 251 .. 400 / #0000ff",
        ] {
            let mut f = sample();
            f.qualifiers[0].1 = Some(format!("This feature has 2 segments: {invalid}"));
            assert_eq!(f.display_segments().len(), 1);
            assert_eq!(f.display_segments()[0].range, f.parts[0]);
        }
    }
    #[test]
    fn origin_order_overlaps_and_unsupported_cleavage_notes_are_not_reinterpreted() {
        let mut feature = sample();
        feature.parts = vec![
            Segment {
                start: 250,
                end: 400,
            },
            Segment {
                start: 100,
                end: 250,
            },
        ];
        assert!(feature.segment_metadata().is_none());
        feature.parts = vec![
            Segment {
                start: 100,
                end: 300,
            },
            Segment {
                start: 200,
                end: 400,
            },
        ];
        assert!(feature.segment_metadata().is_none());
        let mut feature = sample();
        feature.qualifiers[0]
            .1
            .as_mut()
            .unwrap()
            .push_str(" Cleavage site after base 300");
        assert!(feature.segment_metadata().is_none());
    }
    #[test]
    fn reverse_and_insert_keep_names_attached_to_their_bases() {
        let mut d = crate::Document::new("test", &"A".repeat(1000)).unwrap();
        d.features.push(sample());
        let original = d.features[0].display_segments();
        d.revcomp().unwrap();
        let spans = d.features[0].display_segments();
        assert_eq!(spans[0].name.as_deref(), Some("Blue domain"));
        assert_eq!(
            spans[0].range,
            Segment {
                start: 600,
                end: 750
            }
        );
        assert_eq!(
            d.features[0].display_direction(),
            Some(FeatureDirection::Reverse)
        );
        d.revcomp().unwrap();
        assert_eq!(d.features[0].display_segments(), original);
        d.replace(0, 0, "ACGT").unwrap();
        assert_eq!(
            d.features[0].display_segments()[0].range,
            Segment {
                start: 104,
                end: 254
            }
        );
        assert_eq!(
            d.features[0].parts,
            vec![Segment {
                start: 104,
                end: 404
            }]
        );
    }
}
