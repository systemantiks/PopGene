//! Bounded, read-only feature details shared by map and annotation lists.
use gpui_kit::{prelude::*, *};
use popgene_core::{Document, Feature, FeatureDirection, Molecule};
use std::fmt::Write;

pub(crate) fn truncate(value: &str, limit: usize) -> String {
    // Inspect only the preview plus one character. Imported qualifier values
    // can be very large; hovering must not traverse their full contents.
    let mut characters = value.chars();
    let mut preview = characters.by_ref().take(limit).collect::<String>();
    if characters.next().is_some() {
        preview.push('…');
    }
    preview
}
pub(crate) fn number(value: usize) -> String {
    let text = value.to_string();
    text.chars()
        .enumerate()
        .flat_map(|(i, c)| {
            if i > 0 && (text.len() - i).is_multiple_of(3) {
                vec![',', c]
            } else {
                vec![c]
            }
        })
        .collect()
}
pub(crate) fn direction(feature: &Feature, molecule: Molecule) -> &'static str {
    if molecule == Molecule::Protein {
        return "Protein region";
    }
    match feature.display_direction() {
        Some(FeatureDirection::Nondirectional) => "Nondirectional",
        Some(FeatureDirection::Bidirectional) => "Bidirectional",
        Some(FeatureDirection::Forward) => "Forward strand",
        Some(FeatureDirection::Reverse) => "Reverse strand",
        _ if feature.reverse => "Reverse strand",
        _ => "Forward strand",
    }
}
pub(crate) fn wraps_origin(feature: &Feature, length: usize) -> bool {
    feature
        .parts
        .windows(2)
        .any(|parts| parts[0].end == length && parts[1].start == 0)
}
fn segment_window(len: usize, active: Option<usize>) -> std::ops::Range<usize> {
    let start = active
        .unwrap_or(0)
        .saturating_sub(1)
        .min(len.saturating_sub(3));
    start..(start + 3).min(len)
}
pub(crate) fn feature_color(feature: &Feature, index: usize) -> Hsla {
    if let Some(value) = feature.display_color_rgb() {
        return rgb(value).into();
    }
    rgb([0x4cc7aa, 0x65abe1, 0xe2ac64, 0xb985d0, 0x5bc7d2, 0xe38299][index % 6]).into()
}

pub(crate) fn render_feature_card(
    doc: &Document,
    feature: &Feature,
    segment: Option<usize>,
) -> AnyElement {
    let index = doc
        .features
        .iter()
        .position(|f| f.id == feature.id)
        .unwrap_or(0);
    let display_segments = feature.display_segments();
    let color = segment
        .and_then(|index| display_segments.get(index))
        .and_then(|segment| segment.color_rgb)
        .map(|value| rgb(value).into())
        .unwrap_or_else(|| feature_color(feature, index));
    let muted = rgb(0x96a6bc);
    let mut card = div()
        .w(px(342.))
        .p_4()
        .flex()
        .flex_col()
        .gap_2()
        .bg(rgb(0x202a3a))
        .border_1()
        .border_color(rgb(0x3c4b61))
        .rounded_lg()
        .shadow_lg()
        .text_color(rgb(0xe5edf5))
        .text_size(px(12.))
        .child(
            div()
                .font_weight(FontWeight::SEMIBOLD)
                .text_size(px(16.))
                .border_l_4()
                .border_color(color)
                .pl_2()
                .child(truncate(&feature.label, 120)),
        )
        .child(div().text_color(muted).child(format!(
            "{} · {} {} · {}",
            truncate(&feature.kind, 80),
            number(feature.len()),
            doc.unit(),
            direction(feature, doc.molecule)
        )));
    if wraps_origin(feature, doc.sequence.len()) {
        card = card.child(
            div()
                .text_color(rgb(0x5bd8bf))
                .child("Crosses the circular origin"),
        );
    }
    let mut ranges = div()
        .border_t_1()
        .border_color(rgb(0x374357))
        .pt_2()
        .flex()
        .flex_col()
        .gap_1();
    if display_segments.len() == 1 {
        let part = &display_segments[0].range;
        ranges = ranges.child(format!(
            "{}–{} · 1-based, inclusive",
            number(part.start + 1),
            number(part.end)
        ));
    } else {
        ranges = ranges.child(format!(
            "{} segments · 1-based, inclusive",
            display_segments.len()
        ));
        for index in segment_window(display_segments.len(), segment) {
            let display = &display_segments[index];
            let part = &display.range;
            ranges = ranges.child(
                div()
                    .text_color(if segment == Some(index) {
                        rgb(0xf19aa9)
                    } else {
                        rgb(0xdce6f2)
                    })
                    .child(format!(
                        "{}Segment {}{}: {}–{} · {} {}{}",
                        if segment == Some(index) { "• " } else { "" },
                        index + 1,
                        display
                            .name
                            .as_ref()
                            .map(|name| format!(" · {}", truncate(name, 60)))
                            .unwrap_or_default(),
                        number(part.start + 1),
                        number(part.end),
                        number(part.end - part.start),
                        doc.unit(),
                        display
                            .color_rgb
                            .map(|value| format!(" · #{value:06x}"))
                            .unwrap_or_default()
                    )),
            );
        }
        if display_segments.len() > 3 {
            ranges = ranges.child(
                div()
                    .text_color(muted)
                    .child("All segments are available in Features / Primers details."),
            );
        }
    }
    if display_segments.len() == 1
        && let Some(name) = &display_segments[0].name
    {
        ranges = ranges.child(format!("Segment: {}", truncate(name, 90)));
    }
    card = card.child(ranges);
    if feature.edited {
        card = card.child(
            div()
                .text_color(rgb(0xe2ac64))
                .child("Sequence changed — review this annotation."),
        );
    }
    if let Some(sequence) = feature.primer_sequence() {
        card = card
            .child(div().text_color(muted).child("Oligo sequence (5′ → 3′)"))
            .child(
                div()
                    .font_family("monospace")
                    .child(truncate(&sequence, 90)),
            )
            .child(
                div()
                    .text_color(muted)
                    .child("Binding annotation; thermodynamics not calculated."),
            );
    }
    let mut qualifiers = feature.qualifiers.iter().filter(|(key, _)| {
        !matches!(
            key.as_str(),
            "label" | "snapgene_color" | "snapgene_directionality" | "primer_sequence"
        )
    });
    for (key, value) in qualifiers.by_ref().take(4) {
        card = card
            .child(div().text_color(muted).child(truncate(
                if key == "translation" {
                    "Imported translation"
                } else {
                    key
                },
                60,
            )))
            .child(truncate(value.as_deref().unwrap_or("present"), 240));
    }
    if qualifiers.next().is_some() {
        card = card.child(
            div()
                .text_color(muted)
                .child("All qualifiers are available in Features / Primers details."),
        );
    }
    card.child(
        div()
            .border_t_1()
            .border_color(rgb(0x374357))
            .pt_2()
            .text_color(muted)
            .text_size(px(11.))
            .child("Click to select · Double-click to show sequence"),
    )
    .into_any_element()
}

/// Search stored annotation data without synthesizing an oligo from its binding site.
pub(crate) fn matches_query(feature: &Feature, query: &str) -> bool {
    let query = query.trim().to_lowercase();
    query.is_empty()
        || feature.label.to_lowercase().contains(&query)
        || feature.kind.to_lowercase().contains(&query)
        || feature
            .primer_sequence()
            .is_some_and(|sequence| sequence.to_lowercase().contains(&query))
        || feature.qualifiers.iter().any(|(key, value)| {
            key.to_lowercase().contains(&query)
                || value
                    .as_ref()
                    .is_some_and(|value| value.to_lowercase().contains(&query))
        })
}

/// Full, selectable text for persistent annotation inspection and clipboard export.
/// Unlike hover previews, it preserves every segment and qualifier, including
/// repeated qualifiers and valueless flags.
pub(crate) fn detail_text(feature: &Feature, molecule: Molecule) -> String {
    let mut text = format!(
        "{}\n{} · {} · {} {} in annotated segments\n",
        feature.label,
        feature.kind,
        direction(feature, molecule),
        number(feature.len()),
        molecule.unit()
    );
    if feature.edited {
        text.push_str("Sequence changed — review this annotation.\n");
    }
    text.push_str("\nSegments (1-based, inclusive; stored order)\n");
    for (index, part) in feature.parts.iter().enumerate() {
        let _ = writeln!(
            text,
            "{}: {}–{} · {} {}",
            index + 1,
            number(part.start + 1),
            number(part.end),
            number(part.end - part.start),
            molecule.unit()
        );
    }
    let display_segments = feature.display_segments();
    if display_segments
        .iter()
        .any(|segment| segment.name.is_some() || segment.color_rgb.is_some())
    {
        text.push_str("\nDisplay segments (imported names and colors; 1-based inclusive)\n");
        for (index, segment) in display_segments.iter().enumerate() {
            let _ = writeln!(
                text,
                "{}: {}–{} · {} · {}",
                index + 1,
                number(segment.range.start + 1),
                number(segment.range.end),
                segment.name.as_deref().unwrap_or("Unnamed"),
                segment
                    .color_rgb
                    .map(|value| format!("#{value:06x}"))
                    .unwrap_or_else(|| "Color not recorded".into())
            );
        }
    }
    if feature.kind == "primer_bind" {
        text.push_str("\nPrimer\n");
        if let Some(sequence) = feature.primer_sequence() {
            let _ = writeln!(text, "Recorded oligo (5′ → 3′): {sequence}");
            let _ = writeln!(text, "Oligo length: {} nt", number(sequence.len()));
        } else {
            text.push_str("Oligo sequence: not recorded. Binding length does not establish the full oligo length.\n");
        }
        let _ = writeln!(
            text,
            "Binding length: {} {}",
            number(feature.len()),
            molecule.unit()
        );
        let _ = writeln!(
            text,
            "Binding strand: {}",
            if feature.reverse {
                "Reverse"
            } else {
                "Forward"
            }
        );
        if let Some(tm) = feature.primer_melting_temperature() {
            let _ = writeln!(
                text,
                "Imported Tm: {tm} °C (source value; not recalculated)"
            );
        } else {
            text.push_str("Tm: not recorded; thermodynamics not calculated.\n");
        }
        if let Some(date) = feature.primer_added_date() {
            let _ = writeln!(text, "Imported added date: {date}");
        }
    }
    let _ = writeln!(
        text,
        "\nAll stored qualifiers ({})",
        feature.qualifiers.len()
    );
    for (key, value) in &feature.qualifiers {
        match value {
            Some(value) => {
                let _ = writeln!(text, "/{key}={value}");
            }
            None => {
                let _ = writeln!(text, "/{key} (present)");
            }
        }
    }
    text
}

#[cfg(test)]
mod tests {
    use super::{detail_text, direction, matches_query, segment_window, truncate, wraps_origin};
    use popgene_core::{Molecule, Segment};
    fn sample() -> popgene_core::Document {
        popgene_formats::import(
            "annotated.gb",
            include_bytes!("../../../validation/fixtures/annotated.gb"),
        )
        .unwrap()
    }
    #[test]
    fn previews_preserve_unicode_and_bound_large_qualifiers() {
        assert_eq!(truncate("αβγδε", 3), "αβγ…");
        assert_eq!(truncate("αβ", 3), "αβ");
        assert_eq!(truncate("", 0), "");
        assert_eq!(truncate("A", 0), "…");
        assert_eq!(truncate(&"A".repeat(1_000_000), 4), "AAAA…");
    }
    #[test]
    fn compound_card_preserves_segments_and_origin() {
        let document = sample();
        assert!(wraps_origin(&document.features[0], 5000));
        assert_eq!(document.features[0].len(), 220);
        assert_eq!(segment_window(8, Some(6)), 5..8);
        assert_eq!(segment_window(8, None), 0..3);
        assert_eq!(segment_window(2, Some(1)), 0..2);
    }
    #[test]
    fn direction_honors_nondirectional_and_reverse_features() {
        let mut feature = sample().features[2].clone();
        assert_eq!(direction(&feature, Molecule::Dna), "Reverse strand");
        feature
            .qualifiers
            .push(("snapgene_directionality".into(), Some("0".into())));
        assert_eq!(direction(&feature, Molecule::Dna), "Nondirectional");
        assert_eq!(direction(&feature, Molecule::Protein), "Protein region");
    }
    #[test]
    fn persistent_details_preserve_all_segments_qualifiers_and_text() {
        let mut feature = sample().features[0].clone();
        feature.parts = (0..8)
            .map(|i| Segment {
                start: i * 10,
                end: i * 10 + 4,
            })
            .collect();
        let full_text = "Complete qualifier αβγ. ".repeat(30);
        feature.qualifiers = (0..6)
            .map(|i| (format!("field{i}"), Some(format!("value{i}"))))
            .collect();
        feature
            .qualifiers
            .push(("note".into(), Some(full_text.clone())));
        feature
            .qualifiers
            .push(("note".into(), Some("second note".into())));
        feature.qualifiers.push(("pseudo".into(), None));
        let text = detail_text(&feature, Molecule::Dna);
        assert!(text.contains("8: 71–74 · 4 nt"));
        assert!(text.contains("/field5=value5"));
        assert!(text.contains(&format!("/note={full_text}")));
        assert!(text.contains("/note=second note"));
        assert!(text.contains("/pseudo (present)"));
        assert!(text.contains("All stored qualifiers (9)"));
    }
    #[test]
    fn primer_details_distinguish_recorded_oligo_and_binding_length() {
        let mut primer = sample().features[5].clone();
        primer.qualifiers = vec![("primer_sequence".into(), Some("aaACGT".into()))];
        let text = detail_text(&primer, Molecule::Dna);
        assert!(text.contains("Recorded oligo (5′ → 3′): aaACGT"));
        assert!(text.contains("Oligo length: 6 nt"));
        assert!(text.contains("Binding length: 24 nt"));
        assert!(text.contains("Binding strand: Reverse"));
        assert!(text.contains("Tm: not recorded"));
        primer
            .qualifiers
            .push(("binding_meltingTemperature".into(), Some("59.3".into())));
        assert!(
            detail_text(&primer, Molecule::Dna)
                .contains("Imported Tm: 59.3 °C (source value; not recalculated)")
        );
        primer.edited = true;
        assert!(!detail_text(&primer, Molecule::Dna).contains("Imported Tm:"));
        primer.qualifiers.clear();
        assert!(detail_text(&primer, Molecule::Dna).contains("Oligo sequence: not recorded"));
    }
    #[test]
    fn annotation_search_finds_sequence_and_qualifiers_case_insensitively() {
        let document = sample();
        let primer = &document.features[4];
        assert!(matches_query(primer, "cttaggatgatcacgttacagttc"));
        assert!(matches_query(primer, "AUDIT_FORWARD"));
        assert!(matches_query(primer, "2026-10-21"));
        assert!(!matches_query(primer, "missing marker"));
        assert!(matches_query(primer, "  "));
    }
    #[test]
    fn details_distinguish_biological_parts_from_named_display_subdivisions() {
        let doc = popgene_formats::import(
            "segments.gb",
            include_bytes!("../../../validation/fixtures/regression/segment-colours.gb"),
        )
        .unwrap();
        let text = detail_text(&doc.features[0], doc.molecule);
        assert!(text.contains("1: 101–400 · 300 nt"));
        assert!(text.contains("1: 101–250 · Red domain · #ff0000"));
        assert!(text.contains("2: 251–400 · Blue domain · #0000ff"));
        assert!(text.contains("Original explanatory note; retained separately."));
    }
}
