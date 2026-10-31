//! Simulated agarose gel drawing and its plain-text band report.
use gpui_kit::component::tooltip::Tooltip;
use gpui_kit::*;
use popgene_core::gel::{Band, Gel, MODEL_NOTE, format_bp};

const AXIS: f32 = 58.;
const LANE: f32 = 46.;
const LABELS: f32 = 30.;
const WELL: f32 = 10.;
const RUN: f32 = 300.;
const BAND: f32 = 3.;
const REFERENCE_BAND: f32 = 5.;

fn size_label(length: usize) -> String {
    if length >= 1_000 && length.is_multiple_of(100) {
        format!("{} kb", length as f64 / 1_000.)
    } else {
        format!("{} bp", format_bp(length))
    }
}
fn band_lengths(band: &Band) -> String {
    let (largest, smallest) = (band.lengths[0], band.lengths[band.lengths.len() - 1]);
    if band.lengths.len() == 1 {
        format!("{} bp", format_bp(largest))
    } else if largest == smallest {
        format!("{} × {} bp", band.lengths.len(), format_bp(largest))
    } else {
        format!(
            "{}–{} bp ({} fragments)",
            format_bp(smallest),
            format_bp(largest),
            band.lengths.len()
        )
    }
}
fn band_description(band: &Band) -> String {
    let mut text = band_lengths(band);
    text.push_str(&format!(
        " · {:.0}% of the lane's strongest band",
        band.intensity * 100.
    ));
    if band.reference {
        text.push_str(" · reference band");
    }
    if band.circular {
        text.push_str(" · circular, drawn at its linear length");
    }
    text
}

pub fn gel_figure(gel: &Gel) -> impl IntoElement {
    let slab = gel.lanes.len() as f32 * LANE;
    let top = LABELS + WELL;
    let mut figure = div()
        .relative()
        .flex_none()
        .w(px(AXIS + slab))
        .h(px(top + RUN + 12.))
        .child(
            div()
                .absolute()
                .left(px(AXIS))
                .top(px(LABELS))
                .w(px(slab))
                .h(px(WELL + RUN + 12.))
                .rounded_md()
                .bg(rgb(0x05080c))
                .border_1()
                .border_color(rgb(0x2a3646)),
        );
    if let Some(ladder) = gel.lanes.first() {
        let mut previous = f32::NEG_INFINITY;
        for band in &ladder.bands {
            let y = top + band.top as f32 * RUN;
            if y - previous < 12. {
                continue;
            }
            previous = y;
            figure = figure.child(
                div()
                    .absolute()
                    .left(px(0.))
                    .top(px(y - 7.))
                    .w(px(AXIS - 8.))
                    .text_right()
                    .text_size(px(10.))
                    .text_color(rgb(0x95a6ba))
                    .child(size_label(band.lengths[0])),
            );
        }
    }
    for (i, lane) in gel.lanes.iter().enumerate() {
        let x = AXIS + i as f32 * LANE;
        let label = lane.label.clone();
        figure = figure
            .child(
                div()
                    .id(("gel-lane", i))
                    .absolute()
                    .left(px(x + 2.))
                    .top(px(4.))
                    .w(px(LANE - 4.))
                    .text_center()
                    .text_size(px(10.))
                    .text_color(rgb(0xc6d3e0))
                    .truncate()
                    .child(lane.label.clone())
                    .tooltip(move |window, cx| Tooltip::new(label.clone()).build(window, cx)),
            )
            .child(
                div()
                    .absolute()
                    .left(px(x + 7.))
                    .top(px(LABELS + 3.))
                    .w(px(LANE - 14.))
                    .h(px(WELL - 4.))
                    .border_1()
                    .border_color(rgb(0x3a4a5e)),
            );
        for (j, band) in lane.bands.iter().enumerate() {
            let thickness = if band.reference { REFERENCE_BAND } else { BAND };
            let alpha = ((0.3 + 0.7 * band.intensity.sqrt()).min(1.) * 255.).round() as u32;
            let description = band_description(band);
            let drawn = div()
                .id(("gel-band", i * 10_000 + j))
                .absolute()
                .left(px(x + 6.))
                .top(px(top + band.top as f32 * RUN - thickness / 2.))
                .w(px(LANE - 12.))
                .h(px((band.bottom - band.top) as f32 * RUN + thickness))
                .rounded_sm()
                .tooltip(move |window, cx| Tooltip::new(description.clone()).build(window, cx));
            figure = figure.child(if band.circular {
                drawn.border_1().border_dashed().border_color(rgb(0xe5b879))
            } else {
                drawn.bg(rgba(0xeaf4ff00 | alpha))
            });
        }
    }
    figure
}

pub fn report(gel: &Gel) -> String {
    let mut text = format!(
        "Simulated agarose gel: {}% standard agarose; efficient separation {}–{} bp\n{MODEL_NOTE}\n\nLane\tLabel\tBand\tRun position (%)\tRelative intensity (%)\n",
        gel.agarose_percent,
        format_bp(gel.separation.0),
        format_bp(gel.separation.1)
    );
    for (i, lane) in gel.lanes.iter().enumerate() {
        for band in &lane.bands {
            text.push_str(&format!(
                "{}\t{}\t{}{}\t{:.1}\t{:.0}\n",
                i + 1,
                lane.label,
                band_lengths(band),
                if band.circular { " (circular)" } else { "" },
                (band.top + band.bottom) * 50.,
                band.intensity * 100.
            ));
        }
    }
    if !gel.warnings.is_empty() {
        text.push_str("\nNotes\n");
        for warning in &gel.warnings {
            text.push_str(&format!("{warning}\n"));
        }
    }
    text
}

#[cfg(test)]
mod tests {
    use super::{band_lengths, report, size_label};
    use popgene_core::gel::{LADDERS, Lane, LaneFragment, simulate};

    #[test]
    fn axis_and_band_labels_use_readable_units() {
        assert_eq!(size_label(10_000), "10 kb");
        assert_eq!(size_label(1_500), "1.5 kb");
        assert_eq!(size_label(750), "750 bp");
        assert_eq!(size_label(2_322), "2,322 bp");
        let lane = Lane {
            label: "EcoRI".into(),
            fragments: [1_010, 1_000, 500, 500]
                .into_iter()
                .map(|length| LaneFragment {
                    length,
                    circular: false,
                    reference: false,
                })
                .collect(),
            equimolar: true,
        };
        let gel = simulate(&[lane], 1.0).unwrap();
        let labels = gel.lanes[0]
            .bands
            .iter()
            .map(band_lengths)
            .collect::<Vec<_>>();
        assert_eq!(labels, ["1,000–1,010 bp (2 fragments)", "2 × 500 bp"]);
    }

    #[test]
    fn report_lists_every_lane_band_and_assumption() {
        let gel = simulate(&[Lane::ladder(&LADDERS[1])], 0.7).unwrap();
        let text = report(&gel);
        assert!(text.contains("0.7% standard agarose; efficient separation 800–12,000 bp"));
        assert!(text.contains("not measured mobilities"));
        assert_eq!(text.matches("\tλ DNA/HindIII\t").count(), 8);
        assert!(text.contains("\nNotes\n"));
    }
}
