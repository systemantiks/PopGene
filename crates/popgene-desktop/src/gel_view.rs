//! Simulated agarose gel drawing, vector export, plain-text report, shared gel
//! controls and the saved gel defaults.
use gpui_kit::component::{
    Disableable, Selectable, Sizable, StyledExt,
    button::Button,
    input::{Input, InputEvent, InputState},
    tooltip::Tooltip,
};
use gpui_kit::{prelude::*, *};
use popgene_core::gel::{
    self, Band, BandEdge, Buffer, Gel, GelSettings, LADDERS, MODEL_NOTE, RUN_LENGTH, format_bp,
};
use serde::{Deserialize, Serialize};
use std::{fmt::Write as _, rc::Rc};

const AXIS: f32 = 58.;
const LANE: f32 = 46.;
const LABELS: f32 = 30.;
const WELL: f32 = 10.;
const RUN: f32 = 300.;
const BAND: f32 = 3.;
const REFERENCE_BAND: f32 = 5.;
/// Called with (lane, band) or (lane, None) when a lane label is clicked.
pub type GelClick = Rc<dyn Fn(usize, Option<usize>, &mut Window, &mut App)>;

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
pub fn band_description(band: &Band) -> String {
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
/// Top edge and height of a drawn band in gel pixels.
fn band_box(band: &Band, edge: BandEdge) -> (f32, f32) {
    let thickness = if band.reference { REFERENCE_BAND } else { BAND };
    let top = band.top as f32 * RUN
        - match edge {
            BandEdge::Top => 0.,
            BandEdge::Middle => thickness / 2.,
            BandEdge::Bottom => thickness,
        };
    (top, (band.bottom - band.top) as f32 * RUN + thickness)
}
fn alpha(band: &Band) -> f32 {
    (0.3 + 0.7 * band.intensity.sqrt() as f32).min(1.)
}
/// Size labels follow the first ladder lane, skipping labels that would overlap.
fn axis_labels(gel: &Gel) -> Vec<(f32, String)> {
    let mut labels = vec![];
    let mut previous = f32::NEG_INFINITY;
    if let Some(ladder) = gel.lanes.iter().find(|l| l.ladder) {
        for band in &ladder.bands {
            let y = band.top as f32 * RUN;
            if y - previous >= 12. {
                previous = y;
                labels.push((y, size_label(band.lengths[0])));
            }
        }
    }
    labels
}

pub fn gel_figure(
    gel: &Gel,
    selected_lane: Option<usize>,
    highlighted: &dyn Fn(usize, usize) -> bool,
    on_click: Option<GelClick>,
) -> Div {
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
    for (y, label) in axis_labels(gel) {
        figure = figure.child(
            div()
                .absolute()
                .left(px(0.))
                .top(px(top + y - 7.))
                .w(px(AXIS - 8.))
                .text_right()
                .text_size(px(10.))
                .text_color(rgb(0x95a6ba))
                .child(label),
        );
    }
    for (i, lane) in gel.lanes.iter().enumerate() {
        let x = AXIS + i as f32 * LANE;
        let label = lane.label.clone();
        let click = on_click.clone();
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
                    .text_color(rgb(if selected_lane == Some(i) {
                        0x48c9bb
                    } else {
                        0xc6d3e0
                    }))
                    .truncate()
                    .child(lane.label.clone())
                    .tooltip(move |window, cx| Tooltip::new(label.clone()).build(window, cx))
                    .when_some(click, |d, click| {
                        d.cursor_pointer()
                            .on_click(move |_, window, cx| click(i, None, window, cx))
                    }),
            )
            .child(
                div()
                    .absolute()
                    .left(px(x + 7.))
                    .top(px(LABELS + 3.))
                    .w(px(LANE - 14.))
                    .h(px(WELL - 4.))
                    .border_1()
                    .border_color(rgb(if selected_lane == Some(i) {
                        0x48c9bb
                    } else {
                        0x3a4a5e
                    })),
            );
        for (j, band) in lane.bands.iter().enumerate() {
            let (y, height) = band_box(band, gel.settings.band_edge);
            let opacity = (alpha(band) * 255.).round() as u32;
            let description = band_description(band);
            let click = on_click.clone();
            let drawn = div()
                .id(("gel-band", i * 10_000 + j))
                .absolute()
                .left(px(x + 6.))
                .top(px(top + y))
                .w(px(LANE - 12.))
                .h(px(height))
                .rounded_sm()
                .tooltip(move |window, cx| Tooltip::new(description.clone()).build(window, cx))
                .when_some(click, |d, click| {
                    d.cursor_pointer()
                        .on_click(move |_, window, cx| click(i, Some(j), window, cx))
                });
            let drawn = if band.circular {
                drawn.border_1().border_dashed().border_color(rgb(0xe5b879))
            } else {
                drawn.bg(rgba(0xeaf4ff00 | opacity))
            };
            figure = figure.child(if highlighted(i, j) {
                drawn.border_1().border_color(rgb(0x48c9bb))
            } else {
                drawn
            });
        }
    }
    figure
}

fn escape(value: &str) -> String {
    value.chars().fold(String::new(), |mut output, c| {
        match c {
            '&' => output.push_str("&amp;"),
            '<' => output.push_str("&lt;"),
            '>' => output.push_str("&gt;"),
            '"' => output.push_str("&quot;"),
            '\'' => output.push_str("&apos;"),
            c if c.is_control() => output.push('\u{fffd}'),
            c => output.push(c),
        }
        output
    })
}

/// Standalone figure using the on-screen geometry. Inverted draws dark bands on
/// a white gel, which prints well; otherwise light bands on a dark gel.
pub fn gel_svg(gel: &Gel, title: &str, inverted: bool) -> String {
    const MARGIN: f32 = 16.;
    const LANE_LABEL: f32 = 72.;
    let slab = gel.lanes.len() as f32 * LANE;
    let (width, top) = (MARGIN * 2. + AXIS + slab, MARGIN + LANE_LABEL + WELL);
    let height = top + RUN + 12. + MARGIN + 30.;
    let (page, slab_fill, band_fill, ink, muted) = if inverted {
        ("#ffffff", "#ffffff", "#111111", "#111111", "#555555")
    } else {
        ("#101722", "#05080c", "#eaf4ff", "#e4edf4", "#95a6ba")
    };
    let left = MARGIN + AXIS;
    let mut svg = format!(
        "<svg xmlns='http://www.w3.org/2000/svg' width='{width:.0}' height='{height:.0}' viewBox='0 0 {width:.0} {height:.0}'><title>{}</title><desc>{}</desc><rect width='100%' height='100%' fill='{page}'/><g font-family='Noto Sans, sans-serif'>",
        escape(title),
        escape(MODEL_NOTE)
    );
    write!(
        svg,
        "<rect x='{left:.1}' y='{:.1}' width='{slab:.1}' height='{:.1}' rx='4' fill='{slab_fill}' stroke='{muted}' stroke-width='0.8'/>",
        top - WELL,
        WELL + RUN + 12.
    )
    .unwrap();
    for (y, label) in axis_labels(gel) {
        write!(
            svg,
            "<text x='{:.1}' y='{:.1}' font-size='10' text-anchor='end' fill='{muted}'>{}</text>",
            left - 8.,
            top + y + 3.5,
            escape(&label)
        )
        .unwrap();
    }
    for (i, lane) in gel.lanes.iter().enumerate() {
        let x = left + i as f32 * LANE;
        let label: String = lane.label.chars().take(28).collect();
        write!(
            svg,
            "<text transform='translate({:.1} {:.1}) rotate(-50)' font-size='10' fill='{ink}'>{}</text><rect x='{:.1}' y='{:.1}' width='{:.1}' height='{:.1}' fill='none' stroke='{muted}' stroke-width='0.8'/>",
            x + LANE / 2.,
            top - WELL - 4.,
            escape(&label),
            x + 7.,
            top - WELL + 3.,
            LANE - 14.,
            WELL - 4.
        )
        .unwrap();
        for band in &lane.bands {
            let (y, h) = band_box(band, gel.settings.band_edge);
            if band.circular {
                write!(
                    svg,
                    "<rect x='{:.1}' y='{:.2}' width='{:.1}' height='{h:.2}' rx='1' fill='none' stroke='#c8962e' stroke-width='1' stroke-dasharray='3 2'/>",
                    x + 6.,
                    top + y,
                    LANE - 12.
                )
                .unwrap();
            } else {
                write!(
                    svg,
                    "<rect x='{:.1}' y='{:.2}' width='{:.1}' height='{h:.2}' rx='1' fill='{band_fill}' fill-opacity='{:.3}'/>",
                    x + 6.,
                    top + y,
                    LANE - 12.,
                    alpha(band)
                )
                .unwrap();
            }
        }
    }
    write!(
        svg,
        "<text x='{MARGIN:.1}' y='{:.1}' font-size='10' fill='{muted}'>{}% agarose · {} · run length {:.1}× · simulated, not measured</text></g></svg>",
        height - MARGIN,
        gel.settings.agarose_percent,
        gel.settings.buffer.label(),
        gel.settings.run_length
    )
    .unwrap();
    svg
}

pub fn report(gel: &Gel) -> String {
    let mut text = format!(
        "Simulated agarose gel: {}% standard agarose in {}; run length {:.1}×; resolution {}–{} bp\n{MODEL_NOTE}\n\nLane\tLabel\tBand\tGel position (%)\tRelative intensity (%)\n",
        gel.settings.agarose_percent,
        gel.settings.buffer.label(),
        gel.settings.run_length,
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
        for length in &lane.ran_off {
            text.push_str(&format!(
                "{}\t{}\t{} bp\tran off the gel\t\n",
                i + 1,
                lane.label,
                format_bp(*length)
            ));
        }
    }
    let noncutters = gel
        .lanes
        .iter()
        .filter(|l| !l.noncutters.is_empty())
        .map(|l| format!("{}: {}\n", l.label, l.noncutters.join(", ")))
        .collect::<String>();
    if !noncutters.is_empty() {
        text.push_str(&format!("\nEnzymes without a site\n{noncutters}"));
    }
    if !gel.warnings.is_empty() {
        text.push_str("\nNotes\n");
        for warning in &gel.warnings {
            text.push_str(&format!("{warning}\n"));
        }
    }
    text
}

/// Ladder and gel settings that new gels start from.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct GelDefaults {
    pub ladder: String,
    pub settings: GelSettings,
}
impl Default for GelDefaults {
    fn default() -> Self {
        Self {
            ladder: LADDERS[0].name.into(),
            settings: GelSettings::default(),
        }
    }
}
impl GelDefaults {
    fn valid(self) -> Self {
        if gel::ladder(&self.ladder).is_some() && self.settings.validate().is_ok() {
            self
        } else {
            Self::default()
        }
    }
    /// Unreadable or outdated defaults fall back to the built-in ones.
    pub fn load() -> Self {
        defaults_path()
            .and_then(|path| std::fs::read(path).ok())
            .filter(|bytes| bytes.len() < 64 * 1024)
            .and_then(|bytes| serde_json::from_slice::<Self>(&bytes).ok())
            .map_or_else(Self::default, Self::valid)
    }
    pub fn save(&self) -> Result<(), String> {
        let path = defaults_path().ok_or("No application state folder is available")?;
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
        }
        let bytes = serde_json::to_vec_pretty(self).map_err(|e| e.to_string())?;
        popgene_formats::atomic_export(&path, &bytes).map_err(|e| e.to_string())
    }
}
/// Kept beside the recovery database, so POPGENE_SESSION_PATH isolates it too.
fn defaults_path() -> Option<std::path::PathBuf> {
    Some(
        crate::session::default_path()?
            .parent()?
            .join("gel-defaults.json"),
    )
}

pub enum GelControlsEvent {
    Changed,
}
/// Agarose, buffer, run length and band-edge controls shared by gel views.
pub struct GelControls {
    agarose: Entity<InputState>,
    buffer: Buffer,
    run_length: f64,
    band_edge: BandEdge,
    _subscription: Subscription,
}
impl EventEmitter<GelControlsEvent> for GelControls {}
impl GelControls {
    pub fn new(settings: GelSettings, window: &mut Window, cx: &mut Context<Self>) -> Self {
        let agarose = cx.new(|cx| {
            InputState::new(window, cx).default_value(settings.agarose_percent.to_string())
        });
        let subscription = cx.subscribe(&agarose, |_, _, event: &InputEvent, cx| {
            if matches!(event, InputEvent::Change) {
                cx.emit(GelControlsEvent::Changed);
            }
        });
        Self {
            agarose,
            buffer: settings.buffer,
            run_length: settings.run_length,
            band_edge: settings.band_edge,
            _subscription: subscription,
        }
    }
    pub fn settings(&self, cx: &App) -> Result<GelSettings, String> {
        let settings = GelSettings {
            agarose_percent: self
                .agarose
                .read(cx)
                .value()
                .trim()
                .parse::<f64>()
                .map_err(|_| "Enter a number for agarose (%)".to_string())?,
            buffer: self.buffer,
            run_length: self.run_length,
            band_edge: self.band_edge,
        };
        settings.validate()?;
        Ok(settings)
    }
    pub fn set(&mut self, settings: GelSettings, window: &mut Window, cx: &mut Context<Self>) {
        self.agarose.update(cx, |input, cx| {
            input.set_value(settings.agarose_percent.to_string(), window, cx)
        });
        self.buffer = settings.buffer;
        self.run_length = settings.run_length;
        self.band_edge = settings.band_edge;
        cx.emit(GelControlsEvent::Changed);
        cx.notify();
    }
    fn changed(&mut self, cx: &mut Context<Self>) {
        cx.emit(GelControlsEvent::Changed);
        cx.notify();
    }
}
impl Render for GelControls {
    fn render(&mut self, _: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        let field =
            |label: &'static str| div().v_flex().gap_1().child(div().text_xs().child(label));
        let (low, high) = self.buffer.agarose_range();
        div()
            .h_flex()
            .flex_wrap()
            .items_start()
            .gap_4()
            .child(
                field("Agarose (%)")
                    .child(Input::new(&self.agarose).w(px(90.)))
                    .child(
                        div()
                            .text_xs()
                            .text_color(rgb(0x96a6bc))
                            .child(format!("{low}–{high}% in {}", self.buffer.label())),
                    ),
            )
            .child(field("Buffer").child(div().h_flex().gap_1().children(
                [(0usize, Buffer::Tae), (1, Buffer::Tbe)].map(|(i, buffer)| {
                    Button::new(("gel-buffer", i))
                        .small()
                        .selected(self.buffer == buffer)
                        .label(buffer.label())
                        .on_click(cx.listener(move |this, _, _, cx| {
                            this.buffer = buffer;
                            this.changed(cx);
                        }))
                }),
            )))
            .child(
                field("Run length").child(
                    div()
                        .h_flex()
                        .gap_1()
                        .items_center()
                        .child(
                            Button::new("gel-run-shorter")
                                .small()
                                .label("−")
                                .disabled(self.run_length <= *RUN_LENGTH.start() + 1e-9)
                                .on_click(cx.listener(|this, _, _, cx| {
                                    this.run_length = ((this.run_length - 0.1) * 10.).round() / 10.;
                                    this.run_length = this.run_length.max(*RUN_LENGTH.start());
                                    this.changed(cx);
                                })),
                        )
                        .child(
                            div()
                                .w(px(40.))
                                .text_center()
                                .text_sm()
                                .child(format!("{:.1}×", self.run_length)),
                        )
                        .child(
                            Button::new("gel-run-longer")
                                .small()
                                .label("+")
                                .disabled(self.run_length >= *RUN_LENGTH.end() - 1e-9)
                                .on_click(cx.listener(|this, _, _, cx| {
                                    this.run_length = ((this.run_length + 0.1) * 10.).round() / 10.;
                                    this.run_length = this.run_length.min(*RUN_LENGTH.end());
                                    this.changed(cx);
                                })),
                        ),
                ),
            )
            .child(
                field("Position marks band").child(
                    div().h_flex().gap_1().children(
                        [
                            (0usize, BandEdge::Top, "Top"),
                            (1, BandEdge::Middle, "Middle"),
                            (2, BandEdge::Bottom, "Bottom"),
                        ]
                        .map(|(i, edge, label)| {
                            Button::new(("gel-band-edge", i))
                                .small()
                                .selected(self.band_edge == edge)
                                .label(label)
                                .on_click(cx.listener(move |this, _, _, cx| {
                                    this.band_edge = edge;
                                    this.changed(cx);
                                }))
                        }),
                    ),
                ),
            )
    }
}

/// Ladder buttons, with an optional custom entry, for choosing one ladder.
pub fn ladder_buttons(
    selected: Option<usize>,
    custom: bool,
    on_select: impl Fn(Option<usize>, &mut Window, &mut App) + 'static,
) -> Div {
    let on_select = Rc::new(on_select);
    let mut row = div().h_flex().flex_wrap().gap_1();
    for (i, ladder) in LADDERS.iter().enumerate() {
        let on_select = on_select.clone();
        row = row.child(
            Button::new(("gel-ladder", i))
                .small()
                .selected(selected == Some(i))
                .label(ladder.name)
                .on_click(move |_, window, cx| on_select(Some(i), window, cx)),
        );
    }
    if custom {
        row = row.child(
            Button::new("gel-custom-ladder")
                .small()
                .selected(selected.is_none())
                .label("Custom")
                .on_click(move |_, window, cx| on_select(None, window, cx)),
        );
    }
    row
}

pub fn parse_ladder(text: &str) -> Result<Vec<usize>, String> {
    let mut sizes = vec![];
    for size in text
        .split(|c: char| c == ',' || c == ';' || c.is_whitespace())
        .filter(|s| !s.is_empty())
    {
        let size = size
            .parse::<usize>()
            .ok()
            .filter(|s| (1..=popgene_core::MAX_BASES).contains(s))
            .ok_or_else(|| format!("Ladder size {size} must be a whole number of base pairs from 1 to {}, without thousands separators", popgene_core::MAX_BASES))?;
        if !sizes.contains(&size) {
            sizes.push(size);
        }
    }
    if sizes.is_empty() || sizes.len() > gel::MAX_LADDER_SIZES {
        return Err(format!(
            "Enter between 1 and {} ladder sizes",
            gel::MAX_LADDER_SIZES
        ));
    }
    sizes.sort_unstable_by(|a, b| b.cmp(a));
    Ok(sizes)
}

#[cfg(test)]
mod tests {
    use super::{GelDefaults, band_box, band_lengths, gel_svg, parse_ladder, report, size_label};
    use popgene_core::gel::{
        BandEdge, Buffer, GelSettings, LADDERS, Lane, LaneFragment, MODEL_NOTE, simulate,
    };

    fn lane(lengths: &[usize], circular: bool) -> Lane {
        Lane {
            label: "EcoRI <cut> & BamHI".into(),
            fragments: lengths
                .iter()
                .map(|&length| LaneFragment {
                    length,
                    circular,
                    reference: false,
                })
                .collect(),
            equimolar: true,
            noncutters: vec!["NotI".into()],
            ..Default::default()
        }
    }

    #[test]
    fn axis_and_band_labels_use_readable_units() {
        assert_eq!(size_label(10_000), "10 kb");
        assert_eq!(size_label(1_500), "1.5 kb");
        assert_eq!(size_label(750), "750 bp");
        assert_eq!(size_label(2_322), "2,322 bp");
        let gel = simulate(
            &[lane(&[1_010, 1_000, 500, 500], false)],
            &GelSettings::default(),
        )
        .unwrap();
        let labels = gel.lanes[0]
            .bands
            .iter()
            .map(band_lengths)
            .collect::<Vec<_>>();
        assert_eq!(labels, ["1,000–1,010 bp (2 fragments)", "2 × 500 bp"]);
    }

    #[test]
    fn band_edge_moves_the_drawn_band_around_the_calculated_position() {
        let gel = simulate(&[lane(&[1_000], false)], &GelSettings::default()).unwrap();
        let band = &gel.lanes[0].bands[0];
        let y = band.top as f32 * 300.;
        assert_eq!(band_box(band, BandEdge::Top).0, y);
        assert_eq!(band_box(band, BandEdge::Middle).0, y - 1.5);
        assert_eq!(band_box(band, BandEdge::Bottom).0, y - 3.);
    }

    #[test]
    fn report_lists_bands_noncutters_runoff_and_assumptions() {
        let settings = GelSettings {
            agarose_percent: 0.7,
            buffer: Buffer::Tbe,
            run_length: 2.,
            ..Default::default()
        };
        let gel = simulate(
            &[Lane::ladder(&LADDERS[4]), lane(&[3_000, 150], false)],
            &settings,
        )
        .unwrap();
        let text = report(&gel);
        assert!(text.contains("0.7% standard agarose in 1X TBE; run length 2.0×"));
        assert!(text.contains("not measured mobilities"));
        assert!(text.contains("150 bp\tran off the gel"));
        assert!(text.contains("Enzymes without a site\nEcoRI <cut> & BamHI: NotI"));
        assert!(text.contains("\nNotes\n"));
    }

    #[test]
    fn svg_escapes_labels_and_marks_circular_dna() {
        let gel = simulate(
            &[Lane::ladder(&LADDERS[0]), lane(&[3_000], true)],
            &GelSettings::default(),
        )
        .unwrap();
        let svg = gel_svg(&gel, "Gel <1>", true);
        assert!(svg.starts_with("<svg xmlns='http://www.w3.org/2000/svg'"));
        assert!(svg.contains("<title>Gel &lt;1&gt;</title>"));
        assert!(svg.contains("EcoRI &lt;cut&gt; &amp; BamHI"));
        assert!(svg.contains("stroke-dasharray"));
        assert!(svg.contains("10 kb"));
        assert!(svg.contains(&MODEL_NOTE[..40]));
        let pdf = crate::figure::encode(&svg, crate::figure::FigureFormat::Pdf, 1.).unwrap();
        assert!(pdf.starts_with(b"%PDF"));
    }

    #[test]
    fn custom_ladder_sizes_are_whole_unique_and_largest_first() {
        assert_eq!(
            parse_ladder("500, 10000 1000;500").unwrap(),
            [10_000, 1_000, 500]
        );
        assert!(parse_ladder("10,000").is_err());
        assert!(parse_ladder("1.5").is_err());
        assert!(parse_ladder("").is_err());
        assert!(parse_ladder(&"1 ".repeat(65)).is_ok());
        let many = (1..=65)
            .map(|i| i.to_string())
            .collect::<Vec<_>>()
            .join(" ");
        assert!(parse_ladder(&many).is_err());
    }

    #[test]
    fn invalid_defaults_fall_back_to_built_in_ones() {
        let unknown = GelDefaults {
            ladder: "Unknown".into(),
            ..Default::default()
        };
        assert_eq!(unknown.valid(), GelDefaults::default());
        let thin = GelDefaults {
            settings: GelSettings {
                agarose_percent: 9.,
                ..Default::default()
            },
            ..Default::default()
        };
        assert_eq!(thin.valid(), GelDefaults::default());
        let kept = GelDefaults {
            ladder: LADDERS[2].name.into(),
            settings: GelSettings {
                agarose_percent: 2.,
                ..Default::default()
            },
        };
        assert_eq!(kept.clone().valid(), kept);
    }
}
