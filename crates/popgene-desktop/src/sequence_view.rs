//! Virtualized annotated sequence and Sanger trace views, rendered directly by GPUI.
use crate::{
    feature_card,
    sequence_layout::{self, SequenceLayout},
};
use gpui_kit::component::scroll::Scrollbar;
use gpui_kit::component::{
    ActiveTheme, Sizable, StyledExt, VirtualListScrollHandle, button::*, h_virtual_list,
    hover_card::HoverCard, v_virtual_list,
};
use gpui_kit::{prelude::*, *};
use popgene_core::{Document, Molecule, Segment, Trace, translation::frame_codons};
use std::{cell::RefCell, rc::Rc, sync::Arc};

const LEFT: f32 = 58.;
const RIGHT: f32 = 66.;
const INK: u32 = 0xdce5ed;
const MUTED: u32 = 0x8799ab;
const ACCENT: u32 = 0x48c9bb;
const MAX_LANES: usize = 4;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SequenceOptions {
    pub show_complement: bool,
    pub show_features: bool,
    pub show_primers: bool,
    pub show_enzymes: bool,
    pub show_translations: bool,
    /// Bits 0..2: forward +1/+2/+3; bits 3..5: reverse -1/-2/-3.
    pub reading_frames: u8,
    pub colored_bases: bool,
    pub compact: bool,
}
impl Default for SequenceOptions {
    fn default() -> Self {
        Self {
            show_complement: true,
            show_features: true,
            show_primers: true,
            show_enzymes: true,
            show_translations: true,
            reading_frames: 0,
            colored_bases: false,
            compact: false,
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SequenceRestrictionSite {
    pub enzyme_index: usize,
    pub start: usize,
    pub top_cut: Option<usize>,
    pub bottom_cut: Option<usize>,
    pub label: String,
}
#[derive(Clone, Debug)]
pub enum SequenceEvent {
    Select(Vec<Segment>),
    SelectFeature(u64),
    InspectEnzyme { enzyme_index: usize, start: usize },
}
struct EnzymeLabel {
    site: usize,
    x: f32,
    width: f32,
    lane: usize,
}
struct RowPlan {
    enzymes: Vec<EnzymeLabel>,
    cuts: Vec<(usize, bool)>,
    omitted_enzymes: usize,
    top: f32,
    lane_y: Vec<f32>,
    frame_y: Vec<(i8, f32)>,
    height: f32,
}

type VisibleRowBounds = Rc<RefCell<Vec<(usize, Bounds<Pixels>)>>>;

pub struct SequenceView {
    document: Rc<Document>,
    selected: Vec<Segment>,
    options: SequenceOptions,
    sites: Vec<SequenceRestrictionSite>,
    layout: SequenceLayout,
    plans: Vec<RowPlan>,
    offsets: Vec<f32>,
    sizes: Rc<Vec<Size<Pixels>>>,
    per_row: usize,
    viewport_width: f32,
    pending_scroll_to: Option<usize>,
    scroll: VirtualListScrollHandle,
    row_bounds: VisibleRowBounds,
    anchor: Option<usize>,
    dragging: Option<usize>,
    focus: FocusHandle,
}
impl EventEmitter<SequenceEvent> for SequenceView {}
impl SequenceView {
    pub fn new(
        document: Document,
        selected: Vec<Segment>,
        _: &mut Window,
        cx: &mut Context<Self>,
    ) -> Self {
        let layout = SequenceLayout::new(&document);
        let mut view = Self {
            document: Rc::new(document),
            selected,
            options: SequenceOptions::default(),
            sites: vec![],
            layout,
            plans: vec![],
            offsets: vec![],
            sizes: Rc::new(vec![]),
            per_row: 60,
            viewport_width: 800.,
            pending_scroll_to: None,
            scroll: VirtualListScrollHandle::new(),
            row_bounds: Rc::new(RefCell::new(vec![])),
            anchor: None,
            dragging: None,
            focus: cx.focus_handle(),
        };
        view.remeasure();
        view
    }
    fn cell_width(&self) -> f32 {
        if self.options.compact { 9. } else { 11. }
    }
    fn paired(&self) -> bool {
        self.options.show_complement
            && !self.document.single_stranded
            && self.document.molecule != Molecule::Protein
    }
    pub fn set_document(
        &mut self,
        document: Document,
        selected: Vec<Segment>,
        scroll_to: Option<usize>,
        options: SequenceOptions,
        sites: Vec<SequenceRestrictionSite>,
        cx: &mut Context<Self>,
    ) {
        let document_changed = *self.document != document;
        let layout_changed = document_changed || self.options != options || self.sites != sites;
        let selection_changed = self.selected != selected;
        if selection_changed {
            self.anchor = selected.first().map(|segment| segment.start);
            self.dragging = None;
        }
        if document_changed {
            self.layout = SequenceLayout::new(&document);
            self.document = Rc::new(document);
            self.anchor = None;
            self.dragging = None;
        }
        self.options = options;
        self.sites = sites;
        self.selected = selected;
        if layout_changed {
            let base = self.visible_base();
            self.per_row = sequence_layout::bases_per_row(
                self.viewport_width,
                self.cell_width(),
                LEFT + RIGHT,
            );
            self.remeasure();
            self.scroll
                .scroll_to_item(base / self.per_row, ScrollStrategy::Top);
        }
        if let Some(base) = scroll_to {
            self.pending_scroll_to = Some(base.min(self.document.sequence.len().saturating_sub(1)));
        }
        if layout_changed || selection_changed || scroll_to.is_some() {
            cx.notify();
        }
    }
    fn visible_base(&self) -> usize {
        let y = -self.scroll.offset().y.as_f32();
        self.offsets
            .partition_point(|offset| *offset <= y)
            .saturating_sub(1)
            * self.per_row
    }
    fn resize(&mut self, width: f32, cx: &mut Context<Self>) {
        self.viewport_width = width;
        let per_row = sequence_layout::bases_per_row(width, self.cell_width(), LEFT + RIGHT);
        let changed = per_row != self.per_row;
        let base = self.visible_base();
        if changed {
            self.per_row = per_row;
            self.remeasure();
        }
        // Resolve navigation only after the sidebar and viewport have their final width.
        if let Some(target) = self.pending_scroll_to.take() {
            let target = target.min(self.document.sequence.len().saturating_sub(1));
            self.scroll
                .scroll_to_item(target / per_row, ScrollStrategy::Center);
            cx.notify();
        } else if changed {
            self.scroll
                .scroll_to_item(base / per_row, ScrollStrategy::Top);
            cx.notify();
        }
    }
    fn remeasure(&mut self) {
        let n = self.document.sequence.len();
        let rows = n.div_ceil(self.per_row);
        let cw = self.cell_width();
        let mut row_sites = vec![Vec::new(); rows];
        let mut cuts = vec![Vec::new(); rows];
        if self.options.show_enzymes
            && self.document.molecule == Molecule::Dna
            && !self.document.single_stranded
        {
            for (i, site) in self.sites.iter().enumerate() {
                let boundary = site.top_cut.unwrap_or(site.start).min(n.saturating_sub(1));
                if let Some(row) = row_sites.get_mut(boundary / self.per_row) {
                    row.push(i);
                }
                for (boundary, top) in [(site.top_cut, true), (site.bottom_cut, false)] {
                    if let Some(boundary) = boundary.filter(|boundary| *boundary <= n) {
                        let row = boundary.min(n.saturating_sub(1)) / self.per_row;
                        if let Some(values) = cuts.get_mut(row) {
                            values.push((boundary - row * self.per_row, top));
                        }
                    }
                }
            }
        }
        self.plans = (0..rows)
            .map(|row| {
                let start = row * self.per_row;
                let end = (start + self.per_row).min(n);
                let annotations = self.layout.row_summary(
                    start,
                    end,
                    MAX_LANES,
                    self.options.show_features,
                    self.options.show_primers,
                    self.options.show_translations,
                );
                let mut enzymes = vec![];
                let mut lane_end = [0_f32; 2];
                let mut omitted_enzymes = 0;
                row_sites[row]
                    .sort_by_key(|i| self.sites[*i].top_cut.unwrap_or(self.sites[*i].start));
                for &i in &row_sites[row] {
                    let site = &self.sites[i];
                    let x = (site.top_cut.unwrap_or(site.start).saturating_sub(start) as f32 * cw)
                        .min(self.per_row as f32 * cw);
                    let width = (site.label.chars().count() as f32 * 7. + 12.).min(155.);
                    let x = x.min((self.per_row as f32 * cw - width).max(0.));
                    if let Some(lane) = lane_end.iter().position(|end| *end <= x) {
                        lane_end[lane] = x + width + 8.;
                        enzymes.push(EnzymeLabel {
                            site: i,
                            x,
                            width,
                            lane,
                        });
                    } else {
                        omitted_enzymes += 1;
                    }
                }
                let enzyme_lanes = enzymes.iter().map(|e| e.lane + 1).max().unwrap_or(0);
                let top = 34. + enzyme_lanes as f32 * 22.;
                let annotations_y = top + if self.paired() { 61. } else { 31. };
                let mut lane_y = vec![];
                let mut y = annotations_y;
                let mut frame_y = vec![];
                if self.document.molecule != Molecule::Protein {
                    for (bit, frame) in [1, 2, 3, -1, -2, -3].into_iter().enumerate() {
                        if self.options.reading_frames & (1 << bit) != 0 {
                            frame_y.push((frame, y));
                            y += 24.;
                        }
                    }
                }
                for lane in 0..annotations.lanes {
                    lane_y.push(y);
                    let translated = annotations.translation_lanes.contains(&lane);
                    y += if translated { 43. } else { 24. };
                }
                let foot = if annotations.omitted > 0 || omitted_enzymes > 0 {
                    20.
                } else {
                    0.
                };
                RowPlan {
                    enzymes,
                    cuts: std::mem::take(&mut cuts[row]),
                    omitted_enzymes,
                    top,
                    lane_y,
                    frame_y,
                    height: y + 22. + foot,
                }
            })
            .collect();
        self.offsets.clear();
        let mut y = 0.;
        self.sizes = Rc::new(
            self.plans
                .iter()
                .map(|plan| {
                    self.offsets.push(y);
                    y += plan.height;
                    size(px(0.), px(plan.height))
                })
                .collect(),
        );
    }
    fn select_base(&mut self, index: usize, extend: bool, cx: &mut Context<Self>) {
        let anchor = if extend {
            self.anchor
                .or_else(|| self.selected.first().map(|s| s.start))
                .unwrap_or(index)
        } else {
            index
        };
        self.anchor = Some(anchor);
        self.dragging = Some(anchor);
        self.extend_selection(anchor, index, cx);
    }
    fn extend_selection(&mut self, anchor: usize, index: usize, cx: &mut Context<Self>) {
        let selected = vec![Segment {
            start: anchor.min(index),
            end: anchor.max(index) + 1,
        }];
        if self.selected != selected {
            self.selected = selected.clone();
            cx.emit(SequenceEvent::Select(selected));
            cx.notify();
        }
    }
    fn base_cell(
        &self,
        index: usize,
        base: u8,
        y: f32,
        opposite: bool,
        cx: &mut Context<Self>,
    ) -> impl IntoElement {
        let cw = self.cell_width();
        let symbol = if self.document.molecule == Molecule::Protein {
            base as char
        } else {
            (base as char).to_ascii_lowercase()
        };
        div()
            .id((
                if opposite {
                    "complement-base"
                } else {
                    "sequence-base"
                },
                index,
            ))
            .absolute()
            .left(px(LEFT + sequence_layout::base_x(index % self.per_row, cw)))
            .top(px(y))
            .w(px(cw))
            .h(px(22.))
            .flex()
            .items_center()
            .justify_center()
            .cursor_text()
            .text_color(if self.options.colored_bases {
                base_color(base)
            } else {
                rgb(INK)
            })
            .hover(|s| s.bg(rgba(0x7cbbc62a)))
            .on_mouse_down(
                MouseButton::Left,
                cx.listener(move |this, event: &MouseDownEvent, window, cx| {
                    this.focus.focus(window, cx);
                    this.select_base(index, event.modifiers.shift, cx);
                    cx.stop_propagation();
                }),
            )
            .child(symbol.to_string())
    }
    fn render_row(&self, row: usize, cx: &mut Context<Self>) -> AnyElement {
        let plan = &self.plans[row];
        let start = row * self.per_row;
        let end = (start + self.per_row).min(self.document.sequence.len());
        let annotations = self.layout.row(
            start,
            end,
            MAX_LANES,
            self.options.show_features,
            self.options.show_primers,
            self.options.show_translations,
        );
        let cw = self.cell_width();
        let width = self.per_row as f32 * cw;
        let paired = self.paired();
        let mut content = div()
            .id(("sequence-row", row))
            .relative()
            .h(px(plan.height))
            .w_full()
            .font_family(cx.theme().mono_font_family.clone())
            .text_size(px(if self.options.compact { 13. } else { 15. }));
        let rows = self.row_bounds.clone();
        let top = plan.top;
        let cuts = plan.cuts.clone();
        let seq_n = self.document.sequence.len();
        // Measurements stay bounded to the rendered rows; no per-base layout is retained.
        content = content.child(
            canvas(
                move |bounds, _, _| {
                    rows.borrow_mut().push((row, bounds));
                },
                move |bounds, _, window, _| {
                    let left = bounds.left() + px(LEFT);
                    let ruler = bounds.top() + px(top - 8.);
                    let mut line = PathBuilder::stroke(px(1.));
                    line.move_to(point(left, ruler));
                    line.line_to(point(left + px((end - start) as f32 * cw), ruler));
                    for column in 0..=(end - start) {
                        let x = left + px(column as f32 * cw);
                        let major = (start + column).is_multiple_of(10);
                        line.move_to(point(x, ruler));
                        line.line_to(point(x, ruler + px(if major { 6. } else { 3. })));
                    }
                    if let Ok(path) = line.build() {
                        window.paint_path(path, rgb(0x4ca7bc));
                    }
                    if paired {
                        let mid = bounds.top() + px(top + 26.);
                        let mut joins = PathBuilder::stroke(px(1.));
                        for column in 0..(end - start) {
                            let x = left + px((column as f32 + 0.5) * cw);
                            joins.move_to(point(x, mid - px(2.)));
                            joins.line_to(point(x, mid + px(2.)));
                        }
                        if let Ok(path) = joins.build() {
                            window.paint_path(path, rgba(0x8799ab66));
                        }
                    }
                    let mut path = PathBuilder::stroke(px(1.5));
                    for (column, forward) in &cuts {
                        if !forward && !paired {
                            continue;
                        }
                        let x = left + px(*column as f32 * cw);
                        let y = bounds.top() + px(top + if *forward { 0. } else { 30. });
                        path.move_to(point(x, y));
                        path.line_to(point(x, y + px(21.)));
                        let direction = if *forward { 1. } else { -1. };
                        path.line_to(point(x + px(direction * 4.), y + px(21.)));
                    }
                    if !cuts.is_empty()
                        && seq_n > 0
                        && let Ok(path) = path.build()
                    {
                        window.paint_path(path, rgb(ACCENT));
                    }
                },
            )
            .absolute()
            .size_full(),
        );
        for selection in &self.selected {
            let a = selection.start.max(start);
            let b = selection.end.min(end);
            if a < b {
                content = content.child(
                    div()
                        .absolute()
                        .left(px(LEFT + (a - start) as f32 * cw))
                        .top(px(top))
                        .w(px((b - a) as f32 * cw))
                        .h(px(if paired { 52. } else { 22. }))
                        .bg(rgba(0x28665b99)),
                );
            }
        }
        for column in 0..=(end - start) {
            if (start + column).is_multiple_of(10) && column > 0 {
                content = content.child(
                    div()
                        .absolute()
                        .left(px(LEFT + column as f32 * cw - 25.))
                        .top(px(top - 25.))
                        .w(px(25.))
                        .text_right()
                        .text_size(px(10.))
                        .text_color(rgb(0x78b6c7))
                        .child((start + column).to_string()),
                );
            }
        }
        let nucleotide = self.document.molecule != Molecule::Protein;
        let left_label = if nucleotide { "5′" } else { "N" };
        let right_label = if nucleotide { "3′" } else { "C" };
        content = content
            .child(
                div()
                    .absolute()
                    .left(px(20.))
                    .top(px(top + 3.))
                    .text_size(px(11.))
                    .text_color(rgb(MUTED))
                    .child(left_label),
            )
            .child(
                div()
                    .absolute()
                    .left(px(LEFT + width + 10.))
                    .top(px(top + 19.))
                    .text_size(px(10.))
                    .text_color(rgb(MUTED))
                    .child(format!("{end}")),
            )
            .child(
                div()
                    .absolute()
                    .left(px(LEFT + width + 10.))
                    .top(px(top + 3.))
                    .text_size(px(9.))
                    .text_color(rgb(MUTED))
                    .child(right_label),
            );
        if paired {
            content = content
                .child(
                    div()
                        .absolute()
                        .left(px(20.))
                        .top(px(top + 33.))
                        .text_size(px(11.))
                        .text_color(rgb(MUTED))
                        .child("3′"),
                )
                .child(
                    div()
                        .absolute()
                        .left(px(LEFT + width + 10.))
                        .top(px(top + 33.))
                        .text_size(px(9.))
                        .text_color(rgb(MUTED))
                        .child("5′"),
                );
        }
        for index in start..end {
            let base = self.document.sequence.as_bytes()[index];
            content = content.child(self.base_cell(index, base, top, false, cx));
            if paired {
                content = content.child(self.base_cell(
                    index,
                    sequence_layout::complement(base, self.document.molecule),
                    top + 30.,
                    true,
                    cx,
                ));
            }
        }
        for label in &plan.enzymes {
            let site = &self.sites[label.site];
            let enzyme_index = site.enzyme_index;
            let recognition_start = site.start;
            content = content.child(
                div()
                    .absolute()
                    .left(px(LEFT + label.x))
                    .top(px(5. + label.lane as f32 * 22.))
                    .w(px(label.width))
                    .child(
                        Button::new(("sequence-enzyme", row * self.sites.len() + label.site))
                            .ghost()
                            .xsmall()
                            .label(site.label.clone())
                            .accessibility_label(format!(
                                "Inspect {} at recognition base {}",
                                site.label,
                                site.start + 1
                            ))
                            .on_click(cx.listener(move |_, _, _, cx| {
                                cx.emit(SequenceEvent::InspectEnzyme {
                                    enzyme_index,
                                    start: recognition_start,
                                })
                            })),
                    ),
            );
        }
        for (part, span) in annotations.spans.iter().enumerate() {
            let feature = &self.document.features[span.feature_index];
            let x = LEFT + span.start as f32 * cw;
            let length = (span.end - span.start) as f32 * cw;
            let y = plan.lane_y[span.lane];
            let color = span
                .color_rgb
                .map(|color| rgb(color).into())
                .unwrap_or_else(|| feature_card::feature_color(feature, span.feature_index));
            let label_color = annotation_label_color(color);
            let direction = feature_card::direction(feature, self.document.molecule);
            let reverse = direction == "Reverse strand";
            let directional = direction != "Nondirectional" && nucleotide;
            let bidirectional = direction == "Bidirectional";
            let arrow_left = directional && (reverse || bidirectional) && !span.continued_left;
            let arrow_right = directional && (!reverse || bidirectional) && !span.continued_right;
            let id = feature.id;
            let selected = self.selected == feature.parts;
            let doc = self.document.clone();
            let f = feature.clone();
            let trigger = div()
                .id(("sequence-feature", row * 1000 + part))
                .relative()
                .w(px(length))
                .h(px(19.))
                .cursor_pointer()
                .on_mouse_down(
                    MouseButton::Left,
                    cx.listener(move |_, _, _, cx| {
                        cx.emit(SequenceEvent::SelectFeature(id));
                        cx.stop_propagation();
                    }),
                )
                .child(
                    canvas(
                        |_, _, _| (),
                        move |bounds, _, window, _| {
                            let x = bounds.left();
                            let y = bounds.top();
                            let w = bounds.size.width;
                            let h = bounds.size.height;
                            let tip = px(7.).min(w * 0.5);
                            let mut p = PathBuilder::fill();
                            p.move_to(point(x + if arrow_left { tip } else { px(0.) }, y));
                            p.line_to(point(x + w - if arrow_right { tip } else { px(0.) }, y));
                            if arrow_right {
                                p.line_to(point(x + w, y + h * 0.5));
                            }
                            p.line_to(point(x + w - if arrow_right { tip } else { px(0.) }, y + h));
                            p.line_to(point(x + if arrow_left { tip } else { px(0.) }, y + h));
                            if arrow_left {
                                p.line_to(point(x, y + h * 0.5));
                            }
                            p.close();
                            if let Ok(path) = p.build() {
                                window.paint_path(path, color);
                            }
                        },
                    )
                    .absolute()
                    .size_full(),
                )
                .child(
                    div()
                        .absolute()
                        .inset_0()
                        .px_2()
                        .overflow_hidden()
                        .text_ellipsis()
                        .text_center()
                        .text_size(px(11.))
                        .font_family(cx.theme().font_family.clone())
                        .font_weight(FontWeight::SEMIBOLD)
                        .text_color(label_color)
                        .when(selected, |d| d.border_1().border_color(rgb(INK)))
                        .child(if length >= 32. {
                            feature_card::truncate(
                                span.segment_name.as_deref().unwrap_or(&feature.label),
                                100,
                            )
                        } else {
                            String::new()
                        }),
                );
            content = content.child(
                div()
                    .absolute()
                    .left(px(x))
                    .top(px(y))
                    .w(px(length))
                    .h(px(19.))
                    .child(
                        HoverCard::new(("sequence-feature-card", row * 1000 + part))
                            .content(move |_, _, _| {
                                feature_card::render_feature_card(&doc, &f, None)
                            })
                            .trigger(trigger),
                    ),
            );
        }
        for (frame, y) in &plan.frame_y {
            let color = if *frame > 0 { 0x8fd5b3 } else { 0x9bbce8 };
            content = content.child(
                div()
                    .absolute()
                    .left(px(4.))
                    .top(px(*y))
                    .w(px(47.))
                    .h(px(18.))
                    .text_right()
                    .text_size(px(11.))
                    .text_color(rgb(color))
                    .child(format!("{frame:+} {}", if *frame > 0 { "→" } else { "←" })),
            );
            for codon in frame_codons(&self.document.sequence, *frame, start..end) {
                let left = codon.start.max(start) - start;
                let right = codon.end.min(end) - start;
                let a = codon.start;
                let b = codon.end;
                content = content.child(
                    div()
                        .id((
                            "frame-codon",
                            (usize::from(frame.unsigned_abs()) + if *frame < 0 { 3 } else { 0 })
                                * self.document.sequence.len()
                                + codon.center,
                        ))
                        .absolute()
                        .left(px(LEFT + left as f32 * cw))
                        .top(px(*y))
                        .w(px((right - left) as f32 * cw))
                        .h(px(18.))
                        .cursor_pointer()
                        .bg(rgba(if *frame > 0 { 0x8fd5b31a } else { 0x9bbce81a }))
                        .hover(|s| s.bg(rgba(0xb5d4e544)))
                        .on_mouse_down(
                            MouseButton::Left,
                            cx.listener(move |_, _, _, cx| {
                                cx.emit(SequenceEvent::Select(vec![Segment { start: a, end: b }]));
                                cx.stop_propagation();
                            }),
                        )
                        .child(
                            div()
                                .absolute()
                                .left(px((codon.center - start - left) as f32 * cw))
                                .w(px(cw))
                                .text_center()
                                .text_size(px(11.))
                                .text_color(rgb(color))
                                .child((codon.amino_acid as char).to_string()),
                        ),
                );
            }
        }
        for glyph in &annotations.translations {
            let y = plan.lane_y[glyph.lane] + 20.;
            let details = format!(
                "Calculated CDS · genetic code {} · bases {}, {}, {}{}",
                glyph.genetic_code,
                glyph.codon_positions[0] + 1,
                glyph.codon_positions[1] + 1,
                glyph.codon_positions[2] + 1,
                if glyph.split_codon {
                    " · joined codon; only its central base is highlighted"
                } else {
                    ""
                }
            );
            content = content.child(
                div()
                    .absolute()
                    .left(px(LEFT + glyph.start as f32 * cw))
                    .top(px(y))
                    .w(px((glyph.end - glyph.start) as f32 * cw))
                    .h(px(18.))
                    .bg(rgba(0xe2ac6433)),
            );
            content = content.child(
                div()
                    .absolute()
                    .left(px(LEFT + glyph.base as f32 * cw))
                    .top(px(y))
                    .w(px(cw))
                    .h(px(18.))
                    .text_center()
                    .text_size(px(11.))
                    .text_color(rgb(0xf0ca8c))
                    .child(
                        HoverCard::new((
                            "calculated-cds",
                            (glyph.feature_index * self.document.sequence.len()
                                + start
                                + glyph.base)
                                * MAX_LANES
                                + glyph.lane,
                        ))
                        .trigger(div().child(glyph.amino_acid.to_string()))
                        .content(move |_, _, _| {
                            div()
                                .p_3()
                                .max_w(px(360.))
                                .text_xs()
                                .child(details.clone())
                                .into_any_element()
                        }),
                    ),
            );
        }
        let mut omitted = vec![];
        if annotations.omitted > 0 {
            omitted.push(format!("{} more annotation segments", annotations.omitted));
        }
        if plan.omitted_enzymes > 0 {
            omitted.push(format!("{} more enzyme labels", plan.omitted_enzymes));
        }
        if !omitted.is_empty() {
            content = content.child(
                div()
                    .absolute()
                    .left(px(LEFT))
                    .top(px(plan.height - 20.))
                    .text_size(px(10.))
                    .text_color(rgb(MUTED))
                    .child(format!(
                        "{} · use the feature/enzyme list",
                        omitted.join(" · ")
                    )),
            );
        }
        content.into_any_element()
    }
}
impl Render for SequenceView {
    fn render(&mut self, _: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        let rows = self.row_bounds.clone();
        let owner = cx.entity().downgrade();
        let cell_width = self.cell_width();
        let current = self.per_row;
        let measured_width = self.viewport_width;
        let navigation_pending = self.pending_scroll_to.is_some();
        let measure = canvas(
            move |bounds, _, cx| {
                rows.borrow_mut().clear();
                let width = bounds.size.width.as_f32();
                if navigation_pending
                    || (width - measured_width).abs() > 0.5
                    || sequence_layout::bases_per_row(width, cell_width, LEFT + RIGHT) != current
                {
                    let owner = owner.clone();
                    cx.defer(move |cx| {
                        let _ = owner.update(cx, |view, cx| view.resize(width, cx));
                    });
                }
            },
            |_, _, _, _| (),
        )
        .absolute()
        .size_full();
        div()
            .id("annotated-sequence")
            .relative()
            .size_full()
            .min_h_0()
            .bg(rgb(0x111a24))
            .track_focus(&self.focus)
            .on_mouse_move(cx.listener(|this, event: &MouseMoveEvent, _, cx| {
                if !event.dragging() {
                    return;
                }
                let Some(anchor) = this.dragging else {
                    return;
                };
                let hit = this.row_bounds.borrow().iter().find_map(|(row, bounds)| {
                    if event.position.y >= bounds.top() && event.position.y < bounds.bottom() {
                        let count =
                            (this.document.sequence.len() - row * this.per_row).min(this.per_row);
                        let x = (event.position.x - bounds.left()).as_f32() - LEFT;
                        let x = x.clamp(0., (count as f32 * this.cell_width() - 0.01).max(0.));
                        sequence_layout::base_at_x(x, this.cell_width(), count)
                            .map(|column| row * this.per_row + column)
                    } else {
                        None
                    }
                });
                if let Some(index) = hit {
                    this.extend_selection(anchor, index, cx);
                }
            }))
            .on_mouse_up(
                MouseButton::Left,
                cx.listener(|this, _, _, _| this.dragging = None),
            )
            .on_mouse_up_out(
                MouseButton::Left,
                cx.listener(|this, _, _, _| this.dragging = None),
            )
            .child(measure)
            .child(
                v_virtual_list(
                    cx.entity(),
                    "sequence-rows",
                    self.sizes.clone(),
                    |this, range, _, cx| {
                        range
                            .map(|row| this.render_row(row, cx))
                            .collect::<Vec<_>>()
                    },
                )
                .size_full()
                .track_scroll(&self.scroll),
            )
            .child(Scrollbar::vertical(&self.scroll))
    }
}
impl SequenceView {
    pub fn translation_legend(&self) -> AnyElement {
        let show_cds = self.document.molecule != Molecule::Protein
            && self.options.show_translations
            && self
                .document
                .features
                .iter()
                .any(|feature| feature.kind.eq_ignore_ascii_case("CDS"));
        let show_frames =
            self.options.reading_frames != 0 && self.document.molecule != Molecule::Protein;
        div().v_flex().w_full().min_w_0()
            .when(show_cds || show_frames, |d| {
                d.child(div().id("translation-status").v_flex().gap_2().py_2().text_xs().text_color(rgb(MUTED))
                    .when(show_cds, |d| d.child("Calculated CDS · recorded frame and genetic code (1, 2, 11 supported)"))
                    .when(show_frames, |d| d.child("Reading frames · standard code 1 · reverse frames start at the last base"))
                    .child("* stop · X unresolved codon · incomplete end codons omitted · ORFs are not inferred")
                    .when(show_cds && self.layout.omitted_translations > 0, |d| {
                        d.child(div().text_color(rgb(0xe0b472)).child(format!("{} CDS translations omitted:", self.layout.omitted_translations)))
                            .children(self.layout.translation_notices.iter().map(|notice| div().child(notice.clone())))
                    }))
            })
            .into_any_element()
    }
}

pub fn base_color(base: u8) -> Rgba {
    rgb(match base {
        b'A' => 0x61d5a1,
        b'C' => 0x70b6fa,
        b'G' => 0xe5bf6c,
        b'T' | b'U' => 0xee8da1,
        _ => 0xaab8c7,
    })
}

/// Choose the higher-contrast black/white foreground from linear sRGB
/// luminance. Saturated blue is dark despite having the same HSL lightness
/// as bright green, so HSL lightness alone is unsuitable for imported colors.
fn annotation_label_color(background: Hsla) -> Rgba {
    let color = background.to_rgb();
    let linear = |channel: f32| {
        if channel <= 0.04045 {
            channel / 12.92
        } else {
            ((channel + 0.055) / 1.055).powf(2.4)
        }
    };
    let luminance = 0.2126 * linear(color.r) + 0.7152 * linear(color.g) + 0.0722 * linear(color.b);
    let black_contrast = (luminance + 0.05) / 0.05;
    let white_contrast = 1.05 / (luminance + 0.05);
    rgb(if white_contrast > black_contrast {
        0xffffff
    } else {
        0x000000
    })
}

/// Each trace canvas paints a bounded sample interval. The horizontal virtual
/// list creates only visible intervals, including for large sequencing files.
pub struct TraceView {
    trace: Option<Arc<Trace>>,
    sequence: Arc<str>,
    scale: f32,
    max_signal: f32,
    sizes: Rc<Vec<Size<Pixels>>>,
    scroll: VirtualListScrollHandle,
}
const TRACE_CHUNK: usize = 256;
impl TraceView {
    pub fn new(document: Document, scale: f32, _: &mut Context<Self>) -> Self {
        let mut view = Self {
            trace: document.trace.map(Arc::new),
            sequence: document.sequence,
            scale: scale.clamp(0.5, 5.),
            max_signal: 1.,
            sizes: Rc::new(vec![]),
            scroll: VirtualListScrollHandle::new(),
        };
        view.remeasure();
        view
    }
    pub fn set_document(&mut self, document: Document, scale: f32, cx: &mut Context<Self>) {
        let scale = scale.clamp(0.5, 5.);
        let changed = self.trace.as_deref() != document.trace.as_ref()
            || self.sequence != document.sequence
            || self.scale != scale;
        if changed {
            self.trace = document.trace.map(Arc::new);
            self.sequence = document.sequence;
            self.scale = scale;
            self.remeasure();
            cx.notify();
        }
    }
    fn samples(&self) -> usize {
        self.trace
            .as_ref()
            .and_then(|t| t.channels.first())
            .map_or(0, Vec::len)
    }
    fn sample_width(&self) -> f32 {
        (800. / self.samples().max(1) as f32).max(2.) * self.scale
    }
    fn remeasure(&mut self) {
        self.max_signal = self
            .trace
            .as_ref()
            .into_iter()
            .flat_map(|t| &t.channels)
            .flat_map(|c| c.iter())
            .copied()
            .max()
            .unwrap_or(1)
            .max(1) as f32;
        let samples = self.samples();
        let width = self.sample_width();
        self.sizes = Rc::new(
            (0..samples.div_ceil(TRACE_CHUNK))
                .map(|i| {
                    size(
                        px(((samples - i * TRACE_CHUNK).min(TRACE_CHUNK)) as f32 * width),
                        px(420.),
                    )
                })
                .collect(),
        );
    }
    fn render_chunk(&self, chunk: usize) -> impl IntoElement {
        let trace = self.trace.clone();
        let sequence = self.sequence.clone();
        let max_signal = self.max_signal;
        let width = self.sample_width();
        let samples = self.samples();
        let start = chunk * TRACE_CHUNK;
        let end = (start + TRACE_CHUNK).min(samples);
        let cell_width = (end - start) as f32 * width;
        div().w(px(cell_width)).h(px(420.)).flex_shrink_0().child(
            canvas(
                |_, _, _| (),
                move |bounds, _, window, cx| {
                    let Some(trace) = &trace else {
                        return;
                    };
                    let x = |i: usize| bounds.left() + px((i as f32 - start as f32) * width);
                    let baseline = bounds.bottom() - px(38.);
                    let amplitude = bounds.size.height - px(126.);
                    let mut grid = PathBuilder::stroke(px(1.));
                    grid.move_to(point(bounds.left(), baseline));
                    grid.line_to(point(bounds.right(), baseline));
                    if let Ok(path) = grid.build() {
                        window.paint_path(path, rgb(0x334151));
                    }
                    // One overlapping point makes adjacent chunk paths continuous.
                    let line_start = start.saturating_sub(1);
                    let line_end = (end + 1).min(samples);
                    for (channel, base) in trace.channels.iter().zip(trace.channel_order.bytes()) {
                        let mut path = PathBuilder::stroke(px(1.5));
                        for (i, sample) in
                            channel.iter().enumerate().take(line_end).skip(line_start)
                        {
                            let point =
                                point(x(i), baseline - amplitude * (*sample as f32 / max_signal));
                            if i == line_start {
                                path.move_to(point);
                            } else {
                                path.line_to(point);
                            }
                        }
                        if let Ok(path) = path.build() {
                            window.paint_path(path, base_color(base));
                        }
                    }
                    // Peaks are ordered by the ABIF adapter. Skip directly to this
                    // visible interval rather than scanning every call per frame.
                    let first = trace.peaks.partition_point(|&p| p < start);
                    for i in first..trace.peaks.len() {
                        let peak = trace.peaks[i];
                        if peak >= end {
                            break;
                        }
                        let Some(&base) = sequence.as_bytes().get(i) else {
                            continue;
                        };
                        paint_text(
                            (base as char).to_string(),
                            point(x(peak) - px(5.), bounds.top() + px(18.)),
                            16.,
                            base_color(base),
                            window,
                            cx,
                        );
                        if let Some(quality) = trace.qualities.get(i) {
                            paint_text(
                                format!("Q{quality}"),
                                point(x(peak) - px(8.), bounds.top() + px(42.)),
                                10.,
                                rgb(0x8c9bad),
                                window,
                                cx,
                            );
                        }
                        paint_text(
                            (i + 1).to_string(),
                            point(x(peak) - px(5.), baseline + px(10.)),
                            10.,
                            rgb(0x8c9bad),
                            window,
                            cx,
                        );
                    }
                },
            )
            .size_full(),
        )
    }
}
impl Render for TraceView {
    fn render(&mut self, _: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        let label = if let Some(trace) = &self.trace {
            format!(
                "{} base calls  ·  {} signal samples  ·  PBAS{}  ·  read-only",
                self.sequence.len(),
                self.samples(),
                trace.calls_version
            )
        } else {
            "No Sanger trace in this document".into()
        };
        div()
            .v_flex()
            .size_full()
            .gap_4()
            .p_5()
            .bg(rgb(0x151d27))
            .child(div().text_sm().text_color(rgb(0x8c9bad)).child(label))
            .child(
                h_virtual_list(
                    cx.entity(),
                    "trace-chunks",
                    self.sizes.clone(),
                    |this, range, _, _| {
                        range
                            .map(|chunk| this.render_chunk(chunk).into_any_element())
                            .collect()
                    },
                )
                .h(px(430.))
                .w_full()
                .track_scroll(&self.scroll),
            )
    }
}
fn paint_text(
    text: String,
    origin: Point<Pixels>,
    font_size: f32,
    color: Rgba,
    window: &mut Window,
    cx: &mut App,
) {
    let text: SharedString = text.into();
    let run = TextRun {
        len: text.len(),
        font: window.text_style().font(),
        color: color.into(),
        background_color: None,
        underline: None,
        strikethrough: None,
    };
    let line = window
        .text_system()
        .shape_line(text, px(font_size), &[run], None);
    let _ = line.paint(origin, px(font_size), TextAlign::Left, None, window, cx);
}

#[cfg(test)]
mod contrast_tests {
    use super::{annotation_label_color, rgb};
    #[test]
    fn imported_annotation_colors_get_readable_label_foregrounds() {
        for background in [0x0000ff, 0x000000, 0x993366, 0x444444] {
            assert_eq!(
                annotation_label_color(rgb(background).into()),
                rgb(0xffffff)
            );
        }
        for background in [0x00ff00, 0xff0000, 0xffffff, 0x4cc7aa, 0xa6acb3] {
            assert_eq!(
                annotation_label_color(rgb(background).into()),
                rgb(0x000000)
            );
        }
    }
}
