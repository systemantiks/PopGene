//! Virtualized annotated sequence and Sanger trace views, rendered directly by GPUI.
use crate::{
    feature_card,
    sequence_layout::{self, SequenceLayout},
};
use gpui_kit::component::scroll::Scrollbar;
use gpui_kit::component::{
    ActiveTheme, Disableable, Sizable, StyledExt, VirtualListScrollHandle, button::*,
    h_virtual_list, hover_card::HoverCard, v_virtual_list,
};
use gpui_kit::{prelude::*, *};
use popgene_core::{
    Document, Molecule, Segment, Trace,
    translation::{
        FlankingCodon, GeneticCode, OrfOptions, OrfSearch, cds_flanking_codons, find_orfs,
        frame_codons_with_options,
    },
};
use std::{cell::RefCell, rc::Rc, sync::Arc};

const LEFT: f32 = 58.;
const RIGHT: f32 = 66.;
const INK: u32 = 0xdce5ed;
const MUTED: u32 = 0x8799ab;
const ACCENT: u32 = 0x48c9bb;
const MAX_LANES: usize = 4;

#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(default)]
pub struct SequenceOptions {
    pub show_complement: bool,
    pub show_features: bool,
    pub show_primers: bool,
    pub show_enzymes: bool,
    pub show_translations: bool,
    /// Bits 0..2: forward +1/+2/+3; bits 3..5: reverse -1/-2/-3.
    pub reading_frames: u8,
    pub genetic_code: u16,
    pub show_orfs: bool,
    pub orf_min_codons: usize,
    pub orf_alternative_starts: bool,
    pub orf_include_partial: bool,
    pub show_flanking: bool,
    pub flanking_codons: usize,
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
            genetic_code: 1,
            show_orfs: false,
            orf_min_codons: 30,
            orf_alternative_starts: false,
            orf_include_partial: false,
            show_flanking: false,
            flanking_codons: 30,
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
    pub second_top_cut: Option<usize>,
    pub second_bottom_cut: Option<usize>,
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
    orf_y: Vec<(i8, f32)>,
    flank_y: Vec<(bool, f32)>,
    height: f32,
}

type VisibleRowBounds = Rc<RefCell<Vec<(usize, Bounds<Pixels>)>>>;

#[derive(Default)]
struct FlankingTracks {
    // Keep the rendered CDS key separate from the optimistic mouse selection.
    // Non-CDS drag steps all share None and do not rebuild the row layout.
    key: Option<(usize, usize)>,
    codons: Vec<FlankingCodon>,
    label: Option<String>,
}
impl FlankingTracks {
    fn update(
        &mut self,
        document: &Document,
        selected: &[Segment],
        options: SequenceOptions,
        document_changed: bool,
    ) -> bool {
        let key = options
            .show_flanking
            .then(|| {
                document.features.iter().position(|feature| {
                    feature.kind.eq_ignore_ascii_case("CDS") && feature.parts == selected
                })
            })
            .flatten()
            .map(|index| (index, options.flanking_codons));
        if !document_changed && self.key == key {
            return false;
        }
        self.key = key;
        let (codons, label) = key.map_or_else(
            || (vec![], None),
            |(index, count)| {
                let feature = &document.features[index];
                match cds_flanking_codons(document, feature, count) {
                    Ok(codons) => (codons, Some(feature.label.clone())),
                    Err(error) => (vec![], Some(format!("{}: {error}", feature.label))),
                }
            },
        );
        let changed = self.codons != codons || self.label != label;
        self.codons = codons;
        self.label = label;
        changed
    }
}

pub struct SequenceView {
    document: Rc<Document>,
    selected: Vec<Segment>,
    options: SequenceOptions,
    sites: Vec<SequenceRestrictionSite>,
    layout: SequenceLayout,
    custom_colors: sequence_layout::SequenceColorIndex,
    orfs: OrfSearch,
    selected_orf: Option<usize>,
    flanking: FlankingTracks,
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
        let custom_colors = sequence_layout::SequenceColorIndex::new(&document);
        let mut view = Self {
            document: Rc::new(document),
            selected,
            options: SequenceOptions::default(),
            sites: vec![],
            layout,
            custom_colors,
            orfs: OrfSearch::default(),
            selected_orf: None,
            flanking: FlankingTracks::default(),
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
    fn select_orf(&mut self, index: usize, cx: &mut Context<Self>) {
        let Some(orf) = self.orfs.orfs.get(index) else {
            return;
        };
        self.selected_orf = Some(index);
        self.scroll
            .scroll_to_item(orf.start / self.per_row, ScrollStrategy::Top);
        cx.emit(SequenceEvent::Select(
            orf.segments(self.document.sequence.len()),
        ));
        cx.notify();
    }
    fn navigate_orf(&mut self, forward: bool, cx: &mut Context<Self>) {
        let count = self.orfs.orfs.len();
        if count == 0 {
            return;
        }
        let index = match self.selected_orf {
            Some(i) => {
                if forward {
                    (i + 1) % count
                } else {
                    (i + count - 1) % count
                }
            }
            None => {
                if forward {
                    0
                } else {
                    count - 1
                }
            }
        };
        self.select_orf(index, cx);
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
        let flanks_changed = self
            .flanking
            .update(&document, &selected, options, document_changed);
        let layout_changed =
            document_changed || self.options != options || self.sites != sites || flanks_changed;
        let orfs_changed = document_changed
            || self.options.show_orfs != options.show_orfs
            || self.options.genetic_code != options.genetic_code
            || self.options.orf_min_codons != options.orf_min_codons
            || self.options.orf_alternative_starts != options.orf_alternative_starts
            || self.options.orf_include_partial != options.orf_include_partial;
        let selection_changed = self.selected != selected;
        if selection_changed {
            self.anchor = selected.first().map(|segment| segment.start);
            self.dragging = None;
        }
        if document_changed {
            self.layout = SequenceLayout::new(&document);
            self.custom_colors = sequence_layout::SequenceColorIndex::new(&document);
            self.document = Rc::new(document);
            self.anchor = None;
            self.dragging = None;
        }
        self.options = options;
        self.sites = sites;
        self.selected = selected;
        if orfs_changed {
            self.selected_orf = None;
            self.orfs = if options.show_orfs {
                GeneticCode::from_id(options.genetic_code)
                    .map(|genetic_code| {
                        find_orfs(
                            &self.document,
                            OrfOptions {
                                genetic_code,
                                min_codons: options.orf_min_codons,
                                alternative_starts: options.orf_alternative_starts,
                                include_partial: options.orf_include_partial,
                                ..OrfOptions::default()
                            },
                        )
                    })
                    .unwrap_or_default()
            } else {
                OrfSearch::default()
            };
        }
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
                for (boundary, top) in [
                    (site.top_cut, true),
                    (site.bottom_cut, false),
                    (site.second_top_cut, true),
                    (site.second_bottom_cut, false),
                ] {
                    if let Some(boundary) = boundary.filter(|boundary| *boundary <= n) {
                        let row = boundary.min(n.saturating_sub(1)) / self.per_row;
                        if let Some(values) = cuts.get_mut(row) {
                            values.push((boundary - row * self.per_row, top));
                        }
                    }
                }
            }
        }
        let mut orf_rows = vec![0_u8; rows];
        for orf in &self.orfs.orfs {
            let frame = orf.frame(n);
            let bit = usize::from(frame.unsigned_abs() - 1) + if frame < 0 { 3 } else { 0 };
            for part in orf.segments(n) {
                for row in orf_rows
                    .iter_mut()
                    .take((part.end - 1) / self.per_row + 1)
                    .skip(part.start / self.per_row)
                {
                    *row |= 1 << bit;
                }
            }
        }
        let mut flank_rows = vec![0_u8; rows];
        for flank in &self.flanking.codons {
            flank_rows[flank.codon.center / self.per_row] |= if flank.upstream { 1 } else { 2 };
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
                let annotations_y = top
                    + if self.paired() { 61. } else { 31. }
                    + 20.
                        * (usize::from(self.document.sequence_qualities.is_some())
                            + usize::from(
                                self.paired()
                                    && self.document.opposite_sequence_qualities.is_some(),
                            )) as f32;
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
                let mut orf_y = vec![];
                for (bit, frame) in [1, 2, 3, -1, -2, -3].into_iter().enumerate() {
                    if orf_rows[row] & (1 << bit) != 0 {
                        orf_y.push((frame, y));
                        y += 22.;
                    }
                }
                let mut flank_y = vec![];
                for (upstream, bit) in [(true, 1), (false, 2)] {
                    if flank_rows[row] & bit != 0 {
                        flank_y.push((upstream, y));
                        y += 22.;
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
                    orf_y,
                    flank_y,
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
            .text_color(
                if self
                    .selected
                    .iter()
                    .any(|range| range.start <= index && index < range.end)
                {
                    rgb(INK)
                } else if let Some(color) = self.custom_colors.color(index, opposite) {
                    rgb(color)
                } else if self.options.colored_bases {
                    base_color(base)
                } else {
                    rgb(INK)
                },
            )
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
        let opposite_span = sequence_layout::opposite_span(&self.document);
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
                            if !(opposite_span.start..opposite_span.end).contains(&(start + column))
                            {
                                continue;
                            }
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
                        if !forward
                            && (!paired
                                || start + column < opposite_span.start
                                || start + column > opposite_span.end)
                        {
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
            if paired && let Some(opposite) = sequence_layout::opposite_base(&self.document, index)
            {
                content = content.child(self.base_cell(index, opposite, top + 30., true, cx));
            }
        }
        if let Some(qualities) = &self.document.sequence_qualities {
            let y = top + if paired { 61. } else { 31. };
            content = content.child(
                div()
                    .absolute()
                    .left(px(20.))
                    .top(px(y))
                    .text_size(px(10.))
                    .text_color(rgb(MUTED))
                    .child("Q"),
            );
            for (index, &quality) in qualities.iter().enumerate().take(end).skip(start) {
                content = content.child(
                    div()
                        .absolute()
                        .left(px(LEFT + sequence_layout::base_x(index % self.per_row, cw)))
                        .top(px(y))
                        .w(px(cw))
                        .text_center()
                        .text_size(px(8.))
                        .text_color(rgb(if quality.is_some_and(|q| q < 20) {
                            0xf0bd76
                        } else {
                            MUTED
                        }))
                        .child(quality.map_or("—".into(), |q| q.to_string())),
                );
            }
        }
        if paired && self.document.opposite_sequence_qualities.is_some() {
            let y = top
                + 61.
                + if self.document.sequence_qualities.is_some() {
                    20.
                } else {
                    0.
                };
            content = content.child(
                div()
                    .absolute()
                    .left(px(10.))
                    .top(px(y))
                    .text_size(px(10.))
                    .text_color(rgb(MUTED))
                    .child("Q−"),
            );
            for index in start..end {
                if sequence_layout::opposite_base(&self.document, index).is_some() {
                    let quality = sequence_layout::opposite_quality(&self.document, index);
                    content = content.child(
                        div()
                            .absolute()
                            .left(px(LEFT + sequence_layout::base_x(index % self.per_row, cw)))
                            .top(px(y))
                            .w(px(cw))
                            .text_center()
                            .text_size(px(8.))
                            .text_color(rgb(if quality.is_some_and(|q| q < 20) {
                                0xf0bd76
                            } else {
                                MUTED
                            }))
                            .child(quality.map_or("—".into(), |q| q.to_string())),
                    );
                }
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
            let length = ((span.end - span.start) as f32 * cw).max(3.);
            let y = plan.lane_y[span.lane];
            let color = span
                .color_rgb
                .map(|color| rgb(color).into())
                .unwrap_or_else(|| feature_card::feature_color(feature, span.feature_index));
            let label_color = if span.style == popgene_core::FeatureSegmentStyle::Filled {
                annotation_label_color(color)
            } else {
                rgb(INK)
            };
            let style = span.style;
            let cleavages = span
                .cleavage_sites
                .iter()
                .map(|boundary| (*boundary - span.start) as f32 * cw)
                .collect::<Vec<_>>();
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
                            if style == popgene_core::FeatureSegmentStyle::Filled {
                                if let Ok(path) = p.build() {
                                    window.paint_path(path, color);
                                }
                            } else {
                                let mut line = PathBuilder::stroke(px(1.5));
                                line.move_to(point(x, y + h - px(2.)));
                                line.line_to(point(x + w, y + h - px(2.)));
                                if arrow_left {
                                    line.move_to(point(x + tip, y + h - px(6.)));
                                    line.line_to(point(x, y + h - px(2.)));
                                }
                                if arrow_right {
                                    line.move_to(point(x + w - tip, y + h - px(6.)));
                                    line.line_to(point(x + w, y + h - px(2.)));
                                }
                                if let Ok(path) = line.build() {
                                    window.paint_path(path, color);
                                }
                            }
                            for &boundary in &cleavages {
                                let mut tick = PathBuilder::stroke(px(2.));
                                tick.move_to(point(x + px(boundary), y - px(2.)));
                                tick.line_to(point(x + px(boundary), y + h + px(2.)));
                                if let Ok(path) = tick.build() {
                                    window.paint_path(path, label_color);
                                }
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
                        .child(
                            if length >= 32.
                                && (style != popgene_core::FeatureSegmentStyle::Gap
                                    || span.segment_name.is_some())
                            {
                                feature_card::truncate(
                                    span.segment_name.as_deref().unwrap_or(&feature.label),
                                    100,
                                )
                            } else {
                                String::new()
                            },
                        ),
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
            for codon in frame_codons_with_options(
                &self.document.sequence,
                *frame,
                start..end,
                GeneticCode::from_id(self.options.genetic_code).unwrap_or(GeneticCode::STANDARD),
                self.document.circular,
            ) {
                let left = codon.start.max(start) - start;
                let right = codon.end.min(end) - start;
                let selection = codon.selection();
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
                                cx.emit(SequenceEvent::Select(selection.clone()));
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
        for (frame, y) in &plan.orf_y {
            for (index, orf) in self
                .orfs
                .orfs
                .iter()
                .enumerate()
                .filter(|(_, orf)| orf.frame(self.document.sequence.len()) == *frame)
            {
                for (part_index, part) in orf
                    .segments(self.document.sequence.len())
                    .into_iter()
                    .enumerate()
                {
                    let left = part.start.max(start);
                    let right = part.end.min(end);
                    if left >= right {
                        continue;
                    }
                    let label = format!(
                        "ORF {frame:+} · {} aa{}",
                        orf.amino_acid_length,
                        if orf.complete { "" } else { " · open end" }
                    );
                    let details = format!(
                        "Predicted ORF · {} · code {} · {}{}",
                        if orf.reverse {
                            "reverse strand"
                        } else {
                            "forward strand"
                        },
                        orf.genetic_code,
                        if orf.complete {
                            "start to stop"
                        } else {
                            "no terminal stop within the boundary"
                        },
                        if orf.context_dependent_stop {
                            " · terminal stop is context-dependent"
                        } else {
                            ""
                        }
                    );
                    content = content.child(
                        div()
                            .absolute()
                            .left(px(LEFT + (left - start) as f32 * cw))
                            .top(px(*y))
                            .w(px((right - left) as f32 * cw))
                            .h(px(18.))
                            .overflow_hidden()
                            .child(
                                HoverCard::new(("orf", row * 20_000 + index * 2 + part_index))
                                    .content(move |_, _, _| {
                                        div()
                                            .p_3()
                                            .max_w(px(420.))
                                            .text_xs()
                                            .child(details.clone())
                                            .into_any_element()
                                    })
                                    .trigger(
                                        div()
                                            .size_full()
                                            .cursor_pointer()
                                            .bg(rgba(if orf.reverse {
                                                0x6f9ed866
                                            } else {
                                                0xddaa5066
                                            }))
                                            .text_color(rgb(INK))
                                            .text_size(px(10.))
                                            .px_1()
                                            .when(self.selected_orf == Some(index), |d| {
                                                d.border_1().border_color(rgb(ACCENT))
                                            })
                                            .on_mouse_down(
                                                MouseButton::Left,
                                                cx.listener(move |this, _, _, cx| {
                                                    this.select_orf(index, cx);
                                                    cx.stop_propagation();
                                                }),
                                            )
                                            .child(label),
                                    ),
                            ),
                    );
                }
            }
        }
        for &(upstream, y) in &plan.flank_y {
            content = content.child(
                div()
                    .absolute()
                    .left(px(3.))
                    .top(px(y))
                    .text_xs()
                    .text_color(rgb(MUTED))
                    .child(if upstream { "5′ flank" } else { "3′ flank" }),
            );
            for flank in self.flanking.codons.iter().filter(|flank| {
                flank.upstream == upstream && (start..end).contains(&flank.codon.center)
            }) {
                let codon = &flank.codon;
                let left = codon.start.max(start) - start;
                let right = codon.end.min(end) - start;
                let selection = codon.selection();
                let details = format!(
                    "Hypothetical {} CDS flank · bases {}, {}, {}",
                    if flank.upstream {
                        "upstream"
                    } else {
                        "downstream"
                    },
                    codon.positions[0] + 1,
                    codon.positions[1] + 1,
                    codon.positions[2] + 1
                );
                content = content.child(
                    div()
                        .absolute()
                        .left(px(LEFT + left as f32 * cw))
                        .top(px(y))
                        .w(px((right - left) as f32 * cw))
                        .h(px(18.))
                        .bg(rgba(0xc2a6d433))
                        .child(
                            HoverCard::new(("cds-flank", codon.center * 2 + usize::from(upstream)))
                                .content(move |_, _, _| {
                                    div()
                                        .p_3()
                                        .text_xs()
                                        .child(details.clone())
                                        .into_any_element()
                                })
                                .trigger(
                                    div()
                                        .size_full()
                                        .cursor_pointer()
                                        .on_mouse_down(
                                            MouseButton::Left,
                                            cx.listener(move |_, _, _, cx| {
                                                cx.emit(SequenceEvent::Select(selection.clone()));
                                                cx.stop_propagation();
                                            }),
                                        )
                                        .child(
                                            div()
                                                .ml(px((codon.center - start - left) as f32 * cw))
                                                .w(px(cw))
                                                .text_center()
                                                .text_size(px(11.))
                                                .text_color(rgb(0xd5b6ed))
                                                .child((codon.amino_acid as char).to_string()),
                                        ),
                                ),
                        ),
                );
            }
        }
        for glyph in &annotations.translations {
            let y = plan.lane_y[glyph.lane] + 20.;
            let details = format!(
                "Calculated CDS · genetic code {} · bases {}{}{}",
                glyph.genetic_code,
                glyph
                    .codon_positions
                    .iter()
                    .take(glyph.nucleotide_count)
                    .map(|p| (p + 1).to_string())
                    .collect::<Vec<_>>()
                    .join(", "),
                if glyph.exception {
                    " · recorded translation exception applied"
                } else {
                    ""
                },
                if glyph.nucleotide_count < 3 {
                    " · partial stop completed by polyadenylation"
                } else if glyph.split_codon {
                    " · joined codon; central base highlighted"
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
        let body = div()
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
            .child(Scrollbar::vertical(&self.scroll));
        div()
            .v_flex()
            .size_full()
            .min_h_0()
            .when(self.options.show_orfs, |d| {
                d.child(
                    div()
                        .h_flex()
                        .flex_shrink_0()
                        .px_2()
                        .py_1()
                        .gap_2()
                        .items_center()
                        .text_xs()
                        .child(format!(
                            "{} ORFs{}",
                            self.orfs.orfs.len(),
                            if self.orfs.truncated {
                                " (limit reached)"
                            } else {
                                ""
                            }
                        ))
                        .child(
                            Button::new("previous-orf")
                                .label("Previous ORF")
                                .small()
                                .disabled(self.orfs.orfs.is_empty())
                                .on_click(
                                    cx.listener(|this, _, _, cx| this.navigate_orf(false, cx)),
                                ),
                        )
                        .child(
                            Button::new("next-orf")
                                .label("Next ORF")
                                .small()
                                .disabled(self.orfs.orfs.is_empty())
                                .on_click(
                                    cx.listener(|this, _, _, cx| this.navigate_orf(true, cx)),
                                ),
                        )
                        .child(
                            Button::new("copy-orf")
                                .label("Copy protein")
                                .small()
                                .disabled(self.selected_orf.is_none())
                                .on_click(cx.listener(|this, _, _, cx| {
                                    if let Some(orf) =
                                        this.selected_orf.and_then(|i| this.orfs.orfs.get(i))
                                    {
                                        cx.write_to_clipboard(ClipboardItem::new_string(
                                            orf.protein(&this.document.sequence),
                                        ));
                                    }
                                })),
                        ),
                )
            })
            .when_some(
                sequence_layout::duplex_display_notice(&self.document),
                |d, notice| {
                    d.child(
                        div()
                            .flex_shrink_0()
                            .px_2()
                            .py_1()
                            .text_xs()
                            .text_color(rgb(MUTED))
                            .child(notice),
                    )
                },
            )
            .child(body.flex_1())
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
        let code = GeneticCode::from_id(self.options.genetic_code);
        let show_orfs = self.options.show_orfs && self.document.molecule != Molecule::Protein;
        let show_flanking =
            self.options.show_flanking && self.document.molecule != Molecule::Protein;
        div().v_flex().w_full().min_w_0()
            .when(show_cds || show_frames || show_orfs || show_flanking, |d| {
                d.child(div().id("translation-status").v_flex().gap_2().py_2().text_xs().text_color(rgb(MUTED))
                    .when(show_cds, |d| d.child("Calculated CDS · recorded phase, genetic code and exact translation exceptions"))
                    .when(show_frames || show_orfs, |d| d.child(format!("Frames / ORFs · code {} ({}){}",self.options.genetic_code,code.map(|c|c.name()).unwrap_or("invalid code"),if self.document.circular {" · circular origin crossings enabled"} else {" · linear boundaries"})))
                    .when(code.is_some_and(|c|c.has_context_dependent_stops()),|d|d.child("Context-dependent stops: frames show internal residues; CDS ends and ORF candidates use terminal stop flags."))
                    .child("* stop · X unresolved codon · reverse frames start at the last base")
                    .when(show_cds && self.layout.omitted_translations > 0, |d| {
                        d.child(div().text_color(rgb(0xe0b472)).child(format!("{} CDS translations omitted:", self.layout.omitted_translations)))
                            .children(self.layout.translation_notices.iter().map(|notice| div().child(notice.clone())))
                    })
                    .when(show_flanking,|d|d.child(format!("Selected CDS flanks · {} · up to {} codons each side or nearest stop (hypothetical)",self.flanking.label.as_deref().unwrap_or("select a CDS annotation"),self.options.flanking_codons.min(300))))
                    .when(show_orfs,|d|d.child(format!("ORF predictions · ≥{} aa · {} · {} · longest start per stop interval",self.options.orf_min_codons,if self.options.orf_alternative_starts {"code-specific starts"}else{"ATG starts"},if self.options.orf_include_partial {"open ends included"}else{"complete only"})))
                )
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
pub(crate) fn annotation_label_color(background: Hsla) -> Rgba {
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
    ordered_calls: Arc<[usize]>,
}
pub enum TraceEvent {
    SelectCall(usize),
}
impl EventEmitter<TraceEvent> for TraceView {}
const TRACE_CHUNK: usize = 256;
// Leave room for centered bases, quality labels, call numbers and hit targets
// at both physical ends. Interior chunks add no gap to the signal axis.
const TRACE_EDGE_GUTTER: f32 = 40.;
const TRACE_CALL_HIT_RADIUS: f32 = 8.;
#[derive(Clone)]
struct TraceChunkGeometry {
    start: usize,
    end: usize,
    signal_samples: std::ops::Range<usize>,
    sample_width: f32,
    left_gutter: f32,
    width: f32,
}
impl TraceChunkGeometry {
    fn new(visible: std::ops::Range<usize>, chunk: usize, sample_width: f32) -> Self {
        let start = visible.start + chunk * TRACE_CHUNK;
        let end = (start + TRACE_CHUNK).min(visible.end);
        let first = start == visible.start;
        let last = end == visible.end;
        let left_gutter = if first { TRACE_EDGE_GUTTER } else { 0. };
        // Samples are points, so the final chunk has no extra interval after
        // its last point. This also gives a one-sample trace a finite canvas.
        let intervals = end - start - usize::from(last);
        Self {
            start,
            end,
            signal_samples: start.saturating_sub(1).max(visible.start)..(end + 1).min(visible.end),
            sample_width,
            left_gutter,
            width: left_gutter
                + intervals as f32 * sample_width
                + if last { TRACE_EDGE_GUTTER } else { 0. },
        }
    }
    fn sample_x(&self, sample: usize) -> f32 {
        let offset = if sample >= self.start {
            (sample - self.start) as f32
        } else {
            -((self.start - sample) as f32)
        };
        self.left_gutter + offset * self.sample_width
    }
    fn call_at_x(&self, x: f32, peaks: &[usize], ordered_calls: &[usize]) -> Option<usize> {
        if ordered_calls.is_empty() {
            return None;
        }
        let sample =
            self.start as f64 + f64::from(x - self.left_gutter) / f64::from(self.sample_width);
        let after = ordered_calls.partition_point(|&index| (peaks[index] as f64) < sample);
        let peak = match (after.checked_sub(1), ordered_calls.get(after)) {
            (Some(before), Some(&index)) => {
                let left = peaks[ordered_calls[before]];
                let right = peaks[index];
                // At an exact midpoint prefer the earlier signal position.
                if sample - left as f64 <= right as f64 - sample {
                    left
                } else {
                    right
                }
            }
            (Some(before), None) => peaks[ordered_calls[before]],
            (None, Some(&index)) => peaks[index],
            (None, None) => return None,
        };
        if (self.sample_x(peak) - x).abs() > TRACE_CALL_HIT_RADIUS {
            return None;
        }
        // Coincident peaks cannot be separated spatially. Pick their first
        // call index consistently; the First/Last call fields address others.
        let first = ordered_calls.partition_point(|&index| peaks[index] < peak);
        Some(ordered_calls[first])
    }
}
impl TraceView {
    pub fn new(document: Document, scale: f32, _: &mut Context<Self>) -> Self {
        let mut view = Self {
            trace: document.trace.map(Arc::new),
            sequence: document.sequence,
            scale: scale.clamp(0.5, 5.),
            max_signal: 1.,
            sizes: Rc::new(vec![]),
            scroll: VirtualListScrollHandle::new(),
            ordered_calls: Arc::from([]),
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
    fn visible_samples(&self) -> std::ops::Range<usize> {
        self.trace
            .as_ref()
            .and_then(|trace| trace.sample_range.as_ref())
            .map_or(0..self.samples(), |range| range.start..range.end)
    }
    fn sample_width(&self) -> f32 {
        (800. / self.visible_samples().len().max(1) as f32).max(2.) * self.scale
    }
    fn remeasure(&mut self) {
        let mut ordered: Vec<usize> = (0..self.sequence.len()).collect();
        if let Some(trace) = &self.trace {
            ordered.sort_by_key(|&i| (trace.peaks[i], i));
        }
        self.ordered_calls = ordered.into();
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
        let visible = self.visible_samples();
        let width = self.sample_width();
        self.sizes = Rc::new(
            (0..visible.len().div_ceil(TRACE_CHUNK))
                .map(|i| {
                    let geometry = TraceChunkGeometry::new(visible.clone(), i, width);
                    size(px(geometry.width), px(420.))
                })
                .collect(),
        );
    }
    fn render_chunk(&self, chunk: usize, cx: &mut Context<Self>) -> impl IntoElement {
        let trace = self.trace.clone();
        let sequence = self.sequence.clone();
        let ordered_calls = self.ordered_calls.clone();
        let max_signal = self.max_signal;
        let geometry = TraceChunkGeometry::new(self.visible_samples(), chunk, self.sample_width());
        let start = geometry.start;
        let end = geometry.end;
        let hit_bounds = Rc::new(RefCell::new(None::<Bounds<Pixels>>));
        let canvas_bounds = hit_bounds.clone();
        let hit_geometry = geometry.clone();
        let hit_trace = trace.clone();
        let hit_calls = ordered_calls.clone();
        let body = div()
            .relative()
            .w(px(geometry.width))
            .h(px(420.))
            .flex_shrink_0()
            .child(
                div()
                    .id(("trace-call-band", chunk))
                    .absolute()
                    .left_0()
                    .top(px(15.))
                    .w_full()
                    .h(px(40.))
                    .cursor_pointer()
                    .on_mouse_down(
                        MouseButton::Left,
                        cx.listener(move |_, event: &MouseDownEvent, _, cx| {
                            let bounds = *hit_bounds.borrow();
                            if let (Some(bounds), Some(trace)) = (bounds, &hit_trace)
                                && let Some(index) = hit_geometry.call_at_x(
                                    (event.position.x - bounds.left()).as_f32(),
                                    &trace.peaks,
                                    &hit_calls,
                                )
                            {
                                cx.emit(TraceEvent::SelectCall(index));
                            }
                        }),
                    ),
            );
        body.child(
            canvas(
                move |bounds, _, _| *canvas_bounds.borrow_mut() = Some(bounds),
                move |bounds, _, window, cx| {
                    let Some(trace) = &trace else {
                        return;
                    };
                    let x = |i: usize| bounds.left() + px(geometry.sample_x(i));
                    let baseline = bounds.bottom() - px(38.);
                    let amplitude = bounds.size.height - px(126.);
                    let mut grid = PathBuilder::stroke(px(1.));
                    grid.move_to(point(bounds.left(), baseline));
                    grid.line_to(point(bounds.right(), baseline));
                    if let Ok(path) = grid.build() {
                        window.paint_path(path, rgb(0x334151));
                    }
                    // One overlapping point makes adjacent chunk paths continuous.
                    let line_start = geometry.signal_samples.start;
                    let line_end = geometry.signal_samples.end;
                    for (channel, base) in trace.channels.iter().zip(trace.channel_order.bytes()) {
                        let mut path = PathBuilder::stroke(px(1.5));
                        for (i, sample) in
                            channel.iter().enumerate().take(line_end).skip(line_start)
                        {
                            let signal_point =
                                point(x(i), baseline - amplitude * (*sample as f32 / max_signal));
                            if line_end - line_start == 1 {
                                // A single recorded sample has no connecting
                                // segment; retain it as a small point mark.
                                path.move_to(point(signal_point.x - px(1.), signal_point.y));
                                path.line_to(point(signal_point.x + px(1.), signal_point.y));
                            } else if i == line_start {
                                path.move_to(signal_point);
                            } else {
                                path.line_to(signal_point);
                            }
                        }
                        if let Ok(path) = path.build() {
                            window.paint_path(path, base_color(base));
                        }
                    }
                    // Signal-order index preserves calls whose original peak
                    // positions cross while avoiding a full scan per canvas.
                    let first = ordered_calls.partition_point(|&i| trace.peaks[i] < start);
                    for &i in &ordered_calls[first..] {
                        let peak = trace.peaks[i];
                        if peak >= end {
                            break;
                        }
                        let Some(&base) = sequence.as_bytes().get(i) else {
                            continue;
                        };
                        paint_text(
                            (base as char).to_string(),
                            point(x(peak), bounds.top() + px(18.)),
                            16.,
                            base_color(base),
                            window,
                            cx,
                        );
                        if let Some(quality) = trace.qualities.get(i) {
                            paint_text(
                                if trace.manual_calls.get(i) == Some(&true) {
                                    "edit".into()
                                } else if !trace.qualities_known {
                                    "—".into()
                                } else {
                                    format!("Q{quality}")
                                },
                                point(x(peak), bounds.top() + px(42.)),
                                10.,
                                rgb(0x8c9bad),
                                window,
                                cx,
                            );
                        }
                        paint_text(
                            (i + 1).to_string(),
                            point(x(peak), baseline + px(10.)),
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
                "{} base calls  ·  {} signal samples  ·  call family {}{}{}",
                self.sequence.len(),
                self.samples(),
                trace.calls_version,
                if trace.edited { " · edited" } else { "" },
                trace
                    .sample_range
                    .as_ref()
                    .map_or(String::new(), |r| format!(
                        " · sample crop {}–{}",
                        r.start,
                        r.end - 1
                    ))
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
                    |this, range, _, cx| {
                        range
                            .map(|chunk| this.render_chunk(chunk, cx).into_any_element())
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
    center: Point<Pixels>,
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
    let origin = point(center.x - line.width / 2., center.y);
    let _ = line.paint(origin, px(font_size), TextAlign::Left, None, window, cx);
}

#[cfg(test)]
mod trace_geometry_tests {
    use super::{TRACE_CHUNK, TRACE_EDGE_GUTTER, TraceChunkGeometry};

    #[test]
    fn endpoint_calls_and_hit_targets_fit_at_both_zoom_limits() {
        for sample_width in [50., 500.] {
            let geometry = TraceChunkGeometry::new(0..8, 0, sample_width);
            assert_eq!(geometry.sample_x(0), TRACE_EDGE_GUTTER);
            assert_eq!(geometry.width - geometry.sample_x(7), TRACE_EDGE_GUTTER);
            for sample in [0, 7] {
                let x = geometry.sample_x(sample);
                assert!(x - 8. >= 0. && x + 8. <= geometry.width);
            }
        }
    }

    #[test]
    fn one_sample_has_centered_gutters_without_a_zoom_sized_empty_tail() {
        for start in [0, 900] {
            for sample_width in [400., 4000.] {
                let geometry = TraceChunkGeometry::new(start..start + 1, 0, sample_width);
                assert_eq!(geometry.width, 2. * TRACE_EDGE_GUTTER);
                assert_eq!(geometry.sample_x(start), TRACE_EDGE_GUTTER);
                assert_eq!(geometry.signal_samples, start..start + 1);
            }
        }
    }

    #[test]
    fn virtual_chunks_share_one_continuous_sample_axis() {
        let visible = 50..50 + 2 * TRACE_CHUNK + 1;
        for sample_width in [1., 10.] {
            let chunks: Vec<_> = (0..3)
                .map(|chunk| TraceChunkGeometry::new(visible.clone(), chunk, sample_width))
                .collect();
            let mut left = 0.;
            for geometry in &chunks {
                for sample in geometry.signal_samples.clone() {
                    assert_eq!(
                        left + geometry.sample_x(sample),
                        TRACE_EDGE_GUTTER + (sample - visible.start) as f32 * sample_width
                    );
                }
                left += geometry.width;
            }
            assert_eq!(
                left,
                2. * TRACE_EDGE_GUTTER + (visible.len() - 1) as f32 * sample_width
            );
            assert_eq!(chunks[0].signal_samples.start, visible.start);
            assert_eq!(chunks[2].signal_samples.end, visible.end);
        }
    }

    #[test]
    fn dense_peak_clicks_choose_nearest_call_instead_of_later_overlapping_target() {
        let peaks = [0, 6, 11, 19, 27, 35, 100];
        let ordered = [0, 1, 2, 3, 4, 5, 6];
        // At minimum zoom a long trace has one pixel per signal sample.
        let geometry = TraceChunkGeometry::new(0..4000, 0, 1.);
        for (index, &peak) in peaks.iter().enumerate() {
            assert_eq!(
                geometry.call_at_x(geometry.sample_x(peak), &peaks, &ordered),
                Some(index)
            );
        }
        assert_eq!(
            geometry.call_at_x(geometry.sample_x(8) + 0.5, &peaks, &ordered),
            Some(1)
        );
        assert_eq!(
            geometry.call_at_x(geometry.sample_x(8) + 0.6, &peaks, &ordered),
            Some(2)
        );
        assert_eq!(
            geometry.call_at_x(geometry.sample_x(60), &peaks, &ordered),
            None
        );
    }

    #[test]
    fn neighboring_chunks_resolve_hits_against_the_same_peak_order() {
        let peaks = [255, 256, 263];
        let ordered = [0, 1, 2];
        let left = TraceChunkGeometry::new(0..600, 0, 1.);
        let right = TraceChunkGeometry::new(0..600, 1, 1.);
        for sample in [254.9, 255., 255.5, 255.6, 256., 260.] {
            let left_x = TRACE_EDGE_GUTTER + sample;
            let right_x = left_x - left.width;
            assert_eq!(
                left.call_at_x(left_x, &peaks, &ordered),
                right.call_at_x(right_x, &peaks, &ordered)
            );
        }
        assert_eq!(right.call_at_x(0., &peaks, &ordered), Some(1));
        // The nearest call may belong to the adjacent chunk.
        assert_eq!(left.call_at_x(left.width - 0.1, &peaks, &ordered), Some(1));
    }

    #[test]
    fn duplicate_and_crossed_peaks_choose_lowest_call_index_at_that_position() {
        let peaks = [256, 255, 256, 300];
        let ordered = [1, 0, 2, 3];
        let geometry = TraceChunkGeometry::new(0..600, 1, 1.);
        for x in [0., 0.1, 1., 7.9] {
            assert_eq!(geometry.call_at_x(x, &peaks, &ordered), Some(0));
        }
        assert_eq!(geometry.call_at_x(-1., &peaks, &ordered), Some(1));
        assert_eq!(geometry.call_at_x(0., &[], &[]), None);
    }
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

#[cfg(test)]
mod flanking_selection_tests {
    use super::{FlankingTracks, SequenceOptions};
    use popgene_core::{Document, Segment};

    fn fixture() -> Document {
        popgene_formats::import(
            "flanking.gb",
            include_bytes!("../../../validation/fixtures/translation/flanking.gb"),
        )
        .unwrap()
    }

    fn options() -> SequenceOptions {
        SequenceOptions {
            show_flanking: true,
            ..SequenceOptions::default()
        }
    }

    #[test]
    fn mouse_selection_echo_clears_previous_cds_flanks_without_remeasuring_every_drag_step() {
        let doc = fixture();
        let mut tracks = FlankingTracks::default();
        let cds = doc.features[0].parts.clone();
        assert!(tracks.update(&doc, &cds, options(), false));
        assert_eq!(tracks.label.as_deref(), Some("Select_for_flanks"));
        assert_eq!(
            tracks
                .codons
                .iter()
                .map(|flank| (flank.codon.center, flank.codon.amino_acid))
                .collect::<Vec<_>>(),
            vec![(4, b'P'), (1, b'*'), (16, b'G'), (19, b'*')]
        );

        // The mouse handler already changed SequenceView.selected before the
        // model echoes base 1. The rendered CDS key must still invalidate.
        let mouse_selection = [Segment { start: 0, end: 1 }];
        assert!(tracks.update(&doc, &mouse_selection, options(), false));
        assert!(tracks.codons.is_empty());
        assert!(tracks.label.is_none());
        for end in 1..=6 {
            assert!(!tracks.update(&doc, &[Segment { start: 0, end }], options(), false));
        }

        assert!(tracks.update(&doc, &cds, options(), false));
        let codons = tracks.codons.as_ptr();
        assert!(!tracks.update(&doc, &cds, options(), false));
        assert_eq!(tracks.codons.as_ptr(), codons);
    }

    #[test]
    fn flanking_cache_refreshes_for_options_and_edits_with_the_same_selected_cds() {
        let mut doc = fixture();
        let cds = doc.features[0].parts.clone();
        let mut tracks = FlankingTracks::default();
        assert!(tracks.update(&doc, &cds, options(), false));
        let shorter = SequenceOptions {
            flanking_codons: 1,
            ..options()
        };
        assert!(tracks.update(&doc, &cds, shorter, false));
        assert_eq!(tracks.codons.len(), 2);
        doc.replace(3, 6, "GGC").unwrap();
        assert!(tracks.update(&doc, &cds, shorter, true));
        assert_eq!(tracks.codons[0].codon.amino_acid, b'G');
        doc.features[0].label = "Renamed CDS".into();
        assert!(tracks.update(&doc, &cds, shorter, true));
        assert_eq!(tracks.label.as_deref(), Some("Renamed CDS"));
        assert!(tracks.update(&doc, &cds, SequenceOptions::default(), false));
        assert!(tracks.codons.is_empty());
        assert!(tracks.label.is_none());
    }

    #[test]
    fn circular_multi_part_selection_retains_its_biological_order() {
        let mut doc = fixture();
        doc.circular = true;
        let parts = vec![Segment { start: 15, end: 21 }, Segment { start: 0, end: 3 }];
        doc.features[0].parts = parts.clone();
        doc.features[0].location = None;
        let mut tracks = FlankingTracks::default();
        assert!(tracks.update(&doc, &parts, options(), true));
        let codons = tracks.codons.as_ptr();
        assert!(!tracks.update(&doc, &parts, options(), false));
        assert_eq!(tracks.codons.as_ptr(), codons);
        let reversed_order = parts.iter().rev().cloned().collect::<Vec<_>>();
        assert!(tracks.update(&doc, &reversed_order, options(), false));
        assert!(tracks.label.is_none());
        assert_eq!(doc.features[0].parts, parts);
    }
}
