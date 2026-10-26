//! Physical duplex inspection. The union axis never invents top-strand residues.
use crate::{
    feature_card,
    sequence_layout::{SequenceColorIndex, SequenceLayout},
    sequence_view::annotation_label_color,
};
use gpui_kit::component::{
    Disableable, Sizable, StyledExt, VirtualListScrollHandle,
    button::Button,
    input::{Input, InputState, Textarea, TextareaState},
    scroll::Scrollbar,
    v_virtual_list,
};
use gpui_kit::{prelude::*, *};
use popgene_core::{
    Document, Feature, FeatureSegmentStyle, Molecule, Segment, metadata::DuplexGeometry,
};
use std::rc::Rc;

const COLUMNS: usize = 60;
const CELL: f32 = 10.;
const LEFT: f32 = 95.;
const LANES: usize = 6;

pub enum DuplexEvent {
    Close,
    SelectTop(Vec<Segment>),
}

/// Annotation indexing uses a coordinate reference, while `base` reads only
/// the two actual stored strands. The reference is never rendered or copied.
struct DuplexScene {
    source: Document,
    axis: Document,
    geometry: DuplexGeometry,
    layout: SequenceLayout,
    colors: SequenceColorIndex,
}
impl DuplexScene {
    fn new(source: Document) -> Result<Self, String> {
        source.validate().map_err(|e| e.to_string())?;
        let ends = source
            .duplex_ends
            .as_ref()
            .ok_or("This document has no explicit physical duplex ends")?;
        let geometry = ends.geometry(source.sequence.len());
        let mut axis = Document::new(
            &source.name,
            &ends
                .reference_sequence(&source.sequence)
                .map_err(|e| e.to_string())?,
        )
        .map_err(|e| e.to_string())?;
        for original in &source.features {
            let mut feature = original.clone();
            let offset = geometry.top.start;
            if offset != 0 {
                feature.project_display_segments(|part| {
                    vec![Segment {
                        start: part.start + offset,
                        end: part.end + offset,
                    }]
                });
                feature.location = feature
                    .location
                    .as_ref()
                    .map(|location| location.replaced(0, 0, offset))
                    .transpose()
                    .map_err(|e| e.to_string())?
                    .flatten();
                for part in &mut feature.parts {
                    part.start += offset;
                    part.end += offset;
                }
                popgene_core::translation::project_translation_exceptions(&mut feature, |p| {
                    Some(p + offset)
                });
            }
            axis.features.push(feature);
        }
        axis.features
            .extend(source.duplex_annotations.features.iter().cloned());
        axis.sequence_colors = source
            .sequence_colors
            .iter()
            .cloned()
            .map(|mut color| {
                color.range.start += geometry.top.start;
                color.range.end += geometry.top.start;
                color
            })
            .chain(source.duplex_annotations.sequence_colors.iter().cloned())
            .collect();
        let layout = SequenceLayout::annotations_only(&axis);
        let colors = SequenceColorIndex::new(&axis);
        Ok(Self {
            source,
            axis,
            geometry,
            layout,
            colors,
        })
    }
    /// Physical nucleotide and zero-based index in its stored 5′→3′ strand.
    fn base(&self, position: usize, reverse: bool) -> Option<(u8, usize)> {
        if reverse {
            let range = &self.geometry.bottom;
            if position < range.start || position >= range.end {
                return None;
            }
            let index = range.end - 1 - position;
            Some((
                self.source
                    .duplex_ends
                    .as_ref()?
                    .opposite_sequence
                    .as_bytes()[index],
                index,
            ))
        } else {
            let index = position.checked_sub(self.geometry.top.start)?;
            self.source
                .sequence
                .as_bytes()
                .get(index)
                .copied()
                .map(|base| (base, index))
        }
    }
    fn quality(&self, position: usize, reverse: bool) -> Option<u8> {
        let index = self.base(position, reverse)?.1;
        if reverse {
            self.source
                .opposite_sequence_qualities
                .as_ref()?
                .get(index)
                .copied()
                .flatten()
        } else {
            self.source
                .sequence_qualities
                .as_ref()?
                .get(index)
                .copied()
                .flatten()
        }
    }
    fn original_feature(&self, index: usize) -> &Feature {
        if index < self.source.features.len() {
            &self.source.features[index]
        } else {
            &self.source.duplex_annotations.features[index - self.source.features.len()]
        }
    }
    fn selected_top(&self, index: usize) -> Option<Vec<Segment>> {
        let feature = self.axis.features.get(index)?;
        if feature
            .location
            .as_ref()
            .is_some_and(|location| location.uniform_projection().is_none())
        {
            return None;
        }
        let parts = &feature.parts;
        if parts.is_empty()
            || parts
                .iter()
                .any(|p| p.start < self.geometry.top.start || p.end > self.geometry.top.end)
        {
            return None;
        }
        Some(
            parts
                .iter()
                .map(|p| Segment {
                    start: p.start - self.geometry.top.start,
                    end: p.end - self.geometry.top.start,
                })
                .collect(),
        )
    }
    fn details(&self, index: usize) -> String {
        let feature = self.original_feature(index);
        format!(
            "Stored annotation coordinates: {} (1-based inclusive).\nFull duplex axis: {}\n\n{}",
            if index < self.source.features.len() {
                "top strand 5′→3′"
            } else {
                "full duplex extent, left→right"
            },
            self.axis.features[index].location_label(),
            feature_card::detail_text(feature, Molecule::Dna)
        )
    }
    fn row_height(&self, row: usize) -> f32 {
        let start = row * COLUMNS;
        let summary = self.layout.row_summary(
            start,
            (start + COLUMNS).min(self.geometry.extent),
            LANES,
            true,
            true,
            false,
        );
        155. + summary.lanes as f32 * 25. + if summary.omitted > 0 { 20. } else { 0. }
    }
}

pub struct DuplexView {
    scene: Option<DuplexScene>,
    selected: Option<usize>,
    details: Entity<TextareaState>,
    position: Entity<InputState>,
    sizes: Rc<Vec<Size<Pixels>>>,
    scroll: VirtualListScrollHandle,
    error: String,
}
impl EventEmitter<DuplexEvent> for DuplexView {}
impl DuplexView {
    pub fn new(document: Document, window: &mut Window, cx: &mut Context<Self>) -> Self {
        let (scene, error) = match DuplexScene::new(document) {
            Ok(scene) => (Some(scene), String::new()),
            Err(error) => (None, error),
        };
        let sizes = Rc::new(
            scene
                .as_ref()
                .map(|scene| {
                    (0..scene.geometry.extent.div_ceil(COLUMNS))
                        .map(|row| size(px(0.), px(scene.row_height(row))))
                        .collect()
                })
                .unwrap_or_default(),
        );
        let details = cx.new(|cx| {
            TextareaState::new(window, cx)
                .rows(5)
                .default_value("Select an annotation to inspect its complete recorded data.")
        });
        Self {
            scene,
            selected: None,
            details,
            position: cx.new(|cx| InputState::new(window, cx).default_value("1")),
            sizes,
            scroll: VirtualListScrollHandle::new(),
            error,
        }
    }
    fn select(&mut self, index: usize, window: &mut Window, cx: &mut Context<Self>) {
        let Some(scene) = &self.scene else {
            return;
        };
        let Some(feature) = scene.axis.features.get(index) else {
            return;
        };
        self.selected = Some(index);
        let text = scene.details(index);
        self.details
            .update(cx, |s, cx| s.set_value(text, window, cx));
        if let Some(start) = feature.parts.first().map(|p| p.start).or_else(|| {
            feature
                .location
                .as_ref()
                .and_then(|l| l.local_boundaries().first().copied())
        }) {
            self.scroll.scroll_to_item(
                start.min(scene.geometry.extent - 1) / COLUMNS,
                ScrollStrategy::Center,
            );
        }
        cx.notify();
    }
    fn navigate(&mut self, forward: bool, window: &mut Window, cx: &mut Context<Self>) {
        let count = self.scene.as_ref().map_or(0, |s| s.axis.features.len());
        if count == 0 {
            return;
        }
        let index = self
            .selected
            .map_or(if forward { 0 } else { count - 1 }, |index| {
                if forward {
                    (index + 1) % count
                } else {
                    (index + count - 1) % count
                }
            });
        self.select(index, window, cx);
    }
    fn jump(&mut self, cx: &mut Context<Self>) {
        let extent = self.scene.as_ref().map_or(0, |s| s.geometry.extent);
        match self.position.read(cx).value().trim().parse::<usize>() {
            Ok(position) if position > 0 && position <= extent => {
                self.scroll
                    .scroll_to_item((position - 1) / COLUMNS, ScrollStrategy::Center);
                self.error.clear();
            }
            _ => self.error = format!("Enter a full-duplex coordinate from 1 to {extent}"),
        }
        cx.notify();
    }
    fn row(&self, row: usize, cx: &mut Context<Self>) -> AnyElement {
        let scene = self.scene.as_ref().unwrap();
        let start = row * COLUMNS;
        let end = (start + COLUMNS).min(scene.geometry.extent);
        let layout = scene.layout.row(start, end, LANES, true, true, false);
        let mut output = div()
            .relative()
            .w_full()
            .h(px(scene.row_height(row)))
            .font_family("monospace")
            .text_size(px(13.))
            .border_b_1()
            .border_color(rgb(0x344354));
        output = output.child(
            div()
                .absolute()
                .left(px(8.))
                .top(px(5.))
                .text_xs()
                .text_color(rgb(0x95a6ba))
                .child(format!("Axis {}–{end}", start + 1)),
        );
        for (reverse, y) in [(false, 34.), (true, 90.)] {
            let strand = if reverse {
                &scene.geometry.bottom
            } else {
                &scene.geometry.top
            };
            let actual_start = start.max(strand.start);
            let actual_end = end.min(strand.end);
            let label = if actual_start < actual_end {
                let first = scene.base(actual_start, reverse).unwrap().1 + 1;
                let last = scene.base(actual_end - 1, reverse).unwrap().1 + 1;
                format!("{} {first}–{last}", if reverse { "3′ ←" } else { "5′ →" })
            } else {
                if reverse {
                    "Bottom absent"
                } else {
                    "Top absent"
                }
                .into()
            };
            output = output.child(
                div()
                    .absolute()
                    .left(px(4.))
                    .top(px(y + 3.))
                    .w(px(LEFT - 8.))
                    .text_size(px(10.))
                    .text_color(rgb(0x95a6ba))
                    .child(label),
            );
            for position in start..end {
                let x = LEFT + (position - start) as f32 * CELL;
                if let Some((base, index)) = scene.base(position, reverse) {
                    let color = scene.colors.color(position, reverse).unwrap_or(0xdce5ed);
                    let unpaired = scene.base(position, !reverse).is_none();
                    let quality = scene.quality(position, reverse);
                    output = output.child(div().id((if reverse {"physical-bottom"} else {"physical-top"}, position))
                        .absolute().left(px(x)).top(px(y)).w(px(CELL)).h(px(22.)).text_center().text_color(rgb(color))
                        .when(unpaired, |d| d.bg(rgba(0xb17f4033)).border_b_1().border_color(rgb(0xbda261)))
                        .cursor_pointer().on_mouse_down(MouseButton::Left, cx.listener(move |this, _, window, cx| {
                            this.selected = None;
                            this.details.update(cx, |s, cx| s.set_value(format!("{} strand, physical 5′→3′ position {}\nFull duplex axis position {}\nBase: {}{}\nQuality: {}", if reverse { "Bottom" } else { "Top" }, index + 1, position + 1, base as char, if unpaired { "\nUnpaired protruding residue" } else { "\nPaired residue" }, quality.map_or("unknown".into(), |q| q.to_string())), window, cx));
                            cx.notify();
                        })).child((base as char).to_string()));
                }
            }
            let known = if reverse {
                scene.source.opposite_sequence_qualities.is_some()
            } else {
                scene.source.sequence_qualities.is_some()
            };
            output = output.child(
                div()
                    .absolute()
                    .left(px(8.))
                    .top(px(y + 22.))
                    .text_size(px(9.))
                    .text_color(rgb(0x95a6ba))
                    .child(if known { "Quality" } else { "Q unknown" }),
            );
            for position in actual_start..actual_end {
                if scene.base(position, reverse).is_some() {
                    let quality = scene.quality(position, reverse);
                    output = output.child(
                        div()
                            .absolute()
                            .left(px(LEFT + (position - start) as f32 * CELL))
                            .top(px(y + 22.))
                            .w(px(CELL))
                            .text_center()
                            .text_size(px(8.))
                            .text_color(rgb(if quality.is_some_and(|q| q < 20) {
                                0xf0bd76
                            } else {
                                0x95a6ba
                            }))
                            .child(quality.map_or("—".into(), |q| q.to_string())),
                    );
                }
            }
            if actual_start < actual_end {
                output = output.child(
                    div()
                        .absolute()
                        .left(px(LEFT + (actual_end - start) as f32 * CELL + 6.))
                        .top(px(y + 3.))
                        .text_size(px(10.))
                        .text_color(rgb(0x95a6ba))
                        .child(if reverse { "5′" } else { "3′" }),
                );
            }
        }
        // Bond ticks occur only where both actual strands contain a residue.
        for position in start..end {
            if scene.base(position, false).is_some() && scene.base(position, true).is_some() {
                output = output.child(
                    div()
                        .absolute()
                        .left(px(LEFT + (position - start) as f32 * CELL + 4.))
                        .top(px(74.))
                        .w(px(1.))
                        .h(px(6.))
                        .bg(rgb(0x6d8297)),
                );
            }
        }
        for (span_index, span) in layout.spans.iter().enumerate() {
            let feature = &scene.axis.features[span.feature_index];
            let index = span.feature_index;
            let x = LEFT + span.start as f32 * CELL;
            let width = ((span.end - span.start) as f32 * CELL).max(3.);
            let y = 137. + span.lane as f32 * 25.;
            let color: Hsla = span
                .color_rgb
                .map(|v| rgb(v).into())
                .unwrap_or_else(|| feature_card::feature_color(feature, index));
            let filled = span.style == FeatureSegmentStyle::Filled;
            let label_color = if filled {
                annotation_label_color(color)
            } else {
                rgb(0xdce5ed)
            };
            let label = span.segment_name.as_deref().unwrap_or(&feature.label);
            let direction = feature_card::direction(feature, Molecule::Dna);
            let label = if direction == "Reverse strand" {
                format!("← {label}")
            } else if direction == "Bidirectional" {
                format!("↔ {label}")
            } else if direction == "Forward strand" {
                format!("{label} →")
            } else {
                label.to_owned()
            };
            let mut annotation = div()
                .id(("duplex-feature", row * 10000 + span_index))
                .absolute()
                .left(px(x))
                .top(px(y))
                .w(px(width))
                .h(px(19.))
                .overflow_hidden()
                .text_ellipsis()
                .text_center()
                .text_size(px(10.))
                .text_color(label_color)
                .cursor_pointer()
                .when(filled, |d| d.bg(color))
                .when(!filled, |d| d.border_b_1().border_color(color))
                .when(self.selected == Some(index), |d| {
                    d.border_1().border_color(rgb(0xffffff))
                })
                .on_mouse_down(
                    MouseButton::Left,
                    cx.listener(move |this, _, window, cx| this.select(index, window, cx)),
                );
            if width >= 30.
                && (span.style != FeatureSegmentStyle::Gap || span.segment_name.is_some())
            {
                annotation = annotation.child(label);
            }
            output = output.child(annotation);
            for boundary in &span.cleavage_sites {
                output = output.child(
                    div()
                        .absolute()
                        .left(px(LEFT + *boundary as f32 * CELL))
                        .top(px(y - 2.))
                        .w(px(2.))
                        .h(px(23.))
                        .bg(label_color),
                );
            }
        }
        if layout.omitted > 0 {
            output = output.child(
                div()
                    .absolute()
                    .left(px(LEFT))
                    .top(px(137. + layout.lanes as f32 * 25.))
                    .text_xs()
                    .text_color(rgb(0xf0bd76))
                    .child(format!(
                        "{} overlapping annotation parts omitted here; use annotation navigation",
                        layout.omitted
                    )),
            );
        }
        output.into_any_element()
    }
}
impl Render for DuplexView {
    fn render(&mut self, _: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        let mut body = div()
            .v_flex()
            .size_full()
            .min_h_0()
            .gap_2()
            .p_3()
            .bg(rgb(0x151d27))
            .text_color(rgb(0xdce5ed))
            .child(
                div()
                    .h_flex()
                    .justify_between()
                    .child(div().text_lg().child("Physical duplex"))
                    .child(
                        Button::new("close-duplex")
                            .small()
                            .label("Close")
                            .on_click(cx.listener(|_, _, _, cx| cx.emit(DuplexEvent::Close))),
                    ),
            );
        if let Some(scene) = &self.scene {
            let ends = scene.source.duplex_ends.as_ref().unwrap();
            body = body.child(div().text_sm().child(format!("{} · top {} nt · bottom {} nt · union axis {} positions", scene.source.name, scene.source.sequence.len(), ends.opposite_sequence.len(), scene.geometry.extent)))
                .child(div().text_xs().child("Top reads 5′→3′ left to right; bottom reads 5′→3′ right to left. Strand coordinates count actual residues. Gold underlines mark protruding ends; blank cells contain no strand."))
                .child(div().text_xs().child(format!("Left {:?}: {}{} · Right {:?}: {}{}", ends.left.kind, ends.left.overhang, ends.left_modification.as_ref().map_or(String::new(), |m| format!(" · {m}")), ends.right.kind, ends.right.overhang, ends.right_modification.as_ref().map_or(String::new(), |m| format!(" · {m}")))))
                .child(div().h_flex().gap_2().child(Input::new(&self.position).w(px(110.))).child(Button::new("duplex-jump").small().label("Go to axis position").on_click(cx.listener(|this, _, _, cx| this.jump(cx))))
                    .child(Button::new("duplex-copy-top").small().label("Copy top 5′→3′").on_click(cx.listener(|this, _, _, cx| { if let Some(scene) = &this.scene { cx.write_to_clipboard(ClipboardItem::new_string(scene.source.sequence.to_string())); } })))
                    .child(Button::new("duplex-copy-bottom").small().label("Copy bottom 5′→3′").on_click(cx.listener(|this, _, _, cx| { if let Some(scene) = &this.scene { cx.write_to_clipboard(ClipboardItem::new_string(scene.source.duplex_ends.as_ref().unwrap().opposite_sequence.clone())); } }))));
            let count = scene.axis.features.len();
            let can_select_top = self.selected.and_then(|i| scene.selected_top(i)).is_some();
            body = body
                .child(
                    div()
                        .relative()
                        .flex_1()
                        .min_h(px(180.))
                        .child(
                            v_virtual_list(
                                cx.entity(),
                                "physical-duplex-rows",
                                self.sizes.clone(),
                                |this, range, _, cx| {
                                    range.map(|row| this.row(row, cx)).collect::<Vec<_>>()
                                },
                            )
                            .size_full()
                            .track_scroll(&self.scroll),
                        )
                        .child(Scrollbar::vertical(&self.scroll)),
                )
                .child(
                    div()
                        .h_flex()
                        .gap_2()
                        .child(format!(
                            "Annotation {} of {count}",
                            self.selected.map_or(0, |i| i + 1)
                        ))
                        .child(
                            Button::new("duplex-prev-feature")
                                .small()
                                .disabled(count == 0)
                                .label("Previous")
                                .on_click(cx.listener(|this, _, window, cx| {
                                    this.navigate(false, window, cx)
                                })),
                        )
                        .child(
                            Button::new("duplex-next-feature")
                                .small()
                                .disabled(count == 0)
                                .label("Next")
                                .on_click(cx.listener(|this, _, window, cx| {
                                    this.navigate(true, window, cx)
                                })),
                        )
                        .child(
                            Button::new("duplex-select-top")
                                .small()
                                .disabled(!can_select_top)
                                .label("Select on top sequence")
                                .on_click(cx.listener(|this, _, _, cx| {
                                    if let Some(parts) = this
                                        .selected
                                        .and_then(|i| this.scene.as_ref()?.selected_top(i))
                                    {
                                        cx.emit(DuplexEvent::SelectTop(parts));
                                    }
                                })),
                        )
                        .child(
                            Button::new("duplex-copy-details")
                                .small()
                                .disabled(self.selected.is_none())
                                .label("Copy annotation details")
                                .on_click(cx.listener(|this, _, _, cx| {
                                    if let Some((scene, index)) =
                                        this.scene.as_ref().zip(this.selected)
                                    {
                                        cx.write_to_clipboard(ClipboardItem::new_string(
                                            scene.details(index),
                                        ));
                                    }
                                })),
                        ),
                )
                .child(
                    Textarea::new(&self.details)
                        .readonly(true)
                        .font_family("monospace")
                        .h(px(160.))
                        .aria_label("Physical strand or annotation details"),
                );
        }
        if !self.error.is_empty() {
            body = body.child(div().text_color(rgb(0xf2a4a4)).child(self.error.clone()));
        }
        body
    }
}

#[cfg(test)]
mod tests {
    use super::DuplexScene;
    use popgene_core::{
        Document, Feature, Segment,
        metadata::SequenceColor,
        restriction::{EndKind, FragmentEnd},
        workflow::DuplexEnds,
    };
    fn doc() -> Document {
        let mut doc = Document::new("physical", "ACGTAC").unwrap();
        doc.duplex_ends = Some(DuplexEnds {
            left: FragmentEnd {
                kind: EndKind::ThreePrime,
                overhang: "GG".into(),
                enzymes: vec![],
            },
            right: FragmentEnd {
                kind: EndKind::FivePrime,
                overhang: "TTT".into(),
                enzymes: vec![],
            },
            opposite_sequence: "TTTGTACGTGG".into(),
            left_modification: None,
            right_modification: None,
        });
        let feature = |id, start, end| Feature {
            id,
            label: format!("f{id}"),
            kind: "misc_feature".into(),
            parts: vec![Segment { start, end }],
            reverse: true,
            qualifiers: vec![],
            edited: false,
            location: None,
        };
        doc.sequence_qualities = Some(vec![Some(10), Some(20), None, Some(40), Some(50), Some(60)]);
        doc.opposite_sequence_qualities = Some((0..11).map(Some).collect());
        doc.features.push(feature(1, 1, 4));
        doc.duplex_annotations.features.push(feature(2, 0, 2));
        doc.duplex_annotations.features.push(feature(3, 1, 5));
        doc.sequence_colors.push(SequenceColor {
            range: Segment { start: 1, end: 3 },
            reverse: false,
            color_rgb: 0x123456,
        });
        doc.duplex_annotations.sequence_colors.push(SequenceColor {
            range: Segment { start: 0, end: 2 },
            reverse: true,
            color_rgb: 0xabcdef,
        });
        doc
    }
    #[test]
    fn actual_strands_and_colors_do_not_invent_protruding_top_bases() {
        let scene = DuplexScene::new(doc()).unwrap();
        assert_eq!(scene.geometry.extent, 11);
        assert_eq!(scene.geometry.top, Segment { start: 2, end: 8 });
        assert_eq!(scene.geometry.bottom, Segment { start: 0, end: 11 });
        assert_eq!(scene.base(0, false), None);
        assert_eq!(scene.base(8, false), None);
        assert_eq!(scene.base(2, false), Some((b'A', 0)));
        assert_eq!(scene.base(0, true), Some((b'G', 10)));
        assert_eq!(scene.base(10, true), Some((b'T', 0)));
        assert_eq!(scene.quality(0, true), Some(10));
        assert_eq!(scene.quality(10, true), Some(0));
        assert_eq!(scene.quality(2, false), Some(10));
        assert_eq!(
            crate::sequence_layout::opposite_quality(&scene.source, 0),
            Some(8)
        );
        assert_eq!(scene.colors.color(3, false), Some(0x123456));
        assert_eq!(scene.colors.color(0, true), Some(0xabcdef));
        assert_eq!(
            scene.axis.features[0].parts,
            vec![Segment { start: 3, end: 6 }]
        );
        assert_eq!(
            scene.selected_top(0),
            Some(vec![Segment { start: 1, end: 4 }])
        );
        assert_eq!(scene.selected_top(1), None);
        assert_eq!(scene.selected_top(2), None);
        assert!(scene.details(1).contains("full duplex extent"));
        assert!(scene.details(0).contains("top strand 5′→3′"));
        assert!(
            scene
                .layout
                .row(0, 11, 6, true, true, true)
                .translations
                .is_empty()
        );
    }
    #[test]
    fn every_end_polarity_maps_actual_copied_strands_after_reverse_complement() {
        let mut document = doc();
        let original = document.clone();
        for _ in 0..2 {
            let scene = DuplexScene::new(document.clone()).unwrap();
            let top = (0..scene.geometry.extent)
                .filter_map(|p| scene.base(p, false).map(|v| v.0 as char))
                .collect::<String>();
            let bottom = (0..scene.geometry.extent)
                .rev()
                .filter_map(|p| scene.base(p, true).map(|v| v.0 as char))
                .collect::<String>();
            assert_eq!(top, &*document.sequence);
            assert_eq!(
                bottom,
                document.duplex_ends.as_ref().unwrap().opposite_sequence
            );
            document.revcomp().unwrap();
        }
        assert_eq!(document.sequence, original.sequence);
        assert_eq!(document.duplex_annotations, original.duplex_annotations);
    }
    #[test]
    fn all_four_end_polarities_use_physical_coverage_in_standard_sequence_tracks() {
        for left_kind in [EndKind::FivePrime, EndKind::ThreePrime] {
            for right_kind in [EndKind::FivePrime, EndKind::ThreePrime] {
                let top_prefix = if left_kind == EndKind::FivePrime {
                    "AA"
                } else {
                    ""
                };
                let top_suffix = if right_kind == EndKind::ThreePrime {
                    "CC"
                } else {
                    ""
                };
                let bottom_prefix = if right_kind == EndKind::FivePrime {
                    "TTT"
                } else {
                    ""
                };
                let bottom_suffix = if left_kind == EndKind::ThreePrime {
                    "GG"
                } else {
                    ""
                };
                let mut source =
                    Document::new("polarities", &format!("{top_prefix}ACGT{top_suffix}")).unwrap();
                source.duplex_ends = Some(DuplexEnds {
                    left: FragmentEnd {
                        kind: left_kind,
                        overhang: if top_prefix.is_empty() {
                            bottom_suffix
                        } else {
                            top_prefix
                        }
                        .into(),
                        enzymes: vec![],
                    },
                    right: FragmentEnd {
                        kind: right_kind,
                        overhang: if top_suffix.is_empty() {
                            bottom_prefix
                        } else {
                            top_suffix
                        }
                        .into(),
                        enzymes: vec![],
                    },
                    opposite_sequence: format!("{bottom_prefix}ACGT{bottom_suffix}"),
                    left_modification: None,
                    right_modification: None,
                });
                let scene = DuplexScene::new(source.clone()).unwrap();
                for position in 0..source.sequence.len() {
                    assert_eq!(
                        crate::sequence_layout::opposite_base(&source, position),
                        scene
                            .base(position + scene.geometry.top.start, true)
                            .map(|b| b.0)
                    );
                }
                let span = crate::sequence_layout::opposite_span(&source);
                assert_eq!(
                    span,
                    Segment {
                        start: top_prefix.len(),
                        end: top_prefix.len() + 4
                    }
                );
                assert_eq!(
                    (0..source.sequence.len())
                        .filter_map(|p| crate::sequence_layout::opposite_base(&source, p))
                        .collect::<Vec<_>>(),
                    b"TGCA"
                );
            }
        }
    }
}
