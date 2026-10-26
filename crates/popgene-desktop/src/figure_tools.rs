//! Figure settings capture the current document, display options and selection.
use crate::{
    figure::FigureFormat,
    gene_map::{self, MapOptions, MapRestrictionSite},
    sequence_view::{SequenceOptions, SequenceRestrictionSite},
};
use gpui_kit::component::{
    Disableable, Selectable, Sizable, StyledExt,
    button::{Button, ButtonVariants},
    input::{Input, InputState},
};
use gpui_kit::*;
use popgene_core::{Document, Segment};

pub struct FigureContext {
    pub document: Document,
    pub map_options: MapOptions,
    pub map_sites: Vec<MapRestrictionSite>,
    pub sequence_options: SequenceOptions,
    pub sequence_sites: Vec<SequenceRestrictionSite>,
    pub selected: Vec<Segment>,
    pub feature: Option<u64>,
}
pub enum FigureEvent {
    Export {
        svg: String,
        format: FigureFormat,
        scale: f32,
        kind: &'static str,
        page_height: Option<u32>,
    },
    Cancel,
}
pub struct FigureTools {
    context: FigureContext,
    width: Entity<InputState>,
    height: Entity<InputState>,
    rows: Entity<InputState>,
    scale: Entity<InputState>,
    format: FigureFormat,
    sequence: bool,
    selection: bool,
    light: bool,
    error: String,
}
impl EventEmitter<FigureEvent> for FigureTools {}
impl FigureTools {
    pub fn new(context: FigureContext, window: &mut Window, cx: &mut Context<Self>) -> Self {
        Self {
            context,
            width: cx.new(|cx| InputState::new(window, cx).default_value("1200")),
            height: cx.new(|cx| InputState::new(window, cx).default_value("900")),
            rows: cx.new(|cx| InputState::new(window, cx).default_value("100")),
            scale: cx.new(|cx| InputState::new(window, cx).default_value("3.125")),
            format: FigureFormat::Pdf,
            sequence: false,
            selection: false,
            light: true,
            error: String::new(),
        }
    }
    fn export(&mut self, window: &mut Window, cx: &mut Context<Self>) {
        let result: Result<_, String> = (|| {
            let width = self
                .width
                .read(cx)
                .value()
                .parse::<u32>()
                .map_err(|_| "Enter a whole-number width")?;
            let scale = if self.format == FigureFormat::Png {
                let scale = self
                    .scale
                    .read(cx)
                    .value()
                    .parse::<f32>()
                    .map_err(|_| "Enter a raster scale")?;
                if !scale.is_finite() || !(0.25..=8.).contains(&scale) {
                    return Err("Raster scale must be between 0.25 and 8".into());
                }
                scale
            } else {
                1.
            };
            let svg = if self.sequence {
                let rows = self
                    .rows
                    .read(cx)
                    .value()
                    .parse::<usize>()
                    .map_err(|_| "Enter a whole-number row limit")?;
                let regions = if self.selection {
                    self.context.selected.as_slice()
                } else {
                    &[]
                };
                crate::sequence_figure::svg(
                    &self.context.document,
                    &self.context.sequence_options,
                    &self.context.sequence_sites,
                    regions,
                    width as f32,
                    rows,
                    !self.light,
                    self.context.feature,
                )?
            } else {
                let height = self
                    .height
                    .read(cx)
                    .value()
                    .parse::<u32>()
                    .map_err(|_| "Enter a whole-number height")?;
                let font = font("Arial");
                gene_map::map_svg(
                    &self.context.document,
                    self.context.map_options,
                    &self.context.map_sites,
                    width,
                    height,
                    self.light,
                    |text, size| {
                        window
                            .text_system()
                            .shape_line(
                                text.to_owned().into(),
                                px(size),
                                &[TextRun {
                                    len: text.len(),
                                    font: font.clone(),
                                    color: rgb(0x172b3b).into(),
                                    background_color: None,
                                    underline: None,
                                    strikethrough: None,
                                }],
                                None,
                            )
                            .width
                            .as_f32()
                    },
                )?
            };
            let page_height = if self.sequence && self.format == FigureFormat::Pdf {
                let height = self
                    .height
                    .read(cx)
                    .value()
                    .parse::<u32>()
                    .map_err(|_| "Enter a whole-number PDF page height")?;
                if !(400..=8000).contains(&height) {
                    return Err("PDF page height must be 400–8,000 layout pixels".into());
                }
                Some(height)
            } else {
                None
            };
            Ok((svg, scale, page_height))
        })();
        match result {
            Ok((svg, scale, page_height)) => cx.emit(FigureEvent::Export {
                svg,
                format: self.format,
                scale,
                kind: if self.sequence { "sequence" } else { "map" },
                page_height,
            }),
            Err(error) => {
                self.error = error;
                cx.notify();
            }
        }
    }
}
impl Render for FigureTools {
    fn render(&mut self, _: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        let mut body=div().v_flex().gap_3().p_3().child(div().text_xl().child("Export figure"))
            .child(div().text_sm().child(self.context.document.name.clone()))
            .child(div().h_flex().gap_2().children([(false,"Map"),(true,"Sequence")].into_iter().enumerate().map(|(i,(sequence,label))|Button::new(("figure-content",i)).small().selected(self.sequence==sequence).label(label).on_click(cx.listener(move |this,_,_,cx|{this.sequence=sequence;this.error.clear();cx.notify();})))))
            .child(div().text_xs().child(if self.sequence {"Uses current sequence layers, genetic code and selected-CDS flank settings. Original coordinates are retained. Any omitted rows or crowded labels are explicitly marked."} else {"Uses current map layers and enzyme filter, fitted to figure dimensions. Dense-label omissions remain labelled."}))
            .child(div().h_flex().gap_2().children([FigureFormat::Svg,FigureFormat::Pdf,FigureFormat::Png].into_iter().enumerate().map(|(i,format)|Button::new(("figure-format",i)).small().selected(self.format==format).label(format.label()).on_click(cx.listener(move |this,_,_,cx|{this.format=format;cx.notify();})))))
            .child(div().h_flex().gap_3().child(div().v_flex().gap_1().child("Width (px)").child(Input::new(&self.width).w(px(160.)))).child(div().v_flex().gap_1().child(if self.sequence {"Maximum rows"} else {"Height (px)"}).child(Input::new(if self.sequence {&self.rows} else {&self.height}).w(px(160.)))));
        if self.sequence {
            if self.format == FigureFormat::Pdf {
                body=body.child(div().h_flex().gap_3().child("PDF page height (px)").child(Input::new(&self.height).w(px(160.))))
                    .child(div().text_xs().child("PDF pages break between complete sequence rows and retain all exported rows."));
            }
            body = body.child(
                div().h_flex().gap_2().children(
                    [(false, "Whole sequence"), (true, "Current selection")]
                        .into_iter()
                        .enumerate()
                        .map(|(i, (selection, label))| {
                            Button::new(("figure-scope", i))
                                .small()
                                .selected(self.selection == selection)
                                .disabled(selection && self.context.selected.is_empty())
                                .label(label)
                                .on_click(cx.listener(move |this, _, _, cx| {
                                    this.selection = selection;
                                    cx.notify();
                                }))
                        }),
                ),
            );
        }
        body=body.child(div().h_flex().gap_2().children([(true,"White background"),(false,"Dark background")].into_iter().enumerate().map(|(i,(light,label))|Button::new(("figure-background",i)).small().selected(self.light==light).label(label).on_click(cx.listener(move |this,_,_,cx|{this.light=light;cx.notify();})))))
            .child(div().h_flex().gap_3().child("PNG raster scale").child(Input::new(&self.scale).w(px(160.)).disabled(self.format!=FigureFormat::Png)))
            .child(div().text_xs().child("Scale 3.125 produces 300 pixels per nominal inch. SVG and PDF use 96 layout pixels per inch. SVG text is outlined; PDF remains vector."));
        if !self.error.is_empty() {
            body = body.child(div().text_color(rgb(0xef938a)).child(self.error.clone()));
        }
        body.child(
            div()
                .h_flex()
                .justify_end()
                .gap_2()
                .child(
                    Button::new("cancel-figure")
                        .label("Cancel")
                        .on_click(cx.listener(|_, _, _, cx| cx.emit(FigureEvent::Cancel))),
                )
                .child(
                    Button::new("export-figure")
                        .primary()
                        .label("Choose destination…")
                        .on_click(cx.listener(|this, _, window, cx| this.export(window, cx))),
                ),
        )
    }
}
