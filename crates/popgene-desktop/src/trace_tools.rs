//! Non-destructive trace editing preview; Apply is one parent-model undo step.
use crate::sequence_view::{TraceEvent, TraceView};
use gpui_kit::component::{
    Disableable, Sizable, StyledExt,
    button::{Button, ButtonVariants},
    input::{Input, InputState, Textarea, TextareaState},
};
use gpui_kit::*;
use popgene_core::{Document, Segment, trace};

pub enum TraceToolsEvent {
    Apply(Box<Document>),
    Cancel,
    Align,
}
pub struct TraceTools {
    original: Document,
    document: Document,
    undo: Vec<Document>,
    view: Entity<TraceView>,
    fields: Vec<Entity<InputState>>,
    table: Entity<TextareaState>,
    error: String,
    _subscription: Subscription,
}
impl EventEmitter<TraceToolsEvent> for TraceTools {}
impl TraceTools {
    pub fn new(document: Document, window: &mut Window, cx: &mut Context<Self>) -> Self {
        let view = cx.new(|cx| TraceView::new(document.clone(), 1., cx));
        let fields = ["1", "1", "N", "0", "20", "10"]
            .into_iter()
            .map(|value| cx.new(|cx| InputState::new(window, cx).default_value(value)))
            .collect();
        let table = cx.new(|cx| {
            TextareaState::new(window, cx)
                .rows(8)
                .default_value(call_table(&document))
        });
        let subscription = cx.subscribe_in(&view, window, |this, _, event, window, cx| {
            let TraceEvent::SelectCall(index) = event;
            for field in &this.fields[..2] {
                field.update(cx, |s, cx| s.set_value((index + 1).to_string(), window, cx));
            }
            if let Some(trace) = &this.document.trace {
                this.fields[3].update(cx, |s, cx| {
                    s.set_value(trace.peaks[*index].to_string(), window, cx)
                });
            }
            cx.notify();
        });
        Self {
            original: document.clone(),
            document,
            undo: vec![],
            view,
            fields,
            table,
            error: String::new(),
            _subscription: subscription,
        }
    }
    pub fn error(&mut self, message: impl Into<String>, cx: &mut Context<Self>) {
        self.error = message.into();
        cx.notify();
    }
    fn number(&self, index: usize, cx: &Context<Self>) -> Result<usize, String> {
        self.fields[index]
            .read(cx)
            .value()
            .trim()
            .parse()
            .map_err(|_| "Enter whole-number call and sample positions".into())
    }
    fn range(&self, cx: &Context<Self>) -> Result<Segment, String> {
        let start = self
            .number(0, cx)?
            .checked_sub(1)
            .ok_or("Call coordinates start at 1")?;
        let end = self.number(1, cx)?;
        if start >= end || end > self.document.sequence.len() {
            return Err("Call interval is outside this read".into());
        }
        Ok(Segment { start, end })
    }
    fn refresh(&mut self, window: &mut Window, cx: &mut Context<Self>) {
        let length = self.document.sequence.len();
        for field in &self.fields[..2] {
            let current = field
                .read(cx)
                .value()
                .parse::<usize>()
                .unwrap_or(1)
                .clamp(1, length.max(1));
            field.update(cx, |state, cx| {
                state.set_value(current.to_string(), window, cx)
            });
        }
        self.view.update(cx, |view, cx| {
            view.set_document(self.document.clone(), 1., cx)
        });
        self.table.update(cx, |state, cx| {
            state.set_value(call_table(&self.document), window, cx)
        });
        cx.notify();
    }
    fn edit(
        &mut self,
        operation: impl FnOnce(&mut Document) -> Result<(), popgene_core::Error>,
        window: &mut Window,
        cx: &mut Context<Self>,
    ) {
        let previous = self.document.clone();
        match operation(&mut self.document) {
            Ok(()) => {
                self.undo.push(previous);
                while self.undo.len() > 1
                    && (self.undo.len() > 100
                        || self
                            .undo
                            .iter()
                            .map(Document::snapshot_bytes)
                            .sum::<usize>()
                            > 64 * 1024 * 1024)
                {
                    self.undo.remove(0);
                }
                self.error.clear();
                self.refresh(window, cx);
            }
            Err(error) => self.error(error.to_string(), cx),
        }
    }
    fn action(&mut self, kind: usize, window: &mut Window, cx: &mut Context<Self>) {
        let result = (|| -> Result<(), String> {
            if kind == 4 {
                let minimum =
                    u8::try_from(self.number(4, cx)?).map_err(|_| "Quality must be 0–255")?;
                let window_size = self.number(5, cx)?;
                let keep = trace::quality_trim_range(&self.document, minimum, window_size)
                    .map_err(|e| e.to_string())?;
                self.fields[0].update(cx, |s, cx| {
                    s.set_value((keep.start + 1).to_string(), window, cx)
                });
                self.fields[1].update(cx, |s, cx| s.set_value(keep.end.to_string(), window, cx));
                self.error.clear();
                cx.notify();
                return Ok(());
            }
            let range = self.range(cx)?;
            let base = self.fields[2].read(cx).value().to_string();
            match kind {
                0 => self.edit(
                    |doc| trace::replace_call(doc, range.start, &base),
                    window,
                    cx,
                ),
                1 => {
                    let sample = self.number(3, cx)?;
                    self.edit(
                        |doc| trace::insert_call(doc, range.start, &base, sample),
                        window,
                        cx,
                    );
                }
                2 => self.edit(|doc| trace::delete_calls(doc, range), window, cx),
                3 => self.edit(|doc| trace::trim_calls(doc, range), window, cx),
                5 => {
                    let sample = self.number(3, cx)?;
                    self.edit(
                        |doc| trace::insert_call(doc, range.end, &base, sample),
                        window,
                        cx,
                    );
                }
                _ => unreachable!(),
            }
            Ok(())
        })();
        if let Err(error) = result {
            self.error(error, cx);
        }
    }
}
fn call_table(document: &Document) -> String {
    let Some(trace) = &document.trace else {
        return "No chromatogram".into();
    };
    let mut text="Call\tBase\tSignal sample (0-based)\tQuality\tManual\tA signal\tC signal\tG signal\tT signal\n".to_string();
    for (i, base) in document.sequence.bytes().take(5000).enumerate() {
        let peak = trace.peaks[i];
        let signals = b"ACGT"
            .iter()
            .map(|base| {
                trace
                    .channel_order
                    .as_bytes()
                    .iter()
                    .position(|b| b == base)
                    .map_or("—".into(), |channel| {
                        trace.channels[channel][peak].to_string()
                    })
            })
            .collect::<Vec<_>>()
            .join("\t");
        text.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            i + 1,
            base as char,
            peak,
            if trace.qualities_known {
                trace.qualities[i].to_string()
            } else {
                "unknown".into()
            },
            trace.manual_calls.get(i) == Some(&true),
            signals
        ));
    }
    if document.sequence.len() > 5000 {
        text.push_str("Table limited to the first 5,000 calls; the plot and coordinate editing support the full read.\n");
    }
    text
}
impl Render for TraceTools {
    fn render(&mut self, _: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        let mut body=div().id("trace-tools-scroll").v_flex().gap_3().p_4().size_full().min_h_0().overflow_y_scroll().child(div().text_lg().child(format!("Trace · {}",self.document.name)))
            .child(div().text_xs().child("Click a called base in the plot to select its position. Manual calls use zero quality and keep a real sample position. Edits preserve the imported call families and imported signals; Apply records one undo step."));
        let mut families = div().h_flex().gap_2();
        if let Some(trace) = &self.document.trace {
            for (index, calls) in trace.original_calls.iter().enumerate() {
                families = families.child(
                    Button::new(("trace-call-family", index))
                        .small()
                        .label(format!("{} · {} calls", calls.label, calls.sequence.len()))
                        .on_click(cx.listener(move |this, _, window, cx| {
                            this.edit(|doc| trace::choose_call_set(doc, index), window, cx)
                        })),
                );
            }
        }
        body = body
            .child(families)
            .child(div().h(px(500.)).flex_shrink_0().child(self.view.clone()));
        let mut fields = div().h_flex().gap_3();
        for (index, label) in [
            "First call (1-based)",
            "Last call (inclusive)",
            "New base",
            "Sample for insertion (0-based)",
        ]
        .into_iter()
        .enumerate()
        {
            fields = fields.child(
                div()
                    .v_flex()
                    .gap_1()
                    .child(div().text_xs().child(label))
                    .child(Input::new(&self.fields[index]).w(px(if index == 3 {
                        210.
                    } else {
                        130.
                    }))),
            );
        }
        body = body.child(fields);
        let mut operations = div().h_flex().gap_2();
        for (index, label) in [
            (0, "Replace first call"),
            (1, "Insert before first call"),
            (5, "Insert after last call"),
            (2, "Delete interval"),
            (3, "Keep interval / trim"),
        ] {
            operations = operations.child(
                Button::new(("trace-edit-action", index))
                    .small()
                    .label(label)
                    .on_click(
                        cx.listener(move |this, _, window, cx| this.action(index, window, cx)),
                    ),
            );
        }
        body = body
            .child(operations)
            .child(
                div()
                    .h_flex()
                    .gap_2()
                    .child(div().text_sm().child("Quality window mean ≥"))
                    .child(Input::new(&self.fields[4]).w(px(70.)))
                    .child(div().text_sm().child("over calls"))
                    .child(Input::new(&self.fields[5]).w(px(70.)))
                    .child(
                        Button::new("trace-propose-quality-trim")
                            .small()
                            .label("Select longest passing interval")
                            .on_click(
                                cx.listener(|this, _, window, cx| this.action(4, window, cx)),
                            ),
                    ),
            )
            .child(
                Textarea::new(&self.table)
                    .readonly(true)
                    .font_family("monospace")
                    .h(px(180.))
                    .aria_label("Base calls, quality and channel signal values"),
            );
        if !self.error.is_empty() {
            body = body.child(div().text_color(rgb(0xf2a4a4)).child(self.error.clone()));
        }
        body.child(
            div()
                .h_flex()
                .gap_2()
                .child(
                    Button::new("trace-apply")
                        .primary()
                        .label("Apply edits")
                        .disabled(self.document == self.original)
                        .on_click(cx.listener(|this, _, _, cx| {
                            cx.emit(TraceToolsEvent::Apply(Box::new(this.document.clone())))
                        })),
                )
                .child(
                    Button::new("trace-preview-undo")
                        .label("Undo preview edit")
                        .disabled(self.undo.is_empty())
                        .on_click(cx.listener(|this, _, window, cx| {
                            if let Some(previous) = this.undo.pop() {
                                this.document = previous;
                                this.error.clear();
                                this.refresh(window, cx);
                            }
                        })),
                )
                .child(
                    Button::new("trace-align")
                        .label("Align to open reference")
                        .disabled(self.document != self.original)
                        .on_click(cx.listener(|_, _, _, cx| cx.emit(TraceToolsEvent::Align))),
                )
                .child(
                    Button::new("trace-cancel")
                        .ghost()
                        .label("Cancel")
                        .on_click(cx.listener(|_, _, _, cx| cx.emit(TraceToolsEvent::Cancel))),
                ),
        )
    }
}
