//! Editable citations and explicit extraction of original embedded attachments.
use gpui_kit::component::{
    Disableable, Selectable, Sizable, StyledExt,
    button::{Button, ButtonVariants},
    input::{Input, InputState, Textarea, TextareaState},
};
use gpui_kit::*;
use popgene_core::{Document, Source, metadata::Reference};
use popgene_formats::{NativeAttachment, NativeEmbeddedTrace, NativeHistoryEntry};

fn edited_optional(original: &Option<String>, value: &str) -> Option<String> {
    if value == original.as_deref().unwrap_or("") {
        original.clone()
    } else {
        (!value.trim().is_empty()).then(|| value.to_owned())
    }
}

pub enum MetadataEvent {
    Save(Vec<Reference>),
    CreateDocument {
        document: Box<Document>,
        message: &'static str,
    },
    ExportAttachment {
        id: u32,
        name: String,
    },
    Cancel,
}
#[derive(Clone, Copy)]
enum StoredDocument {
    History(u32),
    Trace(usize),
}
pub struct MetadataTools {
    references: Vec<Reference>,
    original_references: Vec<Reference>,
    selected: Option<usize>,
    fields: Vec<(&'static str, Entity<InputState>)>,
    remark: Entity<TextareaState>,
    extras: Vec<(Entity<InputState>, Entity<InputState>)>,
    attachments: Vec<NativeAttachment>,
    embedded_traces: Vec<NativeEmbeddedTrace>,
    trace_page: usize,
    selected_trace: Option<usize>,
    trace_details: Entity<TextareaState>,
    history: Vec<NativeHistoryEntry>,
    history_page: usize,
    selected_history: Option<u32>,
    history_details: Entity<TextareaState>,
    source: Option<Source>,
    restoring: bool,
    loading: bool,
    error: String,
}
impl EventEmitter<MetadataEvent> for MetadataTools {}
impl MetadataTools {
    pub fn new(
        references: Vec<Reference>,
        source: Option<Source>,
        window: &mut Window,
        cx: &mut Context<Self>,
    ) -> Self {
        let fields = [
            "Title",
            "Authors",
            "Consortium",
            "Journal",
            "PubMed ID",
            "Original reference scope",
        ]
        .into_iter()
        .map(|label| (label, cx.new(|cx| InputState::new(window, cx))))
        .collect();
        let mut result = Self {
            original_references: references.clone(),
            references,
            selected: None,
            fields,
            remark: cx.new(|cx| TextareaState::new(window, cx).rows(3)),
            extras: vec![],
            attachments: vec![],
            embedded_traces: vec![],
            trace_page: 0,
            selected_trace: None,
            trace_details: cx.new(|cx| TextareaState::new(window, cx).rows(6)),
            history: vec![],
            history_page: 0,
            selected_history: None,
            history_details: cx.new(|cx| TextareaState::new(window, cx).rows(8)),
            source: source.clone(),
            restoring: false,
            loading: false,
            error: String::new(),
        };
        if !result.references.is_empty() {
            result.load(0, window, cx);
        }
        if let Some(source) = source.filter(|s| s.format == "SnapGene") {
            result.loading = true;
            cx.spawn_in(window, async move |this, cx| {
                let loaded = cx
                    .background_executor()
                    .spawn(async move {
                        (
                            popgene_formats::native_attachments(&source),
                            popgene_formats::native_history(&source),
                            popgene_formats::native_embedded_traces(&source),
                        )
                    })
                    .await;
                let _ = this.update_in(cx, |this, _, cx| {
                    this.loading = false;
                    match loaded.0 {
                        Ok(a) => this.attachments = a,
                        Err(e) => this.error = format!("Attachments could not be read: {e}"),
                    };
                    match loaded.1 {
                        Ok(history) => this.history = history,
                        Err(e) => this
                            .error
                            .push_str(&format!(" History could not be read: {e}")),
                    }
                    match loaded.2 {
                        Ok(traces) => this.embedded_traces = traces,
                        Err(e) => this
                            .error
                            .push_str(&format!(" Embedded traces could not be read: {e}")),
                    }
                    cx.notify();
                });
            })
            .detach();
        }
        result
    }
    pub fn error(&mut self, message: impl Into<String>, cx: &mut Context<Self>) {
        self.error = message.into();
        cx.notify();
    }
    fn inspect_history(&mut self, id: u32, window: &mut Window, cx: &mut Context<Self>) {
        let Some(entry) = self.history.iter().find(|entry| entry.id == id) else {
            return;
        };
        let mut text = format!("{} — {}\n", entry.name, entry.operation);
        for (key, value) in &entry.attributes {
            text.push_str(&format!("Node/@{key}: {value}\n"));
        }
        for (key, value) in &entry.details {
            text.push_str(&format!("{key}: {value}\n"));
        }
        if entry.details_truncated {
            text.push_str("\nDisplay limit reached. Additional metadata remains in the retained original source.\n");
        }
        self.history_details
            .update(cx, |field, cx| field.set_value(text, window, cx));
        self.selected_history = Some(id);
        cx.notify();
    }
    fn inspect_trace(&mut self, index: usize, window: &mut Window, cx: &mut Context<Self>) {
        let Some(trace) = self
            .embedded_traces
            .iter()
            .find(|trace| trace.index == index)
        else {
            return;
        };
        let mut text = format!(
            "{}\nStored container orientation: {}\n{} embedded bytes\n",
            trace.name,
            if trace.reverse { "reverse" } else { "forward" },
            trace.byte_len
        );
        for (key, value) in &trace.properties {
            text.push_str(&format!("{key}: {value}\n"));
        }
        self.trace_details
            .update(cx, |field, cx| field.set_value(text, window, cx));
        self.selected_trace = Some(index);
        cx.notify();
    }
    fn store(&mut self, cx: &App) {
        let Some(index) = self.selected else {
            return;
        };
        let values = self
            .fields
            .iter()
            .map(|(_, field)| field.read(cx).value().to_string())
            .collect::<Vec<_>>();
        let original = self.references[index].clone();
        self.references[index] = Reference {
            title: values[0].clone(),
            authors: edited_optional(&original.authors, &values[1]),
            consortium: edited_optional(&original.consortium, &values[2]),
            journal: edited_optional(&original.journal, &values[3]),
            pubmed: edited_optional(&original.pubmed, &values[4]),
            scope: edited_optional(&original.scope, &values[5]),
            remark: edited_optional(&original.remark, &self.remark.read(cx).value()),
            extra: self
                .extras
                .iter()
                .enumerate()
                .filter_map(|(i, (key, value))| {
                    let input = key.read(cx).value();
                    let key = if original
                        .extra
                        .get(i)
                        .is_some_and(|(key, _)| key == input.as_ref())
                    {
                        input.to_string()
                    } else {
                        input.trim().to_string()
                    };
                    (!key.is_empty()).then(|| (key, value.read(cx).value().to_string()))
                })
                .collect(),
        };
    }
    fn load(&mut self, index: usize, window: &mut Window, cx: &mut Context<Self>) {
        self.selected = Some(index);
        let r = &self.references[index];
        let values = [
            r.title.as_str(),
            r.authors.as_deref().unwrap_or(""),
            r.consortium.as_deref().unwrap_or(""),
            r.journal.as_deref().unwrap_or(""),
            r.pubmed.as_deref().unwrap_or(""),
            r.scope.as_deref().unwrap_or(""),
        ];
        for ((_, field), value) in self.fields.iter().zip(values) {
            field.update(cx, |field, cx| {
                field.set_value(value.to_string(), window, cx)
            });
        }
        self.remark.update(cx, |field, cx| {
            field.set_value(r.remark.clone().unwrap_or_default(), window, cx)
        });
        self.extras = r
            .extra
            .iter()
            .map(|(key, value)| {
                (
                    cx.new(|cx| InputState::new(window, cx).default_value(key.clone())),
                    cx.new(|cx| InputState::new(window, cx).default_value(value.clone())),
                )
            })
            .collect();
    }
    fn open_stored(&mut self, entry: StoredDocument, window: &mut Window, cx: &mut Context<Self>) {
        if self.restoring {
            return;
        }
        self.store(cx);
        if self.references != self.original_references {
            self.error(
                "Apply or cancel the pending citation edits before opening a stored document.",
                cx,
            );
            return;
        }
        let Some(source) = self.source.clone() else {
            return;
        };
        self.restoring = true;
        cx.notify();
        cx.spawn_in(window, async move |this, cx| {
            let restored = cx
                .background_executor()
                .spawn(async move { match entry {
                    StoredDocument::History(id) => popgene_formats::native_history_document(&source, id),
                    StoredDocument::Trace(index) => popgene_formats::native_embedded_trace(&source, index),
                } })
                .await;
            let _ = this.update_in(cx, |this, _, cx| {
                this.restoring = false;
                this.store(cx);
                if this.references != this.original_references {
                    this.error("Citation edits are pending. Apply or cancel them, then reopen the stored document.",cx);
                    return;
                }
                match restored {
                    Ok(mut doc) => {
                        let message = match entry {
                            StoredDocument::History(id) => {
                                doc.name = format!("{} history {}", doc.name, id);
                                "History snapshot restored · review and save the new document"
                            },
                            StoredDocument::Trace(_) => "Embedded trace opened from the original source · review and save the new document",
                        };
                        cx.emit(MetadataEvent::CreateDocument { document: Box::new(doc), message });
                    }
                    Err(e) => this.error = e.to_string(),
                };
                cx.notify();
            });
        })
        .detach();
    }
    fn save(&mut self, cx: &mut Context<Self>) {
        self.store(cx);
        if self.references.len() > 1024 || self.references.iter().any(|r| r.extra.len() > 256) {
            self.error(
                "Use at most 1,024 references and 256 additional fields per reference",
                cx,
            );
            return;
        }
        cx.emit(MetadataEvent::Save(self.references.clone()));
    }
}
impl Render for MetadataTools {
    fn render(&mut self, _: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        let mut body=div().id("metadata-tools-scroll").v_flex().gap_3().p_3().max_h(px(760.)).overflow_y_scroll()
            .child(div().text_xl().child("References, attachments & history"))
            .child(div().text_xs().child("Citations are editable document data. Reference scope retains the source's numbering; review it after sequence edits."));
        let mut list = div()
            .id("reference-list")
            .v_flex()
            .gap_1()
            .h(px((self.references.len().min(5) * 32) as f32))
            .flex_shrink_0()
            .overflow_y_scroll();
        for (i, reference) in self.references.iter().enumerate() {
            let title = if reference.title.is_empty() {
                format!("Reference {}", i + 1)
            } else {
                format!("{}. {}", i + 1, reference.title)
            };
            list = list.child(
                Button::new(("reference", i))
                    .small()
                    .flex_shrink_0()
                    .selected(self.selected == Some(i))
                    .label(title)
                    .on_click(cx.listener(move |this, _, window, cx| {
                        this.store(cx);
                        this.load(i, window, cx);
                        cx.notify();
                    })),
            );
        }
        body = body.child(list).child(
            div()
                .h_flex()
                .gap_2()
                .child(
                    Button::new("new-reference")
                        .small()
                        .label("Add reference")
                        .disabled(self.references.len() >= 1024)
                        .on_click(cx.listener(|this, _, window, cx| {
                            this.store(cx);
                            this.references.push(Reference::default());
                            this.load(this.references.len() - 1, window, cx);
                            cx.notify();
                        })),
                )
                .child(
                    Button::new("remove-reference")
                        .small()
                        .label("Remove selected")
                        .disabled(self.selected.is_none())
                        .on_click(cx.listener(|this, _, window, cx| {
                            if let Some(i) = this.selected.take() {
                                this.references.remove(i);
                                if !this.references.is_empty() {
                                    this.load(i.min(this.references.len() - 1), window, cx);
                                }
                            }
                            cx.notify();
                        })),
                ),
        );
        if self.selected.is_some() {
            body = body
                .children(self.fields.iter().map(|(label, input)| {
                    div()
                        .v_flex()
                        .gap_1()
                        .child(div().text_xs().child(*label))
                        .child(Input::new(input))
                }))
                .child(div().text_xs().child("Remark"))
                .child(Textarea::new(&self.remark).h(px(75.)))
                .child(div().text_xs().child("Additional citation fields"));
            for (i, (key, value)) in self.extras.iter().enumerate() {
                body = body.child(
                    div()
                        .h_flex()
                        .gap_2()
                        .child(Input::new(key).w(px(180.)))
                        .child(Input::new(value).flex_1())
                        .child(
                            Button::new(("remove-reference-field", i))
                                .small()
                                .label("Remove")
                                .on_click(cx.listener(move |this, _, _, cx| {
                                    this.extras.remove(i);
                                    cx.notify();
                                })),
                        ),
                );
            }
            body = body.child(
                Button::new("add-reference-field")
                    .small()
                    .label("Add citation field")
                    .disabled(self.extras.len() >= 256)
                    .on_click(cx.listener(|this, _, window, cx| {
                        let key =
                            cx.new(|cx| InputState::new(window, cx).placeholder("Field name"));
                        let value = cx.new(|cx| InputState::new(window, cx).placeholder("Value"));
                        this.extras.push((key, value));
                        cx.notify();
                    })),
            );
        }
        body = body.child(
            div()
                .mt_3()
                .text_lg()
                .child("Original embedded attachments"),
        );
        if self.loading {
            body = body.child("Reading attachment records…");
        } else if self.attachments.is_empty() {
            body = body.child(
                div()
                    .text_xs()
                    .child("No embedded attachments found in the original source."),
            );
        }
        for attachment in &self.attachments {
            let id = attachment.id;
            let name = attachment.name.clone();
            body = body.child(
                div()
                    .h_flex()
                    .gap_2()
                    .child(div().flex_1().text_sm().child(format!(
                            "{} · {} bytes{}",
                            attachment.name,
                            attachment.size,
                            attachment
                                .modified
                                .as_ref()
                                .map(|d| format!(" · {d}"))
                                .unwrap_or_default()
                        )))
                    .child(
                        Button::new(("export-attachment", id as usize))
                            .small()
                            .label("Save a copy…")
                            .on_click(cx.listener(move |_, _, _, cx| {
                                cx.emit(MetadataEvent::ExportAttachment {
                                    id,
                                    name: name.clone(),
                                })
                            })),
                    ),
            );
        }
        if !self.embedded_traces.is_empty() {
            body = body.child(div().mt_3().text_lg().child("Original embedded traces"))
                .child(div().text_xs().child("Open creates a separate unsaved trace document with its original calls and signal. Container orientation is recorded separately; it does not change the stored calls."));
            for trace in self
                .embedded_traces
                .iter()
                .skip(self.trace_page * 25)
                .take(25)
            {
                let index = trace.index;
                body = body.child(
                    div()
                        .h_flex()
                        .gap_2()
                        .child(div().flex_1().text_sm().child(format!(
                            "{} · {} bytes · {} orientation",
                            trace.name,
                            trace.byte_len,
                            if trace.reverse { "reverse" } else { "forward" }
                        )))
                        .child(
                            Button::new(("embedded-trace-details", index))
                                .small()
                                .label("Details")
                                .selected(self.selected_trace == Some(index))
                                .on_click(cx.listener(move |this, _, window, cx| {
                                    this.inspect_trace(index, window, cx)
                                })),
                        )
                        .child(
                            Button::new(("open-embedded-trace", index))
                                .small()
                                .label("Open trace")
                                .disabled(self.restoring)
                                .on_click(cx.listener(move |this, _, window, cx| {
                                    this.open_stored(StoredDocument::Trace(index), window, cx)
                                })),
                        ),
                );
            }
            if self.embedded_traces.len() > 25 {
                body = body.child(
                    div()
                        .h_flex()
                        .gap_2()
                        .child(
                            Button::new("previous-trace-page")
                                .small()
                                .label("Previous")
                                .disabled(self.trace_page == 0)
                                .on_click(cx.listener(|this, _, _, cx| {
                                    this.trace_page = this.trace_page.saturating_sub(1);
                                    cx.notify();
                                })),
                        )
                        .child(format!(
                            "Trace page {} of {}",
                            self.trace_page + 1,
                            self.embedded_traces.len().div_ceil(25)
                        ))
                        .child(
                            Button::new("next-trace-page")
                                .small()
                                .label("Next")
                                .disabled((self.trace_page + 1) * 25 >= self.embedded_traces.len())
                                .on_click(cx.listener(|this, _, _, cx| {
                                    this.trace_page += 1;
                                    cx.notify();
                                })),
                        ),
                );
            }
            if self.selected_trace.is_some() {
                body = body
                    .child(
                        Button::new("copy-embedded-trace-details")
                            .small()
                            .label("Copy trace details")
                            .on_click(cx.listener(|this, _, _, cx| {
                                cx.write_to_clipboard(ClipboardItem::new_string(
                                    this.trace_details.read(cx).value().to_string(),
                                ))
                            })),
                    )
                    .child(
                        Textarea::new(&self.trace_details)
                            .readonly(true)
                            .h(px(180.))
                            .aria_label("Original embedded trace details"),
                    );
            }
        }
        body=body.child(div().mt_3().text_lg().child("Original procedure history"))
            .child(div().text_xs().child("Restore opens a separate unsaved document. Unsupported operations stay in the original source and cannot be replayed."));
        if self.history.is_empty() && !self.loading {
            body = body.child(
                div()
                    .text_xs()
                    .child("No procedure tree found in the original source."),
            );
        }
        let history_rows = self
            .history
            .len()
            .saturating_sub(self.history_page * 100)
            .min(100);
        let mut history = div()
            .id("procedure-history-list")
            .v_flex()
            .gap_2()
            .h(px((history_rows.min(4) * 76) as f32))
            .max_h(px(300.))
            .flex_shrink_0()
            .overflow_y_scroll();
        for entry in self.history.iter().skip(self.history_page * 100).take(100) {
            let id = entry.id;
            history = history.child(
                div()
                    .v_flex()
                    .flex_shrink_0()
                    .gap_1()
                    .p_2()
                    .border_1()
                    .border_color(rgb(0x344356))
                    .child(
                        div()
                            .h_flex()
                            .gap_2()
                            .child(div().flex_1().text_sm().child(format!(
                                    "#{} {} · {} · recorded extent {} {}{}",
                                    entry.id,
                                    entry.name,
                                    entry.operation,
                                    entry.sequence_length,
                                    entry.molecule.unit(),
                                    entry
                                        .parent
                                        .map(|parent| format!(" · parent #{parent}"))
                                        .unwrap_or_default()
                                )))
                            .child(
                                Button::new(("inspect-history", id as usize))
                                    .small()
                                    .label("Details")
                                    .selected(self.selected_history == Some(id))
                                    .on_click(cx.listener(move |this, _, window, cx| {
                                        this.inspect_history(id, window, cx)
                                    })),
                            )
                            .child(
                                Button::new(("restore-history", id as usize))
                                    .small()
                                    .label("Restore as new")
                                    .disabled(!entry.restorable || self.restoring)
                                    .on_click(cx.listener(move |this, _, window, cx| {
                                        this.open_stored(StoredDocument::History(id), window, cx)
                                    })),
                            ),
                    ),
            );
        }
        body = body.child(history);
        if self.history.len() > 100 {
            body = body.child(
                div()
                    .h_flex()
                    .gap_2()
                    .child(
                        Button::new("previous-history-page")
                            .small()
                            .label("Previous")
                            .disabled(self.history_page == 0)
                            .on_click(cx.listener(|this, _, _, cx| {
                                this.history_page = this.history_page.saturating_sub(1);
                                cx.notify();
                            })),
                    )
                    .child(format!(
                        "Nodes {}–{} of {}",
                        self.history_page * 100 + 1,
                        ((self.history_page + 1) * 100).min(self.history.len()),
                        self.history.len()
                    ))
                    .child(
                        Button::new("next-history-page")
                            .small()
                            .label("Next")
                            .disabled((self.history_page + 1) * 100 >= self.history.len())
                            .on_click(cx.listener(|this, _, _, cx| {
                                this.history_page += 1;
                                cx.notify();
                            })),
                    ),
            );
        }
        if let Some(id) = self.selected_history {
            body = body
                .child(
                    div()
                        .h_flex()
                        .gap_2()
                        .child(format!("History #{id} details"))
                        .child(
                            Button::new("copy-history-details")
                                .small()
                                .label("Copy details")
                                .on_click(cx.listener(|this, _, _, cx| {
                                    cx.write_to_clipboard(ClipboardItem::new_string(
                                        this.history_details.read(cx).value().to_string(),
                                    ))
                                })),
                        ),
                )
                .child(
                    Textarea::new(&self.history_details)
                        .readonly(true)
                        .h(px(240.))
                        .aria_label("Original procedure, oligo and parameter details"),
                );
        }
        if !self.error.is_empty() {
            body = body.child(div().text_color(rgb(0xf29a8b)).child(self.error.clone()));
        }
        body.child(
            div()
                .h_flex()
                .justify_end()
                .gap_2()
                .child(
                    Button::new("cancel-metadata")
                        .label("Cancel")
                        .on_click(cx.listener(|_, _, _, cx| cx.emit(MetadataEvent::Cancel))),
                )
                .child(
                    Button::new("apply-references")
                        .primary()
                        .label("Apply references")
                        .on_click(cx.listener(|this, _, _, cx| this.save(cx))),
                ),
        )
    }
}
