//! Native document creation and sequence editing forms.
use gpui_kit::component::{
    Disableable, Selectable, Sizable, StyledExt,
    button::{Button, ButtonVariants},
    checkbox::Checkbox,
    input::{Input, InputEvent, InputState, Textarea, TextareaState},
};
use gpui_kit::*;
use popgene_core::{
    Document, Molecule, Segment,
    creation::{Mode, Request},
};
use std::io::Read;

const MUTED: u32 = 0x8e9daf;
const LINE: u32 = 0x2c394a;
const SURFACE: u32 = 0x192330;
const ACCENT: u32 = 0x66d7c0;

pub enum CreatorEvent {
    Created(Box<Document>),
    Open,
    Folder,
    Cancel,
    Dirty(bool),
    Dialog(bool),
}

pub struct Creator {
    kind: Option<Molecule>,
    request: Request,
    name: Entity<InputState>,
    text: Entity<TextareaState>,
    count: Entity<InputState>,
    touched: bool,
    suspended: bool,
    import_pending: bool,
    error: String,
    _subscriptions: Vec<Subscription>,
}

#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct CreatorDraft {
    kind: Option<Molecule>,
    request: Request,
    touched: bool,
}

impl EventEmitter<CreatorEvent> for Creator {}

impl Creator {
    pub fn snapshot(&self, cx: &App) -> CreatorDraft {
        let mut request = self.request.clone();
        request.name = self.name.read(cx).value().to_string();
        request.text = self.text.read(cx).value().to_string();
        request.count = self.count.read(cx).value().to_string();
        CreatorDraft {
            kind: self.kind,
            request,
            touched: self.touched,
        }
    }

    pub fn restore(&mut self, draft: CreatorDraft, window: &mut Window, cx: &mut Context<Self>) {
        self.kind = draft.kind;
        self.name.update(cx, |s, cx| {
            s.set_value(draft.request.name.clone(), window, cx)
        });
        self.text.update(cx, |s, cx| {
            s.set_value(draft.request.text.clone(), window, cx)
        });
        self.count.update(cx, |s, cx| {
            s.set_value(draft.request.count.clone(), window, cx)
        });
        self.request = draft.request;
        self.touched = draft.touched;
        cx.notify();
    }

    pub fn new(window: &mut Window, cx: &mut Context<Self>) -> Self {
        let name = cx.new(|cx| {
            InputState::new(window, cx)
                .default_value("Untitled")
                .placeholder("Sequence name")
        });
        let text = cx.new(|cx| {
            TextareaState::new(window, cx)
                .rows(7)
                .placeholder("Paste a sequence or one FASTA record")
        });
        let count = cx.new(|cx| InputState::new(window, cx).default_value("100"));
        let subscriptions = vec![
            cx.subscribe(&name, |this, _, event, cx| {
                if matches!(event, InputEvent::Change) {
                    this.mark_dirty(cx);
                }
            }),
            cx.subscribe(&text, |this, _, event, cx| {
                if matches!(event, InputEvent::Change) {
                    this.mark_dirty(cx);
                }
            }),
            cx.subscribe(&count, |this, _, event, cx| {
                if matches!(event, InputEvent::Change) {
                    this.mark_dirty(cx);
                }
            }),
        ];
        Self {
            kind: None,
            request: Request::new(Molecule::Dna),
            name,
            text,
            count,
            touched: false,
            suspended: false,
            import_pending: false,
            error: String::new(),
            _subscriptions: subscriptions,
        }
    }

    pub fn is_dirty(&self) -> bool {
        self.touched
    }

    pub fn set_suspended(&mut self, suspended: bool, cx: &mut Context<Self>) {
        if self.suspended != suspended {
            self.suspended = suspended;
            cx.notify();
        }
    }

    fn mark_dirty(&mut self, cx: &mut Context<Self>) {
        if !self.touched {
            self.touched = true;
            cx.emit(CreatorEvent::Dirty(true));
        }
        self.error.clear();
        cx.notify();
    }

    fn locked(&self) -> bool {
        self.suspended || self.import_pending
    }

    pub fn choose(&mut self, kind: Molecule, cx: &mut Context<Self>) {
        if self.locked() {
            return;
        }
        choose_molecule(&mut self.request, kind);
        self.kind = Some(kind);
        self.error.clear();
        cx.notify();
    }

    fn create(&mut self, cx: &mut Context<Self>) {
        if self.locked() {
            return;
        }
        let mut request = self.request.clone();
        request.name = self.name.read(cx).value().to_string();
        request.text = self.text.read(cx).value().to_string();
        request.count = self.count.read(cx).value().to_string();
        if request.mode != Mode::Sequence {
            request.reverse = false;
        }
        match request.create() {
            Ok(doc) => cx.emit(CreatorEvent::Created(Box::new(doc))),
            Err(error) => {
                self.error = error.to_string();
                cx.notify();
            }
        }
    }

    fn load_text(&mut self, window: &mut Window, cx: &mut Context<Self>) {
        if self.locked() {
            return;
        }
        self.import_pending = true;
        cx.emit(CreatorEvent::Dialog(true));
        cx.notify();
        let dialog = rfd::AsyncFileDialog::new().set_parent(window);
        // Launch after the input callback has returned its GPUI entity/window leases.
        // The weak entity prevents a late result from reviving a discarded creator.
        cx.spawn_in(window, async move |this, cx| {
            let selected = dialog
                .add_filter(
                    "Sequence text",
                    &["txt", "fa", "fasta", "fna", "faa", "frna"],
                )
                .pick_file()
                .await;
            let result = if let Some(file) = selected {
                let path = file.path().to_path_buf();
                Some(
                    cx.background_executor()
                        .spawn(async move {
                            let result = read_sequence_text(&path);
                            (path, result)
                        })
                        .await,
                )
            } else {
                None
            };
            let _ = this.update_in(cx, |this, window, cx| {
                this.import_pending = false;
                if let Some((path, result)) = result {
                    match result {
                        Ok(text) => {
                            this.text
                                .update(cx, |state, cx| state.set_value(text, window, cx));
                            if this.name.read(cx).value().as_ref() == "Untitled" {
                                let name = path
                                    .file_stem()
                                    .unwrap_or_default()
                                    .to_string_lossy()
                                    .into_owned();
                                this.name
                                    .update(cx, |state, cx| state.set_value(name, window, cx));
                            }
                            this.request.mode = Mode::Sequence;
                            this.mark_dirty(cx);
                        }
                        Err(error) => this.error = error,
                    }
                }
                cx.emit(CreatorEvent::Dialog(false));
                cx.notify();
            });
        })
        .detach();
    }

    fn chooser(&self, cx: &mut Context<Self>) -> AnyElement {
        let cards = [
            (
                "new-dna",
                Molecule::Dna,
                0xeb8192,
                "A · T",
                "DNA",
                "Double-stranded or single-stranded",
            ),
            (
                "new-rna",
                Molecule::Rna,
                0x66d7a7,
                "A · U",
                "RNA",
                "Transcripts and RNA sequences",
            ),
            (
                "new-protein",
                Molecule::Protein,
                0x70bde8,
                "N — C",
                "Protein",
                "One-letter amino acid sequences",
            ),
        ]
        .into_iter()
        .map(|(id, kind, color, emblem, name, detail)| {
            Button::new(id)
                .ghost()
                .accessibility_label(format!("New {name} file"))
                .disabled(self.locked())
                .w(px(224.))
                .h(px(210.))
                .border_1()
                .border_color(rgb(LINE))
                .rounded_xl()
                .child(
                    div()
                        .v_flex()
                        .gap_4()
                        .items_start()
                        .p_4()
                        .w_full()
                        .child(
                            div()
                                .text_2xl()
                                .text_color(rgb(color))
                                .font_family("monospace")
                                .child(emblem),
                        )
                        .child(
                            div()
                                .text_lg()
                                .font_weight(FontWeight::SEMIBOLD)
                                .child(format!("New {name}")),
                        )
                        .child(
                            div()
                                .text_xs()
                                .text_color(rgb(MUTED))
                                .whitespace_normal()
                                .child(detail),
                        ),
                )
                .on_click(cx.listener(move |this, _, _, cx| this.choose(kind, cx)))
        });
        div()
            .size_full()
            .flex()
            .flex_col()
            .items_center()
            .justify_center()
            .gap_6()
            .p_8()
            .child(
                div()
                    .v_flex()
                    .items_center()
                    .gap_2()
                    .child(
                        div()
                            .text_xs()
                            .text_color(rgb(ACCENT))
                            .child("START SOMETHING NEW"),
                    )
                    .child(
                        div()
                            .text_3xl()
                            .font_weight(FontWeight::SEMIBOLD)
                            .child("New sequence"),
                    )
                    .child(
                        div()
                            .text_sm()
                            .text_color(rgb(MUTED))
                            .child("Choose a molecule. Your work stays on this device."),
                    ),
            )
            .child(div().flex().gap_4().children(cards))
            .child(
                div()
                    .flex()
                    .gap_3()
                    .items_center()
                    .child(
                        div()
                            .text_sm()
                            .text_color(rgb(MUTED))
                            .child("Already have a file?"),
                    )
                    .child(
                        Button::new("creator-open")
                            .disabled(self.locked())
                            .ghost()
                            .label("Open file…")
                            .on_click(cx.listener(|_, _, _, cx| cx.emit(CreatorEvent::Open))),
                    )
                    .child(
                        Button::new("creator-folder")
                            .disabled(self.locked())
                            .ghost()
                            .label("Open project folder…")
                            .on_click(cx.listener(|_, _, _, cx| cx.emit(CreatorEvent::Folder))),
                    ),
            )
            .child(
                div()
                    .text_xs()
                    .text_color(rgb(MUTED))
                    .child("SnapGene · GenBank / GenPept · EMBL · FASTA · Sanger traces"),
            )
            .into_any_element()
    }

    fn form(&self, kind: Molecule, cx: &mut Context<Self>) -> AnyElement {
        let mut body =
            div()
                .v_flex()
                .gap_4()
                .w_full()
                .max_w(px(760.))
                .child(
                    div()
                        .flex()
                        .items_center()
                        .justify_between()
                        .child(
                            Button::new("creator-back")
                                .disabled(self.locked())
                                .ghost()
                                .small()
                                .label("← Sequence types")
                                .on_click(cx.listener(|this, _, _, cx| {
                                    this.kind = None;
                                    cx.notify();
                                })),
                        )
                        .child(div().text_xs().text_color(rgb(MUTED)).child("NEW DOCUMENT")),
                )
                .child(
                    div()
                        .v_flex()
                        .gap_1()
                        .child(
                            div()
                                .text_2xl()
                                .font_weight(FontWeight::SEMIBOLD)
                                .child(format!("New {} sequence", kind.label())),
                        )
                        .child(div().text_sm().text_color(rgb(MUTED)).child(
                            "Enter a sequence, set its properties, then create your document.",
                        )),
                )
                .child(field(
                    "Sequence name",
                    Input::new(&self.name)
                        .readonly(self.locked())
                        .aria_label("Sequence name"),
                ))
                .child(
                    div().flex().gap_2().children(
                        [
                            (Mode::Sequence, "Enter sequence"),
                            (
                                Mode::Unknown,
                                if kind == Molecule::Protein {
                                    "Unknown residues"
                                } else {
                                    "Unknown bases"
                                },
                            ),
                            (Mode::Placeholder, "Empty placeholder"),
                        ]
                        .into_iter()
                        .enumerate()
                        .map(|(i, (mode, label))| {
                            Button::new(("creation-mode", i))
                                .disabled(self.locked())
                                .label(label)
                                .small()
                                .selected(self.request.mode == mode)
                                .on_click(cx.listener(move |this, _, _, cx| {
                                    if this.request.mode != mode {
                                        this.request.mode = mode;
                                        this.mark_dirty(cx);
                                    }
                                }))
                        }),
                    ),
                );
        body = match self.request.mode {
            Mode::Sequence => body
                .child(field(if kind == Molecule::Protein { "Amino acids · one-letter code or one FASTA record" } else { "Bases · raw sequence or one FASTA record" }, Textarea::new(&self.text).readonly(self.locked()).aria_label("Sequence text").font_family("monospace").h(px(176.))))
                .child(div().flex().items_center().justify_between()
                    .child(Button::new("read-text").small().disabled(self.locked()).label(if self.import_pending { "Reading text file…" } else { "Read text file…" }).on_click(cx.listener(|this, _, window, cx| this.load_text(window, cx))))
                    .child(Checkbox::new("ignore-numbers").disabled(self.locked()).label("Ignore line numbers").checked(self.request.ignore_numbers).on_change(cx.listener(|this, value, _, cx| { this.request.ignore_numbers = *value; this.mark_dirty(cx); }))))
                .child(hint("Whitespace is ignored. Invalid symbols and multiple FASTA records are reported.")),
            Mode::Unknown => body.child(field(if kind == Molecule::Protein { "Number of X residues" } else { "Number of N bases" }, Input::new(&self.count).readonly(self.locked()).aria_label("Sequence length").w(px(220.)))).child(hint("Choose a length from 1 to 10,000,000 symbols.")),
            Mode::Placeholder => body.child(div().rounded_lg().bg(rgb(SURFACE)).border_1().border_color(rgb(LINE)).p_4().text_sm().child("Create an empty document for notes. You can add the sequence later.")),
        };
        if kind != Molecule::Protein {
            body = body.child(
                div()
                    .flex()
                    .gap_6()
                    .child(field(
                        "Strandedness",
                        div().flex().gap_2().children(
                            [(false, "Double-stranded"), (true, "Single-stranded")]
                                .into_iter()
                                .enumerate()
                                .map(|(i, (value, label))| {
                                    Button::new(("strand", i))
                                        .disabled(self.locked())
                                        .small()
                                        .label(label)
                                        .selected(self.request.single_stranded == value)
                                        .on_click(cx.listener(move |this, _, _, cx| {
                                            if this.request.single_stranded != value {
                                                this.request.single_stranded = value;
                                                this.mark_dirty(cx);
                                            }
                                        }))
                                }),
                        ),
                    ))
                    .child(field(
                        "Topology",
                        div().flex().gap_2().children(
                            [(false, "Linear"), (true, "Circular")]
                                .into_iter()
                                .enumerate()
                                .map(|(i, (value, label))| {
                                    Button::new(("topology", i))
                                        .disabled(self.locked())
                                        .small()
                                        .label(label)
                                        .selected(self.request.circular == value)
                                        .on_click(cx.listener(move |this, _, _, cx| {
                                            if this.request.circular != value {
                                                this.request.circular = value;
                                                this.mark_dirty(cx);
                                            }
                                        }))
                                }),
                        ),
                    )),
            );
            if self.request.mode == Mode::Sequence {
                body = body.child(
                    div()
                        .v_flex()
                        .gap_3()
                        .child(
                            Checkbox::new("reverse")
                                .disabled(self.locked())
                                .label("Reverse complement on creation")
                                .checked(self.request.reverse)
                                .on_change(cx.listener(|this, value, _, cx| {
                                    this.request.reverse = *value;
                                    this.mark_dirty(cx);
                                })),
                        )
                        .child(
                            Checkbox::new("convert")
                                .disabled(self.locked())
                                .label(if kind == Molecule::Rna {
                                    "Convert T to U"
                                } else {
                                    "Convert U to T"
                                })
                                .checked(self.request.convert_tu)
                                .on_change(cx.listener(|this, value, _, cx| {
                                    this.request.convert_tu = *value;
                                    this.mark_dirty(cx);
                                })),
                        ),
                );
            }
        }
        if !self.error.is_empty() {
            body = body.child(error_box(&self.error));
        }
        body = body.child(
            div()
                .flex()
                .gap_3()
                .border_t_1()
                .border_color(rgb(LINE))
                .pt_4()
                .child(
                    Button::new("create-sequence")
                        .disabled(self.locked())
                        .primary()
                        .label("Create sequence")
                        .on_click(cx.listener(|this, _, _, cx| this.create(cx))),
                )
                .child(
                    Button::new("cancel-creation")
                        .disabled(self.locked())
                        .ghost()
                        .label("Cancel")
                        .on_click(cx.listener(|_, _, _, cx| cx.emit(CreatorEvent::Cancel))),
                ),
        );
        div()
            .id("creation-scroll")
            .size_full()
            .overflow_y_scroll()
            .flex()
            .justify_center()
            .p_8()
            .child(body)
            .into_any_element()
    }
}

impl Render for Creator {
    fn render(&mut self, _: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        match self.kind {
            Some(kind) => self.form(kind, cx),
            None => self.chooser(cx),
        }
    }
}

fn choose_molecule(request: &mut Request, molecule: Molecule) {
    if request.molecule != molecule {
        request.molecule = molecule;
        request.single_stranded = molecule == Molecule::Rna;
        request.reverse = false;
        request.convert_tu = false;
        if molecule == Molecule::Protein {
            request.circular = false;
        }
    }
}

fn read_sequence_text(path: &std::path::Path) -> Result<String, String> {
    let file = std::fs::File::open(path).map_err(|e| e.to_string())?;
    let mut bytes = Vec::new();
    file.take(32 * 1024 * 1024 + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| e.to_string())?;
    if bytes.len() > 32 * 1024 * 1024 {
        return Err("Text file exceeds 32 MiB. Use Open for annotated or binary files.".into());
    }
    String::from_utf8(bytes)
        .map_err(|_| "This is not a UTF-8 sequence text file. Use Open for binary files.".into())
}

fn field(label: &'static str, control: impl IntoElement) -> AnyElement {
    div()
        .v_flex()
        .gap_2()
        .child(div().text_sm().font_weight(FontWeight::MEDIUM).child(label))
        .child(control)
        .into_any_element()
}
fn hint(text: &'static str) -> AnyElement {
    div()
        .text_xs()
        .text_color(rgb(MUTED))
        .child(text)
        .into_any_element()
}
fn error_box(text: &str) -> AnyElement {
    div()
        .p_3()
        .rounded_md()
        .bg(rgb(0x3a232e))
        .text_sm()
        .text_color(rgb(0xf0a3b1))
        .child(text.to_owned())
        .into_any_element()
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum EditMode {
    Replace,
    Insert,
    Delete,
}

pub enum EditEvent {
    Apply {
        start: usize,
        end: usize,
        bases: String,
    },
    Cancel,
}
pub struct EditForm {
    document: Document,
    mode: EditMode,
    start: Entity<InputState>,
    end: Entity<InputState>,
    bases: Entity<TextareaState>,
    error: String,
}
impl EventEmitter<EditEvent> for EditForm {}
impl EditForm {
    pub fn error(&mut self, message: impl Into<String>, cx: &mut Context<Self>) {
        self.error = message.into();
        cx.notify();
    }
    pub fn new(
        document: Document,
        selected: Vec<Segment>,
        window: &mut Window,
        cx: &mut Context<Self>,
    ) -> Self {
        let first = selected.first().map_or(0, |p| p.start);
        let last = selected
            .first()
            .map_or(document.sequence.len().min(1), |p| p.end);
        Self {
            mode: if document.placeholder {
                EditMode::Insert
            } else {
                EditMode::Replace
            },
            document,
            start: cx.new(|cx| InputState::new(window, cx).default_value((first + 1).to_string())),
            end: cx.new(|cx| InputState::new(window, cx).default_value(last.to_string())),
            bases: cx.new(|cx| {
                TextareaState::new(window, cx)
                    .rows(5)
                    .placeholder("Sequence to insert")
            }),
            error: if selected.len() > 1 {
                "Choose one contiguous range. Compound features remain separate segments.".into()
            } else {
                String::new()
            },
        }
    }
    fn apply(&mut self, cx: &mut Context<Self>) {
        let result = validated_edit(
            &self.document,
            self.mode,
            &self.start.read(cx).value(),
            &self.end.read(cx).value(),
            &self.bases.read(cx).value(),
        );
        match result {
            Ok((start, end, bases)) => cx.emit(EditEvent::Apply { start, end, bases }),
            Err(error) => {
                self.error = error;
                cx.notify();
            }
        }
    }
}
impl Render for EditForm {
    fn render(&mut self, _: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        let mut body = div()
            .v_flex()
            .gap_4()
            .w_full()
            .child(
                div()
                    .text_xl()
                    .font_weight(FontWeight::SEMIBOLD)
                    .child("Edit sequence"),
            )
            .child(hint(
                "Coordinates are 1-based and inclusive. Edits can be undone.",
            ))
            .child(
                div().flex().gap_2().children(
                    [
                        (EditMode::Replace, "Replace"),
                        (EditMode::Insert, "Insert before"),
                        (EditMode::Delete, "Delete"),
                    ]
                    .into_iter()
                    .enumerate()
                    .map(|(i, (mode, label))| {
                        Button::new(("edit-mode", i))
                            .small()
                            .selected(self.mode == mode)
                            .label(label)
                            .on_click(cx.listener(move |this, _, _, cx| {
                                this.mode = mode;
                                this.error.clear();
                                cx.notify();
                            }))
                    }),
                ),
            )
            .child(
                div()
                    .flex()
                    .gap_4()
                    .child(field(
                        if self.mode == EditMode::Insert {
                            "Before position"
                        } else {
                            "From"
                        },
                        Input::new(&self.start)
                            .aria_label("Start position")
                            .w(px(160.)),
                    ))
                    .children((self.mode != EditMode::Insert).then(|| {
                        field(
                            "To",
                            Input::new(&self.end).aria_label("End position").w(px(160.)),
                        )
                    })),
            );
        if self.mode != EditMode::Delete {
            body = body.child(field(
                "New sequence",
                Textarea::new(&self.bases)
                    .aria_label("Replacement sequence")
                    .font_family("monospace")
                    .h(px(150.)),
            ));
        }
        if !self.error.is_empty() {
            body = body.child(error_box(&self.error));
        }
        body.child(
            div()
                .flex()
                .gap_3()
                .pt_2()
                .child(
                    Button::new("apply-edit")
                        .primary()
                        .label("Apply edit")
                        .on_click(cx.listener(|this, _, _, cx| this.apply(cx))),
                )
                .child(
                    Button::new("cancel-edit")
                        .ghost()
                        .label("Cancel")
                        .on_click(cx.listener(|_, _, _, cx| cx.emit(EditEvent::Cancel))),
                ),
        )
    }
}

fn validated_edit(
    doc: &Document,
    mode: EditMode,
    start: &str,
    end: &str,
    text: &str,
) -> Result<(usize, usize, String), String> {
    let first = start
        .trim()
        .parse::<usize>()
        .map_err(|_| "Enter a valid start position.")?;
    if first == 0 {
        return Err("Positions start at 1.".into());
    }
    let begin = first - 1;
    let end = if mode == EditMode::Insert {
        begin
    } else {
        end.trim()
            .parse::<usize>()
            .map_err(|_| "Enter a valid end position.")?
    };
    if begin > end || end > doc.sequence.len() {
        return Err("Choose a range within the sequence.".into());
    }
    let bases = if mode == EditMode::Delete {
        String::new()
    } else {
        doc.molecule.normalize(text).map_err(|e| e.to_string())?
    };
    if mode != EditMode::Delete && bases.is_empty() {
        return Err("Enter a sequence to insert, or choose Delete.".into());
    }
    let mut check = doc.clone();
    check
        .replace(begin, end, &bases)
        .map_err(|e| e.to_string())?;
    Ok((begin, end, bases))
}

pub enum AnnotationEvent {
    Apply { name: String, reverse: bool },
    Cancel,
}
pub struct AnnotationForm {
    document: Document,
    selected: Vec<Segment>,
    name: Entity<InputState>,
    reverse: bool,
    error: String,
}
impl EventEmitter<AnnotationEvent> for AnnotationForm {}
impl AnnotationForm {
    pub fn new(
        document: Document,
        selected: Vec<Segment>,
        window: &mut Window,
        cx: &mut Context<Self>,
    ) -> Self {
        Self {
            document,
            selected,
            name: cx.new(|cx| InputState::new(window, cx).placeholder("Feature name")),
            reverse: false,
            error: String::new(),
        }
    }
    pub fn error(&mut self, message: impl Into<String>, cx: &mut Context<Self>) {
        self.error = message.into();
        cx.notify();
    }
    fn apply(&mut self, cx: &mut Context<Self>) {
        let name = self.name.read(cx).value().trim().to_owned();
        let result = if name.is_empty() || name.chars().any(char::is_control) {
            Err("Enter a feature name without control characters.".into())
        } else if self.selected.len() != 1 {
            Err("Choose one contiguous sequence range.".into())
        } else {
            self.document
                .clone()
                .annotate(name.clone(), self.selected[0].clone(), self.reverse)
                .map_err(|e| e.to_string())
        };
        match result {
            Ok(()) => cx.emit(AnnotationEvent::Apply {
                name,
                reverse: self.reverse,
            }),
            Err(error) => self.error(error, cx),
        }
    }
}
impl Render for AnnotationForm {
    fn render(&mut self, _: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        let location = self
            .selected
            .iter()
            .map(|p| format!("{}–{}", p.start + 1, p.end))
            .collect::<Vec<_>>()
            .join(", ");
        let mut body = div()
            .v_flex()
            .gap_4()
            .child(
                div()
                    .text_xl()
                    .font_weight(FontWeight::SEMIBOLD)
                    .child("Add annotation"),
            )
            .child(
                div()
                    .text_sm()
                    .text_color(rgb(MUTED))
                    .child(format!("Selected positions {location}")),
            )
            .child(field(
                "Feature name",
                Input::new(&self.name).aria_label("Feature name"),
            ));
        if self.document.molecule != Molecule::Protein {
            body = body.child(
                Checkbox::new("annotation-reverse")
                    .checked(self.reverse)
                    .label("Reverse strand")
                    .on_change(cx.listener(|this, value, _, cx| {
                        this.reverse = *value;
                        cx.notify();
                    })),
            );
        }
        if !self.error.is_empty() {
            body = body.child(error_box(&self.error));
        }
        body.child(
            div()
                .flex()
                .gap_3()
                .child(
                    Button::new("add-feature")
                        .primary()
                        .label("Add feature")
                        .on_click(cx.listener(|this, _, _, cx| this.apply(cx))),
                )
                .child(
                    Button::new("cancel-feature")
                        .ghost()
                        .label("Cancel")
                        .on_click(cx.listener(|_, _, _, cx| cx.emit(AnnotationEvent::Cancel))),
                ),
        )
    }
}

pub enum ExportEvent {
    Choose(popgene_formats::SaveFormat),
    Cancel,
}
pub struct ExportForm {
    document: Document,
    format: popgene_formats::SaveFormat,
    export: bool,
    choose_format: bool,
    supported: Vec<(popgene_formats::SaveFormat, bool)>,
    warnings: Vec<String>,
}
impl EventEmitter<ExportEvent> for ExportForm {}
impl ExportForm {
    pub fn new(
        document: Document,
        format: popgene_formats::SaveFormat,
        export: bool,
        choose_format: bool,
    ) -> Self {
        let supported = popgene_formats::SaveFormat::ALL
            .into_iter()
            .map(|f| (f, f.supports(&document)))
            .collect();
        let warnings = popgene_formats::save_warnings(&document, format);
        Self {
            document,
            format,
            export,
            choose_format,
            supported,
            warnings,
        }
    }
}
impl Render for ExportForm {
    fn render(&mut self, _: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        let supported = self
            .supported
            .iter()
            .any(|(f, valid)| *f == self.format && *valid);
        div().v_flex().gap_4()
            .child(div().text_xl().font_weight(FontWeight::SEMIBOLD)
                .child(if self.export { "Export a copy" } else { "Save sequence" }))
            .child(div().text_sm().text_color(rgb(MUTED)).child(self.document.name.clone()))
            .children(self.supported.iter().copied().enumerate().filter(|(_, (f, _))| self.choose_format || *f == self.format).map(|(i, (format, supported))| {
                Button::new(("file-format", i)).label(format.label())
                    .selected(self.format == format)
                    .disabled(!supported)
                    .on_click(cx.listener(move |this, _, _, cx| { this.format = format; this.warnings = popgene_formats::save_warnings(&this.document, format); cx.notify(); }))
            }))
            .child(div().v_flex().gap_2().p_4().rounded_lg().border_1().border_color(rgb(LINE)).bg(rgb(SURFACE)).text_sm()
                .child(if self.warnings.is_empty() { "The supported contents can be saved in this format." } else { "Review before saving:" })
                .children(self.warnings.iter().map(|warning| div().child(format!("• {warning}"))))
                .child(if self.export { "The open document stays unchanged." } else { "The tab will use the saved file's contents. Undo can restore omitted information during this session. Keep the original file for information this format cannot store." }))
            .children((!supported).then(|| error_box("This format cannot save the current document. Add a sequence to an empty placeholder, or choose a supported format.")))
            .child(div().flex().gap_3()
                .child(Button::new("choose-export").primary()
                    .label(if self.export { "Choose export location…" } else if self.choose_format { "Choose save location…" } else { "Save in this format" })
                    .disabled(!supported)
                    .on_click(cx.listener(|this, _, _, cx| cx.emit(ExportEvent::Choose(this.format)))))
                .child(Button::new("cancel-export").ghost().label("Cancel")
                    .on_click(cx.listener(|_, _, _, cx| cx.emit(ExportEvent::Cancel)))))
    }
}

#[cfg(test)]
mod tests {
    use super::{EditMode, choose_molecule, validated_edit};
    use popgene_core::{Document, Molecule, creation::Request};
    #[test]
    fn changing_molecule_clears_incompatible_transformations() {
        let mut request = Request::new(Molecule::Dna);
        request.circular = true;
        request.reverse = true;
        request.convert_tu = true;
        choose_molecule(&mut request, Molecule::Protein);
        assert!(
            !request.circular
                && !request.reverse
                && !request.convert_tu
                && !request.single_stranded
        );
        choose_molecule(&mut request, Molecule::Rna);
        assert!(request.single_stranded);
    }
    #[test]
    fn edits_validate_positions_and_alphabet_before_emitting() {
        let doc = Document::new("test", "ACGT").unwrap();
        assert_eq!(
            validated_edit(&doc, EditMode::Replace, "2", "3", "TT").unwrap(),
            (1, 3, "TT".into())
        );
        assert_eq!(
            validated_edit(&doc, EditMode::Insert, "5", "ignored", "N").unwrap(),
            (4, 4, "N".into())
        );
        assert!(validated_edit(&doc, EditMode::Insert, "6", "", "N").is_err());
        assert!(validated_edit(&doc, EditMode::Replace, "0", "2", "A").is_err());
        assert!(validated_edit(&doc, EditMode::Replace, "1", "2", "U").is_err());
        assert!(validated_edit(&doc, EditMode::Delete, "1", "4", "").is_err());
        assert_eq!(&*doc.sequence, "ACGT");
    }
}
