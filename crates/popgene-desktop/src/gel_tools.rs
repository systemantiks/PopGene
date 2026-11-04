//! Agarose gel simulation across open sequences, ladders and saved gels.
use crate::{
    digest_tools::initial_methylation,
    figure::FigureFormat,
    gel_view::{self, GelClick, GelControls, GelControlsEvent, GelDefaults},
};
use gpui_kit::component::{
    Disableable, Selectable, Sizable, StyledExt,
    button::{Button, ButtonVariants},
    input::{Input, InputEvent, InputState},
};
use gpui_kit::{prelude::*, *};
use popgene_core::{
    Document,
    gel::{self, GEL_EXTENSION, Gel, GelSpec, LADDERS, Lane, LaneSpec, MAX_LANES, SampleDocument},
    restriction::{ENZYMES, Methylation},
};
use std::rc::Rc;

const MAX_GEL_FILE: u64 = 256 * 1024 * 1024;

pub enum GelEvent {
    Export {
        svg: String,
        format: FigureFormat,
        scale: f32,
        name: String,
    },
    Print {
        svg: String,
    },
    Save {
        json: String,
        name: String,
    },
    /// A file chooser opened by the tool is visible.
    Dialog(bool),
    Cancel,
}
struct Entry {
    spec: LaneSpec,
    /// None while a digestion runs.
    lane: Option<Result<Lane, String>>,
    generation: u64,
}
pub struct GelTools {
    name: String,
    documents: Vec<Document>,
    entries: Vec<Entry>,
    selected: Option<usize>,
    band: Option<(usize, usize)>,
    controls: Entity<GelControls>,
    label: Entity<InputState>,
    enzymes: Entity<InputState>,
    search: Entity<InputState>,
    sizes: Entity<InputState>,
    custom_sizes: Entity<InputState>,
    default_ladder: String,
    inverted: bool,
    generation: u64,
    dialog: bool,
    message: String,
    error: String,
    _subscriptions: Vec<Subscription>,
}
impl EventEmitter<GelEvent> for GelTools {}

/// An uncut lane for a document, with its recorded methylation as the assumption.
pub fn sample_lane(document: &Document) -> Result<LaneSpec, String> {
    Ok(LaneSpec::Sample {
        label: document.name.clone(),
        sample: Box::new(SampleDocument::from_document(document)?),
        enzymes: vec![],
        methylation: initial_methylation(document),
    })
}
pub fn default_lanes(defaults: &GelDefaults, documents: &[&Document]) -> Vec<LaneSpec> {
    std::iter::once(LaneSpec::Ladder {
        name: defaults.ladder.clone(),
    })
    .chain(documents.iter().filter_map(|d| sample_lane(d).ok()))
    .take(MAX_LANES)
    .collect()
}
/// Canonical catalogue names; an empty list loads uncut DNA.
fn enzyme_names(text: &str) -> Result<Vec<String>, String> {
    let names = text
        .split(|c: char| c == ',' || c == ';' || c.is_whitespace())
        .filter(|s| !s.is_empty())
        .map(str::to_string)
        .collect::<Vec<_>>();
    Ok(LaneSpec::enzyme_indices(&names)?
        .into_iter()
        .map(|i| ENZYMES[i].name.to_string())
        .collect())
}

impl GelTools {
    pub fn new(
        name: String,
        documents: Vec<Document>,
        lanes: Vec<LaneSpec>,
        defaults: GelDefaults,
        window: &mut Window,
        cx: &mut Context<Self>,
    ) -> Self {
        let controls = cx.new(|cx| GelControls::new(defaults.settings, window, cx));
        let input = |placeholder: &'static str, window: &mut Window, cx: &mut Context<Self>| {
            cx.new(|cx| InputState::new(window, cx).placeholder(placeholder))
        };
        let label = input("Lane label", window, cx);
        let enzymes = input(
            "Enzymes, comma-separated; empty loads uncut DNA",
            window,
            cx,
        );
        let search = input("Find an enzyme", window, cx);
        let sizes = input("Ladder sizes in bp", window, cx);
        let custom_sizes = input(
            "Custom ladder sizes in bp, e.g. 10000 5000 1000",
            window,
            cx,
        );
        let mut subscriptions =
            vec![cx.subscribe(&controls, |_, _, _: &GelControlsEvent, cx| cx.notify())];
        for input in [&search, &label, &enzymes, &sizes, &custom_sizes] {
            subscriptions.push(cx.subscribe(input, |_, _, _: &InputEvent, cx| cx.notify()));
        }
        let mut tools = Self {
            name,
            documents,
            entries: vec![],
            selected: None,
            band: None,
            controls,
            label,
            enzymes,
            search,
            sizes,
            custom_sizes,
            default_ladder: defaults.ladder,
            inverted: false,
            generation: 0,
            dialog: false,
            message: String::new(),
            error: String::new(),
            _subscriptions: subscriptions,
        };
        for spec in lanes {
            tools.push(spec, cx);
        }
        tools
    }
    fn push(&mut self, spec: LaneSpec, cx: &mut Context<Self>) -> bool {
        if self.entries.len() >= MAX_LANES {
            self.error = format!("A gel holds at most {MAX_LANES} lanes");
            return false;
        }
        self.entries.push(Entry {
            spec,
            lane: None,
            generation: 0,
        });
        self.compute(self.entries.len() - 1, cx);
        true
    }
    fn compute(&mut self, i: usize, cx: &mut Context<Self>) {
        self.generation += 1;
        let generation = self.generation;
        let entry = &mut self.entries[i];
        entry.generation = generation;
        let spec = entry.spec.clone();
        if spec.is_ladder() {
            entry.lane = Some(spec.lane());
            cx.notify();
            return;
        }
        entry.lane = None;
        cx.spawn(async move |this, cx| {
            let lane = cx
                .background_executor()
                .spawn(async move { spec.lane() })
                .await;
            let _ = this.update(cx, |this, cx| {
                // Lanes may be moved or removed while their digestion runs.
                if let Some(entry) = this.entries.iter_mut().find(|e| e.generation == generation) {
                    entry.lane = Some(lane);
                    cx.notify();
                }
            });
        })
        .detach();
        cx.notify();
    }
    fn gel(&self, cx: &App) -> Result<Option<Gel>, String> {
        let mut lanes = vec![];
        for (i, entry) in self.entries.iter().enumerate() {
            match &entry.lane {
                None => return Ok(None),
                Some(Err(e)) => return Err(format!("Lane {}: {e}", i + 1)),
                Some(Ok(lane)) => lanes.push(lane.clone()),
            }
        }
        if lanes.is_empty() {
            return Err("Add a ladder or a sequence lane".into());
        }
        let settings = self.controls.read(cx).settings(cx)?;
        gel::simulate(&lanes, &settings).map(Some)
    }
    fn ready_gel(&mut self, cx: &mut Context<Self>) -> Option<Gel> {
        match self.gel(cx) {
            Ok(Some(gel)) => Some(gel),
            Ok(None) => {
                self.error = "Wait for the lanes to finish calculating".into();
                cx.notify();
                None
            }
            Err(e) => {
                self.error = e;
                cx.notify();
                None
            }
        }
    }
    fn select(&mut self, i: usize, window: &mut Window, cx: &mut Context<Self>) {
        let Some(entry) = self.entries.get(i) else {
            return;
        };
        let (enzymes, sizes) = match &entry.spec {
            LaneSpec::Sample { enzymes, .. } => (enzymes.join(", "), String::new()),
            LaneSpec::CustomLadder { sizes, .. } => (
                String::new(),
                sizes
                    .iter()
                    .map(|s| s.to_string())
                    .collect::<Vec<_>>()
                    .join(" "),
            ),
            LaneSpec::Ladder { .. } => (String::new(), String::new()),
        };
        let label = entry.spec.label().to_string();
        self.label
            .update(cx, |input, cx| input.set_value(label, window, cx));
        self.enzymes
            .update(cx, |input, cx| input.set_value(enzymes, window, cx));
        self.sizes
            .update(cx, |input, cx| input.set_value(sizes, window, cx));
        self.selected = Some(i);
        self.error.clear();
        cx.notify();
    }
    /// Applies the label, enzymes or ladder sizes typed for the selected lane.
    fn apply(&mut self, cx: &mut Context<Self>) {
        let Some(i) = self.selected.filter(|&i| i < self.entries.len()) else {
            return;
        };
        let text = self.label.read(cx).value().trim().to_string();
        let result: Result<LaneSpec, String> = (|| {
            Ok(match &self.entries[i].spec {
                LaneSpec::Ladder { name } => LaneSpec::Ladder { name: name.clone() },
                LaneSpec::CustomLadder { .. } => LaneSpec::CustomLadder {
                    label: if text.is_empty() {
                        "Custom ladder".into()
                    } else {
                        text
                    },
                    sizes: gel_view::parse_ladder(&self.sizes.read(cx).value())?,
                },
                LaneSpec::Sample {
                    sample,
                    methylation,
                    ..
                } => LaneSpec::Sample {
                    label: if text.is_empty() {
                        sample.name.clone()
                    } else {
                        text
                    },
                    sample: sample.clone(),
                    enzymes: enzyme_names(&self.enzymes.read(cx).value())?,
                    methylation: *methylation,
                },
            })
        })();
        match result {
            Ok(spec) => {
                self.entries[i].spec = spec;
                self.error.clear();
                self.band = None;
                self.compute(i, cx);
            }
            Err(e) => {
                self.error = e;
                cx.notify();
            }
        }
    }
    fn apply_to_all(&mut self, cx: &mut Context<Self>) {
        let names = match enzyme_names(&self.enzymes.read(cx).value()) {
            Ok(names) => names,
            Err(e) => {
                self.error = e;
                cx.notify();
                return;
            }
        };
        self.band = None;
        for i in 0..self.entries.len() {
            if let LaneSpec::Sample { enzymes, .. } = &mut self.entries[i].spec {
                *enzymes = names.clone();
                self.compute(i, cx);
            }
        }
        self.error.clear();
    }
    fn toggle_methylation(&mut self, which: usize, cx: &mut Context<Self>) {
        let Some(i) = self.selected else { return };
        if let Some(LaneSpec::Sample { methylation, .. }) =
            self.entries.get_mut(i).map(|e| &mut e.spec)
        {
            match which {
                0 => methylation.dam = !methylation.dam,
                1 => methylation.dcm = !methylation.dcm,
                _ => methylation.cpg = !methylation.cpg,
            }
            self.band = None;
            self.compute(i, cx);
        }
    }
    fn remove(&mut self, cx: &mut Context<Self>) {
        if let Some(i) = self.selected.filter(|&i| i < self.entries.len()) {
            self.entries.remove(i);
            self.selected = None;
            self.band = None;
            cx.notify();
        }
    }
    fn shift(&mut self, right: bool, cx: &mut Context<Self>) {
        let Some(i) = self.selected else { return };
        let j = if right { i + 1 } else { i.wrapping_sub(1) };
        if j < self.entries.len() {
            self.entries.swap(i, j);
            self.selected = Some(j);
            self.band = None;
            cx.notify();
        }
    }
    fn add_custom_ladder(&mut self, cx: &mut Context<Self>) {
        match gel_view::parse_ladder(&self.custom_sizes.read(cx).value()) {
            Ok(sizes) => {
                if self.push(
                    LaneSpec::CustomLadder {
                        label: "Custom ladder".into(),
                        sizes,
                    },
                    cx,
                ) {
                    self.error.clear();
                }
            }
            Err(e) => self.error = e,
        }
        cx.notify();
    }
    fn save(&mut self, cx: &mut Context<Self>) {
        let result = self.controls.read(cx).settings(cx).and_then(|settings| {
            GelSpec::new(
                settings,
                self.entries.iter().map(|e| e.spec.clone()).collect(),
            )
            .to_json()
        });
        match result {
            Ok(json) => cx.emit(GelEvent::Save {
                json,
                name: self.name.clone(),
            }),
            Err(e) => {
                self.error = e;
                cx.notify();
            }
        }
    }
    fn open(&mut self, window: &mut Window, cx: &mut Context<Self>) {
        if self.dialog {
            return;
        }
        self.dialog = true;
        cx.emit(GelEvent::Dialog(true));
        cx.notify();
        let dialog = rfd::AsyncFileDialog::new()
            .set_parent(window)
            .add_filter("PopGene gel", &[GEL_EXTENSION]);
        cx.spawn_in(window, async move |this, cx| {
            let result = match dialog.pick_file().await {
                Some(file) => {
                    let path = file.path().to_path_buf();
                    Some(
                        cx.background_executor()
                            .spawn(async move {
                                let size =
                                    std::fs::metadata(&path).map_err(|e| e.to_string())?.len();
                                if size > MAX_GEL_FILE {
                                    return Err("The gel file exceeds the 256 MiB limit".into());
                                }
                                let text =
                                    std::fs::read_to_string(&path).map_err(|e| e.to_string())?;
                                let name = path
                                    .file_stem()
                                    .unwrap_or_default()
                                    .to_string_lossy()
                                    .into_owned();
                                GelSpec::from_json(&text).map(|spec| (name, spec))
                            })
                            .await,
                    )
                }
                None => None,
            };
            let _ = this.update_in(cx, |this, window, cx| {
                this.dialog = false;
                cx.emit(GelEvent::Dialog(false));
                match result {
                    Some(Ok((name, spec))) => {
                        this.name = name;
                        this.entries.clear();
                        this.selected = None;
                        this.band = None;
                        this.controls
                            .update(cx, |controls, cx| controls.set(spec.settings, window, cx));
                        for lane in spec.lanes {
                            this.push(lane, cx);
                        }
                        this.error.clear();
                        this.message = format!("Opened {}", this.name);
                    }
                    Some(Err(e)) => this.error = e,
                    None => {}
                }
                cx.notify();
            });
        })
        .detach();
    }
    fn export(&mut self, format: FigureFormat, cx: &mut Context<Self>) {
        if let Some(gel) = self.ready_gel(cx) {
            cx.emit(GelEvent::Export {
                svg: gel_view::gel_svg(&gel, &self.name, self.inverted),
                format,
                scale: if format == FigureFormat::Png { 3. } else { 1. },
                name: self.name.clone(),
            });
        }
    }
    fn copy_image(&mut self, cx: &mut Context<Self>) {
        let Some(gel) = self.ready_gel(cx) else {
            return;
        };
        let svg = gel_view::gel_svg(&gel, &self.name, self.inverted);
        self.message = "Rendering gel image…".into();
        cx.spawn(async move |this, cx| {
            let png = cx
                .background_executor()
                .spawn(async move { crate::figure::encode(&svg, FigureFormat::Png, 2.) })
                .await;
            let _ = this.update(cx, |this, cx| {
                match png {
                    Ok(bytes) => {
                        cx.write_to_clipboard(ClipboardItem::new_image(&Image::from_bytes(
                            ImageFormat::Png,
                            bytes,
                        )));
                        this.message = "Gel image copied".into();
                    }
                    Err(e) => this.error = e,
                }
                cx.notify();
            });
        })
        .detach();
        cx.notify();
    }
    fn save_defaults(&mut self, cx: &mut Context<Self>) {
        let settings = match self.controls.read(cx).settings(cx) {
            Ok(settings) => settings,
            Err(e) => {
                self.error = e;
                cx.notify();
                return;
            }
        };
        let ladder = self
            .entries
            .iter()
            .find_map(|e| match &e.spec {
                LaneSpec::Ladder { name } => Some(name.clone()),
                _ => None,
            })
            .unwrap_or_else(|| self.default_ladder.clone());
        self.default_ladder = ladder.clone();
        let defaults = GelDefaults { ladder, settings };
        cx.spawn(async move |this, cx| {
            let result = cx
                .background_executor()
                .spawn(async move { defaults.save() })
                .await;
            let _ = this.update(cx, |this, cx| {
                match result {
                    Ok(()) => this.message = "Saved as the default gel settings".into(),
                    Err(e) => this.error = format!("Could not save the defaults: {e}"),
                }
                cx.notify();
            });
        })
        .detach();
    }

    fn lane_editor(&self, cx: &mut Context<Self>) -> Div {
        let mut editor = div()
            .v_flex()
            .gap_2()
            .child(div().text_sm().child("Selected lane"));
        let Some((i, entry)) = self
            .selected
            .and_then(|i| self.entries.get(i).map(|e| (i, e)))
        else {
            return editor.child(
                div()
                    .text_xs()
                    .text_color(rgb(0x96a6bc))
                    .child("Click a lane label or band to edit that lane."),
            );
        };
        editor = editor
            .child(div().text_xs().child(format!(
                "Lane {} · {}",
                i + 1,
                match &entry.spec {
                    LaneSpec::Ladder { .. } => "catalogue ladder".to_string(),
                    LaneSpec::CustomLadder { .. } => "custom ladder".into(),
                    LaneSpec::Sample { sample, .. } => format!(
                        "{} · {} bp · {}",
                        sample.name,
                        gel::format_bp(sample.sequence.len()),
                        if sample.circular {
                            "circular"
                        } else {
                            "linear"
                        }
                    ),
                }
            )))
            .when(!matches!(entry.spec, LaneSpec::Ladder { .. }), |d| {
                d.child(Input::new(&self.label))
            });
        if let LaneSpec::Sample { methylation, .. } = &entry.spec {
            let query = self.search.read(cx).value().to_ascii_lowercase();
            editor = editor
                .child(Input::new(&self.enzymes))
                .child(Input::new(&self.search))
                .when(!query.is_empty(), |d| {
                    d.child(
                        div().h_flex().flex_wrap().gap_1().children(
                            ENZYMES
                                .iter()
                                .enumerate()
                                .filter(|(_, e)| e.name.to_ascii_lowercase().contains(&query))
                                .take(12)
                                .map(|(k, e)| {
                                    Button::new(("gel-add-enzyme", k))
                                        .small()
                                        .label(e.name)
                                        .on_click(cx.listener(move |this, _, window, cx| {
                                            let mut value =
                                                this.enzymes.read(cx).value().to_string();
                                            if !value.trim().is_empty() {
                                                value.push_str(", ");
                                            }
                                            value.push_str(ENZYMES[k].name);
                                            this.enzymes.update(cx, |input, cx| {
                                                input.set_value(value, window, cx)
                                            });
                                        }))
                                }),
                        ),
                    )
                })
                .child(div().text_xs().child(
                    "Methylation assumptions (complete on both strands; default reaction otherwise)",
                ))
                .child(methylation_buttons(*methylation, cx));
        }
        if matches!(entry.spec, LaneSpec::CustomLadder { .. }) {
            editor = editor.child(Input::new(&self.sizes));
        }
        let sample = matches!(entry.spec, LaneSpec::Sample { .. });
        editor
            .child(
                div()
                    .h_flex()
                    .flex_wrap()
                    .gap_1()
                    .when(!matches!(entry.spec, LaneSpec::Ladder { .. }), |d| {
                        d.child(
                            Button::new("gel-apply-lane")
                                .small()
                                .primary()
                                .label("Apply")
                                .on_click(cx.listener(|this, _, _, cx| this.apply(cx))),
                        )
                    })
                    .when(sample, |d| {
                        d.child(
                            Button::new("gel-apply-all")
                                .small()
                                .label("Apply enzymes to all lanes")
                                .on_click(cx.listener(|this, _, _, cx| this.apply_to_all(cx))),
                        )
                        .child(
                            Button::new("gel-uncut")
                                .small()
                                .label("Uncut")
                                .on_click(cx.listener(|this, _, window, cx| {
                                    this.enzymes
                                        .update(cx, |input, cx| input.set_value("", window, cx));
                                    this.apply(cx);
                                })),
                        )
                    }),
            )
            .child(
                div()
                    .h_flex()
                    .gap_1()
                    .child(
                        Button::new("gel-lane-left")
                            .small()
                            .label("Move left")
                            .disabled(i == 0)
                            .on_click(cx.listener(|this, _, _, cx| this.shift(false, cx))),
                    )
                    .child(
                        Button::new("gel-lane-right")
                            .small()
                            .label("Move right")
                            .disabled(i + 1 >= self.entries.len())
                            .on_click(cx.listener(|this, _, _, cx| this.shift(true, cx))),
                    )
                    .child(
                        Button::new("gel-lane-remove")
                            .small()
                            .label("Remove lane")
                            .on_click(cx.listener(|this, _, _, cx| this.remove(cx))),
                    ),
            )
    }
    fn lane_adder(&self, cx: &mut Context<Self>) -> Div {
        let full = self.entries.len() >= MAX_LANES;
        div()
            .v_flex()
            .gap_2()
            .child(div().text_sm().child("Add lanes"))
            .child(div().text_xs().child("Open sequences (loaded uncut)"))
            .child(
                div()
                    .h_flex()
                    .flex_wrap()
                    .gap_1()
                    .children(self.documents.iter().enumerate().map(|(k, document)| {
                        Button::new(("gel-add-sequence", k))
                            .small()
                            .label(document.name.clone())
                            .disabled(full)
                            .on_click(cx.listener(move |this, _, _, cx| {
                                match sample_lane(&this.documents[k]) {
                                    Ok(spec) => {
                                        this.push(spec, cx);
                                    }
                                    Err(e) => this.error = e,
                                }
                                cx.notify();
                            }))
                    }))
                    .when(self.documents.is_empty(), |d| {
                        d.child(
                            div()
                                .text_xs()
                                .text_color(rgb(0x96a6bc))
                                .child("No double-stranded DNA documents are open."),
                        )
                    }),
            )
            .child(div().text_xs().child("Ladders"))
            .child(
                div()
                    .h_flex()
                    .flex_wrap()
                    .gap_1()
                    .children(LADDERS.iter().enumerate().map(|(k, ladder)| {
                        Button::new(("gel-add-ladder", k))
                            .small()
                            .label(ladder.name)
                            .disabled(full)
                            .on_click(cx.listener(move |this, _, _, cx| {
                                this.push(
                                    LaneSpec::Ladder {
                                        name: LADDERS[k].name.into(),
                                    },
                                    cx,
                                );
                            }))
                    })),
            )
            .child(
                div()
                    .h_flex()
                    .gap_1()
                    .child(div().flex_1().child(Input::new(&self.custom_sizes)))
                    .child(
                        Button::new("gel-add-custom")
                            .small()
                            .label("Add custom ladder")
                            .disabled(full)
                            .on_click(cx.listener(|this, _, _, cx| this.add_custom_ladder(cx))),
                    ),
            )
    }
}
fn methylation_buttons(methylation: Methylation, cx: &mut Context<GelTools>) -> Div {
    div().h_flex().gap_1().children(
        [
            (0usize, "Dam", methylation.dam),
            (1, "Dcm", methylation.dcm),
            (2, "CpG", methylation.cpg),
        ]
        .map(|(k, label, selected)| {
            Button::new(("gel-methylation", k))
                .small()
                .selected(selected)
                .label(label)
                .on_click(cx.listener(move |this, _, _, cx| this.toggle_methylation(k, cx)))
        }),
    )
}

impl Render for GelTools {
    fn render(&mut self, _: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        let mut figure = div().v_flex().gap_2().min_w(px(420.));
        let gel = self.gel(cx);
        match &gel {
            Ok(Some(gel)) => {
                let owner = cx.entity().downgrade();
                let click: GelClick = Rc::new(move |lane, band, window, cx| {
                    let _ = owner.update(cx, |this, cx| {
                        this.select(lane, window, cx);
                        this.band = band.map(|b| (lane, b));
                    });
                });
                let band = self.band;
                figure = figure
                    .child(gel_view::gel_figure(
                        gel,
                        self.selected,
                        &|l, b| band == Some((l, b)),
                        Some(click),
                    ))
                    .children(
                        self.band
                            .and_then(|(l, b)| {
                                gel.lanes.get(l).and_then(|lane| {
                                    lane.bands.get(b).map(|band| {
                                        format!(
                                            "Lane {} · {}: {}",
                                            l + 1,
                                            lane.label,
                                            gel_view::band_description(band)
                                        )
                                    })
                                })
                            })
                            .map(|text| div().text_sm().child(text)),
                    )
                    .child(
                        div()
                            .text_xs()
                            .text_color(rgb(0x96a6bc))
                            .child(gel::MODEL_NOTE),
                    )
                    .children(
                        gel.lanes
                            .iter()
                            .filter(|l| !l.noncutters.is_empty())
                            .map(|l| {
                                div().text_xs().child(format!(
                                    "{}: no site for {}",
                                    l.label,
                                    l.noncutters.join(", ")
                                ))
                            }),
                    )
                    .children(
                        gel.warnings
                            .iter()
                            .map(|w| div().text_xs().text_color(rgb(0xe5b879)).child(w.clone())),
                    );
            }
            Ok(None) => {
                figure = figure.child(div().text_sm().child("Calculating lanes…"));
            }
            Err(e) => {
                figure = figure.child(div().text_sm().text_color(rgb(0xef938a)).child(e.clone()));
            }
        }
        let ladder_notes = self
            .entries
            .iter()
            .filter_map(|e| match &e.spec {
                LaneSpec::Ladder { name } => gel::ladder(name),
                _ => None,
            })
            .fold(Vec::<&gel::Ladder>::new(), |mut seen, ladder| {
                if !seen.contains(&ladder) {
                    seen.push(ladder);
                }
                seen
            });
        figure = figure.children(ladder_notes.into_iter().map(|ladder| {
            div().text_xs().text_color(rgb(0x96a6bc)).child(format!(
                "{}: {} Source: {}",
                ladder.name, ladder.note, ladder.source
            ))
        }));
        let ready = matches!(gel, Ok(Some(_)));
        let actions = div()
            .h_flex()
            .flex_wrap()
            .gap_2()
            .child(
                Button::new("gel-open")
                    .label("Open gel…")
                    .disabled(self.dialog)
                    .on_click(cx.listener(|this, _, window, cx| this.open(window, cx))),
            )
            .child(
                Button::new("gel-save")
                    .label("Save gel…")
                    .disabled(self.entries.is_empty())
                    .on_click(cx.listener(|this, _, _, cx| this.save(cx))),
            )
            .children(
                [
                    (0usize, FigureFormat::Svg, "Export SVG…"),
                    (1, FigureFormat::Pdf, "Export PDF…"),
                    (2, FigureFormat::Png, "Export PNG…"),
                ]
                .map(|(k, format, label)| {
                    Button::new(("gel-export", k))
                        .label(label)
                        .disabled(!ready)
                        .on_click(cx.listener(move |this, _, _, cx| this.export(format, cx)))
                }),
            )
            .child(
                Button::new("gel-inverted")
                    .label("Dark bands on white")
                    .selected(self.inverted)
                    .on_click(cx.listener(|this, _, _, cx| {
                        this.inverted = !this.inverted;
                        cx.notify();
                    })),
            )
            .child(
                Button::new("gel-copy-image")
                    .label("Copy image")
                    .disabled(!ready)
                    .on_click(cx.listener(|this, _, _, cx| this.copy_image(cx))),
            )
            .child(
                Button::new("gel-print")
                    .label("Print…")
                    .disabled(!ready)
                    .on_click(cx.listener(|this, _, _, cx| {
                        if let Some(gel) = this.ready_gel(cx) {
                            cx.emit(GelEvent::Print {
                                svg: gel_view::gel_svg(&gel, &this.name, true),
                            });
                        }
                    })),
            )
            .child(
                Button::new("gel-copy-report")
                    .label("Copy report")
                    .disabled(!ready)
                    .on_click(cx.listener(|this, _, _, cx| {
                        if let Some(gel) = this.ready_gel(cx) {
                            cx.write_to_clipboard(ClipboardItem::new_string(gel_view::report(
                                &gel,
                            )));
                        }
                    })),
            )
            .child(
                Button::new("gel-save-defaults")
                    .label("Save as defaults")
                    .on_click(cx.listener(|this, _, _, cx| this.save_defaults(cx))),
            )
            .child(
                Button::new("gel-close")
                    .label("Close")
                    .on_click(cx.listener(|_, _, _, cx| cx.emit(GelEvent::Cancel))),
            );
        div()
            .id("gel-tools-scroll")
            .v_flex()
            .gap_3()
            .p_5()
            .size_full()
            .overflow_y_scroll()
            .child(div().text_xl().child("Simulate agarose gel"))
            .child(div().text_sm().text_color(rgb(0x96a6bc)).child(format!(
                "{} · {} of {MAX_LANES} lanes",
                self.name,
                self.entries.len()
            )))
            .child(self.controls.clone())
            .child(
                div().h_flex().items_start().gap_5().child(figure).child(
                    div()
                        .v_flex()
                        .gap_4()
                        .flex_1()
                        .min_w(px(300.))
                        .child(self.lane_editor(cx))
                        .child(self.lane_adder(cx)),
                ),
            )
            .when(!self.message.is_empty(), |d| {
                d.child(div().text_sm().child(self.message.clone()))
            })
            .when(!self.error.is_empty(), |d| {
                d.child(
                    div()
                        .text_sm()
                        .text_color(rgb(0xef938a))
                        .child(self.error.clone()),
                )
            })
            .child(actions)
    }
}

#[cfg(test)]
mod tests {
    use super::{default_lanes, enzyme_names, sample_lane};
    use crate::gel_view::GelDefaults;
    use popgene_core::{Document, gel::LaneSpec, metadata::DnaMethylation};

    #[test]
    fn enzyme_lists_use_catalogue_names_and_allow_uncut_lanes() {
        assert_eq!(
            enzyme_names("ecori, BamHI;ecori").unwrap(),
            ["EcoRI", "BamHI"]
        );
        assert!(enzyme_names("").unwrap().is_empty());
        assert!(enzyme_names("EcoR").is_err());
    }

    #[test]
    fn new_gels_start_with_the_default_ladder_and_loadable_documents() {
        let mut plasmid = Document::new("pDemo", "GAATTCAAAA").unwrap();
        plasmid.circular = true;
        plasmid.methylation = Some(DnaMethylation {
            dam: true,
            dcm: false,
            eco_ki: false,
        });
        let mut rna = Document::new("RNA", "ACGT").unwrap();
        rna.single_stranded = true;
        let lanes = default_lanes(&GelDefaults::default(), &[&plasmid, &rna]);
        assert_eq!(lanes.len(), 2);
        assert!(matches!(&lanes[0], LaneSpec::Ladder { name } if name == "GeneRuler 1 kb"));
        let LaneSpec::Sample {
            label,
            sample,
            enzymes,
            methylation,
        } = &lanes[1]
        else {
            panic!("expected a sample lane");
        };
        assert_eq!(label, "pDemo");
        assert!(sample.circular && enzymes.is_empty());
        assert!(methylation.dam && !methylation.dcm && !methylation.cpg);
        assert!(sample_lane(&rna).is_err());
    }
}
