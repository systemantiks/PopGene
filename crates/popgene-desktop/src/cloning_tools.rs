//! Explicit PCR and ordered assembly plans, with a read-only product preview.
use gpui_kit::component::{
    Disableable, Selectable, Sizable, StyledExt,
    button::{Button, ButtonVariants},
    input::{Input, InputEvent, InputState},
};
use gpui_kit::*;
use popgene_core::{
    Document, Molecule,
    cloning::{self, AssemblyInput, LigationOptions, PcrOptions},
    restriction::{DigestConditions, ENZYMES, Methylation},
};
use std::collections::{BTreeMap, BTreeSet};

pub enum CloningEvent {
    Create(Box<Document>),
    Cancel,
}
#[derive(Clone, Copy, PartialEq)]
enum Mode {
    Pcr,
    Ligation,
    Overlap,
    GoldenGate,
}
enum Calculated {
    Product(Box<Document>),
    Fragments(Vec<Vec<Document>>),
}
#[derive(Default)]
struct Calculation {
    generation: u64,
    pending: Option<u64>,
    candidates: Option<Vec<Vec<Document>>>,
    result: Option<Document>,
    error: String,
}
impl Calculation {
    fn invalidate(&mut self, redigest: bool) {
        self.generation = self.generation.wrapping_add(1);
        self.result = None;
        self.error.clear();
        if redigest {
            self.candidates = None;
        }
    }
    fn busy(&self) -> bool {
        self.pending.is_some()
    }
    fn begin(&mut self, redigest: bool) -> Option<u64> {
        if self.busy() {
            return None;
        }
        self.invalidate(redigest);
        self.pending = Some(self.generation);
        self.pending
    }
    fn complete(&mut self, generation: u64, result: Result<Calculated, String>) {
        // A late worker cannot publish obsolete inputs or end a newer operation.
        if self.pending != Some(generation) {
            return;
        }
        self.pending = None;
        if self.generation != generation {
            return;
        }
        match result {
            Ok(Calculated::Product(product)) => self.result = Some(*product),
            Ok(Calculated::Fragments(candidates)) => self.candidates = Some(candidates),
            Err(error) => self.error = error,
        }
    }
    fn current(&self) -> Option<&Document> {
        (!self.busy()).then_some(self.result.as_ref()).flatten()
    }
    fn fasta(&self) -> Option<String> {
        self.current()
            .map(|product| format!(">{}\n{}\n", product.name, product.sequence))
    }
}
pub struct CloningTools {
    sources: Vec<Document>,
    template: usize,
    forward: Option<u64>,
    reverse: Option<u64>,
    mode: Mode,
    order: Vec<usize>,
    orientations: BTreeSet<usize>,
    fields: Vec<(&'static str, Entity<InputState>)>,
    circular: bool,
    kinase: bool,
    allow_nicks: bool,
    forward_phosphate: bool,
    reverse_phosphate: bool,
    methylation: Methylation,
    choices: BTreeMap<usize, usize>,
    calculation: Calculation,
    _subscriptions: Vec<Subscription>,
}
impl EventEmitter<CloningEvent> for CloningTools {}
impl CloningTools {
    pub fn new(documents: Vec<Document>, window: &mut Window, cx: &mut Context<Self>) -> Self {
        let sources = documents
            .into_iter()
            .filter(|d| {
                d.molecule == Molecule::Dna
                    && !d.single_stranded
                    && d.trace.is_none()
                    && !d.sequence.is_empty()
            })
            .collect::<Vec<_>>();
        let fields = [
            ("Product name", "Cloning product"),
            ("Minimum exact overlap (bp)", "15"),
            ("Maximum exact overlap (bp)", "40"),
            ("Type IIS enzyme", "BsaI"),
            ("Temperature (°C)", "37"),
            ("pH", "7.9"),
            ("Magnesium (mM)", "10"),
            ("Glycerol (%)", "5"),
        ]
        .into_iter()
        .map(|(label, value)| {
            (
                label,
                cx.new(|cx| InputState::new(window, cx).default_value(value)),
            )
        })
        .collect::<Vec<_>>();
        let subscriptions = fields
            .iter()
            .enumerate()
            .map(|(i, (_, input))| {
                cx.subscribe(input, move |this: &mut Self, _, event: &InputEvent, cx| {
                    if matches!(event, InputEvent::Change) {
                        this.invalidate(i >= 3, cx);
                    }
                })
            })
            .collect();
        let (forward, reverse) = primers(sources.first());
        let order = if sources.is_empty() { vec![] } else { vec![0] };
        Self {
            sources,
            template: 0,
            forward,
            reverse,
            mode: Mode::Pcr,
            order,
            orientations: BTreeSet::new(),
            fields,
            circular: true,
            kinase: false,
            allow_nicks: false,
            forward_phosphate: false,
            reverse_phosphate: false,
            methylation: Methylation::default(),
            choices: BTreeMap::new(),
            calculation: Calculation::default(),
            _subscriptions: subscriptions,
        }
    }
    pub fn error(&mut self, message: impl Into<String>, cx: &mut Context<Self>) {
        self.calculation.error = message.into();
        cx.notify();
    }
    fn invalidate(&mut self, redigest: bool, cx: &mut Context<Self>) {
        self.calculation.invalidate(redigest);
        if redigest {
            self.choices.clear();
        }
        cx.notify();
    }
    fn begin(&mut self, redigest: bool, cx: &mut Context<Self>) -> Option<u64> {
        let generation = self.calculation.begin(redigest)?;
        if redigest {
            self.choices.clear();
        }
        cx.notify();
        Some(generation)
    }
    fn inputs(&self) -> Vec<AssemblyInput> {
        self.order
            .iter()
            .map(|&i| AssemblyInput {
                document: self.sources[i].clone(),
                reverse: self.orientations.contains(&i),
            })
            .collect()
    }
    fn number(&self, i: usize, cx: &App) -> Result<f64, String> {
        self.fields[i]
            .1
            .read(cx)
            .value()
            .trim()
            .parse::<f64>()
            .ok()
            .filter(|v| v.is_finite())
            .ok_or_else(|| format!("Enter a finite value for {}", self.fields[i].0))
    }
    fn overlap(&self, i: usize, cx: &App) -> Result<usize, String> {
        self.fields[i]
            .1
            .read(cx)
            .value()
            .trim()
            .parse()
            .map_err(|_| format!("Enter a whole number for {}", self.fields[i].0))
    }
    fn enzyme(&self, cx: &App) -> Result<usize, String> {
        let query = self.fields[3].1.read(cx).value();
        ENZYMES
            .iter()
            .position(|e| e.name.eq_ignore_ascii_case(query.trim()))
            .ok_or_else(|| "Enter an exact enzyme name from the catalogue".into())
    }
    fn conditions(&self, cx: &App) -> Result<DigestConditions, String> {
        Ok(DigestConditions {
            temperature_c: self.number(4, cx)?,
            ph: self.number(5, cx)?,
            magnesium_mm: self.number(6, cx)?,
            glycerol_percent: self.number(7, cx)?,
            methylation: self.methylation,
        })
    }
    fn digest(&mut self, window: &mut Window, cx: &mut Context<Self>) {
        let Some(generation) = self.begin(true, cx) else {
            return;
        };
        let parameters = self
            .enzyme(cx)
            .and_then(|e| self.conditions(cx).map(|c| (e, c)));
        let (e, c) = match parameters {
            Ok(p) => p,
            Err(e) => {
                self.calculation.complete(generation, Err(e));
                cx.notify();
                return;
            }
        };
        let inputs = self.inputs();
        cx.spawn_in(window, async move |this, cx| {
            let result = cx
                .background_executor()
                .spawn(async move { cloning::golden_gate_fragments(&inputs, e, c) })
                .await;
            let _ = this.update_in(cx, |this, _, cx| {
                this.calculation.complete(
                    generation,
                    result.map(Calculated::Fragments).map_err(|e| e.to_string()),
                );
                cx.notify();
            });
        })
        .detach();
    }
    fn preview(&mut self, window: &mut Window, cx: &mut Context<Self>) {
        let Some(generation) = self.begin(false, cx) else {
            return;
        };
        let name = self.fields[0].1.read(cx).value().trim().to_string();
        if name.is_empty() {
            self.calculation
                .complete(generation, Err("Enter a product name".into()));
            cx.notify();
            return;
        }
        let mode = self.mode;
        let inputs = self.inputs();
        let ligation = LigationOptions {
            circular: self.circular,
            phosphorylate: self.kinase,
            allow_nicks: self.allow_nicks,
        };
        let pcr_options = PcrOptions {
            forward_phosphorylated: self.forward_phosphate,
            reverse_phosphorylated: self.reverse_phosphate,
        };
        let template = self.sources.get(self.template).cloned();
        let f = self.forward;
        let r = self.reverse;
        let parameters: Result<_, String> = (|| {
            let overlap = if mode == Mode::Overlap {
                (self.overlap(1, cx)?, self.overlap(2, cx)?)
            } else {
                (15, 40)
            };
            let enzyme = if mode == Mode::GoldenGate {
                self.enzyme(cx)?
            } else {
                0
            };
            let conditions = if mode == Mode::GoldenGate {
                self.conditions(cx)?
            } else {
                DigestConditions::default()
            };
            let choices = if mode == Mode::GoldenGate {
                if self.calculation.candidates.is_none() {
                    return Err("Digest inputs and select one fragment per input first".into());
                }
                self.order
                    .iter()
                    .map(|i| {
                        self.choices
                            .get(i)
                            .copied()
                            .ok_or_else(|| "Select one digest fragment for every input".to_string())
                    })
                    .collect::<Result<Vec<_>, _>>()?
            } else {
                vec![]
            };
            Ok((overlap, enzyme, conditions, choices))
        })();
        let (overlap, enzyme, conditions, choices) = match parameters {
            Ok(p) => p,
            Err(e) => {
                self.calculation.complete(generation, Err(e));
                cx.notify();
                return;
            }
        };
        cx.spawn_in(window, async move |this, cx| {
            let result = cx
                .background_executor()
                .spawn(async move {
                    match mode {
                        Mode::Pcr => {
                            let template = template.ok_or_else(|| {
                                popgene_core::Error::Operation("Open a DNA template first".into())
                            })?;
                            let forward = template
                                .features
                                .iter()
                                .find(|feature| Some(feature.id) == f)
                                .ok_or_else(|| {
                                    popgene_core::Error::Operation(
                                        "Select a bound forward primer".into(),
                                    )
                                })?;
                            let reverse = template
                                .features
                                .iter()
                                .find(|feature| Some(feature.id) == r)
                                .ok_or_else(|| {
                                    popgene_core::Error::Operation(
                                        "Select a bound reverse primer".into(),
                                    )
                                })?;
                            cloning::pcr_product(&template, forward, reverse, pcr_options, name)
                        }
                        Mode::Ligation => cloning::ligate(&inputs, ligation, name),
                        Mode::Overlap => cloning::assemble_overlaps(
                            &inputs,
                            overlap.0,
                            overlap.1,
                            ligation.circular,
                            name,
                        ),
                        Mode::GoldenGate => {
                            cloning::golden_gate(&inputs, &choices, enzyme, conditions, name)
                        }
                    }
                })
                .await;
            let _ = this.update_in(cx, |this, _, cx| {
                this.calculation.complete(
                    generation,
                    result
                        .map(|product| Calculated::Product(Box::new(product)))
                        .map_err(|e| e.to_string()),
                );
                cx.notify();
            });
        })
        .detach();
    }
}
fn primers(document: Option<&Document>) -> (Option<u64>, Option<u64>) {
    let find = |reverse| {
        document.and_then(|d| {
            d.features
                .iter()
                .find(|f| f.kind == "primer_bind" && f.reverse == reverse && !f.parts.is_empty())
                .map(|f| f.id)
        })
    };
    (find(false), find(true))
}
fn move_source(order: &mut [usize], source: usize, up: bool) -> bool {
    let Some(position) = order.iter().position(|&i| i == source) else {
        return false;
    };
    let next = if up {
        position.checked_sub(1)
    } else {
        position.checked_add(1).filter(|&i| i < order.len())
    };
    let Some(next) = next else {
        return false;
    };
    order.swap(position, next);
    true
}
fn field(label: &'static str, input: impl IntoElement) -> Div {
    div()
        .v_flex()
        .gap_1()
        .child(div().text_xs().child(label))
        .child(input)
}
impl Render for CloningTools {
    fn render(&mut self, _: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        let mut body=div().id("cloning-tools-scroll").v_flex().gap_3().p_5().max_h(px(760.)).overflow_y_scroll()
            .child(div().text_xl().child("Cloning and PCR"))
            .child(div().text_xs().text_color(rgb(0x96a6bc)).child("Build a sequence product from captured document snapshots. Review the preview before creating a new document. Reaction yield and error rates are not predicted."))
            .child(field("Product name",Input::new(&self.fields[0].1).disabled(self.calculation.busy())))
            .child(div().h_flex().flex_wrap().gap_2().children([(Mode::Pcr,"PCR"),(Mode::Ligation,"Restriction ligation"),(Mode::Overlap,"Overlap assembly"),(Mode::GoldenGate,"Golden Gate")].into_iter().enumerate().map(|(i,(mode,label))|Button::new(("cloning-mode",i)).small().selected(self.mode==mode).disabled(self.calculation.busy()).label(label).on_click(cx.listener(move|this,_,_,cx|{if this.calculation.busy(){return;}this.mode=mode;this.invalidate(false,cx);})))));
        if self.mode == Mode::Pcr {
            body =
                body.child(div().text_sm().child("Template")).child(
                    div().h_flex().flex_wrap().gap_2().children(
                        self.sources.iter().enumerate().map(|(i, d)| {
                            Button::new(("pcr-template", i))
                                .small()
                                .selected(i == self.template)
                                .disabled(self.calculation.busy())
                                .label(d.name.clone())
                                .on_click(cx.listener(move |this, _, _, cx| {
                                    if this.calculation.busy() {
                                        return;
                                    }
                                    this.template = i;
                                    (this.forward, this.reverse) = primers(this.sources.get(i));
                                    this.invalidate(false, cx);
                                }))
                        }),
                    ),
                );
            if let Some(template) = self.sources.get(self.template) {
                let features = template
                    .features
                    .iter()
                    .filter(|f| f.kind == "primer_bind" && !f.parts.is_empty())
                    .collect::<Vec<_>>();
                if features.is_empty() {
                    body=body.child(div().text_color(rgb(0xe7bd7b)).child("Add or design bound primers in the Primers view, then reopen cloning tools."));
                }
                for (i, f) in features.iter().take(200).enumerate() {
                    let id = f.id;
                    let reverse = f.reverse;
                    body = body.child(
                        Button::new(("pcr-primer", i))
                            .small()
                            .selected(if reverse {
                                self.reverse == Some(id)
                            } else {
                                self.forward == Some(id)
                            })
                            .disabled(self.calculation.busy())
                            .label(format!(
                                "{} · {} · {}",
                                if reverse { "R" } else { "F" },
                                f.label,
                                f.location_label()
                            ))
                            .on_click(cx.listener(move |this, _, _, cx| {
                                if this.calculation.busy() {
                                    return;
                                }
                                if reverse {
                                    this.reverse = Some(id)
                                } else {
                                    this.forward = Some(id)
                                }
                                this.invalidate(false, cx);
                            })),
                    );
                }
            }
            body = body.child(
                div().h_flex().gap_2().children(
                    [
                        (
                            0usize,
                            self.forward_phosphate,
                            "Forward primer 5′ phosphate",
                        ),
                        (1, self.reverse_phosphate, "Reverse primer 5′ phosphate"),
                    ]
                    .into_iter()
                    .map(|(i, selected, label)| {
                        Button::new(("pcr-phosphate", i))
                            .small()
                            .selected(selected)
                            .disabled(self.calculation.busy())
                            .label(label)
                            .on_click(cx.listener(move |this, _, _, cx| {
                                if this.calculation.busy() {
                                    return;
                                }
                                if i == 0 {
                                    this.forward_phosphate = !this.forward_phosphate
                                } else {
                                    this.reverse_phosphate = !this.reverse_phosphate
                                }
                                this.invalidate(false, cx);
                            }))
                    }),
                ),
            );
            body=body.child(div().text_xs().text_color(rgb(0x96a6bc)).child("Stored full oligos supply 5′ tails. Their 3′ suffixes must match the declared binding intervals exactly and uniquely. The mature product has blunt ends; polymerase-specific A-tailing is not inferred."));
        } else {
            body = body
                .child(
                    div()
                        .text_sm()
                        .child("Input documents — click to include or remove"),
                )
                .child(div().h_flex().flex_wrap().gap_2().children(
                    self.sources.iter().enumerate().map(|(i, d)| {
                        Button::new(("assembly-source", i))
                            .small()
                            .selected(self.order.contains(&i))
                            .disabled(self.calculation.busy())
                            .label(d.name.clone())
                            .on_click(cx.listener(move |this, _, _, cx| {
                                if this.calculation.busy() {
                                    return;
                                }
                                if this.order.contains(&i) {
                                    this.order.retain(|&j| j != i);
                                } else {
                                    this.order.push(i);
                                }
                                this.invalidate(true, cx);
                            }))
                    }),
                ));
            for (position, &i) in self.order.iter().enumerate() {
                let total = self.order.len();
                body = body.child(
                    div()
                        .h_flex()
                        .gap_2()
                        .items_center()
                        .child(div().flex_1().child(format!(
                            "{}. {} · {} bp",
                            position + 1,
                            self.sources[i].name,
                            self.sources[i].sequence.len()
                        )))
                        .child(
                            Button::new(("assembly-reverse", i))
                                .small()
                                .selected(self.orientations.contains(&i))
                                .disabled(self.calculation.busy())
                                .label("Reverse")
                                .on_click(cx.listener(move |this, _, _, cx| {
                                    if this.calculation.busy() {
                                        return;
                                    }
                                    if !this.orientations.insert(i) {
                                        this.orientations.remove(&i);
                                    }
                                    this.invalidate(false, cx);
                                })),
                        )
                        .child(
                            Button::new(("assembly-up", i))
                                .small()
                                .disabled(self.calculation.busy() || position == 0)
                                .label("↑")
                                .on_click(cx.listener(move |this, _, _, cx| {
                                    if this.calculation.busy() {
                                        return;
                                    }
                                    if move_source(&mut this.order, i, true) {
                                        this.invalidate(true, cx);
                                    }
                                })),
                        )
                        .child(
                            Button::new(("assembly-down", i))
                                .small()
                                .disabled(self.calculation.busy() || position + 1 >= total)
                                .label("↓")
                                .on_click(cx.listener(move |this, _, _, cx| {
                                    if this.calculation.busy() {
                                        return;
                                    }
                                    if move_source(&mut this.order, i, false) {
                                        this.invalidate(true, cx);
                                    }
                                })),
                        ),
                );
            }
            if self.mode != Mode::GoldenGate {
                body = body.child(
                    Button::new("assembly-circular")
                        .small()
                        .selected(self.circular)
                        .disabled(self.calculation.busy())
                        .label("Close circular product")
                        .on_click(cx.listener(|this, _, _, cx| {
                            if this.calculation.busy() {
                                return;
                            }
                            this.circular = !this.circular;
                            this.invalidate(false, cx);
                        })),
                );
            }
            if self.mode == Mode::Ligation {
                body = body.child(
                    div().h_flex().gap_2().children(
                        [
                            (0usize, self.kinase, "Model 5′ kinase treatment"),
                            (1, self.allow_nicks, "Allow one sealed strand per junction"),
                        ]
                        .into_iter()
                        .map(|(i, selected, label)| {
                            Button::new(("ligation-chemistry", i))
                                .small()
                                .selected(selected)
                                .disabled(self.calculation.busy())
                                .label(label)
                                .on_click(cx.listener(move |this, _, _, cx| {
                                    if this.calculation.busy() {
                                        return;
                                    }
                                    if i == 0 {
                                        this.kinase = !this.kinase
                                    } else {
                                        this.allow_nicks = !this.allow_nicks
                                    }
                                    this.invalidate(false, cx);
                                }))
                        }),
                    ),
                );
                body=body.child(div().text_xs().text_color(rgb(0x96a6bc)).child("Use linear fragments, for example products from Digestion. Both cohesive polarity and sequence must match. Native ends with unspecified phosphorylation require explicit kinase treatment; unknown blocking modifications remain unsupported."));
            } else if self.mode == Mode::Overlap {
                body=body.child(div().h_flex().gap_3().children(self.fields.iter().skip(1).take(2).map(|(label,input)|field(label,Input::new(input).disabled(self.calculation.busy())))))
                    .child(div().text_xs().text_color(rgb(0x96a6bc)).child("Exact-overlap Gibson-style product model for fully paired linear fragments. Each junction must have one unambiguous overlap in the requested range; duplicate overlap sequence appears once in the product."));
            } else {
                body = body
                    .child(div().text_xs().child(format!("Supported enzymes: {}. Engineered supplier variants use the parent enzyme's recognition geometry.", cloning::GOLDEN_GATE_ENZYMES.join(", "))))
                    .child(div().h_flex().flex_wrap().gap_3().children(
                        self.fields.iter().skip(3).map(|(label, input)| {
                            field(label, Input::new(input).w(px(140.)).disabled(self.calculation.busy()))
                        }),
                    ))
                    .child(
                        div().h_flex().gap_2().children(
                            [
                                (0usize, self.methylation.dam, "Dam methylation"),
                                (1, self.methylation.dcm, "Dcm methylation"),
                                (2, self.methylation.cpg, "CpG methylation"),
                            ]
                            .into_iter()
                            .map(|(i, selected, label)| {
                                Button::new(("golden-gate-methylation", i))
                                    .small()
                                    .selected(selected)
                                    .disabled(self.calculation.busy())
                                    .label(label)
                                    .on_click(cx.listener(move |this, _, _, cx| {
                                    if this.calculation.busy() { return; }
                                        match i {
                                            0 => this.methylation.dam = !this.methylation.dam,
                                            1 => this.methylation.dcm = !this.methylation.dcm,
                                            _ => this.methylation.cpg = !this.methylation.cpg,
                                        }
                                        this.invalidate(true, cx);
                                    }))
                            }),
                        ),
                    )
                    .child(
                        Button::new("golden-gate-digest")
                            .small()
                            .disabled(self.calculation.busy())
                            .label("Digest inputs and choose fragments")
                            .on_click(cx.listener(|this, _, window, cx| this.digest(window, cx))),
                    );
                if let Some(candidates) = &self.calculation.candidates {
                    for (position, fragments) in candidates.iter().enumerate() {
                        let source = self.order[position];
                        body = body.child(div().text_sm().child(format!(
                            "{} — select one fragment",
                            self.sources[source].name
                        )));
                        for (index, fragment) in fragments.iter().enumerate() {
                            let detail = fragment.duplex_ends.as_ref().map_or(
                                "circular / no cuts".into(),
                                |e| {
                                    format!(
                                        "{:?} {} → {:?} {}",
                                        e.left.kind,
                                        e.left.overhang,
                                        e.right.kind,
                                        e.right.overhang
                                    )
                                },
                            );
                            body = body.child(
                                Button::new(("golden-gate-fragment", position * 128 + index))
                                    .small()
                                    .selected(self.choices.get(&source) == Some(&index))
                                    .disabled(self.calculation.busy())
                                    .label(format!(
                                        "{} · {} bp · {}",
                                        index + 1,
                                        fragment.sequence.len(),
                                        detail
                                    ))
                                    .on_click(cx.listener(move |this, _, _, cx| {
                                        if this.calculation.busy() {
                                            return;
                                        }
                                        this.choices.insert(source, index);
                                        this.invalidate(false, cx);
                                    })),
                            );
                        }
                    }
                }
                body=body.child(div().text_xs().text_color(rgb(0x96a6bc)).child("Fragment choice, order and orientation are explicit. The circular product must have unique complementary overhang partners and no remaining site for the selected enzyme. Enzyme cycling and ligase mismatch fidelity are not predicted."));
            }
        }
        body = body.child(
            Button::new("preview-cloning-product")
                .primary()
                .disabled(self.calculation.busy())
                .label(if self.calculation.busy() {
                    "Calculating…"
                } else {
                    "Preview product"
                })
                .on_click(cx.listener(|this, _, window, cx| this.preview(window, cx))),
        );
        if let Some(product) = self.calculation.current() {
            let procedure = product.comments.last().cloned().unwrap_or_default();
            body = body
                .child(div().text_lg().child(format!(
                    "{} · {} bp · {} · {} annotations",
                    product.name,
                    product.sequence.len(),
                    if product.circular {
                        "circular"
                    } else {
                        "linear"
                    },
                    product.features.len() + product.duplex_annotations.features.len()
                )))
                .child(div().text_xs().child(procedure))
                .child(
                    div()
                        .h_flex()
                        .gap_2()
                        .child(
                            Button::new("copy-cloning-fasta")
                                .small()
                                .label("Copy product FASTA")
                                .on_click(cx.listener(|this, _, _, cx| {
                                    if let Some(fasta) = this.calculation.fasta() {
                                        cx.write_to_clipboard(ClipboardItem::new_string(fasta));
                                    }
                                })),
                        )
                        .child(
                            Button::new("create-cloning-product")
                                .primary()
                                .label("Create product document")
                                .on_click(cx.listener(|this, _, _, cx| {
                                    if let Some(product) = this.calculation.current() {
                                        cx.emit(CloningEvent::Create(Box::new(product.clone())));
                                    }
                                })),
                        ),
                );
            for warning in &product.warnings {
                body = body.child(
                    div()
                        .text_xs()
                        .text_color(rgb(0xe7bd7b))
                        .child(warning.clone()),
                );
            }
        }
        if !self.calculation.error.is_empty() {
            body = body.child(
                div()
                    .text_color(rgb(0xf2a4a4))
                    .child(self.calculation.error.clone()),
            );
        }
        body.child(
            Button::new("close-cloning-tools")
                .ghost()
                .label("Close")
                .on_click(cx.listener(|_, _, _, cx| cx.emit(CloningEvent::Cancel))),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::{Calculated, Calculation, move_source};
    use popgene_core::{
        Document, Feature, Segment,
        cloning::{self, AssemblyInput, PcrOptions},
        metadata::Reference,
        restriction::{DigestConditions, ENZYMES},
    };

    fn pcr(name: &str, tail: &str) -> Document {
        let mut template = Document::new("template", "AACCGTTAGCTAGGTCACGA").unwrap();
        template.references.push(Reference {
            title: format!("{name} publication"),
            ..Default::default()
        });
        template
            .annotate("insert".into(), Segment { start: 5, end: 9 }, false)
            .unwrap();
        template.features[0]
            .qualifiers
            .push(("citation".into(), Some("[1]".into())));
        let primer = |id, reverse, start, sequence: String| Feature {
            id,
            label: if reverse { "reverse" } else { "forward" }.into(),
            kind: "primer_bind".into(),
            parts: vec![Segment {
                start,
                end: start + 4,
            }],
            location: None,
            reverse,
            qualifiers: vec![("primer_sequence".into(), Some(sequence))],
            edited: false,
        };
        cloning::pcr_product(
            &template,
            &primer(2, false, 0, format!("{tail}AACC")),
            &primer(3, true, 16, "TCGT".into()),
            PcrOptions::default(),
            name,
        )
        .unwrap()
    }

    fn golden_gate_inputs() -> (Vec<AssemblyInput>, usize) {
        let inputs = [
            ("vector", "AAAAGGTCTCAAGTCCCGATATACCGGCGTATGAGACCAAAA"),
            ("insert", "AAAAGGTCTCACGTATTAAGGCCATAGTCTGAGACCAAAA"),
        ]
        .into_iter()
        .map(|(name, sequence)| AssemblyInput {
            document: Document::new(name, sequence).unwrap(),
            reverse: false,
        })
        .collect();
        let enzyme = ENZYMES.iter().position(|e| e.name == "BsaI").unwrap();
        (inputs, enzyme)
    }

    fn publish(state: &mut Calculation, product: Document) {
        let token = state.begin(false).unwrap();
        state.complete(token, Ok(Calculated::Product(Box::new(product))));
    }

    #[test]
    fn queued_reorder_resolves_the_source_after_removal_or_another_move() {
        // Rendered source 2 used to be in position 2; source 1 was removed
        // before its Up callback arrived. The old position would panic.
        let mut order = vec![0, 2];
        assert!(move_source(&mut order, 2, true));
        assert_eq!(order, [2, 0]);
        assert!(!move_source(&mut order, 1, true));
        assert!(!move_source(&mut order, 2, true));
        assert!(!move_source(&mut order, 0, false));

        // A previous move put source 1 first; its queued Down action must move
        // source 1, rather than whichever source now occupies its old row.
        let mut order = vec![1, 2, 0];
        assert!(move_source(&mut order, 1, false));
        assert_eq!(order, [2, 1, 0]);
        assert!(!move_source(&mut [], 1, false));
    }

    #[test]
    fn actions_use_only_the_current_product_and_its_source_references() {
        let mut state = Calculation::default();
        publish(&mut state, pcr("first", "GGG"));
        assert_eq!(
            state.fasta().as_deref(),
            Some(">first\nGGGAACCGTTAGCTAGGTCACGA\n")
        );
        // This is the same invalidation used by every input callback. Even an
        // already-rendered Create/Copy callback now reads this cleared state.
        state.invalidate(false);
        assert!(state.current().is_none());
        assert!(state.fasta().is_none());
        publish(&mut state, pcr("second", "TT"));
        let product = state.current().unwrap();
        assert_eq!(&*product.sequence, "TTAACCGTTAGCTAGGTCACGA");
        assert_eq!(product.references.len(), 1);
        assert_eq!(product.references[0].title, "second publication");
        let feature = product
            .features
            .iter()
            .find(|f| f.label == "insert")
            .unwrap();
        assert_eq!(feature.qualifier("citation"), Some("[1]"));
        assert_eq!(feature.parts, [Segment { start: 7, end: 11 }]);
        assert_eq!(
            state.fasta().as_deref(),
            Some(">second\nTTAACCGTTAGCTAGGTCACGA\n")
        );
    }

    #[test]
    fn pending_and_failed_previews_cannot_reuse_a_successful_product() {
        let mut state = Calculation::default();
        publish(&mut state, pcr("previous", ""));
        // begin precedes validation, including an empty product name or a
        // malformed numeric field that never starts a worker.
        let token = state.begin(false).unwrap();
        assert!(state.busy());
        assert!(state.current().is_none());
        assert!(state.fasta().is_none());
        assert!(state.begin(true).is_none());
        state.complete(token, Err("Enter a product name".into()));
        assert!(!state.busy());
        assert!(state.current().is_none());
        assert_eq!(state.error, "Enter a product name");

        publish(&mut state, pcr("another previous", ""));
        let token = state.begin(false).unwrap();
        let (inputs, _) = golden_gate_inputs();
        let failed = cloning::assemble_overlaps(&inputs, 40, 15, false, "invalid overlaps");
        assert!(failed.is_err());
        state.complete(
            token,
            failed
                .map(|d| Calculated::Product(Box::new(d)))
                .map_err(|e| e.to_string()),
        );
        assert!(!state.error.is_empty());
        assert!(state.current().is_none());
        assert!(state.fasta().is_none());
    }

    #[test]
    fn obsolete_product_completion_cannot_publish_or_end_a_newer_digest() {
        let mut state = Calculation::default();
        let old = state.begin(false).unwrap();
        // A delayed InputEvent::Change can invalidate an in-flight operation.
        state.invalidate(true);
        state.complete(old, Ok(Calculated::Product(Box::new(pcr("obsolete", "")))));
        assert!(!state.busy());
        assert!(state.current().is_none());
        let next = state.begin(true).unwrap();
        state.complete(old, Err("late failure".into()));
        assert!(state.busy());
        assert_eq!(state.pending, Some(next));
        assert!(state.error.is_empty());
        let (inputs, enzyme) = golden_gate_inputs();
        let candidates =
            cloning::golden_gate_fragments(&inputs, enzyme, DigestConditions::default()).unwrap();
        state.complete(next, Ok(Calculated::Fragments(candidates)));
        assert!(!state.busy());
        let candidates = state.candidates.as_ref().unwrap();
        assert_eq!(candidates.len(), 2);
        assert_eq!(candidates[0][1].sequence.as_ref(), "AGTCCCGATATACCGG");
        assert_eq!(candidates[1][1].sequence.as_ref(), "CGTATTAAGGCCAT");
        assert!(state.current().is_none());
    }

    #[test]
    fn digest_candidates_survive_product_only_changes_but_not_redigestion() {
        let mut state = Calculation::default();
        let (inputs, enzyme) = golden_gate_inputs();
        let digest = state.begin(true).unwrap();
        state.complete(
            digest,
            Ok(Calculated::Fragments(
                cloning::golden_gate_fragments(&inputs, enzyme, DigestConditions::default())
                    .unwrap(),
            )),
        );
        let candidates = state.candidates.as_ref().unwrap().as_ptr();
        let product = cloning::golden_gate(
            &inputs,
            &[1, 1],
            enzyme,
            DigestConditions::default(),
            "Golden Gate",
        )
        .unwrap();
        publish(&mut state, product);
        assert_eq!(
            state.current().unwrap().sequence.as_ref(),
            "AGTCCCGATATACCGGCGTATTAAGGCCAT"
        );
        state.invalidate(false);
        assert_eq!(state.candidates.as_ref().unwrap().as_ptr(), candidates);
        assert!(state.current().is_none());
        let old = state.begin(true).unwrap();
        assert!(state.candidates.is_none());
        state.invalidate(true);
        state.complete(
            old,
            Ok(Calculated::Fragments(
                cloning::golden_gate_fragments(&inputs, enzyme, DigestConditions::default())
                    .unwrap(),
            )),
        );
        assert!(state.candidates.is_none());
        let next = state.begin(false).unwrap();
        state.complete(old, Err("obsolete digest failure".into()));
        assert_eq!(state.pending, Some(next));
        state.complete(
            next,
            Err("Digest inputs and select one fragment per input first".into()),
        );
        assert!(!state.busy());
        assert!(state.candidates.is_none());
        assert!(state.fasta().is_none());
    }
}
