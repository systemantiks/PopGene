//! Native Primer3 analysis/design and oligo editing. The parent applies events
//! through its undoable document model after checking the captured revision.
use gpui_kit::component::{
    Disableable, Selectable, Sizable, StyledExt,
    button::{Button, ButtonVariants},
    input::{Input, InputEvent, InputState, Textarea, TextareaState},
};
use gpui_kit::*;
use popgene_core::{
    Document, Feature, Molecule, Segment,
    primer::{
        self, DesignOptions, DesignResult, MeltingModel, OligoAnalysis, PrimerPair, Reaction,
        SaltCorrection,
    },
};

pub enum PrimerEvent {
    Save {
        feature_id: Option<u64>,
        name: String,
        sequence: String,
        bindings: Vec<Segment>,
        reverse: bool,
        reaction: Reaction,
        tm_c: Option<f64>,
    },
    AddPair {
        pair: PrimerPair,
        reaction: Reaction,
    },
    Cancel,
}
#[derive(Clone, Copy, PartialEq)]
enum Mode {
    Oligo,
    Design,
}
pub struct PrimerTools {
    document: Document,
    editing: Option<u64>,
    mode: Mode,
    name: Entity<InputState>,
    sequence: Entity<TextareaState>,
    bindings: Entity<InputState>,
    reverse: bool,
    fields: Vec<(&'static str, Entity<InputState>)>,
    melting_model: MeltingModel,
    salt: SaltCorrection,
    busy: bool,
    error: String,
    analysis: Option<OligoAnalysis>,
    design: Option<DesignResult>,
    result_reaction: Option<Reaction>,
    _subscriptions: Vec<Subscription>,
}
impl EventEmitter<PrimerEvent> for PrimerTools {}
impl PrimerTools {
    pub fn new(
        document: Document,
        selected: Vec<Segment>,
        editing: Option<Feature>,
        window: &mut Window,
        cx: &mut Context<Self>,
    ) -> Self {
        let stored_reaction = editing
            .as_ref()
            .and_then(|f| f.qualifier("primer_reaction"));
        let decoded_reaction = stored_reaction
            .map(serde_json::from_str::<Reaction>)
            .transpose();
        let reaction_error = decoded_reaction.as_ref().err().map(|_| "Stored reaction parameters could not be read; review the displayed defaults before applying edits".to_string()).unwrap_or_default();
        let reaction = decoded_reaction.ok().flatten().unwrap_or_default();
        let reaction_values = [
            reaction.dna_nm,
            reaction.monovalent_mm,
            reaction.magnesium_mm,
            reaction.dntp_mm,
            reaction.dmso_percent,
            reaction.formamide_m,
            reaction.temperature_c,
        ];
        let text = editing
            .as_ref()
            .and_then(Feature::primer_sequence)
            .unwrap_or_else(|| {
                selected
                    .iter()
                    .filter_map(|p| document.sequence.get(p.start..p.end))
                    .collect()
            });
        let parts = editing.as_ref().map_or(selected, |f| f.parts.clone());
        let design_target = if editing.is_none() {
            match parts.as_slice() {
                [part] => Some((part.start + 1, part.end)),
                [tail, head]
                    if document.circular
                        && tail.end == document.sequence.len()
                        && head.start == 0 =>
                {
                    Some((tail.start + 1, head.end))
                }
                _ => None,
            }
        } else {
            None
        };
        let binding_text = parts
            .iter()
            .map(|p| format!("{}..{}", p.start + 1, p.end))
            .collect::<Vec<_>>()
            .join(", ");
        let name = cx.new(|cx| {
            InputState::new(window, cx)
                .default_value(editing.as_ref().map_or("Primer", |f| &f.label).to_string())
        });
        let sequence = cx.new(|cx| TextareaState::new(window, cx).rows(3).default_value(text));
        let bindings = cx.new(|cx| InputState::new(window, cx).default_value(binding_text));
        let fields = [
            ("DNA (nM)", "50"),
            ("Monovalent ions (mM)", "50"),
            ("Magnesium (mM)", "1.5"),
            ("dNTPs (mM)", "0.6"),
            ("DMSO (%)", "0"),
            ("Formamide (M)", "0"),
            ("Structure temperature (°C)", "37"),
            ("Primer minimum length", "18"),
            ("Primer optimum length", "20"),
            ("Primer maximum length", "27"),
            ("Minimum Tm (°C)", "57"),
            ("Optimum Tm (°C)", "60"),
            ("Maximum Tm (°C)", "63"),
            ("Minimum GC (%)", "20"),
            ("Maximum GC (%)", "80"),
            ("Maximum Tm difference (°C)", "3"),
            ("Minimum product length", "100"),
            ("Maximum product length", "300"),
            ("Pairs to return", "5"),
            ("Target from (optional)", ""),
            ("Target to (optional)", ""),
        ]
        .into_iter()
        .enumerate()
        .map(|(index, (label, value))| {
            let mut value = reaction_values
                .get(index)
                .map_or_else(|| value.to_string(), |v| v.to_string());
            if let Some((from, to)) = design_target {
                if index == 19 {
                    value = from.to_string();
                } else if index == 20 {
                    value = to.to_string();
                }
            }
            (
                label,
                cx.new(|cx| InputState::new(window, cx).default_value(value)),
            )
        })
        .collect::<Vec<_>>();
        let mut subscriptions = vec![];
        for input in [name.clone(), bindings.clone()]
            .into_iter()
            .chain(fields.iter().map(|(_, input)| input.clone()))
        {
            subscriptions.push(cx.subscribe(&input, |this, _, event: &InputEvent, cx| {
                if matches!(event, InputEvent::Change) {
                    this.invalidate(cx);
                }
            }));
        }
        subscriptions
            .push(cx.subscribe(&sequence, |this, _, _: &InputEvent, cx| this.invalidate(cx)));
        Self {
            document,
            editing: editing.as_ref().map(|f| f.id),
            mode: Mode::Oligo,
            name,
            sequence,
            bindings,
            reverse: editing.is_some_and(|f| f.reverse),
            fields,
            melting_model: reaction.melting_model,
            salt: reaction.salt_correction,
            busy: false,
            error: reaction_error,
            analysis: None,
            design: None,
            result_reaction: None,
            _subscriptions: subscriptions,
        }
    }
    pub fn error(&mut self, message: impl Into<String>, cx: &mut Context<Self>) {
        self.error = message.into();
        cx.notify();
    }
    fn invalidate(&mut self, cx: &mut Context<Self>) {
        self.analysis = None;
        self.design = None;
        self.result_reaction = None;
        self.error.clear();
        cx.notify();
    }
    fn number(&self, index: usize, cx: &App) -> Result<f64, String> {
        self.fields[index]
            .1
            .read(cx)
            .value()
            .trim()
            .parse::<f64>()
            .map_err(|_| format!("Enter a number for {}", self.fields[index].0))
    }
    fn integer(&self, index: usize, cx: &App) -> Result<usize, String> {
        self.fields[index]
            .1
            .read(cx)
            .value()
            .trim()
            .parse::<usize>()
            .map_err(|_| format!("Enter a whole number for {}", self.fields[index].0))
    }
    fn reaction(&self, cx: &App) -> Result<Reaction, String> {
        let r = Reaction {
            dna_nm: self.number(0, cx)?,
            monovalent_mm: self.number(1, cx)?,
            magnesium_mm: self.number(2, cx)?,
            dntp_mm: self.number(3, cx)?,
            dmso_percent: self.number(4, cx)?,
            formamide_m: self.number(5, cx)?,
            temperature_c: self.number(6, cx)?,
            melting_model: self.melting_model,
            salt_correction: self.salt,
        };
        r.validate()?;
        Ok(r)
    }
    fn options(&self, cx: &App) -> Result<DesignOptions, String> {
        let from = self.fields[19].1.read(cx).value();
        let to = self.fields[20].1.read(cx).value();
        let target = if from.trim().is_empty() && to.trim().is_empty() {
            None
        } else {
            let start = self
                .integer(19, cx)?
                .checked_sub(1)
                .ok_or("Target starts at 1 or greater")?;
            let mut end = self.integer(20, cx)?;
            if self.document.circular && (end == 0 || end > self.document.sequence.len()) {
                return Err("Circular target end must be within the molecule (1-based)".into());
            }
            if self.document.circular && end <= start {
                end = end
                    .checked_add(self.document.sequence.len())
                    .ok_or("Target coordinate overflow")?;
            }
            Some(start..end)
        };
        Ok(DesignOptions {
            reaction: self.reaction(cx)?,
            min_length: self.integer(7, cx)?,
            optimal_length: self.integer(8, cx)?,
            max_length: self.integer(9, cx)?,
            min_tm: self.number(10, cx)?,
            optimal_tm: self.number(11, cx)?,
            max_tm: self.number(12, cx)?,
            min_gc: self.number(13, cx)?,
            max_gc: self.number(14, cx)?,
            max_tm_difference: self.number(15, cx)?,
            product_size: self.integer(16, cx)?
                ..self
                    .integer(17, cx)?
                    .checked_add(1)
                    .ok_or("Product length is too large")?,
            num_return: self.integer(18, cx)?,
            target,
        })
    }
    fn analyze(&mut self, window: &mut Window, cx: &mut Context<Self>) {
        if self.busy {
            return;
        }
        let reaction = match self.reaction(cx) {
            Ok(r) => r,
            Err(e) => {
                self.error(e, cx);
                return;
            }
        };
        let sequence = match Molecule::Dna.normalize(&self.sequence.read(cx).value()) {
            Ok(s) => s,
            Err(e) => {
                self.error(e.to_string(), cx);
                return;
            }
        };
        self.busy = true;
        self.error.clear();
        cx.notify();
        cx.spawn_in(window, async move |this, cx| {
            let result = cx
                .background_executor()
                .spawn(async move { primer::analyze(&sequence, reaction) })
                .await;
            let _ = this.update_in(cx, |this, _, cx| {
                this.busy = false;
                match result {
                    Ok(result) => {
                        this.analysis = Some(result);
                        this.result_reaction = Some(reaction);
                    }
                    Err(e) => this.error = e,
                }
                cx.notify();
            });
        })
        .detach();
    }
    fn design(&mut self, window: &mut Window, cx: &mut Context<Self>) {
        if self.busy {
            return;
        }
        let options = match self.options(cx) {
            Ok(o) => o,
            Err(e) => {
                self.error(e, cx);
                return;
            }
        };
        let reaction = options.reaction;
        let sequence = self.document.sequence.clone();
        let circular = self.document.circular;
        self.busy = true;
        self.error.clear();
        cx.notify();
        cx.spawn_in(window, async move |this, cx| {
            let result = cx
                .background_executor()
                .spawn(async move {
                    if circular {
                        primer::design_pairs_circular(&sequence, &options)
                    } else {
                        primer::design_pairs(&sequence, &options)
                    }
                })
                .await;
            let _ = this.update_in(cx, |this, _, cx| {
                this.busy = false;
                match result {
                    Ok(result) => {
                        this.design = Some(result);
                        this.result_reaction = Some(reaction);
                    }
                    Err(e) => this.error = e,
                }
                cx.notify();
            });
        })
        .detach();
    }
    fn save(&mut self, cx: &mut Context<Self>) {
        let result: Result<_, String> = (|| {
            let name = self.name.read(cx).value().trim().to_string();
            if name.is_empty() {
                return Err("Enter a primer name".into());
            }
            let sequence = Molecule::Dna
                .normalize(&self.sequence.read(cx).value())
                .map_err(|e| e.to_string())?;
            if sequence.is_empty() {
                return Err("Enter an oligo sequence".into());
            }
            let bindings = parse_bindings(
                &self.bindings.read(cx).value(),
                self.document.sequence.len(),
            )?;
            Ok((name, sequence, bindings, self.reaction(cx)?))
        })();
        match result {
            Ok((name, sequence, bindings, reaction)) => cx.emit(PrimerEvent::Save {
                feature_id: self.editing,
                name,
                sequence,
                bindings,
                reverse: self.reverse,
                reaction,
                tm_c: self.analysis.as_ref().map(|a| a.tm_c),
            }),
            Err(error) => self.error(error, cx),
        }
    }
}

fn parse_bindings(text: &str, length: usize) -> Result<Vec<Segment>, String> {
    if text.trim().is_empty() {
        return Ok(vec![]);
    }
    text.split(',')
        .map(|part| {
            let (start, end) = part
                .trim()
                .split_once("..")
                .ok_or("Use inclusive ranges such as 10..30, 50..60; leave empty for unbound")?;
            let start = start
                .trim()
                .parse::<usize>()
                .map_err(|_| "Invalid binding start")?
                .checked_sub(1)
                .ok_or("Binding coordinates start at 1")?;
            let end = end
                .trim()
                .parse::<usize>()
                .map_err(|_| "Invalid binding end")?;
            if start >= end || end > length {
                return Err("Binding interval is outside the document".into());
            }
            Ok(Segment { start, end })
        })
        .collect()
}
fn field(label: &'static str, input: impl IntoElement) -> Div {
    div()
        .v_flex()
        .gap_1()
        .child(div().text_xs().child(label))
        .child(input)
}
impl Render for PrimerTools {
    fn render(&mut self, _: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        let mut body = div().id("primer-tools-scroll").v_flex().gap_3().p_5().max_h(px(760.)).overflow_y_scroll()
            .child(div().text_xl().child("Primers"))
            .child(div().h_flex().gap_2().children([(Mode::Oligo, "Analyze / edit oligo"), (Mode::Design, "Design primer pairs")].into_iter().enumerate().map(|(i,(mode,label))| {
                Button::new(("primer-mode", i)).small().selected(self.mode == mode).label(label).disabled(self.busy)
                    .on_click(cx.listener(move |this, _, _, cx| { this.mode = mode; cx.notify(); }))
            })))
            .child(div().text_xs().text_color(rgb(0x96a6bc)).child("Primer3 · full-oligo perfect-match Tm for 2–60 unmodified DNA bases. Tails and binding mismatches are not annealed-region models. Tm is not a recommended PCR annealing temperature."));
        if self.mode == Mode::Oligo {
            body = body
                .child(field(
                    "Primer name",
                    Input::new(&self.name).disabled(self.busy),
                ))
                .child(field(
                    "Oligo 5′ → 3′",
                    Textarea::new(&self.sequence)
                        .disabled(self.busy)
                        .h(px(90.))
                        .font_family("monospace"),
                ))
                .child(field(
                    "Binding intervals (1-based inclusive; empty = unbound)",
                    Input::new(&self.bindings).disabled(self.busy),
                ))
                .child(
                    div().h_flex().gap_2().children(
                        [(false, "Forward"), (true, "Reverse")]
                            .into_iter()
                            .enumerate()
                            .map(|(i, (reverse, label))| {
                                Button::new(("primer-strand", i))
                                    .small()
                                    .selected(self.reverse == reverse)
                                    .disabled(self.busy)
                                    .label(label)
                                    .on_click(cx.listener(move |this, _, _, cx| {
                                        this.reverse = reverse;
                                        cx.notify();
                                    }))
                            }),
                    ),
                );
        }
        body = body
            .child(div().text_sm().child("Reaction parameters"))
            .child(
                div()
                    .h_flex()
                    .flex_wrap()
                    .gap_3()
                    .children(self.fields.iter().take(7).map(|(label, input)| {
                        field(label, Input::new(input).w(px(155.)).disabled(self.busy))
                    })),
            )
            .child(
                div().h_flex().gap_2().children(
                    [
                        (MeltingModel::SantaLucia, "SantaLucia"),
                        (MeltingModel::Breslauer, "Breslauer"),
                    ]
                    .into_iter()
                    .enumerate()
                    .map(|(i, (model, label))| {
                        Button::new(("primer-model", i))
                            .small()
                            .selected(self.melting_model == model)
                            .label(label)
                            .disabled(self.busy)
                            .on_click(cx.listener(move |this, _, _, cx| {
                                this.melting_model = model;
                                this.invalidate(cx);
                            }))
                    }),
                ),
            )
            .child(
                div().h_flex().gap_2().children(
                    [
                        (SaltCorrection::SantaLucia, "SantaLucia salt"),
                        (SaltCorrection::Owczarzy, "Owczarzy salt"),
                        (SaltCorrection::Schildkraut, "Schildkraut salt"),
                    ]
                    .into_iter()
                    .enumerate()
                    .map(|(i, (salt, label))| {
                        Button::new(("primer-salt", i))
                            .small()
                            .selected(self.salt == salt)
                            .label(label)
                            .disabled(self.busy)
                            .on_click(cx.listener(move |this, _, _, cx| {
                                this.salt = salt;
                                this.invalidate(cx);
                            }))
                    }),
                ),
            );
        if self.mode == Mode::Design {
            body =
                body.child(div().text_sm().child("Design constraints"))
                    .child(div().h_flex().flex_wrap().gap_3().children(
                        self.fields.iter().skip(7).map(|(label, input)| {
                            field(label, Input::new(input).w(px(155.)).disabled(self.busy))
                        }),
                    ))
                    .child(
                        Button::new("run-primer-design")
                            .primary()
                            .disabled(self.busy)
                            .label(if self.busy {
                                "Designing…"
                            } else {
                                "Design pairs"
                            })
                            .on_click(cx.listener(|this, _, window, cx| this.design(window, cx))),
                    );
        } else {
            body = body.child(
                div()
                    .h_flex()
                    .gap_2()
                    .child(
                        Button::new("analyze-primer")
                            .primary()
                            .disabled(self.busy)
                            .label("Calculate Tm / structures")
                            .on_click(cx.listener(|this, _, window, cx| this.analyze(window, cx))),
                    )
                    .child(
                        Button::new("save-primer")
                            .disabled(self.busy)
                            .label(if self.editing.is_some() {
                                "Apply primer edits"
                            } else {
                                "Add primer"
                            })
                            .on_click(cx.listener(|this, _, _, cx| this.save(cx))),
                    ),
            );
        }
        if let Some(a) = &self.analysis {
            body = body.child(
                div()
                    .p_3()
                    .v_flex()
                    .gap_1()
                    .bg(rgb(0x203445))
                    .child(format!(
                        "Full-oligo Tm {:.2} °C · GC {:.1}%",
                        a.tm_c, a.gc_percent
                    ))
                    .child(format!(
                        "Hairpin: {} · {:.2} °C · ΔG {:.2} kcal/mol",
                        if a.hairpin.found {
                            "predicted"
                        } else {
                            "none predicted"
                        },
                        a.hairpin.tm_c,
                        a.hairpin.free_energy_cal_mol / 1000.
                    ))
                    .child(format!(
                        "Homodimer: {} · {:.2} °C · ΔG {:.2} kcal/mol",
                        if a.homodimer.found {
                            "predicted"
                        } else {
                            "none predicted"
                        },
                        a.homodimer.tm_c,
                        a.homodimer.free_energy_cal_mol / 1000.
                    ))
                    .children(a.warnings.iter().map(|w| div().text_xs().child(w.clone()))),
            );
        }
        if let Some(result) = &self.design {
            body = body.children(
                result
                    .warnings
                    .iter()
                    .map(|w| div().text_xs().child(w.clone())),
            );
            for (i, pair) in result.pairs.iter().enumerate() {
                let pair = pair.clone();
                let reaction = self.result_reaction.unwrap_or_default();
                body = body.child(
                    div()
                        .v_flex()
                        .gap_1()
                        .p_3()
                        .border_1()
                        .border_color(rgb(0x344558))
                        .child(format!(
                            "Pair {} · product {}–{} ({} bp) · penalty {:.3}",
                            i + 1,
                            pair.product.start + 1,
                            if self.document.circular {
                                (pair.product.end - 1) % self.document.sequence.len() + 1
                            } else {
                                pair.product.end
                            },
                            pair.product.len(),
                            pair.penalty
                        ))
                        .child(div().font_family("monospace").child(format!(
                            "F {} · {:.2} °C",
                            pair.forward.sequence, pair.forward.tm_c
                        )))
                        .child(div().font_family("monospace").child(format!(
                            "R {} · {:.2} °C",
                            pair.reverse.sequence, pair.reverse.tm_c
                        )))
                        .child(
                            Button::new(("add-primer-pair", i))
                                .small()
                                .label("Add this pair")
                                .on_click(cx.listener(move |_, _, _, cx| {
                                    cx.emit(PrimerEvent::AddPair {
                                        pair: pair.clone(),
                                        reaction,
                                    })
                                })),
                        ),
                );
            }
        }
        if !self.error.is_empty() {
            body = body.child(div().text_color(rgb(0xf2a4a4)).child(self.error.clone()));
        }
        body.child(
            Button::new("close-primer-tools")
                .ghost()
                .label("Close")
                .on_click(cx.listener(|_, _, _, cx| cx.emit(PrimerEvent::Cancel))),
        )
    }
}
