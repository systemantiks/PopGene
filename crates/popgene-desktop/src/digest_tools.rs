//! Restriction digestion setup, scientific assumptions and product selection.
use gpui_kit::component::{
    Disableable, Selectable, Sizable, StyledExt,
    button::{Button, ButtonVariants},
    input::{Input, InputEvent, InputState},
};
use gpui_kit::*;
use popgene_core::{
    Document, Molecule,
    restriction::{
        self, CleavageAssessment, DigestConditions, DigestFragment, DigestResult, ENZYMES, EndKind,
        FragmentEnd, Methylation,
    },
};
use std::{collections::BTreeSet, sync::Arc};

pub enum DigestEvent {
    Create {
        fragments: Vec<DigestFragment>,
        conditions: DigestConditions,
        enzymes: Vec<String>,
    },
    Cancel,
}
pub struct DigestTools {
    document: Document,
    enzymes: Entity<InputState>,
    search: Entity<InputState>,
    fields: Vec<(&'static str, Entity<InputState>)>,
    methylation: Methylation,
    result: Option<Arc<DigestResult>>,
    result_conditions: DigestConditions,
    result_enzymes: Vec<String>,
    selected: BTreeSet<usize>,
    busy: bool,
    error: String,
    _subscriptions: Vec<Subscription>,
}
impl EventEmitter<DigestEvent> for DigestTools {}
impl DigestTools {
    pub fn new(
        document: Document,
        selected_enzyme: Option<usize>,
        window: &mut Window,
        cx: &mut Context<Self>,
    ) -> Self {
        let enzymes = cx.new(|cx| {
            InputState::new(window, cx).default_value(
                selected_enzyme
                    .and_then(|i| ENZYMES.get(i))
                    .map_or("EcoRI", |e| e.name),
            )
        });
        let search = cx.new(|cx| InputState::new(window, cx).placeholder("Find an enzyme"));
        let fields = [
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
        let mut subscriptions = vec![];
        for input in std::iter::once(&enzymes).chain(fields.iter().map(|(_, i)| i)) {
            subscriptions.push(cx.subscribe(input, |this, _, event: &InputEvent, cx| {
                if matches!(event, InputEvent::Change) {
                    this.invalidate(cx);
                }
            }));
        }
        subscriptions.push(cx.subscribe(&search, |_, _, _: &InputEvent, cx| cx.notify()));
        let methylation = initial_methylation(&document);
        Self {
            document,
            enzymes,
            search,
            fields,
            methylation,
            result: None,
            result_conditions: DigestConditions::default(),
            result_enzymes: vec![],
            selected: BTreeSet::new(),
            busy: false,
            error: String::new(),
            _subscriptions: subscriptions,
        }
    }
    pub fn error(&mut self, error: impl Into<String>, cx: &mut Context<Self>) {
        self.error = error.into();
        cx.notify();
    }
    fn invalidate(&mut self, cx: &mut Context<Self>) {
        self.result = None;
        self.selected.clear();
        self.error.clear();
        cx.notify();
    }
    fn run(&mut self, window: &mut Window, cx: &mut Context<Self>) {
        if self.busy {
            return;
        }
        let parameters: Result<_, String> = (|| {
            if self.document.molecule != Molecule::Dna
                || self.document.single_stranded
                || self.document.trace.is_some()
            {
                return Err(
                    "Restriction digestion requires double-stranded DNA without trace data".into(),
                );
            }
            let indices = parse_enzymes(&self.enzymes.read(cx).value())?;
            let mut values = vec![];
            for (label, input) in &self.fields {
                values.push(
                    input
                        .read(cx)
                        .value()
                        .trim()
                        .parse::<f64>()
                        .map_err(|_| format!("Enter a number for {label}"))?,
                );
            }
            Ok((
                indices,
                DigestConditions {
                    temperature_c: values[0],
                    ph: values[1],
                    magnesium_mm: values[2],
                    glycerol_percent: values[3],
                    methylation: self.methylation,
                },
            ))
        })();
        let (indices, conditions) = match parameters {
            Ok(p) => p,
            Err(e) => {
                self.error(e, cx);
                return;
            }
        };
        let names = indices
            .iter()
            .map(|&i| ENZYMES[i].name.to_string())
            .collect::<Vec<_>>();
        let document = self.document.clone();
        self.busy = true;
        self.invalidate(cx);
        cx.spawn_in(window, async move |this, cx| {
            let result = cx
                .background_executor()
                .spawn(async move { restriction::digest_document(&document, &indices, conditions) })
                .await;
            let _ = this.update_in(cx, |this, _, cx| {
                this.busy = false;
                match result {
                    Ok(result) => {
                        if result.fragments.len() <= 32 && result.all_sites_resolved {
                            this.selected.extend(0..result.fragments.len());
                        }
                        this.result = Some(Arc::new(result));
                        this.result_conditions = conditions;
                        this.result_enzymes = names;
                    }
                    Err(e) => this.error = e,
                }
                cx.notify();
            });
        })
        .detach();
    }
    fn create(&mut self, cx: &mut Context<Self>) {
        let Some(result) = &self.result else { return };
        if !result.all_sites_resolved || self.selected.is_empty() {
            return;
        }
        cx.emit(DigestEvent::Create {
            fragments: self
                .selected
                .iter()
                .filter_map(|&i| result.fragments.get(i).cloned())
                .collect(),
            conditions: self.result_conditions,
            enzymes: self.result_enzymes.clone(),
        });
    }
}
fn parse_enzymes(text: &str) -> Result<Vec<usize>, String> {
    let mut indices = vec![];
    for name in text
        .split(|c: char| c == ',' || c == ';' || c.is_whitespace())
        .filter(|s| !s.is_empty())
    {
        let i = ENZYMES
            .iter()
            .position(|e| e.name.eq_ignore_ascii_case(name))
            .ok_or_else(|| format!("Unknown enzyme: {name}"))?;
        if !indices.contains(&i) {
            indices.push(i);
        }
    }
    if indices.is_empty() || indices.len() > 64 {
        return Err("Choose between 1 and 64 enzymes".into());
    }
    Ok(indices)
}
fn initial_methylation(document: &Document) -> Methylation {
    let recorded = document.methylation.unwrap_or_default();
    Methylation {
        dam: recorded.dam,
        dcm: recorded.dcm,
        cpg: false,
    }
}
fn methylation_source(document: &Document) -> String {
    match document.methylation {
        Some(recorded) => format!(
            "Initialized from recorded flags: Dam {}, Dcm {}, EcoKI {}. Flags are assumptions, not measured occupancy. CpG starts off; EcoKI sensitivity is not modeled.",
            recorded.dam, recorded.dcm, recorded.eco_ki
        ),
        None => "Template methylation is unknown. All assumptions start off; choose the conditions to simulate. EcoKI sensitivity is not modeled.".into(),
    }
}
fn end_label(end: Option<&FragmentEnd>) -> String {
    match end {
        None => "closed".into(),
        Some(e) => match e.kind {
            EndKind::Native => "original end".into(),
            EndKind::Blunt => "blunt".into(),
            EndKind::FivePrime => format!("5′ {}", e.overhang),
            EndKind::ThreePrime => format!("3′ {}", e.overhang),
        },
    }
}
fn report(
    result: &DigestResult,
    names: &[String],
    conditions: DigestConditions,
    document: &Document,
) -> String {
    let mut text = format!(
        "Restriction digestion: {}\n{} °C; pH {}; Mg {} mM; glycerol {}%; fully methylated Dam={}, Dcm={}, CpG={}\n",
        names.join(", "),
        conditions.temperature_c,
        conditions.ph,
        conditions.magnesium_mm,
        conditions.glycerol_percent,
        conditions.methylation.dam,
        conditions.methylation.dcm,
        conditions.methylation.cpg
    );
    text.push_str(&format!("{}\n", methylation_source(document)));
    for warning in &result.warnings {
        text.push_str(&format!("Assumption: {warning}\n"));
    }
    for assessed in &result.sites {
        let assessment = match &assessed.assessment {
            CleavageAssessment::Permitted => "permitted".into(),
            CleavageAssessment::Blocked(s) => format!("blocked: {s}"),
            CleavageAssessment::Unresolved(s) => format!("unresolved: {s}"),
        };
        text.push_str(&format!(
            "{} at {}: {assessment}\n",
            ENZYMES[assessed.site.enzyme_index].name,
            assessed.site.start + 1
        ));
    }
    text.push_str("\nProduct\tTop strand nt\tOpposite strand nt\tOriginal top coordinates\tLeft end\tRight end\n");
    for (i, f) in result.fragments.iter().enumerate() {
        text.push_str(&format!(
            "{}\t{}\t{}\t{}..{}{}\t{}\t{}\n",
            i + 1,
            f.sequence.len(),
            f.opposite_sequence.len(),
            f.start + 1,
            f.end,
            if f.wraps_origin {
                " (across origin)"
            } else {
                ""
            },
            end_label(f.left_end.as_ref()),
            end_label(f.right_end.as_ref())
        ));
    }
    text
}
impl Render for DigestTools {
    fn render(&mut self, _: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        let query = self.search.read(cx).value().to_ascii_lowercase();
        let mut body = div()
            .id("digest-tools-scroll")
            .v_flex()
            .gap_3()
            .p_5()
            .max_h(px(780.))
            .overflow_y_scroll()
            .child(div().text_xl().child("Restriction digestion"))
            .child(div().text_sm().child(format!(
                "{} · {} bp · {}",
                self.document.name,
                self.document.sequence.len(),
                if self.document.circular {
                    "circular"
                } else {
                    "linear"
                }
            )))
            .child(div().text_xs().child("Enzymes (comma-separated)"))
            .child(Input::new(&self.enzymes).disabled(self.busy))
            .child(Input::new(&self.search).disabled(self.busy));
        if !query.is_empty() {
            body = body.child(
                div().h_flex().flex_wrap().gap_1().children(
                    ENZYMES
                        .iter()
                        .enumerate()
                        .filter(|(_, e)| e.name.to_ascii_lowercase().contains(&query))
                        .take(12)
                        .map(|(i, e)| {
                            Button::new(("digest-add-enzyme", i))
                                .small()
                                .label(e.name)
                                .disabled(self.busy)
                                .on_click(cx.listener(move |this, _, window, cx| {
                                    let mut value = this.enzymes.read(cx).value().to_string();
                                    if !value.trim().is_empty() {
                                        value.push_str(", ");
                                    }
                                    value.push_str(ENZYMES[i].name);
                                    this.enzymes
                                        .update(cx, |input, cx| input.set_value(value, window, cx));
                                    this.invalidate(cx);
                                }))
                        }),
                ),
            );
        }
        body = body
            .child(
                div()
                    .h_flex()
                    .flex_wrap()
                    .gap_3()
                    .children(self.fields.iter().map(|(label, input)| {
                        div()
                            .v_flex()
                            .gap_1()
                            .child(div().text_xs().child(*label))
                            .child(Input::new(input).w(px(155.)).disabled(self.busy))
                    })),
            )
            .child(
                div()
                    .text_xs()
                    .child("Methylation assumptions (complete on both strands)"),
            )
            .child(div().text_xs().child(methylation_source(&self.document)))
            .child(
                div().h_flex().gap_2().children(
                    [
                        (0usize, "Dam", self.methylation.dam),
                        (1, "Dcm", self.methylation.dcm),
                        (2, "CpG", self.methylation.cpg),
                    ]
                    .into_iter()
                    .map(|(i, label, selected)| {
                        Button::new(("digest-methylation", i))
                            .small()
                            .selected(selected)
                            .label(label)
                            .disabled(self.busy)
                            .on_click(cx.listener(move |this, _, _, cx| {
                                match i {
                                    0 => this.methylation.dam = !this.methylation.dam,
                                    1 => this.methylation.dcm = !this.methylation.dcm,
                                    _ => this.methylation.cpg = !this.methylation.cpg,
                                };
                                this.invalidate(cx);
                            }))
                    }),
                ),
            )
            .child(
                Button::new("run-digest")
                    .primary()
                    .disabled(self.busy)
                    .label(if self.busy {
                        "Calculating…"
                    } else {
                        "Simulate digestion"
                    })
                    .on_click(cx.listener(|this, _, window, cx| this.run(window, cx))),
            );
        if let Some(result) = &self.result {
            let result = result.clone();
            let unresolved = result
                .sites
                .iter()
                .filter(|s| matches!(s.assessment, CleavageAssessment::Unresolved(_)))
                .count();
            let blocked = result
                .sites
                .iter()
                .filter(|s| matches!(s.assessment, CleavageAssessment::Blocked(_)))
                .count();
            body = body
                .child(div().text_sm().child(format!(
                    "{} products · {} blocked sites · {} unresolved sites",
                    result.fragments.len(),
                    blocked,
                    unresolved
                )))
                .children(
                    result
                        .warnings
                        .iter()
                        .map(|w| div().text_xs().text_color(rgb(0xe5b879)).child(w.clone())),
                );
            if !result.all_sites_resolved {
                body=body.child(div().text_sm().text_color(rgb(0xe5b879)).child("These products are provisional. Resolve the listed site assumptions before creating documents."));
            }
            let owner = cx.entity().downgrade();
            let selected = self.selected.clone();
            let list_result = result.clone();
            body = body
                .child(
                    uniform_list(
                        "digest-products",
                        result.fragments.len(),
                        move |range, _, _| {
                            range
                                .map(|i| {
                                    let f = &list_result.fragments[i];
                                    let owner = owner.clone();
                                    div().h(px(54.)).px_2().child(
                                        Button::new(("digest-product", i))
                                            .small()
                                            .w_full()
                                            .selected(selected.contains(&i))
                                            .disabled(!list_result.all_sites_resolved)
                                            .label(format!(
                                                "{} · {} / {} nt · {}..{}{} · {} | {}",
                                                i + 1,
                                                f.sequence.len(),
                                                f.opposite_sequence.len(),
                                                f.start + 1,
                                                f.end,
                                                if f.wraps_origin { " ↻" } else { "" },
                                                end_label(f.left_end.as_ref()),
                                                end_label(f.right_end.as_ref())
                                            ))
                                            .on_click(move |_, _, cx| {
                                                let _ = owner.update(cx, |this, cx| {
                                                    if !this.selected.remove(&i) {
                                                        if this.selected.len() >= 64 {
                                                            this.error(
                                                                "Select at most 64 products",
                                                                cx,
                                                            );
                                                            return;
                                                        }
                                                        this.selected.insert(i);
                                                    }
                                                    cx.notify();
                                                });
                                            }),
                                    )
                                })
                                .collect::<Vec<_>>()
                        },
                    )
                    .h(px(240.)),
                )
                .child(
                    div()
                        .h_flex()
                        .gap_2()
                        .child(
                            Button::new("create-digest-products")
                                .primary()
                                .disabled(!result.all_sites_resolved || self.selected.is_empty())
                                .label(format!("Create {} selected products", self.selected.len()))
                                .on_click(cx.listener(|this, _, _, cx| this.create(cx))),
                        )
                        .child(
                            Button::new("copy-digest-report")
                                .label("Copy report")
                                .on_click(cx.listener(|this, _, _, cx| {
                                    if let Some(result) = &this.result {
                                        cx.write_to_clipboard(ClipboardItem::new_string(report(
                                            result,
                                            &this.result_enzymes,
                                            this.result_conditions,
                                            &this.document,
                                        )));
                                    }
                                })),
                        ),
                );
            body = body.children(
                result
                    .sites
                    .iter()
                    .filter_map(|s| match &s.assessment {
                        CleavageAssessment::Permitted => None,
                        CleavageAssessment::Blocked(reason)
                        | CleavageAssessment::Unresolved(reason) => {
                            Some(div().text_xs().child(format!(
                                "{} at {}: {}",
                                ENZYMES[s.site.enzyme_index].name,
                                s.site.start + 1,
                                reason
                            )))
                        }
                    })
                    .take(100),
            );
        }
        if !self.error.is_empty() {
            body = body.child(
                div()
                    .text_sm()
                    .text_color(rgb(0xef938a))
                    .child(self.error.clone()),
            );
        }
        body.child(
            Button::new("close-digest-tools")
                .label("Close")
                .on_click(cx.listener(|_, _, _, cx| cx.emit(DigestEvent::Cancel))),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::{initial_methylation, methylation_source, parse_enzymes};
    use popgene_core::restriction::ENZYMES;
    use popgene_core::{Document, metadata::DnaMethylation};
    #[test]
    fn exact_enzyme_names_and_duplicate_selection() {
        let values = parse_enzymes("ecori, EcoRI; EcoRV").unwrap();
        assert_eq!(values.len(), 2);
        assert_eq!(ENZYMES[values[0]].name, "EcoRI");
        assert!(parse_enzymes("EcoR").is_err());
        assert!(parse_enzymes("").is_err());
    }
    #[test]
    fn imported_flags_initialize_assumptions_without_inventing_cpg_or_occupancy() {
        let mut document = Document::new("template", "GATCCAGCTG").unwrap();
        assert!(!initial_methylation(&document).dam);
        assert!(methylation_source(&document).contains("unknown"));
        document.methylation = Some(DnaMethylation {
            dam: true,
            dcm: false,
            eco_ki: true,
        });
        let assumptions = initial_methylation(&document);
        assert!(assumptions.dam);
        assert!(!assumptions.dcm);
        assert!(!assumptions.cpg);
        assert!(methylation_source(&document).contains("EcoKI true"));
        assert!(methylation_source(&document).contains("not measured occupancy"));
        document.replace(0, 1, "A").unwrap();
        assert!(methylation_source(&document).contains("unknown"));
        assert!(!initial_methylation(&document).dam);
    }
}
