//! Read-only protein composition and explicitly scoped sequence calculations.
use crate::feature_card::number;
use gpui_kit::component::{
    Disableable, Selectable, Sizable, StyledExt,
    button::*,
    input::{Input, InputEvent, InputState},
};
use gpui_kit::{prelude::*, *};
use popgene_core::{
    Document, Molecule, Segment,
    protein::{
        self, DetailedProteinProperties, Modification, PkaModel, ProteinOptions, ProteinProperties,
    },
};
use std::fmt::Write;

fn stat(label: &'static str, value: String) -> Div {
    div()
        .v_flex()
        .gap_1()
        .min_w(px(160.))
        .child(div().text_xs().text_color(rgb(0x96a6bc)).child(label))
        .child(div().text_lg().font_weight(FontWeight::MEDIUM).child(value))
}

fn report(properties: &ProteinProperties) -> String {
    let mut out = format!(
        "Protein properties\n{} sequence symbols; {} non-stop residues; {} stop symbols\n",
        properties.symbol_count, properties.residue_count, properties.stop_count
    );
    if let Some(mass) = properties.mass {
        let _ = writeln!(
            out,
            "Average neutral mass: {:.4} Da\nMonoisotopic neutral mass: {:.6} Da",
            mass.average_da, mass.monoisotopic_da
        );
    }
    if let Some(pi) = properties.theoretical_pi {
        let _ = writeln!(out, "Theoretical pI: {pi:.2} (Bjellqvist pKa set)");
    }
    if let Some(charge) = properties.charge_at_ph7 {
        let _ = writeln!(out, "Net charge at pH 7: {charge:.2} (Bjellqvist pKa set)");
    }
    if let Some(absorbance) = properties.absorbance {
        let _ = writeln!(
            out,
            "Extinction at 280 nm (reduced Cys): {} M^-1 cm^-1\nExtinction at 280 nm (all possible Cys pairs form cystines): {} M^-1 cm^-1\nA280 at 1 mg/mL, 1 cm (reduced Cys): {:.2}",
            absorbance.reduced_extinction,
            absorbance.paired_cysteine_extinction,
            absorbance.reduced_absorbance
        );
    }
    if let Some(reason) = &properties.unavailable_reason {
        let _ = writeln!(out, "{reason}");
    }
    out.push_str("\nComposition (percent of all stored sequence symbols, including stops)\nSymbol\tResidue\tCount\tPercent\n");
    for residue in &properties.composition {
        let _ = writeln!(
            out,
            "{}\t{}\t{}\t{:.2}%",
            residue.symbol, residue.name, residue.count, residue.percent
        );
    }
    out.push_str("\nModel: one linear unmodified chain; free N/C termini, reduced cysteines. Mass uses average/monoisotopic amino-acid masses minus one water per peptide bond. pI/charge uses Bjellqvist pKa values with residue-specific terminal corrections and charge balance over pH 0–14. Extinction uses Pace/Edelhoch Trp 5500, Tyr 1490, cystine 125; paired Cys is a separate hypothetical optical scenario. No disulfide bonds, modifications or processing are inferred for mass or pI.\nMethods: https://web.expasy.org/compute_pi/pi_tool-doc.html\nExtinction: https://web.expasy.org/protparam/protparam-doc.html\nMass tables: https://github.com/biopython/biopython/blob/biopython-188/Bio/Data/IUPACData.py\n");
    out
}

/// Pass a cached analysis for the current protein document revision.
pub(crate) fn render(properties: &ProteinProperties) -> AnyElement {
    let copy = report(properties);
    let mut summary = div()
        .h_flex()
        .flex_wrap()
        .gap_6()
        .child(stat("SEQUENCE SYMBOLS", number(properties.symbol_count)))
        .child(stat("RESIDUES", number(properties.residue_count)))
        .child(stat("STOP SYMBOLS", number(properties.stop_count)));
    if let Some(mass) = properties.mass {
        summary = summary
            .child(stat("AVERAGE MASS", format!("{:.2} Da", mass.average_da)))
            .child(stat(
                "MONOISOTOPIC MASS",
                format!("{:.4} Da", mass.monoisotopic_da),
            ));
    } else {
        summary = summary.child(stat("MOLECULAR MASS", "Unavailable".into()));
    }
    summary = summary.child(stat(
        "THEORETICAL pI · BJELLQVIST",
        properties
            .theoretical_pi
            .map_or_else(|| "Unavailable".into(), |pi| format!("{pi:.2}")),
    ));
    summary = summary.child(stat(
        "CHARGE AT pH 7 · BJELLQVIST",
        properties
            .charge_at_ph7
            .map_or_else(|| "Unavailable".into(), |charge| format!("{charge:+.2}")),
    ));
    if let Some(absorbance) = properties.absorbance {
        summary = summary
            .child(stat(
                "EXTINCTION 280 nm · REDUCED",
                format!("{} M⁻¹cm⁻¹", number(absorbance.reduced_extinction)),
            ))
            .child(stat(
                "EXTINCTION · CYS PAIRS",
                format!("{} M⁻¹cm⁻¹", number(absorbance.paired_cysteine_extinction)),
            ))
            .child(stat(
                "A280 · 1 mg/mL · 1 cm",
                format!("{:.2}", absorbance.reduced_absorbance),
            ));
    } else {
        summary = summary
            .child(stat("EXTINCTION 280 nm", "Unavailable".into()))
            .child(stat("A280 · 1 mg/mL · 1 cm", "Unavailable".into()));
    }
    div().id("protein-properties-scroll").flex_1().min_h_0().overflow_y_scroll().px_6().pb_6()
        .child(div().v_flex().gap_4().p_5().rounded_lg().border_1().border_color(rgb(0x293a4a)).bg(rgb(0x172231))
            .child(div().h_flex().gap_3()
                .child(div().text_lg().font_weight(FontWeight::SEMIBOLD).child("Protein properties"))
                .child(div().flex_1())
                .child(Button::new("copy-protein-properties").small().label("Copy properties")
                    .on_click(move |_, _, cx| cx.write_to_clipboard(ClipboardItem::new_string(copy.clone())))))
            .child(summary)
            .when_some(properties.unavailable_reason.as_ref(), |d, reason| d.child(div().text_sm().text_color(rgb(0xe7bd7b)).child(reason.clone())))
            .child(div().text_sm().text_color(rgb(0x96a6bc)).child("One linear, unmodified chain with free N/C termini and reduced cysteines. Mass and pI require the 20 standard amino acids; ambiguous residues and stop symbols remain in the composition below."))
            .child(div().text_sm().text_color(rgb(0x96a6bc)).child("Theoretical pI and charge use the Bjellqvist pKa set and residue-specific terminal corrections. These estimates can differ from other software and experiment. Modifications, disulfide bonds and protein processing are not inferred for mass or pI."))
            .child(div().text_sm().text_color(rgb(0x96a6bc)).child("Extinction uses Pace/Edelhoch coefficients. The Cys-pairs value assumes every possible cysteine pair forms a cystine; A280 uses reduced cysteines. No other chromophores are modeled."))
            .child(div().h_flex().flex_wrap().gap_2()
                .child(Button::new("protein-methods").small().ghost().label("Method reference")
                    .on_click(|_, _, cx| cx.open_url("https://web.expasy.org/compute_pi/pi_tool-doc.html")))
                .child(Button::new("protein-mass-tables").small().ghost().label("Mass table reference")
                    .on_click(|_, _, cx| cx.open_url("https://github.com/biopython/biopython/blob/biopython-188/Bio/Data/IUPACData.py")))
                .child(Button::new("protein-optics-reference").small().ghost().label("Extinction reference")
                    .on_click(|_, _, cx| cx.open_url("https://web.expasy.org/protparam/protparam-doc.html")))))
        .child(div().pt_5().pb_3().v_flex().gap_1()
            .child(div().text_lg().font_weight(FontWeight::MEDIUM).child("Amino acid composition"))
            .child(div().text_sm().text_color(rgb(0x96a6bc)).child("Percentages use the full sequence length, including stop symbols. Stop symbols are counted separately above.")))
        .child(div().h_flex().gap_4().px_4().py_2().text_xs().text_color(rgb(0x96a6bc))
            .child(div().w(px(54.)).child("CODE"))
            .child(div().flex_1().child("RESIDUE"))
            .child(div().w(px(100.)).text_right().child("COUNT"))
            .child(div().w(px(110.)).text_right().child("PERCENT")))
        .children(properties.composition.iter().map(|residue| div().h_flex().gap_4().px_4().py_2().border_b_1().border_color(rgb(0x293a4a))
            .child(div().w(px(54.)).font_family("monospace").child(residue.symbol.to_string()))
            .child(div().flex_1().child(residue.name))
            .child(div().w(px(100.)).text_right().child(number(residue.count)))
            .child(div().w(px(110.)).text_right().child(format!("{:.2}%", residue.percent)))))
        .into_any_element()
}

pub enum ProteinEvent {
    Cancel,
}
#[derive(Default)]
struct Calculation {
    generation: u64,
    pending: Option<u64>,
    result: Option<DetailedProteinProperties>,
    error: String,
}
impl Calculation {
    fn invalidate(&mut self) {
        self.generation = self.generation.wrapping_add(1);
        self.result = None;
        self.error.clear();
    }
    fn busy(&self) -> bool {
        self.pending.is_some()
    }
    fn begin(&mut self) -> Option<u64> {
        if self.busy() {
            return None;
        }
        self.invalidate();
        self.pending = Some(self.generation);
        self.pending
    }
    fn complete(&mut self, generation: u64, result: Result<DetailedProteinProperties, String>) {
        if self.pending != Some(generation) {
            return;
        }
        self.pending = None;
        if self.generation != generation {
            return;
        }
        match result {
            Ok(result) => self.result = Some(result),
            Err(error) => self.error = error,
        }
    }
    fn current(&self) -> Option<&DetailedProteinProperties> {
        (!self.busy()).then_some(self.result.as_ref()).flatten()
    }
    fn report(&self, name: &str) -> Option<String> {
        self.current().map(|result| detailed_report(name, result))
    }
    fn curve(&self) -> Option<String> {
        let p = self.current().filter(|p| !p.charge_curve.is_empty())?;
        Some(
            format!("pH\tNet charge ({})\n", p.model.label())
                + &p.charge_curve
                    .iter()
                    .map(|(ph, q)| format!("{ph:.1}\t{q:.8}\n"))
                    .collect::<String>(),
        )
    }
}
/// A captured document snapshot. Analysis is read-only and runs on a worker.
pub struct ProteinTools {
    document: Document,
    fields: Vec<(&'static str, Entity<InputState>)>,
    model: PkaModel,
    n_acetyl: bool,
    c_amide: bool,
    secondary: bool,
    calculation: Calculation,
    notice: String,
    _subscriptions: Vec<Subscription>,
}
impl EventEmitter<ProteinEvent> for ProteinTools {}
impl ProteinTools {
    pub fn new(
        document: Document,
        selected: Vec<Segment>,
        window: &mut Window,
        cx: &mut Context<Self>,
    ) -> Self {
        let range = if selected.len() == 1 {
            selected[0].start..selected[0].end
        } else {
            0..document.sequence.len()
        };
        let values = [
            (range.start + 1).to_string(),
            range.end.to_string(),
            "7".into(),
            String::new(),
            String::new(),
            String::new(),
            String::new(),
        ];
        let labels = [
            "First residue (1-based)",
            "Last residue (inclusive)",
            "Charge at pH",
            "Acetylated lysines (positions, comma-separated)",
            "Oxidized methionines (positions)",
            "Phosphorylated S/T/Y (positions)",
            "Disulfides (pairs, e.g. 2-8, 10-14)",
        ];
        let fields = labels
            .into_iter()
            .zip(values)
            .map(|(label, value)| {
                (
                    label,
                    cx.new(|cx| InputState::new(window, cx).default_value(value)),
                )
            })
            .collect::<Vec<_>>();
        let subscriptions = fields
            .iter()
            .map(|(_, input)| {
                cx.subscribe(input, |this: &mut Self, _, event: &InputEvent, cx| {
                    if matches!(event, InputEvent::Change) {
                        this.invalidate(cx);
                    }
                })
            })
            .collect();
        let notice = if selected.len() > 1 {
            "Multiple selected intervals were not joined into an artificial peptide. Choose one contiguous range below.".into()
        } else {
            String::new()
        };
        Self {
            document,
            fields,
            model: PkaModel::Bjellqvist,
            n_acetyl: false,
            c_amide: false,
            secondary: false,
            calculation: Calculation::default(),
            notice,
            _subscriptions: subscriptions,
        }
    }
    fn invalidate(&mut self, cx: &mut Context<Self>) {
        self.calculation.invalidate();
        cx.notify();
    }
    fn options(&self, cx: &App) -> Result<ProteinOptions, String> {
        let value = |i: usize| self.fields[i].1.read(cx).value().to_string();
        let start = value(0)
            .trim()
            .parse::<usize>()
            .map_err(|_| "Enter a valid first residue")?
            .checked_sub(1)
            .ok_or("First residue starts at 1")?;
        let end = value(1)
            .trim()
            .parse::<usize>()
            .map_err(|_| "Enter a valid last residue")?;
        let ph = value(2)
            .trim()
            .parse::<f64>()
            .map_err(|_| "Enter a numeric pH")?;
        Ok(ProteinOptions {
            range: Some(start..end),
            model: self.model,
            ph,
            secondary_structure: self.secondary,
            modifications: parse_modifications(
                self.n_acetyl,
                self.c_amide,
                &value(3),
                &value(4),
                &value(5),
                &value(6),
            )?,
        })
    }
    fn calculate(&mut self, window: &mut Window, cx: &mut Context<Self>) {
        let Some(generation) = self.calculation.begin() else {
            return;
        };
        if self.document.molecule != Molecule::Protein {
            self.calculation.complete(
                generation,
                Err("Protein analysis requires a protein document".into()),
            );
            cx.notify();
            return;
        }
        let options = match self.options(cx) {
            Ok(o) => o,
            Err(error) => {
                self.calculation.complete(generation, Err(error));
                cx.notify();
                return;
            }
        };
        let sequence = self.document.sequence.clone();
        cx.notify();
        cx.spawn_in(window, async move |this, cx| {
            let result = cx
                .background_executor()
                .spawn(async move { protein::analyze_with_options(&sequence, &options) })
                .await;
            let _ = this.update_in(cx, |this, _, cx| {
                this.calculation.complete(generation, result);
                cx.notify();
            });
        })
        .detach();
    }
}
fn positions(text: &str) -> Result<Vec<usize>, String> {
    if text.trim().is_empty() {
        return Ok(vec![]);
    }
    text.split(',')
        .map(|part| {
            part.trim()
                .parse::<usize>()
                .map_err(|_| {
                    "Modification positions must be comma-separated residue numbers".into()
                })
                .and_then(|p| {
                    p.checked_sub(1)
                        .ok_or("Modification positions start at 1".into())
                })
        })
        .collect()
}
fn parse_modifications(
    n: bool,
    c: bool,
    lys: &str,
    met: &str,
    phospho: &str,
    bonds: &str,
) -> Result<Vec<Modification>, String> {
    let mut result = vec![];
    if n {
        result.push(Modification::NAcetyl);
    }
    if c {
        result.push(Modification::CAmide);
    }
    result.extend(positions(lys)?.into_iter().map(Modification::LysAcetyl));
    result.extend(positions(met)?.into_iter().map(Modification::MetOxidation));
    result.extend(
        positions(phospho)?
            .into_iter()
            .map(Modification::Phosphorylation),
    );
    if !bonds.trim().is_empty() {
        for pair in bonds.split(',') {
            let (a, b) = pair
                .trim()
                .split_once('-')
                .ok_or("Disulfides require pairs such as 2-8, 10-14")?;
            let (a, b) = (positions(a)?, positions(b)?);
            if a.len() != 1 || b.len() != 1 {
                return Err("Enter two positions per disulfide".into());
            }
            result.push(Modification::Disulfide(a[0], b[0]));
        }
    }
    Ok(result)
}
fn detailed_report(name: &str, p: &DetailedProteinProperties) -> String {
    let mut out = format!(
        "Protein analysis: {name}\nResidues {}..{} (1-based inclusive), treated as one excised linear chain\nCharge model: {}\npH: {}\n",
        p.range.start + 1,
        p.range.end,
        p.model.label(),
        p.ph
    );
    let _ = writeln!(
        out,
        "{}\nHenderson–Hasselbalch independent groups; pI requires a charge-balance root within pH 0–14.",
        p.model.description()
    );
    for source in p.model.references() {
        let _ = writeln!(out, "{source}");
    }
    if let Some(mass) = p.mass {
        let _ = writeln!(
            out,
            "Average neutral mass: {:.4} Da\nMonoisotopic neutral mass: {:.6} Da",
            mass.average_da, mass.monoisotopic_da
        );
    }
    let _ = writeln!(
        out,
        "Theoretical pI: {}\nNet charge: {}",
        p.theoretical_pi
            .map_or("Unavailable".into(), |v| format!("{v:.6}")),
        p.charge.map_or("Unavailable".into(), |v| format!("{v:.6}"))
    );
    let _ = writeln!(
        out,
        "Extinction 280 nm: {} M^-1 cm^-1\nA280 at 1 mg/mL, 1 cm: {}",
        p.extinction_280
            .map_or("Unavailable".into(), |v| v.to_string()),
        p.absorbance_1mg_ml
            .map_or("Unavailable".into(), |v| format!("{v:.6}"))
    );
    out.push_str("Declared modifications:\n");
    if p.modifications.is_empty() {
        out.push_str("None\n");
    }
    for m in &p.modifications {
        let _ = writeln!(out, "{}", m.label());
    }
    for warning in &p.warnings {
        let _ = writeln!(out, "{warning}");
    }
    if let Some(reason) = &p.unmodified.unavailable_reason {
        let _ = writeln!(out, "{reason}");
    }
    out.push_str("\nComposition (percent of selected symbols, including stops)\n");
    for c in &p.unmodified.composition {
        let _ = writeln!(out, "{}\t{}\t{:.2}%", c.symbol, c.count, c.percent);
    }
    if let Some(s) = &p.secondary {
        let _ = writeln!(
            out,
            "\nSecondary prediction: {}\n{}\nH={} E={} C={}",
            s.method, s.assignments, s.helix_count, s.sheet_count, s.coil_count
        );
    }
    out.push_str("\nMass: Biopython 1.88 amino-acid tables, peptide bonds subtract water. Modification deltas: Unimod 1/2/21/35/2020. Optics: Pace/Edelhoch W5500, Y1490, explicitly declared cystine125.\nhttps://web.expasy.org/compute_pi/pi_tool-doc.html\nhttps://raw.githubusercontent.com/kimrutherford/EMBOSS/master/emboss/data/Epk.dat\nhttps://www.unimod.org/\nhttps://web.expasy.org/protscale/pscale/alpha-helixFasman.html\nhttps://web.expasy.org/protscale/pscale/beta-sheetFasman.html\n");
    out
}
impl Render for ProteinTools {
    fn render(&mut self, _: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        let mut body=div().id("protein-tools-scroll").v_flex().gap_3().p_5().max_h(px(760.)).overflow_y_scroll()
            .child(div().text_xl().child("Protein analysis"))
            .child(div().text_sm().child(self.document.name.clone()))
            .child(div().text_xs().text_color(rgb(0x96a6bc)).child("Analyze a contiguous interval as one peptide. All modification positions use the original document numbering. Terminal modifications apply to the chosen interval. No changes are made to the document."));
        if !self.notice.is_empty() {
            body = body.child(div().text_color(rgb(0xe7bd7b)).child(self.notice.clone()));
        }
        body =
            body.child(div().h_flex().flex_wrap().gap_3().children(
                self.fields.iter().take(3).map(|(label, input)| {
                    div()
                        .v_flex()
                        .gap_1()
                        .w(px(200.))
                        .child(div().text_xs().child(*label))
                        .child(Input::new(input).disabled(self.calculation.busy()))
                }),
            ))
            .child(
                div().h_flex().flex_wrap().gap_2().children(
                    [
                        PkaModel::Bjellqvist,
                        PkaModel::Emboss,
                        PkaModel::Grimsley2009,
                    ]
                    .into_iter()
                    .enumerate()
                    .map(|(i, model)| {
                        Button::new(("protein-pka-model", i))
                            .small()
                            .selected(self.model == model)
                            .disabled(self.calculation.busy())
                            .label(model.label())
                            .on_click(cx.listener(move |this, _, _, cx| {
                                if this.calculation.busy() {
                                    return;
                                }
                                this.model = model;
                                this.invalidate(cx);
                            }))
                    }),
                ),
            )
            .child(
                div()
                    .text_xs()
                    .text_color(rgb(0x96a6bc))
                    .child(self.model.description()),
            )
            .child(
                div().h_flex().flex_wrap().gap_2().children(
                    [
                        (0, self.n_acetyl, "N-terminal acetylation"),
                        (1, self.c_amide, "C-terminal amidation"),
                        (2, self.secondary, "Predict secondary structure"),
                    ]
                    .into_iter()
                    .map(|(i, selected, label)| {
                        Button::new(("protein-model-toggle", i as usize))
                            .small()
                            .selected(selected)
                            .disabled(self.calculation.busy())
                            .label(label)
                            .on_click(cx.listener(move |this, _, _, cx| {
                                if this.calculation.busy() {
                                    return;
                                }
                                match i {
                                    0 => this.n_acetyl = !this.n_acetyl,
                                    1 => this.c_amide = !this.c_amide,
                                    _ => this.secondary = !this.secondary,
                                }
                                this.invalidate(cx);
                            }))
                    }),
                ),
            );
        for (label, input) in self.fields.iter().skip(3) {
            body = body.child(
                div()
                    .v_flex()
                    .gap_1()
                    .child(div().text_xs().child(*label))
                    .child(Input::new(input).disabled(self.calculation.busy())),
            );
        }
        body = body.child(
            Button::new("calculate-protein-analysis")
                .label(if self.calculation.busy() {
                    "Calculating…"
                } else {
                    "Calculate"
                })
                .disabled(self.calculation.busy())
                .on_click(cx.listener(|this, _, window, cx| this.calculate(window, cx))),
        );
        if let Some(p) = self.calculation.current() {
            body = body
                .child(
                    div()
                        .h_flex()
                        .flex_wrap()
                        .gap_4()
                        .child(stat(
                            "AVERAGE MASS",
                            p.mass.map_or("Unavailable".into(), |m| {
                                format!("{:.2} Da", m.average_da)
                            }),
                        ))
                        .child(stat(
                            "MONOISOTOPIC MASS",
                            p.mass.map_or("Unavailable".into(), |m| {
                                format!("{:.4} Da", m.monoisotopic_da)
                            }),
                        ))
                        .child(stat(
                            "THEORETICAL pI",
                            p.theoretical_pi
                                .map_or("Unavailable".into(), |v| format!("{v:.3}")),
                        ))
                        .child(stat(
                            "NET CHARGE AT CHOSEN pH",
                            p.charge
                                .map_or("Unavailable".into(), |v| format!("{v:+.3}")),
                        ))
                        .child(stat(
                            "EXTINCTION 280 nm",
                            p.extinction_280
                                .map_or("Unavailable".into(), |v| format!("{v} M⁻¹cm⁻¹")),
                        ))
                        .child(stat(
                            "A280 · 1 mg/mL · 1 cm",
                            p.absorbance_1mg_ml
                                .map_or("Unavailable".into(), |v| format!("{v:.3}")),
                        )),
                )
                .child(
                    Button::new("copy-detailed-protein")
                        .small()
                        .label("Copy full report")
                        .on_click(cx.listener(|this, _, _, cx| {
                            if let Some(text) = this.calculation.report(&this.document.name) {
                                cx.write_to_clipboard(ClipboardItem::new_string(text));
                            }
                        })),
                );
            if !p.charge_curve.is_empty() {
                body = body.child(
                    Button::new("copy-protein-charge-curve")
                        .small()
                        .label("Copy pH/charge table")
                        .on_click(cx.listener(|this, _, _, cx| {
                            if let Some(curve) = this.calculation.curve() {
                                cx.write_to_clipboard(ClipboardItem::new_string(curve));
                            }
                        })),
                );
            }
            if let Some(reason) = &p.unmodified.unavailable_reason {
                body = body.child(div().text_color(rgb(0xe7bd7b)).child(reason.clone()));
            }
            for warning in &p.warnings {
                body = body.child(
                    div()
                        .text_xs()
                        .text_color(rgb(0x96a6bc))
                        .child(warning.clone()),
                );
            }
            if let Some(s) = &p.secondary {
                body=body.child(div().text_sm().child(format!("Secondary estimate: H {} · E {} · C {}",s.helix_count,s.sheet_count,s.coil_count)))
                    .child(div().text_xs().text_color(rgb(0xe7bd7b)).child(s.method))
                    .child(div().text_xs().text_color(rgb(0xe7bd7b)).child("This legacy sequence heuristic has no calibrated confidence or 3D coordinates. C includes unassigned turns; use experimental evidence or a validated modern predictor for structural conclusions."));
                for (i, chunk) in s.assignments.as_bytes().chunks(80).take(25).enumerate() {
                    body = body.child(div().font_family("monospace").text_xs().child(format!(
                        "{:>7}  {}",
                        p.range.start + i * 80 + 1,
                        String::from_utf8_lossy(chunk)
                    )));
                }
                if s.assignments.len() > 2000 {
                    body=body.child(div().text_xs().child("Showing the first 2,000 assignments. Copy full report for the complete result."));
                }
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
            Button::new("close-protein-tools")
                .ghost()
                .label("Close")
                .on_click(cx.listener(|_, _, _, cx| cx.emit(ProteinEvent::Cancel))),
        )
    }
}

#[cfg(test)]
mod detailed_tests {
    use super::{
        Calculation, Modification, PkaModel, ProteinOptions, detailed_report, parse_modifications,
    };
    fn result(model: PkaModel) -> super::DetailedProteinProperties {
        popgene_core::protein::analyze_with_options(
            "KGGD",
            &ProteinOptions {
                range: Some(1..3),
                model,
                ..ProteinOptions::default()
            },
        )
        .unwrap()
    }
    fn ready() -> Calculation {
        let mut state = Calculation::default();
        let token = state.begin().unwrap();
        state.complete(token, Ok(result(PkaModel::Grimsley2009)));
        state
    }

    #[test]
    fn copy_actions_resolve_current_inputs_instead_of_rendered_output() {
        let mut state = ready();
        assert!(
            state
                .report("scope")
                .unwrap()
                .contains("Theoretical pI: 5.500000")
        );
        assert!(
            state
                .curve()
                .unwrap()
                .starts_with("pH\tNet charge (Grimsley 2009 (IPC Arg))")
        );
        state.invalidate(); // Same path as all input Change events and parameter toggles.
        assert!(state.report("scope").is_none() && state.curve().is_none());
        let token = state.begin().unwrap();
        assert!(state.report("scope").is_none() && state.curve().is_none());
        state.complete(token, Ok(result(PkaModel::Emboss)));
        assert!(
            state
                .report("scope")
                .unwrap()
                .contains("Theoretical pI: 5.550000")
        );
        assert!(
            state
                .curve()
                .unwrap()
                .starts_with("pH\tNet charge (Original EMBOSS)")
        );
    }

    #[test]
    fn invalid_or_pending_calculations_cannot_keep_a_previous_report() {
        let mut state = ready();
        let token = state.begin().unwrap();
        assert!(state.busy());
        assert!(state.begin().is_none());
        assert!(state.current().is_none());
        state.complete(token, Err("Enter a numeric pH".into()));
        assert!(!state.busy());
        assert_eq!(state.error, "Enter a numeric pH");
        assert!(state.report("scope").is_none() && state.curve().is_none());
        state.invalidate();
        assert!(state.error.is_empty());
    }

    #[test]
    fn delayed_worker_completion_cannot_restore_invalidated_values_or_clear_new_work() {
        let mut state = Calculation::default();
        let old = state.begin().unwrap();
        state.invalidate(); // A delayed Change event while the worker is running.
        state.complete(old, Ok(result(PkaModel::Grimsley2009)));
        assert!(!state.busy());
        assert!(state.current().is_none());
        let next = state.begin().unwrap();
        state.complete(old, Err("obsolete failure".into()));
        assert_eq!(state.pending, Some(next));
        assert!(state.error.is_empty());
        state.complete(next, Ok(result(PkaModel::Bjellqvist)));
        assert_eq!(state.current().unwrap().model, PkaModel::Bjellqvist);
        assert!(!state.busy());
    }

    #[test]
    fn grimsley_report_keeps_model_provenance_selected_scope_and_numeric_result() {
        let p = popgene_core::protein::analyze_with_options(
            "KGGD",
            &ProteinOptions {
                range: Some(1..3),
                model: PkaModel::Grimsley2009,
                ..ProteinOptions::default()
            },
        )
        .unwrap();
        let report = detailed_report("Selected peptide", &p);
        for expected in [
            "Residues 2..3",
            "Grimsley 2009 (IPC Arg)",
            "Arg 12.04",
            "10.1002/pro.19",
            "10.1186/s13062-016-0159-9",
            "Theoretical pI: 5.500000",
            "Net charge: -0.166138",
            "root within pH 0–14",
            "new N/C termini",
        ] {
            assert!(report.contains(expected), "missing {expected}");
        }
    }
    #[test]
    fn modification_inputs_preserve_document_numbering_and_reject_bad_pairs() {
        assert_eq!(
            parse_modifications(true, false, "4", "2", "6", "3-9").unwrap(),
            vec![
                Modification::NAcetyl,
                Modification::LysAcetyl(3),
                Modification::MetOxidation(1),
                Modification::Phosphorylation(5),
                Modification::Disulfide(2, 8)
            ]
        );
        for invalid in ["0", "abc", "3-"] {
            assert!(parse_modifications(false, false, "", "", "", invalid).is_err());
        }
    }
}
