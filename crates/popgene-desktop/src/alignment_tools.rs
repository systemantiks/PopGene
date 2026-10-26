//! Bounded, asynchronous native alignment workspace. Input documents are snapshots;
//! a consensus is emitted as a new document and never overwrites an input read.
use gpui_kit::component::{
    Disableable, Selectable, Sizable, StyledExt,
    button::{Button, ButtonVariants},
    input::{Input, InputEvent, InputState, Textarea, TextareaState},
};
use gpui_kit::*;
use popgene_core::{
    Document, Molecule,
    alignment::{self, AlignmentMode, AlignmentOptions},
};

pub enum AlignmentEvent {
    CreateConsensus(Box<Document>),
    Cancel,
}
#[derive(Clone, Copy, PartialEq)]
enum Output {
    Alignment,
    Differences,
    Fasta,
}
struct ResultText {
    summary: String,
    alignment: String,
    differences: String,
    fasta: String,
    consensus: Option<Document>,
}
#[derive(Default)]
struct Calculation {
    generation: u64,
    pending: Option<u64>,
    result: Option<ResultText>,
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
    fn complete(&mut self, generation: u64, result: Result<ResultText, String>) {
        // An obsolete completion must not clear a newer calculation's busy state.
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
    fn current(&self) -> Option<&ResultText> {
        (!self.busy()).then_some(self.result.as_ref()).flatten()
    }
    fn text(&self, output: Output) -> Option<&str> {
        self.current().map(|result| match output {
            Output::Alignment => result.alignment.as_str(),
            Output::Differences => result.differences.as_str(),
            Output::Fasta => result.fasta.as_str(),
        })
    }
    fn consensus(&self) -> Option<&Document> {
        self.current().and_then(|result| result.consensus.as_ref())
    }
}
pub struct AlignmentTools {
    documents: Vec<Document>,
    included: Vec<bool>,
    reference: usize,
    options: AlignmentOptions,
    scores: Vec<Entity<InputState>>,
    text: Entity<TextareaState>,
    output: Output,
    calculation: Calculation,
    _subscriptions: Vec<Subscription>,
}
impl EventEmitter<AlignmentEvent> for AlignmentTools {}
impl AlignmentTools {
    pub fn error(&mut self, message: impl Into<String>, cx: &mut Context<Self>) {
        self.calculation.error = message.into();
        cx.notify();
    }
    pub fn new(documents: Vec<Document>, window: &mut Window, cx: &mut Context<Self>) -> Self {
        let included = vec![true; documents.len()];
        let scores = ["2", "-3", "-5", "-1"]
            .into_iter()
            .map(|value| cx.new(|cx| InputState::new(window, cx).default_value(value)))
            .collect::<Vec<_>>();
        let subscriptions = scores
            .iter()
            .map(|input| {
                cx.subscribe_in(input, window, |this: &mut Self, _, event, window, cx| {
                    if matches!(event, InputEvent::Change) {
                        this.invalidate(window, cx);
                    }
                })
            })
            .collect();
        let text = cx.new(|cx| TextareaState::new(window, cx).rows(20));
        Self {
            documents,
            included,
            reference: 0,
            options: AlignmentOptions::default(),
            scores,
            text,
            output: Output::Alignment,
            calculation: Calculation::default(),
            _subscriptions: subscriptions,
        }
    }
    fn invalidate(&mut self, window: &mut Window, cx: &mut Context<Self>) {
        self.calculation.invalidate();
        self.update_text(window, cx);
    }
    fn update_text(&mut self, window: &mut Window, cx: &mut Context<Self>) {
        let value = self
            .calculation
            .text(self.output)
            .unwrap_or_default()
            .to_owned();
        self.text
            .update(cx, |state, cx| state.set_value(value, window, cx));
        cx.notify();
    }
    fn run(&mut self, poa: bool, window: &mut Window, cx: &mut Context<Self>) {
        let Some(generation) = self.calculation.begin() else {
            return;
        };
        // Clear previous output before validation as well as asynchronous work.
        self.update_text(window, cx);
        let parsed = self
            .scores
            .iter()
            .map(|field| field.read(cx).value().trim().parse::<i32>())
            .collect::<Result<Vec<_>, _>>();
        let Ok(scores) = parsed else {
            self.calculation
                .complete(generation, Err("Scores must be whole numbers".into()));
            cx.notify();
            return;
        };
        let mut options = self.options;
        options.match_score = scores[0];
        options.mismatch_score = scores[1];
        options.gap_open = scores[2];
        options.gap_extend = scores[3];
        let mut documents = Vec::new();
        if let Some(reference) = self.documents.get(self.reference) {
            documents.push(reference.clone());
        }
        documents.extend(
            self.documents
                .iter()
                .enumerate()
                .filter(|(i, _)| *i != self.reference && self.included[*i])
                .map(|(_, d)| d.clone()),
        );
        if documents.len() < 2 {
            self.calculation.complete(
                generation,
                Err("Include at least two open documents".into()),
            );
            cx.notify();
            return;
        }
        cx.spawn_in(window, async move |this, cx| {
            let result = cx
                .background_executor()
                .spawn(async move { calculate(documents, options, poa) })
                .await;
            let _ = this.update_in(cx, |this, window, cx| {
                this.calculation.complete(generation, result);
                this.update_text(window, cx);
            });
        })
        .detach();
    }
}
fn calculate(
    documents: Vec<Document>,
    options: AlignmentOptions,
    poa: bool,
) -> Result<ResultText, String> {
    let molecule = documents[0].molecule;
    let names = documents
        .iter()
        .map(|d| d.name.clone())
        .collect::<Vec<_>>()
        .join(", ");
    let consensus_document = |sequence: &str, method: &str| -> Result<Document, String> {
        let mut doc = Document::new_typed(
            format!("{} consensus", documents[0].name),
            sequence,
            molecule,
        )
        .map_err(|e| e.to_string())?;
        doc.description = format!(
            "{method}. Inputs: {names}. Calculated from document snapshots; alignment columns and read coordinates remain in the alignment result."
        );
        Ok(doc)
    };
    if poa {
        let consensus = alignment::poa_consensus(&documents, options).map_err(|e| e.to_string())?;
        let summary = format!(
            "Partial-order consensus · {} sequences · {} residues · linear gap penalty {} per column",
            documents.len(),
            consensus.len(),
            options.gap_open
        );
        return Ok(ResultText {
            alignment: format!(
                "{summary}\n\n{consensus}\n\nProgressive graph consensus is input-order dependent. Query orientations use the selected pairwise mode; POA graph alignment is global. No quality weighting is applied.\nInputs: {names}"
            ),
            differences: "Use Align to inspect residue differences and read quality values.".into(),
            fasta: format!(">consensus\n{consensus}\n"),
            consensus: Some(consensus_document(&consensus, &summary)?),
            summary,
        });
    }
    if documents.len() == 2 {
        let pair = alignment::align_pair(&documents[0], &documents[1], options)
            .map_err(|e| e.to_string())?;
        let summary = format!(
            "{:?} · score {} · {:.1}% exact identity over {} paired residues · {} gap columns{}",
            options.mode,
            pair.score,
            pair.identity_percent(),
            pair.aligned_residues,
            pair.gaps,
            if pair.query.reverse {
                " · query reverse complement"
            } else {
                ""
            }
        );
        let reference_notice = if pair.crosses_origin {
            "\nThis alignment crosses the circular origin. Reference coordinates wrap to base 1; no reference base is used for more than one lap.\n"
        } else if documents[0].circular && options.mode == AlignmentMode::Global {
            "\nGlobal alignment compares complete sequences at their stored origins. Query in reference and Local modes search across the circular seam.\n"
        } else if documents[0].circular {
            "\nCircular reference searched across the origin, limited to one circumference.\n"
        } else {
            ""
        };
        Ok(ResultText {
            alignment: format!("{summary}{reference_notice}\n\n{}", pair.pretty(70)),
            differences: pair.differences(),
            fasta: pair.fasta(),
            consensus: None,
            summary,
        })
    } else {
        let result = alignment::align_multiple(&documents, options).map_err(|e| e.to_string())?;
        let summary = format!(
            "Reference-guided alignment · {} rows · {} columns · reference {}",
            result.rows.len(),
            result.rows[0].sequence.len(),
            documents[0].name
        );
        let mut differences = String::new();
        let mut linear_reference = documents[0].clone();
        linear_reference.circular = false;
        for query in &documents[1..] {
            let pair = alignment::align_pair(&linear_reference, query, options)
                .map_err(|e| e.to_string())?;
            differences.push_str(&format!(
                "{} vs {}\n{}\n\n",
                documents[0].name,
                query.name,
                pair.differences()
            ));
        }
        Ok(ResultText {
            alignment: format!(
                "{summary}\nShared insertions are left aligned; this reference-guided heuristic does not optimize all pairs. Consensus uses majority residues and IUPAC ties, without quality weighting. Circular inputs use their current origins.\n\n{}\nConsensus\n{}",
                result.pretty(70),
                result.consensus
            ),
            differences,
            fasta: result.fasta(),
            consensus: Some(consensus_document(&result.consensus, &summary)?),
            summary,
        })
    }
}
impl Render for AlignmentTools {
    fn render(&mut self, _: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        let mut controls =
            div()
                .id("alignment-controls")
                .v_flex()
                .gap_3()
                .w(px(280.))
                .min_h_0()
                .flex_shrink_0()
                .overflow_y_scroll()
                .child(div().text_lg().child("Align open sequences"))
                .child(div().text_xs().child(
                    "Choose a reference and included reads. The reference is always included.",
                ));
        for (index, document) in self.documents.iter().enumerate() {
            controls = controls.child(
                div()
                    .v_flex()
                    .gap_1()
                    .child(div().text_sm().child(format!(
                        "{} · {} {}",
                        document.name,
                        document.sequence.len(),
                        if document.molecule == Molecule::Protein {
                            "aa"
                        } else {
                            "nt"
                        }
                    )))
                    .child(
                        div()
                            .h_flex()
                            .gap_1()
                            .child(
                                Button::new(("alignment-reference", index))
                                    .small()
                                    .selected(self.reference == index)
                                    .label("Reference")
                                    .disabled(self.calculation.busy())
                                    .on_click(cx.listener(move |this, _, window, cx| {
                                        if !this.calculation.busy() && this.reference != index {
                                            this.reference = index;
                                            this.invalidate(window, cx);
                                        }
                                    })),
                            )
                            .child(
                                Button::new(("alignment-include", index))
                                    .small()
                                    .selected(self.included[index] || self.reference == index)
                                    .label("Include")
                                    .disabled(self.calculation.busy() || self.reference == index)
                                    .on_click(cx.listener(move |this, _, window, cx| {
                                        if !this.calculation.busy() && this.reference != index {
                                            this.included[index] = !this.included[index];
                                            this.invalidate(window, cx);
                                        }
                                    })),
                            ),
                    ),
            );
        }
        for (i, (mode, label)) in [
            (AlignmentMode::QueryInReference, "Query in reference"),
            (AlignmentMode::Global, "Global"),
            (AlignmentMode::Local, "Local (pair only)"),
        ]
        .into_iter()
        .enumerate()
        {
            controls = controls.child(
                Button::new(("alignment-mode", i))
                    .small()
                    .selected(self.options.mode == mode)
                    .label(label)
                    .disabled(self.calculation.busy())
                    .on_click(cx.listener(move |this, _, window, cx| {
                        if !this.calculation.busy() && this.options.mode != mode {
                            this.options.mode = mode;
                            this.invalidate(window, cx);
                        }
                    })),
            );
        }
        controls = controls
            .child(
                Button::new("alignment-both-strands")
                    .small()
                    .selected(self.options.both_strands)
                    .label("Search both nucleotide strands")
                    .disabled(self.calculation.busy())
                    .on_click(cx.listener(|this, _, window, cx| {
                        if !this.calculation.busy() {
                            this.options.both_strands = !this.options.both_strands;
                            this.invalidate(window, cx);
                        }
                    })),
            )
            .child(
                Button::new("alignment-blosum")
                    .small()
                    .selected(self.options.protein_blosum62)
                    .label("Protein BLOSUM62")
                    .disabled(self.calculation.busy())
                    .on_click(cx.listener(|this, _, window, cx| {
                        if !this.calculation.busy() {
                            this.options.protein_blosum62 = !this.options.protein_blosum62;
                            this.invalidate(window, cx);
                        }
                    })),
            );
        for (index, label) in ["Match", "Mismatch", "Gap opening", "Gap extension"]
            .into_iter()
            .enumerate()
        {
            controls = controls.child(
                div()
                    .h_flex()
                    .gap_2()
                    .child(div().w(px(120.)).text_sm().child(label))
                    .child(
                        Input::new(&self.scores[index])
                            .disabled(self.calculation.busy())
                            .w(px(110.)),
                    ),
            );
        }
        controls=controls.child(div().text_xs().child("Affine gap cost = opening + (length − 1) × extension. BLOSUM62 replaces match/mismatch scores for proteins."))
            .child(Button::new("alignment-run").primary().label(if self.calculation.busy() {"Calculating…"}else{"Align"}).disabled(self.calculation.busy()).on_click(cx.listener(|this,_,window,cx|this.run(false,window,cx))))
            .child(Button::new("alignment-poa").label("Calculate POA consensus").disabled(self.calculation.busy()).on_click(cx.listener(|this,_,window,cx|this.run(true,window,cx))))
            .child(Button::new("alignment-close").ghost().label("Close").on_click(cx.listener(|_,_,_,cx|cx.emit(AlignmentEvent::Cancel))));
        let mut output = div().v_flex().flex_1().min_w_0().min_h_0().gap_2();
        let mut toolbar = div().h_flex().gap_2();
        for (i, (mode, label)) in [
            (Output::Alignment, "Alignment"),
            (Output::Differences, "Differences / quality"),
            (Output::Fasta, "Gapped FASTA"),
        ]
        .into_iter()
        .enumerate()
        {
            toolbar = toolbar.child(
                Button::new(("alignment-output", i))
                    .small()
                    .selected(self.output == mode)
                    .label(label)
                    .on_click(cx.listener(move |this, _, window, cx| {
                        this.output = mode;
                        this.update_text(window, cx);
                    })),
            )
        }
        toolbar = toolbar.child(
            Button::new("alignment-copy")
                .small()
                .label("Copy")
                .disabled(self.calculation.current().is_none())
                .on_click(cx.listener(|this, _, _, cx| {
                    if let Some(text) = this.calculation.text(this.output) {
                        cx.write_to_clipboard(ClipboardItem::new_string(text.to_owned()));
                    }
                })),
        );
        if let Some(result) = self.calculation.current() {
            output = output.child(div().text_sm().child(result.summary.clone()));
            if result.consensus.is_some() {
                toolbar = toolbar.child(
                    Button::new("alignment-new-consensus")
                        .small()
                        .label("Open consensus")
                        .on_click(cx.listener(|this, _, _, cx| {
                            if let Some(document) = this.calculation.consensus() {
                                cx.emit(AlignmentEvent::CreateConsensus(Box::new(
                                    document.clone(),
                                )));
                            }
                        })),
                );
            }
        }
        if !self.calculation.error.is_empty() {
            output = output.child(
                div()
                    .text_color(rgb(0xf2a4a4))
                    .child(self.calculation.error.clone()),
            );
        }
        output = output.child(toolbar).child(
            Textarea::new(&self.text)
                .readonly(true)
                .font_family("monospace")
                .flex_1()
                .min_h_0()
                .aria_label("Alignment result with original residue coordinates"),
        );
        div()
            .h_flex()
            .items_stretch()
            .size_full()
            .min_h_0()
            .gap_4()
            .p_4()
            .child(controls)
            .child(output)
    }
}

#[cfg(test)]
mod calculation_tests {
    use super::{Calculation, Output, ResultText, calculate};
    use popgene_core::{
        Document,
        alignment::{AlignmentMode, AlignmentOptions},
    };

    fn documents() -> Vec<Document> {
        [
            ("reference", "ACGT"),
            ("read one", "ACGT"),
            ("read two", "ACCT"),
        ]
        .into_iter()
        .map(|(name, sequence)| Document::new(name, sequence).unwrap())
        .collect()
    }

    fn options() -> AlignmentOptions {
        AlignmentOptions {
            mode: AlignmentMode::Global,
            both_strands: false,
            ..AlignmentOptions::default()
        }
    }

    fn consensus() -> ResultText {
        calculate(documents(), options(), false).unwrap()
    }

    fn ready() -> Calculation {
        let mut state = Calculation::default();
        let generation = state.begin().unwrap();
        state.complete(generation, Ok(consensus()));
        assert_eq!(state.consensus().unwrap().sequence.as_ref(), "ACGT");
        state
    }

    #[test]
    fn invalidation_removes_every_output_and_consensus_before_recalculation() {
        let mut state = ready();
        state.error = "Prior action failed".into();
        state.invalidate();
        assert!(state.current().is_none());
        assert!(state.consensus().is_none());
        assert!(state.error.is_empty());
        for output in [Output::Alignment, Output::Differences, Output::Fasta] {
            assert!(state.text(output).is_none());
        }

        // Excluding a read changes the calculation from MSA to a pair: the
        // previous three-row consensus must no longer be an available action.
        let generation = state.begin().unwrap();
        state.complete(
            generation,
            calculate(documents()[..2].to_vec(), options(), false),
        );
        assert!(state.current().is_some());
        assert!(state.consensus().is_none());
        assert!(!state.text(Output::Fasta).unwrap().contains("read two"));

        // A changed reference creates a newly named consensus from that run.
        state.invalidate();
        let mut changed_reference = documents();
        changed_reference.swap(0, 1);
        let generation = state.begin().unwrap();
        state.complete(generation, calculate(changed_reference, options(), false));
        assert_eq!(state.consensus().unwrap().name, "read one consensus");
    }

    #[test]
    fn pending_and_failed_runs_never_expose_an_older_consensus() {
        let mut state = ready();
        let generation = state.begin().unwrap();
        assert!(state.busy());
        assert!(state.begin().is_none());
        assert!(state.current().is_none());
        assert!(state.consensus().is_none());
        assert!(state.text(Output::Fasta).is_none());
        state.complete(generation, Err("Scores must be whole numbers".into()));
        assert!(!state.busy());
        assert_eq!(state.error, "Scores must be whole numbers");
        assert!(state.consensus().is_none());

        let generation = state.begin().unwrap();
        assert!(state.error.is_empty());
        let invalid_options = AlignmentOptions {
            gap_open: 1,
            ..options()
        };
        let failed = calculate(documents(), invalid_options, false);
        assert!(failed.is_err());
        state.complete(generation, failed);
        assert!(!state.error.is_empty());
        assert!(state.current().is_none());
        assert!(state.consensus().is_none());

        let generation = state.begin().unwrap();
        state.complete(generation, calculate(documents(), options(), true));
        assert!(state.error.is_empty());
        assert_eq!(state.consensus().unwrap().sequence.as_ref(), "ACGT");
    }

    #[test]
    fn obsolete_completion_cannot_restore_output_or_clear_a_newer_run() {
        let mut state = Calculation::default();
        let obsolete = state.begin().unwrap();
        // Input events can invalidate a pending calculation even if controls
        // normally prevent mouse/keyboard edits while the worker is running.
        state.invalidate();
        assert!(state.busy());
        state.complete(obsolete, Ok(consensus()));
        assert!(!state.busy());
        assert!(state.current().is_none());
        assert!(state.consensus().is_none());

        let current = state.begin().unwrap();
        state.complete(obsolete, Err("Late failure".into()));
        assert!(state.busy());
        assert!(state.error.is_empty());
        state.complete(current, Ok(consensus()));
        assert!(!state.busy());
        assert!(state.consensus().is_some());
        state.complete(obsolete, Err("Still obsolete".into()));
        assert!(state.error.is_empty());
        assert!(state.consensus().is_some());
    }

    #[test]
    fn presentation_tabs_keep_the_same_current_result_and_consensus() {
        let state = ready();
        let generation = state.generation;
        let current = state.current().unwrap() as *const ResultText;
        assert!(state.text(Output::Alignment).unwrap().contains("3 rows"));
        assert!(
            state
                .text(Output::Differences)
                .unwrap()
                .contains("Query quality")
        );
        assert!(
            state
                .text(Output::Fasta)
                .unwrap()
                .contains(">read two\nACCT\n")
        );
        assert_eq!(state.generation, generation);
        assert_eq!(state.current().unwrap() as *const ResultText, current);
        assert_eq!(state.consensus().unwrap().sequence.as_ref(), "ACGT");
    }
}
