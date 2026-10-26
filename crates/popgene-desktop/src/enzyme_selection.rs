//! Explicit per-tab display lists and manually chosen original-source panels.
use crate::model::{EnzymeFilter, parse_enzyme_names, validate_enzyme_names};
use gpui_kit::component::{
    Sizable, StyledExt,
    button::{Button, ButtonVariants},
    input::{Input, InputEvent, InputState, Textarea, TextareaState},
};
use gpui_kit::{prelude::*, *};
use popgene_core::restriction::ENZYMES;
use popgene_formats::NativeEnzymePanels;
use std::collections::BTreeSet;

#[derive(Default)]
struct StoredPanels {
    source_name: Option<String>,
    loading: bool,
    contents: NativeEnzymePanels,
    error: Option<String>,
}
impl StoredPanels {
    fn new(source_name: Option<&str>) -> Self {
        Self {
            source_name: source_name.map(str::to_owned),
            loading: source_name.is_some(),
            ..Default::default()
        }
    }
    fn finish(&mut self, result: Result<NativeEnzymePanels, String>) {
        self.loading = false;
        match result {
            Ok(contents) => {
                self.contents = contents;
                self.error = None;
            }
            Err(error) => self.error = Some(error),
        }
    }
    /// Preserve every exact member, including unknown names and duplicates.
    /// An explicit empty panel is a blank draft, not "all catalogue enzymes".
    fn draft(&self, index: usize) -> Option<String> {
        self.contents
            .panels
            .get(index)
            .map(|panel| panel.enzyme_names.join("\n"))
    }
    fn parse_draft(
        &self,
        selected: Option<usize>,
        input: &str,
    ) -> Result<BTreeSet<String>, String> {
        let names = parse_enzyme_names(input)?;
        if let Some(panel) = selected.and_then(|index| self.contents.panels.get(index)) {
            // A source member such as "EcoRI,BamHI" is one unknown name.
            // Adding another enzyme must not reinterpret that unchanged member
            // as two valid names. Correcting or removing the member is explicit.
            // Native members use ASCII whitespace separators. Unicode spaces
            // inside a source name must not become valid draft separators.
            let tokens = input.split_ascii_whitespace().collect::<BTreeSet<_>>();
            let unchanged_members = panel
                .enzyme_names
                .iter()
                .filter(|name| tokens.contains(name.as_str()))
                .cloned()
                .collect();
            validate_enzyme_names(&unchanged_members)?;
        }
        Ok(names)
    }
}

pub enum EnzymeSelectionEvent {
    Apply(Option<BTreeSet<String>>),
    Cancel,
}
pub struct EnzymeSelection {
    names: Entity<TextareaState>,
    search: Entity<InputState>,
    filter: EnzymeFilter,
    stored: StoredPanels,
    selected_panel: Option<usize>,
    error: String,
    _subscriptions: Vec<Subscription>,
}
impl EventEmitter<EnzymeSelectionEvent> for EnzymeSelection {}
impl EnzymeSelection {
    pub fn new(
        selected: Option<&BTreeSet<String>>,
        filter: EnzymeFilter,
        original_source_name: Option<&str>,
        window: &mut Window,
        cx: &mut Context<Self>,
    ) -> Self {
        let value = selected
            .map(|names| names.iter().cloned().collect::<Vec<_>>().join("\n"))
            .unwrap_or_default();
        let names = cx.new(|cx| {
            TextareaState::new(window, cx)
                .rows(6)
                .placeholder("EcoRI\nBamHI")
                .default_value(value)
        });
        let search = cx
            .new(|cx| InputState::new(window, cx).placeholder("Search catalogue to add an enzyme"));
        let subscriptions = vec![
            cx.subscribe(&names, |this, _, event: &InputEvent, cx| {
                if matches!(event, InputEvent::Change) {
                    this.error.clear();
                    cx.notify();
                }
            }),
            cx.subscribe(&search, |_, _, _: &InputEvent, cx| cx.notify()),
        ];
        Self {
            names,
            search,
            filter,
            stored: StoredPanels::new(original_source_name),
            selected_panel: None,
            error: String::new(),
            _subscriptions: subscriptions,
        }
    }
    pub fn error(&mut self, message: impl Into<String>, cx: &mut Context<Self>) {
        self.error = message.into();
        cx.notify();
    }
    pub fn stored_panels_loaded(
        &mut self,
        result: Result<NativeEnzymePanels, String>,
        cx: &mut Context<Self>,
    ) {
        // Background completion never changes the current user draft or its
        // validation error. Only an explicit panel click replaces the names.
        self.stored.finish(result);
        cx.notify();
    }
    fn use_panel(&mut self, index: usize, window: &mut Window, cx: &mut Context<Self>) {
        if let Some(value) = self.stored.draft(index) {
            self.selected_panel = Some(index);
            self.names
                .update(cx, |state, cx| state.set_value(value, window, cx));
            self.error.clear();
            cx.notify();
        }
    }
    fn apply(&mut self, cx: &mut Context<Self>) {
        match self
            .stored
            .parse_draft(self.selected_panel, &self.names.read(cx).value())
        {
            Ok(names) => cx.emit(EnzymeSelectionEvent::Apply(Some(names))),
            Err(error) => self.error(error, cx),
        }
    }
    fn add(&mut self, name: &'static str, window: &mut Window, cx: &mut Context<Self>) {
        let current = self.names.read(cx).value();
        if !current
            .split(|c: char| c.is_whitespace() || matches!(c, ',' | ';'))
            .any(|s| s == name)
        {
            let value = if current.trim().is_empty() {
                name.to_string()
            } else {
                format!("{}\n{name}", current.trim_end())
            };
            self.names
                .update(cx, |state, cx| state.set_value(value, window, cx));
        }
        self.error.clear();
        cx.notify();
    }
}
impl Render for EnzymeSelection {
    fn render(&mut self, _: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        let query = self.search.read(cx).value().trim().to_ascii_lowercase();
        let matches = ENZYMES
            .iter()
            .filter(|enzyme| !query.is_empty() && enzyme.name.to_ascii_lowercase().contains(&query))
            .take(8)
            .collect::<Vec<_>>();
        let parsed = self
            .stored
            .parse_draft(self.selected_panel, &self.names.read(cx).value());
        let summary = match &parsed {
            Ok(names) if names.is_empty() => "Empty list: no enzymes will be displayed.".to_owned(),
            Ok(names) => format!("{} catalogue enzymes selected.", names.len()),
            Err(_) => "Unrecognized names must be corrected before applying.".to_owned(),
        };
        div().v_flex().gap_3()
            .child(div().text_xl().child("Choose displayed enzymes"))
            .child(div().text_sm().child("Use exact catalogue names separated by spaces, commas or new lines. This list applies to the map, sequence, enzyme table and figures."))
            .child(div().text_sm().child(format!("Additional filter: {}. Change it in the enzyme filter menu.", self.filter.label())))
            .when(self.stored.source_name.is_some(), |d| d.child(self.stored_panels(cx)))
            .child(Textarea::new(&self.names).aria_label("Displayed enzyme names").font_family("monospace"))
            .child(div().text_sm().child(summary))
            .child(Input::new(&self.search).small().aria_label("Search enzyme catalogue"))
            .child(div().h_flex().flex_wrap().gap_2().children(matches.iter().map(|enzyme| {
                let name = enzyme.name;
                Button::new(SharedString::from(format!("add-displayed-enzyme-{name}"))).small().label(name)
                    .on_click(cx.listener(move |this, _, window, cx| this.add(name, window, cx)))
            })))
            .when(!query.is_empty() && matches.is_empty(), |d| d.child(div().text_sm().child("No catalogue names match this search.")))
            .when(!self.error.is_empty(), |d| d.child(div().text_sm().text_color(rgb(0xf08c8c)).child(self.error.clone())))
            .child(div().h_flex().gap_2().justify_end()
                .child(Button::new("cancel-enzyme-selection").ghost().label("Cancel").on_click(cx.listener(|_,_,_,cx| cx.emit(EnzymeSelectionEvent::Cancel))))
                .child(Button::new("all-enzyme-selection").label("Use all catalogue enzymes").on_click(cx.listener(|_,_,_,cx| cx.emit(EnzymeSelectionEvent::Apply(None)))))
                .child(Button::new("apply-enzyme-selection").primary().label(if parsed.as_ref().is_ok_and(|names| names.is_empty()) { "Show no enzymes" } else { "Apply list" }).on_click(cx.listener(|this,_,_,cx| this.apply(cx)))))
    }
}

impl EnzymeSelection {
    fn stored_panels(&self, cx: &Context<Self>) -> impl IntoElement {
        let source = self.stored.source_name.as_deref().unwrap_or_default();
        div().v_flex().gap_2().min_w_0().flex_shrink_0()
            .child(div().text_sm().font_weight(FontWeight::MEDIUM).child("Stored panels from the original source"))
            .child(div().text_xs().min_w_0().truncate().child(format!("Original source: {source}")))
            .child(div().text_xs().child("Choose a panel to replace the draft below, then review and Apply. Stored panels do not identify the originally active selection."))
            .when(self.stored.loading, |d| d.child(div().text_sm().child("Reading stored panels…")))
            .when_some(self.stored.error.as_ref(), |d, error| d.child(div().text_sm().text_color(rgb(0xf08c8c)).child(format!("Stored panels could not be read: {error}"))))
            .when(!self.stored.loading && self.stored.error.is_none() && self.stored.contents.panels.is_empty(), |d| d.child(div().text_sm().child("No supported named panels were found in the original source.")))
            .when(!self.stored.contents.panels.is_empty() || !self.stored.contents.warnings.is_empty(), |d| d.child(
                div().id("stored-enzyme-panels-scroll").v_flex().gap_2().h(px(180.)).flex_shrink_0().overflow_y_scroll()
                    .children(self.stored.contents.panels.iter().enumerate().map(|(index, panel)| {
                        div().h_flex().gap_3().min_w_0().flex_shrink_0()
                            .child(div().v_flex().flex_1().min_w_0()
                                .child(div().text_sm().truncate().child(panel.name.clone()))
                                .child(div().text_xs().child(format!("{} members · source packet {}, panel {}", panel.enzyme_names.len(), panel.packet_index + 1, panel.panel_index + 1))))
                            .child(Button::new(("use-stored-enzyme-panel", index)).small().label("Use panel").tooltip(panel.name.clone())
                                .on_click(cx.listener(move |this, _, window, cx| this.use_panel(index, window, cx))))
                    }))
                    .children(self.stored.contents.warnings.iter().map(|warning| div().text_xs().flex_shrink_0().child(warning.clone())))
            ))
    }
}

#[cfg(test)]
mod tests {
    use super::{NativeEnzymePanels, StoredPanels, parse_enzyme_names};
    use popgene_core::Source;

    fn panels(xml: &str) -> NativeEnzymePanels {
        let bytes = [
            vec![14],
            (xml.len() as u32).to_be_bytes().to_vec(),
            xml.as_bytes().to_vec(),
        ]
        .concat();
        popgene_formats::native_enzyme_panels(&Source {
            file_name: "original.dna".into(),
            format: "SnapGene".into(),
            bytes: bytes.into(),
            retained_blocks: vec![],
        })
        .unwrap()
    }

    #[test]
    fn source_panel_drafts_preserve_unknown_members_for_apply_validation() {
        let mut stored = StoredPanels::new(Some("original.dna"));
        stored.finish(Ok(panels(r#"<CustomEnzymeSets><CustomEnzymeSet type="0" name="Needs correction" enzymeNames="EcoRI UnknownEnzyme EcoRI"/><CustomEnzymeSet type="0" name="Explicitly empty" enzymeNames=""/></CustomEnzymeSets>"#)));
        let draft = stored.draft(0).unwrap();
        assert_eq!(draft, "EcoRI\nUnknownEnzyme\nEcoRI");
        assert!(
            parse_enzyme_names(&draft)
                .unwrap_err()
                .contains("UnknownEnzyme")
        );
        assert_eq!(stored.draft(1).as_deref(), Some(""));
        assert!(
            parse_enzyme_names(&stored.draft(1).unwrap())
                .unwrap()
                .is_empty()
        );
        assert_eq!(stored.draft(2), None);
    }

    #[test]
    fn source_read_failure_is_separate_from_available_drafts() {
        let mut stored = StoredPanels::new(Some("original.dna"));
        assert!(stored.loading);
        stored.finish(Ok(panels(r#"<CustomEnzymeSets><CustomEnzymeSet type="0" name="Panel" enzymeNames="EcoRI BamHI"/></CustomEnzymeSets>"#)));
        let draft = stored.draft(0);
        stored.finish(Err("Malformed stored panel".into()));
        assert!(!stored.loading);
        assert_eq!(stored.error.as_deref(), Some("Malformed stored panel"));
        assert_eq!(stored.draft(0), draft);
        assert_eq!(stored.source_name.as_deref(), Some("original.dna"));
    }

    #[test]
    fn source_member_delimiters_require_explicit_correction() {
        let mut stored = StoredPanels::new(Some("original.dna"));
        stored.finish(Ok(panels(r#"<CustomEnzymeSets><CustomEnzymeSet type="0" name="Unknown exact name" enzymeNames="EcoRI,BamHI"/></CustomEnzymeSets>"#)));
        let unchanged = stored.draft(0).unwrap();
        assert!(
            stored
                .parse_draft(Some(0), &unchanged)
                .unwrap_err()
                .contains("EcoRI,BamHI")
        );
        assert!(stored.parse_draft(Some(0), "EcoRI,BamHI\nHindIII").is_err());
        assert_eq!(
            stored.parse_draft(Some(0), "EcoRI\nBamHI").unwrap().len(),
            2
        );
        // Manually entered comma-separated lists keep their established behavior.
        assert_eq!(stored.parse_draft(None, "EcoRI,BamHI").unwrap().len(), 2);
    }

    #[test]
    fn unicode_whitespace_inside_source_member_is_not_a_draft_separator() {
        let mut stored = StoredPanels::new(Some("original.dna"));
        stored.finish(Ok(panels(r#"<CustomEnzymeSets><CustomEnzymeSet type="0" name="Unicode member" enzymeNames="EcoRI&#160;BamHI"/></CustomEnzymeSets>"#)));
        let unchanged = stored.draft(0).unwrap();
        assert_eq!(unchanged, "EcoRI\u{a0}BamHI");
        assert!(stored.parse_draft(Some(0), &unchanged).is_err());
        assert!(
            stored
                .parse_draft(Some(0), &format!("{unchanged}\nHindIII"))
                .is_err()
        );
        assert_eq!(
            stored.parse_draft(Some(0), "EcoRI\nBamHI").unwrap().len(),
            2
        );
    }
}
