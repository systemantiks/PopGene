//! Native GPUI workspace. Scientific state lives in `model`, independent of the UI.
use crate::{
    enzyme_sidebar::{Details, EnzymeSidebar, EnzymeSidebarEvent},
    feature_card,
    forms::{
        AnnotationEvent, AnnotationForm, Creator, CreatorEvent, EditEvent, EditForm, ExportEvent,
        ExportForm,
    },
    gene_map::{GeneMap, MapEvent, MapOptions, MapRestrictionSite},
    model::{Close, EnzymeFilter, Model, View},
    sequence_view::{
        SequenceEvent, SequenceOptions, SequenceRestrictionSite, SequenceView, TraceView,
    },
};
use gpui_kit::component::{
    button::*,
    hover_card::HoverCard,
    input::{Input, InputEvent, InputState, Textarea, TextareaState},
    menu::{DropdownMenu, PopupMenuItem},
    tab::{Tab as ViewTab, TabBar},
    *,
};
use gpui_kit::{assets::IconName as AssetIcon, prelude::*, *};
use popgene_core::{Document, Feature, Molecule, Segment, restriction::ENZYMES};
use popgene_formats::{SaveFormat, save_warnings};
use std::{
    collections::{HashMap, HashSet},
    path::PathBuf,
    time::Duration,
};

pub const BG: u32 = 0x101722;
pub const PANEL: u32 = 0x17212e;
pub const BORDER: u32 = 0x2a3646;
pub const INK: u32 = 0xe4edf4;
pub const MUTED: u32 = 0x8999ad;
pub const ACCENT: u32 = 0x48c9bb;

actions!(
    popgene,
    [
        NewSequence,
        OpenFiles,
        OpenFolder,
        Save,
        SaveAs,
        Export,
        CloseTab,
        Quit,
        CheckForUpdates,
        Undo,
        Redo,
        NextTab,
        PreviousTab,
        FindSequence,
        Escape
    ]
);

struct TabUi {
    creator: Option<Entity<Creator>>,
    map: Entity<GeneMap>,
    sequence: Entity<SequenceView>,
    sequence_options: SequenceOptions,
    trace: Entity<TraceView>,
    split: Entity<ResizableState>,
    start: Entity<InputState>,
    end: Entity<InputState>,
    search: Entity<InputState>,
    filter: Entity<InputState>,
    description: Entity<TextareaState>,
    description_saved: String,
    feature_details: Entity<TextareaState>,
    feature_details_saved: Option<(Feature, Molecule)>,
    protein_properties: Option<(
        std::sync::Arc<str>,
        popgene_core::protein::ProteinProperties,
    )>,
    source_metadata: Option<Entity<TextareaState>>,
    source_metadata_loaded: bool,
    source_metadata_error: Option<String>,
    selected: Vec<Segment>,
    range_synced: (String, String),
}
enum Modal {
    Edit(Entity<EditForm>),
    Annotation(Entity<AnnotationForm>),
    Export(Entity<ExportForm>),
    Report,
    About,
    UpdateInfo(String),
}
enum SaveIntent {
    Save {
        tab_id: u64,
        close: bool,
        format: SaveFormat,
    },
    Export {
        tab_id: u64,
        format: SaveFormat,
    },
    Original {
        tab_id: u64,
    },
}

pub struct PopGene {
    pub model: Model,
    ui: HashMap<u64, TabUi>,
    focus: FocusHandle,
    modal_focus: FocusHandle,
    modal_was_open: bool,
    project_filter: Entity<InputState>,
    collapsed: HashSet<PathBuf>,
    navigator: bool,
    modal: Option<Modal>,
    polling: bool,
    file_dialog_open: bool,
    creator_dialog_draft: Option<(u64, bool)>,
    tabs_scroll: ScrollHandle,
    shown_tab: Option<u64>,
    enzyme_sidebar: Option<Entity<EnzymeSidebar>>,
    enzyme_source: Option<(u64, usize, Option<usize>)>,
    enzyme_revision: Option<(
        std::sync::Arc<popgene_core::restriction::RestrictionAnalysis>,
        String,
    )>,
}

impl PopGene {
    pub fn new(paths: Vec<PathBuf>, window: &mut Window, cx: &mut Context<Self>) -> Self {
        let filter = cx.new(|cx| InputState::new(window, cx).placeholder("Filter project files…"));
        cx.subscribe(&filter, |_, _, _: &InputEvent, cx| cx.notify())
            .detach();
        let focus = cx.focus_handle();
        focus.focus(window, cx);
        let mut this = Self {
            model: Model::new(),
            ui: HashMap::new(),
            focus,
            modal_focus: cx.focus_handle(),
            modal_was_open: false,
            project_filter: filter,
            collapsed: HashSet::new(),
            navigator: true,
            modal: None,
            polling: false,
            file_dialog_open: false,
            creator_dialog_draft: None,
            tabs_scroll: ScrollHandle::new(),
            shown_tab: None,
            enzyme_sidebar: None,
            enzyme_source: None,
            enzyme_revision: None,
        };
        this.model.open_paths(paths);
        this.wake(cx);
        #[cfg(target_os = "macos")]
        this.watch_native_quit(cx);
        this
    }
    #[cfg(target_os = "macos")]
    fn watch_native_quit(&mut self, cx: &mut Context<Self>) {
        if !crate::updater::install_quit_guard() {
            return;
        }
        cx.spawn(async move |this, cx| {
            loop {
                cx.background_executor()
                    .timer(Duration::from_millis(250))
                    .await;
                let Some(requested) = crate::updater::take_quit_request() else {
                    if this.upgrade().is_none() {
                        break;
                    }
                    continue;
                };
                if this
                    .update(cx, |this, cx| {
                        if requested {
                            cx.activate(true);
                            this.request_quit(cx);
                            // Sparkle's status window may still be frontmost.
                            // Raise the GPUI workspace so its save prompt is visible.
                            cx.defer(|cx| {
                                for handle in cx.windows() {
                                    let _ = cx.update_window(handle, |_, window, _| {
                                        window.activate_window();
                                    });
                                }
                            });
                        } else {
                            // Installation can fail after approval (for example,
                            // authorization is denied). Future edits need a fresh
                            // save/discard decision before any subsequent quit.
                            this.model.allow_close = false;
                            this.cancel_close();
                            cx.notify();
                        }
                    })
                    .is_err()
                {
                    crate::updater::reply_to_quit(false);
                    break;
                }
            }
        })
        .detach();
    }
    fn check_for_updates(&mut self, cx: &mut Context<Self>) {
        if self.blocked() {
            self.model.status = "Finish the open dialog before checking for updates.".into();
            cx.notify();
            return;
        }
        let this = cx.entity().downgrade();
        // Sparkle can create native windows. Do so outside a GPUI entity borrow.
        cx.defer(move |cx| {
            if let Err(message) = crate::updater::check_for_updates() {
                let _ = this.update(cx, |this, cx| {
                    this.modal = Some(Modal::UpdateInfo(message));
                    cx.notify();
                });
            }
        });
    }
    fn finish_quit(&self, cx: &mut Context<Self>) {
        if self.model.allow_close && !crate::updater::reply_to_quit(true) {
            cx.quit();
        }
    }
    fn cancel_close(&mut self) {
        self.model.cancel_close();
        crate::updater::reply_to_quit(false);
    }
    pub fn request_quit(&mut self, cx: &mut Context<Self>) {
        if self.blocked() {
            crate::updater::reply_to_quit(false);
            return;
        }
        if !self.model.request_close(Close::Quit) {
            crate::updater::reply_to_quit(false);
        }
        self.finish_quit(cx);
        cx.notify();
    }
    fn wake(&mut self, cx: &mut Context<Self>) {
        cx.notify();
        if self.polling || !self.model.has_background_work() {
            return;
        }
        self.polling = true;
        cx.spawn(async move |this, cx| {
            loop {
                cx.background_executor()
                    .timer(Duration::from_millis(50))
                    .await;
                let again = this
                    .update(cx, |this, cx| {
                        if this.model.poll() {
                            cx.notify();
                        }
                        this.finish_quit(cx);
                        let busy = this.model.has_background_work();
                        if !busy {
                            this.polling = false;
                        }
                        busy
                    })
                    .unwrap_or(false);
                if !again {
                    break;
                }
            }
        })
        .detach();
    }
    fn blocked(&self) -> bool {
        self.file_dialog_open || self.modal.is_some() || self.model.pending_close.is_some()
    }
    fn mutable(&self) -> bool {
        !self.blocked() && self.model.current().is_some_and(|t| !t.busy() && !t.draft)
    }
    fn open(&mut self, cx: &mut Context<Self>) {
        if self.blocked() {
            return;
        }
        self.file_dialog_open = true;
        cx.notify();
        let response = crate::open_dialog::prompt(cx);
        cx.spawn(async move |this, cx| {
            let result = response.await;
            let _ = this.update(cx, |this, cx| {
                this.file_dialog_open = false;
                match result {
                    Ok(Some(paths)) => this.model.open_paths(paths),
                    Ok(None) => {}
                    Err(error) => {
                        this.model.status = error.clone();
                        this.model.issues.push(error);
                    }
                }
                this.wake(cx);
            });
        })
        .detach();
    }
    fn folder(&mut self, cx: &mut Context<Self>) {
        if self.blocked() {
            return;
        }
        self.file_dialog_open = true;
        cx.notify();
        cx.spawn(async move |this, cx| {
            let path = rfd::AsyncFileDialog::new().pick_folder().await;
            let _ = this.update(cx, |this, cx| {
                this.file_dialog_open = false;
                if let Some(path) = path {
                    this.model.open_folder(path.path().to_owned());
                }
                this.wake(cx);
            });
        })
        .detach();
    }
    fn save(&mut self, save_as: bool, cx: &mut Context<Self>) {
        if !self.mutable() {
            return;
        }
        let i = self.model.active;
        let t = &self.model.tabs[i];
        let format = t
            .save_format
            .unwrap_or_else(|| SaveFormat::recommended(&t.editor.document));
        let path = (!save_as).then(|| t.path.clone()).flatten();
        if let Some(path) = &path
            && format.supports(&t.editor.document)
            && save_warnings(&t.editor.document, format).is_empty()
        {
            self.model.save_to(i, path.clone(), format);
            self.wake(cx);
            return;
        }
        self.format_dialog(i, format, false, false, path, cx);
    }
    fn save_before_close(&mut self, cx: &mut Context<Self>) {
        if self.file_dialog_open || self.model.waiting_save.is_some() {
            return;
        }
        let Some(i) = self.model.close_target() else {
            return;
        };
        let t = &self.model.tabs[i];
        if t.draft || t.busy() {
            return;
        }
        let format = t
            .save_format
            .unwrap_or_else(|| SaveFormat::recommended(&t.editor.document));
        let path = t.path.clone();
        if let Some(path) = &path
            && format.supports(&t.editor.document)
            && save_warnings(&t.editor.document, format).is_empty()
        {
            self.model.save_and_close(path.clone(), format);
            self.wake(cx);
            return;
        }
        self.format_dialog(i, format, false, true, path, cx);
    }
    fn format_dialog(
        &mut self,
        index: usize,
        format: SaveFormat,
        export: bool,
        close: bool,
        path: Option<PathBuf>,
        cx: &mut Context<Self>,
    ) {
        let tab = &self.model.tabs[index];
        let id = tab.id;
        // Unsupported current formats must offer conversion instead of a disabled Save.
        let path = path.filter(|_| format.supports(&tab.editor.document));
        let form = cx
            .new(|_| ExportForm::new(tab.editor.document.clone(), format, export, path.is_none()));
        cx.subscribe(&form, move |this, _, event, cx| {
            if this.file_dialog_open {
                return;
            }
            match event {
                ExportEvent::Cancel => this.modal = None,
                ExportEvent::Choose(format) => {
                    let Some(i) = this.model.tabs.iter().position(|t| t.id == id && !t.busy())
                    else {
                        return;
                    };
                    if let Some(path) = &path {
                        this.modal = None;
                        if close {
                            if this.model.close_target() == Some(i) {
                                this.model.save_and_close(path.clone(), *format);
                            }
                        } else {
                            this.model.save_to(i, path.clone(), *format);
                        }
                        this.wake(cx);
                    } else {
                        let doc = &this.model.tabs[i].editor.document;
                        let ext = format.extension_for(doc);
                        let filename = format!("{}.{}", safe_name(&doc.name), ext);
                        let intent = if export {
                            SaveIntent::Export {
                                tab_id: id,
                                format: *format,
                            }
                        } else {
                            SaveIntent::Save {
                                tab_id: id,
                                close,
                                format: *format,
                            }
                        };
                        this.choose_save_path(intent, filename, Some(ext), cx);
                    }
                }
            }
            cx.notify();
        })
        .detach();
        self.modal = Some(Modal::Export(form));
        cx.notify();
    }
    fn recover_source(&mut self, cx: &mut Context<Self>) {
        if !self.mutable() {
            return;
        }
        let t = self.model.current().unwrap();
        let Some(source) = &t.original_source else {
            return;
        };
        self.choose_save_path(
            SaveIntent::Original { tab_id: t.id },
            source.file_name.clone(),
            None,
            cx,
        );
    }
    fn choose_save_path(
        &mut self,
        intent: SaveIntent,
        filename: String,
        extension: Option<&'static str>,
        cx: &mut Context<Self>,
    ) {
        if self.file_dialog_open {
            return;
        }
        self.file_dialog_open = true;
        cx.notify();
        cx.spawn(async move |this, cx| {
            let mut dialog = rfd::AsyncFileDialog::new().set_file_name(filename);
            if let Some(extension) = extension {
                dialog = dialog.add_filter("Sequence document", &[extension]);
            }
            let path = dialog.save_file().await;
            let _ = this.update(cx, |this, cx| {
                this.file_dialog_open = false;
                if let Some(path) = path {
                    let tab_id = match intent {
                        SaveIntent::Save { tab_id, .. }
                        | SaveIntent::Export { tab_id, .. }
                        | SaveIntent::Original { tab_id } => tab_id,
                    };
                    // Imports may finish while the sheet is open, changing active/index.
                    // Always resolve the original document identity on completion.
                    if let Some(i) = this
                        .model
                        .tabs
                        .iter()
                        .position(|tab| tab.id == tab_id && !tab.busy())
                    {
                        let path = path.path().to_owned();
                        match intent {
                            SaveIntent::Save {
                                close: false,
                                format,
                                ..
                            } => {
                                this.modal = None;
                                this.model.save_to(i, path, format);
                            }
                            SaveIntent::Save {
                                close: true,
                                format,
                                ..
                            } => {
                                this.modal = None;
                                if this.model.close_target() == Some(i) {
                                    this.model.save_and_close(path, format);
                                }
                            }
                            SaveIntent::Export { format, .. } => {
                                this.model.export_to(i, path, format);
                                this.modal = None;
                            }
                            SaveIntent::Original { .. } => {
                                this.model.export_original_to(i, path);
                            }
                        }
                    }
                }
                this.wake(cx);
            });
        })
        .detach();
    }
    fn close_tab(&mut self, id: u64, cx: &mut Context<Self>) {
        if self.blocked() {
            return;
        }
        self.model.request_close(Close::Tab(id));
        cx.notify();
    }
    fn new_sequence(&mut self, cx: &mut Context<Self>) {
        if !self.blocked() {
            self.model.new_sequence();
            cx.notify();
        }
    }
    fn switch(&mut self, index: usize, window: &mut Window, cx: &mut Context<Self>) {
        if !self.blocked() {
            self.model.active = index;
            self.focus.focus(window, cx);
            cx.notify();
        }
    }
    fn inspect_enzyme(
        &mut self,
        tab_id: u64,
        enzyme_index: usize,
        selected_start: Option<usize>,
        cx: &mut Context<Self>,
    ) {
        if self.blocked() {
            return;
        }
        if enzyme_index >= ENZYMES.len()
            || !self
                .model
                .current()
                .is_some_and(|tab| tab.id == tab_id && tab.restriction.is_some())
        {
            return;
        }
        self.enzyme_source = Some((tab_id, enzyme_index, selected_start));
        self.enzyme_revision = None;
        cx.notify();
    }
    fn close_enzyme_sidebar(&mut self) {
        self.enzyme_sidebar = None;
        self.enzyme_source = None;
        self.enzyme_revision = None;
    }
    fn sync_enzyme_sidebar(&mut self, window: &mut Window, cx: &mut Context<Self>) {
        let Some((tab_id, enzyme_index, start)) = self.enzyme_source else {
            return;
        };
        let Some(tab) = self
            .model
            .current()
            .filter(|tab| tab.id == tab_id && !tab.draft)
        else {
            self.close_enzyme_sidebar();
            if !self.blocked() {
                self.focus.focus(window, cx);
            }
            return;
        };
        let had_details = self.enzyme_revision.is_some();
        // A prior site is no longer selected after an edit or a different selection.
        let start = start.filter(|start| {
            tab.selected_enzyme == Some(enzyme_index)
                && tab
                    .selected
                    .first()
                    .is_some_and(|part| part.start == *start)
        });
        if self.enzyme_source != Some((tab_id, enzyme_index, start)) {
            self.enzyme_source = Some((tab_id, enzyme_index, start));
            self.enzyme_revision = None;
        }
        // Hide old positions while edits trigger fresh background analysis.
        let Some(analysis) = &tab.restriction else {
            self.enzyme_revision = None;
            if had_details && !self.blocked() {
                self.focus.focus(window, cx);
            }
            return;
        };
        if self
            .enzyme_revision
            .as_ref()
            .is_some_and(|(previous, name)| {
                std::sync::Arc::ptr_eq(previous, analysis) && *name == tab.editor.document.name
            })
        {
            return;
        }
        let Some(details) = Details::from_document(
            &tab.editor.document,
            analysis,
            enzyme_index,
            start,
            tab.restriction_ambiguous,
        ) else {
            return;
        };
        self.enzyme_revision = Some((analysis.clone(), tab.editor.document.name.clone()));
        if let Some(sidebar) = &self.enzyme_sidebar {
            sidebar.update(cx, |sidebar, cx| sidebar.replace(details, cx));
        } else {
            let sidebar = cx.new(|cx| EnzymeSidebar::new(details, window, cx));
            cx.subscribe_in(&sidebar, window, |this, _, event, window, cx| {
                if this.blocked() {
                    return;
                }
                match event {
                    EnzymeSidebarEvent::Close => {
                        this.close_enzyme_sidebar();
                        this.focus.focus(window, cx);
                    }
                    EnzymeSidebarEvent::Inspect(index) => {
                        if let Some((id, _, _)) = this.enzyme_source {
                            this.inspect_enzyme(id, *index, None, cx);
                        }
                    }
                    EnzymeSidebarEvent::SelectSite(start) => {
                        if let Some((id, index, _)) = this.enzyme_source
                            && let Some(tab) = this.model.current_mut().filter(|tab| tab.id == id)
                            && let Some(site) = tab.restriction.as_ref().and_then(|a| {
                                a.sites
                                    .iter()
                                    .find(|site| site.enzyme_index == index && site.start == *start)
                            })
                        {
                            let parts = tab.site_parts(site);
                            tab.select(parts, None, true);
                            tab.selected_enzyme = Some(index);
                            tab.show_sequence();
                            this.enzyme_source = Some((id, index, Some(*start)));
                            this.enzyme_revision = None;
                        }
                    }
                }
                cx.notify();
            })
            .detach();
            self.enzyme_sidebar = Some(sidebar);
        }
    }
    fn cycle(&mut self, forward: bool, window: &mut Window, cx: &mut Context<Self>) {
        let n = self.model.tabs.len();
        if n > 0 && !self.blocked() {
            let i = (self.model.active + if forward { 1 } else { n - 1 }) % n;
            self.switch(i, window, cx);
        }
    }
    fn sync_ui(&mut self, window: &mut Window, cx: &mut Context<Self>) {
        self.sync_enzyme_sidebar(window, cx);
        let blocked = self.blocked();
        if blocked != self.modal_was_open {
            if blocked {
                self.modal_focus.focus(window, cx);
            } else {
                self.focus.focus(window, cx);
            }
            self.modal_was_open = blocked;
        }
        self.ui
            .retain(|id, _| self.model.tabs.iter().any(|t| t.id == *id));
        for i in 0..self.model.tabs.len() {
            let t = &self.model.tabs[i];
            let id = t.id;
            if let std::collections::hash_map::Entry::Vacant(entry) = self.ui.entry(id) {
                let doc = t.editor.document.clone();
                let creator = t.draft.then(|| cx.new(|cx| Creator::new(window, cx)));
                if let Some(creator) = &creator {
                    cx.subscribe_in(creator, window, move |this, _, event, window, cx| {
                        let Some(i) = this.model.tabs.iter().position(|t| t.id == id) else {
                            return;
                        };
                        match event {
                            CreatorEvent::Dirty(dirty) => {
                                if let Some((pending_id, saved_dirty)) =
                                    &mut this.creator_dialog_draft
                                    && *pending_id == id
                                {
                                    *saved_dirty = *dirty;
                                } else {
                                    this.model.tabs[i].draft_dirty = *dirty;
                                }
                            }
                            CreatorEvent::Dialog(open) => {
                                this.file_dialog_open = *open;
                                if *open {
                                    // Preserve an untouched chooser if an earlier queued import
                                    // completes while its sequence-text sheet is open.
                                    this.creator_dialog_draft =
                                        Some((id, this.model.tabs[i].draft_dirty));
                                    this.model.tabs[i].draft_dirty = true;
                                } else if let Some((pending_id, dirty)) =
                                    this.creator_dialog_draft.take()
                                    && pending_id == id
                                {
                                    this.model.tabs[i].draft_dirty = dirty;
                                }
                            }
                            CreatorEvent::Created(doc) => {
                                this.model.active = i;
                                this.model.create_in_active((**doc).clone());
                                this.focus.focus(window, cx);
                            }
                            CreatorEvent::Open => this.open(cx),
                            CreatorEvent::Folder => this.folder(cx),
                            CreatorEvent::Cancel => this.close_tab(id, cx),
                        }
                        cx.notify();
                    })
                    .detach();
                }
                let map = cx.new(|cx| GeneMap::new(doc.clone(), window, cx));
                cx.subscribe(&map, move |this, _, event, cx| {
                    if this.blocked() {
                        return;
                    }
                    let mut inspect = None;
                    if let Some(t) = this.model.tabs.iter_mut().find(|t| t.id == id) {
                        match event {
                            MapEvent::Select {
                                parts,
                                feature,
                                show_sequence,
                            } => {
                                t.select(parts.clone(), Some(*feature), true);
                                if *show_sequence {
                                    t.show_sequence();
                                }
                            }
                            MapEvent::SelectRestriction {
                                parts,
                                enzyme_index,
                            } => {
                                t.select(parts.clone(), None, true);
                                t.selected_enzyme = Some(*enzyme_index);
                                inspect =
                                    Some((*enzyme_index, parts.first().map(|part| part.start)));
                            }
                        }
                    }
                    if let Some((enzyme_index, start)) = inspect {
                        this.inspect_enzyme(id, enzyme_index, start, cx);
                    }
                    cx.notify();
                })
                .detach();
                let sequence = cx.new(|cx| SequenceView::new(doc.clone(), vec![], window, cx));
                cx.subscribe(&sequence, move |this, _, event, cx| {
                    if this.blocked() {
                        return;
                    }
                    let mut inspect = None;
                    if let Some(t) = this.model.tabs.iter_mut().find(|t| t.id == id) {
                        match event {
                            SequenceEvent::Select(parts) => t.select(parts.clone(), None, false),
                            SequenceEvent::SelectFeature(feature_id) => {
                                if let Some(feature) = t
                                    .editor
                                    .document
                                    .features
                                    .iter()
                                    .find(|f| f.id == *feature_id)
                                {
                                    t.select(feature.parts.clone(), Some(*feature_id), false);
                                }
                            }
                            SequenceEvent::InspectEnzyme {
                                enzyme_index,
                                start,
                            } => {
                                if let Some(site) = t.restriction.as_ref().and_then(|analysis| {
                                    analysis.sites.iter().find(|site| {
                                        site.enzyme_index == *enzyme_index && site.start == *start
                                    })
                                }) {
                                    let parts = t.site_parts(site);
                                    t.select(parts, None, true);
                                    t.selected_enzyme = Some(*enzyme_index);
                                    inspect = Some((*enzyme_index, Some(*start)));
                                }
                            }
                        }
                    }
                    if let Some((enzyme_index, start)) = inspect {
                        this.inspect_enzyme(id, enzyme_index, start, cx);
                    }
                    cx.notify();
                })
                .detach();
                let trace = cx.new(|cx| TraceView::new(doc.clone(), 1.0, cx));
                let start = input(&t.range_start, "From", window, cx);
                let end = input(&t.range_end, "To", window, cx);
                let search = input("", "Find exact sequence…", window, cx);
                cx.subscribe(&search, move |this, state, event, cx| {
                    if this.blocked() {
                        return;
                    }
                    if matches!(event, InputEvent::Change | InputEvent::PressEnter { .. }) {
                        let query = state.read(cx).value().to_string();
                        if let Some(t) = this.model.tabs.iter_mut().find(|t| t.id == id) {
                            if t.query == query && matches!(event, InputEvent::PressEnter { .. }) {
                                if matches!(event, InputEvent::PressEnter { shift: true, .. }) {
                                    t.previous_hit();
                                } else {
                                    t.next_hit();
                                }
                            } else {
                                t.find(query);
                            }
                        }
                        cx.notify();
                    }
                })
                .detach();
                let filter = input("", "Filter annotations…", window, cx);
                cx.subscribe(&filter, |_, _, _: &InputEvent, cx| cx.notify())
                    .detach();
                let description = cx.new(|cx| {
                    TextareaState::new(window, cx).default_value(doc.description.clone())
                });
                let feature_details = cx.new(|cx| TextareaState::new(window, cx));
                cx.subscribe(&description, move |this, state, event, cx| {
                    if this.blocked() || !matches!(event, InputEvent::Change) {
                        return;
                    }
                    let value = state.read(cx).value().to_string();
                    if let Some(t) = this.model.tabs.iter_mut().find(|t| t.id == id)
                        && !t.busy()
                    {
                        let _ = t.editor.apply(|d| {
                            d.description = value.clone();
                            Ok(())
                        });
                        if let Some(ui) = this.ui.get_mut(&id) {
                            ui.description_saved = value;
                        }
                        t.message("Description updated", false);
                        cx.notify();
                    }
                })
                .detach();
                entry.insert(TabUi {
                    creator,
                    map,
                    sequence,
                    sequence_options: SequenceOptions::default(),
                    trace,
                    split: cx.new(|_| ResizableState::default()),
                    start,
                    end,
                    search,
                    filter,
                    description,
                    description_saved: doc.description,
                    feature_details,
                    feature_details_saved: None,
                    protein_properties: None,
                    source_metadata: None,
                    source_metadata_loaded: false,
                    source_metadata_error: None,
                    selected: vec![],
                    range_synced: (t.range_start.clone(), t.range_end.clone()),
                });
            }
        }
        if let Some(t) = self.model.current_mut() {
            let u = self.ui.get_mut(&t.id).unwrap();
            if !t.draft {
                u.creator = None;
            }
            if let Some(creator) = &u.creator {
                creator.update(cx, |creator, cx| creator.set_suspended(blocked, cx));
                t.draft_dirty = creator.read(cx).is_dirty();
            }
            let range = (t.range_start.clone(), t.range_end.clone());
            if u.selected != t.selected || u.range_synced != range {
                u.selected = t.selected.clone();
                u.range_synced = range;
                u.start
                    .update(cx, |s, cx| s.set_value(t.range_start.clone(), window, cx));
                u.end
                    .update(cx, |s, cx| s.set_value(t.range_end.clone(), window, cx));
            }
            if u.search.read(cx).value().as_ref() != t.query {
                u.search
                    .update(cx, |state, cx| state.set_value(t.query.clone(), window, cx));
            }
            if u.description_saved != t.editor.document.description {
                u.description_saved = t.editor.document.description.clone();
                if u.description.read(cx).value().as_ref() != u.description_saved {
                    u.description.update(cx, |s, cx| {
                        s.set_value(u.description_saved.clone(), window, cx)
                    });
                }
            }
            if matches!(t.view, View::Features | View::Primers) {
                let selected = t.editor.document.features.iter().find(|feature| {
                    Some(feature.id) == t.feature
                        && (t.view != View::Primers || feature.kind == "primer_bind")
                });
                let changed = match (&u.feature_details_saved, selected) {
                    (Some((previous, molecule)), Some(feature)) => {
                        previous != feature || *molecule != t.editor.document.molecule
                    }
                    (None, None) => false,
                    _ => true,
                };
                if changed {
                    let value = selected.map_or_else(String::new, |feature| {
                        feature_card::detail_text(feature, t.editor.document.molecule)
                    });
                    u.feature_details
                        .update(cx, |state, cx| state.set_value(value, window, cx));
                    u.feature_details_saved = selected
                        .cloned()
                        .map(|feature| (feature, t.editor.document.molecule));
                }
            }
            if t.view == View::Properties
                && t.editor.document.molecule == Molecule::Protein
                && u.protein_properties.as_ref().is_none_or(|(sequence, _)| {
                    !std::sync::Arc::ptr_eq(sequence, &t.editor.document.sequence)
                })
            {
                u.protein_properties = Some((
                    t.editor.document.sequence.clone(),
                    popgene_core::protein::analyze(&t.editor.document.sequence),
                ));
            }
            if matches!(t.view, View::History | View::Parameters) && !u.source_metadata_loaded {
                u.source_metadata_loaded = true;
                if let Some(source) = &t.original_source {
                    match popgene_formats::source_metadata(source) {
                        Ok(fields) if !fields.is_empty() => {
                            let value = fields
                                .into_iter()
                                .map(|(label, value)| format!("{label}\n{value}"))
                                .collect::<Vec<_>>()
                                .join("\n\n");
                            u.source_metadata = Some(
                                cx.new(|cx| TextareaState::new(window, cx).default_value(value)),
                            );
                        }
                        Ok(_) => {}
                        Err(error) => u.source_metadata_error = Some(error.to_string()),
                    }
                }
            }
            if t.view.has_map() && !t.draft {
                u.map.update(cx, |m, cx| {
                    m.set_document(
                        t.editor.document.clone(),
                        t.selected.clone(),
                        MapOptions {
                            linear: t.map_linear,
                            labels: t.show_labels,
                            features: t.show_features,
                            primers: t.show_primers,
                            positions: t.show_positions,
                            sites: t.show_sites,
                            zoom: t.map_zoom,
                        },
                        map_sites(t),
                        cx,
                    )
                });
            }
            if t.view.has_sequence() && !t.draft {
                u.sequence.update(cx, |s, cx| {
                    s.set_document(
                        t.editor.document.clone(),
                        t.selected.clone(),
                        t.scroll_to.take(),
                        u.sequence_options,
                        sequence_sites(t, u.sequence_options.show_enzymes),
                        cx,
                    )
                });
            }
            if t.view == View::Trace && !t.draft {
                u.trace.update(cx, |s, cx| {
                    s.set_document(t.editor.document.clone(), t.trace_scale, cx)
                });
            }
        }
    }
    fn edit(&mut self, window: &mut Window, cx: &mut Context<Self>) {
        if !self.mutable() {
            return;
        }
        let t = self.model.current().unwrap();
        if t.editor.document.trace.is_some() {
            return;
        }
        if t.selected.len() != 1 && !t.editor.document.placeholder {
            self.model
                .current_mut()
                .unwrap()
                .message("Select one contiguous range to edit.", true);
            cx.notify();
            return;
        }
        let id = t.id;
        let form =
            cx.new(|cx| EditForm::new(t.editor.document.clone(), t.selected.clone(), window, cx));
        cx.subscribe_in(&form, window, move |this, form, event, window, cx| {
            match event {
                EditEvent::Cancel => this.modal = None,
                EditEvent::Apply { start, end, bases } => {
                    if let Some(t) = this.model.tabs.iter_mut().find(|t| t.id == id) {
                        if t.apply(|doc| doc.replace(*start, *end, bases)) {
                            this.modal = None;
                            this.focus.focus(window, cx);
                        } else {
                            form.update(cx, |f, cx| f.error(t.status.clone(), cx));
                        }
                    }
                }
            }
            cx.notify();
        })
        .detach();
        self.modal = Some(Modal::Edit(form));
        cx.notify();
    }
    fn annotate(&mut self, window: &mut Window, cx: &mut Context<Self>) {
        if !self.mutable() {
            return;
        }
        let t = self.model.current().unwrap();
        if t.selected.len() != 1 || t.editor.document.trace.is_some() {
            return;
        }
        let part = t.selected[0].clone();
        let id = t.id;
        let form = cx.new(|cx| {
            AnnotationForm::new(t.editor.document.clone(), t.selected.clone(), window, cx)
        });
        cx.subscribe(&form, move |this, form, event, cx| {
            match event {
                AnnotationEvent::Cancel => this.modal = None,
                AnnotationEvent::Apply { name, reverse } => {
                    if let Some(t) = this.model.tabs.iter_mut().find(|t| t.id == id) {
                        if t.apply(|d| d.annotate(name.clone(), part.clone(), *reverse)) {
                            this.modal = None;
                        } else {
                            form.update(cx, |f, cx| f.error(t.status.clone(), cx));
                        }
                    }
                }
            }
            cx.notify();
        })
        .detach();
        self.modal = Some(Modal::Annotation(form));
        cx.notify();
    }
    fn export(&mut self, _window: &mut Window, cx: &mut Context<Self>) {
        if !self.mutable() {
            return;
        }
        let i = self.model.active;
        let format = SaveFormat::recommended(&self.model.tabs[i].editor.document);
        self.format_dialog(i, format, true, false, None, cx);
    }
    fn history_controls(&self, cx: &mut Context<Self>) -> AnyElement {
        let t = self.model.current();
        let editable = self.mutable();
        let undo = editable && t.is_some_and(|t| t.editor.can_undo());
        let redo = editable && t.is_some_and(|t| t.editor.can_redo());
        div()
            .h_flex()
            .gap_1()
            .child(
                Button::new("undo")
                    .accessibility_label("Undo")
                    .ghost()
                    .small()
                    .icon(IconName::Undo)
                    .tooltip("Undo document edit · ⌘Z")
                    .disabled(!undo)
                    .on_click(cx.listener(|s, _, _, cx| {
                        if let Some(t) = s.model.current_mut() {
                            t.undo();
                        }
                        cx.notify();
                    })),
            )
            .child(
                Button::new("redo")
                    .accessibility_label("Redo")
                    .ghost()
                    .small()
                    .icon(IconName::Redo)
                    .tooltip("Redo · ⇧⌘Z")
                    .disabled(!redo)
                    .on_click(cx.listener(|s, _, _, cx| {
                        if let Some(t) = s.model.current_mut() {
                            t.redo();
                        }
                        cx.notify();
                    })),
            )
            .into_any_element()
    }
    fn open_button(&self, cx: &mut Context<Self>) -> AnyElement {
        let button = Button::new("open-files")
            .accessibility_label("Open file or project folder")
            .ghost()
            .small()
            .icon(IconName::FolderOpen)
            .label("Open")
            .tooltip("Open files or a project folder · ⌘O");
        #[cfg(target_os = "macos")]
        {
            button
                .on_click(cx.listener(|s, _, _, cx| s.open(cx)))
                .into_any_element()
        }
        #[cfg(not(target_os = "macos"))]
        {
            let owner = cx.entity().downgrade();
            button
                .dropdown_menu(move |menu, _, _| {
                    let files = owner.clone();
                    let folders = owner.clone();
                    menu.item(PopupMenuItem::new("Open Files…").on_click(move |_, _, cx| {
                        let _ = files.update(cx, |s, cx| s.open(cx));
                    }))
                    .item(
                        PopupMenuItem::new("Open Project Folder…").on_click(move |_, _, cx| {
                            let _ = folders.update(cx, |s, cx| s.folder(cx));
                        }),
                    )
                })
                .into_any_element()
        }
    }
    fn toolbar(&self, cx: &mut Context<Self>) -> AnyElement {
        let t = self.model.current();
        let editable = self.mutable();
        let nucleotide = editable
            && t.is_some_and(|t| {
                t.editor.document.molecule != Molecule::Protein
                    && t.editor.document.trace.is_none()
                    && !t.editor.document.placeholder
            });
        div()
            .h_flex()
            .h(px(52.))
            .px_4()
            .gap_2()
            .bg(rgb(PANEL))
            .border_b_1()
            .border_color(rgb(BORDER))
            .child(self.open_button(cx))
            .child(
                Button::new("save")
                    .primary()
                    .small()
                    .icon(gpui_kit::assets::IconName::Save)
                    .label("Save")
                    .disabled(!editable)
                    .on_click(cx.listener(|s, _, _, cx| s.save(false, cx))),
            )
            .child(
                Button::new("save-as")
                    .ghost()
                    .small()
                    .label("Save as…")
                    .disabled(!editable)
                    .on_click(cx.listener(|s, _, _, cx| s.save(true, cx))),
            )
            .child(div().w(px(1.)).h(px(22.)).mx_2().bg(rgb(BORDER)))
            .child(self.history_controls(cx))
            .child(
                Button::new("reverse")
                    .ghost()
                    .small()
                    .label("Reverse complement")
                    .disabled(!nucleotide)
                    .on_click(cx.listener(|s, _, _, cx| {
                        if let Some(t) = s.model.current_mut() {
                            t.apply(Document::revcomp);
                        }
                        cx.notify();
                    })),
            )
            .child(
                Button::new("topology")
                    .ghost()
                    .small()
                    .label(if t.is_some_and(|t| t.editor.document.circular) {
                        "Circular"
                    } else {
                        "Linear"
                    })
                    .disabled(!nucleotide)
                    .on_click(cx.listener(|s, _, _, cx| {
                        if let Some(t) = s.model.current_mut() {
                            t.apply(|d| {
                                d.circular = !d.circular;
                                Ok(())
                            });
                        }
                        cx.notify();
                    })),
            )
            .child(div().flex_1())
            .child(
                Button::new("export")
                    .ghost()
                    .small()
                    .icon(gpui_kit::assets::IconName::Download)
                    .label("Export")
                    .disabled(!editable || t.is_some_and(|t| t.editor.document.placeholder))
                    .on_click(cx.listener(|s, _, w, cx| s.export(w, cx))),
            )
            .into_any_element()
    }
    fn tabs(&mut self, cx: &mut Context<Self>) -> AnyElement {
        let active = self.model.current().map(|t| t.id);
        if self.shown_tab != active {
            self.tabs_scroll.scroll_to_item(self.model.active);
            self.shown_tab = active;
        }
        let tabs = self
            .model
            .tabs
            .iter()
            .enumerate()
            .map(|(i, t)| {
                let id = t.id;
                let selected = i == self.model.active;
                let color = match t.editor.document.molecule {
                    Molecule::Dna => ACCENT,
                    Molecule::Rna => 0xa9cf76,
                    Molecule::Protein => 0xa39be9,
                };
                div()
                    .id(("tab", id))
                    .h_flex()
                    .h_full()
                    .min_w(px(174.))
                    .max_w(px(255.))
                    .px_3()
                    .gap_2()
                    .border_r_1()
                    .border_b_2()
                    .border_color(if selected { rgb(ACCENT) } else { rgb(BORDER) })
                    .bg(rgb(if selected { PANEL } else { BG }))
                    .cursor_pointer()
                    .hover(|s| s.bg(rgb(0x223040)))
                    .on_click(cx.listener(move |s, _, w, cx| s.switch(i, w, cx)))
                    .child(div().text_color(rgb(color)).text_sm().child(if t.draft {
                        "＋"
                    } else if t.editor.document.trace.is_some() {
                        "⌁"
                    } else {
                        "◌"
                    }))
                    .child(
                        div()
                            .flex_1()
                            .min_w_0()
                            .text_sm()
                            .truncate()
                            .child(t.title()),
                    )
                    .when(t.dirty(), |d| {
                        d.child(div().text_color(rgb(ACCENT)).child("•"))
                    })
                    .child(
                        Button::new(("close-tab", id))
                            .accessibility_label("Close tab")
                            .ghost()
                            .xsmall()
                            .icon(IconName::Close)
                            .tooltip("Close tab")
                            .on_click(cx.listener(move |s, _, _, cx| {
                                cx.stop_propagation();
                                s.close_tab(id, cx);
                            })),
                    )
            })
            .collect::<Vec<_>>();
        div()
            .h_flex()
            .h(px(42.))
            .bg(rgb(BG))
            .border_b_1()
            .border_color(rgb(BORDER))
            .child(
                div()
                    .id("file-tabs")
                    .h_flex()
                    .h_full()
                    .flex_1()
                    .min_w_0()
                    .overflow_x_scroll()
                    .track_scroll(&self.tabs_scroll)
                    .children(tabs),
            )
            .child(
                Button::new("new-sequence")
                    .accessibility_label("New sequence")
                    .ghost()
                    .small()
                    .icon(IconName::Plus)
                    .tooltip("New sequence · ⌘N")
                    .on_click(cx.listener(|s, _, _, cx| s.new_sequence(cx))),
            )
            .into_any_element()
    }
    fn navigator(&self, cx: &mut Context<Self>) -> AnyElement {
        let query = self.project_filter.read(cx).value().to_lowercase();
        let mut tree = div().v_flex().gap_1();
        if let Some(project) = &self.model.project {
            for entry in &project.entries {
                tree = tree.children(self.tree(entry, 0, &query, cx));
            }
        } else {
            tree = tree.child(div().p_3().text_sm().text_color(rgb(MUTED)).child(
                "Open a folder to browse a SnapGene project or a collection of sequence files.",
            ));
        }
        div()
            .v_flex()
            .size_full()
            .text_left()
            .bg(rgb(PANEL))
            .border_r_1()
            .border_color(rgb(BORDER))
            .child(
                div()
                    .h_flex()
                    .px_4()
                    .pt_5()
                    .pb_3()
                    .gap_2()
                    .child(eyebrow("PROJECT"))
                    .child(div().flex_1())
                    .child(
                        Button::new("project-folder")
                            .accessibility_label("Open project folder")
                            .ghost()
                            .xsmall()
                            .icon(IconName::FolderOpen)
                            .tooltip("Open project folder")
                            .on_click(cx.listener(|s, _, _, cx| s.folder(cx))),
                    )
                    .child(
                        Button::new("refresh")
                            .accessibility_label("Refresh project")
                            .ghost()
                            .xsmall()
                            .icon(IconName::RefreshCw)
                            .tooltip("Refresh project")
                            .disabled(self.model.project.is_none())
                            .on_click(cx.listener(|s, _, _, cx| {
                                if let Some(p) = &s.model.project {
                                    s.model.open_folder(p.root.clone());
                                    s.wake(cx);
                                }
                            })),
                    ),
            )
            .when_some(self.model.project.as_ref(), |d, p| {
                d.child(
                    div()
                        .px_4()
                        .pb_3()
                        .font_weight(FontWeight::SEMIBOLD)
                        .truncate()
                        .child(
                            p.root
                                .file_name()
                                .unwrap_or_default()
                                .to_string_lossy()
                                .into_owned(),
                        ),
                )
            })
            .child(
                div().px_3().pb_3().child(
                    Input::new(&self.project_filter)
                        .small()
                        .prefix(Icon::new(IconName::Search).small())
                        .aria_label("Filter project files"),
                ),
            )
            .child(
                div()
                    .id("project-tree")
                    .min_h_0()
                    .flex_1()
                    .overflow_y_scroll()
                    .px_2()
                    .child(tree),
            )
            .when(self.model.project.is_some(), |d| {
                d.child(
                    div().p_3().child(
                        Button::new("open-project")
                            .small()
                            .w_full()
                            .accessibility_label("Open supported files")
                            .child(navigation_label("Open supported files", None))
                            .on_click(cx.listener(|s, _, _, cx| {
                                if let Some(p) = &s.model.project {
                                    s.model.open_paths(p.supported_files());
                                    s.wake(cx);
                                }
                            })),
                    ),
                )
            })
            .child(
                div()
                    .border_t_1()
                    .border_color(rgb(BORDER))
                    .px_3()
                    .py_3()
                    .flex_shrink_0()
                    .v_flex()
                    .gap_2()
                    .child(eyebrow("OPEN DOCUMENTS"))
                    .child(
                        div()
                            .id("open-documents-scroll")
                            .max_h(px(240.))
                            .overflow_y_scroll()
                            .v_flex()
                            .gap_1()
                            .children(self.model.tabs.iter().enumerate().map(|(i, t)| {
                                Button::new(("open-doc", t.id))
                                    .ghost()
                                    .small()
                                    .w_full()
                                    .flex_shrink_0()
                                    .justify_start()
                                    .selected(i == self.model.active)
                                    .accessibility_label(t.title())
                                    .child(navigation_label(
                                        format!(
                                            "{}{}",
                                            if t.dirty() { "• " } else { "" },
                                            t.title()
                                        ),
                                        None,
                                    ))
                                    .on_click(cx.listener(move |s, _, w, cx| s.switch(i, w, cx)))
                            })),
                    ),
            )
            .child(
                div().p_3().border_t_1().border_color(rgb(BORDER)).child(
                    Button::new("import-report")
                        .ghost()
                        .small()
                        .w_full()
                        .accessibility_label(format!(
                            "Import report · {} issues",
                            self.model.issues.len()
                        ))
                        .child(navigation_label(
                            format!("Import report · {} issues", self.model.issues.len()),
                            None,
                        ))
                        .on_click(cx.listener(|s, _, _, cx| {
                            s.modal = Some(Modal::Report);
                            cx.notify();
                        })),
                ),
            )
            .into_any_element()
    }
    fn tree(
        &self,
        e: &popgene_formats::project::Entry,
        depth: usize,
        query: &str,
        cx: &mut Context<Self>,
    ) -> Vec<AnyElement> {
        fn matches(e: &popgene_formats::project::Entry, q: &str) -> bool {
            q.is_empty()
                || e.name.to_lowercase().contains(q)
                || e.children.iter().any(|e| matches(e, q))
        }
        if !matches(e, query) {
            return vec![];
        }
        let path = e.path.clone();
        let directory = e.directory;
        let supported = e.supported;
        let collapsed = self.collapsed.contains(&path) && query.is_empty();
        let selected = self
            .model
            .current()
            .is_some_and(|t| t.source_path.as_ref() == Some(&path));
        let mut rows = vec![
            div()
                .pl(px(depth as f32 * 12.))
                .child(
                    Button::new(SharedString::from(format!("project:{}", path.display())))
                        .ghost()
                        .small()
                        .w_full()
                        .justify_start()
                        .selected(selected)
                        .accessibility_label(e.name.clone())
                        .child(navigation_label(
                            e.name.clone(),
                            Some(if directory {
                                if collapsed {
                                    IconName::ChevronRight
                                } else {
                                    IconName::ChevronDown
                                }
                            } else {
                                IconName::File
                            }),
                        ))
                        .disabled(!directory && !supported)
                        .tooltip(path.display().to_string())
                        .on_click(cx.listener(move |s, _, _, cx| {
                            if s.blocked() {
                                return;
                            }
                            if directory {
                                if !s.collapsed.remove(&path) {
                                    s.collapsed.insert(path.clone());
                                }
                                cx.notify();
                            } else {
                                s.model.open_paths(vec![path.clone()]);
                                s.wake(cx);
                            }
                        })),
                )
                .into_any_element(),
        ];
        if directory && !collapsed {
            for child in &e.children {
                rows.extend(self.tree(child, depth + 1, query, cx));
            }
        }
        rows
    }
    fn content(&mut self, cx: &mut Context<Self>) -> AnyElement {
        let Some(t) = self.model.current() else {
            return div()
                .v_flex()
                .size_full()
                .items_center()
                .justify_center()
                .gap_4()
                .child(div().text_2xl().child("Your next discovery starts here"))
                .child(
                    div()
                        .text_color(rgb(MUTED))
                        .child("Create a sequence or open an existing file."),
                )
                .child(
                    Button::new("empty-new")
                        .primary()
                        .label("New sequence")
                        .on_click(cx.listener(|s, _, _, cx| s.new_sequence(cx))),
                )
                .child(
                    Button::new("empty-open")
                        .ghost()
                        .label("Open files…")
                        .on_click(cx.listener(|s, _, _, cx| s.open(cx))),
                )
                .into_any_element();
        };
        let u = &self.ui[&t.id];
        if t.draft {
            return div()
                .size_full()
                .children(u.creator.clone())
                .into_any_element();
        }
        let doc = &t.editor.document;
        let views: Vec<_> = [
            (View::Map, "Map"),
            (View::Sequence, "Sequence"),
            (View::Enzymes, "Enzymes"),
            (View::Features, "Features"),
            (View::Primers, "Primers"),
            (View::Trace, "Trace"),
            (View::Parameters, "Parameters"),
            (View::Properties, "Properties"),
            (View::History, "History"),
        ]
        .into_iter()
        .filter(|(v, _)| {
            (!matches!(v, View::Trace | View::Parameters) || doc.trace.is_some())
                && (*v != View::Properties || doc.molecule == Molecule::Protein)
                && (*v != View::Primers || doc.molecule != Molecule::Protein)
                && (*v != View::Enzymes || doc.molecule == Molecule::Dna)
        })
        .collect();
        let selected = views
            .iter()
            .position(|(view, _)| *view == t.view || (*view == View::Map && t.view == View::Split))
            .unwrap_or(0);
        let nav = TabBar::new("document-view-tabs")
            .underline()
            .h(px(46.))
            .items_end()
            .flex_shrink_0()
            .px_4()
            .selected_index(selected)
            .children(views.into_iter().map(|(view, label)| {
                ViewTab::new()
                    .label(label)
                    .on_click(cx.listener(move |s, _, w, cx| {
                        if !s.blocked() {
                            if let Some(t) = s.model.current_mut() {
                                t.view = view;
                            }
                            s.focus.focus(w, cx);
                            cx.notify();
                        }
                    }))
            }));
        let title = div()
            .h_flex()
            .px_6()
            .pt_5()
            .pb_4()
            .gap_3()
            .child(
                div()
                    .v_flex()
                    .gap_1()
                    .min_w_0()
                    .flex_1()
                    .child(
                        div()
                            .text_xl()
                            .font_weight(FontWeight::SEMIBOLD)
                            .truncate()
                            .child(doc.name.clone()),
                    )
                    .child(div().text_xs().text_color(rgb(MUTED)).child(format!(
                        "{}  ·  {} annotations  ·  {}",
                        doc.summary(),
                        doc.features.len(),
                        if t.dirty() {
                            "Unsaved changes"
                        } else {
                            "Up to date"
                        }
                    ))),
            )
            .child(
                div()
                    .px_2()
                    .py_1()
                    .rounded_md()
                    .bg(rgb(0x203a3d))
                    .text_xs()
                    .text_color(rgb(ACCENT))
                    .child(doc.molecule.label()),
            );
        let body = match t.view {
            View::Map => self.map_content(cx),
            View::Split => self.split_content(cx),
            View::Enzymes => self.enzyme_content(cx),
            View::Sequence => self.sequence_content(cx),
            View::Features | View::Primers => self.feature_table(t.view == View::Primers, cx),
            View::Trace => div()
                .v_flex()
                .flex_1()
                .min_h_0()
                .child(
                    div()
                        .h_flex()
                        .px_6()
                        .gap_2()
                        .child(eyebrow("SANGER TRACE"))
                        .child(
                            Button::new("trace-minus")
                                .ghost()
                                .small()
                                .label("−")
                                .on_click(cx.listener(|s, _, _, cx| {
                                    let t = s.model.current_mut().unwrap();
                                    t.trace_scale = (t.trace_scale / 1.25).max(0.5);
                                    cx.notify();
                                })),
                        )
                        .child(format!("{:.1}×", t.trace_scale))
                        .child(
                            Button::new("trace-plus")
                                .ghost()
                                .small()
                                .label("+")
                                .on_click(cx.listener(|s, _, _, cx| {
                                    let t = s.model.current_mut().unwrap();
                                    t.trace_scale = (t.trace_scale * 1.25).min(5.);
                                    cx.notify();
                                })),
                        )
                        .child(
                            div()
                                .text_xs()
                                .text_color(rgb(MUTED))
                                .child("Imported calls and quality · read only"),
                        ),
                )
                .child(div().flex_1().min_h_0().child(u.trace.clone()))
                .into_any_element(),
            View::Properties => self.ui[&t.id].protein_properties.as_ref().map_or_else(
                || {
                    div()
                        .child("Calculating protein properties…")
                        .into_any_element()
                },
                |(_, properties)| crate::protein_properties::render(properties),
            ),
            View::Parameters => self.trace_parameters(cx),
            View::History => self.history(cx),
        };
        div()
            .v_flex()
            .size_full()
            .bg(rgb(BG))
            .child(nav)
            .when(!t.view.has_map() && t.view != View::Sequence, |d| {
                d.child(title)
            })
            .child(body)
            .into_any_element()
    }
    fn split_button(&self, cx: &mut Context<Self>) -> Button {
        let t = self.model.current().unwrap();
        Button::new("map-sequence-split")
            .ghost()
            .small()
            .icon(AssetIcon::Rows2)
            .tooltip(if t.view == View::Split {
                "Show map only"
            } else {
                "Split view · map above sequence"
            })
            .accessibility_label("Split map and sequence view")
            .selected(t.view == View::Split)
            .toggled(t.view == View::Split)
            .disabled(t.editor.document.placeholder)
            .on_click(cx.listener(|s, _, w, cx| {
                if !s.blocked() {
                    if let Some(t) = s.model.current_mut() {
                        t.view = if t.view == View::Split {
                            View::Map
                        } else {
                            View::Split
                        };
                    }
                    s.focus.focus(w, cx);
                    cx.notify();
                }
            }))
    }
    fn split_content(&self, cx: &mut Context<Self>) -> AnyElement {
        let t = self.model.current().unwrap();
        let u = &self.ui[&t.id];
        div()
            .flex_1()
            .min_h_0()
            .child(
                v_resizable(("map-sequence-panes", t.id))
                    .with_state(&u.split)
                    .child(
                        resizable_panel()
                            .size(px(420.))
                            .size_range(px(280.)..Pixels::MAX)
                            .child(
                                div()
                                    .v_flex()
                                    .size_full()
                                    .overflow_hidden()
                                    .child(self.map_content(cx)),
                            ),
                    )
                    .child(
                        resizable_panel().size_range(px(160.)..Pixels::MAX).child(
                            div()
                                .v_flex()
                                .size_full()
                                .overflow_hidden()
                                .child(self.sequence_content(cx)),
                        ),
                    ),
            )
            .into_any_element()
    }
    fn cutter_filter(&self, control_id: &'static str, cx: &mut Context<Self>) -> AnyElement {
        let t = self.model.current().unwrap();
        let id = t.id;
        let selected = t.enzyme_filter;
        let owner = cx.entity().downgrade();
        Button::new(control_id)
            .ghost()
            .small()
            .when(control_id == "enzyme-cutter-filter", |b| {
                b.label(format!("{} ▾", selected.label()))
            })
            .when(control_id != "enzyme-cutter-filter", |b| {
                b.icon(AssetIcon::ListFilter)
            })
            .tooltip(format!("Restriction enzymes · {}", selected.label()))
            .accessibility_label("Restriction enzyme filter")
            .disabled(!t.restriction_eligible())
            .dropdown_menu(move |mut menu, _, _| {
                menu = menu.label(format!("Included catalogue · {} enzymes", ENZYMES.len()));
                for filter in EnzymeFilter::ALL {
                    let owner = owner.clone();
                    menu = menu.item(
                        PopupMenuItem::new(filter.label())
                            .checked(filter == selected)
                            .on_click(move |_, _, cx| {
                                let _ = owner.update(cx, |s, cx| {
                                    if s.blocked() {
                                        return;
                                    }
                                    if let Some(t) = s.model.tabs.iter_mut().find(|t| t.id == id) {
                                        t.enzyme_filter = filter;
                                    }
                                    cx.notify();
                                });
                            }),
                    );
                }
                menu
            })
            .into_any_element()
    }
    fn display_menu(&self, cx: &mut Context<Self>) -> AnyElement {
        let t = self.model.current().unwrap();
        let id = t.id;
        let settings = [
            ("Restriction sites", t.show_sites),
            ("Features", t.show_features),
            ("Primers", t.show_primers),
            ("Ruler", t.show_positions),
            ("Labels", t.show_labels),
            ("Group equivalent cut sites", t.nonredundant),
        ];
        let owner = cx.entity().downgrade();
        Button::new("map-display")
            .ghost()
            .small()
            .label("Display ▾")
            .accessibility_label("Map display options")
            .dropdown_menu(move |mut menu, _, _| {
                menu = menu.label("Map layers");
                for (i, (label, checked)) in settings.iter().copied().enumerate() {
                    let owner = owner.clone();
                    menu = menu.item(PopupMenuItem::new(label).checked(checked).on_click(
                        move |_, _, cx| {
                            let _ = owner.update(cx, |s, cx| {
                                if s.blocked() {
                                    return;
                                }
                                if let Some(t) = s.model.tabs.iter_mut().find(|t| t.id == id) {
                                    match i {
                                        0 => t.show_sites = !t.show_sites,
                                        1 => t.show_features = !t.show_features,
                                        2 => t.show_primers = !t.show_primers,
                                        3 => t.show_positions = !t.show_positions,
                                        4 => t.show_labels = !t.show_labels,
                                        _ => t.nonredundant = !t.nonredundant,
                                    }
                                }
                                cx.notify();
                            });
                        },
                    ));
                }
                menu
            })
            .into_any_element()
    }
    fn map_content(&self, cx: &mut Context<Self>) -> AnyElement {
        let t = self.model.current().unwrap();
        let u = &self.ui[&t.id];
        let doc = &t.editor.document;
        let layers = [
            ("Restriction sites", AssetIcon::Scissors, t.show_sites),
            ("Features", AssetIcon::ArrowBigRight, t.show_features),
            ("Primers", AssetIcon::ArrowRightLeft, t.show_primers),
            ("Ruler", AssetIcon::Ruler, t.show_positions),
            ("Labels", AssetIcon::Tags, t.show_labels),
        ];
        div()
            .id("map-viewer")
            .h_flex()
            .items_stretch()
            .flex_1()
            .min_h_0()
            .min_w_0()
            .overflow_hidden()
            .child(
                viewer_rail("map-layers", false)
                    .when(doc.circular, |d| {
                        d.child(
                            Button::new("circle")
                                .ghost()
                                .small()
                                .icon(AssetIcon::Circle)
                                .tooltip("Circular map")
                                .accessibility_label("Circular map")
                                .selected(!t.map_linear)
                                .toggled(!t.map_linear)
                                .on_click(cx.listener(|s, _, _, cx| {
                                    s.model.current_mut().unwrap().map_linear = false;
                                    cx.notify();
                                })),
                        )
                        .child(
                            Button::new("linear")
                                .ghost()
                                .small()
                                .icon(AssetIcon::MoveHorizontal)
                                .tooltip("Linear map")
                                .accessibility_label("Linear map")
                                .selected(t.map_linear)
                                .toggled(t.map_linear)
                                .on_click(cx.listener(|s, _, _, cx| {
                                    s.model.current_mut().unwrap().map_linear = true;
                                    cx.notify();
                                })),
                        )
                    })
                    .child(self.split_button(cx))
                    .children(layers.into_iter().enumerate().map(
                        |(i, (label, icon, selected))| {
                            let button = Button::new(("map-layer", i))
                                .ghost()
                                .small()
                                .icon(icon)
                                .tooltip(if i == 0 && !t.restriction_eligible() {
                                    "Restriction sites require double-stranded DNA".to_owned()
                                } else {
                                    format!(
                                        "{} {}",
                                        if selected { "Hide" } else { "Show" },
                                        label.to_lowercase()
                                    )
                                })
                                .accessibility_label(label)
                                .selected(selected)
                                .toggled(selected)
                                .disabled(i == 0 && !t.restriction_eligible())
                                .on_click(cx.listener(move |s, _, _, cx| {
                                    let t = s.model.current_mut().unwrap();
                                    match i {
                                        0 => t.show_sites = !t.show_sites,
                                        1 => t.show_features = !t.show_features,
                                        2 => t.show_primers = !t.show_primers,
                                        3 => t.show_positions = !t.show_positions,
                                        _ => t.show_labels = !t.show_labels,
                                    }
                                    cx.notify();
                                }));
                            div()
                                .v_flex()
                                .gap_1()
                                .child(button)
                                .when(i == 0, |d| d.child(self.cutter_filter("cutter-filter", cx)))
                        },
                    ))
            )
            .child(div().flex_1().min_w_0().min_h_0().child(u.map.clone()))
            .child(
                viewer_rail("map-tools", true)
                    .child(eyebrow("VIEW"))
                    .child(div().h_flex().justify_center().gap_1()
                    .child(
                        Button::new("map-zoom-out")
                            .ghost()
                            .small()
                            .label("−")
                            .accessibility_label("Zoom out map")
                            .disabled(t.map_zoom <= 1.)
                            .on_click(cx.listener(|s, _, _, cx| {
                                let t = s.model.current_mut().unwrap();
                                t.map_zoom = (t.map_zoom / 1.25).max(1.);
                                cx.notify();
                            })),
                    )
                    .child(
                        div()
                            .text_xs()
                            .w(px(40.))
                            .text_center()
                            .child(format!("{:.0}%", t.map_zoom * 100.)),
                    )
                    .child(
                        Button::new("map-zoom-in")
                            .ghost()
                            .small()
                            .label("+")
                            .accessibility_label("Zoom in map")
                            .disabled(t.map_zoom >= 3.)
                            .on_click(cx.listener(|s, _, _, cx| {
                                let t = s.model.current_mut().unwrap();
                                t.map_zoom = (t.map_zoom * 1.25).min(3.);
                                cx.notify();
                            })),
                    )
                    )
                    .child(
                        Button::new("map-fit")
                            .ghost()
                            .small()
                            .label("Fit")
                            .accessibility_label("Fit map to window")
                            .on_click(cx.listener(|s, _, _, cx| {
                                s.model.current_mut().unwrap().map_zoom = 1.;
                                cx.notify();
                            })),
                    )
                    .child(
                        Button::new("map-find")
                            .ghost()
                            .small()
                            .label("Find")
                            .selected(t.map_find && t.view != View::Split)
                            .on_click(cx.listener(|s, _, w, cx| {
                                let t = s.model.current_mut().unwrap();
                                t.map_find = true;
                                s.ui[&t.id].search.read(cx).focus_handle(cx).focus(w, cx);
                                cx.notify();
                            })),
                    )
                    .child(self.display_menu(cx))

                    .when(t.view != View::Split, |d| d.child(
                    div()
                        .v_flex()
                        .w_full()
                        .pt_2()
                        .gap_2()
                        .child(
                            div()
                                .w_full()
                                .child(Input::new(&u.search).small().aria_label("Find DNA on map")),
                        )
                        .child(
                            Button::new("map-previous-match")
                                .ghost()
                                .small()
                                .label("Previous")
                                .disabled(t.hits.is_empty())
                                .on_click(cx.listener(|s, _, _, cx| {
                                    s.model.current_mut().unwrap().previous_hit();
                                    cx.notify();
                                })),
                        )
                        .child(
                            Button::new("map-next-match")
                                .ghost()
                                .small()
                                .label("Next match")
                                .disabled(t.hits.is_empty())
                                .on_click(cx.listener(|s, _, _, cx| {
                                    s.model.current_mut().unwrap().next_hit();
                                    cx.notify();
                                })),
                        )
                        .child(div().text_xs().text_color(rgb(MUTED)).child(
                            if t.search_error.is_empty() {
                                format!(
                                    "{} matches{}",
                                    t.hits.len(),
                                    if t.hits.len() == 1000 {
                                        " (capped)"
                                    } else {
                                        ""
                                    }
                                )
                            } else {
                                t.search_error.clone()
                            },
                        )),
                    ))
            .when(t.show_sites && t.restriction_eligible() && (t.analyzing() || t.restriction_error.is_some() || t.restriction_ambiguous || t.restriction.as_ref().is_some_and(|a| a.truncated)), |d| {
                d.child(div().py_1().text_xs().text_color(rgb(0xe5b879)).child(
                    if let Some(error) = &t.restriction_error { error.clone() }
                    else if t.analyzing() { "Finding restriction sites…".into() }
                    else if t.restriction_ambiguous { "Ambiguous-base matches omitted; restriction counts may be incomplete.".into() }
                    else { "Restriction positions are capped; inspect an enzyme for retained and total counts.".into() }
                ))
            })

            )
            .into_any_element()
    }

    fn enzyme_content(&self, cx: &mut Context<Self>) -> AnyElement {
        let t = self.model.current().unwrap();
        let tab_id = t.id;
        let Some(analysis) = &t.restriction else {
            return div()
                .p_6()
                .v_flex()
                .gap_3()
                .child(if !t.restriction_eligible() {
                    "Restriction analysis requires a nonempty double-stranded DNA sequence."
                        .to_owned()
                } else if let Some(error) = &t.restriction_error {
                    error.clone()
                } else {
                    "Finding restriction sites…".to_owned()
                })
                .into_any_element();
        };
        let mut rows = Vec::new();
        for (i, e) in ENZYMES.iter().enumerate() {
            let count = analysis.counts[i];
            if !t.enzyme_filter.includes(count, e.recognition.len()) {
                continue;
            }
            let sites: Vec<_> = analysis
                .sites
                .iter()
                .filter(|s| s.enzyme_index == i)
                .collect();
            let retained = sites.len();
            let mut row = div()
                .v_flex()
                .py_3()
                .px_4()
                .gap_2()
                .border_b_1()
                .border_color(rgb(BORDER))
                .child(
                    div()
                        .h_flex()
                        .flex_wrap()
                        .gap_3()
                        .child(
                            Button::new(("enzyme-details", i))
                                .ghost()
                                .small()
                                .w(px(110.))
                                .font_weight(FontWeight::SEMIBOLD)
                                .accessibility_label(format!("Inspect {}", e.name))
                                .child(navigation_label(e.name, None))
                                .on_click(cx.listener(move |s, _, _, cx| {
                                    s.inspect_enzyme(tab_id, i, None, cx)
                                })),
                        )
                        .child(
                            div()
                                .w(px(130.))
                                .font_family("monospace")
                                .text_color(rgb(ACCENT))
                                .child(e.recognition),
                        )
                        .child(div().text_xs().text_color(rgb(MUTED)).child(format!(
                            "{count} recognition site{}",
                            if count == 1 { "" } else { "s" }
                        ))),
                );
            if e.name == "NmeAIII" {
                row =
                    row.child(div().text_xs().text_color(rgb(MUTED)).child(
                        "Nominal cut: NmeAIII can shift by one base with sequence context.",
                    ));
            }
            let mut buttons = Vec::new();
            for (j, site) in sites.iter().take(64).enumerate() {
                let parts = t.site_parts(site);
                let start = site.start;
                let label = match (site.top_cut, site.bottom_cut) {
                    (Some(top), Some(bottom)) => format!(
                        "{} / {}{}",
                        cut_label(
                            top,
                            t.editor.document.sequence.len(),
                            t.editor.document.circular
                        ),
                        cut_label(
                            bottom,
                            t.editor.document.sequence.len(),
                            t.editor.document.circular
                        ),
                        if site.reverse { " ←" } else { " →" }
                    ),
                    _ => format!("{} · cut outside sequence", site.start + 1),
                };
                buttons.push(
                    Button::new((SharedString::from(format!("enzyme-{i}")), j))
                        .ghost()
                        .small()
                        .accessibility_label(format!(
                            "{} cut after top / bottom base {label}",
                            e.name
                        ))
                        .label(label)
                        .on_click(cx.listener(move |s, _, _, cx| {
                            let t = s.model.current_mut().unwrap();
                            t.select(parts.clone(), None, true);
                            t.selected_enzyme = Some(i);
                            s.inspect_enzyme(tab_id, i, Some(start), cx);
                            cx.notify();
                        })),
                );
            }
            row = row.child(div().h_flex().flex_wrap().gap_1().children(buttons));
            if count > 64 {
                row = row.child(div().text_xs().text_color(rgb(MUTED)).child(format!(
                    "First {} positions shown; {retained} retained; {count} total sites",
                    retained.min(64)
                )));
            }
            rows.push(row);
        }
        let empty = rows.is_empty();
        div().v_flex().flex_1().min_h_0()
            .child(div().h_flex().px_6().pb_3().gap_3().child(self.cutter_filter("enzyme-cutter-filter", cx))
                .child(div().text_xs().text_color(rgb(MUTED)).child(format!("{} enzymes · sequence matches only",ENZYMES.len()))))
            .child(div().px_6().pb_3().text_xs().text_color(rgb(MUTED)).child("Positions: cut after top / bottom base. Click an enzyme or position to open its details sidebar. Methylation and reaction conditions are not modeled."))
            .when(t.restriction_ambiguous,|d|d.child(div().px_6().pb_3().text_xs().text_color(rgb(0xe5b879)).child("Ambiguous-base recognition windows omitted; cutter counts may be incomplete.")))
            .child(div().id("enzyme-table").flex_1().min_h_0().overflow_y_scroll().px_4()
                .when(empty,|d|d.child(div().p_4().child("No enzymes match this filter."))).children(rows))
            .into_any_element()
    }
    fn toggle_sequence_layer(&mut self, id: u64, layer: usize, cx: &mut Context<Self>) {
        if self.blocked() {
            return;
        }
        if let Some(ui) = self.ui.get_mut(&id) {
            let options = &mut ui.sequence_options;
            match layer {
                0 => options.show_complement = !options.show_complement,
                1 => options.show_features = !options.show_features,
                2 => options.show_primers = !options.show_primers,
                3 => options.show_enzymes = !options.show_enzymes,
                4 => options.show_translations = !options.show_translations,
                5 => options.colored_bases = !options.colored_bases,
                6 => options.compact = !options.compact,
                7..=12 => options.reading_frames ^= 1 << (layer - 7),
                13 => options.reading_frames = if options.reading_frames == 63 { 0 } else { 63 },
                _ => {}
            }
            cx.notify();
        }
    }
    fn sequence_display_menu(&self, cx: &mut Context<Self>) -> AnyElement {
        let tab = self.model.current().unwrap();
        let id = tab.id;
        let options = self.ui[&id].sequence_options;
        let nucleotide = tab.editor.document.molecule != Molecule::Protein;
        let duplex = nucleotide && !tab.editor.document.single_stranded;
        let settings = [
            (0, "Complementary strand", options.show_complement, duplex),
            (1, "Features", options.show_features, true),
            (2, "Primers", options.show_primers, nucleotide),
            (
                3,
                "Restriction enzymes",
                options.show_enzymes,
                tab.restriction_eligible(),
            ),
            (
                4,
                "Calculated CDS translations",
                options.show_translations,
                nucleotide,
            ),
            (5, "Colour bases", options.colored_bases, true),
            (6, "Compact spacing", options.compact, true),
            (
                13,
                "All six reading frames",
                options.reading_frames == 63,
                nucleotide,
            ),
            (7, "Frame +1", options.reading_frames & 1 != 0, nucleotide),
            (8, "Frame +2", options.reading_frames & 2 != 0, nucleotide),
            (9, "Frame +3", options.reading_frames & 4 != 0, nucleotide),
            (10, "Frame −1", options.reading_frames & 8 != 0, nucleotide),
            (11, "Frame −2", options.reading_frames & 16 != 0, nucleotide),
            (12, "Frame −3", options.reading_frames & 32 != 0, nucleotide),
        ];
        let owner = cx.entity().downgrade();
        Button::new("sequence-display")
            .ghost()
            .small()
            .label("Display ▾")
            .accessibility_label("Sequence display options")
            .dropdown_menu(move |mut menu, _, _| {
                menu = menu.label("Sequence display");
                for (layer, label, checked, available) in settings {
                    if !available {
                        continue;
                    }
                    let owner = owner.clone();
                    menu = menu.item(PopupMenuItem::new(label).checked(checked).on_click(
                        move |_, _, cx| {
                            let _ =
                                owner.update(cx, |s, cx| s.toggle_sequence_layer(id, layer, cx));
                        },
                    ));
                }
                menu
            })
            .into_any_element()
    }
    fn sequence_content(&self, cx: &mut Context<Self>) -> AnyElement {
        let t = self.model.current().unwrap();
        let u = &self.ui[&t.id];
        if t.editor.document.placeholder {
            return div()
                .v_flex()
                .flex_1()
                .items_center()
                .justify_center()
                .gap_3()
                .child(div().text_xl().child("Empty placeholder"))
                .child(
                    div()
                        .text_color(rgb(MUTED))
                        .child("Add a sequence when you are ready."),
                )
                .child(
                    Button::new("add-sequence")
                        .primary()
                        .label("Add sequence…")
                        .on_click(cx.listener(|s, _, w, cx| s.edit(w, cx))),
                )
                .into_any_element();
        }
        let id = t.id;
        let options = u.sequence_options;
        let nucleotide = t.editor.document.molecule != Molecule::Protein;
        let layers = [
            (
                1,
                "Features",
                AssetIcon::ArrowBigRight,
                options.show_features,
                true,
            ),
            (
                2,
                "Primers",
                AssetIcon::ArrowRightLeft,
                options.show_primers,
                nucleotide,
            ),
            (
                3,
                "Restriction enzymes",
                AssetIcon::Scissors,
                options.show_enzymes,
                t.restriction_eligible(),
            ),
        ];
        let search_status = if t.query.is_empty() {
            if t.selection_summary.is_empty() {
                format!(
                    "{} {} · Click a base; Shift-click to extend",
                    number(t.editor.document.sequence.len()),
                    t.editor.document.unit()
                )
            } else {
                t.selection_summary.clone()
            }
        } else if !t.search_error.is_empty() {
            t.search_error.clone()
        } else if t.hits.is_empty() {
            "No matches".into()
        } else {
            format!(
                "{} of {} matches{}",
                t.hit_index + 1,
                t.hits.len(),
                if t.hits.len() == 1000 {
                    " (capped)"
                } else {
                    ""
                }
            )
        };
        let analysis_notice = if !options.show_enzymes || !t.restriction_eligible() {
            None
        } else if let Some(error) = &t.restriction_error {
            Some(error.clone())
        } else if t.analyzing() {
            Some("Finding restriction sites…".into())
        } else if t.restriction_ambiguous {
            Some("Ambiguous-base matches omitted; restriction counts may be incomplete.".into())
        } else if t
            .restriction
            .as_ref()
            .is_some_and(|analysis| analysis.truncated)
        {
            Some(
                "Restriction positions are capped; enzyme details show retained and total counts."
                    .into(),
            )
        } else {
            None
        };
        div()
            .id("sequence-viewer")
            .h_flex()
            .items_stretch()
            .flex_1()
            .min_h_0()
            .min_w_0()
            .overflow_hidden()
            .child(
                viewer_rail("sequence-layers", false)
                    .child(self.split_button(cx))
                    .children(
                        layers
                            .into_iter()
                            .filter(|(_, _, _, _, available)| *available)
                            .map(|(layer, label, icon, selected, _)| {
                                Button::new(("sequence-layer", layer))
                                    .ghost()
                                    .small()
                                    .icon(icon)
                                    .accessibility_label(format!("Sequence {label}"))
                                    .tooltip(format!(
                                        "{} {}",
                                        if selected { "Hide" } else { "Show" },
                                        label.to_lowercase()
                                    ))
                                    .selected(selected)
                                    .toggled(selected)
                                    .on_click(cx.listener(move |s, _, _, cx| {
                                        s.toggle_sequence_layer(id, layer, cx)
                                    }))
                            }),
                    )
                    .when(t.restriction_eligible(), |d| {
                        d.child(self.cutter_filter("sequence-cutter-filter", cx))
                    }),
            )
            .child(
                div()
                    .v_flex()
                    .flex_1()
                    .min_w_0()
                    .min_h_0()
                    .child(div().flex_1().min_h_0().child(u.sequence.clone()))
                    .child(
                        div()
                            .h_flex()
                            .h(px(30.))
                            .flex_shrink_0()
                            .px_3()
                            .py_1()
                            .gap_2()
                            .border_t_1()
                            .border_color(rgb(BORDER))
                            .child(
                                div()
                                    .min_w_0()
                                    .flex_1()
                                    .text_xs()
                                    .text_color(rgb(MUTED))
                                    .truncate()
                                    .child(search_status),
                            ),
                    ),
            )
            .child(
                viewer_rail("sequence-tools", true)
                    .child(self.sequence_display_menu(cx))
                    .child(eyebrow("FROM"))
                    .child(
                        div()
                            .w_full()
                            .child(Input::new(&u.start).small().aria_label("Range start")),
                    )
                    .child(eyebrow("TO"))
                    .child(
                        div()
                            .w_full()
                            .child(Input::new(&u.end).small().aria_label("Range end")),
                    )
                    .child(
                        Button::new("select-range")
                            .ghost()
                            .small()
                            .label("Select")
                            .on_click(cx.listener(|s, _, _, cx| {
                                if s.blocked() {
                                    return;
                                }
                                let id = s.model.current().unwrap().id;
                                let u = &s.ui[&id];
                                let start = u.start.read(cx).value().to_string();
                                let end = u.end.read(cx).value().to_string();
                                s.model.current_mut().unwrap().select_range(&start, &end);
                                cx.notify();
                            })),
                    )
                    .child(
                        Button::new("edit-sequence")
                            .ghost()
                            .small()
                            .label("Edit…")
                            .accessibility_label("Edit sequence")
                            .disabled(t.editor.document.trace.is_some() || t.busy())
                            .on_click(cx.listener(|s, _, w, cx| s.edit(w, cx))),
                    )
                    .child(
                        Button::new("annotate-selection")
                            .ghost()
                            .small()
                            .label("Annotate…")
                            .disabled(
                                t.selected.len() != 1
                                    || t.editor.document.trace.is_some()
                                    || t.busy(),
                            )
                            .on_click(cx.listener(|s, _, w, cx| s.annotate(w, cx))),
                    )
                    .child(
                        Button::new("copy-selection")
                            .ghost()
                            .small()
                            .icon(IconName::Copy)
                            .label("Copy")
                            .disabled(t.selected.is_empty())
                            .on_click(cx.listener(|s, _, _, cx| {
                                if !s.blocked()
                                    && let Some(t) = s.model.current()
                                {
                                    cx.write_to_clipboard(ClipboardItem::new_string(
                                        t.selected_sequence(),
                                    ));
                                }
                            })),
                    )
                    .child(
                        div()
                            .v_flex()
                            .w_full()
                            .gap_2()
                            .pt_2()
                            .child(eyebrow("FIND"))
                            .child(
                                div().w_full().child(
                                    Input::new(&u.search)
                                        .small()
                                        .prefix(Icon::new(IconName::Search).small())
                                        .aria_label("Find exact sequence"),
                                ),
                            )
                            .child(
                                Button::new("previous-hit")
                                    .ghost()
                                    .small()
                                    .label("Previous")
                                    .accessibility_label("Previous sequence match")
                                    .disabled(t.hits.is_empty())
                                    .on_click(cx.listener(|s, _, _, cx| {
                                        if s.blocked() {
                                            return;
                                        }
                                        s.model.current_mut().unwrap().previous_hit();
                                        cx.notify();
                                    })),
                            )
                            .child(
                                Button::new("next-hit")
                                    .ghost()
                                    .small()
                                    .label("Next")
                                    .accessibility_label("Next sequence match")
                                    .disabled(t.hits.is_empty())
                                    .on_click(cx.listener(|s, _, _, cx| {
                                        if s.blocked() {
                                            return;
                                        }
                                        s.model.current_mut().unwrap().next_hit();
                                        cx.notify();
                                    })),
                            ),
                    )
                    .when_some(analysis_notice, |d, notice| {
                        d.child(
                            div()
                                .pb_1()
                                .text_xs()
                                .text_color(rgb(0xe5b879))
                                .child(notice),
                        )
                    })
                    .child(u.sequence.read(cx).translation_legend()),
            )
            .into_any_element()
    }

    fn feature_table(&self, primers: bool, cx: &mut Context<Self>) -> AnyElement {
        let t = self.model.current().unwrap();
        let u = &self.ui[&t.id];
        let query = u.filter.read(cx).value().to_string();
        let doc = std::rc::Rc::new(t.editor.document.clone());
        let features: Vec<_> = doc
            .features
            .iter()
            .enumerate()
            .filter(|(_, f)| {
                (!primers || f.kind == "primer_bind") && feature_card::matches_query(f, &query)
            })
            .map(|(i, _)| i)
            .collect();
        let count = features.len();
        let selected = t.feature;
        let weak = cx.entity().downgrade();
        let list = uniform_list(
            if primers {
                "primers-list"
            } else {
                "features-list"
            },
            count,
            move |range, _, _| {
                range
                    .map(|row| {
                        let index = features[row];
                        let f = doc.features[index].clone();
                        let color = feature_card::feature_color(&f, index);
                        let name = f.label.clone();
                        let location = f.location_label();
                        let oligo = primers.then(|| f.primer_sequence()).flatten();
                        let secondary = if primers {
                            oligo.as_ref().map_or_else(
                                || "Oligo sequence not recorded".to_string(),
                                |sequence| format!("5′ {sequence} 3′"),
                            )
                        } else {
                            f.kind.clone()
                        };
                        let length = if primers {
                            oligo.as_ref().map_or_else(
                                || "Unknown oligo".to_string(),
                                |sequence| format!("{} nt oligo", sequence.len()),
                            )
                        } else {
                            format!("{} {}", f.len(), doc.unit())
                        };
                        let bound_length = format!("{} {} bound", f.len(), doc.unit());
                        let w = weak.clone();
                        let id = f.id;
                        let parts = f.parts.clone();
                        let d = doc.clone();
                        let card = f;
                        let trigger = div()
                            .id(("feature-row", id))
                            .h_flex()
                            .h(px(if primers { 62. } else { 52. }))
                            .px_4()
                            .gap_3()
                            .border_b_1()
                            .border_color(rgb(BORDER))
                            .cursor_pointer()
                            .bg(rgb(if selected == Some(id) { 0x213c43 } else { BG }))
                            .hover(|s| s.bg(rgb(0x223040)))
                            .on_click(move |event, wnd, cx| {
                                let parts = parts.clone();
                                let double = event.click_count() == 2;
                                let _ = w.update(cx, |s, cx| {
                                    if s.blocked() {
                                        return;
                                    }
                                    if let Some(t) = s.model.current_mut() {
                                        t.select(parts, Some(id), true);
                                        if double {
                                            t.show_sequence();
                                        }
                                    }
                                    s.focus.focus(wnd, cx);
                                    cx.notify();
                                });
                            })
                            .child(div().w(px(4.)).h(px(24.)).rounded_sm().bg(color))
                            .child(
                                div()
                                    .flex_1()
                                    .min_w_0()
                                    .v_flex()
                                    .child(
                                        div()
                                            .font_weight(FontWeight::MEDIUM)
                                            .truncate()
                                            .child(name),
                                    )
                                    .child(
                                        div()
                                            .text_xs()
                                            .text_color(rgb(MUTED))
                                            .truncate()
                                            .when(primers, |d| d.font_family("monospace"))
                                            .child(secondary),
                                    ),
                            )
                            .child(
                                div()
                                    .w(px(185.))
                                    .text_xs()
                                    .text_color(rgb(MUTED))
                                    .truncate()
                                    .child(location),
                            )
                            .child(
                                div()
                                    .w(px(if primers { 115. } else { 74. }))
                                    .v_flex()
                                    .text_right()
                                    .text_xs()
                                    .child(length)
                                    .when(primers, |d| {
                                        d.child(div().text_color(rgb(MUTED)).child(bound_length))
                                    }),
                            )
                            .when(card.edited, |d| {
                                d.child(div().text_xs().text_color(rgb(0xe0b472)).child("Review"))
                            });
                        HoverCard::new(("feature-preview", id))
                            .open_delay(Duration::from_millis(300))
                            .trigger(trigger)
                            .content(move |_, _, _| {
                                feature_card::render_feature_card(&d, &card, None)
                            })
                            .into_any_element()
                    })
                    .collect::<Vec<_>>()
            },
        )
        .size_full();
        div()
            .v_flex()
            .flex_1()
            .min_h_0()
            .px_6()
            .pb_5()
            .gap_3()
            .child(
                div()
                    .h_flex()
                    .gap_3()
                    .child(
                        div().w(px(270.)).child(
                            Input::new(&u.filter)
                                .small()
                                .prefix(Icon::new(IconName::Search).small())
                                .aria_label("Filter annotations by name, type, sequence or qualifier"),
                        ),
                    )
                    .child(div().text_xs().text_color(rgb(MUTED)).child(format!(
                        "{count} {}",
                        if primers {
                            "primer binding sites"
                        } else {
                            "annotations"
                        }
                    ))),
            )
            .child(
                div()
                    .h_flex()
                    .px_4()
                    .text_xs()
                    .text_color(rgb(MUTED))
                    .child(div().flex_1().child(if primers { "NAME / OLIGO (5′ → 3′)" } else { "NAME / TYPE" }))
                    .child(div().w(px(185.)).child(if primers { "BINDING / STRAND" } else { "LOCATION" }))
                    .child(
                        div()
                            .w(px(if primers { 115. } else { 74. }))
                            .text_right()
                            .child(if primers { "OLIGO / BINDING" } else { "LENGTH" }),
                    ),
            )
            .child(
                div()
                    .flex_1()
                    .min_h_0()
                    .border_1()
                    .border_color(rgb(BORDER))
                    .rounded_lg()
                    .overflow_hidden()
                    .child(list),
            )
            .when(count == 0, |d| {
                d.child(
                    div()
                        .text_sm()
                        .text_color(rgb(MUTED))
                        .child("No matching annotations."),
                )
            })
            .when_some(u.feature_details_saved.as_ref(), |d, (feature, _)| {
                let oligo = feature.primer_sequence().map(|s| s.into_owned());
                d.child(
                    div()
                        .v_flex()
                        .gap_2()
                        .flex_shrink_0()
                        .child(
                            div()
                                .h_flex()
                                .gap_2()
                                .child(
                                    div()
                                        .flex_1()
                                        .min_w_0()
                                        .text_sm()
                                        .font_weight(FontWeight::MEDIUM)
                                        .truncate()
                                        .child(format!("Selected annotation · {}", feature.label)),
                                )
                                .when_some(oligo, |d, oligo| {
                                    d.child(
                                        Button::new("copy-primer-oligo")
                                            .ghost()
                                            .small()
                                            .label("Copy oligo")
                                            .disabled(self.blocked())
                                            .on_click(move |_, _, cx| {
                                                cx.write_to_clipboard(ClipboardItem::new_string(oligo.clone()));
                                            }),
                                    )
                                })
                                .child(
                                    Button::new("copy-feature-details")
                                        .ghost()
                                        .small()
                                        .label("Copy details")
                                        .disabled(self.blocked())
                                        .on_click(cx.listener(|s, _, _, cx| {
                                            if let Some(tab) = s.model.current()
                                                && let Some(ui) = s.ui.get(&tab.id)
                                            {
                                                cx.write_to_clipboard(ClipboardItem::new_string(
                                                    ui.feature_details.read(cx).value().to_string(),
                                                ));
                                            }
                                        })),
                                ),
                        )
                        .child(
                            Textarea::new(&u.feature_details)
                                .readonly(true)
                                .h(px(200.))
                                .aria_label("Complete selected annotation details"),
                        ),
                )
            })
            .child(div().text_xs().text_color(rgb(MUTED)).child(if primers {
                "Click a row for the complete oligo, segments and stored metadata · Tm is not calculated"
            } else {
                "Click for all segments and qualifiers · double-click to show the annotated sequence"
            }))
            .into_any_element()
    }
    fn trace_parameters(&self, _cx: &mut Context<Self>) -> AnyElement {
        let tab = self.model.current().unwrap();
        let ui = &self.ui[&tab.id];
        let trace = tab.editor.document.trace.as_ref().unwrap();
        div().v_flex().flex_1().min_h_0().px_6().pb_6().gap_3()
            .child(div().text_lg().child("Trace parameters"))
            .child(div().text_sm().text_color(rgb(MUTED)).child("Read-only values from the original ABIF file. Missing instrument fields remain absent."))
            .child(format!("{} calls · {} samples/channel · channel order {} · PBAS/PCON/PLOC{}", tab.editor.document.sequence.len(), trace.channels.first().map_or(0, Vec::len), trace.channel_order, trace.calls_version))
            .when_some(ui.source_metadata.as_ref(), |d, metadata| d.child(
                Textarea::new(metadata).readonly(true).flex_1().min_h_0().aria_label("Original ABIF trace parameters")
            ))
            .when(ui.source_metadata.is_none() && ui.source_metadata_error.is_none(), |d| d.child("No supported instrument parameters are stored in this file."))
            .when_some(ui.source_metadata_error.as_ref(), |d, error| d.child(error.clone()))
            .into_any_element()
    }
    fn history(&self, cx: &mut Context<Self>) -> AnyElement {
        let tab = self.model.current().unwrap();
        let ui = &self.ui[&tab.id];
        let doc = &tab.editor.document;
        let details = div()
            .v_flex()
            .gap_3()
            .p_6()
            .child(
                div()
                    .text_lg()
                    .font_weight(FontWeight::SEMIBOLD)
                    .child("Document details"),
            )
            .child(doc.summary())
            .child(
                div()
                    .h_flex()
                    .gap_4()
                    .child(stat(
                        "LENGTH",
                        format!("{} {}", number(doc.sequence.len()), doc.unit()),
                    ))
                    .when(doc.molecule != Molecule::Protein, |d| {
                        d.child(stat("GC CONTENT", format!("{:.1}%", tab.gc)))
                    }),
            )
            .child(eyebrow("DESCRIPTION"))
            .child(
                Textarea::new(&self.ui[&tab.id].description)
                    .h(px(100.))
                    .disabled(tab.busy())
                    .readonly(self.blocked())
                    .aria_label("Document description"),
            )
            .when_some(doc.accession.as_ref(), |d, accession| {
                d.child(stat("ACCESSION", accession.clone()))
            })
            .when_some(ui.source_metadata.as_ref(), |d, metadata| {
                d.child(
                    div()
                        .h_flex()
                        .gap_3()
                        .child(div().flex_1().child(eyebrow("ORIGINAL IMPORTED METADATA")))
                        .child(
                            Button::new("copy-source-metadata")
                                .ghost()
                                .small()
                                .label("Copy metadata")
                                .disabled(self.blocked())
                                .on_click(cx.listener(|s, _, _, cx| {
                                    if let Some(tab) = s.model.current()
                                        && let Some(metadata) = s.ui.get(&tab.id).and_then(|ui| ui.source_metadata.as_ref())
                                    {
                                        cx.write_to_clipboard(ClipboardItem::new_string(metadata.read(cx).value().to_string()));
                                    }
                                })),
                        ),
                )
                .child(div().text_xs().text_color(rgb(MUTED)).child("From the original imported file. These fields and references are not updated by edits or format conversions."))
                .child(Textarea::new(metadata).readonly(true).h(px(240.)).aria_label("Original imported document metadata and references"))
            })
            .when_some(ui.source_metadata_error.as_ref(), |d, error| {
                d.child(div().text_sm().text_color(rgb(MUTED)).child(format!("Original metadata could not be read: {error}")))
            })
            .when(!doc.warnings.is_empty(), |d| {
                d.child(eyebrow("COMPATIBILITY NOTES"))
                    .children(doc.warnings.iter().map(|warning| {
                        div()
                            .text_sm()
                            .text_color(rgb(MUTED))
                            .child(warning.clone())
                    }))
            })
            .when(!doc.comments.is_empty(), |d| {
                d.child(eyebrow("COMMENTS")).children(
                    doc.comments
                        .iter()
                        .map(|comment| div().text_sm().child(comment.clone())),
                )
            });
        let mut panel=div().v_flex().gap_4().p_6().rounded_lg().border_1().border_color(rgb(BORDER)).bg(rgb(PANEL))
            .child(div().text_lg().font_weight(FontWeight::SEMIBOLD).child("Source & history"))
            .child(div().text_sm().text_color(rgb(MUTED)).child("Original imported bytes are available during this session. Standard-format saves contain only what that format supports; keep the original file for unsupported history and metadata."));
        if let Some(source) = &tab.original_source {
            panel=panel.child(div().font_weight(FontWeight::MEDIUM).child(source.file_name.clone()))
                .child(div().text_sm().text_color(rgb(MUTED)).child(format!("{} · {} bytes retained",source.format,number(source.bytes.len()))))
                .child(div().text_sm().child(if source.retained_blocks.iter().any(|(k,_)|matches!(k,7|11)){"SnapGene procedure history is present in the original source file. Its operations and graph are not yet interpreted."}else{"No interpreted procedure history is available for this document."}))
                .child(Button::new("recover-source").small().icon(gpui_kit::assets::IconName::Download).label("Save a copy of the original…").disabled(self.model.current().unwrap().busy()).on_click(cx.listener(|s,_,_,cx|{
                    s.recover_source(cx);
                })))
                .children(source.retained_blocks.iter().map(|(k,n)|div().text_xs().text_color(rgb(MUTED)).child(format!("Block {k:02X} · {} bytes",number(*n)))));
        } else {
            panel = panel.child(
                "This sequence was created in PopGene and has no imported source attachment.",
            );
        }
        div().id("history-scroll").flex_1().min_h_0().overflow_y_scroll().px_6().pb_6().child(details).child(panel)
            .child(div().py_4().text_sm().text_color(rgb(MUTED)).child("Edits do not alter the retained original. A recovered source copy contains the imported state. Undo is limited to 32 steps / 32 MiB of sequence snapshots in this session."))
            .into_any_element()
    }
    fn modal(&self, cx: &mut Context<Self>) -> Option<AnyElement> {
        let content = if let Some(Modal::Export(form)) = &self.modal {
            div().child(form.clone()).into_any_element()
        } else if let Some(index) = self.model.close_target() {
            let t = &self.model.tabs[index];
            let busy = self.model.waiting_save.is_some() || self.file_dialog_open;
            div().v_flex().gap_5().p_6().w(px(500.)).child(div().text_xl().font_weight(FontWeight::SEMIBOLD).child("Save changes before closing?"))
                .child(format!("{} has unsaved {}.",t.title(),if t.draft{"creation settings"}else{"changes"}))
                .when(t.draft,|d|d.child(div().text_sm().text_color(rgb(MUTED)).child("Create the sequence first to save it as a document. Cancel returns to the draft.")))
                .when(t.error,|d|d.child(div().text_sm().text_color(rgb(0xf08c8c)).child(t.status.clone())))
                .child(div().h_flex().gap_2().justify_end()
                    .child(Button::new("cancel-close").ghost().label("Cancel").disabled(busy).on_click(cx.listener(|s,_,_,cx|{s.cancel_close();cx.notify();})))
                    .child(Button::new("discard-close").label("Discard").disabled(busy).on_click(cx.listener(|s,_,_,cx|{s.model.discard_close();s.finish_quit(cx);cx.notify();})))
                    .child(Button::new("save-close").primary().label(if busy{"Saving…"}else{"Save & close"}).disabled(busy||t.draft).on_click(cx.listener(|s,_,_,cx|{
                        s.save_before_close(cx);
                    })))).into_any_element()
        } else {
            match self.modal.as_ref()?{
            Modal::Edit(e)=>div().child(e.clone()).into_any_element(),
            Modal::Annotation(e)=>div().child(e.clone()).into_any_element(),
            Modal::Export(e)=>div().child(e.clone()).into_any_element(),
            Modal::UpdateInfo(message)=>div().v_flex().gap_4().w(px(460.)).p_6()
                .child(div().text_xl().child("Check for Updates"))
                .child(message.clone())
                .child(Button::new("close-update-info").primary().label("Done").on_click(cx.listener(|s,_,_,cx|{s.modal=None;cx.notify();}))).into_any_element(),
            Modal::Report=>div().v_flex().gap_4().w(px(640.)).p_6().child(div().text_xl().child("Project import report"))
                .child(div().text_sm().text_color(rgb(MUTED)).child("Files that failed to open remain untouched. Other documents are still available."))
                .child(div().id("report-scroll").max_h(px(350.)).overflow_y_scroll().v_flex().gap_3().children(self.model.issues.iter().map(|issue|div().text_sm().child(issue.clone())))
                    .when(self.model.issues.is_empty(),|d|d.child("No file-open failures recorded.")))
                .child(Button::new("close-report").primary().label("Done").on_click(cx.listener(|s,_,_,cx|{s.modal=None;cx.notify();}))).into_any_element(),
            Modal::About=>div().v_flex().gap_4().w(px(460.)).p_6().child(div().text_2xl().child("PopGene"))
                .child(concat!(env!("CARGO_PKG_VERSION"), " · Native sequence workspace"))
                .child(div().text_sm().text_color(rgb(MUTED)).child("Built in Rust with GPUI and GPUI Kit. Scientific workflows and file storage run locally. No account required."))
                .child("GPL-3.0-or-later · Open source")
                .child(div().text_sm().text_color(rgb(MUTED)).child("github.com/systemantiks/PopGene"))
                .child(Button::new("close-about").primary().label("Done").on_click(cx.listener(|s,_,_,cx|{s.modal=None;cx.notify();}))).into_any_element(),
        }
        };
        Some(
            div()
                .absolute()
                .inset_0()
                .size_full()
                .flex()
                .items_center()
                .justify_center()
                .bg(rgba(0x050a13bb))
                .occlude()
                .on_mouse_down(MouseButton::Left, |_, _, cx| cx.stop_propagation())
                .child(
                    div()
                        .id("modal-content")
                        .w(px(740.))
                        .max_w(relative(0.94))
                        .max_h(relative(0.9))
                        .overflow_y_scroll()
                        .p_6()
                        .rounded_xl()
                        .border_1()
                        .border_color(rgb(BORDER))
                        .bg(rgb(PANEL))
                        .shadow_lg()
                        .child(content)
                        .focus_trap("workspace-modal", &self.modal_focus),
                )
                .into_any_element(),
        )
    }
}
impl Render for PopGene {
    fn render(&mut self, window: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        let sequence_sites_tab = self
            .model
            .current()
            .filter(|tab| {
                tab.view.has_sequence()
                    && self
                        .ui
                        .get(&tab.id)
                        .is_none_or(|ui| ui.sequence_options.show_enzymes)
            })
            .map(|tab| tab.id);
        self.model
            .ensure_current_analysis(self.enzyme_source.map(|(id, _, _)| id), sequence_sites_tab);
        if self.model.has_background_work() && !self.polling {
            self.wake(cx);
        }
        self.sync_ui(window, cx);
        let tabs = self.tabs(cx);
        let toolbar = self.toolbar(cx);
        let navigator = self.navigator(cx);
        let content = self.content(cx);
        let enzyme_visible = self.enzyme_source.is_some();
        let enzyme_panel = if self.enzyme_revision.is_some() {
            self.enzyme_sidebar
                .as_ref()
                .map(|sidebar| div().size_full().child(sidebar.clone()).into_any_element())
                .unwrap_or_else(|| div().into_any_element())
        } else {
            div()
                .size_full()
                .p_4()
                .v_flex()
                .gap_3()
                .child(
                    self.model
                        .current()
                        .map_or("Restriction details unavailable", |tab| {
                            if !tab.restriction_eligible() {
                                "Restriction analysis requires nonempty double-stranded DNA"
                            } else {
                                tab.restriction_error
                                    .as_deref()
                                    .unwrap_or("Updating restriction sites…")
                            }
                        })
                        .to_owned(),
                )
                .child(
                    Button::new("close-pending-enzyme")
                        .small()
                        .label("Close details")
                        .on_click(cx.listener(|s, _, window, cx| {
                            s.close_enzyme_sidebar();
                            s.focus.focus(window, cx);
                            cx.notify();
                        })),
                )
                .into_any_element()
        };
        let status = self
            .model
            .current()
            .map_or_else(|| self.model.status.clone(), |t| t.status.clone());
        let error = self.model.current().is_some_and(|t| t.error);
        let panel = h_resizable("workspace-panels")
            .child(
                resizable_panel()
                    .size(px(238.))
                    .size_range(px(180.)..px(360.))
                    .visible(self.navigator)
                    .child(navigator),
            )
            .child(
                resizable_panel()
                    .size_range(px(420.)..px(10000.))
                    .child(content),
            )
            .child(
                resizable_panel()
                    .size(px(390.))
                    .size_range(px(340.)..px(520.))
                    .visible(enzyme_visible)
                    .child(enzyme_panel),
            );
        let modal = self.modal(cx);
        div()
            .id("popgene-workspace")
            .relative()
            .size_full()
            .v_flex()
            .bg(rgb(BG))
            .text_color(rgb(INK))
            .text_sm()
            .font_family(cx.theme().font_family.clone())
            .key_context("PopGene")
            .track_focus(&self.focus)
            .on_action(cx.listener(|s, _: &NewSequence, _, cx| s.new_sequence(cx)))
            .on_action(cx.listener(|s, _: &OpenFiles, _, cx| s.open(cx)))
            .on_action(cx.listener(|s, _: &OpenFolder, _, cx| s.folder(cx)))
            .on_action(cx.listener(|s, _: &Save, _, cx| s.save(false, cx)))
            .on_action(cx.listener(|s, _: &SaveAs, _, cx| s.save(true, cx)))
            .on_action(cx.listener(|s, _: &Export, w, cx| s.export(w, cx)))
            .on_action(cx.listener(|s, _: &CloseTab, _, cx| {
                if let Some(t) = s.model.current() {
                    s.close_tab(t.id, cx);
                }
            }))
            .on_action(cx.listener(|s, _: &Quit, _, cx| s.request_quit(cx)))
            .on_action(cx.listener(|s, _: &CheckForUpdates, _, cx| s.check_for_updates(cx)))
            .on_action(cx.listener(|s, _: &Undo, _, cx| {
                if s.mutable() {
                    s.model.current_mut().unwrap().undo();
                    cx.notify();
                }
            }))
            .on_action(cx.listener(|s, _: &Redo, _, cx| {
                if s.mutable() {
                    s.model.current_mut().unwrap().redo();
                    cx.notify();
                }
            }))
            .on_action(cx.listener(|s, _: &NextTab, w, cx| s.cycle(true, w, cx)))
            .on_action(cx.listener(|s, _: &PreviousTab, w, cx| s.cycle(false, w, cx)))
            .on_action(cx.listener(|s, _: &FindSequence, w, cx| {
                if s.blocked() {
                    return;
                }
                if let Some(t) = s.model.current_mut()
                    && !t.draft
                {
                    if t.view == View::Map {
                        t.map_find = true;
                    } else {
                        t.show_sequence();
                    }
                    s.ui[&t.id].search.read(cx).focus_handle(cx).focus(w, cx);
                    cx.notify();
                }
            }))
            .on_action(cx.listener(|s, _: &Escape, w, cx| {
                if s.model.waiting_save.is_none() && !s.file_dialog_open {
                    let had_overlay = s.modal.is_some() || s.model.pending_close.is_some();
                    s.cancel_close();
                    if !had_overlay {
                        s.close_enzyme_sidebar();
                    }
                    s.modal = None;
                    s.focus.focus(w, cx);
                    cx.notify();
                }
            }))
            .on_drop(cx.listener(|s, paths: &ExternalPaths, _, cx| {
                if !s.blocked() {
                    s.model.open_paths(paths.paths().to_vec());
                    s.wake(cx);
                }
            }))
            .child(
                TitleBar::new().bg(rgb(BG)).child(
                    div()
                        .h_flex()
                        .w_full()
                        .gap_2()
                        .px_3()
                        .child(
                            div()
                                .font_weight(FontWeight::BOLD)
                                .text_color(rgb(ACCENT))
                                .child("PopGene"),
                        )
                        .child(
                            Button::new("navigator")
                                .accessibility_label("Toggle project navigator")
                                .ghost()
                                .small()
                                .icon(IconName::PanelLeft)
                                .tooltip("Toggle project navigator")
                                .on_click(cx.listener(|s, _, _, cx| {
                                    s.navigator = !s.navigator;
                                    cx.notify();
                                })),
                        )
                        .child(div().flex_1())
                        .child(
                            Button::new("about")
                                .accessibility_label("About PopGene")
                                .ghost()
                                .xsmall()
                                .icon(IconName::Info)
                                .tooltip("About PopGene")
                                .on_click(cx.listener(|s, _, _, cx| {
                                    if !s.blocked() {
                                        s.modal = Some(Modal::About);
                                        cx.notify();
                                    }
                                })),
                        ),
                ),
            )
            .child(tabs)
            .child(toolbar)
            .child(div().flex_1().min_h_0().child(panel))
            .child(
                div()
                    .h_flex()
                    .h(px(28.))
                    .px_4()
                    .gap_3()
                    .bg(rgb(PANEL))
                    .border_t_1()
                    .border_color(rgb(BORDER))
                    .text_xs()
                    .child(div().size(px(5.)).rounded_full().bg(rgb(if error {
                        0xf08c8c
                    } else {
                        ACCENT
                    })))
                    .child(
                        div()
                            .flex_1()
                            .min_w_0()
                            .truncate()
                            .text_color(rgb(if error { 0xf08c8c } else { MUTED }))
                            .child(status),
                    )
                    .when(self.model.is_busy(), |d| d.child("Working…"))
                    .child(div().text_color(rgb(MUTED)).child(format!(
                        concat!("{} documents · PopGene ", env!("CARGO_PKG_VERSION")),
                        self.model.tabs.len()
                    ))),
            )
            .children(modal)
    }
}
/// Fixed rails keep filters and asynchronous notices from resizing the canvas.
fn viewer_rail(id: &'static str, tools: bool) -> Stateful<Div> {
    div()
        .id(id)
        .v_flex()
        .w(px(if tools { 144. } else { 44. }))
        .flex_shrink_0()
        .min_h_0()
        .overflow_y_scroll()
        .p_2()
        .gap_2()
        .bg(rgb(PANEL))
        .border_color(rgb(BORDER))
        .when(tools, |d| d.border_l_1())
        .when(!tools, |d| d.border_r_1().items_center())
}

fn input(
    text: &str,
    placeholder: &'static str,
    window: &mut Window,
    cx: &mut App,
) -> Entity<InputState> {
    cx.new(|cx| {
        InputState::new(window, cx)
            .placeholder(placeholder)
            .default_value(text.to_owned())
    })
}
fn eyebrow(text: impl Into<SharedString>) -> Div {
    div()
        .text_xs()
        .font_weight(FontWeight::SEMIBOLD)
        .text_color(rgb(MUTED))
        .child(text.into())
}

// Kit buttons center their inner content independently of the outer style.
// A full-width child keeps directory and document labels on a common left edge.
fn navigation_label(text: impl Into<SharedString>, icon: Option<IconName>) -> Div {
    div()
        .h_flex()
        .w_full()
        .min_w_0()
        .gap_2()
        .text_left()
        .when_some(icon, |d, icon| {
            d.child(Icon::new(icon).small().flex_shrink_0())
        })
        .child(div().flex_1().min_w_0().truncate().child(text.into()))
}
fn stat(label: &'static str, value: String) -> Div {
    div().v_flex().gap_1().child(eyebrow(label)).child(value)
}
fn safe_name(name: &str) -> String {
    name.replace(['/', '\\'], "_")
}
fn cut_label(position: usize, length: usize, circular: bool) -> String {
    if position == 0 {
        if circular {
            number(length)
        } else {
            "before 1".into()
        }
    } else {
        number(position)
    }
}
fn sequence_sites(t: &crate::model::Tab, visible: bool) -> Vec<SequenceRestrictionSite> {
    if !visible {
        return vec![];
    }
    let mut grouped: HashMap<(usize, usize, Option<usize>, Option<usize>), usize> = HashMap::new();
    let mut sites: Vec<SequenceRestrictionSite> = Vec::new();
    for site in t.displayed_sites() {
        let enzyme = &ENZYMES[site.enzyme_index];
        let key = (
            site.start,
            enzyme.recognition.len(),
            site.top_cut,
            site.bottom_cut,
        );
        if t.nonredundant
            && let Some(index) = grouped.get(&key)
        {
            sites[*index].label.push_str(" / ");
            sites[*index].label.push_str(enzyme.name);
            continue;
        }
        grouped.insert(key, sites.len());
        sites.push(SequenceRestrictionSite {
            enzyme_index: site.enzyme_index,
            start: site.start,
            top_cut: site.top_cut,
            bottom_cut: site.bottom_cut,
            label: enzyme.name.into(),
        });
    }
    sites.sort_by_key(|site| (site.top_cut, site.start, site.enzyme_index));
    sites
}
fn map_sites(t: &crate::model::Tab) -> Vec<MapRestrictionSite> {
    if !t.show_sites {
        return vec![];
    }
    let mut grouped: HashMap<(usize, usize, usize, usize), usize> = HashMap::new();
    let mut sites: Vec<MapRestrictionSite> = Vec::new();
    for site in t.displayed_sites() {
        let e = &ENZYMES[site.enzyme_index];
        let key = (
            site.start,
            e.recognition.len(),
            site.top_cut.unwrap(),
            site.bottom_cut.unwrap(),
        );
        if t.nonredundant
            && let Some(index) = grouped.get(&key)
        {
            sites[*index].label.push_str(" / ");
            sites[*index].label.push_str(e.name);
            continue;
        }
        grouped.insert(key, sites.len());
        sites.push(MapRestrictionSite {
            label: e.name.into(),
            position: site.top_cut.unwrap(),
            start: site.start,
            length: e.recognition.len(),
            reverse: site.reverse,
            enzyme_index: site.enzyme_index,
        });
    }
    sites.sort_by_key(|site| site.position);
    sites
}
pub fn number(n: usize) -> String {
    let s = n.to_string();
    let mut out = String::new();
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

#[cfg(test)]
mod map_control_tests {
    use super::{Document, cut_label, map_sites, sequence_sites};
    #[test]
    fn grouping_combines_equivalent_cuts_but_preserves_distinct_cleavage_offsets() {
        let doc = Document::new("cuts", "AAAACTCGAGTTTTGGTACCAAAA").unwrap();
        let mut tab = crate::model::Tab::from_document(1, doc);
        tab.restriction = Some(std::sync::Arc::new(popgene_core::restriction::find_sites(
            &tab.editor.document.sequence,
            false,
        )));
        let grouped = map_sites(&tab);
        assert!(
            grouped
                .iter()
                .any(|s| s.label.contains("PaeR7I") && s.label.contains("XhoI"))
        );
        assert!(grouped.iter().any(|s| s.label == "Acc65I"));
        assert!(grouped.iter().any(|s| s.label == "KpnI"));
        tab.nonredundant = false;
        let separate = map_sites(&tab);
        assert!(separate.iter().any(|s| s.label == "PaeR7I"));
        assert!(separate.iter().any(|s| s.label == "XhoI"));
        assert!(!tab.dirty());
    }
    #[test]
    fn sequence_sites_share_cutter_grouping_but_not_map_visibility() {
        let doc = Document::new("Sequence labels", "CTCGAGGGTACC").unwrap();
        let mut tab = crate::model::Tab::from_document(1, doc);
        tab.restriction = Some(std::sync::Arc::new(popgene_core::restriction::find_sites(
            &tab.editor.document.sequence,
            false,
        )));
        tab.show_sites = false;
        assert!(map_sites(&tab).is_empty());
        let sites = sequence_sites(&tab, true);
        assert!(
            sites
                .iter()
                .any(|site| site.label.contains("PaeR7I") && site.label.contains("XhoI"))
        );
        let acc = sites.iter().find(|site| site.label == "Acc65I").unwrap();
        let kpn = sites.iter().find(|site| site.label == "KpnI").unwrap();
        assert_eq!(acc.start, kpn.start);
        assert_ne!(acc.top_cut, kpn.top_cut);
        assert_ne!(acc.bottom_cut, kpn.bottom_cut);
        assert!(sequence_sites(&tab, false).is_empty());
        assert!(!tab.dirty());
    }

    #[test]
    fn cleavage_boundary_labels_do_not_add_a_residue_offset() {
        assert_eq!(cut_label(2345, 5000, true), "2,345");
        assert_eq!(cut_label(0, 5000, true), "5,000");
        assert_eq!(cut_label(0, 5000, false), "before 1");
    }

    #[test]
    fn imported_fixture_maps_keep_independently_verified_cleavage_boundaries() {
        // Supplier offsets, checked independently against Biopython's next-base
        // coordinates: FspAI RTGC|GCAY; PaqCI CACCTGC(4/8); NmeAIII GCCGAG(21/19).
        // NmeAIII's context-dependent alternative is not inferred from sequence.
        for (name, bytes) in [
            (
                "annotated.gb",
                include_bytes!("../../../validation/fixtures/annotated.gb").as_slice(),
            ),
            (
                "annotated.dna",
                include_bytes!("../../../validation/fixtures/annotated.dna").as_slice(),
            ),
        ] {
            let doc = popgene_formats::import(name, bytes).unwrap();
            let mut tab = crate::model::Tab::from_document(1, doc);
            tab.restriction = Some(std::sync::Arc::new(popgene_core::restriction::find_sites(
                &tab.editor.document.sequence,
                tab.editor.document.circular,
            )));
            let sites = map_sites(&tab);
            for (enzyme, position, label) in [
                ("EcoRI", 2345, "2,345"),
                ("FspAI", 3270, "3,270"),
                ("PaqCI", 3444, "3,444"),
                ("NmeAIII", 1583, "1,583"),
            ] {
                let site = sites.iter().find(|site| site.label == enzyme).unwrap();
                assert_eq!(site.position, position, "{name}: {enzyme}");
                assert_eq!(cut_label(site.position, 5000, true), label);
            }
        }
    }
}
