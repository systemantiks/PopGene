//! Toolkit-independent document and workspace state.
//!
//! File operations run on workers; the frontend polls for completion and redraws only
//! when something changed. A tab's editor remains locked while its snapshot is saved.
use popgene_core::restriction::{self, ENZYMES, RestrictionAnalysis, RestrictionSite};
use popgene_core::{Document, Editor, Error, Molecule, Segment, Source};
use popgene_formats::{self as formats, SaveFormat, project::Project};
use std::{
    collections::{BTreeSet, VecDeque},
    path::{Path, PathBuf},
    sync::mpsc::{self, Receiver, TryRecvError},
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};

const SAMPLE: &[u8] = include_bytes!("../../../validation/fixtures/annotated.gb");
pub const MAX_TABS: usize = 64;
const MAX_SOURCE_BYTES: usize = 256 * 1024 * 1024;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum View {
    Map,
    Sequence,
    Split,
    Enzymes,
    Features,
    Primers,
    Trace,
    Parameters,
    Properties,
    History,
}
impl View {
    pub fn has_map(self) -> bool {
        matches!(self, Self::Map | Self::Split)
    }
    pub fn has_sequence(self) -> bool {
        matches!(self, Self::Sequence | Self::Split)
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum EnzymeFilter {
    #[default]
    UniqueSix,
    Unique,
    OneOrTwo,
    AllSix,
    All,
}
impl EnzymeFilter {
    pub const ALL: [Self; 5] = [
        Self::UniqueSix,
        Self::Unique,
        Self::OneOrTwo,
        Self::AllSix,
        Self::All,
    ];
    pub fn label(self) -> &'static str {
        match self {
            Self::UniqueSix => "Unique 6+ cutters",
            Self::Unique => "Unique cutters",
            Self::OneOrTwo => "1–2 cutters",
            Self::AllSix => "All 6+ cutters",
            Self::All => "All cutters",
        }
    }
    pub fn includes(self, count: usize, motif_len: usize) -> bool {
        match self {
            Self::UniqueSix => count == 1 && motif_len >= 6,
            Self::Unique => count == 1,
            Self::OneOrTwo => (1..=2).contains(&count),
            Self::AllSix => count > 0 && motif_len >= 6,
            Self::All => count > 0,
        }
    }
}

struct AnalysisJob {
    receiver: Receiver<Result<RestrictionAnalysis, String>>,
    cancel: Arc<AtomicBool>,
}
impl Drop for AnalysisJob {
    fn drop(&mut self) {
        self.cancel.store(true, Ordering::Relaxed);
    }
}

pub struct Tab {
    pub id: u64,
    pub draft: bool,
    pub draft_dirty: bool,
    pub unsaved: bool,
    pub editor: Editor,
    pub saved: Document,
    pub source_path: Option<PathBuf>,
    pub path: Option<PathBuf>,
    pub save_format: Option<SaveFormat>,
    pub disk_hash: Option<[u8; 32]>,
    /// The first imported bytes remain recoverable during this session even after Save As.
    pub original_source: Option<Source>,
    pub view: View,
    pub map_linear: bool,
    pub show_labels: bool,
    pub show_features: bool,
    pub show_primers: bool,
    pub show_positions: bool,
    pub show_sites: bool,
    pub map_zoom: f32,
    pub map_find: bool,
    pub enzyme_filter: EnzymeFilter,
    pub nonredundant: bool,
    pub restriction: Option<Arc<RestrictionAnalysis>>,
    pub restriction_error: Option<String>,
    pub restriction_ambiguous: bool,
    pub selected_enzyme: Option<usize>,
    analysis_job: Option<AnalysisJob>,
    pub selected: Vec<Segment>,
    pub selection_summary: String,
    pub feature: Option<u64>,
    pub range_start: String,
    pub range_end: String,
    /// Zero-based base coordinate; the responsive sequence view chooses the row.
    pub scroll_to: Option<usize>,
    pub query: String,
    pub hits: Vec<Segment>,
    pub hit_index: usize,
    pub search_error: String,
    pub trace_scale: f32,
    pub gc: f32,
    pub status: String,
    pub error: bool,
    job: Option<Receiver<Result<Completed, String>>>,
    job_path: Option<PathBuf>,
}

enum Completed {
    Saved(Box<Document>, PathBuf, SaveFormat, [u8; 32]),
    Exported(PathBuf),
}
impl Tab {
    /// Reveal sequence targets without collapsing an existing split layout.
    pub fn show_sequence(&mut self) {
        if !self.view.has_sequence() {
            self.view = View::Sequence;
        }
    }

    pub fn from_document(id: u64, doc: Document) -> Self {
        let view = if doc.trace.is_some() {
            View::Trace
        } else if doc.molecule == Molecule::Protein || doc.placeholder {
            View::Sequence
        } else {
            View::Map
        };
        let gc = doc.gc();
        let end = doc.sequence.len().min(100).to_string();
        let restriction_ambiguous = doc.sequence.bytes().any(|b| !b"ACGT".contains(&b));
        let original_source = doc.source.clone();
        Self {
            id,
            draft: false,
            draft_dirty: false,
            unsaved: false,
            saved: doc.clone(),
            editor: Editor::new(doc),
            source_path: None,
            path: None,
            save_format: None,
            disk_hash: None,
            original_source,
            view,
            map_linear: false,
            show_labels: true,
            show_features: true,
            show_primers: true,
            show_positions: true,
            show_sites: true,
            map_zoom: 1.0,
            map_find: false,
            enzyme_filter: EnzymeFilter::default(),
            nonredundant: true,
            restriction: None,
            restriction_error: None,
            restriction_ambiguous,
            selected_enzyme: None,
            analysis_job: None,
            selected: vec![],
            selection_summary: String::new(),
            feature: None,
            range_start: "1".into(),
            range_end: end,
            scroll_to: None,
            query: String::new(),
            hits: vec![],
            hit_index: 0,
            search_error: String::new(),
            trace_scale: 1.0,
            gc,
            status: "Ready".into(),
            error: false,
            job: None,
            job_path: None,
        }
    }
    pub fn title(&self) -> String {
        if self.draft {
            return "New Sequence".into();
        }
        self.source_path
            .as_ref()
            .and_then(|p| p.file_name())
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| self.editor.document.name.clone())
    }
    pub fn dirty(&self) -> bool {
        if self.draft {
            self.draft_dirty
        } else {
            self.unsaved || self.editor.document != self.saved
        }
    }
    pub fn busy(&self) -> bool {
        self.job.is_some()
    }
    pub fn message(&mut self, message: impl Into<String>, error: bool) {
        self.status = message.into();
        self.error = error;
    }
    pub fn changed(&mut self) {
        self.analysis_job = None;
        self.restriction = None;
        self.restriction_error = None;
        self.restriction_ambiguous = self
            .editor
            .document
            .sequence
            .bytes()
            .any(|b| !b"ACGT".contains(&b));
        self.selected_enzyme = None;
        self.gc = self.editor.document.gc();
        self.selected.clear();
        self.selection_summary.clear();
        self.feature = None;
        self.hits.clear();
        self.query.clear();
        self.search_error.clear();
        self.scroll_to = Some(0);
        self.range_start = "1".into();
        self.range_end = self.editor.document.sequence.len().min(100).to_string();
    }
    pub fn apply(&mut self, operation: impl FnOnce(&mut Document) -> Result<(), Error>) -> bool {
        if self.busy() || self.draft {
            return false;
        }
        match self.editor.apply(operation) {
            Ok(()) => {
                self.changed();
                self.message("Changes applied", false);
                true
            }
            Err(error) => {
                self.message(error.to_string(), true);
                false
            }
        }
    }
    pub fn undo(&mut self) {
        if self.busy() || self.draft || !self.editor.can_undo() {
            return;
        }
        self.editor.undo();
        self.changed();
        self.message("Undone", false);
    }
    pub fn redo(&mut self) {
        if self.busy() || self.draft || !self.editor.can_redo() {
            return;
        }
        self.editor.redo();
        self.changed();
        self.message("Redone", false);
    }
    pub fn select(&mut self, parts: Vec<Segment>, feature: Option<u64>, scroll: bool) {
        self.selected_enzyme = None;
        let len = self.editor.document.sequence.len();
        if parts.iter().any(|p| p.start >= p.end || p.end > len) {
            self.message("Selection is outside the sequence.", true);
            return;
        }
        if let Some(part) = parts.first() {
            self.range_start = (part.start + 1).to_string();
            self.range_end = part.end.to_string();
            if scroll {
                self.scroll_to = Some(part.start);
            }
        }
        let length: usize = parts.iter().map(|p| p.end - p.start).sum();
        let gc: usize = parts
            .iter()
            .map(|p| {
                self.editor.document.sequence[p.start..p.end]
                    .bytes()
                    .filter(|b| *b == b'G' || *b == b'C')
                    .count()
            })
            .sum();
        let locations = parts
            .iter()
            .map(|p| format!("{}–{}", p.start + 1, p.end))
            .collect::<Vec<_>>()
            .join(", ");
        let name = feature
            .and_then(|id| self.editor.document.features.iter().find(|f| f.id == id))
            .map(|f| format!("{} · ", f.label))
            .unwrap_or_default();
        self.selection_summary = if length == 0 {
            String::new()
        } else {
            format!(
                "{name}{locations} · {length} {}{}",
                if self.editor.document.molecule == Molecule::Dna
                    && !self.editor.document.single_stranded
                {
                    "bp"
                } else {
                    self.editor.document.molecule.unit()
                },
                if self.editor.document.molecule == Molecule::Protein {
                    String::new()
                } else {
                    format!(" · {:.1}% GC", gc as f32 * 100. / length as f32)
                }
            )
        };
        self.selected = parts;
        self.feature = feature;
    }
    pub fn restriction_eligible(&self) -> bool {
        let doc = &self.editor.document;
        !self.draft
            && doc.molecule == Molecule::Dna
            && !doc.single_stranded
            && !doc.placeholder
            && !doc.sequence.is_empty()
    }
    pub fn analyzing(&self) -> bool {
        self.analysis_job.is_some()
    }
    fn ensure_analysis(&mut self) {
        if !self.restriction_eligible()
            || self.restriction.is_some()
            || self.restriction_error.is_some()
            || self.analyzing()
        {
            return;
        }
        let sequence = self.editor.document.sequence.clone();
        let circular = self.editor.document.circular;
        let cancel = Arc::new(AtomicBool::new(false));
        let token = cancel.clone();
        self.analysis_job = Some(AnalysisJob {
            cancel,
            receiver: worker(move || {
                restriction::find_sites_cancellable(&sequence, circular, || {
                    token.load(Ordering::Relaxed)
                })
                .ok_or_else(|| "Analysis canceled".into())
            }),
        });
    }
    fn poll_analysis(&mut self) -> bool {
        let Some(job) = &self.analysis_job else {
            return false;
        };
        let result = match job.receiver.try_recv() {
            Ok(result) => result,
            Err(TryRecvError::Empty) => return false,
            Err(TryRecvError::Disconnected) => {
                Err("Restriction analysis ended unexpectedly.".into())
            }
        };
        self.analysis_job = None;
        match result {
            Ok(result) => self.restriction = Some(Arc::new(result)),
            Err(error) => self.restriction_error = Some(error),
        }
        true
    }
    pub fn displayed_sites(&self) -> Vec<&RestrictionSite> {
        let Some(analysis) = &self.restriction else {
            return vec![];
        };
        analysis
            .sites
            .iter()
            .filter(|site| {
                let enzyme = &ENZYMES[site.enzyme_index];
                self.enzyme_filter
                    .includes(analysis.counts[site.enzyme_index], enzyme.recognition.len())
                    && site.top_cut.is_some()
                    && site.bottom_cut.is_some()
            })
            .collect()
    }
    pub fn site_parts(&self, site: &RestrictionSite) -> Vec<Segment> {
        let n = self.editor.document.sequence.len();
        let end = site.start + ENZYMES[site.enzyme_index].recognition.len();
        if end > n {
            vec![
                Segment {
                    start: site.start,
                    end: n,
                },
                Segment {
                    start: 0,
                    end: end - n,
                },
            ]
        } else {
            vec![Segment {
                start: site.start,
                end,
            }]
        }
    }
    pub fn select_range(&mut self, start: &str, end: &str) -> bool {
        match (start.trim().parse::<usize>(), end.trim().parse::<usize>()) {
            (Ok(start), Ok(end))
                if start > 0 && start <= end && end <= self.editor.document.sequence.len() =>
            {
                self.select(
                    vec![Segment {
                        start: start - 1,
                        end,
                    }],
                    None,
                    true,
                );
                true
            }
            _ => {
                self.message("Enter a valid 1-based inclusive range.", true);
                false
            }
        }
    }
    pub fn selected_sequence(&self) -> String {
        self.selected
            .iter()
            .map(|s| &self.editor.document.sequence[s.start..s.end])
            .collect()
    }
    pub fn find(&mut self, query: String) {
        self.query = query;
        match self.editor.document.find(&self.query) {
            Ok(hits) => {
                self.hits = hits;
                self.hit_index = 0;
                self.search_error.clear();
                self.select_hit();
            }
            Err(error) => {
                self.hits.clear();
                self.search_error = error.to_string();
            }
        }
    }
    pub fn next_hit(&mut self) {
        if !self.hits.is_empty() {
            self.hit_index = (self.hit_index + 1) % self.hits.len();
            self.select_hit();
        }
    }
    pub fn previous_hit(&mut self) {
        if !self.hits.is_empty() {
            self.hit_index = (self.hit_index + self.hits.len() - 1) % self.hits.len();
            self.select_hit();
        }
    }
    fn select_hit(&mut self) {
        if let Some(hit) = self.hits.get(self.hit_index).cloned() {
            let n = self.editor.document.sequence.len();
            let parts = if hit.end > n {
                vec![
                    Segment {
                        start: hit.start,
                        end: n,
                    },
                    Segment {
                        start: 0,
                        end: hit.end - n,
                    },
                ]
            } else {
                vec![hit]
            };
            self.select(parts, None, true);
        }
    }
    fn start_job(
        &mut self,
        path: PathBuf,
        task: impl FnOnce() -> Result<Completed, String> + Send + 'static,
    ) {
        self.job_path = Some(path);
        self.job = Some(worker(task));
    }
    fn poll(&mut self) -> bool {
        let Some(result) = receive(&self.job) else {
            return false;
        };
        self.job = None;
        self.job_path = None;
        match result {
            Err(error) => self.message(error, true),
            Ok(Completed::Saved(doc, path, format, hash)) => {
                // The working document must describe what is actually on disk. Keep the
                // pre-conversion version in Undo so an explicitly accepted lossy Save As
                // can still be recovered during this session.
                let canonical = *doc;
                let content_changed = !same_working_data(&self.editor.document, &canonical);
                if !content_changed {
                    // A new source attachment is bookkeeping, not an undoable edit.
                    self.editor.document = canonical.clone();
                } else {
                    self.editor
                        .apply(|document| {
                            *document = canonical.clone();
                            Ok(())
                        })
                        .expect("saved documents were validated before writing");
                }
                self.saved = canonical;
                if content_changed {
                    self.changed();
                }
                if self.editor.document.trace.is_none() && self.view == View::Trace {
                    self.view = View::Sequence;
                }
                self.unsaved = false;
                self.source_path = Some(path.clone());
                self.path = Some(path);
                self.save_format = Some(format);
                self.disk_hash = Some(hash);
                self.message(format!("Saved · {}", format.label()), false);
            }
            Ok(Completed::Exported(path)) => self.message(
                format!(
                    "Exported {}",
                    path.file_name().unwrap_or_default().to_string_lossy()
                ),
                false,
            ),
        }
        true
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Close {
    Tab(u64),
    Quit,
}
enum Loaded {
    File(Box<Document>, PathBuf, Option<[u8; 32]>),
    Folder(Project),
}
pub struct Model {
    pub tabs: Vec<Tab>,
    pub active: usize,
    pub project: Option<Project>,
    pub status: String,
    pub issues: Vec<String>,
    pub pending_close: Option<Close>,
    pub waiting_save: Option<u64>,
    pub allow_close: bool,
    next_id: u64,
    queue: VecDeque<PathBuf>,
    opening: Option<Receiver<Result<Loaded, String>>>,
    opening_path: Option<PathBuf>,
    discarded: BTreeSet<u64>,
}
impl Model {
    pub fn new() -> Self {
        let mut doc = formats::import("annotated.gb", SAMPLE).expect("checked-in sample is valid");
        doc.name = "Synthetic circular DNA".into();
        let mut tab = Tab::from_document(1, doc);
        tab.status = "Original synthetic sample · ready to explore".into();
        Self {
            tabs: vec![tab],
            active: 0,
            project: None,
            status: "Ready · local files only".into(),
            issues: vec![],
            pending_close: None,
            waiting_save: None,
            allow_close: false,
            next_id: 2,
            queue: VecDeque::new(),
            opening: None,
            opening_path: None,
            discarded: BTreeSet::new(),
        }
    }
    pub fn current(&self) -> Option<&Tab> {
        self.tabs.get(self.active)
    }
    pub fn current_mut(&mut self) -> Option<&mut Tab> {
        self.tabs.get_mut(self.active)
    }
    pub fn is_busy(&self) -> bool {
        self.opening.is_some() || !self.queue.is_empty() || self.tabs.iter().any(Tab::busy)
    }
    pub fn has_background_work(&self) -> bool {
        self.is_busy() || self.tabs.iter().any(Tab::analyzing)
    }
    pub fn ensure_current_analysis(
        &mut self,
        inspected_tab: Option<u64>,
        sequence_sites_tab: Option<u64>,
    ) {
        for (i, tab) in self.tabs.iter_mut().enumerate() {
            if i == self.active
                && (tab.view == View::Enzymes
                    || (tab.view.has_map() && tab.show_sites)
                    || (tab.view.has_sequence() && sequence_sites_tab == Some(tab.id))
                    || inspected_tab == Some(tab.id))
            {
                tab.ensure_analysis();
            } else {
                // Keep cached results, but cancel unnecessary work when switching documents/views.
                tab.analysis_job = None;
            }
        }
    }
    pub fn new_sequence(&mut self) -> Option<u64> {
        if self.pending_close.is_some() {
            return None;
        }
        if self.tabs.len() >= MAX_TABS {
            self.status = "Close a tab before creating another sequence (64-tab limit).".into();
            return None;
        }
        let mut request = popgene_core::creation::Request::new(Molecule::Dna);
        request.mode = popgene_core::creation::Mode::Placeholder;
        let mut tab =
            Tab::from_document(self.next_id, request.create().expect("valid empty draft"));
        self.next_id += 1;
        tab.draft = true;
        tab.status = "New sequence · not created yet".into();
        let id = tab.id;
        self.tabs.push(tab);
        self.active = self.tabs.len() - 1;
        Some(id)
    }
    pub fn create_in_active(&mut self, doc: Document) -> bool {
        let Some(tab) = self.current_mut() else {
            return false;
        };
        if !tab.draft || tab.busy() {
            return false;
        }
        if let Err(error) = doc.validate() {
            tab.message(error.to_string(), true);
            return false;
        }
        let id = tab.id;
        *tab = Tab::from_document(id, doc);
        tab.unsaved = true;
        tab.message("Created · not saved yet", false);
        true
    }
    pub fn open_folder(&mut self, path: PathBuf) {
        self.open_paths(vec![path]);
    }
    pub fn open_paths(&mut self, paths: Vec<PathBuf>) {
        if self.pending_close.is_some() {
            return;
        }
        for path in paths {
            let path = resolved(&path);
            if let Some(i) = self
                .tabs
                .iter()
                .position(|t| t.source_path.as_ref() == Some(&path))
            {
                self.active = i;
                continue;
            }
            if self.opening_path.as_ref() != Some(&path) && !self.queue.contains(&path) {
                self.queue.push_back(path);
            }
        }
        self.start_next();
    }
    fn start_next(&mut self) {
        if self.opening.is_some() {
            return;
        }
        while let Some(path) = self.queue.pop_front() {
            if let Some(i) = self
                .tabs
                .iter()
                .position(|t| t.source_path.as_ref() == Some(&path))
            {
                self.active = i;
                continue;
            }
            self.status = format!(
                "Opening {}…",
                path.file_name().unwrap_or_default().to_string_lossy()
            );
            self.opening_path = Some(path.clone());
            self.opening = Some(worker(move || {
                if path.is_dir() {
                    return formats::project::scan(&path)
                        .map(Loaded::Folder)
                        .map_err(|e| format!("{}: {e}", path.display()));
                }
                // Detect concurrent edits even while loading, rather than associating a
                // newer file hash with an older parsed document.
                let before = formats::digest(&path).map_err(|e| e.to_string())?;
                let doc = formats::open(&path).map_err(|e| format!("{}: {e}", path.display()))?;
                let after = formats::digest(&path).map_err(|e| e.to_string())?;
                if before != after {
                    return Err(format!(
                        "{} changed while opening. Open it again to load the current version.",
                        path.display()
                    ));
                }
                let hash = Some(after);
                Ok(Loaded::File(Box::new(doc), path, hash))
            }));
            break;
        }
    }
    pub fn poll(&mut self) -> bool {
        let mut changed = false;
        for tab in &mut self.tabs {
            changed |= tab.poll();
            changed |= tab.poll_analysis();
        }
        if let Some(result) = receive(&self.opening) {
            changed = true;
            self.opening = None;
            self.opening_path = None;
            match result {
                Ok(Loaded::Folder(project)) => {
                    self.issues.extend(project.warnings.clone());
                    self.status = format!(
                        "Project: {} files; {} supported",
                        project.files,
                        project.supported_files().len()
                    );
                    self.project = Some(project);
                }
                Ok(Loaded::File(doc, path, hash)) => {
                    self.add_document(*doc, Some(path), hash);
                }
                Err(error) => {
                    self.status = "Some files could not be opened — see Import report".into();
                    self.issues.push(error);
                }
            }
            self.start_next();
        }
        if let Some(id) = self.waiting_save {
            if let Some(tab) = self.tabs.iter().find(|t| t.id == id) {
                if !tab.busy() {
                    changed = true;
                    self.waiting_save = None;
                    if !tab.dirty() {
                        self.advance_close();
                    } else {
                        self.status = "Save did not complete; document is still open.".into();
                    }
                }
            } else {
                self.waiting_save = None;
                changed = true;
            }
        }
        changed
    }
    pub fn add_document(
        &mut self,
        doc: Document,
        path: Option<PathBuf>,
        hash: Option<[u8; 32]>,
    ) -> bool {
        let path = path.map(|p| resolved(&p));
        if let Some(ref path) = path
            && let Some(i) = self
                .tabs
                .iter()
                .position(|t| t.source_path.as_ref() == Some(path))
        {
            self.active = i;
            return true;
        }
        let replace = self
            .current()
            .filter(|t| t.draft && !t.dirty())
            .map(|t| t.id);
        let used: usize = self
            .tabs
            .iter()
            .map(|tab| {
                let current = tab.editor.document.source.as_ref();
                let original = tab.original_source.as_ref();
                current.map_or(0, |s| s.bytes.len())
                    + original.map_or(0, |s| {
                        if current.is_some_and(|c| Arc::ptr_eq(&c.bytes, &s.bytes)) {
                            0
                        } else {
                            s.bytes.len()
                        }
                    })
            })
            .sum();
        let incoming = doc.source.as_ref().map_or(0, |s| s.bytes.len());
        if (self.tabs.len() >= MAX_TABS && replace.is_none())
            || used.saturating_add(incoming) > MAX_SOURCE_BYTES
        {
            self.issues.push(format!("{}: open-tab limit reached (64 tabs / 256 MiB source attachments). Close files and try again.", path.as_ref().map_or_else(|| doc.name.clone(), |p| p.display().to_string())));
            self.status = "Open-tab limit reached — see Import report".into();
            return false;
        }
        if let Some(id) = replace {
            self.remove_tab(id);
        }
        let mut tab = Tab::from_document(self.next_id, doc);
        self.next_id += 1;
        tab.source_path = path.clone();
        tab.save_format = path
            .as_deref()
            .and_then(SaveFormat::from_path)
            // Content detection wins over a misleading extension. A GenBank record
            // named .dna must use Save As, not silently become a binary SnapGene file.
            .filter(|format| {
                tab.editor.document.source.as_ref().is_none_or(|source| {
                    matches!(
                        (format, source.format.as_str()),
                        (SaveFormat::SnapGene, "SnapGene")
                            | (SaveFormat::GenBank, "GenBank")
                            | (SaveFormat::Fasta, "FASTA")
                            | (SaveFormat::Embl, "EMBL")
                            | (SaveFormat::Abif, "ABIF")
                    )
                })
            })
            .filter(|format| format.supports(&tab.editor.document));
        tab.path = if hash.is_some() && tab.save_format.is_some() {
            path
        } else {
            None
        };
        tab.disk_hash = hash;
        self.tabs.push(tab);
        self.active = self.tabs.len() - 1;
        self.status = "File opened · original source retained".into();
        true
    }
    fn collision(&self, index: usize, path: &Path, export: bool) -> bool {
        self.tabs.iter().enumerate().any(|(i, tab)| {
            (i != index || export)
                && (tab.source_path.as_deref() == Some(path) || tab.path.as_deref() == Some(path))
                || tab.job_path.as_deref() == Some(path)
        }) || self.opening_path.as_deref() == Some(path)
            || self.queue.iter().any(|p| p == path)
    }
    pub fn save_to(&mut self, index: usize, path: PathBuf, format: SaveFormat) -> bool {
        if !self.tabs.get(index).is_some_and(|t| !t.draft && !t.busy()) {
            return false;
        }
        if SaveFormat::from_path(&path) != Some(format) {
            let extension = format.extension_for(&self.tabs[index].editor.document);
            self.tabs[index].message(
                format!(
                    "Choose a {} filename, such as .{extension}.",
                    format.label()
                ),
                true,
            );
            return false;
        }
        if !format.supports(&self.tabs[index].editor.document) {
            self.tabs[index].message(
                format!(
                    "{} cannot save this document. Choose another format or add a sequence first.",
                    format.label()
                ),
                true,
            );
            return false;
        }
        let path = resolved(&path);
        if self.collision(index, &path, false) {
            self.tabs[index].message("That file is open or being written in another tab. Choose a different save location.", true);
            return false;
        }
        let tab = &mut self.tabs[index];
        let expected = if tab.source_path.as_ref() == Some(&path) {
            tab.disk_hash
        } else {
            None
        };
        let doc = tab.editor.document.clone();
        tab.message(format!("Saving {}…", format.label()), false);
        tab.start_job(path.clone(), move || {
            // Re-read the encoded bytes before replacing any file, so a codec error or
            // misleading extension cannot turn a successful save into unreadable data.
            let canonical = saved_projection(&path, &doc, format)?;
            let hash = formats::save(&path, &doc, format, expected).map_err(|e| e.to_string())?;
            Ok(Completed::Saved(Box::new(canonical), path, format, hash))
        });
        true
    }
    pub fn export_to(&mut self, index: usize, path: PathBuf, format: SaveFormat) -> bool {
        if !self.tabs.get(index).is_some_and(|t| !t.draft && !t.busy()) {
            return false;
        }
        if SaveFormat::from_path(&path) != Some(format) {
            let extension = format.extension_for(&self.tabs[index].editor.document);
            self.tabs[index].message(
                format!(
                    "Choose a {} filename, such as .{extension}.",
                    format.label()
                ),
                true,
            );
            return false;
        }
        if !format.supports(&self.tabs[index].editor.document) {
            self.tabs[index].message(
                format!("{} cannot export this document. Choose another format or add a sequence first.", format.label()),
                true,
            );
            return false;
        }
        let path = resolved(&path);
        if self.collision(index, &path, true) {
            self.tabs[index].message(
                "That file is already open or being written. Export to a different location.",
                true,
            );
            return false;
        }
        let tab = &mut self.tabs[index];
        let doc = tab.editor.document.clone();
        tab.message(format!("Exporting {}…", format.label()), false);
        tab.start_job(path.clone(), move || {
            saved_projection(&path, &doc, format)?;
            formats::save(&path, &doc, format, None).map_err(|e| e.to_string())?;
            Ok(Completed::Exported(path))
        });
        true
    }
    pub fn export_original_to(&mut self, index: usize, path: PathBuf) -> bool {
        let Some(tab) = self.tabs.get(index).filter(|t| !t.draft && !t.busy()) else {
            return false;
        };
        let Some(source) = tab.original_source.as_ref() else {
            return false;
        };
        let bytes = source.bytes.clone();
        let path = resolved(&path);
        if self.collision(index, &path, true) {
            self.tabs[index].message("Choose a destination that is not open in a tab.", true);
            return false;
        }
        let tab = &mut self.tabs[index];
        tab.message("Recovering original source…", false);
        tab.start_job(path.clone(), move || {
            formats::atomic_export(&path, &bytes).map_err(|e| e.to_string())?;
            Ok(Completed::Exported(path))
        });
        true
    }
    pub fn request_close(&mut self, request: Close) -> bool {
        if self.is_busy() {
            self.status = "Wait for file operations to finish before closing.".into();
            return false;
        }
        self.discarded.clear();
        self.pending_close = Some(request);
        self.advance_close();
        true
    }
    pub fn close_target(&self) -> Option<usize> {
        match self.pending_close {
            Some(Close::Tab(id)) => self.tabs.iter().position(|t| t.id == id && t.dirty()),
            Some(Close::Quit) => self
                .tabs
                .iter()
                .position(|t| t.dirty() && !self.discarded.contains(&t.id)),
            None => None,
        }
    }
    fn advance_close(&mut self) {
        if let Some(i) = self.close_target() {
            self.active = i;
            return;
        }
        match self.pending_close.take() {
            Some(Close::Tab(id)) => self.remove_tab(id),
            Some(Close::Quit) => self.allow_close = true,
            None => {}
        }
    }
    pub fn discard_close(&mut self) {
        if self.waiting_save.is_some() {
            return;
        }
        if let Some(i) = self.close_target() {
            let id = self.tabs[i].id;
            match self.pending_close {
                Some(Close::Tab(_)) => {
                    self.remove_tab(id);
                    self.pending_close = None;
                }
                Some(Close::Quit) => {
                    self.discarded.insert(id);
                    self.advance_close();
                }
                None => {}
            }
        }
    }
    pub fn cancel_close(&mut self) {
        self.pending_close = None;
        self.waiting_save = None;
        self.discarded.clear();
    }
    pub fn save_and_close(&mut self, path: PathBuf, format: SaveFormat) -> bool {
        let Some(i) = self.close_target() else {
            return false;
        };
        if self.save_to(i, path, format) {
            self.waiting_save = Some(self.tabs[i].id);
            true
        } else {
            false
        }
    }
    fn remove_tab(&mut self, id: u64) {
        if let Some(i) = self.tabs.iter().position(|t| t.id == id) {
            self.tabs.remove(i);
            if i < self.active {
                self.active -= 1;
            }
            self.active = self.active.min(self.tabs.len().saturating_sub(1));
        }
    }
}
/// Source attachments and importer notices are not sequence edits.
fn same_working_data(left: &Document, right: &Document) -> bool {
    let mut left = left.clone();
    let mut right = right.clone();
    left.source = None;
    right.source = None;
    left.warnings.clear();
    right.warnings.clear();
    left == right
}

/// Validate the public file representation before any destination is replaced.
fn saved_projection(path: &Path, doc: &Document, format: SaveFormat) -> Result<Document, String> {
    let bytes = formats::encode(doc, format).map_err(|e| e.to_string())?;
    let name = path.file_name().unwrap_or_default().to_string_lossy();
    let canonical = formats::import(&name, &bytes).map_err(|e| e.to_string())?;
    if canonical.molecule != doc.molecule || canonical.sequence != doc.sequence {
        return Err(format!(
            "The selected filename would not preserve the {} sequence. Use a .{} filename.",
            doc.molecule.label(),
            format.extension_for(doc)
        ));
    }
    Ok(canonical)
}

fn worker<T: Send + 'static>(
    task: impl FnOnce() -> Result<T, String> + Send + 'static,
) -> Receiver<Result<T, String>> {
    let (tx, rx) = mpsc::channel();
    std::thread::spawn(move || {
        let result =
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(task)).unwrap_or_else(|_| {
                Err("The operation failed unexpectedly. Existing documents are unchanged.".into())
            });
        let _ = tx.send(result);
    });
    rx
}
fn receive<T>(receiver: &Option<Receiver<Result<T, String>>>) -> Option<Result<T, String>> {
    match receiver.as_ref()?.try_recv() {
        Ok(value) => Some(value),
        Err(TryRecvError::Empty) => None,
        Err(TryRecvError::Disconnected) => Some(Err(
            "The file worker ended unexpectedly. Your document is still open.".into(),
        )),
    }
}
fn resolved(path: &Path) -> PathBuf {
    if let Ok(path) = std::fs::canonicalize(path) {
        return path;
    }
    // Resolve parent symlinks for new destinations too, so aliases cannot bypass
    // open-document and in-flight save collision checks.
    if let (Some(parent), Some(name)) = (path.parent(), path.file_name())
        && let Ok(parent) = std::fs::canonicalize(parent)
    {
        return parent.join(name);
    }
    if path.is_absolute() {
        path.into()
    } else {
        std::env::current_dir().unwrap_or_default().join(path)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{
        sync::atomic::{AtomicU64, Ordering},
        time::{Duration, Instant},
    };
    static NEXT_TEMP: AtomicU64 = AtomicU64::new(0);

    #[test]
    fn restriction_filters_keep_recognition_counts_distinct_from_lengths() {
        assert!(EnzymeFilter::UniqueSix.includes(1, 6));
        assert!(!EnzymeFilter::UniqueSix.includes(2, 6));
        assert!(!EnzymeFilter::UniqueSix.includes(1, 4));
        assert!(EnzymeFilter::Unique.includes(1, 4));
        assert!(EnzymeFilter::OneOrTwo.includes(2, 4));
        assert!(!EnzymeFilter::All.includes(0, 6));
    }

    #[test]
    fn sequence_edit_cancels_analysis_without_locking_document_or_reusing_old_results() {
        let mut model = Model::new();
        model.ensure_current_analysis(None, None);
        let tab = model.current_mut().unwrap();
        assert!(tab.analyzing());
        assert!(!tab.busy());
        let cancel = tab.analysis_job.as_ref().unwrap().cancel.clone();
        assert!(tab.apply(|doc| doc.replace(0, 1, "A")));
        assert!(cancel.load(Ordering::Relaxed));
        assert!(!tab.analyzing());
        assert!(tab.restriction.is_none());
        assert!(!tab.poll_analysis());
    }

    #[test]
    fn changing_tab_cancels_unfinished_analysis_and_map_controls_are_not_edits() {
        let mut model = Model::new();
        model.ensure_current_analysis(None, None);
        let cancel = model.tabs[0].analysis_job.as_ref().unwrap().cancel.clone();
        model.tabs[0].show_labels = false;
        model.tabs[0].map_zoom = 2.;
        model.tabs[0].enzyme_filter = EnzymeFilter::All;
        assert!(!model.tabs[0].dirty());
        model.new_sequence();
        model.ensure_current_analysis(None, None);
        assert!(cancel.load(Ordering::Relaxed));
        assert!(!model.tabs[0].analyzing());
        assert!(!model.current().unwrap().analyzing());
    }

    #[test]
    fn sequence_layer_requests_analysis_without_enabling_map_sites() {
        let mut model = Model::new();
        let tab = model.current_mut().unwrap();
        tab.view = View::Sequence;
        tab.show_sites = false;
        let id = tab.id;
        model.ensure_current_analysis(None, Some(id));
        let tab = model.current().unwrap();
        assert!(tab.analyzing());
        assert!(!tab.show_sites);
        assert!(!tab.dirty());
        let cancel = tab.analysis_job.as_ref().unwrap().cancel.clone();
        model.ensure_current_analysis(None, None);
        assert!(cancel.load(Ordering::Relaxed));
        assert!(!model.current().unwrap().analyzing());

        // A request for a different tab never starts the active tab's analysis.
        model.ensure_current_analysis(None, Some(id + 100));
        assert!(!model.current().unwrap().analyzing());
        model.current_mut().unwrap().view = View::History;
        model.ensure_current_analysis(None, Some(id));
        assert!(!model.current().unwrap().analyzing());
    }

    #[test]
    fn sequence_scroll_targets_are_base_coordinates_independent_of_wrapping() {
        let mut tab = Model::new().tabs.remove(0);
        tab.select(
            vec![Segment {
                start: 137,
                end: 149,
            }],
            None,
            true,
        );
        assert_eq!(tab.scroll_to.take(), Some(137));
        assert_eq!(tab.range_start, "138");
        assert_eq!(tab.range_end, "149");
        tab.select(
            vec![Segment {
                start: 299,
                end: 300,
            }],
            None,
            false,
        );
        assert_eq!(tab.scroll_to, None);
    }

    #[test]
    fn enzyme_inspector_refreshes_after_edits_outside_map_and_cancels_when_closed() {
        for view in [View::Sequence, View::History] {
            let mut model = Model::new();
            let tab = model.current_mut().unwrap();
            tab.view = view;
            let id = tab.id;
            model.ensure_current_analysis(None, None);
            assert!(!model.current().unwrap().analyzing());

            model.ensure_current_analysis(Some(id), None);
            assert!(model.current().unwrap().analyzing());
            let deadline = Instant::now() + Duration::from_secs(10);
            while model.current().unwrap().analyzing() {
                model.poll();
                assert!(Instant::now() < deadline, "restriction analysis timed out");
                std::thread::sleep(Duration::from_millis(1));
            }
            let tab = model.current_mut().unwrap();
            assert!(tab.restriction.is_some());
            assert!(tab.apply(|doc| {
                let replacement = if doc.sequence.starts_with('A') {
                    "T"
                } else {
                    "A"
                };
                doc.replace(0, 1, replacement)
            }));
            assert!(tab.restriction.is_none());
            assert!(!tab.analyzing());

            model.ensure_current_analysis(Some(id), None);
            let tab = model.current().unwrap();
            assert!(tab.analyzing());
            let cancel = tab.analysis_job.as_ref().unwrap().cancel.clone();
            model.ensure_current_analysis(None, None);
            assert!(cancel.load(Ordering::Relaxed));
            assert!(!model.current().unwrap().analyzing());
            assert_eq!(model.current().unwrap().view, view);
        }
    }

    #[test]
    fn split_view_keeps_analysis_and_sequence_navigation_active() {
        let mut model = Model::new();
        model.current_mut().unwrap().view = View::Split;
        model.ensure_current_analysis(None, None);
        let tab = model.current_mut().unwrap();
        assert!(tab.analyzing());
        let cancel = tab.analysis_job.as_ref().unwrap().cancel.clone();

        // A map feature can reveal a wrapped sequence without leaving split view.
        let parts = vec![
            Segment {
                start: 4900,
                end: 5000,
            },
            Segment { start: 0, end: 120 },
        ];
        tab.select(parts.clone(), None, true);
        tab.show_sequence();
        assert_eq!(tab.view, View::Split);
        assert_eq!(tab.selected, parts);
        assert_eq!(tab.scroll_to.take(), Some(4900));
        assert!(!tab.dirty());

        // Selecting bases in the lower pane uses the same selection as the map.
        tab.select(vec![Segment { start: 10, end: 20 }], None, false);
        assert_eq!(
            tab.selected_sequence(),
            &tab.editor.document.sequence[10..20]
        );
        assert!(tab.scroll_to.is_none());
        model.ensure_current_analysis(None, None);
        assert!(!cancel.load(Ordering::Relaxed));

        model.current_mut().unwrap().show_sites = false;
        model.ensure_current_analysis(None, None);
        assert!(cancel.load(Ordering::Relaxed));
        assert!(!model.current().unwrap().analyzing());
    }

    #[test]
    fn origin_selection_summary_keeps_both_parts_and_length() {
        let mut model = Model::new();
        let tab = model.current_mut().unwrap();
        tab.select(
            vec![
                Segment {
                    start: 4900,
                    end: 5000,
                },
                Segment { start: 0, end: 120 },
            ],
            None,
            true,
        );
        assert!(tab.selection_summary.contains("4901–5000, 1–120"));
        assert!(tab.selection_summary.contains("220 bp"));
        assert_eq!(tab.scroll_to, Some(4900));
    }
    struct Temp(PathBuf);
    impl Temp {
        fn new() -> Self {
            let path = std::env::temp_dir().join(format!(
                "popgene-model-{}-{}",
                std::process::id(),
                NEXT_TEMP.fetch_add(1, Ordering::Relaxed)
            ));
            std::fs::create_dir_all(&path).unwrap();
            Self(path)
        }
    }
    impl Drop for Temp {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }
    fn settle(model: &mut Model) {
        let deadline = Instant::now() + Duration::from_secs(10);
        while model.is_busy() {
            model.poll();
            assert!(Instant::now() < deadline, "file operation timed out");
            std::thread::sleep(Duration::from_millis(5));
        }
        model.poll();
    }
    fn two_tabs() -> Model {
        let mut model = Model::new();
        model.add_document(
            Document::new("second", "ACGTACGT").unwrap(),
            Some(PathBuf::from("/example/second.gb")),
            None,
        );
        model
    }
    #[test]
    fn creation_draft_is_independent_and_background_close_targets_it() {
        let mut model = Model::new();
        let original = model.tabs[0].editor.document.clone();
        let id = model.new_sequence().unwrap();
        assert_eq!(model.tabs[1].title(), "New Sequence");
        assert!(!model.tabs[1].dirty());
        model.tabs[1].draft_dirty = true;
        model.active = 0;
        model.request_close(Close::Tab(id));
        assert_eq!(model.active, 1);
        assert_eq!(model.close_target(), Some(1));
        assert_eq!(model.tabs[0].editor.document, original);
        model.cancel_close();
        assert!(model.tabs[1].draft_dirty);
    }
    #[test]
    fn newly_created_document_is_dirty_without_an_edit() {
        let mut model = Model::new();
        let id = model.new_sequence().unwrap();
        let doc = Document::new_typed("protein", "MKWVTF*X", Molecule::Protein).unwrap();
        assert!(model.create_in_active(doc));
        let tab = &model.tabs[1];
        assert_eq!(tab.id, id);
        assert!(!tab.draft);
        assert!(tab.dirty());
        assert_eq!(tab.editor.document, tab.saved);
        assert_eq!(tab.view, View::Sequence);
        model.request_close(Close::Tab(id));
        assert_eq!(model.close_target(), Some(1));
    }
    #[test]
    fn blank_chooser_closes_without_prompt() {
        let mut model = Model::new();
        let id = model.new_sequence().unwrap();
        model.request_close(Close::Tab(id));
        assert_eq!(model.tabs.len(), 1);
        assert!(model.pending_close.is_none());
    }
    #[test]
    fn opened_document_replaces_only_untouched_chooser() {
        let mut model = Model::new();
        model.new_sequence();
        model.add_document(Document::new("opened", "ACGT").unwrap(), None, None);
        assert_eq!(model.tabs.len(), 2);
        assert!(!model.current().unwrap().draft);
        model.new_sequence();
        model.current_mut().unwrap().draft_dirty = true;
        model.add_document(Document::new("another", "ACGT").unwrap(), None, None);
        assert_eq!(model.tabs.len(), 4);
        assert!(model.tabs[2].draft_dirty);
    }
    #[test]
    fn tabs_keep_independent_edits_selection_and_undo() {
        let mut model = two_tabs();
        let original = model.tabs[0].editor.document.clone();
        assert!(model.tabs[1].apply(|d| d.replace(0, 1, "T")));
        model.tabs[1].select(vec![Segment { start: 2, end: 4 }], None, true);
        assert!(model.tabs[1].dirty());
        assert!(!model.tabs[0].dirty());
        assert_eq!(model.tabs[0].editor.document, original);
        assert!(model.tabs[0].selected.is_empty());
        assert!(!model.tabs[0].editor.can_undo());
        model.tabs[1].undo();
        assert!(!model.tabs[1].dirty());
        model.tabs[1].redo();
        assert!(model.tabs[1].dirty());
    }
    #[test]
    fn closing_clean_tab_keeps_other_dirty_document() {
        let mut model = two_tabs();
        model.tabs[1].apply(|d| d.replace(0, 1, "T"));
        model.request_close(Close::Tab(model.tabs[0].id));
        assert_eq!(model.tabs.len(), 1);
        assert!(model.tabs[0].dirty());
        assert_eq!(&*model.tabs[0].editor.document.sequence, "TCGTACGT");
    }
    #[test]
    fn quit_visits_every_dirty_tab_and_cancellation_preserves_all() {
        let mut model = two_tabs();
        for tab in &mut model.tabs {
            tab.apply(|d| {
                d.description = "edited".into();
                Ok(())
            });
        }
        model.request_close(Close::Quit);
        assert_eq!(model.close_target(), Some(0));
        model.discard_close();
        assert_eq!(model.close_target(), Some(1));
        model.cancel_close();
        assert_eq!(model.tabs.len(), 2);
        assert!(model.tabs.iter().all(Tab::dirty));
        assert!(!model.allow_close);
    }
    #[test]
    fn repeated_open_activates_existing_tab() {
        let mut model = two_tabs();
        model.active = 0;
        model.open_paths(vec![PathBuf::from("/example/second.gb")]);
        assert_eq!(model.tabs.len(), 2);
        assert_eq!(model.active, 1);
        assert!(!model.is_busy());
    }
    #[test]
    fn open_accepts_a_project_folder_and_file_in_one_selection() {
        let temp = Temp::new();
        let folder = temp.0.join("project");
        let nested = folder.join("sequences");
        std::fs::create_dir_all(&nested).unwrap();
        let file = nested.join("one.fasta");
        std::fs::write(&file, ">one\nACGT\n").unwrap();
        let mut model = Model::new();
        model.tabs[0].apply(|d| {
            d.description = "Unsaved work".into();
            Ok(())
        });

        model.open_paths(vec![folder.clone(), file.clone(), file.clone()]);
        settle(&mut model);

        let project = model.project.as_ref().unwrap();
        assert_eq!(project.root, resolved(&folder));
        assert_eq!(project.supported_files(), vec![resolved(&file)]);
        assert_eq!(model.tabs.len(), 2);
        assert_eq!(
            model.current().unwrap().editor.document.sequence.as_ref(),
            "ACGT"
        );
        assert_eq!(model.tabs[0].editor.document.description, "Unsaved work");
        assert!(model.tabs[0].dirty());
        assert!(model.issues.is_empty());
    }
    #[test]
    fn close_confirmation_targets_requested_document() {
        let mut model = two_tabs();
        model.tabs[0].apply(|d| {
            d.description = "edited".into();
            Ok(())
        });
        model.active = 1;
        model.request_close(Close::Tab(model.tabs[0].id));
        assert_eq!(model.active, 0);
        assert_eq!(model.close_target(), Some(0));
    }
    #[test]
    fn concurrent_saves_cannot_choose_same_destination() {
        let temp = Temp::new();
        let path = temp.0.join("one.gb");
        let mut model = two_tabs();
        assert!(model.save_to(0, path.clone(), SaveFormat::GenBank));
        assert!(!model.save_to(1, path.clone(), SaveFormat::GenBank));
        settle(&mut model);
        assert!(!model.save_to(1, path.clone(), SaveFormat::GenBank));
        assert_eq!(
            formats::open(&path).unwrap().name,
            model.tabs[0].editor.document.name
        );
        assert!(model.tabs[1].error);
    }
    #[test]
    fn failed_save_keeps_dirty_tab_and_close_confirmation() {
        let temp = Temp::new();
        let path = temp.0.join("one.gb");
        let mut model = two_tabs();
        assert!(model.save_to(1, path.clone(), SaveFormat::GenBank));
        settle(&mut model);
        model.tabs[1].apply(|d| d.replace(0, 1, "T"));
        let id = model.tabs[1].id;
        model.request_close(Close::Tab(id));
        std::fs::write(&path, b"external modification").unwrap();
        assert!(model.save_and_close(path.clone(), SaveFormat::GenBank));
        settle(&mut model);
        assert_eq!(model.tabs.len(), 2);
        assert!(model.tabs[1].dirty());
        assert!(model.tabs[1].error);
        assert_eq!(model.close_target(), Some(1));
        assert_eq!(std::fs::read(path).unwrap(), b"external modification");
    }
    #[test]
    fn successful_save_finishes_pending_close() {
        let temp = Temp::new();
        let mut model = two_tabs();
        model.tabs[1].apply(|d| d.replace(0, 1, "T"));
        model.request_close(Close::Tab(model.tabs[1].id));
        let path = temp.0.join("saved.gb");
        assert!(model.save_and_close(path.clone(), SaveFormat::GenBank));
        settle(&mut model);
        assert_eq!(model.tabs.len(), 1);
        assert!(model.pending_close.is_none());
        assert_eq!(&*formats::open(&path).unwrap().sequence, "TCGTACGT");
    }
    #[test]
    fn imported_standard_file_saves_in_place_and_reopens_cleanly() {
        let temp = Temp::new();
        let path = temp.0.join("annotated.gb");
        std::fs::write(&path, SAMPLE).unwrap();
        let mut model = Model::new();
        model.open_paths(vec![path.clone()]);
        settle(&mut model);
        assert_eq!(model.tabs[1].path.as_ref(), Some(&resolved(&path)));
        assert_eq!(model.tabs[1].save_format, Some(SaveFormat::GenBank));
        assert!(model.tabs[1].disk_hash.is_some());
        model.tabs[1].apply(|d| {
            d.description = "Edited description".into();
            Ok(())
        });
        assert!(model.save_to(1, path.clone(), SaveFormat::GenBank));
        settle(&mut model);
        assert!(!model.tabs[1].dirty());
        assert_eq!(model.tabs[1].editor.document, formats::open(&path).unwrap());
        assert_eq!(
            model.tabs[1].editor.document.description,
            "Edited description"
        );
        assert!(std::fs::read(path).unwrap().starts_with(b"LOCUS"));
    }
    #[test]
    fn misleading_extension_requires_save_as_without_changing_file_representation() {
        let temp = Temp::new();
        let path = temp.0.join("genbank-disguised.dna");
        std::fs::write(&path, SAMPLE).unwrap();
        let mut model = Model::new();
        model.open_paths(vec![path.clone()]);
        settle(&mut model);
        assert!(model.tabs[1].path.is_none());
        assert!(model.tabs[1].save_format.is_none());
        assert_eq!(
            SaveFormat::recommended(&model.tabs[1].editor.document),
            SaveFormat::GenBank
        );
        assert_eq!(std::fs::read(path).unwrap(), SAMPLE);
    }
    #[test]
    fn ordinary_save_keeps_selection_and_does_not_add_undo_for_source_bookkeeping() {
        let temp = Temp::new();
        let path = temp.0.join("dna.gb");
        let mut model = Model::new();
        model.add_document(Document::new("dna", "ACGT").unwrap(), None, None);
        model.tabs[1].select(vec![Segment { start: 1, end: 3 }], None, false);
        assert!(model.save_to(1, path, SaveFormat::GenBank));
        settle(&mut model);
        assert!(!model.tabs[1].editor.can_undo());
        assert_eq!(model.tabs[1].selected, [Segment { start: 1, end: 3 }]);
    }
    #[test]
    fn save_as_projects_to_disk_but_undo_and_original_recovery_keep_richer_input() {
        let temp = Temp::new();
        let path = temp.0.join("converted.fasta");
        let recovered = temp.0.join("original.gb");
        let mut model = Model::new();
        let original = model.tabs[0].editor.document.clone();
        assert!(!original.features.is_empty());
        assert!(model.save_to(0, path.clone(), SaveFormat::Fasta));
        settle(&mut model);
        assert!(!model.tabs[0].dirty());
        assert!(model.tabs[0].editor.document.features.is_empty());
        assert!(!model.tabs[0].editor.document.circular);
        assert_eq!(model.tabs[0].editor.document, formats::open(&path).unwrap());
        assert_eq!(model.tabs[0].path.as_ref(), Some(&resolved(&path)));
        assert_eq!(model.tabs[0].save_format, Some(SaveFormat::Fasta));
        assert!(model.export_original_to(0, recovered.clone()));
        settle(&mut model);
        assert_eq!(std::fs::read(recovered).unwrap(), SAMPLE);
        model.tabs[0].undo();
        assert_eq!(model.tabs[0].editor.document, original);
        assert!(model.tabs[0].dirty());
    }
    #[test]
    fn private_format_is_read_only_and_cannot_complete_save_and_close() {
        let temp = Temp::new();
        let path = temp.0.join("saved.popgene");
        let legacy = temp.0.join("legacy.popgene");
        let document = Document::new("Legacy synthetic fixture", "ACGTACGT").unwrap();
        {
            let db = rusqlite::Connection::open(&legacy).unwrap();
            db.execute_batch("PRAGMA application_id=1347375175; PRAGMA user_version=3; CREATE TABLE document (id INTEGER PRIMARY KEY, payload TEXT NOT NULL, original BLOB);").unwrap();
            db.execute(
                "INSERT INTO document VALUES (1, ?1, NULL)",
                [serde_json::to_string(&document).unwrap()],
            )
            .unwrap();
        }
        let mut model = Model::new();
        model.open_paths(vec![legacy]);
        settle(&mut model);
        assert_eq!(model.tabs.len(), 2);
        assert_eq!(model.tabs[1].editor.document.sequence.as_ref(), "ACGTACGT");
        assert!(model.tabs[1].path.is_none());
        assert!(model.tabs[1].save_format.is_none());
        model.tabs[1].apply(|d| {
            d.description = "unsaved".into();
            Ok(())
        });
        model.request_close(Close::Tab(model.tabs[1].id));
        assert!(!model.save_and_close(path.clone(), SaveFormat::GenBank));
        assert!(!path.exists());
        assert!(model.tabs[1].dirty());
        assert_eq!(model.close_target(), Some(1));
    }
    #[test]
    fn protein_fasta_extension_cannot_silently_change_molecule_type() {
        let temp = Temp::new();
        let mut model = Model::new();
        model.add_document(
            Document::new_typed("protein", "ACG", Molecule::Protein).unwrap(),
            None,
            None,
        );
        let misleading = temp.0.join("protein.fasta");
        assert!(model.save_to(1, misleading.clone(), SaveFormat::Fasta));
        settle(&mut model);
        assert!(model.tabs[1].error);
        assert!(!misleading.exists());
        let path = temp.0.join("protein.faa");
        assert!(model.save_to(1, path.clone(), SaveFormat::Fasta));
        settle(&mut model);
        assert!(!model.tabs[1].error);
        assert_eq!(formats::open(&path).unwrap().molecule, Molecule::Protein);
    }
    #[test]
    fn malformed_file_is_isolated_and_queue_deduplicates() {
        let temp = Temp::new();
        let good = temp.0.join("good.fasta");
        let bad = temp.0.join("bad.dna");
        std::fs::write(&good, ">good\nACGT\n").unwrap();
        std::fs::write(&bad, [9, 0, 0]).unwrap();
        let mut model = Model::new();
        model.open_paths(vec![good.clone(), bad, good]);
        settle(&mut model);
        assert_eq!(model.tabs.len(), 2);
        assert_eq!(model.issues.len(), 1);
        assert_eq!(model.tabs[1].editor.document.sequence.as_ref(), "ACGT");
    }
    #[test]
    fn origin_search_keeps_separate_segments_and_exact_copy() {
        let mut tab = Tab::from_document(1, Document::new("circle", "ACGTTGCA").unwrap());
        tab.editor.document.circular = true;
        tab.find("CAAC".into());
        assert_eq!(
            tab.selected,
            [Segment { start: 6, end: 8 }, Segment { start: 0, end: 2 }]
        );
        assert_eq!(tab.selected_sequence(), "CAAC");
    }
    #[test]
    fn source_file_cannot_be_overwritten_by_export() {
        let temp = Temp::new();
        let path = temp.0.join("source.gb");
        std::fs::write(&path, b"original data").unwrap();
        let mut model = Model::new();
        model.tabs[0].source_path = Some(resolved(&path));
        assert!(!model.export_to(0, path.clone(), SaveFormat::GenBank));
        assert_eq!(std::fs::read(path).unwrap(), b"original data");
    }
    #[test]
    fn invalid_ranges_and_edits_do_not_change_selection_or_document() {
        let mut tab = Tab::from_document(1, Document::new("DNA", "ACGT").unwrap());
        tab.select(vec![Segment { start: 1, end: 3 }], None, false);
        assert!(!tab.select_range("0", "100"));
        assert_eq!(tab.selected, [Segment { start: 1, end: 3 }]);
        assert!(!tab.apply(|d| d.replace(0, 2, "Q")));
        assert_eq!(tab.editor.document.sequence.as_ref(), "ACGT");
        assert!(!tab.dirty());
    }
}
