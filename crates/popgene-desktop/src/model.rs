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

#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
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

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
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

/// User display selection, separate from source records and biological edits.
pub fn parse_enzyme_names(input: &str) -> Result<BTreeSet<String>, String> {
    if input.len() > 64 * 1024 {
        return Err("The enzyme list exceeds 64 KiB.".into());
    }
    let names = input
        .split(|c: char| c.is_whitespace() || matches!(c, ',' | ';'))
        .filter(|name| !name.is_empty())
        .map(str::to_owned)
        .collect();
    validate_enzyme_names(&names)?;
    Ok(names)
}

pub fn validate_enzyme_names(names: &BTreeSet<String>) -> Result<(), String> {
    if names.len() > ENZYMES.len() {
        return Err(format!(
            "Select at most {} catalogue enzymes.",
            ENZYMES.len()
        ));
    }
    let unknown = names
        .iter()
        .filter(|name| !ENZYMES.iter().any(|enzyme| enzyme.name == name.as_str()))
        .take(9)
        .collect::<Vec<_>>();
    if !unknown.is_empty() {
        let listed = unknown
            .iter()
            .take(8)
            .map(|name| name.chars().take(60).collect::<String>())
            .collect::<Vec<_>>()
            .join(", ");
        return Err(format!(
            "Unknown catalogue enzyme{}: {listed}{}. Use exact names, including case (for example, EcoRI).",
            if unknown.len() == 1 { "" } else { "s" },
            if unknown.len() > 8 { ", …" } else { "" }
        ));
    }
    Ok(())
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
    /// None uses the catalogue; an explicitly empty set displays no enzymes.
    pub custom_enzyme_names: Option<BTreeSet<String>>,
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

pub struct PrimerEdit {
    pub feature_id: Option<u64>,
    pub name: String,
    pub sequence: String,
    pub bindings: Vec<Segment>,
    pub reverse: bool,
    pub reaction: popgene_core::primer::Reaction,
    pub tm_c: Option<f64>,
}
impl Tab {
    pub fn save_primers(&mut self, edits: &[PrimerEdit]) -> Result<Vec<u64>, String> {
        if self.busy() || self.draft {
            return Err("Finish the current operation first.".into());
        }
        let mut ids = vec![];
        let success = self.apply(|doc| {
            if doc.molecule == Molecule::Protein {
                return Err(Error::Operation(
                    "Primers require a nucleotide template.".into(),
                ));
            }
            for edit in edits {
                edit.reaction.validate().map_err(Error::Operation)?;
                let sequence = Molecule::Dna.normalize(&edit.sequence)?;
                if sequence.is_empty() || edit.name.trim().is_empty() {
                    return Err(Error::Operation(
                        "Enter a primer name and oligo sequence.".into(),
                    ));
                }
                let binding_length: usize = edit
                    .bindings
                    .iter()
                    .map(|p| p.end.saturating_sub(p.start))
                    .sum();
                if edit
                    .bindings
                    .iter()
                    .any(|p| p.start >= p.end || p.end > doc.sequence.len())
                    || binding_length > sequence.len()
                {
                    return Err(Error::Operation(
                        "Primer bindings must fit the template and not exceed the oligo length."
                            .into(),
                    ));
                }
                let id = match edit.feature_id {
                    Some(id) => {
                        if !doc
                            .features
                            .iter()
                            .any(|f| f.id == id && f.kind == "primer_bind")
                        {
                            return Err(Error::Operation(
                                "The primer being edited is no longer present.".into(),
                            ));
                        }
                        id
                    }
                    None => doc
                        .features
                        .iter()
                        .chain(&doc.duplex_annotations.features)
                        .map(|f| f.id)
                        .max()
                        .unwrap_or(0)
                        .checked_add(1)
                        .ok_or_else(|| {
                            Error::Operation("Annotation identifier limit reached.".into())
                        })?,
                };
                let mut feature = doc.features.iter().find(|f| f.id == id).cloned().unwrap_or(
                    popgene_core::Feature {
                        id,
                        label: String::new(),
                        kind: "primer_bind".into(),
                        parts: vec![],
                        location: None,
                        reverse: false,
                        qualifiers: vec![],
                        edited: false,
                    },
                );
                feature.label = edit.name.trim().to_string();
                feature.parts = edit.bindings.clone();
                feature.location = None;
                feature.reverse = edit.reverse;
                feature.edited = false;
                feature.qualifiers.retain(|(key, _)| {
                    !matches!(
                        key.as_str(),
                        "primer_sequence"
                            | "binding_meltingTemperature"
                            | "calculated_melting_temperature"
                            | "primer_reaction"
                            | "primer_calculation_method"
                    )
                });
                feature
                    .qualifiers
                    .push(("primer_sequence".into(), Some(sequence)));
                feature.qualifiers.push((
                    "primer_reaction".into(),
                    Some(
                        serde_json::to_string(&edit.reaction)
                            .map_err(|e| Error::Operation(e.to_string()))?,
                    ),
                ));
                if let Some(tm) = edit.tm_c.filter(|v| v.is_finite()) {
                    feature.qualifiers.push((
                        "calculated_melting_temperature".into(),
                        Some(format!("{tm:.6}")),
                    ));
                    feature.qualifiers.push((
                        "primer_calculation_method".into(),
                        Some(popgene_core::primer::METHOD.into()),
                    ));
                }
                if let Some(existing) = doc.features.iter_mut().find(|f| f.id == id) {
                    *existing = feature;
                } else {
                    doc.features.push(feature);
                }
                ids.push(id);
            }
            Ok(())
        });
        if success {
            self.message(
                format!(
                    "Saved {} primer{}",
                    ids.len(),
                    if ids.len() == 1 { "" } else { "s" }
                ),
                false,
            );
            Ok(ids)
        } else {
            Err(self.status.clone())
        }
    }
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
            custom_enzyme_names: None,
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
        } else if let Some(annotation) =
            feature.and_then(|id| self.editor.document.features.iter().find(|f| f.id == id))
        {
            self.range_start.clear();
            self.range_end.clear();
            if scroll
                && let Some(boundary) = annotation
                    .location
                    .as_ref()
                    .and_then(|l| l.local_boundaries().first().copied())
            {
                self.scroll_to = Some(boundary.min(len.saturating_sub(1)));
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
            feature
                .and_then(|id| self.editor.document.features.iter().find(|f| f.id == id))
                .map(|f| {
                    format!(
                        "{} · {} · no sequence bases selected",
                        f.label,
                        f.location_label()
                    )
                })
                .unwrap_or_default()
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
    pub fn enzyme_selection_label(&self) -> String {
        match &self.custom_enzyme_names {
            None => "All catalogue enzymes".into(),
            Some(names) if names.is_empty() => "Custom list: no enzymes".into(),
            Some(names) => format!("Custom list: {} enzymes", names.len()),
        }
    }
    pub fn set_enzyme_selection(&mut self, names: Option<BTreeSet<String>>) -> Result<(), String> {
        if let Some(names) = &names {
            validate_enzyme_names(names)?;
        }
        self.custom_enzyme_names = names;
        Ok(())
    }
    /// Shared by map/sequence sites, the enzyme table, and figure inputs.
    pub fn includes_enzyme(&self, index: usize, count: usize) -> bool {
        ENZYMES.get(index).is_some_and(|enzyme| {
            self.enzyme_filter
                .includes(count, enzyme.recognition_length())
                && self
                    .custom_enzyme_names
                    .as_ref()
                    .is_none_or(|names| names.contains(enzyme.name))
        })
    }
    pub fn displayed_enzyme_indices(&self) -> Vec<usize> {
        let Some(analysis) = &self.restriction else {
            return vec![];
        };
        analysis
            .counts
            .iter()
            .enumerate()
            .filter_map(|(index, &count)| self.includes_enzyme(index, count).then_some(index))
            .collect()
    }
    pub fn displayed_sites(&self) -> Vec<&RestrictionSite> {
        let Some(analysis) = &self.restriction else {
            return vec![];
        };
        analysis
            .sites
            .iter()
            .filter(|site| {
                self.includes_enzyme(site.enzyme_index, analysis.counts[site.enzyme_index])
            })
            .collect()
    }
    pub fn site_parts(&self, site: &RestrictionSite) -> Vec<Segment> {
        let n = self.editor.document.sequence.len();
        let end = site.start + site.recognition_length;
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
    File(Vec<Document>, PathBuf, Option<[u8; 32]>),
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
    pub(crate) next_id: u64,
    queue: VecDeque<PathBuf>,
    opening: Option<Receiver<Result<Loaded, String>>>,
    opening_path: Option<PathBuf>,
    pub(crate) discarded: BTreeSet<u64>,
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
                let docs =
                    formats::open_records(&path).map_err(|e| format!("{}: {e}", path.display()))?;
                let after = formats::digest(&path).map_err(|e| e.to_string())?;
                if before != after {
                    return Err(format!(
                        "{} changed while opening. Open it again to load the current version.",
                        path.display()
                    ));
                }
                let hash = Some(after);
                Ok(Loaded::File(docs, path, hash))
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
                Ok(Loaded::File(docs, path, hash)) => {
                    self.add_records(docs, path, hash);
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
                            | (SaveFormat::Scf, "SCF")
                            | (SaveFormat::Ztr, "ZTR")
                            | (SaveFormat::Fastq, "FASTQ")
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
    /// Generated products open atomically and are unsaved from the outset.
    pub fn add_generated_documents(
        &mut self,
        documents: Vec<Document>,
        status: &str,
    ) -> Result<(), String> {
        if documents.is_empty() {
            return Err("No products were selected".into());
        }
        let replace = self
            .current()
            .filter(|t| t.draft && !t.dirty())
            .map(|t| t.id);
        if self.tabs.len() + documents.len() - usize::from(replace.is_some()) > MAX_TABS {
            return Err(
                "The products would exceed the 64-tab limit. Select fewer products or close files."
                    .into(),
            );
        }
        let used: usize = self
            .tabs
            .iter()
            .map(|tab| {
                let source = tab.editor.document.source.as_ref();
                source.map_or(0, |s| s.bytes.len())
                    + tab.original_source.as_ref().map_or(0, |original| {
                        if source.is_some_and(|s| Arc::ptr_eq(&s.bytes, &original.bytes)) {
                            0
                        } else {
                            original.bytes.len()
                        }
                    })
            })
            .sum();
        let incoming: usize = documents
            .iter()
            .map(|doc| doc.source.as_ref().map_or(0, |s| s.bytes.len()))
            .sum();
        if used.saturating_add(incoming) > MAX_SOURCE_BYTES {
            return Err("These products exceed the retained-source memory limit. Close files before creating them.".into());
        }
        let count = u64::try_from(documents.len()).map_err(|_| "Too many products")?;
        let next = self
            .next_id
            .checked_add(count)
            .ok_or("Document identifier limit reached")?;
        for doc in &documents {
            doc.validate().map_err(|e| e.to_string())?;
        }
        let tabs = documents
            .into_iter()
            .enumerate()
            .map(|(i, doc)| {
                let mut tab = Tab::from_document(self.next_id + i as u64, doc);
                tab.unsaved = true;
                tab.message(status, false);
                tab
            })
            .collect::<Vec<_>>();
        if let Some(id) = replace {
            self.remove_tab(id);
        }
        self.tabs.extend(tabs);
        self.next_id = next;
        self.active = self.tabs.len() - 1;
        self.status = status.into();
        Ok(())
    }

    /// A multi-record source opens atomically as independent documents. Each
    /// record uses Save As; saving one must never replace its sibling records.
    pub fn add_records(
        &mut self,
        mut docs: Vec<Document>,
        path: PathBuf,
        hash: Option<[u8; 32]>,
    ) -> bool {
        if docs.len() == 1 {
            return self.add_document(docs.remove(0), Some(path), hash);
        }
        if docs.is_empty() {
            return false;
        }
        let path = resolved(&path);
        if let Some(i) = self
            .tabs
            .iter()
            .position(|t| t.source_path.as_ref() == Some(&path))
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
                current.map_or(0, |s| s.bytes.len())
                    + tab.original_source.as_ref().map_or(0, |s| {
                        if current.is_some_and(|c| Arc::ptr_eq(&c.bytes, &s.bytes)) {
                            0
                        } else {
                            s.bytes.len()
                        }
                    })
            })
            .sum();
        let incoming: usize = docs
            .iter()
            .map(|d| d.source.as_ref().map_or(0, |s| s.bytes.len()))
            .sum();
        if self.tabs.len() + docs.len() - usize::from(replace.is_some()) > MAX_TABS
            || used.saturating_add(incoming) > MAX_SOURCE_BYTES
        {
            self.issues.push(format!("{}: {} records would exceed the open-tab limit (64 tabs / 256 MiB source attachments). No records were opened; close tabs and try again.", path.display(), docs.len()));
            self.status = "Multi-record file exceeds workspace limits — see report".into();
            return false;
        }
        if let Some(error) = docs.iter().find_map(|d| d.validate().err()) {
            self.issues.push(format!("{}: {error}", path.display()));
            return false;
        }
        if let Some(id) = replace {
            self.remove_tab(id);
        }
        let first = self.tabs.len();
        let count = docs.len();
        for (index, doc) in docs.into_iter().enumerate() {
            let mut tab = Tab::from_document(self.next_id, doc);
            self.next_id += 1;
            tab.source_path = Some(path.clone());
            tab.disk_hash = hash;
            tab.message(
                format!(
                    "Record {} of {count} · use Save As to save this record separately",
                    index + 1
                ),
                false,
            );
            self.tabs.push(tab);
        }
        self.active = first;
        self.status = format!("Opened {count} records · original file retained");
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
    pub fn export_attachment_to(
        &mut self,
        index: usize,
        attachment_id: u32,
        path: PathBuf,
    ) -> bool {
        let Some(tab) = self.tabs.get(index).filter(|t| !t.draft && !t.busy()) else {
            return false;
        };
        let Some(source) = tab.original_source.clone() else {
            return false;
        };
        let path = resolved(&path);
        if self.collision(index, &path, true) {
            self.tabs[index].message(
                "Choose an attachment destination that is not open or being written.",
                true,
            );
            return false;
        }
        let tab = &mut self.tabs[index];
        tab.message("Extracting embedded attachment…", false);
        tab.start_job(path.clone(), move || {
            let bytes = formats::native_attachment_bytes(&source, attachment_id)
                .map_err(|e| e.to_string())?;
            formats::atomic_export(&path, &bytes).map_err(|e| e.to_string())?;
            Ok(Completed::Exported(path))
        });
        true
    }
    pub fn export_figure_to(
        &mut self,
        index: usize,
        path: PathBuf,
        svg: String,
        format: crate::figure::FigureFormat,
        scale: f32,
        page_height: Option<u32>,
    ) -> bool {
        if self
            .tabs
            .get(index)
            .is_none_or(|tab| tab.draft || tab.busy())
        {
            return false;
        }
        if path
            .extension()
            .and_then(|s| s.to_str())
            .is_none_or(|s| !s.eq_ignore_ascii_case(format.extension()))
        {
            self.tabs[index].message(
                format!(
                    "Use the .{} extension for this figure format.",
                    format.extension()
                ),
                true,
            );
            return false;
        }
        let path = resolved(&path);
        if self.collision(index, &path, true) {
            self.tabs[index].message(
                "Export figures to a destination that is not open or being written.",
                true,
            );
            return false;
        }
        let tab = &mut self.tabs[index];
        tab.message("Rendering publication figure…", false);
        tab.start_job(path.clone(), move || {
            let bytes = if page_height.is_some() {
                crate::figure::encode_with_pages(&svg, format, scale, page_height)?
            } else {
                crate::figure::encode(&svg, format, scale)?
            };
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
pub(crate) mod tests {
    use super::*;
    use std::{
        sync::atomic::{AtomicU64, Ordering},
        time::{Duration, Instant},
    };
    static NEXT_TEMP: AtomicU64 = AtomicU64::new(0);

    pub(crate) fn physical_phase_two_fragments() -> Vec<Document> {
        use popgene_core::{metadata, restriction, workflow};
        let mut source = Document::new("physical phase two", "ATGGAATTCAAGTAA").unwrap();
        source
            .annotate("coding α".into(), Segment { start: 0, end: 15 }, false)
            .unwrap();
        source.features[0].kind = "CDS".into();
        source.features[0].qualifiers.extend([
            ("codon_start".into(), Some("2".into())),
            ("citation".into(), Some("1".into())),
            ("snapgene_color".into(), Some("#336699".into())),
        ]);
        source.references.push(metadata::Reference {
            title: "Physical fragment integration α".into(),
            authors: Some("Synthetic, A.".into()),
            ..Default::default()
        });
        source.sequence_qualities = Some((0..15).map(|q| (q != 1).then_some(q)).collect());
        source.opposite_sequence_qualities =
            Some((30..45).map(|q| (q != 32).then_some(q)).collect());
        source.sequence_colors = vec![
            metadata::SequenceColor {
                range: Segment { start: 0, end: 3 },
                reverse: false,
                color_rgb: 0x112233,
            },
            metadata::SequenceColor {
                range: Segment { start: 9, end: 12 },
                reverse: true,
                color_rgb: 0x993366,
            },
        ];
        let enzyme = restriction::ENZYMES
            .iter()
            .position(|enzyme| enzyme.name == "EcoRI")
            .unwrap();
        let conditions = restriction::DigestConditions::default();
        let digest = restriction::digest_document(&source, &[enzyme], conditions).unwrap();
        let fragments =
            workflow::digest_products(&source, &digest.fragments, &["EcoRI".into()], conditions)
                .unwrap();
        assert_eq!(fragments[0].sequence.as_ref(), "ATGG");
        assert_eq!(
            fragments[0].duplex_ends.as_ref().unwrap().opposite_sequence,
            "AATTCCAT"
        );
        assert_eq!(fragments[1].sequence.as_ref(), "AATTCAAGTAA");
        assert_eq!(
            fragments[1].duplex_ends.as_ref().unwrap().opposite_sequence,
            "TTACTTG"
        );
        assert_eq!(
            fragments[0].duplex_annotations.features[0].qualifier("codon_start"),
            Some("2")
        );
        assert_eq!(fragments[1].features[0].qualifier("codon_start"), Some("1"));
        fragments
    }

    pub(crate) fn assert_physical_fragment(actual: &Document, original: &Document) {
        assert_eq!(actual.sequence, original.sequence);
        assert_eq!(
            actual.duplex_ends, original.duplex_ends,
            "{:?}",
            actual.warnings
        );
        assert_eq!(actual.sequence_qualities, original.sequence_qualities);
        assert_eq!(
            actual.opposite_sequence_qualities,
            original.opposite_sequence_qualities
        );
        assert_eq!(actual.sequence_colors, original.sequence_colors);
        assert_eq!(
            actual.duplex_annotations.sequence_colors,
            original.duplex_annotations.sequence_colors
        );
        assert_eq!(
            actual.annotation_lineage.origins,
            original.annotation_lineage.origins
        );
        assert_eq!(actual.annotation_lineage.pieces.len(), 1);
        assert!(
            actual
                .annotation_lineage
                .pieces
                .iter()
                .all(|piece| actual.annotation_lineage.actionable_piece(actual, piece))
        );
        let features: Vec<_> = actual
            .features
            .iter()
            .chain(&actual.duplex_annotations.features)
            .collect();
        let expected: Vec<_> = original
            .features
            .iter()
            .chain(&original.duplex_annotations.features)
            .collect();
        assert_eq!(features.len(), 1);
        assert_eq!(features[0].location_label(), expected[0].location_label());
        assert_eq!(
            features[0].qualifier("codon_start"),
            expected[0].qualifier("codon_start")
        );
        assert_eq!(features[0].display_color_rgb(), Some(0x336699));
        let citation: usize = features[0]
            .qualifier("citation")
            .unwrap()
            .trim()
            .trim_start_matches('[')
            .trim_end_matches(']')
            .parse()
            .unwrap();
        assert_eq!(
            actual.references[citation - 1].title,
            "Physical fragment integration α"
        );
        assert_eq!(
            actual.references[citation - 1].authors.as_deref(),
            Some("Synthetic, A.")
        );
    }

    pub(crate) fn assert_physical_phase_two_rejoin(left: &Document, right: &Document) {
        use popgene_core::{cloning, translation};
        let product = cloning::ligate(
            &[left, right].map(|document| cloning::AssemblyInput {
                document: document.clone(),
                reverse: false,
            }),
            cloning::LigationOptions::default(),
            "saved physical rejoin",
        )
        .unwrap();
        assert_eq!(product.sequence.as_ref(), "ATGGAATTCAAGTAA");
        assert_eq!(
            product.duplex_ends.as_ref().unwrap().opposite_sequence,
            "TTACTTGAATTCCAT"
        );
        assert_eq!(product.features.len(), 1);
        assert!(product.duplex_annotations.features.is_empty());
        let feature = &product.features[0];
        assert_eq!(feature.location_label(), "1..15");
        assert_eq!(feature.qualifier("codon_start"), Some("2"));
        assert_eq!(feature.display_color_rgb(), Some(0x336699));
        assert_eq!(
            translation::translate_cds(&product, feature)
                .unwrap()
                .amino_acids,
            b"WNSS"
        );
        let citation: usize = feature
            .qualifier("citation")
            .unwrap()
            .trim()
            .trim_start_matches('[')
            .trim_end_matches(']')
            .parse()
            .unwrap();
        assert_eq!(
            product.references[citation - 1].title,
            "Physical fragment integration α"
        );
        assert_eq!(
            product.sequence_qualities,
            Some((0..15).map(|q| (q != 1).then_some(q)).collect())
        );
        assert_eq!(
            product.opposite_sequence_qualities,
            Some((30..45).map(|q| (q != 32).then_some(q)).collect())
        );
        assert!(
            product
                .sequence_colors
                .iter()
                .any(|color| color.range == Segment { start: 0, end: 3 }
                    && !color.reverse
                    && color.color_rgb == 0x112233)
        );
        assert!(
            product
                .sequence_colors
                .iter()
                .any(|color| color.range == Segment { start: 9, end: 12 }
                    && color.reverse
                    && color.color_rgb == 0x993366)
        );
    }

    #[test]
    fn physical_fragment_saved_projection_reopen_and_undo_redo_preserve_reassembly() {
        for format in [SaveFormat::GenBank, SaveFormat::Embl] {
            let temp = Temp::new();
            let fragments = physical_phase_two_fragments();
            let mut model = Model::new();
            model.tabs = fragments
                .iter()
                .enumerate()
                .map(|(i, document)| Tab::from_document(i as u64 + 1, document.clone()))
                .collect();
            let mut reopened = Vec::new();
            for (index, original) in fragments.iter().enumerate() {
                let warnings = formats::save_warnings(original, format);
                assert!(
                    !warnings
                        .iter()
                        .any(|warning| warning.contains("supplement is omitted")),
                    "{warnings:?}"
                );
                let path = temp
                    .0
                    .join(format!("fragment-{index}.{}", format.extension()));
                assert!(model.save_to(index, path.clone(), format));
                settle(&mut model);
                let tab = &model.tabs[index];
                assert!(!tab.error, "{}: {}", format.label(), tab.status);
                assert!(!tab.dirty());
                assert!(!tab.unsaved);
                assert_eq!(tab.path.as_ref(), Some(&path.canonicalize().unwrap()));
                assert_eq!(tab.editor.document, tab.saved);
                let bytes = std::fs::read(&path).unwrap();
                assert_eq!(
                    tab.saved.source.as_ref().unwrap().bytes.as_ref(),
                    bytes.as_slice()
                );
                let from_disk = formats::open(&path).unwrap();
                assert_physical_fragment(&tab.editor.document, original);
                assert_physical_fragment(&from_disk, original);
                assert_eq!(from_disk.annotation_lineage, tab.saved.annotation_lineage);
                reopened.push(from_disk);
            }
            // Reassembly uses only separately decoded public files, not the
            // original source or the model's pre-save documents.
            assert_physical_phase_two_rejoin(&reopened[0], &reopened[1]);
            let canonical = model.tabs[1].editor.document.clone();
            assert!(model.tabs[1].apply(|document| document.replace(5, 6, "T")));
            let edited = model.tabs[1].editor.document.clone();
            assert!(model.tabs[1].dirty());
            assert!(edited.annotation_lineage.pieces.is_empty());
            model.tabs[1].undo();
            assert_eq!(model.tabs[1].editor.document, canonical);
            assert!(!model.tabs[1].dirty());
            assert_physical_phase_two_rejoin(&reopened[0], &model.tabs[1].editor.document);
            model.tabs[1].redo();
            assert_eq!(model.tabs[1].editor.document, edited);
            assert!(model.tabs[1].dirty());
            assert!(
                model.tabs[1]
                    .editor
                    .document
                    .annotation_lineage
                    .pieces
                    .is_empty()
            );
        }
    }

    fn lineage_document() -> Document {
        let mut source = Document::new("saved ancestry", "ATGGAATTCAAGTAA").unwrap();
        source
            .annotate("whole CDS".into(), Segment { start: 0, end: 15 }, false)
            .unwrap();
        source.features[0].kind = "CDS".into();
        source.features[0]
            .qualifiers
            .push(("codon_start".into(), Some("1".into())));
        // This is a fully paired, exact-location fixture: flat-file physical-end
        // loss and native fuzzy-location policy are separate codec constraints.
        popgene_core::workflow::extract_region(&source, 0, 15, "saved ancestry").unwrap()
    }

    #[test]
    fn lineage_supported_save_canonical_projection_and_reopen_remain_actionable() {
        for format in [SaveFormat::SnapGene, SaveFormat::GenBank, SaveFormat::Embl] {
            let temp = Temp::new();
            let path = temp.0.join(format!("ancestry.{}", format.extension()));
            let original = lineage_document();
            assert!(!original.annotation_lineage.pieces.is_empty());
            let origin_ids = original
                .annotation_lineage
                .origins
                .iter()
                .map(|origin| origin.id.clone())
                .collect::<Vec<_>>();
            let mut model = Model::new();
            model.tabs = vec![Tab::from_document(1, original.clone())];
            assert!(model.save_to(0, path.clone(), format));
            settle(&mut model);
            let tab = &model.tabs[0];
            assert!(!tab.error, "{}: {}", format.label(), tab.status);
            assert!(!tab.dirty());
            let canonical = &tab.editor.document;
            assert_eq!(
                canonical.annotation_lineage.instance_id,
                original.annotation_lineage.instance_id
            );
            assert_eq!(
                canonical
                    .annotation_lineage
                    .origins
                    .iter()
                    .map(|origin| origin.id.clone())
                    .collect::<Vec<_>>(),
                origin_ids
            );
            assert_eq!(canonical.annotation_lineage.pieces.len(), 1);
            assert!(canonical.annotation_lineage.pieces.iter().all(|piece| {
                canonical
                    .annotation_lineage
                    .actionable_piece(canonical, piece)
            }));
            let reopened = formats::open(&path).unwrap();
            assert_eq!(reopened.annotation_lineage, canonical.annotation_lineage);
            assert_eq!(tab.saved.annotation_lineage, canonical.annotation_lineage);
            assert_eq!(
                reopened.feature_sequence(&reopened.features[0]).unwrap(),
                "ATGGAATTCAAGTAA"
            );
        }
    }

    #[test]
    fn lineage_lossy_save_is_undoable_while_export_retains_live_evidence() {
        let temp = Temp::new();
        let export = temp.0.join("export.fasta");
        let saved = temp.0.join("saved.fasta");
        let original = lineage_document();
        assert!(
            formats::save_warnings(&original, SaveFormat::Fasta)
                .iter()
                .any(|warning| warning.contains("ancestry"))
        );
        let mut model = Model::new();
        model.tabs = vec![Tab::from_document(1, original.clone())];
        assert!(model.export_to(0, export.clone(), SaveFormat::Fasta));
        settle(&mut model);
        assert!(!model.tabs[0].error, "{}", model.tabs[0].status);
        assert_eq!(model.tabs[0].editor.document, original);
        assert!(
            formats::open(&export)
                .unwrap()
                .annotation_lineage
                .is_empty()
        );
        assert!(model.save_to(0, saved.clone(), SaveFormat::Fasta));
        settle(&mut model);
        let tab = &mut model.tabs[0];
        assert!(!tab.error, "{}", tab.status);
        assert!(!tab.dirty());
        assert!(tab.editor.document.annotation_lineage.is_empty());
        assert!(tab.saved.annotation_lineage.is_empty());
        assert!(formats::open(&saved).unwrap().annotation_lineage.is_empty());
        tab.undo();
        assert_eq!(tab.editor.document, original);
        assert!(tab.dirty());
        assert!(
            tab.editor
                .document
                .annotation_lineage
                .pieces
                .iter()
                .all(|piece| tab
                    .editor
                    .document
                    .annotation_lineage
                    .actionable_piece(&tab.editor.document, piece))
        );
        tab.redo();
        assert!(tab.editor.document.annotation_lineage.is_empty());
        assert!(!tab.dirty());
    }

    #[test]
    fn generated_products_open_atomically_and_are_unsaved() {
        let mut model = Model::new();
        let before = model
            .tabs
            .iter()
            .map(|t| (t.id, t.editor.document.clone()))
            .collect::<Vec<_>>();
        let first = Document::new("first product", "ACGT").unwrap();
        let mut invalid = Document::new("invalid product", "ACGT").unwrap();
        invalid.sequence = "?".into();
        assert!(
            model
                .add_generated_documents(vec![first.clone(), invalid], "created")
                .is_err()
        );
        assert_eq!(
            model
                .tabs
                .iter()
                .map(|t| (t.id, t.editor.document.clone()))
                .collect::<Vec<_>>(),
            before
        );
        model
            .add_generated_documents(
                vec![first, Document::new("second product", "TGCA").unwrap()],
                "created",
            )
            .unwrap();
        assert_eq!(model.tabs.len(), before.len() + 2);
        for tab in &model.tabs[before.len()..] {
            assert!(tab.dirty());
            assert!(tab.path.is_none());
            assert!(tab.original_source.is_none());
            assert_eq!(tab.status, "created");
        }
        assert_eq!(
            model
                .tabs
                .iter()
                .take(before.len())
                .map(|t| (t.id, t.editor.document.clone()))
                .collect::<Vec<_>>(),
            before
        );
    }

    fn primer_edit(name: &str, bindings: Vec<Segment>) -> PrimerEdit {
        PrimerEdit {
            feature_id: None,
            name: name.into(),
            sequence: "ACGTACGTACGTACGTACGT".into(),
            bindings,
            reverse: false,
            reaction: Default::default(),
            tm_c: Some(55.25),
        }
    }

    #[test]
    fn primer_pair_save_is_one_atomic_undoable_edit() {
        let mut model = Model::new();
        let tab = &mut model.tabs[0];
        let before = tab.editor.document.clone();
        let edits = [
            primer_edit("forward", vec![Segment { start: 0, end: 20 }]),
            primer_edit("reverse", vec![Segment { start: 20, end: 40 }]),
        ];
        let ids = tab.save_primers(&edits).unwrap();
        assert_eq!(ids.len(), 2);
        assert_ne!(ids[0], ids[1]);
        let feature = tab
            .editor
            .document
            .features
            .iter()
            .find(|f| f.id == ids[0])
            .unwrap();
        assert_eq!(
            feature.qualifier("calculated_melting_temperature"),
            Some("55.250000")
        );
        let reaction: popgene_core::primer::Reaction =
            serde_json::from_str(feature.qualifier("primer_reaction").unwrap()).unwrap();
        assert_eq!(reaction, Default::default());
        tab.undo();
        assert_eq!(tab.editor.document, before);
        let invalid = [
            primer_edit("valid", vec![]),
            primer_edit("invalid", vec![Segment { start: 0, end: 21 }]),
        ];
        assert!(tab.save_primers(&invalid).is_err());
        assert_eq!(tab.editor.document, before);
    }

    #[test]
    fn unbound_oligo_survives_template_edit_and_can_be_bound_later() {
        let mut model = Model::new();
        let tab = &mut model.tabs[0];
        let ids = tab.save_primers(&[primer_edit("Unbound", vec![])]).unwrap();
        assert!(tab.apply(|d| d.replace(0, 1, "A")));
        let feature = tab
            .editor
            .document
            .features
            .iter()
            .find(|f| f.id == ids[0])
            .unwrap();
        assert!(feature.parts.is_empty());
        assert_eq!(feature.location_label(), "Unbound oligo");
        let mut bound = primer_edit("Bound", vec![Segment { start: 1, end: 21 }]);
        bound.feature_id = Some(ids[0]);
        tab.save_primers(&[bound]).unwrap();
        assert_eq!(
            tab.editor
                .document
                .features
                .iter()
                .filter(|f| f.id == ids[0])
                .count(),
            1
        );
        assert_eq!(
            tab.editor
                .document
                .features
                .iter()
                .find(|f| f.id == ids[0])
                .unwrap()
                .parts,
            vec![Segment { start: 1, end: 21 }]
        );
    }

    #[test]
    fn multi_record_open_preserves_siblings_and_requires_independent_save() {
        let mut model = Model::new();
        let docs = formats::import_records("two.fa", b">first\nACGT\n>second\nTTAA\n").unwrap();
        let path = std::env::temp_dir().join("popgene-two-records.fa");
        assert!(model.add_records(docs, path.clone(), None));
        assert_eq!(model.tabs.len(), 3);
        assert_eq!(model.active, 1);
        for tab in &model.tabs[1..] {
            assert_eq!(tab.source_path.as_ref(), Some(&resolved(&path)));
            assert!(tab.path.is_none());
            assert!(tab.save_format.is_none());
        }
        model.tabs[1].apply(|d| d.replace(0, 1, "T"));
        assert_eq!(&*model.tabs[2].editor.document.sequence, "TTAA");
        let docs = formats::import_records("two.fa", b">first\nACGT\n>second\nTTAA\n").unwrap();
        assert!(model.add_records(docs, path, None));
        assert_eq!(model.tabs.len(), 3);
    }

    #[test]
    fn over_capacity_multi_record_import_has_no_partial_tabs() {
        let mut model = Model::new();
        let before = model.tabs[0].editor.document.clone();
        let docs = (0..MAX_TABS)
            .map(|i| Document::new(format!("record-{i}"), "ACGT").unwrap())
            .collect();
        assert!(!model.add_records(docs, PathBuf::from("many.fa"), None));
        assert_eq!(model.tabs.len(), 1);
        assert_eq!(model.tabs[0].editor.document, before);
        assert!(
            model
                .issues
                .last()
                .unwrap()
                .contains("No records were opened")
        );
    }

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
    fn custom_enzyme_names_require_canonical_catalogue_entries() {
        assert_eq!(
            parse_enzyme_names("EcoRI, BamHI\nEcoRI; BamHI").unwrap(),
            BTreeSet::from(["BamHI".into(), "EcoRI".into()])
        );
        assert!(parse_enzyme_names(" ; ,\n").unwrap().is_empty());
        let error = parse_enzyme_names("ecori, ImaginaryI").unwrap_err();
        assert!(error.contains("ecori") && error.contains("ImaginaryI"));
        assert!(error.contains("including case"));
        assert!(parse_enzyme_names(&"x".repeat(65537)).is_err());
        // A rejected replacement must not discard the previous valid selection.
        let mut tab = Tab::from_document(1, Document::new("list", "GAATTC").unwrap());
        let selected = Some(parse_enzyme_names("EcoRI").unwrap());
        tab.set_enzyme_selection(selected.clone()).unwrap();
        assert!(
            tab.set_enzyme_selection(Some(BTreeSet::from(["ecori".into()])))
                .is_err()
        );
        assert_eq!(tab.custom_enzyme_names, selected);
        assert!(!tab.dirty());
    }

    #[test]
    fn custom_enzyme_selection_intersects_count_filters_and_survives_edits() {
        let mut tab = Tab::from_document(1, Document::new("cuts", "GAATTCGGATCCGAATTC").unwrap());
        tab.set_enzyme_selection(Some(parse_enzyme_names("EcoRI BamHI").unwrap()))
            .unwrap();
        tab.restriction = Some(Arc::new(restriction::find_sites(
            &tab.editor.document.sequence,
            false,
        )));
        let names = |tab: &Tab| {
            tab.displayed_enzyme_indices()
                .into_iter()
                .map(|i| ENZYMES[i].name)
                .collect::<BTreeSet<_>>()
        };
        assert_eq!(names(&tab), BTreeSet::from(["BamHI"]));
        tab.enzyme_filter = EnzymeFilter::All;
        assert_eq!(names(&tab), BTreeSet::from(["BamHI", "EcoRI"]));
        assert_eq!(tab.displayed_sites().len(), 3);
        assert!(!tab.dirty());
        assert!(tab.apply(|doc| doc.replace(0, 6, "AAAAAA")));
        assert_eq!(tab.custom_enzyme_names.as_ref().unwrap().len(), 2);
        assert!(tab.restriction.is_none());
        tab.restriction = Some(Arc::new(restriction::find_sites(
            &tab.editor.document.sequence,
            false,
        )));
        tab.enzyme_filter = EnzymeFilter::UniqueSix;
        assert_eq!(names(&tab), BTreeSet::from(["BamHI", "EcoRI"]));
        assert_eq!(tab.displayed_sites().len(), 2);
        tab.set_enzyme_selection(Some(BTreeSet::new())).unwrap();
        assert!(tab.displayed_enzyme_indices().is_empty());
        assert!(tab.displayed_sites().is_empty());
        assert_eq!(tab.enzyme_selection_label(), "Custom list: no enzymes");
        tab.set_enzyme_selection(None).unwrap();
        assert!(tab.displayed_enzyme_indices().len() > 2);
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
    pub(crate) fn settle(model: &mut Model) {
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
