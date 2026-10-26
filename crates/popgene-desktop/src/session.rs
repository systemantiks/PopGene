//! Transactional application recovery; never writes to an imported sequence file.
//!
//! Source bytes live in checksummed SQLite blobs rather than serde's skipped
//! Source.bytes field. The writer holds an exclusive connection so a second
//! running application cannot replace the first application's recovery state.
use crate::{
    model::{EnzymeFilter, MAX_TABS, Model, Tab, View, validate_enzyme_names},
    sequence_view::SequenceOptions,
};
use popgene_core::{Document, Editor, Segment, Source};
use popgene_formats::SaveFormat;
use rusqlite::{Connection, OptionalExtension, params};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    collections::HashMap,
    path::{Path, PathBuf},
    sync::{Arc, mpsc},
    time::Duration,
};

const SCHEMA: u32 = 1;
const MAX_JSON: usize = 256 * 1024 * 1024;
const MAX_BLOBS: usize = 256 * 1024 * 1024;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Snapshot {
    schema: u32,
    tabs: Vec<SavedTab>,
    active_id: Option<u64>,
    project: Option<PathBuf>,
    pub options: HashMap<u64, SequenceOptions>,
    #[serde(default)]
    pub drafts: HashMap<u64, crate::forms::CreatorDraft>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
struct SavedTab {
    id: u64,
    document: Document,
    saved: Document,
    original_source: Option<Source>,
    unsaved: bool,
    draft: bool,
    draft_dirty: bool,
    path: Option<PathBuf>,
    source_path: Option<PathBuf>,
    disk_hash: Option<[u8; 32]>,
    view: View,
    map_linear: bool,
    show_labels: bool,
    show_features: bool,
    show_primers: bool,
    show_positions: bool,
    show_sites: bool,
    map_zoom: f32,
    enzyme_filter: EnzymeFilter,
    #[serde(default)]
    custom_enzyme_names: Option<std::collections::BTreeSet<String>>,
    nonredundant: bool,
    selected: Vec<Segment>,
    feature: Option<u64>,
    query: String,
    trace_scale: f32,
}

impl Snapshot {
    pub fn capture(model: &Model, mut options: HashMap<u64, SequenceOptions>) -> Self {
        let tabs: Vec<_> = model
            .tabs
            .iter()
            .filter_map(|t| {
                let discarded = model.allow_close && model.discarded.contains(&t.id);
                // Discarding a never-saved document means it must not reappear.
                if discarded && (t.unsaved || t.draft) {
                    return None;
                }
                Some(SavedTab {
                    id: t.id,
                    document: if discarded {
                        t.saved.clone()
                    } else {
                        t.editor.document.clone()
                    },
                    saved: t.saved.clone(),
                    original_source: t.original_source.clone(),
                    unsaved: t.unsaved,
                    draft: t.draft,
                    draft_dirty: t.draft_dirty,
                    path: t.path.clone(),
                    source_path: t.source_path.clone(),
                    disk_hash: t.disk_hash,
                    view: t.view,
                    map_linear: t.map_linear,
                    show_labels: t.show_labels,
                    show_features: t.show_features,
                    show_primers: t.show_primers,
                    show_positions: t.show_positions,
                    show_sites: t.show_sites,
                    map_zoom: t.map_zoom,
                    enzyme_filter: t.enzyme_filter,
                    custom_enzyme_names: t.custom_enzyme_names.clone(),
                    nonredundant: t.nonredundant,
                    selected: t.selected.clone(),
                    feature: t.feature,
                    query: t.query.clone(),
                    trace_scale: t.trace_scale,
                })
            })
            .collect();
        options.retain(|id, _| tabs.iter().any(|t| t.id == *id));
        Self {
            schema: SCHEMA,
            tabs,
            active_id: model.current().map(|t| t.id),
            project: model.project.as_ref().map(|p| p.root.clone()),
            options,
            drafts: HashMap::new(),
        }
    }

    fn validate(&self) -> Result<(), String> {
        if self.schema != SCHEMA {
            return Err("Recovery data uses an unsupported version.".into());
        }
        if self.tabs.len() > MAX_TABS {
            return Err("Recovery data exceeds the tab limit.".into());
        }
        let mut ids = std::collections::HashSet::new();
        for t in &self.tabs {
            if t.id == 0 || t.id == u64::MAX || !ids.insert(t.id) {
                return Err("Invalid recovery tab identity.".into());
            }
            t.document.validate().map_err(|e| e.to_string())?;
            t.saved.validate().map_err(|e| e.to_string())?;
            if let Some(names) = &t.custom_enzyme_names {
                validate_enzyme_names(names)?;
            }
            if !t.map_zoom.is_finite() || !t.trace_scale.is_finite() {
                return Err("Invalid recovery display scale.".into());
            }
        }
        Ok(())
    }

    pub fn restore(self, model: &mut Model) -> Result<HashMap<u64, SequenceOptions>, String> {
        self.validate()?;
        let mut tabs = Vec::with_capacity(self.tabs.len());
        for saved in self.tabs {
            let mut tab = Tab::from_document(saved.id, saved.document);
            tab.saved = saved.saved;
            tab.original_source = saved.original_source;
            tab.unsaved = saved.unsaved;
            tab.draft = saved.draft;
            tab.draft_dirty = saved.draft_dirty;
            tab.path = saved.path;
            tab.source_path = saved.source_path;
            tab.disk_hash = saved.disk_hash;
            tab.save_format = tab
                .path
                .as_deref()
                .and_then(SaveFormat::from_path)
                .filter(|f| f.supports(&tab.editor.document));
            tab.view = saved.view;
            tab.map_linear = saved.map_linear;
            tab.show_labels = saved.show_labels;
            tab.show_features = saved.show_features;
            tab.show_primers = saved.show_primers;
            tab.show_positions = saved.show_positions;
            tab.show_sites = saved.show_sites;
            tab.map_zoom = saved.map_zoom.clamp(0.25, 20.0);
            tab.enzyme_filter = saved.enzyme_filter;
            tab.custom_enzyme_names = saved.custom_enzyme_names;
            tab.nonredundant = saved.nonredundant;
            tab.trace_scale = saved.trace_scale.clamp(0.1, 20.0);
            tab.find(saved.query);
            let selection = saved
                .selected
                .into_iter()
                .filter(|s| s.start < s.end && s.end <= tab.editor.document.sequence.len())
                .collect();
            let feature = saved
                .feature
                .filter(|id| tab.editor.document.features.iter().any(|f| f.id == *id));
            tab.select(selection, feature, true);
            if matches!(tab.view, View::Trace | View::Parameters)
                && tab.editor.document.trace.is_none()
            {
                tab.view = View::Sequence;
            }
            // Never associate the old document with a new hash. A changed/missing
            // source must use Save As, preserving both the recovered and disk data.
            if let Some(path) = &tab.source_path {
                let unchanged = tab
                    .disk_hash
                    .is_some_and(|expected| popgene_formats::digest(path).ok() == Some(expected));
                if !unchanged {
                    tab.path = None;
                    tab.save_format = None;
                    tab.unsaved = true;
                    tab.message(
                        "Recovered snapshot · source changed or is unavailable; use Save As.",
                        false,
                    );
                } else {
                    tab.message("Session restored · recovery snapshot loaded", false);
                }
            } else {
                tab.message("Recovered unsaved document", false);
            }
            // Preserve the previous saved baseline as an undoable recovery action.
            if tab.dirty() && tab.saved != tab.editor.document {
                let recovered = tab.editor.document.clone();
                tab.editor = Editor::new(tab.saved.clone());
                tab.editor
                    .apply(|d| {
                        *d = recovered;
                        Ok(())
                    })
                    .map_err(|e| e.to_string())?;
            }
            tabs.push(tab);
        }
        model.active = self
            .active_id
            .and_then(|id| tabs.iter().position(|t| t.id == id))
            .unwrap_or(0);
        model.next_id = tabs.iter().map(|t| t.id).max().unwrap_or(0) + 1;
        model.tabs = tabs;
        model.status = "Session restored · automatic recovery enabled".into();
        if let Some(root) = self.project {
            match popgene_formats::project::scan(&root) {
                Ok(project) => model.project = Some(project),
                Err(e) => model
                    .issues
                    .push(format!("Restored project is unavailable: {e}")),
            }
        }
        Ok(self.options)
    }
}

enum Request {
    Save(Box<Snapshot>, Option<mpsc::Sender<Result<(), String>>>),
}

pub struct Session {
    sender: mpsc::SyncSender<Request>,
    results: mpsc::Receiver<Result<(), String>>,
    submitted: Option<Snapshot>,
}

impl Session {
    pub fn open(path: PathBuf) -> Result<(Self, Option<Snapshot>), String> {
        let (sender, requests) = mpsc::sync_channel(1);
        let (ready, initial) = mpsc::sync_channel(1);
        let (completed, results) = mpsc::channel();
        std::thread::spawn(move || {
            let mut store = match Database::open(&path) {
                Ok(store) => store,
                Err(e) => {
                    let _ = ready.send(Err(e));
                    return;
                }
            };
            match store.load() {
                Ok(snapshot) => {
                    if ready.send(Ok(snapshot)).is_err() {
                        return;
                    }
                }
                Err(e) => {
                    let _ = ready.send(Err(e));
                    return;
                }
            }
            while let Ok(Request::Save(snapshot, acknowledged)) = requests.recv() {
                let result = store.save(&snapshot);
                if let Some(ack) = acknowledged {
                    let _ = ack.send(result.clone());
                }
                let _ = completed.send(result);
            }
        });
        let snapshot = initial
            .recv_timeout(Duration::from_secs(10))
            .map_err(|e| format!("Recovery storage did not start: {e}"))??;
        Ok((
            Self {
                sender,
                results,
                submitted: snapshot.clone(),
            },
            snapshot,
        ))
    }

    pub fn autosave(&mut self, snapshot: Snapshot) {
        if self.submitted.as_ref() == Some(&snapshot) {
            return;
        }
        if self
            .sender
            .try_send(Request::Save(Box::new(snapshot.clone()), None))
            .is_ok()
        {
            self.submitted = Some(snapshot);
        }
    }

    pub fn poll(&mut self) -> Option<Result<(), String>> {
        let mut last = None;
        while let Ok(result) = self.results.try_recv() {
            if result.is_err() {
                self.submitted = None;
            }
            last = Some(result);
        }
        last
    }

    /// A normal quit must finish writing the final Save/Discard choices.
    pub fn flush(&mut self, snapshot: Snapshot) -> Result<(), String> {
        let (ack, result) = mpsc::channel();
        self.sender
            .send(Request::Save(Box::new(snapshot.clone()), Some(ack)))
            .map_err(|e| format!("Recovery writer is unavailable: {e}"))?;
        result
            .recv_timeout(Duration::from_secs(30))
            .map_err(|e| format!("Recovery write did not finish: {e}"))??;
        self.submitted = Some(snapshot);
        Ok(())
    }
}

pub fn default_path() -> Option<PathBuf> {
    if let Some(path) = std::env::var_os("POPGENE_SESSION_PATH") {
        return Some(path.into());
    }
    #[cfg(target_os = "macos")]
    {
        std::env::var_os("HOME")
            .map(|p| PathBuf::from(p).join("Library/Application Support/PopGene/recovery.sqlite"))
    }
    #[cfg(target_os = "windows")]
    {
        std::env::var_os("LOCALAPPDATA").map(|p| PathBuf::from(p).join("PopGene/recovery.sqlite"))
    }
    #[cfg(not(any(target_os = "macos", target_os = "windows")))]
    {
        std::env::var_os("XDG_STATE_HOME")
            .map(PathBuf::from)
            .or_else(|| std::env::var_os("HOME").map(|p| PathBuf::from(p).join(".local/state")))
            .map(|p| p.join("popgene/recovery.sqlite"))
    }
}

struct Database {
    connection: Connection,
}
impl Database {
    fn open(path: &Path) -> Result<Self, String> {
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
        }
        let connection = Connection::open(path).map_err(|e| e.to_string())?;
        connection
            .busy_timeout(Duration::from_millis(200))
            .map_err(|e| e.to_string())?;
        connection.execute_batch("PRAGMA locking_mode=EXCLUSIVE; PRAGMA journal_mode=DELETE; PRAGMA synchronous=FULL;
            BEGIN EXCLUSIVE;
            CREATE TABLE IF NOT EXISTS workspace (id INTEGER PRIMARY KEY CHECK(id=1), json BLOB NOT NULL);
            CREATE TABLE IF NOT EXISTS sources (tab INTEGER NOT NULL, slot INTEGER NOT NULL, hash BLOB NOT NULL, PRIMARY KEY(tab,slot));
            CREATE TABLE IF NOT EXISTS blobs (hash BLOB PRIMARY KEY, bytes BLOB NOT NULL);
            COMMIT;").map_err(|e| format!("Cannot open recovery storage (another app instance may be using it): {e}"))?;
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            std::fs::set_permissions(path, std::fs::Permissions::from_mode(0o600))
                .map_err(|e| e.to_string())?;
        }
        Ok(Self { connection })
    }
    fn load(&self) -> Result<Option<Snapshot>, String> {
        let length: Option<i64> = self
            .connection
            .query_row("SELECT length(json) FROM workspace WHERE id=1", [], |r| {
                r.get(0)
            })
            .optional()
            .map_err(|e| e.to_string())?;
        let Some(length) = length else {
            return Ok(None);
        };
        if length < 0 || length as usize > MAX_JSON {
            return Err("Recovery document exceeds the size limit.".into());
        }
        let total: i64 = self
            .connection
            .query_row(
                "SELECT coalesce(sum(length(bytes)),0) FROM blobs",
                [],
                |r| r.get(0),
            )
            .map_err(|e| e.to_string())?;
        if total < 0 || total as usize > MAX_BLOBS {
            return Err("Recovery sources exceed the size limit.".into());
        }
        let json: Vec<u8> = self
            .connection
            .query_row("SELECT json FROM workspace WHERE id=1", [], |r| r.get(0))
            .map_err(|e| e.to_string())?;
        let mut snapshot: Snapshot =
            serde_json::from_slice(&json).map_err(|e| format!("Invalid recovery document: {e}"))?;
        snapshot.validate()?;
        let mut blobs: HashMap<Vec<u8>, Arc<[u8]>> = HashMap::new();
        for tab in &mut snapshot.tabs {
            for (slot, source) in [
                &mut tab.document.source,
                &mut tab.saved.source,
                &mut tab.original_source,
            ]
            .into_iter()
            .enumerate()
            {
                let Some(source) = source else {
                    continue;
                };
                let hash: Vec<u8> = self
                    .connection
                    .query_row(
                        "SELECT hash FROM sources WHERE tab=?1 AND slot=?2",
                        params![tab.id as i64, slot],
                        |r| r.get(0),
                    )
                    .map_err(|e| format!("Recovery source is missing: {e}"))?;
                if !blobs.contains_key(&hash) {
                    let bytes: Vec<u8> = self
                        .connection
                        .query_row("SELECT bytes FROM blobs WHERE hash=?1", [&hash], |r| {
                            r.get(0)
                        })
                        .map_err(|e| format!("Recovery source is missing: {e}"))?;
                    if Sha256::digest(&bytes).as_slice() != hash {
                        return Err("Recovery source checksum mismatch.".into());
                    }
                    blobs.insert(hash.clone(), bytes.into());
                }
                source.bytes = blobs[&hash].clone();
            }
        }
        Ok(Some(snapshot))
    }
    fn save(&mut self, snapshot: &Snapshot) -> Result<(), String> {
        snapshot.validate()?;
        let json = serde_json::to_vec(snapshot).map_err(|e| e.to_string())?;
        if json.len() > MAX_JSON {
            return Err("Recovery document exceeds the size limit.".into());
        }
        let transaction = self.connection.transaction().map_err(|e| e.to_string())?;
        transaction
            .execute("DELETE FROM sources", [])
            .map_err(|e| e.to_string())?;
        let mut total = 0usize;
        let mut hashes = std::collections::HashSet::new();
        for tab in &snapshot.tabs {
            for (slot, source) in [
                &tab.document.source,
                &tab.saved.source,
                &tab.original_source,
            ]
            .into_iter()
            .enumerate()
            {
                let Some(source) = source else {
                    continue;
                };
                let hash = Sha256::digest(&source.bytes).to_vec();
                if hashes.insert(hash.clone()) {
                    total = total.saturating_add(source.bytes.len());
                    if total > MAX_BLOBS {
                        return Err("Recovery sources exceed the size limit.".into());
                    }
                    transaction
                        .execute(
                            "INSERT OR IGNORE INTO blobs(hash,bytes) VALUES(?1,?2)",
                            params![hash, source.bytes.as_ref()],
                        )
                        .map_err(|e| e.to_string())?;
                }
                transaction
                    .execute(
                        "INSERT INTO sources(tab,slot,hash) VALUES(?1,?2,?3)",
                        params![tab.id as i64, slot, hash],
                    )
                    .map_err(|e| e.to_string())?;
            }
        }
        transaction
            .execute(
                "DELETE FROM blobs WHERE hash NOT IN (SELECT hash FROM sources)",
                [],
            )
            .map_err(|e| e.to_string())?;
        transaction
            .execute(
                "INSERT OR REPLACE INTO workspace(id,json) VALUES(1,?1)",
                [json],
            )
            .map_err(|e| e.to_string())?;
        transaction.commit().map_err(|e| e.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn path(name: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "popgene-recovery-{}-{}-{name}.sqlite",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ))
    }
    fn lineage_source() -> Document {
        let source = Document::new("recovery lineage", "ATGGAATTCAAGTAA").unwrap();
        let bytes = popgene_formats::encode(&source, SaveFormat::GenBank).unwrap();
        let mut source = popgene_formats::import("recovery-lineage.gb", &bytes).unwrap();
        // User-created feature epochs and the source seed both need recovery.
        source
            .annotate("whole CDS".into(), Segment { start: 0, end: 15 }, false)
            .unwrap();
        source.features[0].kind = "CDS".into();
        source.features[0]
            .qualifiers
            .push(("codon_start".into(), Some("1".into())));
        source
    }
    fn lineage_fragments(source: &Document) -> Vec<Document> {
        use popgene_core::{restriction, workflow};
        let enzyme = restriction::ENZYMES
            .iter()
            .position(|enzyme| enzyme.name == "EcoRI")
            .unwrap();
        let conditions = restriction::DigestConditions::default();
        let digest = restriction::digest_document(source, &[enzyme], conditions).unwrap();
        workflow::digest_products(source, &digest.fragments, &["EcoRI".into()], conditions).unwrap()
    }
    fn assert_lineage_rejoins(left: &Document, right: &Document) {
        use popgene_core::{cloning, translation};
        let product = cloning::ligate(
            &[
                cloning::AssemblyInput {
                    document: left.clone(),
                    reverse: false,
                },
                cloning::AssemblyInput {
                    document: right.clone(),
                    reverse: false,
                },
            ],
            cloning::LigationOptions::default(),
            "recovered rejoin",
        )
        .unwrap();
        assert_eq!(product.sequence.as_ref(), "ATGGAATTCAAGTAA");
        assert_eq!(product.features.len(), 1);
        assert_eq!(product.features[0].location_label(), "1..15");
        assert_eq!(
            translation::translate_cds(&product, &product.features[0])
                .unwrap()
                .amino_acids
                .as_slice(),
            b"MEFK*"
        );
    }
    #[test]
    fn lineage_source_seed_and_fragment_evidence_survive_sqlite_restart() {
        let path = path("lineage-roundtrip");
        let source = lineage_source();
        let original_bytes = source.source.as_ref().unwrap().bytes.clone();
        let seed = source.annotation_lineage.instance_id.clone();
        let epochs = source.annotation_lineage.feature_epochs.clone();
        assert!(!epochs.is_empty());
        let fragments = lineage_fragments(&source);
        let mut model = Model::new();
        model.tabs = vec![
            Tab::from_document(1, source),
            Tab::from_document(2, fragments[0].clone()),
            Tab::from_document(3, fragments[1].clone()),
        ];
        model.active = 2;
        let snapshot = Snapshot::capture(&model, HashMap::new());
        {
            let mut db = Database::open(&path).unwrap();
            db.save(&snapshot).unwrap();
        }
        let db = Database::open(&path).unwrap();
        let loaded = db.load().unwrap().unwrap();
        let mut restored = Model::new();
        loaded.restore(&mut restored).unwrap();
        assert_eq!(restored.active, 2);
        assert_eq!(restored.tabs.len(), 3);
        assert_eq!(
            restored.tabs[0]
                .editor
                .document
                .annotation_lineage
                .instance_id,
            seed
        );
        assert_eq!(
            restored.tabs[0]
                .editor
                .document
                .annotation_lineage
                .feature_epochs,
            epochs
        );
        assert_eq!(
            restored.tabs[0].original_source.as_ref().unwrap().bytes,
            original_bytes
        );
        for (index, fragment) in fragments.iter().enumerate() {
            let recovered = &restored.tabs[index + 1].editor.document;
            assert_eq!(recovered.annotation_lineage, fragment.annotation_lineage);
            assert!(recovered.annotation_lineage.pieces.iter().all(|piece| {
                recovered
                    .annotation_lineage
                    .actionable_piece(recovered, piece)
            }));
        }
        assert_lineage_rejoins(
            &restored.tabs[1].editor.document,
            &restored.tabs[2].editor.document,
        );
        // The otherwise empty source registry's seed must survive too: new
        // digestion after recovery can join a fragment from before recovery.
        let redigested = lineage_fragments(&restored.tabs[0].editor.document);
        assert_lineage_rejoins(&restored.tabs[1].editor.document, &redigested[1]);
        drop(db);
        std::fs::remove_file(path).unwrap();
    }
    #[test]
    fn recovered_lineage_edit_undo_redo_restores_only_the_matching_baseline() {
        let path = path("lineage-undo");
        let fragments = lineage_fragments(&lineage_source());
        let mut model = Model::new();
        model.tabs = vec![Tab::from_document(1, fragments[1].clone())];
        assert!(model.tabs[0].apply(|document| document.replace(5, 6, "T")));
        let edited = model.tabs[0].editor.document.clone();
        assert_ne!(edited.sequence, fragments[1].sequence);
        assert!(edited.annotation_lineage.pieces.is_empty());
        {
            let mut db = Database::open(&path).unwrap();
            db.save(&Snapshot::capture(&model, HashMap::new())).unwrap();
        }
        let db = Database::open(&path).unwrap();
        let mut restored = Model::new();
        db.load().unwrap().unwrap().restore(&mut restored).unwrap();
        let tab = &mut restored.tabs[0];
        assert_eq!(tab.editor.document, edited);
        assert!(tab.dirty());
        assert!(tab.editor.can_undo());
        tab.undo();
        assert_eq!(tab.editor.document, fragments[1]);
        assert!(!tab.dirty());
        assert_lineage_rejoins(&fragments[0], &tab.editor.document);
        tab.redo();
        assert_eq!(tab.editor.document, edited);
        assert!(tab.dirty());
        assert!(tab.editor.document.annotation_lineage.pieces.is_empty());
        drop(db);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn saved_physical_fragment_files_recover_clean_and_restore_undo_redo_baselines() {
        use crate::model::tests::{
            assert_physical_fragment, assert_physical_phase_two_rejoin,
            physical_phase_two_fragments, settle,
        };
        for format in [SaveFormat::GenBank, SaveFormat::Embl] {
            let directory = path("physical-fragment-files").with_extension("files");
            std::fs::create_dir(&directory).unwrap();
            let directory = directory.canonicalize().unwrap();
            let recovery = directory.join("recovery.sqlite");
            let fragments = physical_phase_two_fragments();
            let mut model = Model::new();
            model.tabs = fragments
                .iter()
                .enumerate()
                .map(|(index, document)| Tab::from_document(index as u64 + 1, document.clone()))
                .collect();
            model.active = 1;
            let files: Vec<_> = (0..2)
                .map(|index| directory.join(format!("fragment-{index}.{}", format.extension())))
                .collect();
            for (index, file) in files.iter().enumerate() {
                assert!(model.save_to(index, file.clone(), format));
                settle(&mut model);
                assert!(!model.tabs[index].error, "{}", model.tabs[index].status);
                assert!(!model.tabs[index].dirty());
            }
            let file_bytes: Vec<_> = files
                .iter()
                .map(|file| std::fs::read(file).unwrap())
                .collect();
            {
                let mut database = Database::open(&recovery).unwrap();
                database
                    .save(&Snapshot::capture(&model, HashMap::new()))
                    .unwrap();
            }
            drop(model);
            let mut restored = Model::new();
            {
                let database = Database::open(&recovery).unwrap();
                database
                    .load()
                    .unwrap()
                    .unwrap()
                    .restore(&mut restored)
                    .unwrap();
            }
            assert_eq!(restored.tabs.len(), 2);
            assert_eq!(restored.active, 1);
            for (index, tab) in restored.tabs.iter().enumerate() {
                assert!(!tab.dirty());
                assert!(!tab.unsaved);
                assert_eq!(tab.path.as_ref(), Some(&files[index]));
                assert_eq!(tab.source_path.as_ref(), Some(&files[index]));
                assert_eq!(tab.save_format, Some(format));
                assert_eq!(tab.editor.document, tab.saved);
                assert_eq!(
                    tab.saved.source.as_ref().unwrap().bytes.as_ref(),
                    file_bytes[index]
                );
                assert_physical_fragment(&tab.editor.document, &fragments[index]);
            }
            assert_physical_phase_two_rejoin(
                &restored.tabs[0].editor.document,
                &restored.tabs[1].editor.document,
            );

            // Save a later edit to recovery only. The on-disk file and saved
            // canonical baseline must remain independently rejoinable.
            assert!(restored.tabs[1].apply(|document| document.replace(5, 6, "T")));
            let edited = restored.tabs[1].editor.document.clone();
            assert!(edited.annotation_lineage.pieces.is_empty());
            {
                let mut database = Database::open(&recovery).unwrap();
                database
                    .save(&Snapshot::capture(&restored, HashMap::new()))
                    .unwrap();
            }
            drop(restored);
            let mut restored = Model::new();
            {
                let database = Database::open(&recovery).unwrap();
                database
                    .load()
                    .unwrap()
                    .unwrap()
                    .restore(&mut restored)
                    .unwrap();
            }
            assert_eq!(restored.tabs[1].editor.document, edited);
            assert!(restored.tabs[1].dirty());
            assert!(restored.tabs[1].editor.can_undo());
            restored.tabs[1].undo();
            assert!(!restored.tabs[1].dirty());
            assert_physical_fragment(&restored.tabs[1].editor.document, &fragments[1]);
            assert_physical_phase_two_rejoin(
                &restored.tabs[0].editor.document,
                &restored.tabs[1].editor.document,
            );
            restored.tabs[1].redo();
            assert_eq!(restored.tabs[1].editor.document, edited);
            assert!(restored.tabs[1].dirty());
            assert!(
                restored.tabs[1]
                    .editor
                    .document
                    .annotation_lineage
                    .pieces
                    .is_empty()
            );
            for (index, file) in files.iter().enumerate() {
                assert_eq!(std::fs::read(file).unwrap(), file_bytes[index]);
            }
            let reopened: Vec<_> = files
                .iter()
                .map(|file| popgene_formats::open(file).unwrap())
                .collect();
            assert_physical_phase_two_rejoin(&reopened[0], &reopened[1]);
            std::fs::remove_dir_all(directory).unwrap();
        }
    }
    #[test]
    fn sqlite_snapshots_without_lineage_get_fresh_seeds_without_inventing_ancestry() {
        let path = path("lineage-legacy");
        let fragments = lineage_fragments(&lineage_source());
        let mut model = Model::new();
        model.tabs = vec![Tab::from_document(1, fragments[1].clone())];
        let snapshot = Snapshot::capture(&model, HashMap::new());
        {
            let mut db = Database::open(&path).unwrap();
            db.save(&snapshot).unwrap();
            let mut legacy = serde_json::to_value(&snapshot).unwrap();
            for tab in legacy["tabs"].as_array_mut().unwrap() {
                for field in ["document", "saved"] {
                    tab[field]
                        .as_object_mut()
                        .unwrap()
                        .remove("annotation_lineage");
                }
            }
            db.connection
                .execute(
                    "UPDATE workspace SET json=?1 WHERE id=1",
                    [serde_json::to_vec(&legacy).unwrap()],
                )
                .unwrap();
        }
        let db = Database::open(&path).unwrap();
        let mut restored = Model::new();
        db.load().unwrap().unwrap().restore(&mut restored).unwrap();
        let tab = &restored.tabs[0];
        assert!(tab.editor.document.annotation_lineage.is_empty());
        assert_ne!(
            tab.editor.document.annotation_lineage.instance_id,
            fragments[1].annotation_lineage.instance_id
        );
        assert_eq!(tab.editor.document.sequence, fragments[1].sequence);
        assert_eq!(tab.editor.document.duplex_ends, fragments[1].duplex_ends);
        assert!(!tab.dirty());
        use popgene_core::cloning;
        let product = cloning::ligate(
            &[
                cloning::AssemblyInput {
                    document: fragments[0].clone(),
                    reverse: false,
                },
                cloning::AssemblyInput {
                    document: tab.editor.document.clone(),
                    reverse: false,
                },
            ],
            cloning::LigationOptions::default(),
            "missing old evidence",
        )
        .unwrap();
        assert_eq!(product.sequence.as_ref(), "ATGGAATTCAAGTAA");
        assert_eq!(product.features.len(), 2);
        drop(db);
        std::fs::remove_file(path).unwrap();
    }
    #[test]
    fn unsaved_edits_and_original_bytes_survive_a_restart() {
        let path = path("roundtrip");
        let mut model = Model::new();
        let original = model.tabs[0]
            .editor
            .document
            .source
            .as_ref()
            .unwrap()
            .bytes
            .clone();
        model.tabs[0].apply(|d| d.replace(0, 3, "AAA"));
        model.tabs[0].select(vec![Segment { start: 2, end: 8 }], None, true);
        model.tabs[0].view = View::Split;
        let options = HashMap::from([(
            1,
            SequenceOptions {
                reading_frames: 63,
                ..Default::default()
            },
        )]);
        let snapshot = Snapshot::capture(&model, options.clone());
        {
            let mut db = Database::open(&path).unwrap();
            db.save(&snapshot).unwrap();
        }
        let db = Database::open(&path).unwrap();
        let restored = db.load().unwrap().unwrap();
        let mut result = Model::new();
        assert_eq!(restored.restore(&mut result).unwrap(), options);
        assert_eq!(
            result.tabs[0].editor.document,
            model.tabs[0].editor.document
        );
        assert_eq!(
            result.tabs[0].original_source.as_ref().unwrap().bytes,
            original
        );
        assert!(result.tabs[0].dirty());
        assert!(result.tabs[0].editor.can_undo());
        assert_eq!(result.tabs[0].view, View::Split);
        assert_eq!(result.tabs[0].selected, model.tabs[0].selected);
        drop(db);
        std::fs::remove_file(path).unwrap();
    }
    #[test]
    fn custom_enzyme_lists_recover_per_tab_including_explicit_empty_and_legacy_defaults() {
        let path = path("enzyme-lists");
        let mut model = Model::new();
        model.tabs[0]
            .set_enzyme_selection(Some(
                crate::model::parse_enzyme_names("EcoRI BamHI").unwrap(),
            ))
            .unwrap();
        model.tabs[0].enzyme_filter = EnzymeFilter::All;
        model
            .add_generated_documents(vec![Document::new("empty list", "GAATTC").unwrap()], "test")
            .unwrap();
        model.tabs[1]
            .set_enzyme_selection(Some(Default::default()))
            .unwrap();
        let snapshot = Snapshot::capture(&model, HashMap::new());
        {
            let mut db = Database::open(&path).unwrap();
            db.save(&snapshot).unwrap();
        }
        let db = Database::open(&path).unwrap();
        let saved = db.load().unwrap().unwrap();
        let mut recovered = Model::new();
        saved.restore(&mut recovered).unwrap();
        assert_eq!(
            recovered.tabs[0].custom_enzyme_names,
            model.tabs[0].custom_enzyme_names
        );
        assert_eq!(recovered.tabs[0].enzyme_filter, EnzymeFilter::All);
        assert_eq!(
            recovered.tabs[1].custom_enzyme_names,
            Some(Default::default())
        );
        assert!(!recovered.tabs[0].dirty());

        let mut legacy = serde_json::to_value(&snapshot).unwrap();
        for tab in legacy["tabs"].as_array_mut().unwrap() {
            tab.as_object_mut().unwrap().remove("custom_enzyme_names");
        }
        let legacy: Snapshot = serde_json::from_value(legacy).unwrap();
        assert!(legacy.tabs.iter().all(|t| t.custom_enzyme_names.is_none()));
        let mut invalid = snapshot;
        invalid.tabs[0].custom_enzyme_names =
            Some(std::collections::BTreeSet::from(["ecori".into()]));
        assert!(invalid.validate().unwrap_err().contains("ecori"));
        drop(db);
        std::fs::remove_file(path).unwrap();
    }
    #[test]
    #[ignore = "subprocess helper for abrupt termination test"]
    fn recovery_crash_writer_child() {
        let Some(path) = std::env::var_os("POPGENE_TEST_CRASH_DATABASE") else {
            return;
        };
        let path = PathBuf::from(path);
        let mut model = Model::new();
        model.tabs[0].apply(|d| d.replace(0, 3, "AAA"));
        model
            .add_generated_documents(lineage_fragments(&lineage_source()), "crash lineage")
            .unwrap();
        let mut db = Database::open(&path).unwrap();
        db.save(&Snapshot::capture(&model, HashMap::new())).unwrap();
        db.connection.execute_batch("BEGIN EXCLUSIVE; UPDATE workspace SET json=x'00'; DELETE FROM sources; DELETE FROM blobs;").unwrap();
        std::fs::write(path.with_extension("ready"), b"transaction in progress").unwrap();
        // Parent kills this actual process, without Rust or SQLite cleanup.
        std::thread::sleep(Duration::from_secs(30));
        std::process::exit(91);
    }
    #[test]
    fn process_crash_rolls_back_partial_recovery_and_retains_committed_source() {
        let path = path("crash");
        let ready = path.with_extension("ready");
        let mut child = std::process::Command::new(std::env::current_exe().unwrap())
            .args([
                "--exact",
                "session::tests::recovery_crash_writer_child",
                "--ignored",
            ])
            .env("POPGENE_TEST_CRASH_DATABASE", &path)
            .stdout(std::process::Stdio::null())
            .spawn()
            .unwrap();
        for _ in 0..1000 {
            if ready.exists() || child.try_wait().unwrap().is_some() {
                break;
            }
            std::thread::sleep(Duration::from_millis(10));
        }
        let reached_transaction = ready.exists();
        let _ = child.kill();
        let status = child.wait().unwrap();
        assert!(
            reached_transaction,
            "Child did not reach an open transaction: {status}"
        );
        assert!(!status.success());
        let db = Database::open(&path).unwrap();
        let snapshot = db.load().unwrap().unwrap();
        let mut restored = Model::new();
        snapshot.restore(&mut restored).unwrap();
        assert!(restored.tabs[0].editor.document.sequence.starts_with("AAA"));
        assert!(restored.tabs[0].dirty());
        let baseline = Model::new();
        assert_eq!(
            restored.tabs[0].original_source.as_ref().unwrap().bytes,
            baseline.tabs[0].original_source.as_ref().unwrap().bytes
        );
        assert_eq!(restored.tabs.len(), 3);
        assert_lineage_rejoins(
            &restored.tabs[1].editor.document,
            &restored.tabs[2].editor.document,
        );
        drop(db);
        std::fs::remove_file(path).unwrap();
        std::fs::remove_file(ready).unwrap();
    }
    #[test]
    fn external_changes_cannot_be_overwritten_by_a_recovered_snapshot() {
        let path = path("external").with_extension("gb");
        let mut model = Model::new();
        std::fs::write(&path, b"before").unwrap();
        let tab = &mut model.tabs[0];
        tab.path = Some(path.clone());
        tab.source_path = Some(path.clone());
        tab.disk_hash = Some(popgene_formats::digest(&path).unwrap());
        let snapshot = Snapshot::capture(&model, HashMap::new());
        std::fs::write(&path, b"external edit").unwrap();
        snapshot.restore(&mut model).unwrap();
        assert!(model.tabs[0].path.is_none());
        assert!(model.tabs[0].unsaved);
        assert_eq!(std::fs::read(&path).unwrap(), b"external edit");
        std::fs::remove_file(path).unwrap();
    }
    #[test]
    fn explicit_discard_is_respected_and_closed_tabs_do_not_return() {
        let mut model = Model::new();
        model.tabs[0].apply(|d| d.replace(0, 3, "AAA"));
        model.request_close(crate::model::Close::Quit);
        model.discard_close();
        let expected = model.tabs[0].saved.clone();
        let snapshot = Snapshot::capture(&model, HashMap::new());
        snapshot.restore(&mut model).unwrap();
        assert_eq!(model.tabs[0].editor.document, expected);
        assert!(!model.tabs[0].dirty());
        model.allow_close = false;
        model.request_close(crate::model::Close::Tab(model.tabs[0].id));
        let snapshot = Snapshot::capture(&model, HashMap::new());
        snapshot.restore(&mut model).unwrap();
        assert!(model.tabs.is_empty());
    }
    #[test]
    fn second_instance_cannot_replace_active_recovery_state() {
        let path = path("exclusive");
        let first = Database::open(&path).unwrap();
        assert!(Database::open(&path).is_err());
        drop(first);
        assert!(Database::open(&path).is_ok());
        std::fs::remove_file(path).unwrap();
    }
    #[test]
    fn corrupt_source_is_rejected_without_changing_recovery_state() {
        let path = path("corrupt");
        let mut db = Database::open(&path).unwrap();
        db.save(&Snapshot::capture(&Model::new(), HashMap::new()))
            .unwrap();
        db.connection
            .execute("UPDATE blobs SET bytes=x'00'", [])
            .unwrap();
        assert!(db.load().unwrap_err().contains("checksum"));
        assert!(db.load().is_err());
        drop(db);
        std::fs::remove_file(path).unwrap();
    }
}
