//! Public file formats only: no private document container or opaque application payloads.
use super::{Error, Result, abif, digest, embl, fasta_bytes, genbank_bytes, import, snapgene};
use popgene_core::{Document, Molecule};
use sha2::{Digest, Sha256};
use std::{io::Write, path::Path};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SaveFormat {
    SnapGene,
    GenBank,
    Fasta,
    Embl,
    Abif,
}
impl SaveFormat {
    pub const ALL: [Self; 5] = [
        Self::SnapGene,
        Self::GenBank,
        Self::Fasta,
        Self::Embl,
        Self::Abif,
    ];
    pub fn label(self) -> &'static str {
        match self {
            Self::SnapGene => "SnapGene DNA",
            Self::GenBank => "GenBank / GenPept",
            Self::Fasta => "FASTA",
            Self::Embl => "EMBL",
            Self::Abif => "ABIF trace (original)",
        }
    }
    pub fn extension(self) -> &'static str {
        match self {
            Self::SnapGene => "dna",
            Self::GenBank => "gb",
            Self::Fasta => "fasta",
            Self::Embl => "embl",
            Self::Abif => "ab1",
        }
    }
    pub fn extension_for(self, doc: &Document) -> &'static str {
        match (self, doc.molecule) {
            (Self::Fasta, Molecule::Protein) => "faa",
            (Self::Fasta, Molecule::Rna) => "frna",
            (Self::GenBank, Molecule::Protein) => "gp",
            _ => self.extension(),
        }
    }
    pub fn from_path(path: &Path) -> Option<Self> {
        let ext = path.extension()?.to_str()?.to_ascii_lowercase();
        match ext.as_str() {
            "dna" => Some(Self::SnapGene),
            "gb" | "gbk" | "genbank" | "gbff" | "gbank" | "gbs" | "gen" | "gp" | "genpept"
            | "genpep" | "gpff" | "ape" => Some(Self::GenBank),
            "fa" | "fas" | "fasta" | "fna" | "ffn" | "faa" | "frna" | "pep" => Some(Self::Fasta),
            "emb" | "embl" | "ebl" | "dat" => Some(Self::Embl),
            "ab1" | "abi" | "abd" | "abf" => Some(Self::Abif),
            _ => None,
        }
    }
    pub fn supports(self, doc: &Document) -> bool {
        if doc.placeholder || doc.sequence.is_empty() {
            return false;
        }
        match self {
            Self::SnapGene => doc.molecule == Molecule::Dna,
            Self::Embl => doc.molecule != Molecule::Protein,
            Self::Abif => original_matches(doc, "ABIF"),
            _ => true,
        }
    }
    pub fn recommended(doc: &Document) -> Self {
        if let Some(source) = &doc.source {
            let format = match source.format.as_str() {
                "SnapGene" => Some(Self::SnapGene),
                "GenBank" => Some(Self::GenBank),
                "FASTA" => Some(Self::Fasta),
                "EMBL" => Some(Self::Embl),
                "ABIF" => Some(Self::Abif),
                _ => None,
            };
            if let Some(format) = format
                && format.supports(doc)
            {
                return format;
            }
        }
        Self::GenBank
    }
}

fn equivalent(a: &Document, b: &Document) -> bool {
    a.name == b.name
        && a.molecule == b.molecule
        && a.single_stranded == b.single_stranded
        && a.placeholder == b.placeholder
        && a.sequence == b.sequence
        && a.circular == b.circular
        && a.description == b.description
        && a.accession == b.accession
        && a.comments == b.comments
        && a.features == b.features
        && a.trace == b.trace
}
fn original_matches(doc: &Document, format: &str) -> bool {
    let Some(source) = &doc.source else {
        return false;
    };
    if source.format != format || source.bytes.is_empty() {
        return false;
    }
    let parsed = if format == "ABIF" {
        abif::read(&source.file_name, &source.bytes)
    } else {
        import(&source.file_name, &source.bytes)
    };
    parsed.is_ok_and(|original| equivalent(doc, &original))
}
fn source_format(format: SaveFormat) -> &'static str {
    match format {
        SaveFormat::SnapGene => "SnapGene",
        SaveFormat::GenBank => "GenBank",
        SaveFormat::Fasta => "FASTA",
        SaveFormat::Embl => "EMBL",
        SaveFormat::Abif => "ABIF",
    }
}

/// Encode a documented third-party format. Unchanged same-format inputs retain their original bytes.
pub fn encode(doc: &Document, format: SaveFormat) -> Result<Vec<u8>> {
    doc.validate()?;
    if !format.supports(doc) {
        return Err(Error::Message(if doc.placeholder {
            "Add sequence before saving this empty placeholder.".into()
        } else if format == SaveFormat::Abif {
            "ABIF writing is limited to an unchanged original trace. Save sequence/annotations as GenBank or FASTA instead.".into()
        } else {
            format!("{} cannot store this molecule type.", format.label())
        }));
    }
    if format != SaveFormat::SnapGene && original_matches(doc, source_format(format)) {
        return Ok(doc.source.as_ref().unwrap().bytes.to_vec());
    }
    match format {
        SaveFormat::SnapGene => snapgene::write(doc),
        SaveFormat::GenBank => genbank_bytes(doc),
        SaveFormat::Fasta => Ok(fasta_bytes(doc)),
        SaveFormat::Embl => embl::write(doc),
        SaveFormat::Abif => unreachable!("supports verifies the unchanged original"),
    }
}

/// Present these messages before writing; they describe information absent from the saved document.
pub fn save_warnings(doc: &Document, format: SaveFormat) -> Vec<String> {
    if !format.supports(doc) {
        return vec![if doc.placeholder {
            "Add sequence before saving this empty placeholder.".into()
        } else {
            format!(
                "{} cannot save this document in its current state.",
                format.label()
            )
        }];
    }
    if format == SaveFormat::SnapGene {
        return snapgene::write_warnings(doc);
    }
    if original_matches(doc, source_format(format)) {
        return vec![];
    }
    let mut messages = vec![];
    if format == SaveFormat::Fasta {
        if doc.name.contains(['\r', '\n']) {
            messages.push("Line breaks in the record name are replaced by spaces.".into());
        }
        if !doc.features.is_empty() {
            messages.push("FASTA omits all annotations and primer binding sites.".into());
        }
        if doc.circular || doc.single_stranded != (doc.molecule == Molecule::Rna) {
            messages.push(
                "FASTA omits topology and strand-count information; it reopens as linear.".into(),
            );
        }
        if !doc.description.is_empty() || doc.accession.is_some() || !doc.comments.is_empty() {
            messages.push("FASTA omits the description, accession and comments; only the record name and sequence are saved.".into());
        }
        if doc.molecule != Molecule::Dna {
            messages.push(format!("Use .{} to retain the molecule type on reopening; FASTA has no universal molecule-type field.", format.extension_for(doc)));
        }
    }
    if matches!(format, SaveFormat::GenBank | SaveFormat::Embl) {
        let normalized = |text: &str| text.split_whitespace().collect::<Vec<_>>().join(" ");
        if normalized(&doc.description) != doc.description
            || doc.comments.iter().any(|s| normalized(s) != *s)
            || doc
                .features
                .iter()
                .flat_map(|f| &f.qualifiers)
                .any(|(_, value)| value.as_ref().is_some_and(|s| normalized(s) != *s))
        {
            messages.push("Whitespace and line breaks in descriptions, comments and qualifier text are normalized on reopening.".into());
        }
        let normalized = safe_name(&doc.name);
        if normalized != doc.name {
            messages.push(format!(
                "The record name will be saved as {normalized:?} to meet file-format naming rules."
            ));
        }
        if doc.features.iter().any(|f| f.edited) {
            messages.push("Annotations changed by sequence edits carry a review note; the in-memory review flag is not stored.".into());
        }
    }
    if doc.trace.is_some() {
        messages.push("Trace signal channels, peak locations and quality scores are omitted. Keep the original ABIF file for the chromatogram.".into());
    }
    if let Some(source) = &doc.source {
        if source.format == "GenBank"
            && format == SaveFormat::GenBank
            && import(&source.file_name, &source.bytes)
                .is_ok_and(|original| original.sequence != doc.sequence)
        {
            messages.push("Source references and unmodeled GenBank header fields are retained; review any coordinates or identifiers affected by sequence edits.".into());
        }
        if source.format == "SnapGene" {
            messages.push("SnapGene history, alignments, end chemistry, display settings, unbound primers and unmodeled metadata are omitted. Keep the original .dna file for those records.".into());
        } else if format == SaveFormat::Embl
            || (source.format == "EMBL" && format != SaveFormat::Embl)
        {
            messages.push("Only the sequence, topology, supported features, description, accession and comments are written; other source header fields and references are omitted.".into());
        }
    }
    messages
}

pub(super) fn safe_name(name: &str) -> String {
    let name: String = name
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '_' {
                c
            } else {
                '_'
            }
        })
        .take(40)
        .collect();
    if name.is_empty() {
        "Untitled".into()
    } else {
        name
    }
}

/// Atomic replacement with conflict checks before serialization and immediately before rename.
pub fn save(
    path: &Path,
    doc: &Document,
    format: SaveFormat,
    expected: Option<[u8; 32]>,
) -> Result<[u8; 32]> {
    if SaveFormat::from_path(path) != Some(format) {
        return Err(Error::Message(format!(
            "Use a .{} filename for {}. PopGene does not write a private file format.",
            format.extension_for(doc),
            format.label()
        )));
    }
    if format == SaveFormat::Fasta {
        let ext = path
            .extension()
            .and_then(|e| e.to_str())
            .unwrap_or("")
            .to_ascii_lowercase();
        if (doc.molecule == Molecule::Protein && !matches!(ext.as_str(), "faa" | "pep"))
            || (doc.molecule == Molecule::Rna && ext != "frna")
        {
            return Err(Error::Message(format!(
                "Use .{} for this molecule type, or choose GenBank / GenPept.",
                format.extension_for(doc)
            )));
        }
    }
    let check = || -> Result<()> {
        if let Some(hash) = expected
            && digest(path)? != hash
        {
            return Err(Error::Message(
                "The file changed outside PopGene. Use Save As to keep both versions.".into(),
            ));
        }
        Ok(())
    };
    check()?;
    let bytes = encode(doc, format)?;
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let mut temp = tempfile::NamedTempFile::new_in(parent)?;
    temp.write_all(&bytes)?;
    temp.as_file().sync_all()?;
    check()?;
    temp.persist(path).map_err(|e| Error::Io(e.error))?;
    #[cfg(unix)]
    std::fs::File::open(parent)?.sync_all()?;
    Ok(Sha256::digest(bytes).into())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn no_private_extension_can_be_written() {
        assert_eq!(SaveFormat::from_path(Path::new("old.popgene")), None);
        let doc = Document::new("dna", "ACGT").unwrap();
        let tmp = tempfile::tempdir().unwrap();
        assert!(
            save(
                &tmp.path().join("wrong.popgene"),
                &doc,
                SaveFormat::GenBank,
                None
            )
            .is_err()
        );
        assert!(!tmp.path().join("wrong.popgene").exists());
    }
    #[test]
    fn fasta_is_plain_standard_text_with_explicit_losses() {
        let mut doc = Document::new("named", "ACGT").unwrap();
        doc.circular = true;
        doc.comments.push("Comment".into());
        assert_eq!(encode(&doc, SaveFormat::Fasta).unwrap(), b">named\nACGT\n");
        let warnings = save_warnings(&doc, SaveFormat::Fasta).join(" ");
        assert!(warnings.contains("topology"));
        assert!(warnings.contains("comments"));
    }
    #[test]
    fn unchanged_abif_copy_retains_every_byte_and_rejects_edits() {
        let bytes = include_bytes!("../../../validation/fixtures/synthetic-complete.ab1");
        let mut doc = import("trace.ab1", bytes).unwrap();
        assert_eq!(encode(&doc, SaveFormat::Abif).unwrap(), bytes);
        assert!(save_warnings(&doc, SaveFormat::Abif).is_empty());
        doc.description.push_str(" changed");
        assert!(!SaveFormat::Abif.supports(&doc));
        assert!(encode(&doc, SaveFormat::Abif).is_err());
    }
    #[test]
    fn genpept_roundtrip_retains_protein_and_stop_symbols() {
        let doc = Document::new_typed("protein", "ACGMKWVTF*XZ", Molecule::Protein).unwrap();
        let bytes = encode(&doc, SaveFormat::GenBank).unwrap();
        assert!(
            String::from_utf8_lossy(&bytes)
                .lines()
                .next()
                .unwrap()
                .contains(" aa ")
        );
        let read = import("protein.gp", &bytes).unwrap();
        assert_eq!(read.sequence, doc.sequence);
        assert_eq!(read.molecule, Molecule::Protein);
    }
    #[test]
    fn unsupported_empty_placeholder_has_no_save_format() {
        let mut doc = Document::new("empty", "A").unwrap();
        doc.sequence = "".into();
        doc.placeholder = true;
        for format in SaveFormat::ALL {
            assert!(!format.supports(&doc));
            assert!(encode(&doc, format).is_err());
        }
    }
}
