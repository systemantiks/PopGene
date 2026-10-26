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
    Scf,
    Ztr,
    Fastq,
}
impl SaveFormat {
    pub const ALL: [Self; 8] = [
        Self::SnapGene,
        Self::GenBank,
        Self::Fasta,
        Self::Embl,
        Self::Abif,
        Self::Scf,
        Self::Ztr,
        Self::Fastq,
    ];
    pub fn label(self) -> &'static str {
        match self {
            Self::SnapGene => "SnapGene DNA / RNA / Protein",
            Self::GenBank => "GenBank / GenPept",
            Self::Fasta => "FASTA",
            Self::Embl => "EMBL",
            Self::Abif => "ABIF trace",
            Self::Scf => "SCF trace",
            Self::Ztr => "ZTR trace",
            Self::Fastq => "FASTQ (Phred+33)",
        }
    }
    pub fn extension(self) -> &'static str {
        match self {
            Self::SnapGene => "dna",
            Self::GenBank => "gb",
            Self::Fasta => "fasta",
            Self::Embl => "embl",
            Self::Abif => "ab1",
            Self::Scf => "scf",
            Self::Ztr => "ztr",
            Self::Fastq => "fastq",
        }
    }
    pub fn extension_for(self, doc: &Document) -> &'static str {
        match (self, doc.molecule) {
            (Self::SnapGene, Molecule::Protein) => "prot",
            (Self::SnapGene, Molecule::Rna) => "rna",
            (Self::Fasta, Molecule::Protein) => "faa",
            (Self::Fasta, Molecule::Rna) => "frna",
            (Self::GenBank, Molecule::Protein) => "gp",
            _ => self.extension(),
        }
    }
    pub fn from_path(path: &Path) -> Option<Self> {
        let ext = path.extension()?.to_str()?.to_ascii_lowercase();
        match ext.as_str() {
            "dna" | "rna" | "prot" => Some(Self::SnapGene),
            "gb" | "gbk" | "genbank" | "gbff" | "gbank" | "gbs" | "gen" | "gp" | "genpept"
            | "genpep" | "gpff" | "ape" => Some(Self::GenBank),
            "fa" | "fas" | "fasta" | "fna" | "ffn" | "faa" | "frna" | "pep" => Some(Self::Fasta),
            "emb" | "embl" | "ebl" | "dat" => Some(Self::Embl),
            "ab1" | "abi" | "abd" | "abf" => Some(Self::Abif),
            "scf" => Some(Self::Scf),
            "ztr" => Some(Self::Ztr),
            "fastq" | "fq" => Some(Self::Fastq),
            _ => None,
        }
    }
    pub fn supports(self, doc: &Document) -> bool {
        if doc.placeholder || doc.sequence.is_empty() {
            return false;
        }
        match self {
            Self::Embl => doc.molecule != Molecule::Protein,
            Self::Scf => doc
                .trace
                .as_ref()
                .is_some_and(|trace| trace.qualities_known),
            Self::Ztr => doc.trace.is_some(),
            Self::Fastq => super::fastq::supports(doc),
            Self::Abif => {
                doc.trace.is_some()
                    && doc
                        .source
                        .as_ref()
                        .is_some_and(|s| s.format == "ABIF" && !s.bytes.is_empty())
            }
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
                "SCF" => Some(Self::Scf),
                "ZTR" => Some(Self::Ztr),
                "FASTQ" => Some(Self::Fastq),
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
        && a.methylation == b.methylation
        && a.placeholder == b.placeholder
        && a.sequence == b.sequence
        && a.sequence_qualities == b.sequence_qualities
        && a.opposite_sequence_qualities == b.opposite_sequence_qualities
        && a.sequence_colors == b.sequence_colors
        && a.circular == b.circular
        && a.description == b.description
        && a.accession == b.accession
        && a.comments == b.comments
        && a.features == b.features
        && a.references == b.references
        && a.duplex_ends == b.duplex_ends
        && a.duplex_annotations == b.duplex_annotations
        && a.annotation_lineage == b.annotation_lineage
        && a.trace == b.trace
}
pub(super) fn original_matches(doc: &Document, format: &str) -> bool {
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
pub(super) fn source_format(format: SaveFormat) -> &'static str {
    match format {
        SaveFormat::SnapGene => "SnapGene",
        SaveFormat::GenBank => "GenBank",
        SaveFormat::Fasta => "FASTA",
        SaveFormat::Embl => "EMBL",
        SaveFormat::Abif => "ABIF",
        SaveFormat::Scf => "SCF",
        SaveFormat::Ztr => "ZTR",
        SaveFormat::Fastq => "FASTQ",
    }
}

/// Encode a documented third-party format. Unchanged same-format inputs retain their original bytes.
pub fn encode(doc: &Document, format: SaveFormat) -> Result<Vec<u8>> {
    doc.validate()?;
    if !format.supports(doc) {
        return Err(Error::Message(if doc.placeholder {
            "Add sequence before saving this empty placeholder.".into()
        } else if format == SaveFormat::Fastq {
            "FASTQ requires a recorded quality for every nucleotide. Unknown confidence cannot be replaced by an invented Phred zero.".into()
        } else if format == SaveFormat::Abif {
            "ABIF writing requires an imported ABIF chromatogram. Use SCF for a portable trace or FASTQ for calls and qualities.".into()
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
        SaveFormat::Abif => abif::write(doc),
        SaveFormat::Scf => super::scf::write(doc),
        SaveFormat::Ztr => super::ztr::write(doc),
        SaveFormat::Fastq => super::fastq::write(doc),
    }
}

/// Present these messages before writing; they describe information absent from the saved document.
pub fn save_warnings(doc: &Document, format: SaveFormat) -> Vec<String> {
    let mut warnings = base_warnings(doc, format);
    warnings.extend(super::lineage::warnings(doc, format));
    warnings
}

fn base_warnings(doc: &Document, format: SaveFormat) -> Vec<String> {
    if !format.supports(doc) {
        return vec![if doc.placeholder {
            "Add sequence before saving this empty placeholder.".into()
        } else if format == SaveFormat::Fastq {
            "FASTQ requires a recorded quality for every nucleotide; unknown confidence cannot be encoded as Phred zero.".into()
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
    let (physical_preserved, mut messages) = super::duplex::preservation(doc, format);
    if doc.methylation.is_some() && !physical_preserved {
        messages.push("Recorded whole-template Dam/Dcm/EcoKI methylation flags are omitted by this format. These flags describe source assumptions, not measured site occupancy.".into());
    }
    if doc.opposite_sequence_qualities.is_some() && !physical_preserved {
        messages.push("Opposite-strand quality scores are omitted by this format. FASTQ/trace exports carry only the active top-strand calls and qualities.".into());
    }
    if !doc.sequence_colors.is_empty() && !physical_preserved {
        messages.push("Custom per-strand nucleotide colors are omitted by this format.".into());
    }
    if doc.duplex_ends.is_some() && !physical_preserved {
        messages.push("Physical duplex ends, overhang geometry, opposite-strand-only residues and terminal modifications are omitted; only the top-strand sequence is saved. Use native DNA format to preserve physical ends.".into());
    }
    if !doc.duplex_annotations.is_empty() && !physical_preserved {
        messages.push("Annotations and colors spanning opposite-strand-only residues require the full duplex extent and are omitted by this top-strand-only format.".into());
    }
    if matches!(
        format,
        SaveFormat::Scf | SaveFormat::Ztr | SaveFormat::Fastq
    ) {
        messages.push("This format saves the active calls and their qualities. Document annotations, references, topology, alternate call families and per-call manual-edit flags are omitted.".into());
        if format == SaveFormat::Ztr {
            messages.push("ZTR retains displayed A/C/G/T signals, active peaks, called confidence and comments. Unmodeled original chunks, other confidence channels and display crop are omitted; preserve the source file for these records.".into());
        }
        if format == SaveFormat::Scf {
            messages.push("SCF retains all displayed signal channels and peak coordinates. The display crop is not stored. Original SCF comments/private bytes are retained; ABIF instrument records are not portable to SCF.".into());
        }
        if format == SaveFormat::Fastq && doc.trace.is_some() {
            messages.push("FASTQ omits chromatogram signal, peak positions and instrument metadata. Keep a trace-format copy.".into());
        }
        if format == SaveFormat::Fastq && doc.molecule == Molecule::Rna {
            messages.push("FASTQ writes nucleotide U as T and reopens as DNA.".into());
        }
        return messages;
    }
    if doc.sequence_qualities.is_some() && !physical_preserved {
        messages.push(
            "Per-symbol sequence quality scores are omitted; use FASTQ to preserve them.".into(),
        );
    }

    if format == SaveFormat::Abif {
        messages.push("Current calls, peaks and qualities replace edited family PBAS/PCON/PLOC1. Other instrument records, analysed/raw signals and original basecaller family 2 retain their original bytes. The former edited family 1 remains in the source file's unused bytes but is no longer an active call family.".into());
        if doc
            .trace
            .as_ref()
            .is_some_and(|t| t.sample_range.is_some() || t.manual_calls.iter().any(|m| *m))
        {
            messages.push("ABIF stores the edited calls and zero qualities for manual edits. The display crop and per-call manual flags are not standard ABIF fields and will not reopen; all raw signal samples remain.".into());
        }
        if !doc.features.is_empty() || !doc.references.is_empty() {
            messages.push("ABIF does not store PopGene annotations or document references.".into());
        }
        messages.push("ABIF retains the original instrument metadata; changes to the document name, description, comments and accession are not written.".into());
        return messages;
    }
    if format == SaveFormat::Fasta {
        if doc.name.contains(['\r', '\n']) {
            messages.push("Line breaks in the record name are replaced by spaces.".into());
        }
        if !doc.features.is_empty() {
            messages.push("FASTA omits all annotations, primer binding sites, segment colors/names, line styles and cleavage markers.".into());
        }
        if !doc.references.is_empty() {
            messages.push("FASTA omits all bibliographic references.".into());
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
        let extended_genbank =
            format == SaveFormat::GenBank && super::genbank_display::uses_extension(doc);
        if extended_genbank {
            messages.push("GenBank display extensions use an Exported LOCUS and retain the document name in KEYWORDS. A clearly labeled software interoperability reference is appended only when the required final marker is absent; existing scientific references remain intact.".into());
        }
        if doc.references.iter().any(|r| !r.extra.is_empty()) {
            messages.push("Bibliographic attributes without dedicated target fields are preserved as readable reference remarks/comments; their original separate attribute structure is not retained.".into());
        }
        if doc
            .features
            .iter()
            .any(|f| f.parts.is_empty() && f.location.is_none())
        {
            messages.push("Unbound oligos have no INSDC feature location; this conversion is refused until they are explicitly removed or saved in native format.".into());
        }
        if doc.features.iter().any(|f| {
            f.display_color_rgb().is_some() || f.display_segments().iter().any(|s| s.name.is_some())
        }) {
            messages.push("Annotation display metadata remains in conventional qualifiers/notes. Readers that do not implement this extension may show default colors, arrows or segment names.".into());
        }
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
        if normalized != doc.name && !extended_genbank {
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
            messages.push("Native procedure history, embedded attachments, alignments, display settings and unmodeled native metadata are omitted from this target format. Structured references and supported feature qualifiers are written where the target supports them. Keep the native source for the omitted records.".into());
        } else if source.format == "EMBL" && format != SaveFormat::Embl {
            messages.push("Additional unmodeled EMBL header lines are omitted when converting formats; structured references are retained by GenBank.".into());
        } else if source.format == "GenBank" && format == SaveFormat::Embl {
            messages.push("Additional unmodeled GenBank header fields are omitted in EMBL conversion; structured references and feature locations are retained.".into());
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
    if format == SaveFormat::SnapGene
        && !path
            .extension()
            .and_then(|ext| ext.to_str())
            .is_some_and(|ext| ext.eq_ignore_ascii_case(format.extension_for(doc)))
    {
        return Err(Error::Message(format!(
            "Use .{} for this native molecule type.",
            format.extension_for(doc)
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
    fn methylation_assumptions_are_explicit_on_lossy_format_conversions() {
        let mut doc = Document::new("recorded flags", "ACGT").unwrap();
        assert!(
            save_warnings(&doc, SaveFormat::SnapGene)
                .iter()
                .any(|warning| warning.contains("methylation flags are unknown"))
        );
        let bytes = encode(&doc, SaveFormat::GenBank).unwrap();
        doc = import("flags.gb", &bytes).unwrap();
        assert!(save_warnings(&doc, SaveFormat::GenBank).is_empty());
        doc.methylation = Some(popgene_core::metadata::DnaMethylation::default());
        for format in [SaveFormat::GenBank, SaveFormat::Embl, SaveFormat::Fasta] {
            assert!(
                save_warnings(&doc, format)
                    .iter()
                    .any(|warning| warning.contains("methylation flags are omitted"))
            );
            assert_eq!(
                import("converted", &encode(&doc, format).unwrap())
                    .unwrap()
                    .methylation,
                None
            );
        }
        assert!(
            !save_warnings(&doc, SaveFormat::SnapGene)
                .iter()
                .any(|warning| warning.contains("methylation flags are unknown"))
        );
    }
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
    fn abif_copy_is_exact_and_edited_family_roundtrips_without_changing_other_records() {
        let bytes = include_bytes!("../../../validation/fixtures/synthetic-complete.ab1");
        let mut doc = import("trace.ab1", bytes).unwrap();
        assert_eq!(encode(&doc, SaveFormat::Abif).unwrap(), bytes);
        assert!(save_warnings(&doc, SaveFormat::Abif).is_empty());
        let old_trace = doc.trace.clone().unwrap();
        popgene_core::trace::replace_call(&mut doc, 1, "N").unwrap();
        assert!(SaveFormat::Abif.supports(&doc));
        let encoded = encode(&doc, SaveFormat::Abif).unwrap();
        let reopened = import("edited.ab1", &encoded).unwrap();
        assert_eq!(reopened.sequence, doc.sequence);
        let new_trace = reopened.trace.unwrap();
        assert_eq!(new_trace.qualities, vec![40, 0, 30, 20, 5]);
        assert_eq!(new_trace.channels, old_trace.channels);
        assert_eq!(new_trace.peaks, old_trace.peaks);
        assert_eq!(&encoded[..18], &bytes[..18]);
        assert_eq!(&encoded[30..bytes.len()], &bytes[30..]);
        assert!(
            save_warnings(&doc, SaveFormat::Abif)
                .join(" ")
                .contains("manual flags")
        );
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
