//! Bounded file adapters. Unsupported records fail explicitly instead of losing data.
use gb_io::{
    reader::SeqReader,
    seq::{Location, Seq, Topology},
};
#[cfg(test)]
use popgene_core::Segment;
use popgene_core::{Document, Feature, Molecule, Source};
use rusqlite::{Connection, OpenFlags};
use sha2::{Digest, Sha256};
use std::{fs, io::Write, path::Path, sync::Arc};
use thiserror::Error;
mod abif;
mod canonical_json;
mod duplex;
mod embl;
mod fastq;
mod genbank_display;
mod lineage;
mod location;
mod records;
mod references;
mod save;
mod scf;
mod source_metadata;
mod ztr;
pub use records::{import_records, open_records};
pub use save::{SaveFormat, encode, save, save_warnings};
pub use source_metadata::source_metadata;
pub mod project;
mod snapgene;
pub use snapgene::{NativeAttachment, native_attachment_bytes, native_attachments};
pub use snapgene::{NativeEmbeddedTrace, native_embedded_trace, native_embedded_traces};
pub use snapgene::{NativeEnzymePanel, NativeEnzymePanels, native_enzyme_panels};
pub use snapgene::{NativeHistoryEntry, native_history, native_history_document};
mod snapgene_preserve;

const MAX_FILE: u64 = 128 * 1024 * 1024;
const APP_ID: i32 = 0x504f5047;

#[derive(Debug, Error)]
pub enum Error {
    #[error("{0}")]
    Message(String),
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error(transparent)]
    Core(#[from] popgene_core::Error),
    #[error("Invalid PopGene container: {0}")]
    Sql(#[from] rusqlite::Error),
    #[error("Invalid PopGene data: {0}")]
    Json(#[from] serde_json::Error),
}
type Result<T> = std::result::Result<T, Error>;

pub fn digest(path: &Path) -> Result<[u8; 32]> {
    Ok(Sha256::digest(read_bounded(path)?).into())
}
fn read_bounded(path: &Path) -> Result<Vec<u8>> {
    use std::io::Read;
    let file = fs::File::open(path)?;
    if file.metadata()?.len() > MAX_FILE {
        return Err(Error::Message(
            "File exceeds this prototype's 128 MiB limit.".into(),
        ));
    }
    let mut bytes = vec![];
    file.take(MAX_FILE + 1).read_to_end(&mut bytes)?;
    if bytes.len() as u64 > MAX_FILE {
        return Err(Error::Message("File is too large.".into()));
    }
    Ok(bytes)
}

pub fn open(path: &Path) -> Result<Document> {
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    if ext == "popgene" {
        return load_native(path);
    }
    let bytes = read_bounded(path)?;
    let name = path.file_name().unwrap_or_default().to_string_lossy();
    import(&name, &bytes)
}

pub fn import(name: &str, bytes: &[u8]) -> Result<Document> {
    let mut records = import_records(name, bytes)?;
    if records.len() != 1 {
        return Err(Error::Message(format!(
            "Multiple records found ({}). Use the multi-record import interface to open every sequence.",
            records.len()
        )));
    }
    Ok(records.remove(0))
}

fn import_one(name: &str, bytes: &[u8]) -> Result<Document> {
    if bytes.len() as u64 > MAX_FILE {
        return Err(Error::Message("File too large.".into()));
    }
    if bytes.first() == Some(&9) {
        return snapgene::read(name, bytes);
    }
    if bytes.starts_with(b"\xaeZTR\r\n\x1a\n") {
        return ztr::read(name, bytes);
    }
    if bytes.starts_with(b".scf") {
        return scf::read(name, bytes);
    }
    if bytes.starts_with(b"ABIF") {
        return abif::read(name, bytes);
    }
    let text=std::str::from_utf8(bytes).map_err(|_|Error::Message("Unsupported binary file. PopGene opens SnapGene DNA/RNA/protein, Sanger ABIF, GenBank/GenPept, EMBL, FASTA and legacy .popgene documents.".into()))?;
    let text = text.trim_start_matches('\u{feff}').trim_start();
    let (mut doc, format) = if text.starts_with('>') {
        (fasta(name, text)?, "FASTA")
    } else if text.starts_with("LOCUS") {
        (genbank(text.as_bytes())?, "GenBank")
    } else if text.starts_with("ID   ") {
        (embl::read(text)?, "EMBL")
    } else {
        return Err(Error::Message(
            "Unrecognized format. Choose a SnapGene DNA/RNA/protein, GenBank/GenPept, EMBL, FASTA or ABIF file. Legacy .popgene files can be opened for conversion.".into(),
        ));
    };
    doc.source = Some(Source {
        file_name: name.into(),
        format: format.into(),
        bytes: Arc::from(bytes),
        retained_blocks: vec![],
    });
    if let Some(target) = match format {
        "GenBank" => Some(SaveFormat::GenBank),
        "EMBL" => Some(SaveFormat::Embl),
        _ => None,
    } {
        duplex::read(&mut doc, bytes, target);
        lineage::read(&mut doc, bytes, target);
    }
    doc.validate()?;
    Ok(doc)
}

fn fasta(name: &str, text: &str) -> Result<Document> {
    let mut lines = text.lines();
    let header = lines
        .next()
        .unwrap_or_default()
        .trim_start_matches('>')
        .trim();
    if header.is_empty() {
        return Err(Error::Message("FASTA needs a nonempty record name.".into()));
    }
    let mut sequence = String::new();
    for line in lines {
        if line.starts_with('>') {
            return Err(Error::Message("Multiple records found. This prototype opens one sequence per file; no records have been silently dropped.".into()));
        }
        sequence.push_str(line);
        if sequence.len() > popgene_core::MAX_BASES * 2 {
            return Err(Error::Message("Sequence exceeds prototype limit.".into()));
        }
    }
    let ext = Path::new(name)
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    let (header, molecule) = if let Some(h) = header.strip_suffix(" [molecule=DNA]") {
        (h, Molecule::Dna)
    } else if let Some(h) = header.strip_suffix(" [molecule=RNA]") {
        (h, Molecule::Rna)
    } else if let Some(h) = header.strip_suffix(" [molecule=Protein]") {
        (h, Molecule::Protein)
    } else {
        (
            header,
            match ext.as_str() {
                "faa" | "pep" => Molecule::Protein,
                "frna" => Molecule::Rna,
                _ => Molecule::Dna,
            },
        )
    };
    // Generic FASTA has no molecule tag. Use a protein alphabet only when the
    // symbols cannot represent nucleotides; ambiguous A/C/G-only files remain DNA.
    let generic = !matches!(ext.as_str(), "fna" | "ffn" | "frna" | "faa" | "pep")
        && !text
            .lines()
            .next()
            .unwrap_or_default()
            .contains(" [molecule=");
    let molecule = if molecule == Molecule::Dna
        && generic
        && sequence.chars().any(|c| {
            !c.is_ascii_whitespace() && !"ACGTURYSWKMBDHVN".contains(c.to_ascii_uppercase())
        })
        && Molecule::Protein.normalize(&sequence).is_ok()
    {
        Molecule::Protein
    } else {
        molecule
    };
    let molecule = if molecule == Molecule::Dna
        && generic
        && sequence.bytes().any(|b| matches!(b, b'U' | b'u'))
        && !sequence.bytes().any(|b| matches!(b, b'T' | b't'))
    {
        Molecule::Rna
    } else {
        molecule
    };
    let mut doc = Document::new_typed(header, &sequence, molecule)?;
    if molecule == Molecule::Dna {
        doc.warnings.push("Untyped FASTA defaults to DNA. For RNA or protein, use New Sequence → Import text, or a typed .frna/.faa file.".into());
    }
    doc.warnings.push("FASTA contains no circularity or annotations. Imported as linear; set circularity explicitly if appropriate.".into());
    Ok(doc)
}

fn one_genbank(bytes: &[u8]) -> Result<Seq> {
    // gb-io preallocates from LOCUS length; reject excessive declared sizes before parsing.
    let text = std::str::from_utf8(bytes)
        .map_err(|_| Error::Message("GenBank must be UTF-8 text.".into()))?;
    for line in text.lines().filter(|l| l.starts_with("LOCUS")) {
        let parts: Vec<_> = line.split_whitespace().collect();
        for pair in parts.windows(2) {
            if matches!(pair[1], "bp" | "aa")
                && pair[0].parse::<usize>().unwrap_or(usize::MAX) > popgene_core::MAX_BASES
            {
                return Err(Error::Message(
                    "Declared GenBank length exceeds 10 million bases.".into(),
                ));
            }
        }
    }
    // Avoid unbounded location recursion in upstream parsing.
    let mut depth = 0i32;
    for b in bytes {
        if *b == b'(' {
            depth += 1;
            if depth > 64 {
                return Err(Error::Message("GenBank nesting is too deep.".into()));
            }
        } else if *b == b')' {
            depth = (depth - 1).max(0);
        }
    }
    // gb-io handles nucleotide LOCUS lines only. Normalize a GenPept header for
    // its shared INSDC parser, then restore amino-acid stop symbols afterward.
    let first = text.lines().next().unwrap_or_default();
    let fields: Vec<_> = first.split_whitespace().collect();
    let protein = fields
        .windows(2)
        .any(|p| p[1] == "aa" && p[0].parse::<usize>().is_ok());
    let mut adapted = String::new();
    let mut protein_sequence = None;
    if protein {
        let count_at = fields.iter().position(|f| *f == "aa").unwrap() - 1;
        if count_at < 2 {
            return Err(Error::Message("Invalid GenPept LOCUS line.".into()));
        }
        let topology = if fields.contains(&"circular") {
            "circular"
        } else {
            "linear"
        };
        adapted = format!(
            "LOCUS       {} {} bp PROTEIN {topology} UNK 01-JAN-1970\n",
            fields[1], fields[count_at]
        );
        let mut sequence = String::new();
        let mut origin = false;
        for line in text.lines().skip(1) {
            if line.starts_with("ORIGIN") {
                origin = true;
            }
            if line.starts_with("//") {
                origin = false;
            }
            if origin && !line.starts_with("ORIGIN") {
                for b in line.bytes() {
                    if b.is_ascii_alphabetic() || b == b'*' {
                        sequence.push(b as char);
                    } else if !b.is_ascii_whitespace() && !b.is_ascii_digit() {
                        return Err(Error::Message("Invalid GenPept sequence symbol.".into()));
                    }
                }
                adapted.push_str(&line.replace('*', "X"));
            } else {
                adapted.push_str(line);
            }
            adapted.push('\n');
        }
        protein_sequence = Some(sequence.into_bytes());
    }
    let parse_bytes = if protein { adapted.as_bytes() } else { bytes };
    let mut records = SeqReader::new(parse_bytes);
    let mut record = records
        .next()
        .ok_or_else(|| Error::Message("No GenBank record found.".into()))?
        .map_err(|e| Error::Message(format!("GenBank parse error: {e}")))?;
    if records.next().is_some() {
        return Err(Error::Message(
            "Multiple GenBank records are not yet supported. Split them into separate files."
                .into(),
        ));
    }
    if let Some(sequence) = protein_sequence {
        record.seq = sequence;
    }
    Ok(record)
}

fn genbank(bytes: &[u8]) -> Result<Document> {
    let record = one_genbank(bytes)?;
    let molecule_type = record
        .molecule_type
        .as_deref()
        .unwrap_or("DNA")
        .to_ascii_uppercase();
    let molecule = if molecule_type.contains("RNA") {
        Molecule::Rna
    } else if molecule_type.contains("DNA") {
        Molecule::Dna
    } else if molecule_type == "PROTEIN" {
        Molecule::Protein
    } else {
        return Err(Error::Message(format!(
            "Unsupported GenBank molecule type {molecule_type}."
        )));
    };
    if record.features.len() > 20_000 {
        return Err(Error::Message(
            "More than 20,000 annotations are not supported yet.".into(),
        ));
    }
    let raw = std::str::from_utf8(&record.seq)
        .map_err(|_| Error::Message("Invalid sequence encoding.".into()))?;
    let sequence = if molecule == Molecule::Rna {
        raw.replace(['T', 't'], "U")
    } else {
        raw.into()
    };
    let label = if genbank_display::advertised(&record) {
        record
            .keywords
            .as_ref()
            .filter(|value| !value.trim().is_empty() && value.trim() != ".")
    } else {
        None
    };
    let mut doc = Document::new_typed(
        label
            .map(|label| label.split_whitespace().collect::<Vec<_>>().join(" "))
            .or_else(|| record.name.clone())
            .unwrap_or_else(|| "Untitled".into()),
        &sequence,
        molecule,
    )?;
    doc.single_stranded = molecule_type.starts_with("SS-")
        || (molecule == Molecule::Rna && !molecule_type.starts_with("DS-"));
    if molecule == Molecule::Rna && raw.bytes().any(|b| b == b'T' || b == b't') {
        doc.warnings.push(
            "GenBank declares RNA: T was converted to U; original bytes are retained.".into(),
        );
    }
    doc.circular = record.is_circular();
    doc.description = record.definition.unwrap_or_default();
    doc.accession = record.accession;
    doc.comments = record.comments;
    doc.references = record
        .references
        .iter()
        .map(references::from_genbank)
        .collect();
    for (index, f) in record.features.iter().enumerate() {
        let explicit = location::from_genbank(&f.location)?;
        explicit
            .validate(doc.sequence.len(), doc.circular)
            .map_err(|e| Error::Message(e.to_string()))?;
        let (parts, reverse) = explicit
            .uniform_projection()
            .unwrap_or_else(|| (explicit.local_parts(), false));
        let mut canonical = if parts.len() == 1 {
            Location::simple_range(parts[0].start as i64, parts[0].end as i64)
        } else {
            Location::Join(
                parts
                    .iter()
                    .map(|p| Location::simple_range(p.start as i64, p.end as i64))
                    .collect(),
            )
        };
        if reverse {
            canonical = Location::Complement(Box::new(canonical));
        }
        let feature_location = (canonical != f.location).then_some(explicit);
        let label = f
            .qualifier_values("label")
            .next()
            .or_else(|| f.qualifier_values("gene").next())
            .or_else(|| f.qualifier_values("product").next())
            .unwrap_or(&f.kind)
            .to_owned();
        doc.features.push(Feature {
            location: feature_location,
            id: index as u64 + 1,
            label,
            kind: f.kind.to_string(),
            parts,
            reverse,
            qualifiers: f
                .qualifiers
                .iter()
                .map(|(k, v)| {
                    (
                        k.to_string(),
                        v.as_ref()
                            .map(|s| s.split_whitespace().collect::<Vec<_>>().join(" ")),
                    )
                })
                .collect(),
            edited: false,
        });
    }
    doc.warnings.push("GenBank metadata and original bytes are retained; qualifier whitespace is normalized. Primer oligos are read only from explicit sequence metadata, including any tails and original letter case.".into());
    Ok(doc)
}

pub fn genbank_bytes(doc: &Document) -> Result<Vec<u8>> {
    lineage::write(doc, SaveFormat::GenBank, genbank_bytes_physical)
}

fn genbank_bytes_physical(doc: &Document) -> Result<Vec<u8>> {
    duplex::write(doc, SaveFormat::GenBank, genbank_bytes_plain)
}

fn genbank_bytes_plain(doc: &Document) -> Result<Vec<u8>> {
    doc.validate()?;
    let valid_key =
        |key: &str| !key.is_empty() && key.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_');
    if doc
        .features
        .iter()
        .any(|f| !valid_key(&f.kind) || f.qualifiers.iter().any(|(key, _)| !valid_key(key)))
    {
        return Err(Error::Message(
            "Feature and qualifier keys must be simple INSDC identifiers.".into(),
        ));
    }
    if doc.placeholder {
        return Err(Error::Message(
            "Add sequence before saving this empty placeholder.".into(),
        ));
    }
    let mut record = if let Some(source) = &doc.source {
        if source.format == "GenBank" {
            one_genbank(&source.bytes)?
        } else {
            Seq::empty()
        }
    } else {
        Seq::empty()
    };
    let source_advertised = genbank_display::advertised(&record);
    let extended = genbank_display::uses_extension(doc);
    record.name = Some(save::safe_name(&doc.name));
    record.seq = doc.sequence.as_bytes().to_vec();
    record.len = Some(doc.sequence.len());
    record.topology = if doc.circular {
        Topology::Circular
    } else {
        Topology::Linear
    };
    record.molecule_type = Some(
        match (doc.molecule, doc.single_stranded) {
            (Molecule::Dna, true) => "ss-DNA",
            (Molecule::Dna, false) => "DNA",
            (Molecule::Rna, false) => "ds-RNA",
            (Molecule::Protein, _) => "PROTEIN",
            _ => "RNA",
        }
        .into(),
    );
    record.definition = Some(doc.description.clone());
    record.contig = None;
    record.accession = doc.accession.clone();
    record.comments = doc.comments.clone();
    record.references = doc
        .references
        .iter()
        .enumerate()
        .map(|(i, r)| references::to_genbank(r, i))
        .collect();
    if extended {
        genbank_display::apply_header(doc, &mut record, source_advertised);
    }
    if doc
        .features
        .iter()
        .any(|feature| feature.parts.is_empty() && feature.location.is_none())
    {
        return Err(Error::Message("GenBank has no unbound-primer location. Save native SnapGene to preserve independent oligos, or remove unbound oligos explicitly before conversion.".into()));
    }
    record.features = doc
        .features
        .iter()
        .enumerate()
        .map(|(feature_index, f)| {
            let mut location = if let Some(explicit) = &f.location {
                location::to_genbank(explicit)
            } else if f.parts.len() == 1 {
                Location::simple_range(f.parts[0].start as i64, f.parts[0].end as i64)
            } else {
                Location::Join(
                    f.parts
                        .iter()
                        .map(|p| Location::simple_range(p.start as i64, p.end as i64))
                        .collect(),
                )
            };
            if f.reverse && f.location.is_none() {
                location = Location::Complement(Box::new(location));
            }
            let mut qualifiers: Vec<_> = f
                .qualifiers
                .iter()
                .map(|(k, v)| (k.clone().into(), v.clone()))
                .collect();
            let inferred = f
                .qualifiers
                .iter()
                .find_map(|(key, value)| (key == "label").then_some(value.as_deref()).flatten())
                .or_else(|| {
                    f.qualifiers.iter().find_map(|(key, value)| {
                        (key == "gene").then_some(value.as_deref()).flatten()
                    })
                })
                .or_else(|| {
                    f.qualifiers.iter().find_map(|(key, value)| {
                        (key == "product").then_some(value.as_deref()).flatten()
                    })
                })
                .unwrap_or(&f.kind);
            if inferred != f.label {
                qualifiers.retain(|(key, _)| key != "label");
                qualifiers.insert(0, ("label".into(), Some(f.label.clone())));
            }
            if f.edited {
                let index = f
                    .terminal_display_note()
                    .and_then(|note| {
                        qualifiers.iter().rposition(|(key, value)| {
                            key == "note" && value.as_deref() == Some(note)
                        })
                    })
                    .unwrap_or(qualifiers.len());
                qualifiers.insert(
                    index,
                    (
                        "note".into(),
                        Some("Sequence changed in PopGene; annotation needs review.".into()),
                    ),
                );
            }
            if let Some(note) = f.interchange_display_note_at(feature_index)
                && qualifiers
                    .iter()
                    .rev()
                    .find(|(key, _)| key == "note")
                    .and_then(|(_, value)| value.as_deref())
                    != Some(note.as_str())
            {
                qualifiers.push(("note".into(), Some(note)));
            }
            gb_io::seq::Feature {
                kind: f.kind.clone().into(),
                location,
                qualifiers,
            }
        })
        .collect();
    let mut output = vec![];
    record.write(&mut output)?;
    // gb-io writes a minimally spaced topology and only bp units. Rebuild its
    // LOCUS line with INSDC column widths, including proper GenPept aa units.
    let body = output.iter().position(|b| *b == b'\n').unwrap() + 1;
    let name = record.name.as_deref().unwrap_or("Untitled");
    let length = doc.sequence.len().to_string();
    let name_length = if name.len() + 1 + length.len() <= 28 {
        format!(
            "{name}{}{length}",
            " ".repeat(28 - name.len() - length.len())
        )
    } else {
        format!("{name} {length}")
    };
    let protein = doc.molecule == Molecule::Protein;
    let units = if protein { "aa" } else { "bp" };
    let molecule = if protein {
        ""
    } else {
        record.molecule_type.as_deref().unwrap_or("DNA")
    };
    let topology = if protein {
        ""
    } else if doc.circular {
        "circular"
    } else {
        "linear"
    };
    let date = record
        .date
        .map(|d| d.to_string())
        .unwrap_or_else(|| "01-JAN-1970".into());
    let mut formatted = format!(
        "LOCUS       {name_length} {units}    {molecule:<7} {topology:<8} {} {date}\n",
        record.division
    )
    .into_bytes();
    formatted.extend_from_slice(&output[body..]);
    output = formatted;
    Ok(output)
}

pub fn fasta_bytes(doc: &Document) -> Vec<u8> {
    let mut output = format!(">{}\n", doc.name.replace(['\r', '\n'], " ")).into_bytes();
    for chunk in doc.sequence.as_bytes().chunks(80) {
        output.extend_from_slice(chunk);
        output.push(b'\n');
    }
    output
}

pub fn atomic_export(path: &Path, bytes: &[u8]) -> Result<()> {
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let mut temp = tempfile::NamedTempFile::new_in(parent)?;
    temp.write_all(bytes)?;
    temp.as_file().sync_all()?;
    temp.persist(path).map_err(|e| Error::Io(e.error))?;
    #[cfg(unix)]
    fs::File::open(parent)?.sync_all()?;
    Ok(())
}

fn load_native(path: &Path) -> Result<Document> {
    if fs::metadata(path)?.len() > MAX_FILE {
        return Err(Error::Message("Document exceeds 128 MiB.".into()));
    }
    let db = Connection::open_with_flags(
        path,
        OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX,
    )?;
    let app: i32 = db.pragma_query_value(None, "application_id", |r| r.get(0))?;
    let version: i32 = db.pragma_query_value(None, "user_version", |r| r.get(0))?;
    if app != APP_ID || !(1..=3).contains(&version) {
        return Err(Error::Message(
            "Unrecognized PopGene container or newer schema version.".into(),
        ));
    }
    let (json, original): (String, Option<Vec<u8>>) = db.query_row(
        "SELECT payload,original FROM document WHERE id=1",
        [],
        |row| Ok((row.get(0)?, row.get(1)?)),
    )?;
    let mut doc: Document = serde_json::from_str(&json)?;
    if let Some(source) = &mut doc.source {
        source.bytes = original
            .ok_or_else(|| Error::Message("Missing original source bytes.".into()))?
            .into();
    }
    doc.validate()?;
    Ok(doc)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rusqlite::params;
    fn legacy_fixture(path: &Path, doc: &Document) {
        let db = Connection::open(path).unwrap();
        db.execute_batch(&format!("PRAGMA application_id={APP_ID}; PRAGMA user_version=3; CREATE TABLE document (id INTEGER PRIMARY KEY, payload TEXT NOT NULL, original BLOB);")).unwrap();
        db.execute(
            "INSERT INTO document VALUES (1, ?1, ?2)",
            params![
                serde_json::to_string(doc).unwrap(),
                doc.source.as_ref().map(|s| s.bytes.as_ref())
            ],
        )
        .unwrap();
    }
    #[test]
    fn legacy_molecules_and_placeholders_remain_readable() {
        use popgene_core::creation::{Mode, Request};
        let tmp = tempfile::tempdir().unwrap();
        for molecule in [Molecule::Dna, Molecule::Rna, Molecule::Protein] {
            for mode in [Mode::Sequence, Mode::Unknown, Mode::Placeholder] {
                let mut r = Request::new(molecule);
                r.mode = mode;
                r.text = match molecule {
                    Molecule::Dna => "ACGTNR",
                    Molecule::Rna => "ACGUNR",
                    Molecule::Protein => "MKWVTF*XZ",
                }
                .into();
                let d = r.create().unwrap();
                let path = tmp
                    .path()
                    .join(format!("typed-{molecule:?}-{mode:?}.popgene"));
                legacy_fixture(&path, &d);
                assert_eq!(d, open(&path).unwrap());
                if mode == Mode::Placeholder {
                    assert!(genbank_bytes(&d).is_err());
                }
            }
        }
    }
    #[test]
    fn typed_fasta_and_rna_genbank_keep_molecule_semantics() {
        for molecule in [Molecule::Rna, Molecule::Protein] {
            // ACG alone is ambiguous between all molecule alphabets: type must survive export.
            let d = Document::new_typed("typed", "ACG", molecule).unwrap();
            let filename = if molecule == Molecule::Rna {
                "typed.frna"
            } else {
                "typed.faa"
            };
            let reread = import(filename, &fasta_bytes(&d)).unwrap();
            assert_eq!(reread.molecule, molecule);
            assert_eq!(reread.sequence, d.sequence);
            assert_eq!(reread.name, d.name);
        }
        let d = Document::new_typed("RNA", "ACGUUNRY", Molecule::Rna).unwrap();
        let reread = import("RNA.gb", &genbank_bytes(&d).unwrap()).unwrap();
        assert_eq!(reread.molecule, Molecule::Rna);
        assert_eq!(reread.sequence, d.sequence);
        assert!(reread.single_stranded);
        assert_eq!(
            import("p.faa", b">protein\nMKWVTF*X").unwrap().molecule,
            Molecule::Protein
        );
        assert_eq!(
            import("r.frna", b">RNA\nAUGC").unwrap().molecule,
            Molecule::Rna
        );
        let protein = Document::new_typed("protein", "ACG", Molecule::Protein).unwrap();
        assert_eq!(
            import("protein.gp", &genbank_bytes(&protein).unwrap())
                .unwrap()
                .molecule,
            Molecule::Protein
        );
    }
    const GB: &[u8] = include_bytes!("../../../validation/fixtures/annotated.gb");
    #[test]
    fn import_real_fixture_has_ordered_origin_and_reverse_features() {
        let d = import("annotated.gb", GB).unwrap();
        assert_eq!(d.sequence.len(), 5000);
        assert!(d.circular);
        assert_eq!(d.features.len(), 6);
        assert_eq!(
            d.features[0].parts,
            vec![
                Segment {
                    start: 4900,
                    end: 5000
                },
                Segment { start: 0, end: 120 }
            ]
        );
        assert!(d.features[2].reverse);
        assert_eq!(
            d.features[2].parts,
            vec![Segment {
                start: 220,
                end: 500
            }]
        );
        assert_eq!(
            d.features[0].extract(&d.sequence),
            format!("{}{}", &d.sequence[4900..], &d.sequence[..120])
        );
        assert_eq!(
            format!("{:x}", Sha256::digest(d.sequence.as_bytes())),
            "194483b119b6096a9b16b901a874325f723759c5e05e82c86b4f0113716f184f"
        );
    }
    #[test]
    fn genbank_export_preserves_supported_fields() {
        let d = import("source.gb", GB).unwrap();
        let out = genbank_bytes(&d).unwrap();
        let reread = import("out.gb", &out).unwrap();
        assert_eq!(d.sequence, reread.sequence);
        assert_eq!(d.features, reread.features);
        assert_eq!(d.comments, reread.comments);
        assert_eq!(d.accession, reread.accession);
    }
    #[test]
    fn genbank_and_snapgene_fixture_display_metadata_agrees() {
        let genbank = import("annotated.gb", GB).unwrap();
        let native = import(
            "annotated.dna",
            include_bytes!("../../../validation/fixtures/annotated.dna"),
        )
        .unwrap();
        for label in ["origin_wrap", "audit_forward", "audit_reverse"] {
            let gb = genbank.features.iter().find(|f| f.label == label).unwrap();
            let dna = native.features.iter().find(|f| f.label == label).unwrap();
            assert_eq!(gb.display_color_rgb(), dna.display_color_rgb(), "{label}");
            assert_eq!(gb.primer_sequence(), dna.primer_sequence(), "{label}");
            assert_eq!(gb.parts, dna.parts, "{label}");
            assert_eq!(gb.reverse, dna.reverse, "{label}");
        }
        let forward = native
            .features
            .iter()
            .find(|f| f.label == "overlap_forward")
            .unwrap();
        assert_eq!(
            forward.display_direction(),
            Some(popgene_core::FeatureDirection::Nondirectional)
        );
        let primer = native
            .features
            .iter()
            .find(|f| f.label == "audit_forward")
            .unwrap();
        assert_eq!(primer.primer_melting_temperature(), Some("56"));
        // Accessors do not insert private qualifiers into ordinary GenBank output.
        let exported = genbank_bytes(&genbank).unwrap();
        assert!(!String::from_utf8_lossy(&exported).contains("/primer_sequence="));
    }
    #[test]
    fn legacy_reader_retains_original_bytes() {
        let d = import("source.gb", GB).unwrap();
        let tmp = tempfile::tempdir().unwrap();
        let path = tmp.path().join("test.popgene");
        legacy_fixture(&path, &d);
        assert_eq!(d, open(&path).unwrap());
    }
    #[test]
    fn legacy_trace_preserves_calls_signals_and_original_bytes() {
        let d = import(
            "trace.ab1",
            include_bytes!("../../../validation/fixtures/synthetic-complete.ab1"),
        )
        .unwrap();
        let tmp = tempfile::tempdir().unwrap();
        let path = tmp.path().join("trace.popgene");
        legacy_fixture(&path, &d);
        assert_eq!(d, open(&path).unwrap());
        let db = Connection::open(path).unwrap();
        assert_eq!(
            db.pragma_query_value(None, "user_version", |r| r.get::<_, i32>(0))
                .unwrap(),
            3
        );
    }
    #[test]
    fn version_one_documents_remain_readable() {
        let d = import("source.gb", GB).unwrap();
        let tmp = tempfile::tempdir().unwrap();
        let path = tmp.path().join("old.popgene");
        legacy_fixture(&path, &d);
        let db = Connection::open(&path).unwrap();
        let mut old = serde_json::to_value(&d).unwrap();
        old.as_object_mut().unwrap().remove("trace");
        for field in ["molecule", "single_stranded", "placeholder"] {
            old.as_object_mut().unwrap().remove(field);
        }
        old["source"]
            .as_object_mut()
            .unwrap()
            .remove("retained_blocks");
        db.execute(
            "UPDATE document SET payload=?1",
            params![serde_json::to_string(&old).unwrap()],
        )
        .unwrap();
        db.pragma_update(None, "user_version", 1).unwrap();
        db.close().unwrap();
        assert_eq!(d, open(&path).unwrap());
    }
    #[test]
    fn external_change_blocks_save_without_overwrite() {
        let d = import("source.gb", GB).unwrap();
        let tmp = tempfile::tempdir().unwrap();
        let path = tmp.path().join("test.gb");
        let hash = save(&path, &d, SaveFormat::GenBank, None).unwrap();
        fs::write(&path, b"external change").unwrap();
        assert!(save(&path, &d, SaveFormat::GenBank, Some(hash)).is_err());
        assert_eq!(fs::read(path).unwrap(), b"external change");
    }
    #[test]
    fn fasta_loss_is_explicit_and_bases_survive() {
        let d = import("source.gb", GB).unwrap();
        let reread = import("test.fa", &fasta_bytes(&d)).unwrap();
        assert_eq!(d.sequence, reread.sequence);
        assert!(!reread.circular);
        assert!(reread.features.is_empty());
        assert!(!reread.warnings.is_empty());
    }
    #[test]
    fn multirecord_does_not_silently_drop_records() {
        assert!(
            import("multi.fa", b">a\nACGT\n>b\nAAAA\n")
                .unwrap_err()
                .to_string()
                .contains("Multiple")
        );
    }
    #[test]
    fn invalid_input_and_native_fail() {
        assert!(import("bad", b">a\n$").is_err());
        assert!(import("dna.fna", b">a\nU").is_err());
        assert!(import("dna.fasta", b">a [molecule=DNA]\nZ").is_err());
        assert_eq!(
            import("protein.fasta", b">a\nMKWVTF*XZ").unwrap().molecule,
            Molecule::Protein
        );
        let tmp = tempfile::NamedTempFile::new().unwrap();
        assert!(load_native(tmp.path()).is_err());
    }
    #[test]
    fn oversized_declared_length_rejected_before_parser() {
        assert!(
            import("huge.gb", b"LOCUS       X 999999999999 bp DNA linear\n")
                .unwrap_err()
                .to_string()
                .contains("Declared")
        );
    }
}
