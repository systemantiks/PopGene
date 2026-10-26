//! All-or-nothing multi-record imports. Each document retains its exact record
//! bytes; callers must save records separately rather than overwrite the input.
use super::{Error, Result, import_one, load_native, read_bounded};
use popgene_core::Document;
use std::{ops::Range, path::Path};

const MAX_RECORDS: usize = 1024;

pub fn open_records(path: &Path) -> Result<Vec<Document>> {
    if path
        .extension()
        .is_some_and(|extension| extension.eq_ignore_ascii_case("popgene"))
    {
        return load_native(path).map(|document| vec![document]);
    }
    let bytes = read_bounded(path)?;
    import_records(
        &path.file_name().unwrap_or_default().to_string_lossy(),
        &bytes,
    )
}

/// Open every record, or return an error without a partial result. Binary native
/// files are single documents. FASTA, GenBank/GenPept and EMBL may contain records
/// with independent molecule types and annotations. The retained source of each
/// returned document is its original record slice, including its own metadata.
pub fn import_records(name: &str, bytes: &[u8]) -> Result<Vec<Document>> {
    if bytes.len() as u64 > super::MAX_FILE {
        return Err(Error::Message(
            "File exceeds the 128 MiB import limit.".into(),
        ));
    }
    if bytes.starts_with(b"@") {
        return super::fastq::read_records(name, bytes);
    }
    let ranges = record_ranges(bytes)?;
    let count = ranges.len();
    let mut documents = Vec::with_capacity(count);
    for (index, range) in ranges.into_iter().enumerate() {
        let mut document = import_one(name, &bytes[range])
            .map_err(|error| Error::Message(format!("Record {} of {count}: {error}", index + 1)))?;
        if count > 1 {
            document.warnings.push(format!(
                "Record {} of {count} in {name}. This document retains its own original record; save it to a separate file to preserve the other records in the input.",
                index + 1
            ));
        }
        documents.push(document);
    }
    Ok(documents)
}

fn limit(count: usize) -> Result<()> {
    if count > MAX_RECORDS {
        return Err(Error::Message(format!(
            "The file contains more than {MAX_RECORDS} records. Split it into smaller files; no records were imported."
        )));
    }
    Ok(())
}

fn record_ranges(bytes: &[u8]) -> Result<Vec<Range<usize>>> {
    if bytes.first() == Some(&9) || bytes.starts_with(b"ABIF") {
        return Ok(std::iter::once(0..bytes.len()).collect());
    }
    let Ok(original) = std::str::from_utf8(bytes) else {
        return Ok(std::iter::once(0..bytes.len()).collect());
    };
    let text = original.trim_start_matches('\u{feff}').trim_start();
    let offset = original.len() - text.len();
    if text.starts_with('>') {
        let mut starts = vec![0];
        let mut cursor = offset;
        for (index, line) in text.split_inclusive('\n').enumerate() {
            if index > 0 && line.starts_with('>') {
                starts.push(cursor);
                limit(starts.len())?;
            }
            cursor += line.len();
        }
        let mut result = Vec::with_capacity(starts.len());
        for (index, start) in starts.iter().enumerate() {
            result.push(*start..starts.get(index + 1).copied().unwrap_or(bytes.len()));
        }
        return Ok(result);
    }
    let prefix = if text.starts_with("LOCUS") {
        "LOCUS"
    } else if text.starts_with("ID   ") {
        "ID   "
    } else {
        return Ok(std::iter::once(0..bytes.len()).collect());
    };
    let mut result = vec![];
    let mut start = 0;
    let mut cursor = offset;
    let mut terminated = false;
    for (index, line) in text.split_inclusive('\n').enumerate() {
        let content = line.trim_end_matches(['\r', '\n']);
        if content.starts_with("LOCUS") {
            let fields: Vec<_> = content.split_whitespace().collect();
            if fields.windows(2).any(|pair| {
                matches!(pair[1], "bp" | "aa")
                    && pair[0].parse::<usize>().unwrap_or(usize::MAX) > popgene_core::MAX_BASES
            }) {
                return Err(Error::Message(
                    "Declared GenBank length exceeds 10 million bases.".into(),
                ));
            }
        }
        if terminated && !content.trim().is_empty() {
            if !content.starts_with(prefix) {
                return Err(Error::Message(
                    "Unexpected content after a record terminator; no records were imported."
                        .into(),
                ));
            }
            result.push(start..cursor);
            limit(result.len() + 1)?;
            start = cursor;
            terminated = false;
        } else if index > 0 && !terminated && content.starts_with(prefix) {
            return Err(Error::Message(
                "A record is missing its // terminator; no records were imported.".into(),
            ));
        }
        if content == "//" {
            terminated = true;
        }
        cursor += line.len();
    }
    if !terminated {
        return Err(Error::Message(
            "A record is missing its // terminator; no records were imported.".into(),
        ));
    }
    result.push(start..bytes.len());
    limit(result.len())?;
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{SaveFormat, encode, import, source_metadata};
    use popgene_core::Molecule;

    #[test]
    fn fasta_records_keep_names_molecule_types_and_exact_source_slices() {
        let first = b"\xef\xbb\xbf>dna [molecule=DNA]\r\nACGT\r\n";
        let second = b">rna [molecule=RNA]\nACGU\n";
        let third = b">protein [molecule=Protein]\nMK*\n";
        let bytes = [first.as_slice(), second, third].concat();
        let docs = import_records("mixed.fasta", &bytes).unwrap();
        assert_eq!(docs.len(), 3);
        for (index, (name, molecule, source)) in [
            ("dna", Molecule::Dna, first.as_slice()),
            ("rna", Molecule::Rna, second.as_slice()),
            ("protein", Molecule::Protein, third.as_slice()),
        ]
        .into_iter()
        .enumerate()
        {
            assert_eq!(docs[index].name, name);
            assert_eq!(docs[index].molecule, molecule);
            assert_eq!(docs[index].source.as_ref().unwrap().bytes.as_ref(), source);
            assert_eq!(encode(&docs[index], SaveFormat::Fasta).unwrap(), source);
            assert!(
                docs[index]
                    .warnings
                    .iter()
                    .any(|warning| warning.contains("separate file"))
            );
        }
        assert!(import("mixed.fasta", &bytes).is_err());
    }

    #[test]
    fn genbank_records_keep_independent_references_and_annotations() {
        let record = |name: &str, author: &str, sequence: &str| {
            format!(
                "LOCUS       {name} 6 bp DNA linear UNK 01-JAN-2000\nDEFINITION  {name}.\nREFERENCE   1 (bases 1 to 6)\n  AUTHORS   {author}\n  TITLE     Direct Submission\nFEATURES             Location/Qualifiers\n     misc_feature    2..4\n                     /label=\"{name} annotation\"\nORIGIN\n        1 {sequence}\n//\n"
            )
        };
        let first = record("first", "Alpha", "acgtac");
        let second = record("second", "Beta", "tttaaa");
        let docs = import_records("multi.gb", format!("{first}{second}").as_bytes()).unwrap();
        assert_eq!(
            docs.iter()
                .map(|doc| doc.sequence.as_ref())
                .collect::<Vec<_>>(),
            ["ACGTAC", "TTTAAA"]
        );
        for (index, author) in ["Alpha", "Beta"].iter().enumerate() {
            assert_eq!(docs[index].features.len(), 1);
            let metadata = source_metadata(docs[index].source.as_ref().unwrap()).unwrap();
            assert!(metadata.iter().any(|(_, text)| text.contains(author)));
            assert!(
                !metadata
                    .iter()
                    .any(|(_, text)| text.contains(if index == 0 { "Beta" } else { "Alpha" }))
            );
        }
        assert_eq!(
            encode(&docs[0], SaveFormat::GenBank).unwrap(),
            first.as_bytes()
        );
        assert_eq!(
            encode(&docs[1], SaveFormat::GenBank).unwrap(),
            second.as_bytes()
        );
    }

    #[test]
    fn embl_records_are_atomic_and_allow_whitespace_between_records() {
        let record = |name: &str, sequence: &str| {
            format!(
                "ID   {name}; SV 1; linear; DNA; STD; UNC; 4 BP.\nDE   {name}.\nSQ   Sequence 4 BP;\n     {sequence} 4\n//\n"
            )
        };
        let first = record("first", "acgt");
        let second = record("second", "tgca");
        let bytes = format!("{first}\n{second}");
        let docs = import_records("multi.embl", bytes.as_bytes()).unwrap();
        assert_eq!(docs.len(), 2);
        assert_eq!(docs[1].sequence.as_ref(), "TGCA");
        assert!(
            import_records(
                "bad.embl",
                format!("{first}{}", second.replace("tgca", "t?ca")).as_bytes()
            )
            .is_err()
        );
        assert!(
            import_records("bad.embl", format!("{first}trailing garbage\n").as_bytes()).is_err()
        );
        assert!(
            import_records(
                "bad.embl",
                format!("{}{}", first.replace("//\n", ""), second).as_bytes()
            )
            .is_err()
        );
        let folder = tempfile::tempdir().unwrap();
        let path = folder.path().join("multi.embl");
        std::fs::write(&path, bytes).unwrap();
        assert_eq!(open_records(&path).unwrap().len(), 2);
    }

    #[test]
    fn record_limit_and_invalid_later_fasta_do_not_return_partial_results() {
        assert!(import_records("bad.fa", b">one\nACGT\n>two\n??\n").is_err());
        let many = ">one\nA\n".repeat(MAX_RECORDS + 1);
        assert!(
            import_records("many.fa", many.as_bytes())
                .unwrap_err()
                .to_string()
                .contains("1024")
        );
    }
}
