//! Original-import metadata for inspection only. Never merged into edited state.
use crate::{Error, MAX_FILE, Result, snapgene};
use popgene_core::Source;

/// Headers/references from the retained original source, not current document
/// metadata. UI callers must label this provenance, especially after edits.
pub fn source_metadata(source: &Source) -> Result<Vec<(String, String)>> {
    if source.bytes.len() as u64 > MAX_FILE {
        return Err(Error::Message(
            "Original source exceeds the metadata inspection limit.".into(),
        ));
    }
    match source.format.as_str() {
        "ABIF" => crate::abif::source_metadata(&source.bytes),
        "SnapGene" => snapgene::source_metadata(&source.bytes),
        "GenBank" => text_headers(&source.bytes, false),
        "EMBL" => text_headers(&source.bytes, true),
        _ => Ok(vec![]),
    }
}

fn text_headers(bytes: &[u8], embl: bool) -> Result<Vec<(String, String)>> {
    let text = std::str::from_utf8(bytes)
        .map_err(|_| Error::Message("Original header is not UTF-8 text.".into()))?;
    let mut fields: Vec<(String, String)> = vec![];
    let mut total = 0;
    for line in text.trim_start_matches('\u{feff}').trim_start().lines() {
        if (!embl
            && (line.starts_with("FEATURES")
                || line.starts_with("ORIGIN")
                || line.starts_with("CONTIG")))
            || (embl && line.starts_with("SQ"))
            || line.trim() == "//"
        {
            break;
        }
        if line.trim().is_empty()
            || (embl && ["XX", "FH", "FT"].iter().any(|key| line.starts_with(key)))
        {
            continue;
        }
        total += line.len();
        if total > 1024 * 1024 || fields.len() >= 10_000 {
            fields.push(("Display limit".into(), "Additional header data remains in the retained original file (metadata display limited to 1 MiB / 10,000 fields).".into()));
            break;
        }
        if embl {
            let Some(key) = line.get(..2) else {
                continue;
            };
            let Some(value) = line.get(5..) else {
                continue;
            };
            let label = match key {
                "ID" => "Record",
                "AC" => "Accession",
                "SV" => "Version",
                "DT" => "Date",
                "DE" => "Definition",
                "KW" => "Keywords",
                "OS" => "Organism",
                "OC" => "Taxonomy",
                "OG" => "Organelle",
                "RN" => "Reference",
                "RP" => "Reference positions",
                "RA" => "Reference authors",
                "RT" => "Reference title",
                "RL" => "Reference journal",
                "RX" => "Reference identifier",
                "RC" => "Reference comment",
                "DR" => "Database cross-reference",
                "CC" => "Comment",
                _ => key,
            };
            if let Some((previous, content)) = fields.last_mut()
                && previous == label
            {
                content.push('\n');
                content.push_str(value.trim_end());
            } else {
                fields.push((label.into(), value.trim_end().into()));
            }
        } else if line.starts_with(' ') {
            if let Some((_, content)) = fields.last_mut() {
                content.push('\n');
                content.push_str(line.trim());
            }
        } else {
            let boundary = line.find(char::is_whitespace).unwrap_or(line.len());
            fields.push((line[..boundary].into(), line[boundary..].trim().into()));
        }
    }
    Ok(fields)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn genbank_references_and_headers_remain_original_after_edits() {
        let mut doc = crate::import(
            "annotated.gb",
            include_bytes!("../../../validation/fixtures/annotated.gb"),
        )
        .unwrap();
        let fields = source_metadata(doc.source.as_ref().unwrap()).unwrap();
        assert!(fields.iter().any(|(k, v)| k == "REFERENCE"
            && v.contains("PopGene audit.")
            && v.contains("SnapGene GenBank format")));
        assert!(
            fields
                .iter()
                .any(|(k, v)| k == "VERSION" && v == "POPGENE823.1")
        );
        assert!(!fields.iter().any(|(_, v)| v.contains("/label")));
        doc.description = "Edited description".into();
        assert_eq!(
            fields,
            source_metadata(doc.source.as_ref().unwrap()).unwrap()
        );
    }

    #[test]
    fn embl_header_reference_lines_are_not_silently_dropped() {
        let fields = text_headers(b"ID   sample; SV 1; linear; DNA; STD; UNC; 4 BP.\nXX\nRN   [1]\nRA   Original author;\nRT   Original title;\nRL   Example journal.\nXX\nFT   misc_feature    1..4\nSQ   Sequence 4 BP;\n     acgt 4\n//\n", true).unwrap();
        assert!(
            fields
                .iter()
                .any(|(k, v)| k == "Reference authors" && v == "Original author;")
        );
        assert!(
            fields
                .iter()
                .any(|(k, v)| k == "Reference journal" && v == "Example journal.")
        );
        assert!(
            !fields
                .iter()
                .any(|(_, v)| v.contains("acgt") || v.contains("misc_feature"))
        );
    }

    #[test]
    fn snapgene_references_and_sequence_class_are_exposed() {
        let doc = crate::import(
            "fixture.dna",
            include_bytes!("../../../validation/fixtures/annotated.dna"),
        )
        .unwrap();
        let fields = source_metadata(doc.source.as_ref().unwrap()).unwrap();
        assert!(fields.iter().any(|(k, v)| k.contains("Reference")
            && v.contains("PopGene audit.")
            && v.contains("Direct Submission")));
        assert!(
            fields
                .iter()
                .any(|(k, v)| k == "SequenceClass" && v == "UNA")
        );
    }
}
