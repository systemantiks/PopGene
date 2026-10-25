//! Single-record EMBL/ENA nucleotide flat files using the shared INSDC feature table.
//! See https://ena-docs.readthedocs.io/en/latest/submit/fileprep/flat-file-example.html
use super::{Error, Result, genbank, genbank_bytes, save::safe_name};
use popgene_core::{Document, Molecule};

fn bad(message: &str) -> Error {
    Error::Message(format!("EMBL: {message}"))
}

pub(super) fn read(text: &str) -> Result<Document> {
    let mut lines = text.lines();
    let first = lines.next().ok_or_else(|| bad("Missing ID line."))?;
    let id = first
        .strip_prefix("ID   ")
        .ok_or_else(|| bad("Missing ID line."))?;
    let fields: Vec<_> = id.split(';').map(str::trim).collect();
    if fields.len() < 4 {
        return Err(bad("Unsupported ID layout."));
    }
    let name = fields[0].split_whitespace().next().unwrap_or("Untitled");
    let length_field = fields.last().unwrap().trim_end_matches('.');
    let mut length_fields = length_field.split_whitespace();
    let length: usize = length_fields
        .next()
        .and_then(|v| v.parse().ok())
        .ok_or_else(|| bad("Invalid declared sequence length."))?;
    if length == 0 || length > popgene_core::MAX_BASES {
        return Err(bad("Declared length must be 1–10,000,000 bases."));
    }
    if length_fields.next().map(str::to_ascii_uppercase).as_deref() != Some("BP") {
        return Err(bad("Only nucleotide EMBL records are supported."));
    }
    let molecule = fields
        .iter()
        .find(|field| {
            let field = field.to_ascii_uppercase();
            field.contains("DNA") || field.contains("RNA")
        })
        .ok_or_else(|| bad("Missing or unsupported molecule type."))?;
    let molecule = molecule
        .split_whitespace()
        .last()
        .unwrap_or(molecule)
        .to_ascii_uppercase();
    if fields.len() >= 7
        && !matches!(
            fields[2].to_ascii_lowercase().as_str(),
            "linear" | "circular"
        )
    {
        return Err(bad("Unsupported or missing topology in the ID line."));
    }
    let circular = fields
        .iter()
        .flat_map(|field| field.split_whitespace())
        .any(|word| word.eq_ignore_ascii_case("circular"));
    let topology = if circular { "circular" } else { "linear" };
    let mut header = format!(
        "LOCUS       {} {length} bp {molecule} {topology} UNK 01-JAN-1970\n",
        safe_name(name)
    );
    let mut features = String::new();
    let mut sequence = String::new();
    let mut description = vec![];
    let mut accession = None;
    let mut comments: Vec<String> = vec![];
    let mut in_comment = false;
    let mut in_sequence = false;
    let mut terminated = false;
    let mut omitted = false;
    for line in lines {
        if !line.starts_with("CC   ") {
            in_comment = false;
        }
        if terminated {
            if !line.trim().is_empty() {
                return Err(bad(
                    "Multiple records or trailing content are not supported; split records into separate files.",
                ));
            }
            continue;
        }
        if line.trim() == "//" {
            terminated = true;
            continue;
        }
        if line.starts_with("ID   ") {
            return Err(bad(
                "Multiple records are not supported; split records into separate files.",
            ));
        }
        if in_sequence {
            for b in line.bytes() {
                if b.is_ascii_alphabetic() {
                    sequence.push(b as char);
                } else if !b.is_ascii_whitespace() && !b.is_ascii_digit() {
                    return Err(bad("Invalid sequence symbol."));
                }
            }
            if sequence.len() > popgene_core::MAX_BASES {
                return Err(bad("Sequence exceeds 10 million bases."));
            }
        } else if line.starts_with("SQ   ") {
            in_sequence = true;
        } else if let Some(value) = line.strip_prefix("DE   ") {
            description.push(value);
        } else if let Some(value) = line.strip_prefix("AC   ") {
            if accession.is_none() {
                accession = value
                    .split(';')
                    .next()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .map(str::to_owned);
            }
        } else if let Some(value) = line.strip_prefix("CC   ") {
            if in_comment {
                let comment = comments.last_mut().unwrap();
                comment.push(' ');
                comment.push_str(value);
            } else {
                comments.push(value.to_owned());
            }
            in_comment = true;
        } else if let Some(value) = line.strip_prefix("FT") {
            if line.len() < 5 || !line.is_char_boundary(2) {
                return Err(bad("Malformed feature table line."));
            }
            features.push_str("  ");
            features.push_str(value);
            features.push('\n');
        } else if !line.trim().is_empty() && !line.starts_with("XX") && !line.starts_with("FH") {
            omitted = true;
        }
    }
    if !terminated || !in_sequence {
        return Err(bad("Missing SQ sequence section or // terminator."));
    }
    if sequence.len() != length {
        return Err(bad("Declared length does not match the sequence."));
    }
    header.push_str(&format!("DEFINITION  {}\n", description.join(" ")));
    if let Some(value) = accession {
        header.push_str(&format!("ACCESSION   {value}\n"));
    }
    for value in comments {
        header.push_str(&format!("COMMENT     {value}\n"));
    }
    if !features.is_empty() {
        header.push_str("FEATURES             Location/Qualifiers\n");
        header.push_str(&features);
    }
    header.push_str("ORIGIN\n");
    header.push_str(&sequence);
    header.push_str("\n//\n");
    let mut doc = genbank(header.as_bytes())?;
    doc.warnings.clear();
    doc.warnings.push("EMBL: exact feature ranges, strand, topology, qualifiers and common header fields are imported. Fuzzy, remote and between-base locations are rejected. Original bytes are retained.".into());
    if omitted {
        doc.warnings.push("Additional EMBL header fields and references are retained only in the original bytes; they are omitted when an edited record is saved.".into());
    }
    Ok(doc)
}

fn field(output: &mut String, key: &str, value: &str) {
    let value = value.split_whitespace().collect::<Vec<_>>().join(" ");
    let mut width = 5;
    output.push_str(key);
    output.push_str("   ");
    for word in value.split_whitespace() {
        if width > 5 && width + word.len() + 1 > 79 {
            output.push('\n');
            output.push_str(key);
            output.push_str("   ");
            width = 5;
        }
        if width > 5 {
            output.push(' ');
            width += 1;
        }
        output.push_str(word);
        width += word.len();
    }
    output.push('\n');
}

pub(super) fn write(doc: &Document) -> Result<Vec<u8>> {
    if doc.molecule == Molecule::Protein || doc.placeholder {
        return Err(bad("Only nonempty nucleotide sequences can be written."));
    }
    let gb =
        String::from_utf8(genbank_bytes(doc)?).map_err(|_| bad("Invalid sequence encoding."))?;
    let molecule = match (doc.molecule, doc.single_stranded) {
        (Molecule::Dna, true) => "ss-DNA",
        (Molecule::Dna, false) => "DNA",
        (Molecule::Rna, false) => "ds-RNA",
        _ => "RNA",
    };
    let topology = if doc.circular { "circular" } else { "linear" };
    let mut output = format!(
        "ID   {}; SV 1; {topology}; {molecule}; STD; UNC; {} BP.\nXX\n",
        safe_name(&doc.name),
        doc.sequence.len()
    );
    if let Some(accession) = &doc.accession {
        field(
            &mut output,
            "AC",
            &format!("{};", accession.replace(['\r', '\n'], " ")),
        );
        output.push_str("XX\n");
    }
    field(&mut output, "DE", &doc.description);
    output.push_str("XX\n");
    for comment in &doc.comments {
        field(&mut output, "CC", comment);
        output.push_str("XX\n");
    }
    let mut in_features = false;
    for line in gb.lines() {
        if line.starts_with("FEATURES") {
            in_features = true;
            output.push_str("FH   Key             Location/Qualifiers\nFH\n");
            continue;
        }
        if line.starts_with("ORIGIN") {
            break;
        }
        if in_features {
            output.push_str("FT");
            output.push_str(&line[2..]);
            output.push('\n');
        }
    }
    if in_features {
        output.push_str("XX\n");
    }
    let count = |base| doc.sequence.bytes().filter(|b| *b == base).count();
    let a = count(b'A');
    let c = count(b'C');
    let g = count(b'G');
    let t = count(b'T');
    output.push_str(&format!(
        "SQ   Sequence {} BP; {a} A; {c} C; {g} G; {t} T; {} other;\n",
        doc.sequence.len(),
        doc.sequence.len() - a - c - g - t
    ));
    let mut total = 0;
    for line in doc.sequence.as_bytes().chunks(60) {
        output.push_str("    ");
        for chunk in line.chunks(10) {
            output.push(' ');
            output.push_str(&String::from_utf8_lossy(chunk).to_ascii_lowercase());
        }
        total += line.len();
        output.push_str(&format!(" {total:>9}\n"));
    }
    output.push_str("//\n");
    Ok(output.into_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{SaveFormat, encode, import};
    #[test]
    fn annotated_circular_embl_retains_ordered_locations_and_qualifiers() {
        let doc = import(
            "source.gb",
            include_bytes!("../../../validation/fixtures/annotated.gb"),
        )
        .unwrap();
        let bytes = encode(&doc, SaveFormat::Embl).unwrap();
        let read = import("copy.embl", &bytes).unwrap();
        assert_eq!(read.sequence, doc.sequence);
        assert_eq!(read.circular, doc.circular);
        assert_eq!(read.features, doc.features);
        assert_eq!(read.accession, doc.accession);
        assert_eq!(encode(&read, SaveFormat::Embl).unwrap(), bytes);
    }
    #[test]
    fn wrapped_comments_preserve_paragraph_boundaries() {
        let mut doc = Document::new("comments", "ACGT").unwrap();
        doc.comments = vec!["A long comment with many words that wraps across several EMBL lines while retaining its original paragraph boundary and complete meaning.".into(), "A separate second comment.".into()];
        let read = import("comments.embl", &write(&doc).unwrap()).unwrap();
        assert_eq!(read.comments, doc.comments);
    }
    #[test]
    fn rna_and_malformed_embl() {
        let doc = Document::new_typed("rna", "ACGUNRY", Molecule::Rna).unwrap();
        let bytes = write(&doc).unwrap();
        let read = import("r.embl", &bytes).unwrap();
        assert_eq!(read.sequence, doc.sequence);
        assert_eq!(read.molecule, Molecule::Rna);
        let text = String::from_utf8(bytes).unwrap();
        assert!(import("bad.embl", text.replace("7 BP.", "8 BP.").as_bytes()).is_err());
        assert!(import("multi.embl", format!("{text}{text}").as_bytes()).is_err());
        assert!(import("truncated.embl", text.trim_end_matches("//\n").as_bytes()).is_err());
        assert!(
            import(
                "huge.embl",
                b"ID   test; SV 1; linear; DNA; STD; UNC; 10000001 BP.\nSQ\n//\n"
            )
            .is_err()
        );
    }
}
