//! Standard Sanger FASTQ (Phred+33), including wrapped sequence/quality lines.
//! https://doi.org/10.1093/nar/gkp1137 — no automatic legacy-offset guessing.
use super::{Error, Result};
use popgene_core::{Document, Molecule, Source};
fn bad(message: &str) -> Error {
    Error::Message(format!("FASTQ (Phred+33): {message}"))
}
pub fn read_records(name: &str, bytes: &[u8]) -> Result<Vec<Document>> {
    let text = std::str::from_utf8(bytes).map_err(|_| bad("Non-UTF8 input"))?;
    let mut lines = text.split_inclusive('\n').peekable();
    let mut offset = 0;
    let mut documents = vec![];
    while let Some(line) = lines.next() {
        let start = offset;
        offset += line.len();
        let header = line
            .trim_end_matches(['\r', '\n'])
            .strip_prefix('@')
            .filter(|s| !s.is_empty())
            .ok_or_else(|| bad("Expected a nonempty @ read header"))?;
        let mut sequence = String::new();
        let mut separator = false;
        for line in lines.by_ref() {
            offset += line.len();
            let line = line.trim_end_matches(['\r', '\n']);
            if let Some(repeated) = line.strip_prefix('+') {
                if !repeated.is_empty() && repeated != header {
                    return Err(bad("Repeated + header differs from @ header"));
                }
                separator = true;
                break;
            }
            if line.bytes().any(|b| b.is_ascii_whitespace()) {
                return Err(bad("Whitespace inside a sequence line"));
            }
            sequence.push_str(line);
            if sequence.len() > popgene_core::MAX_BASES {
                return Err(bad("Read exceeds ten million bases"));
            }
        }
        if !separator || sequence.is_empty() {
            return Err(bad("Missing sequence or + separator"));
        }
        let mut qualities = Vec::with_capacity(sequence.len());
        while qualities.len() < sequence.len() {
            let line = lines
                .next()
                .ok_or_else(|| bad("Quality string is shorter than the sequence"))?;
            offset += line.len();
            for byte in line.trim_end_matches(['\r', '\n']).bytes() {
                if !(33..=126).contains(&byte) {
                    return Err(bad("Quality symbols must be printable ASCII 33–126"));
                }
                qualities.push(Some(byte - 33));
            }
        }
        if qualities.len() != sequence.len() {
            return Err(bad("Quality string is longer than the sequence"));
        }
        let mut document = Document::new(header, &sequence)?;
        document.sequence_qualities = Some(qualities);
        document.source = Some(Source {
            file_name: name.into(),
            format: "FASTQ".into(),
            bytes: bytes[start..offset].into(),
            retained_blocks: vec![],
        });
        document.warnings.push("FASTQ quality encoding is standard Phred+33. Legacy Solexa/Illumina+64 files require explicit conversion; their offset cannot be inferred reliably.".into());
        document.validate()?;
        documents.push(document);
        if documents.len() > 1024 {
            return Err(bad("More than 1,024 reads; split into smaller input files"));
        }
    }
    if documents.is_empty() {
        return Err(bad("No reads"));
    }
    if documents.len() > 1 {
        let count = documents.len();
        for (index, doc) in documents.iter_mut().enumerate() {
            doc.warnings.push(format!("Read {} of {count}; save each open read to a separate file to retain the other input reads.",index+1));
        }
    }
    Ok(documents)
}
pub fn qualities(document: &Document) -> Option<Vec<u8>> {
    document
        .trace
        .as_ref()
        .filter(|t| t.qualities_known)
        .map(|t| t.qualities.clone())
        .or_else(|| {
            document
                .sequence_qualities
                .as_ref()
                .and_then(|qualities| qualities.iter().copied().collect())
        })
}
pub fn supports(document: &Document) -> bool {
    document.molecule != Molecule::Protein
        && (document
            .trace
            .as_ref()
            .filter(|trace| trace.qualities_known)
            .is_some_and(|trace| {
                trace.qualities.len() == document.sequence.len()
                    && trace.qualities.iter().all(|q| *q <= 93)
            })
            || document
                .sequence_qualities
                .as_ref()
                .is_some_and(|qualities| {
                    qualities.len() == document.sequence.len()
                        && qualities.iter().all(|q| q.is_some_and(|q| q <= 93))
                }))
}
pub fn write(document: &Document) -> Result<Vec<u8>> {
    document.validate()?;
    if !supports(document) {
        return Err(bad(
            "Every nucleotide needs a recorded quality in the 0–93 range; unknown confidence cannot be represented as Phred zero",
        ));
    }
    let header = document.name.replace(['\n', '\r'], " ");
    let quality: String = qualities(document)
        .unwrap()
        .iter()
        .map(|q| (q + 33) as char)
        .collect();
    Ok(format!(
        "@{header}\n{}\n+\n{quality}\n",
        document.sequence.replace('U', "T")
    )
    .into_bytes())
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn wrapped_and_special_quality_symbols_preserve_records() {
        let input = b"@one\nAC\nGT\n+one\n!@\n+~\n@two\nN\n+\nI\n";
        let docs = read_records("reads.fq", input).unwrap();
        assert_eq!(docs.len(), 2);
        assert_eq!(docs[0].sequence.as_ref(), "ACGT");
        assert_eq!(
            docs[0].sequence_qualities.as_deref(),
            Some([Some(0), Some(31), Some(10), Some(93)].as_slice())
        );
        assert_eq!(
            docs[0].source.as_ref().unwrap().bytes.as_ref(),
            &input[..22]
        );
        let reopened = read_records("one.fastq", &write(&docs[0]).unwrap()).unwrap();
        assert_eq!(reopened[0].sequence_qualities, docs[0].sequence_qualities);
        assert_eq!(docs[1].sequence_qualities, Some(vec![Some(40)]));
    }
    #[test]
    fn edited_unknown_confidence_is_not_exported_as_phred_zero() {
        let mut document = read_records("read.fastq", b"@r\nACGT\n+\n!III\n")
            .unwrap()
            .remove(0);
        assert_eq!(document.sequence_qualities.as_ref().unwrap()[0], Some(0));
        document.replace(2, 2, "A").unwrap();
        assert_eq!(document.sequence_qualities.as_ref().unwrap()[2], None);
        assert!(!supports(&document));
        assert!(!crate::SaveFormat::Fastq.supports(&document));
        assert!(
            write(&document)
                .unwrap_err()
                .to_string()
                .contains("unknown confidence")
        );
    }
    #[test]
    fn malformed_input_never_silently_drops_reads() {
        for data in [
            b"@x\nA\n+\n".as_slice(),
            b"@x\nA\n+y\nI\n",
            b"@x\nA\n+\nII\n",
            b"@x\nA\n+\n \n",
            b"@x\nA\n+\nI\n@bad\nT\n+\n",
        ] {
            assert!(read_records("bad.fq", data).is_err());
        }
    }
}
