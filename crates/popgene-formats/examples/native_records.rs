//! Validate locally observed native records and generate edited files for native UI review.
//! Usage: native_records INPUT_DIR OUTPUT_DIR. Existing different output is refused.
use popgene_formats::{
    SaveFormat, encode, native_attachment_bytes, native_attachments, native_history,
    native_history_document, open,
};
use std::{error::Error, path::PathBuf};
fn main() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args_os().skip(1);
    let inputs = PathBuf::from(args.next().ok_or("INPUT_DIR required")?);
    let output = PathBuf::from(args.next().ok_or("OUTPUT_DIR required")?);
    std::fs::create_dir_all(&output)?;
    for name in [
        "reference-details.dna",
        "reference-after-replace.dna",
        "reference-left-threeprime.dna",
        "reference-both-threeprime.dna",
        "reference-left-fiveprime.dna",
    ] {
        let path = inputs.join(name);
        let mut doc = open(&path)?;
        let source = doc.source.as_ref().unwrap();
        assert_eq!(encode(&doc, SaveFormat::SnapGene)?, std::fs::read(&path)?);
        let attachments = native_attachments(source)?;
        assert_eq!(attachments.len(), 1);
        if let Some(a) = attachments.first() {
            assert_eq!(
                native_attachment_bytes(source, a.id)?,
                std::fs::read(inputs.join("attachment-note.txt"))?
            );
        }
        assert_eq!(doc.references[0].title, "Synthetic reference title");
        let history = native_history(source)?;
        println!(
            "{name}: {} actual top residues, {} references, {} attachment bytes, {} history nodes; ends {:?}",
            doc.sequence.len(),
            doc.references.len(),
            attachments.first().map_or(0, |a| a.size),
            history.len(),
            doc.duplex_ends
        );
        if name == "reference-after-replace.dna" {
            assert_eq!(
                native_history_document(source, 0)?.sequence.as_ref(),
                "AAACCCGGGTTT"
            );
        }
        let index = 5;
        doc.replace(index, index + 1, "A")?;
        let bytes = encode(&doc, SaveFormat::SnapGene)?;
        let destination = output.join(format!("popgene-{name}"));
        if destination.exists() && std::fs::read(&destination)? != bytes {
            return Err(format!("Refusing overwrite of {}", destination.display()).into());
        }
        std::fs::write(&destination, &bytes)?;
        let restored = open(&destination)?;
        assert_eq!(restored.sequence, doc.sequence);
        assert_eq!(restored.duplex_ends, doc.duplex_ends);
        let entries = native_history(restored.source.as_ref().unwrap())?;
        let old = entries
            .iter()
            .find(|e| e.parent == Some(entries[0].id))
            .unwrap();
        let recovered = native_history_document(restored.source.as_ref().unwrap(), old.id)?;
        assert_eq!(recovered.sequence, open(&path)?.sequence);
        println!("  validated edited output {}", destination.display());
    }
    Ok(())
}
