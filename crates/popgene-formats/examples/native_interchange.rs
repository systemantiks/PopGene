//! Regenerate synthetic native molecule fixtures; optionally verify SnapGene-authored inputs.
//! Usage: cargo run -p popgene-formats --example native_interchange -- OUTPUT_DIR [INPUT_DIR]
use popgene_core::{Document, Feature, Molecule, Segment};
use popgene_formats::{SaveFormat, encode, open};
use std::{error::Error, path::PathBuf};
fn main() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args_os().skip(1);
    let output = PathBuf::from(args.next().ok_or("Provide output directory")?);
    let inputs = args.next().map(PathBuf::from);
    std::fs::create_dir_all(&output)?;
    for (name, molecule, sequence) in [
        ("synthetic-native-rna", Molecule::Rna, "ACGUAUGCUUAGNNRY"),
        (
            "synthetic-native-protein",
            Molecule::Protein,
            "MKWVTFISLLFLFSSAYSRG*X",
        ),
    ] {
        let mut doc = Document::new_typed(name, sequence, molecule)?;
        doc.circular = molecule == Molecule::Rna;
        doc.single_stranded = false;
        doc.description = "Synthetic native interoperability fixture; no biological claim.".into();
        doc.features.push(Feature {
            location: None,
            id: 1,
            label: "Synthetic domain".into(),
            kind: "misc_feature".into(),
            parts: vec![Segment { start: 1, end: 8 }],
            reverse: false,
            qualifiers: vec![("note".into(), Some("Synthetic domain metadata".into()))],
            edited: false,
        });
        let extension = SaveFormat::SnapGene.extension_for(&doc);
        let path = output.join(format!("{name}.{extension}"));
        let bytes = encode(&doc, SaveFormat::SnapGene)?;
        if path.exists() && std::fs::read(&path)? != bytes {
            return Err(format!("Refusing to overwrite different {}", path.display()).into());
        }
        std::fs::write(&path, bytes)?;
        let restored = open(&path)?;
        assert_eq!(restored.molecule, molecule);
        assert_eq!(restored.sequence, doc.sequence);
        println!(
            "{}: {} {}",
            path.display(),
            restored.sequence.len(),
            molecule.label()
        );
    }
    if let Some(inputs) = inputs {
        for name in ["native-rna.rna", "native-protein.prot"] {
            let input = inputs.join(name);
            let mut doc = open(&input)?;
            assert_eq!(encode(&doc, SaveFormat::SnapGene)?, std::fs::read(&input)?);
            doc.description = "Edited in PopGene; original imported annotations retained.".into();
            let path = output.join(format!("edited-{name}"));
            let bytes = encode(&doc, SaveFormat::SnapGene)?;
            if path.exists() && std::fs::read(&path)? != bytes {
                return Err(format!("Refusing to overwrite different {}", path.display()).into());
            }
            std::fs::write(&path, bytes)?;
            println!(
                "{}: {} {}",
                path.display(),
                doc.sequence.len(),
                doc.molecule.label()
            );
        }
    }
    Ok(())
}
