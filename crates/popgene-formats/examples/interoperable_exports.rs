//! Generate deterministic exports for independent readers; no private output format.
use popgene_core::{Document, Feature, Molecule, Segment};
use popgene_formats::{SaveFormat, import, save};
use std::{fs, path::PathBuf};
fn main() {
    let directory = PathBuf::from(std::env::args().nth(1).expect("output directory"));
    fs::create_dir_all(&directory).unwrap();
    let mut doc = import(
        "annotated.gb",
        include_bytes!("../../../validation/fixtures/annotated.gb"),
    )
    .unwrap();
    doc.name = "interoperable_circular".into();
    for format in [
        SaveFormat::GenBank,
        SaveFormat::Embl,
        SaveFormat::Fasta,
        SaveFormat::SnapGene,
    ] {
        save(
            &directory.join(format!("circular.{}", format.extension())),
            &doc,
            format,
            None,
        )
        .unwrap();
    }
    let mut protein = Document::new_typed("protein", "ACGMKWVTF*XZ", Molecule::Protein).unwrap();
    protein.features.push(Feature {
        location: None,
        id: 1,
        label: "region".into(),
        kind: "Region".into(),
        parts: vec![Segment { start: 2, end: 8 }],
        reverse: false,
        qualifiers: vec![("note".into(), Some("Synthetic protein region".into()))],
        edited: false,
    });
    save(
        &directory.join("protein.gp"),
        &protein,
        SaveFormat::GenBank,
        None,
    )
    .unwrap();
    save(
        &directory.join("protein.faa"),
        &protein,
        SaveFormat::Fasta,
        None,
    )
    .unwrap();
    let rna = Document::new_typed("rna", "ACGUNRY", Molecule::Rna).unwrap();
    save(&directory.join("rna.gb"), &rna, SaveFormat::GenBank, None).unwrap();
    save(&directory.join("rna.embl"), &rna, SaveFormat::Embl, None).unwrap();
    save(&directory.join("rna.frna"), &rna, SaveFormat::Fasta, None).unwrap();
    let trace = import(
        "trace.ab1",
        include_bytes!("../../../validation/fixtures/synthetic-complete.ab1"),
    )
    .unwrap();
    save(&directory.join("trace.ab1"), &trace, SaveFormat::Abif, None).unwrap();
    println!(
        "Generated interoperable fixtures in {}",
        directory.display()
    );
}
