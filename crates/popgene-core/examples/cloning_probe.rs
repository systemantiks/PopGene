//! JSON-lines bridge to actual cloning APIs for independent validation.
use popgene_core::{
    Document, Feature, FeatureLocation, Segment,
    cloning::{self, AssemblyInput, LigationOptions, PcrOptions},
    restriction::{self, DigestConditions, ENZYMES},
    translation, workflow,
};
use serde_json::{Value, json};
use std::io::{self, BufRead};
fn number(value: &Value, key: &str) -> usize {
    value[key].as_u64().unwrap() as usize
}
fn document(value: &Value) -> Document {
    let mut d = Document::new(
        value["name"].as_str().unwrap_or("input"),
        value["sequence"].as_str().unwrap(),
    )
    .unwrap();
    d.circular = value["circular"].as_bool().unwrap_or(false);
    if let Some(annotations) = value["annotations"].as_array() {
        for f in annotations {
            d.annotate(
                f["label"].as_str().unwrap().into(),
                Segment {
                    start: number(f, "start"),
                    end: number(f, "end"),
                },
                f["reverse"].as_bool().unwrap_or(false),
            )
            .unwrap();
            let feature = d.features.last_mut().unwrap();
            if let Some(kind) = f["kind"].as_str() {
                feature.kind = kind.into();
            }
            if let Some(qualifiers) = f.get("qualifiers") {
                feature.qualifiers = serde_json::from_value(qualifiers.clone()).unwrap();
            }
            if let Some(location) = f.get("location") {
                let location: FeatureLocation = serde_json::from_value(location.clone()).unwrap();
                (feature.parts, feature.reverse) = location.uniform_projection().unwrap();
                feature.location = Some(location);
            }
        }
    }
    d.validate().unwrap();
    d
}
fn primer(value: &Value, reverse: bool, n: usize) -> Feature {
    let start = number(value, "start");
    let end = start + number(value, "length");
    let parts = if end > n {
        vec![
            Segment { start, end: n },
            Segment {
                start: 0,
                end: end - n,
            },
        ]
    } else {
        vec![Segment { start, end }]
    };
    Feature {
        id: if reverse { 1001 } else { 1000 },
        label: if reverse { "reverse" } else { "forward" }.into(),
        kind: "primer_bind".into(),
        parts,
        location: None,
        reverse,
        qualifiers: vec![(
            "primer_sequence".into(),
            Some(value["sequence"].as_str().unwrap().into()),
        )],
        edited: false,
    }
}
fn run(v: &Value) -> Result<Document, popgene_core::Error> {
    let conditions = DigestConditions::default();
    let enzyme = || {
        ENZYMES
            .iter()
            .position(|e| e.name == v["enzyme"].as_str().unwrap())
            .unwrap()
    };
    let inputs = || {
        v["inputs"]
            .as_array()
            .unwrap()
            .iter()
            .map(|d| AssemblyInput {
                document: document(d),
                reverse: d["reverse"].as_bool().unwrap_or(false),
            })
            .collect::<Vec<_>>()
    };
    match v["method"].as_str().unwrap() {
        "pcr" => {
            let d = document(&v["template"]);
            cloning::pcr_product(
                &d,
                &primer(&v["forward"], false, d.sequence.len()),
                &primer(&v["reverse"], true, d.sequence.len()),
                PcrOptions::default(),
                "PCR",
            )
        }
        "digest_ligation" => {
            let source = document(&v["template"]);
            let enzyme = enzyme();
            let digest = restriction::digest_document(&source, &[enzyme], conditions).unwrap();
            let products = workflow::digest_products(
                &source,
                &digest.fragments,
                &[ENZYMES[enzyme].name.into()],
                conditions,
            )?;
            let indices = v["fragments"].as_array().unwrap();
            let inputs = indices
                .iter()
                .map(|i| AssemblyInput {
                    document: products[i.as_u64().unwrap() as usize].clone(),
                    reverse: v["reverse"].as_bool().unwrap_or(false),
                })
                .collect::<Vec<_>>();
            cloning::ligate(
                &inputs,
                LigationOptions {
                    circular: v["circular"].as_bool().unwrap_or(false),
                    ..Default::default()
                },
                "ligation",
            )
        }
        "overlap" => cloning::assemble_overlaps(
            &inputs(),
            number(v, "minimum"),
            number(v, "maximum"),
            v["circular"].as_bool().unwrap_or(false),
            "overlap",
        ),
        "golden_gate" => cloning::golden_gate(
            &inputs(),
            &v["fragments"]
                .as_array()
                .unwrap()
                .iter()
                .map(|i| i.as_u64().unwrap() as usize)
                .collect::<Vec<_>>(),
            enzyme(),
            conditions,
            "Golden Gate",
        ),
        _ => panic!("unknown oracle method"),
    }
}
fn main() {
    for line in io::stdin().lock().lines() {
        let v: Value = serde_json::from_str(&line.unwrap()).unwrap();
        match run(&v) {
            Ok(d) => println!(
                "{}",
                json!({"sequence":d.sequence,"circular":d.circular,"opposite":d.duplex_ends.as_ref().map(|e|&e.opposite_sequence),"annotations":d.features.iter().map(|f| {
                    let cds = translation::translate_cds(&d, f).ok();
                    json!({"label":f.label,"parts":f.parts,"reverse":f.reverse,"sequence":d.feature_sequence(f).unwrap(),"location":f.location.as_ref().map(FeatureLocation::expression),"qualifiers":f.qualifiers,"phase":cds.as_ref().map(|t|t.phase+1),"translation":cds.map(|t|String::from_utf8(t.amino_acids).unwrap())})
                }).collect::<Vec<_>>() })
            ),
            Err(e) => println!("{}", json!({"error":e.to_string()})),
        }
    }
}
