use super::super::{
    Document, Feature, Segment, native_primer, packet, read, serialize_node, write, write_warnings,
    xml,
};
use super::REFUSAL;

const OBSERVED: &[u8] = include_bytes!("../../../../../validation/fixtures/annotated.dna");
fn blocks(bytes: &[u8]) -> Vec<(u8, Vec<u8>)> {
    crate::snapgene_preserve::packets(bytes)
        .unwrap()
        .into_iter()
        .map(|(kind, body)| (kind, body.to_vec()))
        .collect()
}
fn primer_packets(bytes: &[u8]) -> Vec<Vec<u8>> {
    blocks(bytes)
        .into_iter()
        .filter_map(|(kind, body)| (kind == 5).then_some(body))
        .collect()
}
fn alter_packet(transform: impl Fn(&mut super::Node)) -> Vec<u8> {
    let mut bytes = vec![];
    for (kind, body) in blocks(OBSERVED) {
        if kind == 5 {
            let mut root = xml(&body).unwrap();
            transform(&mut root);
            let mut body = String::new();
            serialize_node(&root, &mut body);
            packet(&mut bytes, kind, body.as_bytes()).unwrap();
        } else {
            packet(&mut bytes, kind, &body).unwrap();
        }
    }
    bytes
}
fn add_feature(doc: &mut Document) {
    doc.annotate(
        "unrelated annotation".into(),
        Segment { start: 12, end: 16 },
        false,
    )
    .unwrap();
}
fn add_primer(doc: &mut Document) {
    let id = doc
        .features
        .iter()
        .chain(&doc.duplex_annotations.features)
        .map(|f| f.id)
        .max()
        .unwrap_or(0)
        + 1;
    doc.features.push(Feature {
        id,
        label: "new oligo".into(),
        kind: "primer_bind".into(),
        parts: vec![Segment { start: 0, end: 12 }],
        reverse: false,
        location: None,
        qualifiers: vec![("primer_sequence".into(), Some(doc.sequence[..12].into()))],
        edited: false,
    });
}
fn settings(bytes: &[u8]) -> super::Node {
    xml(&primer_packets(bytes)[0])
        .unwrap()
        .child("HybridizationParams")
        .unwrap()
        .clone()
}
fn serialized(node: &super::Node) -> String {
    let mut text = String::new();
    serialize_node(node, &mut text);
    text
}
fn assert_refused(doc: &Document) {
    assert!(write_warnings(doc).iter().any(|w| w == REFUSAL));
    assert!(write(doc).unwrap_err().to_string().contains(REFUSAL));
}

#[test]
fn actual_observed_settings_and_complete_packet_survive_unrelated_feature_changes() {
    let source = read("annotated.dna", OBSERVED).unwrap();
    let params = settings(OBSERVED);
    assert_eq!(params.attrs.len(), 5);
    assert_eq!(params.attr("minContinuousMatchLen"), "10");
    assert_eq!(params.attr("allowMismatch"), "1");
    assert_eq!(params.attr("minMeltingTemperature"), "40");
    assert_eq!(params.attr("showAdditionalFivePrimeMatches"), "1");
    assert_eq!(params.attr("minimumFivePrimeAnnealing"), "15");
    for mode in 0..3 {
        let mut doc = source.clone();
        match mode {
            0 => add_feature(&mut doc),
            1 => {
                doc.features
                    .iter_mut()
                    .find(|f| f.kind != "primer_bind")
                    .unwrap()
                    .label = "changed feature name".into()
            }
            _ => {
                let index = doc
                    .features
                    .iter()
                    .position(|f| f.kind != "primer_bind")
                    .unwrap();
                doc.features.remove(index);
            }
        }
        let output = write(&doc).unwrap();
        assert_eq!(primer_packets(&output), primer_packets(OBSERVED));
        let restored = read("updated.dna", &output).unwrap();
        assert_eq!(
            restored
                .features
                .iter()
                .filter(|f| f.kind == "primer_bind")
                .count(),
            source
                .features
                .iter()
                .filter(|f| f.kind == "primer_bind")
                .count()
        );
        assert_eq!(write(&restored).unwrap(), output);
    }
}

#[test]
fn simple_settings_survive_added_changed_deleted_and_last_deleted_primers() {
    let source = read("annotated.dna", OBSERVED).unwrap();
    let expected = serialized(&settings(OBSERVED));
    for mode in 0..4 {
        let mut doc = source.clone();
        match mode {
            0 => add_primer(&mut doc),
            1 => {
                let index = doc
                    .features
                    .iter()
                    .position(|f| native_primer(f, &doc).is_some())
                    .unwrap();
                doc.features[index].label = "changed oligo name".into();
            }
            2 => {
                let index = doc
                    .features
                    .iter()
                    .position(|f| f.kind == "primer_bind")
                    .unwrap();
                doc.features.remove(index);
            }
            _ => doc.features.retain(|f| f.kind != "primer_bind"),
        }
        let bytes = write(&doc).unwrap();
        let packets = primer_packets(&bytes);
        assert_eq!(packets.len(), 1);
        assert_eq!(serialized(&settings(&bytes)), expected);
        let root = xml(&packets[0]).unwrap();
        assert_eq!(root.children[0].name, "HybridizationParams");
        assert!(root.attr("nextValidID").parse::<u64>().unwrap() >= 2);
        if mode == 3 {
            assert!(root.children.iter().all(|node| node.name != "Primer"));
        }
        let restored = read("updated.dna", &bytes).unwrap();
        assert_eq!(
            restored
                .features
                .iter()
                .filter(|f| f.kind == "primer_bind")
                .count(),
            doc.features
                .iter()
                .filter(|f| f.kind == "primer_bind")
                .count()
        );
        assert_eq!(write(&restored).unwrap(), bytes);
    }
}

#[test]
fn original_counter_lower_bound_prevents_reusing_retired_native_primer_ids() {
    let bytes = alter_packet(|root| {
        root.attrs.insert("nextValidID".into(), "400".into());
    });
    let mut doc = read("counter.dna", &bytes).unwrap();
    add_primer(&mut doc);
    let output = write(&doc).unwrap();
    let root = xml(&primer_packets(&output)[0]).unwrap();
    assert_eq!(root.attr("nextValidID"), "401");
    assert_eq!(
        root.children
            .iter()
            .find(|node| node.attr("name") == "new oligo")
            .unwrap()
            .attr("recentID"),
        "400"
    );
    assert_eq!(
        serialized(root.child("HybridizationParams").unwrap()),
        serialized(&settings(&bytes))
    );
}

#[test]
fn unknown_and_interspersed_envelopes_are_kept_exact_until_primer_set_changes() {
    for mode in 0..6 {
        let bytes = alter_packet(|root| match mode {
            0 => {
                root.attrs.insert("recycledIDs".into(), "9,11".into());
            }
            1 => root.children.insert(
                0,
                super::Node {
                    name: "UnknownSettings".into(),
                    attrs: [("primerID".into(), "1".into())].into(),
                    ..Default::default()
                },
            ),
            2 => {
                root.children[0]
                    .attrs
                    .insert("unknownBindingPolicy".into(), "retain".into());
            }
            3 => {
                let params = root.children.remove(0);
                root.children.insert(1, params);
            }
            4 => root.children.insert(0, root.children[0].clone()),
            _ => {
                root.attrs
                    .insert("nextValidID".into(), "invalid-source-counter".into());
            }
        });
        let mut doc = read("opaque.dna", &bytes).unwrap();
        add_feature(&mut doc);
        assert_eq!(
            primer_packets(&write(&doc).unwrap()),
            primer_packets(&bytes)
        );
        add_primer(&mut doc);
        assert_refused(&doc);
        doc.features.retain(|feature| feature.kind != "primer_bind");
        assert_refused(&doc);
    }
}

#[test]
fn existing_repeated_packets_and_duplicate_native_ids_preserve_original_record_identity() {
    let mut bytes = vec![];
    for (kind, body) in blocks(OBSERVED) {
        packet(&mut bytes, kind, &body).unwrap();
        if kind == 5 {
            packet(&mut bytes, kind, &body).unwrap();
        }
    }
    let mut doc = read("repeated.dna", &bytes).unwrap();
    add_feature(&mut doc);
    assert_eq!(
        primer_packets(&write(&doc).unwrap()),
        primer_packets(&bytes)
    );
    let primer_indices = doc
        .features
        .iter()
        .enumerate()
        .filter_map(|(i, f)| (f.kind == "primer_bind").then_some(i))
        .collect::<Vec<_>>();
    assert_eq!(primer_indices.len(), 4);
    // The two packet bodies/IDs are identical, but they are distinct records.
    // Content equality alone must not mistake reversing their occurrence order
    // for preserving the source envelope associations.
    doc.features.swap(primer_indices[0], primer_indices[2]);
    assert_refused(&doc);
    doc.features.swap(primer_indices[0], primer_indices[2]);
    add_primer(&mut doc);
    assert_refused(&doc);
}

#[test]
fn empty_original_settings_packet_remains_after_feature_or_primer_edits() {
    let bytes = alter_packet(|root| root.children.retain(|node| node.name != "Primer"));
    let mut doc = read("empty.dna", &bytes).unwrap();
    add_feature(&mut doc);
    assert_eq!(
        primer_packets(&write(&doc).unwrap()),
        primer_packets(&bytes)
    );
    add_primer(&mut doc);
    let output = write(&doc).unwrap();
    assert_eq!(
        serialized(&settings(&output)),
        serialized(&settings(&bytes))
    );
    assert_eq!(
        xml(&primer_packets(&output)[0])
            .unwrap()
            .children
            .iter()
            .filter(|node| node.name == "Primer")
            .count(),
        1
    );
}

#[test]
fn envelope_comments_remain_raw_and_prevent_unsupported_rebuild() {
    let mut bytes = vec![];
    for (kind, body) in blocks(OBSERVED) {
        if kind == 5 {
            let text = String::from_utf8(body).unwrap().replace(
                "</Primers>",
                "<!-- original packet-scoped metadata --><?policy retain?></Primers>",
            );
            packet(&mut bytes, kind, text.as_bytes()).unwrap();
        } else {
            packet(&mut bytes, kind, &body).unwrap();
        }
    }
    let mut doc = read("comments.dna", &bytes).unwrap();
    add_feature(&mut doc);
    assert_eq!(
        primer_packets(&write(&doc).unwrap()),
        primer_packets(&bytes)
    );
    add_primer(&mut doc);
    assert_refused(&doc);
}
