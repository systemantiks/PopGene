use super::*;

fn primer(id: u64, label: &str, start: usize, reverse: bool) -> Feature {
    let mut feature = Feature {
        id,
        label: label.into(),
        kind: "primer_bind".into(),
        parts: vec![Segment {
            start,
            end: start + 20,
        }],
        location: None,
        reverse,
        qualifiers: vec![],
        edited: false,
    };
    feature.qualifiers.push((
        "primer_sequence".into(),
        Some(feature.extract(&"ACGT".repeat(125))),
    ));
    feature.qualifiers.push((
        "primer_reaction".into(),
        Some(serde_json::to_string(&popgene_core::primer::Reaction::default()).unwrap()),
    ));
    feature
}

#[test]
fn generated_primer_feature_colors_survive_native_save_and_reopen() {
    let mut doc = Document::new("designed primers", &"ACGT".repeat(125)).unwrap();
    doc.features = vec![
        primer(1, "Forward", 139, false),
        primer(2, "Reverse", 394, true),
    ];
    assert!(doc.features.iter().all(|f| f.display_color_rgb().is_none()));
    let colors = [0xb985d0, 0xb985d0];
    let encoded = write(&doc).unwrap();
    let reopened = read("primers.dna", &encoded).unwrap();
    assert!(
        !crate::snapgene_preserve::packets(&encoded)
            .unwrap()
            .iter()
            .any(|(kind, _)| *kind == 5)
    );
    for (index, (original, restored)) in doc.features.iter().zip(&reopened.features).enumerate() {
        assert_eq!(original.effective_color_rgb(index), colors[index]);
        assert_eq!(restored.display_color_rgb(), Some(colors[index]));
        assert_eq!(restored.effective_color_rgb(index), colors[index]);
        assert_eq!(restored.parts, original.parts);
        assert_eq!(restored.reverse, original.reverse);
        assert_eq!(restored.primer_sequence(), original.primer_sequence());
        assert_eq!(
            restored.qualifier("primer_reaction"),
            original.qualifier("primer_reaction")
        );
    }
    assert_eq!(write(&reopened).unwrap(), encoded);
}

#[test]
fn rebuilt_defaults_and_imported_segment_colors_survive_native_packet_grouping() {
    let mut doc = Document::new("mixed colors", &"ACGT".repeat(125)).unwrap();
    let mut native_primer = primer(2, "Native blue primer", 24, false);
    native_primer
        .qualifiers
        .retain(|(key, _)| key != "primer_reaction");
    native_primer
        .qualifiers
        .push(("snapgene_color".into(), Some("blue".into())));
    doc.features = vec![
        primer(1, "Uncolored forward", 139, false),
        native_primer,
        primer(3, "Uncolored reverse", 394, true),
    ];
    let mut imported_bytes = vec![];
    packet(&mut imported_bytes, 9, b"SnapGene\0\x01\0\x0f\0\x14").unwrap();
    packet(&mut imported_bytes, 0, b"\x02ACGTACGTACGTACGTACGT").unwrap();
    packet(&mut imported_bytes, 10, b"<Features><Feature name=\"Imported colored domains\" type=\"misc_feature\"><Segment range=\"1-10\" type=\"standard\" color=\"#aabbcc\" name=\"First\"/><Segment range=\"11-20\" type=\"standard\" color=\"#dd3311\" name=\"Second\"/></Feature></Features>").unwrap();
    let mut imported = read("imported.dna", &imported_bytes)
        .unwrap()
        .features
        .remove(0);
    imported.id = 4;
    doc.features.push(imported);
    let expected = doc
        .features
        .iter()
        .enumerate()
        .map(|(index, f)| {
            (
                f.label.clone(),
                f.effective_color_rgb(index),
                f.display_segments()
                    .iter()
                    .map(|s| s.color_rgb.unwrap_or_else(|| f.effective_color_rgb(index)))
                    .collect::<Vec<_>>(),
            )
        })
        .collect::<Vec<_>>();
    let encoded = write(&doc).unwrap();
    let reopened = read("mixed.dna", &encoded).unwrap();
    assert_eq!(
        reopened.features.last().unwrap().label,
        "Native blue primer"
    );
    let packets = crate::snapgene_preserve::packets(&encoded).unwrap();
    assert_eq!(packets.iter().filter(|(kind, _)| *kind == 10).count(), 1);
    assert_eq!(packets.iter().filter(|(kind, _)| *kind == 5).count(), 1);
    for (label, color, segment_colors) in expected {
        let (index, restored) = reopened
            .features
            .iter()
            .enumerate()
            .find(|(_, f)| f.label == label)
            .unwrap();
        assert_eq!(restored.effective_color_rgb(index), color, "{label}");
        assert_eq!(
            restored
                .display_segments()
                .iter()
                .map(|s| s
                    .color_rgb
                    .unwrap_or_else(|| restored.effective_color_rgb(index)))
                .collect::<Vec<_>>(),
            segment_colors,
            "{label}"
        );
    }
}
