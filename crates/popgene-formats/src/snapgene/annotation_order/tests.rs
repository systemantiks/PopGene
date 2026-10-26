use super::super::{Document, Feature, Segment, packet, read, write, write_warnings, xml};

fn feature(id: u64, primer: bool, start: usize) -> Feature {
    Feature {
        id,
        label: format!("annotation {id}"),
        kind: if primer {
            "primer_bind"
        } else {
            "misc_feature"
        }
        .into(),
        parts: vec![Segment {
            start,
            end: start + 4,
        }],
        reverse: false,
        location: None,
        qualifiers: if primer {
            vec![("primer_sequence".into(), Some("ACGT".into()))]
        } else {
            vec![]
        },
        edited: false,
    }
}
fn warnings(doc: &Document) -> Vec<String> {
    write_warnings(doc)
        .into_iter()
        .filter(|warning| warning.contains("implicit display color"))
        .collect()
}
fn assert_no_palette_warning(doc: &Document) {
    let messages = warnings(doc);
    assert!(messages.is_empty(), "{messages:?}");
}
fn assert_no_encoded_primer_color(bytes: &[u8]) {
    let packets = crate::snapgene_preserve::packets(bytes).unwrap();
    let mut count = 0;
    for (_, body) in packets.iter().filter(|(kind, _)| *kind == 5) {
        for primer in xml(body)
            .unwrap()
            .children
            .iter()
            .filter(|node| node.name == "Primer")
        {
            assert!(!primer.attrs.contains_key("color"));
            count += 1;
        }
    }
    assert!(count > 0);
}
fn all_features(doc: &Document) -> impl Iterator<Item = &Feature> {
    doc.features.iter().chain(&doc.duplex_annotations.features)
}
/// Independent observed result: compare actual import colors, matching these
/// synthetic fixtures by unique location/strand (not rewritten native IDs).
fn changed_uncolored_primers(before: &Document, after: &Document) -> usize {
    all_features(before)
        .enumerate()
        .filter(|(_, feature)| {
            feature.kind == "primer_bind" && feature.display_color_rgb().is_none()
        })
        .filter(|(old_index, feature)| {
            let (new_index, restored) = all_features(after)
                .enumerate()
                .find(|(_, restored)| {
                    restored.kind == feature.kind
                        && restored.parts == feature.parts
                        && restored.reverse == feature.reverse
                })
                .unwrap();
            feature.effective_color_rgb(*old_index) != restored.effective_color_rgb(new_index)
        })
        .count()
}
fn reopen(doc: &Document) -> Document {
    read("reopened.dna", &write(doc).unwrap()).unwrap()
}
fn source(duplex: bool) -> Document {
    let mut bytes = vec![];
    packet(&mut bytes, 9, b"SnapGene\0\x01\0\x0f\0\x14").unwrap();
    packet(&mut bytes, 0, b"\x02ACGTACGTACGTACGTACGTACGTACGTACGT").unwrap();
    // The first record has three stored sites but two distinct imported models.
    packet(&mut bytes, 5, b"<Primers><Primer recentID=\"7\" name=\"multi\" sequence=\"ACGT\" chemistry=\"untouched\"><BindingSite location=\"0-3\" boundStrand=\"0\"/><BindingSite location=\"8-11\" boundStrand=\"0\"/><BindingSite location=\"0-3\" boundStrand=\"0\"/></Primer></Primers>").unwrap();
    packet(&mut bytes, 10, b"<Features><Feature recentID=\"3\" name=\"neighbor\" type=\"misc_feature\"><Segment range=\"17-20\" color=\"#aabbcc\"/></Feature></Features>").unwrap();
    if duplex {
        packet(&mut bytes, 8, b"<AdditionalSequenceProperties><UpstreamStickiness>-4</UpstreamStickiness><DownstreamStickiness>0</DownstreamStickiness></AdditionalSequenceProperties>").unwrap();
    }
    read("original.dna", &bytes).unwrap()
}

#[test]
fn interleaved_primer_retains_default_after_native_regrouping_without_color_attribute() {
    let mut doc = Document::new("palette", &"ACGT".repeat(8)).unwrap();
    doc.features = vec![
        feature(1, false, 0),
        feature(2, true, 4),
        feature(3, false, 8),
    ];
    assert_no_palette_warning(&doc);
    let bytes = write(&doc).unwrap();
    let reopened = read("regrouped.dna", &bytes).unwrap();
    assert_eq!(changed_uncolored_primers(&doc, &reopened), 0);
    assert_eq!(doc.features[1].effective_color_rgb(1), 0xb985d0);
    assert_eq!(reopened.features[2].effective_color_rgb(2), 0xb985d0);
    assert!(reopened.features[2].qualifier("snapgene_color").is_none());
    assert_no_encoded_primer_color(&bytes);
}

#[test]
fn primer_default_and_explicit_colors_remain_stable_across_different_index_shifts() {
    let mut doc = Document::new("cycle", &"ACGT".repeat(8)).unwrap();
    doc.features = std::iter::once(feature(1, true, 0))
        .chain((2..=7).map(|id| feature(id, false, (id as usize - 1) * 4)))
        .collect();
    assert_no_palette_warning(&doc);
    assert_eq!(changed_uncolored_primers(&doc, &reopen(&doc)), 0);
    doc.features.pop(); // A five-place move also leaves the primer default stable.
    assert_no_palette_warning(&doc);
    assert_eq!(changed_uncolored_primers(&doc, &reopen(&doc)), 0);
    for color in [
        "black",
        "orange",
        "blue",
        "#123456",
        "unverified-native-color",
    ] {
        doc.features[0]
            .qualifiers
            .push(("snapgene_color".into(), Some(color.into())));
        assert_no_palette_warning(&doc); // unknown values remain uninterpreted
        assert_eq!(
            reopen(&doc)
                .features
                .last()
                .unwrap()
                .qualifier("snapgene_color"),
            Some(color)
        );
        doc.features[0].qualifiers.pop();
    }
    doc.features[0]
        .qualifiers
        .push(("snapgene_color".into(), Some(String::new())));
    assert_no_palette_warning(&doc); // empty color remains unrecorded for display
    assert_eq!(changed_uncolored_primers(&doc, &reopen(&doc)), 0);
    assert_eq!(
        reopen(&doc)
            .features
            .last()
            .unwrap()
            .qualifier("snapgene_color"),
        Some("")
    );
    doc.features[0].qualifiers.pop();
    doc.features[0]
        .qualifiers
        .push(("note".into(), Some("reaction details retained".into())));
    assert_no_palette_warning(&doc); // additional metadata routes it through Feature
    let written = write(&doc).unwrap();
    assert!(
        !crate::snapgene_preserve::packets(&written)
            .unwrap()
            .iter()
            .any(|(kind, _)| *kind == 5)
    );
    assert_eq!(
        changed_uncolored_primers(&doc, &read("features.dna", &written).unwrap()),
        0
    );
    doc.features = vec![feature(1, false, 0), feature(2, true, 4)];
    assert_no_palette_warning(&doc); // already grouped
}

#[test]
fn unchanged_and_metadata_only_saves_keep_packet_order_without_warning() {
    let mut doc = source(false);
    let original = doc.source.as_ref().unwrap().bytes.to_vec();
    assert!(write_warnings(&doc).is_empty());
    assert_eq!(write(&doc).unwrap(), original);
    doc.description = "metadata edit only".into();
    assert_no_palette_warning(&doc);
    let updated = write(&doc).unwrap();
    let annotation_packets = |bytes: &[u8]| {
        crate::snapgene_preserve::packets(bytes)
            .unwrap()
            .into_iter()
            .filter(|(kind, _)| matches!(kind, 5 | 10))
            .map(|(kind, body)| (kind, body.to_vec()))
            .collect::<Vec<_>>()
    };
    assert_eq!(annotation_packets(&original), annotation_packets(&updated));
    assert_eq!(
        changed_uncolored_primers(&doc, &read("metadata.dna", &updated).unwrap()),
        0
    );
}

#[test]
fn preserved_multisite_primer_uses_source_site_order_after_model_reordering() {
    let mut doc = source(false);
    doc.features.swap(0, 1);
    doc.features[2].label = "changed neighbor".into();
    // Reopened order is neighbor, original first site, original second site.
    // Both binding-site models keep one default despite restored source order.
    assert_no_palette_warning(&doc);
    let bytes = write(&doc).unwrap();
    let reopened = read("multisite.dna", &bytes).unwrap();
    assert_eq!(changed_uncolored_primers(&doc, &reopened), 0);
    assert_no_encoded_primer_color(&bytes);
    assert_eq!(reopened.features[1].parts[0].start, 0);
    assert_eq!(reopened.features[2].parts[0].start, 8);
    let blocks = crate::snapgene_preserve::packets(&bytes).unwrap();
    let primer = xml(blocks.iter().find(|(kind, _)| *kind == 5).unwrap().1).unwrap();
    assert_eq!(primer.children.len(), 1);
    assert_eq!(primer.children[0].children.len(), 3);
    assert_eq!(primer.children[0].attr("chemistry"), "untouched");
}

#[test]
fn duplex_partition_preserves_multisite_associations_and_effective_indices() {
    let mut doc = source(true);
    assert_eq!(doc.features.len(), 2);
    assert_eq!(doc.duplex_annotations.features.len(), 1);
    doc.features[1].label = "changed top neighbor".into();
    assert_no_palette_warning(&doc);
    let bytes = write(&doc).unwrap();
    let reopened = read("duplex.dna", &bytes).unwrap();
    assert_eq!(changed_uncolored_primers(&doc, &reopened), 0);
    assert_no_encoded_primer_color(&bytes);
    assert_eq!(reopened.features[1].kind, "primer_bind");
    assert_eq!(reopened.duplex_annotations.features.len(), 1);
    assert_eq!(reopened.duplex_annotations.features[0].parts[0].start, 0);
}

#[test]
fn duplex_partition_keeps_both_top_and_opposite_only_primer_defaults() {
    let mut doc = source(true);
    doc.source = None;
    doc.features = vec![feature(1, true, 4)];
    doc.duplex_annotations.features = vec![feature(2, false, 0)];
    // Writer emits the overhang Feature first, but reader partitions it after
    // the top Primer again. The displayed Primer remains at index zero.
    assert_no_palette_warning(&doc);
    assert_eq!(changed_uncolored_primers(&doc, &reopen(&doc)), 0);
    doc.features.push(feature(3, false, 12));
    doc.duplex_annotations.features = vec![feature(4, true, 0), feature(5, false, 0)];
    assert_no_palette_warning(&doc);
    let bytes = write(&doc).unwrap();
    let reopened = read("duplex-reordered.dna", &bytes).unwrap();
    assert_no_encoded_primer_color(&bytes);
    assert_eq!(changed_uncolored_primers(&doc, &reopened), 0);
    assert_eq!(reopened.features[1].effective_color_rgb(1), 0xb985d0);
    assert_eq!(
        reopened.duplex_annotations.features[1].effective_color_rgb(3),
        0xb985d0
    );
}
