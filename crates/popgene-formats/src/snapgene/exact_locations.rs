use super::*;
use popgene_core::translation::translate_cds;

fn range(start: usize, end: usize) -> FeatureLocation {
    FeatureLocation::Range {
        start,
        end,
        before: false,
        after: false,
    }
}

fn complement(location: FeatureLocation) -> FeatureLocation {
    FeatureLocation::Complement(Box::new(location))
}

fn feature(id: u64, location: FeatureLocation) -> Feature {
    let (parts, reverse) = location
        .uniform_projection()
        .unwrap_or_else(|| (location.local_parts(), false));
    Feature {
        id,
        label: format!("location {id}"),
        kind: "misc_feature".into(),
        parts,
        location: Some(location),
        reverse,
        qualifiers: vec![],
        edited: false,
    }
}

#[test]
fn extracted_complete_cds_keep_native_sequence_translation_and_named_colors() {
    let mut source = Document::new("template", "CCCATGAAACCCTAAGGGTTACCCTTTCATCCC").unwrap();
    for (index, (parts, reverse, note)) in [
        (
            vec![Segment { start: 3, end: 9 }, Segment { start: 12, end: 15 }],
            false,
            "This forward directional feature has 2 segments: 1: 4 .. 9 / #00aa66 / Start 2: 13 .. 15 / #9933cc / Stop",
        ),
        (
            vec![Segment { start: 18, end: 21 }, Segment { start: 24, end: 30 }],
            true,
            "This reverse directional feature has 2 segments: 1: 19 .. 21 / #ff6600 / Stop 2: 25 .. 30 / #3366cc / Start",
        ),
    ]
    .into_iter()
    .enumerate()
    {
        source.features.push(Feature {
            id: index as u64 + 1,
            label: format!("CDS {index}"),
            kind: "CDS".into(),
            parts,
            location: None,
            reverse,
            qualifiers: vec![
                ("note".into(), Some(note.into())),
                ("codon_start".into(), Some("1".into())),
                ("transl_table".into(), Some("11".into())),
            ],
            edited: false,
        });
    }
    let product = popgene_core::workflow::extract_region(&source, 3, 27, "complete CDSs").unwrap();
    assert!(product.features.iter().all(|f| f.location.is_some()));
    let encoded = write(&product).unwrap();
    let restored = read("product.dna", &encoded).unwrap();
    assert_eq!(restored.sequence, product.sequence);
    for (original, reopened) in product.features.iter().zip(&restored.features) {
        assert_eq!(product.feature_sequence(original).unwrap(), "ATGAAATAA");
        assert_eq!(restored.feature_sequence(reopened).unwrap(), "ATGAAATAA");
        assert_eq!(reopened.parts, original.parts);
        assert_eq!(reopened.reverse, original.reverse);
        assert_eq!(reopened.display_segments(), original.display_segments());
        let translation = translate_cds(&product, original).unwrap();
        assert_eq!(translation.amino_acids, b"MK*");
        assert_eq!(translate_cds(&restored, reopened).unwrap(), translation);
    }
    assert_eq!(write(&restored).unwrap(), encoded);

    // Cutting either CDS boundary remains a real partial location, not a reason
    // to reuse the exact native segment encoding.
    let clipped = popgene_core::workflow::extract_region(&source, 4, 25, "partial CDSs").unwrap();
    assert!(write(&clipped).is_err());
}

#[test]
fn exact_nested_reverse_overlapping_and_origin_joins_keep_order() {
    use FeatureLocation::Join;
    let cases = [
        (range(0, 3), "AAA", vec!["1-3"], false),
        (complement(range(0, 3)), "TTT", vec!["1-3"], true),
        (
            Join(vec![range(9, 12), Join(vec![range(0, 3), range(3, 6)])]),
            "TTTAAACCC",
            vec!["10-12", "1-3", "4-6"],
            false,
        ),
        (
            complement(Join(vec![range(0, 3), range(6, 9)])),
            "CCCTTT",
            vec!["1-3", "7-9"],
            true,
        ),
        (
            Join(vec![complement(range(0, 3)), complement(range(6, 9))]),
            "TTTCCC",
            vec!["7-9", "1-3"],
            true,
        ),
        (
            Join(vec![range(0, 6), range(3, 9)]),
            "AAACCCCCCGGG",
            vec!["1-6", "4-9"],
            false,
        ),
        (
            complement(Join(vec![range(9, 12), range(0, 6)])),
            "GGGTTTAAA",
            vec!["10-12", "1-6"],
            true,
        ),
        (
            Join(vec![range(0, 3), range(0, 3)]),
            "AAAAAA",
            vec!["1-3", "1-3"],
            false,
        ),
    ];
    for (location, expected, ranges, reverse) in cases {
        let mut doc = Document::new("exact location", "AAACCCGGGTTT").unwrap();
        doc.circular = true;
        doc.features.push(feature(1, location));
        assert_eq!(doc.feature_sequence(&doc.features[0]).unwrap(), expected);
        let encoded = write(&doc).unwrap();
        let packets = crate::snapgene_preserve::packets(&encoded).unwrap();
        let features = xml(packets.iter().find(|(kind, _)| *kind == 10).unwrap().1).unwrap();
        let node = &features.children[0];
        assert_eq!(node.attr("directionality"), if reverse { "2" } else { "1" });
        assert_eq!(
            node.children
                .iter()
                .filter(|node| node.name == "Segment")
                .map(|node| node.attr("range"))
                .collect::<Vec<_>>(),
            ranges
        );
        let restored = read("exact.dna", &encoded).unwrap();
        assert_eq!(restored.features[0].parts, doc.features[0].parts);
        assert_eq!(restored.features[0].reverse, reverse);
        assert_eq!(
            restored.feature_sequence(&restored.features[0]).unwrap(),
            expected
        );
    }
}

#[test]
fn exact_rna_and_protein_locations_preserve_their_alphabet() {
    for (molecule, sequence, location, expected) in [
        (
            Molecule::Rna,
            "AAACCCGGGUUU",
            complement(FeatureLocation::Join(vec![range(0, 3), range(6, 9)])),
            "CCCUUU",
        ),
        (
            Molecule::Protein,
            "MKWVTFISL",
            FeatureLocation::Join(vec![range(0, 3), range(6, 9)]),
            "MKWISL",
        ),
    ] {
        let mut doc = Document::new_typed("exact", sequence, molecule).unwrap();
        doc.features.push(feature(1, location));
        let restored = read("exact", &write(&doc).unwrap()).unwrap();
        assert_eq!(restored.molecule, molecule);
        assert_eq!(
            restored.feature_sequence(&restored.features[0]).unwrap(),
            expected
        );
    }
    // Even a cancelled pair of complements is not a valid protein operation.
    let mut protein = Document::new_typed("protein", "MKWVTFISL", Molecule::Protein).unwrap();
    protein
        .features
        .push(feature(1, complement(complement(range(0, 3)))));
    assert!(write(&protein).is_err());
}

#[test]
fn exact_locations_keep_full_duplex_axis_and_bottom_only_annotations() {
    let mut doc = Document::new("duplex extent", "TAACCCGGGTTT").unwrap();
    doc.features = vec![
        feature(1, range(0, 4)),
        feature(2, complement(range(4, 7))),
        feature(3, complement(range(10, 12))),
    ];
    ends::read(&mut doc, &[(8, b"<AdditionalSequenceProperties><UpstreamStickiness>-1</UpstreamStickiness><DownstreamStickiness>2</DownstreamStickiness><UpstreamModification>Unmodified</UpstreamModification><DownstreamModification>Unmodified</DownstreamModification></AdditionalSequenceProperties>")]).unwrap();
    assert_eq!(doc.features.len(), 1);
    assert_eq!(doc.duplex_annotations.features.len(), 2);
    let encoded = write(&doc).unwrap();
    let restored = read("duplex.dna", &encoded).unwrap();
    assert_eq!(restored.sequence, doc.sequence);
    assert_eq!(restored.duplex_ends, doc.duplex_ends);
    assert_eq!(restored.features[0].parts, doc.features[0].parts);
    assert_eq!(restored.features[0].reverse, doc.features[0].reverse);
    for (original, reopened) in doc
        .duplex_annotations
        .features
        .iter()
        .zip(&restored.duplex_annotations.features)
    {
        assert_eq!(reopened.parts, original.parts);
        assert_eq!(reopened.reverse, original.reverse);
    }
}

#[test]
fn unrepresented_locations_remain_refused_even_with_a_local_projection() {
    use FeatureLocation::*;
    let cases = [
        Range {
            start: 0,
            end: 3,
            before: true,
            after: false,
        },
        complement(Range {
            start: 0,
            end: 3,
            before: false,
            after: true,
        }),
        // A fuzzy interior child is invisible to outer partial_ends(), but
        // still cannot become an exact native segment.
        Join(vec![
            range(0, 3),
            Range {
                start: 3,
                end: 6,
                before: true,
                after: false,
            },
            range(6, 9),
        ]),
        Order(vec![range(0, 3), range(6, 9)]),
        Bond(vec![range(0, 1), range(6, 7)]),
        OneOf(vec![range(0, 3), range(6, 9)]),
        Remote {
            accession: "OTHER.1".into(),
            location: Some(Box::new(range(0, 3))),
        },
        Join(vec![
            range(0, 3),
            Gap {
                length: Some(3),
                unknown100: false,
            },
            range(6, 9),
        ]),
        Join(vec![range(0, 3), complement(range(6, 9))]),
        Join(vec![range(0, 3), Between { left: 5, right: 6 }]),
    ];
    for location in cases {
        let mut doc = Document::new("unsupported", "AAACCCGGGTTT").unwrap();
        doc.features.push(feature(1, location));
        doc.validate().unwrap();
        let error = write(&doc).unwrap_err().to_string();
        assert!(error.contains("without changing its meaning"), "{error}");
    }
}
