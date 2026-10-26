use super::*;
use popgene_core::{FeatureLocation as L, workflow};

fn annotation(id: u64, label: &str, start: usize, end: usize) -> Feature {
    Feature {
        id,
        label: label.into(),
        kind: "misc_feature".into(),
        parts: vec![Segment { start, end }],
        reverse: false,
        location: None,
        qualifiers: vec![("label".into(), Some(label.into()))],
        edited: false,
    }
}
fn tracked() -> Document {
    let mut source = Document::new("ancestry", "ATGGAATTCAAGTAA").unwrap();
    let mut feature = annotation(17, "α coding region", 0, 15);
    feature.kind = "CDS".into();
    feature
        .qualifiers
        .push(("codon_start".into(), Some("1".into())));
    source.features.push(feature);
    workflow::extract_region(&source, 0, 15, "fragment").unwrap()
}
fn roundtrip(doc: &Document, format: SaveFormat) -> Document {
    let bytes = crate::encode(doc, format)
        .unwrap_or_else(|error| panic!("{} roundtrip for {}: {error}", format.label(), doc.name));
    if let Some(directory) = std::env::var_os("POPGENE_LINEAGE_QA_DIR") {
        use std::io::Write as _;
        let directory = std::path::PathBuf::from(directory);
        std::fs::create_dir_all(&directory).unwrap();
        let name = doc
            .name
            .chars()
            .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
            .collect::<String>();
        let path = directory.join(format!("{name}.{}", format.extension_for(doc)));
        std::fs::OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(path)
            .unwrap()
            .write_all(&bytes)
            .unwrap();
    }
    let reread = crate::import(&format!("case.{}", format.extension_for(doc)), &bytes).unwrap();
    assert!(
        !reread.annotation_lineage.is_empty(),
        "{:?}",
        reread.warnings
    );
    assert_eq!(
        reread.annotation_lineage.instance_id,
        doc.annotation_lineage.instance_id
    );
    assert_eq!(
        reread.annotation_lineage.origins,
        doc.annotation_lineage.origins,
        "{}: warnings={:?}; before={:?}; after={:?}",
        format.label(),
        crate::save_warnings(doc, format),
        doc.features,
        reread.features
    );
    assert_eq!(
        reread.annotation_lineage.pieces.len(),
        doc.annotation_lineage.pieces.len()
    );
    for piece in &reread.annotation_lineage.pieces {
        assert!(reread.annotation_lineage.actionable_piece(&reread, piece));
    }
    assert_eq!(crate::encode(&reread, format).unwrap(), bytes);
    reread
}

#[test]
fn ancestry_survives_native_genbank_embl_canonical_save_and_exact_reencode() {
    let doc = tracked();
    assert_eq!(doc.annotation_lineage.pieces.len(), 1);
    for format in [SaveFormat::SnapGene, SaveFormat::GenBank, SaveFormat::Embl] {
        let reread = roundtrip(&doc, format);
        assert_eq!(reread.sequence, doc.sequence);
        assert_eq!(
            reread.feature_sequence(&reread.features[0]).unwrap(),
            "ATGGAATTCAAGTAA"
        );
        assert!(reread.comments.iter().all(|c| !c.contains(START_PREFIX)));
    }
}

#[test]
fn context_golden_matches_nested_desktop_record_without_map_feature_dependence() {
    // Synthetic reconstruction of the exact typed context independently checked
    // from a source0016 desktop Save/reopen. Python compact JSON produced both
    // golden hashes; nested Feature/Reference ordering matters for legacy reads.
    let mut doc = Document::new("unused name", "ATGGAATTCAAGTAA").unwrap();
    let mut feature = annotation(1, "coding", 0, 15);
    feature.kind = "CDS".into();
    feature.qualifiers = [
        ("snapgene_color", "#2255cc"),
        ("codon_start", "2"),
        ("label", "coding"),
        ("note", "color: #2255cc; direction: RIGHT"),
        ("transl_table", "1"),
        ("snapgene_directionality", "1"),
        ("snapgene_readingFrame", "2"),
    ]
    .into_iter()
    .map(|(key, value)| (key.into(), Some(value.into())))
    .collect();
    doc.features.push(feature);
    doc.references.push(popgene_core::metadata::Reference {
        title: "Direct Submission".into(),
        authors: Some(".".into()),
        journal: Some("PopGene export using SnapGene GenBank format".into()),
        extra: vec![("type".into(), "Journal Article".into())],
        ..Default::default()
    });
    let hashes = context(&doc).unwrap();
    assert_eq!(
        hashes.canonical,
        "86c7eef66fe86fd5b7c824016793e91172a5abf4fb39f647e10df827f2c4f4bf"
    );
    assert_eq!(
        hashes.legacy_insertion(),
        "57d566e0387786fd788c808d2c02e4c6f4e1dd2d67e2923a3d372b970af2dd61"
    );
}

#[test]
fn canonical_and_exact_legacy_contexts_reopen_but_changed_context_stays_inactive() {
    let doc = tracked();
    for format in [SaveFormat::SnapGene, SaveFormat::GenBank, SaveFormat::Embl] {
        let plain: fn(&Document) -> Result<Vec<u8>> = match format {
            SaveFormat::SnapGene => crate::snapgene::write_plain,
            SaveFormat::GenBank => crate::genbank_bytes_plain,
            SaveFormat::Embl => crate::embl::write_plain,
            _ => unreachable!(),
        };
        let first = plain(&doc).unwrap();
        let canonical = crate::import(&format!("context.{}", format.extension()), &first).unwrap();
        let registry = projected_registry(&doc, &canonical, format).unwrap();
        let hashes = context(&canonical).unwrap();
        assert_ne!(hashes.canonical, hashes.legacy_insertion());
        for expected in [
            &hashes.canonical,
            hashes.legacy_insertion(),
            &"0".repeat(64),
        ] {
            let valid = expected != "0".repeat(64);
            let value = serde_json::json!({
                "format": format_name(format),
                "records_sha256": records_hash(&first, format).unwrap(),
                "context_sha256": expected,
                "registry": registry,
            });
            let mut with_carrier = doc.clone();
            with_carrier.comments.push(block(&value).unwrap());
            let bytes = plain(&with_carrier).unwrap();
            assert_eq!(
                records_hash(&first, format).unwrap(),
                records_hash(&bytes, format).unwrap()
            );
            let reopened =
                crate::import(&format!("legacy.{}", format.extension()), &bytes).unwrap();
            if valid {
                assert_eq!(reopened.annotation_lineage, registry);
                assert!(
                    reopened
                        .annotation_lineage
                        .actionable_piece(&reopened, &reopened.annotation_lineage.pieces[0])
                );
                assert_eq!(crate::encode(&reopened, format).unwrap(), bytes);
            } else {
                assert!(reopened.annotation_lineage.is_empty());
                assert!(reopened.comments.iter().any(|c| c.contains(START_PREFIX)));
                assert!(
                    reopened
                        .warnings
                        .iter()
                        .any(|w| w.contains("canonical context changed"))
                );
            }
        }
    }
}

#[test]
fn actual_emission_order_binds_native_primers_without_converting_them() {
    let mut source = Document::new("primer first", "ACGTACGTACGT").unwrap();
    let mut primer = annotation(91, "bound primer", 0, 4);
    primer.kind = "primer_bind".into();
    primer.qualifiers.extend([
        ("primer_sequence".into(), Some("ACGT".into())),
        ("snapgene_color".into(), Some("blue".into())),
    ]);
    source.features = vec![primer, annotation(72, "neighbor", 5, 8)];
    let doc = workflow::extract_region(&source, 0, 12, "all").unwrap();
    let bytes = crate::encode(&doc, SaveFormat::SnapGene).unwrap();
    let packets = crate::snapgene_preserve::packets(&bytes).unwrap();
    assert_eq!(packets.iter().filter(|(kind, _)| *kind == 5).count(), 1);
    let reread = roundtrip(&doc, SaveFormat::SnapGene);
    assert_eq!(reread.features[0].label, "neighbor");
    assert_eq!(reread.features[1].label, "bound primer");
    let old = doc
        .annotation_lineage
        .pieces
        .iter()
        .find(|p| p.feature_id == 91)
        .unwrap();
    let new = reread
        .annotation_lineage
        .pieces
        .iter()
        .find(|p| p.feature_id == 2)
        .unwrap();
    assert_eq!(old.origin_id, new.origin_id);
    assert_eq!(old.runs, new.runs);
}

#[test]
fn duplicate_features_and_repeated_join_leaves_remain_distinct() {
    let mut source = Document::new("duplicates", "ACGTACGTACGT").unwrap();
    let mut a = annotation(2, "same", 0, 4);
    a.parts.push(Segment { start: 0, end: 4 });
    a.location = Some(L::Join(vec![
        L::Range {
            start: 0,
            end: 4,
            before: false,
            after: false,
        },
        L::Range {
            start: 0,
            end: 4,
            before: false,
            after: false,
        },
    ]));
    a.qualifiers
        .push(("snapgene_color".into(), Some("red".into())));
    a.qualifiers
        .push(("snapgene_directionality".into(), Some("1".into())));
    let mut b = a.clone();
    b.id = 8;
    source.features = vec![a, b];
    let doc = workflow::extract_region(&source, 0, 12, "duplicates").unwrap();
    assert_ne!(
        doc.annotation_lineage.pieces[0].origin_id,
        doc.annotation_lineage.pieces[1].origin_id
    );
    for format in [SaveFormat::SnapGene, SaveFormat::GenBank, SaveFormat::Embl] {
        let reread = roundtrip(&doc, format);
        assert_eq!(reread.features.len(), 2);
        assert_eq!(
            reread.feature_sequence(&reread.features[0]).unwrap(),
            "ACGTACGT"
        );
        let leaves = reread.annotation_lineage.pieces[0]
            .runs
            .iter()
            .map(|r| r.leaf)
            .collect::<BTreeSet<_>>();
        assert_eq!(leaves, BTreeSet::from([0, 1]));
    }
}

#[test]
fn changed_qualifiers_are_not_laundered_and_fasta_explicitly_drops_ancestry() {
    let mut doc = tracked();
    doc.features[0]
        .qualifiers
        .push(("note".into(), Some("user changed this annotation".into())));
    assert!(
        !doc.annotation_lineage
            .actionable_piece(&doc, &doc.annotation_lineage.pieces[0])
    );
    assert!(
        crate::save_warnings(&doc, SaveFormat::SnapGene)
            .iter()
            .any(|w| w.contains("not carried as active"))
    );
    let bytes = crate::encode(&doc, SaveFormat::SnapGene).unwrap();
    let reread = crate::import("inactive.dna", &bytes).unwrap();
    assert!(reread.annotation_lineage.pieces.is_empty());
    assert!(reread.annotation_lineage.tracking_limited);
    assert_eq!(
        reread.annotation_lineage.instance_id,
        doc.annotation_lineage.instance_id
    );
    let mut renamed = reread.clone();
    renamed.name = "new name".into();
    let encoded = crate::encode(&renamed, SaveFormat::SnapGene).unwrap();
    assert_eq!(
        crate::import("seed.dna", &encoded)
            .unwrap()
            .annotation_lineage
            .instance_id,
        doc.annotation_lineage.instance_id
    );
    let tracked = tracked();
    assert!(
        crate::save_warnings(&tracked, SaveFormat::Fasta)
            .iter()
            .any(|w| w.contains("ancestry is omitted"))
    );
    let fasta = crate::encode(&tracked, SaveFormat::Fasta).unwrap();
    assert!(
        crate::import("lost.fasta", &fasta)
            .unwrap()
            .annotation_lineage
            .is_empty()
    );
}

#[test]
fn unknown_or_conflicting_display_values_are_not_laundered_by_a_valid_value() {
    let cases = [
        vec![
            ("snapgene_color", Some("nonsense")),
            ("note", Some("color: #ff0000; direction: RIGHT")),
        ],
        vec![
            ("snapgene_color", Some("red")),
            ("snapgene_color", Some("blue")),
        ],
        vec![
            ("snapgene_color", Some("#ff0000")),
            ("snapgene_color", Some("nonsense")),
        ],
        vec![("snapgene_color", Some("red")), ("snapgene_color", None)],
        vec![
            ("snapgene_directionality", Some("1")),
            ("snapgene_directionality", Some("2")),
        ],
        vec![
            ("snapgene_directionality", Some("1")),
            ("snapgene_directionality", Some("nonsense")),
        ],
        vec![
            ("snapgene_directionality", Some("1")),
            ("snapgene_directionality", None),
        ],
    ];
    for values in cases {
        let mut source = Document::new("display metadata", "ACGTACGT").unwrap();
        let mut feature = annotation(1, "display", 0, 8);
        feature.qualifiers.extend(
            values
                .iter()
                .map(|(key, value)| ((*key).into(), value.map(str::to_owned))),
        );
        source.features.push(feature);
        let doc = workflow::extract_region(&source, 0, 8, "display fragment").unwrap();
        assert!(
            doc.annotation_lineage
                .actionable_piece(&doc, &doc.annotation_lineage.pieces[0])
        );
        let bytes = crate::encode(&doc, SaveFormat::SnapGene).unwrap();
        let reread = crate::import("display.dna", &bytes).unwrap();
        assert_eq!(reread.sequence, doc.sequence);
        assert_eq!(reread.features[0].parts, doc.features[0].parts);
        assert_eq!(
            reread.features[0].effective_color_rgb(0),
            doc.features[0].effective_color_rgb(0)
        );
        assert!(reread.annotation_lineage.pieces.is_empty(), "{values:?}");
        assert!(reread.annotation_lineage.tracking_limited, "{values:?}");
        let warnings = crate::save_warnings(&doc, SaveFormat::SnapGene).join(" ");
        assert!(
            warnings.contains("Active annotation ancestry is omitted")
                && warnings.contains("qualifiers")
                && warnings.contains("ancestry-limited marker"),
            "{values:?}: {warnings}"
        );
    }
}

#[test]
fn individually_equivalent_display_values_keep_actionable_ancestry() {
    let mut source = Document::new("equivalent display", "ACGTACGT").unwrap();
    let mut feature = annotation(1, "display", 0, 8);
    feature.qualifiers.extend([
        ("snapgene_color".into(), Some("red".into())),
        ("snapgene_color".into(), Some("#FF0000".into())),
        ("snapgene_directionality".into(), Some("1".into())),
        ("snapgene_directionality".into(), Some("1".into())),
    ]);
    source.features.push(feature);
    let doc = workflow::extract_region(&source, 0, 8, "equivalent display").unwrap();
    let bytes = crate::encode(&doc, SaveFormat::SnapGene).unwrap();
    let reread = crate::import("equivalent.dna", &bytes).unwrap();
    assert_eq!(reread.annotation_lineage.pieces.len(), 1);
    assert!(
        reread
            .annotation_lineage
            .actionable_piece(&reread, &reread.annotation_lineage.pieces[0])
    );
    assert_eq!(reread.features[0].display_color_rgb(), Some(0xff0000));
    assert_eq!(
        reread.features[0].display_direction(),
        Some(FeatureDirection::Forward)
    );
    assert!(
        !crate::save_warnings(&doc, SaveFormat::SnapGene)
            .iter()
            .any(|warning| warning.contains("ancestry is omitted"))
    );
}

#[test]
fn corrupted_unknown_and_changed_record_carriers_remain_readable_inactive() {
    let doc = tracked();
    let bytes = crate::encode(&doc, SaveFormat::GenBank).unwrap();
    let text = String::from_utf8(bytes).unwrap();
    for modified in [
        text.replacen("v1 begin", "v2 begin", 1),
        text.replacen("/codon_start=1", "/codon_start=2", 1),
        text.replacen(";sha256=", ";sha256=f", 1),
    ] {
        assert_ne!(modified, text);
        let reread = crate::import("inactive.gb", modified.as_bytes()).unwrap();
        assert!(reread.annotation_lineage.is_empty());
        assert!(reread.comments.iter().any(|c| c.contains(START_PREFIX)));
        assert!(
            reread
                .warnings
                .iter()
                .any(|w| w.contains("no ancestry was activated"))
        );
    }
}

#[test]
fn split_empty_seed_cannot_activate_after_ordinary_comment_folding() {
    let source = Document::new("empty seed", "AACCGGTT").unwrap();
    for format in [SaveFormat::GenBank, SaveFormat::Embl] {
        let plain: fn(&Document) -> Result<Vec<u8>> = match format {
            SaveFormat::GenBank => crate::genbank_bytes_plain,
            SaveFormat::Embl => crate::embl::write_plain,
            _ => unreachable!(),
        };
        let bytes = plain(&source).unwrap();
        let canonical = crate::import(&format!("seed.{}", format.extension()), &bytes).unwrap();
        let seeded = write_registry(
            &source,
            &canonical,
            &bytes,
            format,
            plain,
            &source.annotation_lineage,
        )
        .unwrap();
        let mut inactive = if format == SaveFormat::GenBank {
            crate::genbank(&seeded).unwrap()
        } else {
            crate::embl::read(std::str::from_utf8(&seeded).unwrap()).unwrap()
        };
        for comment in &mut inactive.comments {
            *comment = comment
                .replace("lineage v1", "lineage\nv1")
                .replace("lineage end]", "lineage\nend]");
        }
        inactive
            .comments
            .push("Unrelated [PopGene example prose".into());
        assert!(inactive.annotation_lineage.is_empty());
        assert_eq!(marker_candidates(&inactive.comments, START_PREFIX, END), 2);
        if format == SaveFormat::Embl {
            let error = crate::encode(&inactive, format).unwrap_err().to_string();
            assert!(error.contains("previously inactive") && error.contains("would become active"));
        } else {
            let bytes = crate::encode(&inactive, format).unwrap();
            let reread = crate::import("still-inactive.gb", &bytes).unwrap();
            assert!(reread.annotation_lineage.is_empty());
            assert_ne!(
                reread.annotation_lineage.instance_id,
                source.annotation_lineage.instance_id
            );
            assert_eq!(marker_candidates(&reread.comments, START_PREFIX, END), 2);
        }
    }
}

#[test]
fn full_duplex_carrier_retains_axis_and_does_not_replace_lost_geometry_or_fuzziness() {
    use popgene_core::restriction::{EndKind, FragmentEnd};
    let end = |kind, overhang: &str| FragmentEnd {
        kind,
        overhang: overhang.into(),
        enzymes: vec![],
    };
    let mut source = Document::new("physical", "AACCCGGGT").unwrap();
    source.duplex_ends = Some(workflow::DuplexEnds {
        left: end(EndKind::ThreePrime, "A"),
        right: end(EndKind::FivePrime, "AA"),
        opposite_sequence: "AAACCCGGGTTA".into(),
        left_modification: None,
        right_modification: None,
    });
    source.features.push(annotation(21, "top", 3, 6));
    let mut overhang = annotation(9, "bottom only", 0, 1);
    overhang.reverse = true;
    source.duplex_annotations.features.push(overhang);
    source.validate().unwrap();
    let mut doc = workflow::duplex_reference(&source).unwrap();
    doc.sequence = source.sequence.clone();
    doc.duplex_ends = source.duplex_ends.clone();
    workflow::partition_duplex_annotations(&mut doc).unwrap();
    let reread = roundtrip(&doc, SaveFormat::SnapGene);
    assert_eq!(reread.duplex_annotations.features.len(), 1);
    assert!(
        reread
            .annotation_lineage
            .pieces
            .iter()
            .any(|p| p.axis == Axis::Duplex)
    );
    let warnings = crate::save_warnings(&doc, SaveFormat::GenBank).join(" ");
    assert!(
        warnings.contains("readable PopGene comment extension"),
        "{warnings}"
    );
    let converted = crate::import(
        "ordinary.gb",
        &crate::encode(&doc, SaveFormat::GenBank).unwrap(),
    )
    .unwrap();
    assert_eq!(converted.sequence, doc.sequence);
    assert_eq!(converted.duplex_ends, doc.duplex_ends);
    assert_eq!(
        converted.annotation_lineage.pieces.len(),
        doc.annotation_lineage.pieces.len()
    );
    assert!(!converted.annotation_lineage.tracking_limited);
    let partial = workflow::extract_region(&tracked(), 0, 8, "partial").unwrap();
    let error = crate::encode(&partial, SaveFormat::SnapGene)
        .unwrap_err()
        .to_string();
    assert!(error.contains("fuzzy"), "{error}");
    roundtrip(&partial, SaveFormat::GenBank);
}

#[test]
fn payload_limit_keeps_ordinary_native_records_with_explicit_tracking_loss() {
    let mut source = Document::new("large annotation metadata", "ACGT").unwrap();
    let mut feature = annotation(1, "large metadata", 0, 4);
    let note = "metadata".repeat(270_000);
    feature.qualifiers.push(("note".into(), Some(note.clone())));
    source.features.push(feature);
    let doc = workflow::extract_region(&source, 0, 4, "bounded carrier").unwrap();
    assert!(!doc.annotation_lineage.pieces.is_empty());
    assert!(doc.annotation_lineage.snapshot_bytes() < popgene_core::lineage::MAX_BYTES);
    let warnings = crate::save_warnings(&doc, SaveFormat::SnapGene).join(" ");
    assert!(
        warnings.contains("Active annotation ancestry is omitted") && warnings.contains("4 MiB")
    );
    let bytes = crate::encode(&doc, SaveFormat::SnapGene).unwrap();
    let reread = crate::import("bounded.dna", &bytes).unwrap();
    assert_eq!(reread.sequence.as_ref(), "ACGT");
    assert_eq!(reread.features[0].qualifier("note"), Some(note.as_str()));
    assert!(reread.annotation_lineage.pieces.is_empty());
    assert!(reread.annotation_lineage.tracking_limited);
}

#[test]
fn native_notes_capacity_keeps_ordinary_output_when_even_limited_marker_cannot_fit() {
    const XML_LIMIT: usize = 16 * 1024 * 1024;
    let mut doc = tracked();
    doc.comments = vec!["x".into()];
    let probe = crate::snapgene::write_plain(&doc).unwrap();
    let overhead = crate::snapgene_preserve::packets(&probe)
        .unwrap()
        .into_iter()
        .find(|(kind, _)| *kind == 6)
        .unwrap()
        .1
        .len()
        - 1;
    doc.comments[0] = "x".repeat(XML_LIMIT - overhead - 16);
    let ordinary = crate::snapgene::write_plain(&doc).unwrap();
    assert_eq!(
        crate::snapgene_preserve::packets(&ordinary)
            .unwrap()
            .into_iter()
            .find(|(kind, _)| *kind == 6)
            .unwrap()
            .1
            .len(),
        XML_LIMIT - 16
    );
    let bytes = crate::encode(&doc, SaveFormat::SnapGene).unwrap();
    assert_eq!(bytes, ordinary);
    let reread = crate::import("no-marker.dna", &bytes).unwrap();
    assert_eq!(reread.sequence, doc.sequence);
    assert_eq!(reread.comments, doc.comments);
    assert_eq!(reread.features[0].parts, doc.features[0].parts);
    assert!(reread.annotation_lineage.is_empty());
    assert!(!reread.annotation_lineage.tracking_limited);
    assert_ne!(
        reread.annotation_lineage.instance_id,
        doc.annotation_lineage.instance_id
    );
    let warnings = crate::save_warnings(&doc, SaveFormat::SnapGene).join(" ");
    assert!(
        warnings.contains("no ancestry marker could be stored"),
        "{warnings}"
    );
    assert!(warnings.contains("16 MiB"), "{warnings}");
    assert!(
        warnings.contains("cannot retain this tracking limitation or the source identity"),
        "{warnings}"
    );
    assert!(
        !warnings.contains("with an explicit ancestry-limited marker"),
        "{warnings}"
    );
}

#[test]
fn external_reorder_of_identical_native_models_is_not_reassociated() {
    let mut source = Document::new("same records", "ACGTACGT").unwrap();
    let mut a = annotation(1, "identical", 0, 4);
    a.qualifiers
        .push(("snapgene_color".into(), Some("red".into())));
    let mut b = a.clone();
    b.id = 2;
    source.features = vec![a, b];
    let doc = workflow::extract_region(&source, 0, 8, "identical records").unwrap();
    let bytes = crate::encode(&doc, SaveFormat::SnapGene).unwrap();
    let mut altered = Vec::new();
    for (kind, body) in crate::snapgene_preserve::packets(&bytes).unwrap() {
        let changed;
        let body = if kind == 10 {
            let text = std::str::from_utf8(body).unwrap();
            let first = text.find("<Feature ").unwrap();
            let first_end = first + text[first..].find("</Feature>").unwrap() + "</Feature>".len();
            let second = first_end + text[first_end..].find("<Feature ").unwrap();
            let second_end =
                second + text[second..].find("</Feature>").unwrap() + "</Feature>".len();
            changed = format!(
                "{}{}{}{}{}",
                &text[..first],
                &text[second..second_end],
                &text[first_end..second],
                &text[first..first_end],
                &text[second_end..]
            );
            assert_ne!(changed, text);
            changed.as_bytes()
        } else {
            body
        };
        altered.push(kind);
        altered.extend_from_slice(&(body.len() as u32).to_be_bytes());
        altered.extend_from_slice(body);
    }
    let reread = crate::import("reordered.dna", &altered).unwrap();
    assert_eq!(reread.features[0].label, reread.features[1].label);
    assert_eq!(reread.features[0].parts, reread.features[1].parts);
    assert!(reread.annotation_lineage.is_empty());
    assert!(reread.comments.iter().any(|c| c.contains(START_PREFIX)));
    assert!(
        reread
            .warnings
            .iter()
            .any(|w| w.contains("record association"))
    );
}
