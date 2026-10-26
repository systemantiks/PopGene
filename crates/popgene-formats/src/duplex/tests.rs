use super::*;
use popgene_core::{
    Feature, FeatureLocation, Segment,
    cloning::{self, AssemblyInput, LigationOptions},
    metadata::Reference,
    restriction::{self, DigestConditions, ENZYMES},
    translation, workflow,
};

fn source(phase: usize) -> Document {
    let mut doc = Document::new("source", "ATGGAATTCAAGTAA").unwrap();
    doc.annotate("coding α".into(), Segment { start: 0, end: 15 }, false)
        .unwrap();
    doc.features[0].kind = "CDS".into();
    doc.features[0]
        .qualifiers
        .push(("codon_start".into(), Some(phase.to_string())));
    doc
}
fn digest(doc: &Document) -> Vec<Document> {
    let enzyme = ENZYMES.iter().position(|e| e.name == "EcoRI").unwrap();
    let conditions = DigestConditions::default();
    let result = restriction::digest_document(doc, &[enzyme], conditions).unwrap();
    workflow::digest_products(doc, &result.fragments, &["EcoRI".into()], conditions).unwrap()
}
fn ligate(docs: &[Document], reverse: bool, circular: bool) -> Document {
    cloning::ligate(
        &docs
            .iter()
            .map(|d| AssemblyInput {
                document: d.clone(),
                reverse,
            })
            .collect::<Vec<_>>(),
        LigationOptions {
            circular,
            phosphorylate: true,
            allow_nicks: false,
        },
        "rejoined",
    )
    .unwrap()
}
fn reopen(doc: &Document, format: SaveFormat) -> Document {
    let (preserved, warnings) = preservation(doc, format);
    assert!(preserved, "{}: {warnings:?}", format.label());
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join(format!("saved.{}", format.extension()));
    crate::save(&path, doc, format, None).unwrap();
    let bytes = std::fs::read(&path).unwrap();
    let result = crate::open(&path).unwrap();
    assert_eq!(crate::encode(&result, format).unwrap(), bytes);
    assert_eq!(result.sequence, doc.sequence);
    assert_eq!(result.duplex_ends, doc.duplex_ends, "{:?}", result.warnings);
    assert!(result.comments.iter().all(|c| !c.contains(START)));
    result
}
fn aa(doc: &Document, f: &Feature) -> String {
    String::from_utf8(translation::translate_cds(doc, f).unwrap().amino_acids).unwrap()
}
fn emit_qa(doc: &Document, format: SaveFormat, label: &str) {
    if let Some(directory) = std::env::var_os("POPGENE_PHYSICAL_QA_DIR") {
        let directory = std::path::PathBuf::from(directory);
        std::fs::create_dir_all(&directory).unwrap();
        let path = directory.join(format!("{label}.{}", format.extension()));
        assert!(!path.exists(), "QA output must be fresh");
        crate::save(&path, doc, format, None).unwrap();
    }
}

#[test]
fn real_save_reopen_of_cohesive_partial_cds_rejoins_all_three_phases() {
    for format in [SaveFormat::GenBank, SaveFormat::Embl] {
        for (phase, expected) in [(1, "MEFK*"), (2, "WNSS"), (3, "GIQV")] {
            let doc = source(phase);
            let fragments = digest(&doc);
            for (index, fragment) in fragments.iter().enumerate() {
                emit_qa(
                    fragment,
                    format,
                    &format!("ecori-phase-{phase}-fragment-{}", index + 1),
                );
            }
            assert_eq!(fragments[0].sequence.as_ref(), "ATGG");
            assert_eq!(
                fragments[0].duplex_ends.as_ref().unwrap().opposite_sequence,
                "AATTCCAT"
            );
            assert_eq!(fragments[1].sequence.as_ref(), "AATTCAAGTAA");
            assert_eq!(
                fragments[1].duplex_ends.as_ref().unwrap().opposite_sequence,
                "TTACTTG"
            );
            let saved = fragments
                .iter()
                .map(|f| reopen(f, format))
                .collect::<Vec<_>>();
            for (original, reread) in fragments.iter().zip(&saved) {
                assert_eq!(
                    reread.annotation_lineage.origins,
                    original.annotation_lineage.origins
                );
                assert_eq!(
                    reread.annotation_lineage.pieces.len(),
                    original.annotation_lineage.pieces.len(),
                    "{:?}",
                    reread.warnings
                );
                assert!(
                    reread
                        .annotation_lineage
                        .pieces
                        .iter()
                        .all(|p| reread.annotation_lineage.actionable_piece(reread, p))
                );
            }
            let product = ligate(&saved, false, false);
            assert_eq!(product.sequence, doc.sequence);
            assert_eq!(product.features.len(), 1);
            assert_eq!(product.features[0].location_label(), "1..15");
            assert_eq!(aa(&product, &product.features[0]), expected);
            assert!(!product.features[0].edited);
            emit_qa(&product, format, &format!("ecori-phase-{phase}-rejoined"));
        }
    }
}

#[test]
fn colored_cited_fragments_keep_physical_state_and_verified_appended_display_reference() {
    let mut doc = source(2);
    doc.features[0].qualifiers.extend([
        ("snapgene_color".into(), Some("#336699".into())),
        ("citation".into(), Some("1".into())),
    ]);
    doc.references.push(Reference {
        title: "Physical fragment integration α".into(),
        authors: Some("Synthetic, A.".into()),
        ..Default::default()
    });
    doc.sequence_qualities = Some((0..15).map(|q| (q != 1).then_some(q)).collect());
    doc.opposite_sequence_qualities = Some((30..45).map(|q| (q != 32).then_some(q)).collect());
    doc.sequence_colors.push(SequenceColor {
        range: Segment { start: 9, end: 12 },
        reverse: true,
        color_rgb: 0x993366,
    });
    let fragments = digest(&doc);
    for format in [SaveFormat::GenBank, SaveFormat::Embl] {
        let reopened: Vec<_> = fragments.iter().map(|f| reopen(f, format)).collect();
        for (index, (original, saved)) in fragments.iter().zip(&reopened).enumerate() {
            assert_eq!(saved.references[0], doc.references[0]);
            assert_eq!(saved.sequence_qualities, original.sequence_qualities);
            assert_eq!(
                saved.opposite_sequence_qualities,
                original.opposite_sequence_qualities
            );
            assert_eq!(saved.sequence_colors, original.sequence_colors);
            assert_eq!(
                saved.duplex_annotations.sequence_colors,
                original.duplex_annotations.sequence_colors
            );
            assert_eq!(
                saved.annotation_lineage.origins,
                original.annotation_lineage.origins
            );
            assert_eq!(saved.annotation_lineage.pieces.len(), 1);
            assert!(
                saved
                    .annotation_lineage
                    .actionable_piece(saved, &saved.annotation_lineage.pieces[0])
            );
            if format == SaveFormat::GenBank && index == 1 {
                assert_eq!(saved.references.len(), 2);
                assert!(crate::genbank_display::is_generated_reference(
                    &saved.references[1],
                    1
                ));
            } else {
                assert_eq!(saved.references.len(), 1);
            }
        }
        let rejoined = ligate(&reopened, false, false);
        assert_eq!(rejoined.features.len(), 1);
        assert_eq!(aa(&rejoined, &rejoined.features[0]), "WNSS");
        assert_eq!(rejoined.features[0].display_color_rgb(), Some(0x336699));
        let slot = rejoined.features[0]
            .qualifier("citation")
            .unwrap()
            .trim_matches(['[', ']'])
            .parse::<usize>()
            .unwrap();
        assert_eq!(rejoined.references[slot - 1], doc.references[0]);
    }
    // An arbitrary extra record, a changed marker, or an invalid original
    // citation cannot borrow the writer's precisely recognized exception.
    let original = &fragments[1];
    let canonical = ordinary(&plain(original, SaveFormat::GenBank), SaveFormat::GenBank);
    assert!(references_for_codec(original, &canonical, SaveFormat::GenBank).is_ok());
    for field in ["title", "journal", "remark", "scope", "extra"] {
        let mut forged = canonical.clone();
        let reference = forged.references.last_mut().unwrap();
        match field {
            "title" => reference.title = "Unrelated reference".into(),
            "journal" => reference.journal = Some("Other SnapGene-like text".into()),
            "remark" => reference.remark = Some("Different meaning".into()),
            "scope" => reference.scope = Some("99".into()),
            "extra" => reference.extra.push(("unknown".into(), "record".into())),
            _ => unreachable!(),
        }
        assert!(
            references_for_codec(original, &forged, SaveFormat::GenBank).is_err(),
            "{field}"
        );
    }
    let mut invalid = original.clone();
    invalid.features[0]
        .qualifiers
        .iter_mut()
        .find(|(key, _)| key == "citation")
        .unwrap()
        .1 = Some("2".into());
    assert!(references_for_codec(&invalid, &canonical, SaveFormat::GenBank).is_err());
}

#[test]
fn reverse_and_circular_reassembly_survive_public_fragment_files() {
    for format in [SaveFormat::GenBank, SaveFormat::Embl] {
        let doc = source(1);
        let mut parts = digest(&doc)
            .iter()
            .map(|f| reopen(f, format))
            .collect::<Vec<_>>();
        parts.reverse();
        let reversed = ligate(&parts, true, false);
        assert_eq!(reversed.sequence.as_ref(), "TTACTTGAATTCCAT");
        assert_eq!(reversed.features.len(), 1);
        assert_eq!(reversed.features[0].location_label(), "complement(1..15)");
        assert_eq!(aa(&reversed, &reversed.features[0]), "MEFK*");
        let mut circle = doc;
        circle.circular = true;
        let parts = digest(&circle)
            .iter()
            .map(|f| reopen(f, format))
            .collect::<Vec<_>>();
        let closed = ligate(&parts, false, true);
        assert_eq!(closed.sequence.as_ref(), "AATTCAAGTAAATGG");
        assert_eq!(closed.features.len(), 1);
        assert_eq!(closed.features[0].location_label(), "join(12..15,1..11)");
        assert_eq!(aa(&closed, &closed.features[0]), "MEFK*");
    }
}

#[test]
fn current_physical_scores_chemistry_colors_and_duplicate_occurrences_are_exact() {
    let mut doc = source(1);
    doc.sequence_qualities = Some((0..15).map(|i| (i != 1).then_some(i)).collect());
    doc.opposite_sequence_qualities = Some((30..45).map(|i| (i != 32).then_some(i)).collect());
    doc.methylation = Some(DnaMethylation {
        dam: true,
        dcm: false,
        eco_ki: true,
    });
    doc.references.push(Reference {
        title: "Current source α".into(),
        ..Reference::default()
    });
    doc.features[0]
        .qualifiers
        .push(("citation".into(), Some("1".into())));
    doc.sequence_colors.push(SequenceColor {
        range: Segment { start: 4, end: 8 },
        reverse: true,
        color_rgb: 0xff0000,
    });
    let mut left = digest(&doc).remove(0);
    left.duplex_ends.as_mut().unwrap().left_modification = Some("5′ hydroxyl".into());
    let mut second = left.duplex_annotations.features[0].clone();
    second.id = 100;
    left.duplex_annotations.features.push(second);
    // Current supplemental metadata is useful independently of ancestry. This
    // deliberately added feature is not assigned an invented ancestor.
    left.annotation_lineage = Default::default();
    for format in [SaveFormat::GenBank, SaveFormat::Embl] {
        let saved = reopen(&left, format);
        assert_eq!(saved.sequence_qualities, left.sequence_qualities);
        assert_eq!(
            saved.opposite_sequence_qualities,
            left.opposite_sequence_qualities
        );
        assert_eq!(saved.methylation, left.methylation);
        assert_eq!(
            saved.duplex_annotations.sequence_colors,
            left.duplex_annotations.sequence_colors
        );
        assert_eq!(saved.duplex_annotations.features.len(), 2);
        assert_eq!(saved.duplex_annotations.features[0].id, 1);
        assert_eq!(saved.duplex_annotations.features[1].id, 2);
        assert_eq!(
            saved.duplex_annotations.features[0].label,
            saved.duplex_annotations.features[1].label
        );
        assert_eq!(saved.references[0].title, "Current source α");
        assert_eq!(saved.sequence_qualities.as_ref().unwrap()[0], Some(0));
        assert_eq!(saved.sequence_qualities.as_ref().unwrap()[1], None);
        let warnings = crate::save_warnings(&left, format).join(" ");
        assert!(
            warnings.contains("actual top strand") && !warnings.contains("qualities are omitted"),
            "{warnings}"
        );
    }
}

fn ordinary(bytes: &[u8], format: SaveFormat) -> Document {
    if format == SaveFormat::GenBank {
        crate::genbank(bytes).unwrap()
    } else {
        crate::embl::read(std::str::from_utf8(bytes).unwrap()).unwrap()
    }
}
fn plain(doc: &Document, format: SaveFormat) -> Vec<u8> {
    if format == SaveFormat::GenBank {
        crate::genbank_bytes_plain(doc).unwrap()
    } else {
        crate::embl::write_plain(doc).unwrap()
    }
}
fn altered_payload(bytes: &[u8], format: SaveFormat, change: impl FnOnce(&mut Value)) -> Vec<u8> {
    let mut doc = ordinary(bytes, format);
    let parsed = parse(&doc.comments).unwrap().unwrap();
    let mut value = parsed.value;
    change(&mut value);
    doc.comments[parsed.comment].replace_range(parsed.start..parsed.end, &block(&value).unwrap());
    let changed = plain(&doc, format);
    assert_eq!(
        crate::lineage::records_hash(bytes, format).unwrap(),
        crate::lineage::records_hash(&changed, format).unwrap()
    );
    changed
}

#[test]
fn physical_context_is_canonical_across_map_features_and_reads_exact_legacy() {
    let doc = Document::new("context", "ACGT").unwrap();
    let hashes = context(&doc).unwrap();
    assert_eq!(
        hashes.canonical,
        "183f836fc2cdb2d8beaee38ef101aa3e703d21d9dbd500c0c1379f6426127979"
    );
    assert_eq!(
        hashes.legacy_insertion(),
        "abaf39b9a5858d39dd40e5179a9d066741acbb56a0e3cfaad2ec1a4ffce0e781"
    );
    let fragment = digest(&source(2)).remove(0);
    for format in [SaveFormat::GenBank, SaveFormat::Embl] {
        let bytes = crate::encode(&fragment, format).unwrap();
        let ordinary = ordinary(&bytes, format);
        let hashes = context(&ordinary).unwrap();
        let parsed = parse(&ordinary.comments).unwrap().unwrap();
        assert_eq!(parsed.value["context_sha256"], hashes.canonical);
        let legacy = altered_payload(&bytes, format, |value| {
            value["context_sha256"] = Value::from(hashes.legacy_insertion());
        });
        let reread = crate::import(&format!("legacy.{}", format.extension()), &legacy).unwrap();
        assert_eq!(reread.duplex_ends, fragment.duplex_ends);
        assert_eq!(reread.duplex_annotations, fragment.duplex_annotations);
        assert!(!reread.annotation_lineage.is_empty());
        assert_eq!(crate::encode(&reread, format).unwrap(), legacy);
        let corrupt = altered_payload(&bytes, format, |value| {
            value["context_sha256"] = Value::from("0".repeat(64));
        });
        inactive(&corrupt, format);
    }
}
fn inactive(bytes: &[u8], format: SaveFormat) {
    let doc = crate::import(&format!("edited.{}", format.extension()), bytes).unwrap();
    assert!(doc.duplex_ends.is_none());
    assert!(doc.duplex_annotations.is_empty());
    assert!(doc.annotation_lineage.pieces.is_empty());
    assert!(doc.comments.iter().any(|c| c.contains(START)));
    assert!(
        doc.warnings
            .iter()
            .any(|w| w.contains("no physical metadata was activated"))
    );
}

#[test]
fn forged_structures_fail_atomic_validation_even_with_a_recomputed_checksum() {
    let fragment = digest(&source(1)).remove(0);
    for format in [SaveFormat::GenBank, SaveFormat::Embl] {
        let bytes = crate::encode(&fragment, format).unwrap();
        let changes: &[fn(&mut Value)] = &[
            |v| v["fragment"]["ends"]["opposite_sequence"] = Value::from("TATTCCAT"),
            |v| v["fragment"]["ends"]["right"]["kind"] = Value::from("ThreePrime"),
            |v| v["fragment"]["ends"]["extra"] = Value::from(true),
            |v| v["fragment"]["duplex_annotations"]["features"][0]["id"] = Value::from(99),
            |v| v["fragment"]["duplex_annotations"]["features"][0]["extra"] = Value::from(true),
            |v| {
                v["fragment"]["duplex_annotations"]["features"][0]["parts"][0]["end"] =
                    Value::from(99)
            },
            |v| {
                v["fragment"]["duplex_annotations"]["features"][0]["qualifiers"]
                    .as_array_mut()
                    .unwrap()
                    .push(serde_json::json!(["citation", "99"]))
            },
            |v| v["fragment"]["sequence_qualities"] = serde_json::json!([0]),
            |v| {
                v["fragment"]["opposite_sequence_qualities"] =
                    serde_json::json!([99, 0, 0, 0, 0, 0, 0, 0])
            },
        ];
        for change in changes {
            inactive(&altered_payload(&bytes, format, change), format);
        }
    }
}

#[test]
fn stale_top_records_metadata_unknown_versions_and_collisions_stay_readable_inactive() {
    let mut fragment = digest(&source(1)).remove(1);
    fragment.references.push(Reference {
        title: "Bound reference".into(),
        ..Reference::default()
    });
    // An ordinary top annotation ensures raw record changes are independently
    // bound even when the changed annotation is unrelated to the extent CDS.
    fragment.features.push(Feature {
        id: 99,
        label: "top mark".into(),
        kind: "misc_feature".into(),
        parts: vec![Segment { start: 4, end: 5 }],
        reverse: false,
        qualifiers: vec![],
        edited: false,
        location: None,
    });
    for format in [SaveFormat::GenBank, SaveFormat::Embl] {
        let bytes = crate::encode(&fragment, format).unwrap();
        let edits: &[fn(&mut Document)] = &[
            |d| d.sequence = "TATTCAAGTAA".into(),
            |d| d.circular = true,
            |d| d.single_stranded = true,
            |d| d.description.push_str(" external edit"),
            |d| d.references[0].title = "Changed reference".into(),
            |d| d.features[0].label = "Changed annotation".into(),
            |d| {
                d.features.reverse();
            },
            |d| {
                let c = d.comments.iter_mut().find(|c| c.contains(START)).unwrap();
                *c = c.replace("physical fragment v1", "physical fragment v99");
            },
            |d| {
                let c = d
                    .comments
                    .iter()
                    .find(|c| c.contains(START))
                    .unwrap()
                    .clone();
                d.comments.push(c);
            },
        ];
        for edit in edits {
            let mut changed = ordinary(&bytes, format);
            edit(&mut changed);
            inactive(&plain(&changed, format), format);
        }
        let mut stripped = ordinary(&bytes, format);
        stripped.comments.retain(|c| !c.contains(START));
        let stripped = crate::import(
            &format!("stripped.{}", format.extension()),
            &plain(&stripped, format),
        )
        .unwrap();
        assert!(stripped.duplex_ends.is_none());
        assert!(stripped.annotation_lineage.pieces.is_empty());
    }
}

#[test]
fn optional_budget_fallback_keeps_ordinary_output_and_unextended_quality_limits() {
    let mut fragment = digest(&source(1)).remove(0);
    fragment.duplex_annotations.features[0].label = "x".repeat(MAX_PAYLOAD);
    // The annotation change makes its old ancestry inactive independently. A
    // fresh ordinary fragment with no lineage also follows the same fallback.
    fragment.annotation_lineage = Default::default();
    for format in [SaveFormat::GenBank, SaveFormat::Embl] {
        let bytes = crate::encode(&fragment, format).unwrap();
        assert_eq!(bytes, plain(&fragment, format));
        let saved = crate::import(&format!("large.{}", format.extension()), &bytes).unwrap();
        assert_eq!(saved.sequence, fragment.sequence);
        assert!(saved.duplex_ends.is_none());
        assert!(saved.duplex_annotations.is_empty());
        assert!(
            crate::save_warnings(&fragment, format)
                .join(" ")
                .contains("4 MiB")
        );
    }
    let ordinary = source(1);
    assert!(
        !crate::encode(&ordinary, SaveFormat::GenBank)
            .unwrap()
            .windows(START.len())
            .any(|w| w == START.as_bytes())
    );
    let mut paired_quality = ordinary;
    paired_quality.opposite_sequence_qualities = Some(vec![Some(0); 15]);
    assert!(
        crate::save_warnings(&paired_quality, SaveFormat::GenBank)
            .join(" ")
            .contains("Opposite-strand quality scores are omitted")
    );
}

#[test]
fn a_large_quality_vector_drops_only_optional_physical_and_active_ancestry_data() {
    let mut doc = source(1);
    doc.sequence = format!("{}{}", doc.sequence, "A".repeat(850_000)).into();
    doc.sequence_qualities = Some(vec![None; doc.sequence.len()]);
    let fragment = digest(&doc).remove(1);
    assert!(
        fragment
            .annotation_lineage
            .pieces
            .iter()
            .any(|p| fragment.annotation_lineage.actionable_piece(&fragment, p))
    );
    assert!(
        preflight(&fragment)
            .unwrap_err()
            .to_string()
            .contains("4 MiB")
    );
    for format in [SaveFormat::GenBank, SaveFormat::Embl] {
        let bytes = crate::encode(&fragment, format).unwrap();
        let saved = crate::import(&format!("large.{}", format.extension()), &bytes).unwrap();
        assert_eq!(saved.sequence, fragment.sequence);
        assert!(saved.duplex_ends.is_none());
        assert!(saved.annotation_lineage.pieces.is_empty());
        assert!(saved.annotation_lineage.tracking_limited);
        let warnings = crate::save_warnings(&fragment, format).join(" ");
        assert!(
            warnings.contains("4 MiB")
                && warnings.contains("Active annotation ancestry is omitted"),
            "{warnings}"
        );
    }
}

#[test]
fn a_stale_valid_carrier_cannot_override_current_edited_end_chemistry() {
    let doc = digest(&source(1)).remove(0);
    for format in [SaveFormat::GenBank, SaveFormat::Embl] {
        let bytes = crate::encode(&doc, format).unwrap();
        let parsed = ordinary(&bytes, format);
        let old = parsed
            .comments
            .iter()
            .find(|c| c.contains(START))
            .unwrap()
            .clone();
        let mut edited = doc.clone();
        edited.duplex_ends.as_mut().unwrap().left_modification = Some("5′ hydroxyl".into());
        edited.comments.push(old);
        let error = crate::encode(&edited, format).unwrap_err().to_string();
        assert!(
            error.contains("existing supplement text cannot be replaced safely"),
            "{error}"
        );
        assert!(
            crate::save_warnings(&edited, format)
                .join(" ")
                .contains("cannot be written")
        );
        for comment in &mut edited.comments {
            *comment = comment
                .replace("physical fragment v1", "physical fragment\nv1")
                .replace("physical fragment end]", "physical fragment\nend]");
        }
        assert!(
            crate::encode(&edited, format)
                .unwrap_err()
                .to_string()
                .contains("existing supplement text cannot be replaced safely")
        );
    }
}

#[test]
fn canonical_serialization_cannot_reactivate_previously_inactive_carriers() {
    for physical in [false, true] {
        let doc = if physical {
            digest(&source(1)).remove(1)
        } else {
            workflow::extract_region(&source(1), 0, 15, "whole").unwrap()
        };
        for format in [SaveFormat::GenBank, SaveFormat::Embl] {
            let text = String::from_utf8(crate::encode(&doc, format).unwrap()).unwrap();
            let phase = if physical { 3 } else { 1 };
            // This changes raw binding bytes while retaining precisely the
            // parsed biological annotation. The ordinary writer removes quotes
            // around codon_start and would otherwise make the old hash valid.
            let stale = text.replace(
                &format!("/codon_start={phase}"),
                &format!("/codon_start=\"{phase}\""),
            );
            assert_ne!(stale, text);
            let mut inactive =
                crate::import(&format!("stale.{}", format.extension()), stale.as_bytes()).unwrap();
            assert!(inactive.duplex_ends.is_none());
            assert!(inactive.annotation_lineage.is_empty());
            inactive
                .comments
                .push("Ordinary unrelated comment change".into());
            let error = crate::encode(&inactive, format).unwrap_err().to_string();
            assert!(
                error.contains("previously inactive") && error.contains("would become active"),
                "{error}"
            );
            // Unsupported versions remain readable, inactive text rather than
            // making ordinary saves needlessly impossible.
            for comment in &mut inactive.comments {
                *comment = comment
                    .replace("physical fragment v1", "physical fragment v99")
                    .replace("lineage v1", "lineage v99");
            }
            let bytes = crate::encode(&inactive, format).unwrap();
            let reread = crate::import(&format!("unknown.{}", format.extension()), &bytes).unwrap();
            assert!(reread.duplex_ends.is_none());
            assert!(reread.annotation_lineage.is_empty());
            assert!(reread.comments.iter().any(|c| c.contains("v99 begin]")));

            let mut split_prefix = ordinary(text.as_bytes(), format);
            for comment in &mut split_prefix.comments {
                *comment = comment
                    .replace("physical fragment v1", "physical fragment\nv1")
                    .replace("physical fragment end]", "physical fragment\nend]")
                    .replace("lineage v1", "lineage\nv1")
                    .replace("lineage end]", "lineage\nend]");
            }
            split_prefix
                .comments
                .push("Unrelated [PopGene example prose".into());
            if format == SaveFormat::Embl {
                // EMBL folds the embedded newlines into marker-separating spaces.
                let error = crate::encode(&split_prefix, format)
                    .unwrap_err()
                    .to_string();
                assert!(
                    error.contains("previously inactive") && error.contains("would become active"),
                    "{error}"
                );
            } else {
                // GenBank retains these embedded newlines, so the text stays
                // inactive and must remain savable without unnecessary refusal.
                let bytes = crate::encode(&split_prefix, format).unwrap();
                let reread = crate::import("split.gb", &bytes).unwrap();
                assert!(reread.duplex_ends.is_none() && reread.annotation_lineage.is_empty());
                assert!(reread.comments.iter().any(|c| c.contains("lineage\nv1")));
            }

            let mut unmatched = source(1);
            unmatched.comments = vec![END.into(), "[PopGene lineage end]".into()];
            let bytes = crate::encode(&unmatched, format).unwrap();
            let reread =
                crate::import(&format!("unmatched.{}", format.extension()), &bytes).unwrap();
            assert!(reread.duplex_ends.is_none() && reread.annotation_lineage.is_empty());
            assert!(reread.comments.iter().any(|c| c.contains(END)));
        }
    }
}

#[test]
fn all_four_terminal_polarities_and_repeated_join_occurrences_survive_reassembly() {
    use popgene_core::restriction::{EndKind, FragmentEnd};
    let axis = "TTAACCGGAATTCTGGAACCGGTT";
    for left in [EndKind::FivePrime, EndKind::ThreePrime] {
        for right in [EndKind::FivePrime, EndKind::ThreePrime] {
            let ts = usize::from(left == EndKind::ThreePrime) * 2;
            let te = axis.len() - usize::from(right == EndKind::FivePrime) * 3;
            let bs = usize::from(left == EndKind::FivePrime) * 2;
            let be = axis.len() - usize::from(right == EndKind::ThreePrime) * 3;
            let mut doc = Document::new("asymmetric", &axis[ts..te]).unwrap();
            doc.duplex_ends = Some(DuplexEnds {
                left: FragmentEnd {
                    kind: left,
                    overhang: if left == EndKind::FivePrime {
                        axis[..2].into()
                    } else {
                        popgene_core::reverse_complement(&axis[..2])
                    },
                    enzymes: vec![],
                },
                right: FragmentEnd {
                    kind: right,
                    overhang: if right == EndKind::ThreePrime {
                        axis[axis.len() - 3..].into()
                    } else {
                        popgene_core::reverse_complement(&axis[axis.len() - 3..])
                    },
                    enzymes: vec![],
                },
                opposite_sequence: popgene_core::reverse_complement(&axis[bs..be]),
                left_modification: Some("5′ phosphate".into()),
                right_modification: Some("5′ phosphate".into()),
            });
            doc.duplex_annotations.features.push(Feature {
                id: 8,
                label: "whole reverse".into(),
                kind: "misc_feature".into(),
                parts: vec![Segment {
                    start: 0,
                    end: axis.len(),
                }],
                reverse: true,
                location: None,
                qualifiers: vec![],
                edited: false,
            });
            for format in [SaveFormat::GenBank, SaveFormat::Embl] {
                let parts = digest(&doc)
                    .iter()
                    .map(|f| reopen(f, format))
                    .collect::<Vec<_>>();
                let product = ligate(&parts, false, false);
                assert_eq!(product.sequence, doc.sequence);
                assert_eq!(product.duplex_ends, doc.duplex_ends);
                let reference = workflow::duplex_reference(&product).unwrap();
                assert_eq!(reference.features.len(), 1);
                assert_eq!(
                    reference.feature_sequence(&reference.features[0]).unwrap(),
                    popgene_core::reverse_complement(axis)
                );
            }
        }
    }
    let mut repeated = source(1);
    let range = |start, end| FeatureLocation::Range {
        start,
        end,
        before: false,
        after: false,
    };
    let location = FeatureLocation::Join(vec![range(0, 8), range(0, 8), range(4, 12)]);
    repeated.features[0].parts = location.local_parts();
    repeated.features[0].location = Some(location.clone());
    repeated.features[0].kind = "misc_feature".into();
    for format in [SaveFormat::GenBank, SaveFormat::Embl] {
        let parts = digest(&repeated)
            .iter()
            .map(|f| reopen(f, format))
            .collect::<Vec<_>>();
        let result = ligate(&parts, false, false);
        assert_eq!(result.features.len(), 1);
        assert_eq!(result.features[0].location, Some(location.clone()));
        assert_eq!(
            result.feature_sequence(&result.features[0]).unwrap(),
            "ATGGAATTATGGAATTAATTCAAG"
        );
    }
}

#[test]
fn native_fuzzy_refusal_is_unchanged_and_genuine_partial_endpoints_are_not_completed() {
    let mut doc = source(1);
    doc.features[0].location = Some(FeatureLocation::Range {
        start: 0,
        end: 15,
        before: true,
        after: true,
    });
    let fragments = digest(&doc);
    for format in [SaveFormat::GenBank, SaveFormat::Embl] {
        let saved = fragments
            .iter()
            .map(|f| reopen(f, format))
            .collect::<Vec<_>>();
        let product = ligate(&saved, false, false);
        assert_eq!(product.features.len(), 1);
        assert_eq!(product.features[0].location_label(), "<1..>15");
        assert!(
            crate::encode(&product, SaveFormat::SnapGene)
                .unwrap_err()
                .to_string()
                .contains("fuzzy")
        );
    }
}
