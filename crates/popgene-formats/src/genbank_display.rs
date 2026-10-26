//! The published extended GenBank convention requires an advertised header as
//! well as terminal notes. Never replace biological references with that marker.
//! https://support.snapgene.com/hc/en-us/articles/10242682237588
use gb_io::seq::{Reference, Seq};
use popgene_core::Document;

fn export_reference(index: usize) -> Reference {
    Reference {
        description: (index + 1).to_string(),
        authors: Some(".".into()),
        consortium: None,
        title: "Direct Submission".into(),
        journal: Some("PopGene export using SnapGene GenBank format".into()),
        pubmed: None,
        remark: Some("Software interoperability marker for display annotations; not a scientific citation or a database submission.".into()),
    }
}

/// Recognize only our complete generated terminal record, not the broader
/// third-party display advertisement accepted by `marker`.
pub(super) fn is_generated_reference(
    actual: &popgene_core::metadata::Reference,
    index: usize,
) -> bool {
    let expected = super::references::from_genbank(&export_reference(index));
    let words = |value: Option<&str>| {
        value.map(|text| text.split_whitespace().collect::<Vec<_>>().join(" "))
    };
    actual.title == expected.title
        && actual.authors == expected.authors
        && actual.consortium == expected.consortium
        && actual.pubmed == expected.pubmed
        && actual.scope == expected.scope
        && actual.extra.is_empty()
        && words(actual.journal.as_deref()) == words(expected.journal.as_deref())
        && words(actual.remark.as_deref()) == words(expected.remark.as_deref())
}

pub(super) fn advertised(record: &Seq) -> bool {
    record
        .name
        .as_deref()
        .is_some_and(|name| name.contains("Exported"))
        && marker(record.references.last())
}
fn marker(reference: Option<&Reference>) -> bool {
    reference.is_some_and(|reference| {
        reference.title == "Direct Submission"
            && reference
                .journal
                .as_deref()
                .is_some_and(|journal| journal.contains("SnapGene"))
    })
}
pub(super) fn uses_extension(doc: &Document) -> bool {
    doc.features
        .iter()
        .any(|feature| feature.interchange_display_note().is_some())
        || doc.source.as_ref().is_some_and(|source| {
            source.format == "GenBank"
                && super::one_genbank(&source.bytes).is_ok_and(|record| advertised(&record))
        })
}
pub(super) fn apply_header(doc: &Document, record: &mut Seq, source_advertised: bool) {
    if !source_advertised
        && let Some(keywords) = record
            .keywords
            .as_deref()
            .filter(|value| !value.is_empty() && *value != "." && *value != doc.name)
    {
        let comment = format!("Original GenBank KEYWORDS: {keywords}");
        if !record.comments.contains(&comment) {
            record.comments.push(comment);
        }
    }
    record.name = Some("Exported".into());
    record.keywords = Some(doc.name.split_whitespace().collect::<Vec<_>>().join(" "));
    if !marker(record.references.last()) {
        record
            .references
            .push(export_reference(record.references.len()));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popgene_core::{Feature, Segment, metadata::Reference as Citation};
    fn normalized(value: &str) -> String {
        value.split_whitespace().collect::<Vec<_>>().join(" ")
    }
    #[test]
    fn generated_display_notes_keep_default_arrows_and_explicit_nondirectionality() {
        use popgene_core::FeatureDirection;
        for reverse in [false, true] {
            for named in [false, true] {
                for nondirectional in [false, true] {
                    let mut doc = Document::new("default direction", "AACCGGTTAACC").unwrap();
                    let mut qualifiers = vec![("snapgene_color".into(), Some("#ff0000".into()))];
                    if named {
                        qualifiers.push(("segment_name".into(), Some("domain".into())));
                    }
                    if nondirectional {
                        qualifiers.push(("snapgene_directionality".into(), Some("0".into())));
                    }
                    doc.features.push(Feature {
                        id: 1,
                        label: "colored domain".into(),
                        kind: "misc_feature".into(),
                        parts: vec![Segment { start: 2, end: 8 }],
                        location: None,
                        reverse,
                        qualifiers,
                        edited: false,
                    });
                    let expected = if nondirectional {
                        FeatureDirection::Nondirectional
                    } else if reverse {
                        FeatureDirection::Reverse
                    } else {
                        FeatureDirection::Forward
                    };
                    for format in [crate::SaveFormat::GenBank, crate::SaveFormat::Embl] {
                        let bytes = crate::encode(&doc, format).unwrap();
                        let reopened = crate::import("direction", &bytes).unwrap();
                        let feature = &reopened.features[0];
                        assert_eq!(feature.display_direction(), Some(expected));
                        assert_eq!(feature.display_color_rgb(), Some(0xff0000));
                        assert_eq!(feature.parts, doc.features[0].parts);
                        assert_eq!(feature.reverse, reverse);
                        assert_eq!(
                            reopened.feature_sequence(feature).unwrap(),
                            doc.feature_sequence(&doc.features[0]).unwrap()
                        );
                        if named {
                            assert_eq!(
                                feature.display_segments()[0].name.as_deref(),
                                Some("domain")
                            );
                        }
                    }
                }
            }
        }
    }
    #[test]
    fn native_style_export_advertises_grammar_without_changing_biology_or_citations() {
        let mut doc = Document::new("Custom map label", "AACCGGTTAACC").unwrap();
        doc.references.push(Citation {
            title: "Actual scientific citation".into(),
            authors: Some("Example,A.".into()),
            journal: Some("Example journal".into()),
            ..Default::default()
        });
        doc.features.push(Feature {
            id: 1,
            label: "Directional color".into(),
            kind: "misc_feature".into(),
            parts: vec![Segment { start: 2, end: 8 }],
            reverse: true,
            location: None,
            qualifiers: vec![
                ("snapgene_color".into(), Some("#aabbcc".into())),
                ("snapgene_directionality".into(), Some("3".into())),
                ("note".into(), Some("Biological note.".into())),
            ],
            edited: true,
        });
        let original = doc.clone();
        let bytes = crate::genbank_bytes(&doc).unwrap();
        let record = super::super::one_genbank(&bytes).unwrap();
        assert!(advertised(&record));
        assert_eq!(record.keywords.as_deref(), Some("Custom map label"));
        assert_eq!(record.references.len(), 2);
        assert_eq!(record.references[0].title, "Actual scientific citation");
        assert_eq!(record.references[0].authors.as_deref(), Some("Example,A."));
        assert!(
            normalized(record.references[1].remark.as_deref().unwrap())
                .contains("not a scientific citation")
        );
        assert_eq!(
            record.features[0].qualifier_values("note").last(),
            Some("color: #aabbcc; direction: BOTH")
        );
        let reopened = crate::import("styles.gb", &bytes).unwrap();
        assert_eq!(reopened.name, doc.name);
        assert_eq!(reopened.features[0].parts, doc.features[0].parts);
        assert_eq!(reopened.features[0].reverse, doc.features[0].reverse);
        assert_eq!(reopened.features[0].display_color_rgb(), Some(0xaabbcc));
        assert_eq!(
            reopened.features[0].display_direction(),
            Some(popgene_core::FeatureDirection::Bidirectional)
        );
        let again = crate::genbank_bytes(&reopened).unwrap();
        assert_eq!(
            super::super::one_genbank(&again).unwrap().references.len(),
            2
        );
        assert_eq!(doc, original);
    }
    #[test]
    fn primer_oligo_tail_case_and_review_note_keep_terminal_metadata() {
        let mut doc = Document::new("primer case", "AACCGGTTAACC").unwrap();
        doc.features.push(Feature {
            id: 1,
            label: "oligo".into(),
            kind: "primer_bind".into(),
            parts: vec![Segment { start: 0, end: 6 }],
            reverse: false,
            location: None,
            qualifiers: vec![
                ("primer_sequence".into(), Some("ttAACCGG".into())),
                (
                    "note".into(),
                    Some("color: blue; added: 2026-10-20; 5' phosphorylated".into()),
                ),
            ],
            edited: true,
        });
        let bytes = crate::genbank_bytes(&doc).unwrap();
        let record = super::super::one_genbank(&bytes).unwrap();
        assert_eq!(
            record.features[0]
                .qualifier_values("note")
                .last()
                .map(normalized)
                .as_deref(),
            Some("color: blue; sequence: ttAACCGG; added: 2026-10-20; 5' phosphorylated")
        );
        let restored = crate::import("primer.gb", &bytes).unwrap();
        assert_eq!(
            restored.features[0].primer_sequence().as_deref(),
            Some("ttAACCGG")
        );
        assert_eq!(restored.features[0].primer_added_date(), Some("2026-10-20"));
    }
    #[test]
    fn ordinary_keywords_remain_provenance_not_an_unadvertised_map_label() {
        let input=b"LOCUS       original 12 bp DNA linear UNK 01-JAN-2000\nKEYWORDS    biological; classification.\nORIGIN\n        1 aaccggttaacc\n//\n";
        let doc = crate::import("plain.gb", input).unwrap();
        assert_eq!(doc.name, "original");
        let bytes = crate::genbank_bytes(&doc).unwrap();
        let record = super::super::one_genbank(&bytes).unwrap();
        assert!(!advertised(&record));
        assert_eq!(
            record.keywords.as_deref(),
            Some("biological; classification.")
        );
        assert!(record.references.is_empty());
    }

    #[test]
    fn named_segment_uses_its_document_palette_color_in_standard_exports() {
        let mut doc = Document::new("named default color", "AACCGGTTAACC").unwrap();
        for (id, label, start, end) in [(1, "First", 0, 2), (2, "Named", 2, 8)] {
            doc.features.push(Feature {
                id,
                label: label.into(),
                kind: "misc_feature".into(),
                parts: vec![Segment { start, end }],
                location: None,
                reverse: false,
                qualifiers: if id == 2 {
                    vec![("segment_name".into(), Some("Uncolored domain".into()))]
                } else {
                    vec![]
                },
                edited: false,
            });
        }
        assert_eq!(doc.features[1].effective_color_rgb(1), 0x65abe1);
        for format in [crate::SaveFormat::GenBank, crate::SaveFormat::Embl] {
            let bytes = crate::encode(&doc, format).unwrap();
            let reopened = crate::import("named", &bytes).unwrap();
            assert_eq!(reopened.features[1].display_color_rgb(), Some(0x65abe1));
            assert_eq!(
                reopened.features[1].display_segments()[0].name.as_deref(),
                Some("Uncolored domain")
            );
            assert_eq!(reopened.features[1].parts, doc.features[1].parts);
            assert_eq!(
                reopened.features[0].effective_color_rgb(0),
                doc.features[0].effective_color_rgb(0)
            );
        }
    }
}
