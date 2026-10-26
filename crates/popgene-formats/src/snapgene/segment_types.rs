//! Preserve unsupported native segment type values as source provenance, not
//! an interpretation of their appearance or biological meaning.
use super::{Document, Feature, record_preservation::Records};
use std::collections::BTreeSet;

const QUALIFIER: &str = "snapgene_segment_type";
pub(super) const REBUILD_REFUSAL: &str = "Native export cannot rebuild an annotation associated with an unrecognized native segment type without losing its unsupported semantics. Keep its original XML unchanged or remove the affected annotation. After deleting one, save and reopen before adding a replacement. Unrelated metadata edits can retain the original record.";

fn unsupported(value: &str) -> bool {
    // An absent type is accepted by existing readers and observed fixtures.
    // Only these established values are interpreted; no unknown value is
    // treated as evidence for line, gap, or another native display convention.
    !matches!(value, "" | "standard" | "gap")
}

#[derive(Default)]
pub(super) struct SourceTypes<'a> {
    seen: BTreeSet<&'a str>,
}
impl<'a> SourceTypes<'a> {
    pub(super) fn retain(
        &mut self,
        value: &'a str,
        qualifiers: &mut Vec<(String, Option<String>)>,
    ) -> bool {
        if !unsupported(value) {
            return false;
        }
        // Membership checks do not scan the growing qualifier list. Borrow
        // the source strings and retain first-occurrence order in the model.
        if self.seen.insert(value) {
            qualifiers.push((QUALIFIER.into(), Some(value.into())));
        }
        true
    }
}

pub(super) fn import_warning(count: usize) -> String {
    // One aggregate warning, without copying unbounded source names/values.
    format!(
        "{count} feature segment(s) have unrecognized native segment types. Their ranges are displayed using a fallback; the types' appearance and biological semantics are not interpreted. Original type values remain source provenance in /{QUALIFIER}, and original XML is retained. Native saves that must rebuild an affected annotation are refused."
    )
}

fn has_unsupported_type(feature: &Feature) -> bool {
    // Check every value: a same-named user qualifier containing 'standard'
    // must not mask an unsupported value captured from the native attribute.
    feature
        .qualifiers
        .iter()
        .any(|(key, value)| key == QUALIFIER && value.as_deref().is_some_and(unsupported))
}

pub(super) fn requires_rebuild(
    current: &Document,
    original: Option<&Document>,
    preserved: &Records,
) -> bool {
    let original_unsupported = original
        .into_iter()
        .flat_map(|doc| doc.features.iter().chain(&doc.duplex_annotations.features))
        .filter(|feature| has_unsupported_type(feature))
        .map(|feature| feature.id)
        .collect::<BTreeSet<_>>();
    current
        .features
        .iter()
        .chain(&current.duplex_annotations.features)
        .any(|feature| {
            !preserved.contains(feature.id)
                && (has_unsupported_type(feature) || original_unsupported.contains(&feature.id))
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::snapgene::{escaped, packet, read, write, write_warnings, xml};

    fn fixture(features: &str, properties: Option<&[u8]>) -> Vec<u8> {
        let mut bytes = Vec::new();
        packet(&mut bytes, 9, b"SnapGene\0\x01\0\x0f\0\x14").unwrap();
        packet(&mut bytes, 0, b"\x02AAACCCGGGTTT").unwrap();
        packet(&mut bytes, 10, features.as_bytes()).unwrap();
        if let Some(properties) = properties {
            packet(&mut bytes, 8, properties).unwrap();
        }
        bytes
    }

    fn unknown_fixture() -> Vec<u8> {
        fixture(
            "<Features><Feature name=\"opaque\" type=\"misc_feature\"><Segment range=\"1-3\" type=\"future-display\" color=\"#123456\" auxiliaryStyle=\"retain exactly\"/><Q name=\"snapgene_segment_type\"><V text=\"standard\"/></Q></Feature><Feature name=\"neighbor\" type=\"misc_feature\"><Segment range=\"7-9\" type=\"standard\" color=\"#654321\"/></Feature></Features>",
            None,
        )
    }

    fn feature_packet(bytes: &[u8]) -> Vec<u8> {
        crate::snapgene_preserve::packets(bytes)
            .unwrap()
            .into_iter()
            .find(|(kind, _)| *kind == 10)
            .unwrap()
            .1
            .to_vec()
    }

    fn assert_refused(doc: &Document) {
        assert_eq!(write(doc).unwrap_err().to_string(), REBUILD_REFUSAL);
        assert!(
            write_warnings(doc)
                .iter()
                .any(|warning| warning == REBUILD_REFUSAL)
        );
    }

    #[test]
    fn warning_is_aggregate_and_bounded_without_truncating_source_provenance() {
        let value = format!("future-{}<&\"", "界".repeat(2048));
        let segment = format!(
            "<Segment range=\"1-3\" type=\"{}\" color=\"#123456\"/>",
            escaped(&value)
        );
        let bytes = fixture(
            &format!(
                "<Features><Feature name=\"opaque\" type=\"misc_feature\">{}</Feature></Features>",
                segment.repeat(64)
            ),
            None,
        );
        let doc = read("unknown.dna", &bytes).unwrap();
        let warnings = doc
            .warnings
            .iter()
            .filter(|w| w.contains("unrecognized native segment types"))
            .collect::<Vec<_>>();
        assert_eq!(warnings.len(), 1);
        assert!(warnings[0].starts_with("64 feature segment(s)"));
        assert!(warnings[0].len() < 512);
        assert!(warnings[0].contains("semantics are not interpreted"));
        assert_eq!(
            doc.features[0]
                .qualifiers
                .iter()
                .filter(|(key, _)| key == QUALIFIER)
                .count(),
            1
        );
        assert_eq!(doc.features[0].qualifier(QUALIFIER), Some(value.as_str()));
        assert_eq!(write(&doc).unwrap(), bytes);

        let segments = ["zeta", "alpha", "zeta", "beta", "alpha"]
            .into_iter()
            .map(|value| format!("<Segment range=\"1-3\" type=\"{value}\"/>"))
            .collect::<String>();
        let bytes = fixture(
            &format!(
                "<Features><Feature name=\"ordered provenance\">{segments}</Feature></Features>"
            ),
            None,
        );
        let doc = read("ordered.dna", &bytes).unwrap();
        assert_eq!(
            doc.features[0]
                .qualifiers
                .iter()
                .filter(|(key, _)| key == QUALIFIER)
                .map(|(_, value)| value.as_deref().unwrap())
                .collect::<Vec<_>>(),
            ["zeta", "alpha", "beta"]
        );
    }

    #[test]
    fn untouched_metadata_and_neighbor_edits_retain_opaque_xml() {
        let bytes = unknown_fixture();
        let mut doc = read("unknown.dna", &bytes).unwrap();
        assert!(has_unsupported_type(&doc.features[0]));
        // A same-named ordinary qualifier cannot override captured provenance.
        assert_eq!(
            doc.features[0]
                .qualifiers
                .iter()
                .filter(|(key, _)| key == QUALIFIER)
                .count(),
            2
        );
        assert_eq!(write(&doc).unwrap(), bytes);
        doc.description = "Metadata only".into();
        let metadata = write(&doc).unwrap();
        assert_eq!(feature_packet(&metadata), feature_packet(&bytes));
        assert!(!write_warnings(&doc).iter().any(|w| w == REBUILD_REFUSAL));

        doc.features[1].label = "renamed neighbor".into();
        let neighbor = write(&doc).unwrap();
        let root = xml(&feature_packet(&neighbor)).unwrap();
        let segment = root.children[0]
            .children
            .iter()
            .find(|n| n.name == "Segment")
            .unwrap();
        assert_eq!(segment.attr("type"), "future-display");
        assert_eq!(segment.attr("auxiliaryStyle"), "retain exactly");
        assert_eq!(segment.attr("range"), "1-3");
        assert_eq!(root.children[1].attr("name"), "renamed neighbor");
        assert!(has_unsupported_type(
            &read("neighbor.dna", &neighbor).unwrap().features[0]
        ));
    }

    #[test]
    fn affected_rebuild_is_refused_but_intentional_feature_deletion_is_allowed() {
        let original = read("unknown.dna", &unknown_fixture()).unwrap();
        let mut renamed = original.clone();
        renamed.features[0].label = "changed opaque annotation".into();
        assert_refused(&renamed);

        let mut removed_provenance = original.clone();
        removed_provenance.features[0]
            .qualifiers
            .retain(|(key, _)| key != QUALIFIER);
        assert_refused(&removed_provenance);

        let mut changed_molecule = original.clone();
        changed_molecule.replace(10, 11, "A").unwrap();
        assert_refused(&changed_molecule);

        let mut deleted = original;
        deleted.features.remove(0);
        let reopened = read("deleted.dna", &write(&deleted).unwrap()).unwrap();
        assert_eq!(reopened.features.len(), 1);
        assert_eq!(reopened.features[0].label, "neighbor");
        assert!(!has_unsupported_type(&reopened.features[0]));
    }

    #[test]
    fn source_provenance_survives_workflow_reversal_detachment_and_standard_conversion() {
        let original = read("unknown.dna", &unknown_fixture()).unwrap();
        let product = popgene_core::workflow::extract_region(&original, 0, 12, "product").unwrap();
        assert!(product.source.is_none());
        assert!(has_unsupported_type(&product.features[0]));
        assert_refused(&product);

        let mut reversed = original.clone();
        reversed.revcomp().unwrap();
        assert!(has_unsupported_type(&reversed.features[0]));
        assert_refused(&reversed);

        let mut detached = original.clone();
        detached.source = None;
        assert_refused(&detached);

        let gb = crate::encode(&original, crate::SaveFormat::GenBank).unwrap();
        let converted = crate::import("converted.gb", &gb).unwrap();
        assert!(has_unsupported_type(&converted.features[0]));
        assert_refused(&converted);
    }

    #[test]
    fn recycled_identifier_is_conservative_until_deletion_is_saved_and_reopened() {
        let bytes = fixture(
            "<Features><Feature name=\"opaque\" type=\"misc_feature\"><Segment range=\"1-3\" type=\"future-display\"/></Feature></Features>",
            None,
        );
        let mut deleted = read("unknown.dna", &bytes).unwrap();
        let original_id = deleted.features[0].id;
        deleted.features.clear();
        let mut replaced = deleted.clone();
        replaced
            .annotate(
                "new annotation".into(),
                popgene_core::Segment { start: 3, end: 6 },
                false,
            )
            .unwrap();
        assert_eq!(replaced.features[0].id, original_id);
        assert!(!has_unsupported_type(&replaced.features[0]));
        // The current/source models cannot distinguish reused identity from
        // clearing provenance while editing the old record. Keep the guard.
        assert_refused(&replaced);

        let saved_deletion = write(&deleted).unwrap();
        let mut reopened = read("deleted.dna", &saved_deletion).unwrap();
        reopened
            .annotate(
                "new annotation".into(),
                popgene_core::Segment { start: 3, end: 6 },
                false,
            )
            .unwrap();
        let restored = read("replacement.dna", &write(&reopened).unwrap()).unwrap();
        assert_eq!(restored.features.len(), 1);
        assert_eq!(restored.features[0].label, "new annotation");
        assert_eq!(
            restored.feature_sequence(&restored.features[0]).unwrap(),
            "CCC"
        );
    }

    #[test]
    fn unsupported_full_duplex_annotation_also_requires_preserved_xml() {
        let original = unknown_fixture();
        let bytes = fixture(
            std::str::from_utf8(&feature_packet(&original)).unwrap(),
            Some(b"<AdditionalSequenceProperties><UpstreamStickiness>-3</UpstreamStickiness><DownstreamStickiness>0</DownstreamStickiness></AdditionalSequenceProperties>"),
        );
        let mut doc = read("duplex.dna", &bytes).unwrap();
        assert_eq!(doc.duplex_annotations.features.len(), 1);
        assert!(has_unsupported_type(&doc.duplex_annotations.features[0]));
        doc.description = "Metadata only".into();
        assert_eq!(
            feature_packet(&write(&doc).unwrap()),
            feature_packet(&bytes)
        );
        doc.duplex_annotations.features[0].label = "changed".into();
        assert_refused(&doc);
    }

    #[test]
    fn established_standard_gap_and_missing_types_keep_existing_behavior() {
        let bytes = fixture(
            "<Features><Feature name=\"known\" type=\"misc_feature\"><Segment range=\"1-3\" color=\"#123456\"/><Segment range=\"4-6\" type=\"gap\"/><Segment range=\"7-9\" type=\"standard\" color=\"#654321\"/><Segment range=\"10-12\" type=\"\" color=\"#abcdef\"/></Feature></Features>",
            None,
        );
        let mut doc = read("known.dna", &bytes).unwrap();
        assert!(
            !doc.warnings
                .iter()
                .any(|w| w.contains("unrecognized native segment types"))
        );
        assert!(!has_unsupported_type(&doc.features[0]));
        assert_eq!(doc.feature_sequence(&doc.features[0]).unwrap(), "AAAGGGTTT");
        assert_eq!(doc.features[0].display_gaps().len(), 1);
        doc.features[0].label = "rebuilt known annotation".into();
        let reopened = read("rebuilt.dna", &write(&doc).unwrap()).unwrap();
        assert_eq!(
            reopened.feature_sequence(&reopened.features[0]).unwrap(),
            "AAAGGGTTT"
        );
        assert_eq!(
            reopened.features[0].display_gaps(),
            doc.features[0].display_gaps()
        );
    }
}
