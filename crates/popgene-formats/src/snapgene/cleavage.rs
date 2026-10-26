//! Native `cleavageArrows` is a comma-separated list of interbase positions.
//! AutoSnapGene's published reader and original pUC19 example establish the
//! encoding: reverse AmpR has an arrow at 2417 between ranges 1626-2417 and
//! 2418-2486. Positions are on the sequence axis, independent of feature strand.
//! Reference: vendor-notices/autosnapgene-reference/SOURCE.json.
use super::{Feature, Node, Result, bad};
use std::collections::BTreeSet;

pub(super) fn read_attribute(node: &Node, length: usize) -> Result<Vec<usize>> {
    let value = node.attr("cleavageArrows");
    if value.is_empty() {
        return Ok(vec![]);
    }
    if value.len() > 1024 * 1024 {
        return Err(bad("Cleavage-arrow attribute exceeds 1 MiB"));
    }
    let mut sites = vec![];
    let mut seen = BTreeSet::new();
    for item in value.split(',') {
        let item = item.trim();
        if item.is_empty() || !item.bytes().all(|b| b.is_ascii_digit()) {
            return Err(bad("Invalid cleavage-arrow position"));
        }
        let site = item
            .parse::<usize>()
            .ok()
            .filter(|site| *site <= length)
            .ok_or_else(|| bad("Cleavage-arrow position is outside the sequence"))?;
        if !seen.insert(site) || sites.len() == 20_000 {
            return Err(bad("Duplicate or excessive cleavage-arrow positions"));
        }
        sites.push(site);
    }
    Ok(sites)
}

pub(super) fn append_note(note: &mut String, sites: &[usize]) {
    if sites.is_empty() {
        return;
    }
    if !note.is_empty() {
        note.push(' ');
    }
    note.push_str(if sites.len() == 1 {
        "Cleavage site after base "
    } else {
        "Cleavage sites after bases "
    });
    note.push_str(
        &sites
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join(", "),
    );
}

pub(super) fn apply(feature: &mut Feature, sites: &[usize]) -> Result<()> {
    if sites.is_empty() {
        return Ok(());
    }
    if !sites.iter().all(|site| {
        feature
            .parts
            .iter()
            .any(|part| part.start <= *site && *site <= part.end)
    }) {
        return Err(bad("Cleavage-arrow position is outside its feature"));
    }
    if feature.display_cleavage_sites() == sites {
        return Ok(());
    }
    // Keep every original qualifier. The canonical terminal note joins the
    // native arrow positions to any recognized display styling, so subsequent
    // coordinate edits use the same mapping as imported GenBank metadata.
    let mut note = feature.terminal_display_note().unwrap_or("").to_owned();
    for marker in ["Cleavage site after base ", "Cleavage sites after bases "] {
        if let Some(index) = note.rfind(marker) {
            note.truncate(index);
            note.truncate(note.trim_end().len());
        }
    }
    append_note(&mut note, sites);
    feature.qualifiers.push(("note".into(), Some(note)));
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::snapgene::{packet, read, write};

    fn fixture(arrows: &str, direction: &str, segments: &str, notes: &str) -> Vec<u8> {
        let mut bytes = vec![];
        packet(&mut bytes, 9, b"SnapGene\0\x01\0\x0f\0\x14").unwrap();
        packet(&mut bytes, 0, b"\x02ACGTACGTACGTACGTACGT").unwrap();
        packet(
            &mut bytes,
            10,
            format!(
                "<Features><Feature name=\"cleavage test\" type=\"misc_feature\" directionality=\"{direction}\" cleavageArrows=\"{arrows}\">{segments}{notes}</Feature></Features>"
            ).as_bytes(),
        ).unwrap();
        bytes
    }

    const SEGMENTS: &str = "<Segment range=\"1-10\" color=\"#ccffcc\" type=\"standard\"/><Segment range=\"11-20\" color=\"#ff0000\" type=\"standard\" name=\"signal\"/>";

    #[test]
    fn native_arrows_preserve_strand_segments_and_qualifiers() {
        for direction in ["1", "2"] {
            let bytes = fixture(
                "10",
                direction,
                SEGMENTS,
                "<Q name=\"note\"><V text=\"Original source prose\"/></Q>",
            );
            let mut doc = read("arrows.dna", &bytes).unwrap();
            assert_eq!(write(&doc).unwrap(), bytes);
            let original = doc.features[0].clone();
            assert_eq!(original.display_cleavage_sites(), [10]);
            assert_eq!(original.reverse, direction == "2");
            assert_eq!(
                original.display_segments()[1].name.as_deref(),
                Some("signal")
            );
            assert!(
                original
                    .qualifiers
                    .contains(&("note".into(), Some("Original source prose".into())))
            );
            doc.features[0].label = "renamed to force annotation reconstruction".into();
            for pass in 0..3 {
                doc.features[0].label = format!("renamed {pass}");
                let bytes = write(&doc).unwrap();
                let xml = crate::snapgene_preserve::packets(&bytes)
                    .unwrap()
                    .into_iter()
                    .find(|(kind, _)| *kind == 10)
                    .unwrap()
                    .1;
                assert!(
                    std::str::from_utf8(xml)
                        .unwrap()
                        .contains("cleavageArrows=\"10\"")
                );
                doc = read("arrows.dna", &bytes).unwrap();
                assert_eq!(doc.features[0].parts, original.parts);
                let mut expected = original.qualifiers.clone();
                let mut actual = doc.features[0].qualifiers.clone();
                expected.sort();
                actual.sort();
                assert_eq!(actual, expected);
                assert_eq!(doc.features[0].display_cleavage_sites(), [10]);
            }
        }
    }

    #[test]
    fn native_arrows_remap_after_insertion_and_reverse_complement() {
        let mut doc = read("arrows.dna", &fixture("0,10,20", "2", SEGMENTS, "")).unwrap();
        doc.replace(4, 4, "AAA").unwrap();
        assert_eq!(doc.features[0].display_cleavage_sites(), [0, 13, 23]);
        doc.revcomp().unwrap();
        assert_eq!(doc.features[0].display_cleavage_sites(), [23, 10, 0]);
        let restored = read("edited.dna", &write(&doc).unwrap()).unwrap();
        assert_eq!(restored.features[0].display_cleavage_sites(), [23, 10, 0]);
        assert_eq!(restored.features[0].parts, doc.features[0].parts);
        assert_eq!(restored.features[0].reverse, doc.features[0].reverse);
    }

    #[test]
    fn invalid_native_arrows_are_not_silently_dropped() {
        for arrows in [
            "-1",
            "21",
            "1,x",
            "1,",
            "10,10",
            "999999999999999999999999999999",
        ] {
            assert!(
                read("invalid.dna", &fixture(arrows, "1", SEGMENTS, "")).is_err(),
                "{arrows}"
            );
        }
        assert!(
            read(
                "outside.dna",
                &fixture("5", "1", "<Segment range=\"10-20\"/>", "")
            )
            .is_err()
        );
    }

    #[test]
    fn native_attribute_overrides_conflicting_note_without_erasing_it() {
        let notes = "<Q name=\"note\"><V text=\"color: #ccffcc; direction: LEFT Cleavage site after base 4\"/></Q>";
        let doc = read(
            "arrows.dna",
            &fixture(
                "10",
                "2",
                "<Segment range=\"1-20\" color=\"#ccffcc\"/>",
                notes,
            ),
        )
        .unwrap();
        assert_eq!(doc.features[0].display_cleavage_sites(), [10]);
        assert!(
            doc.features[0]
                .qualifiers
                .iter()
                .any(|(key, value)| key == "note"
                    && value.as_ref().is_some_and(|text| text.ends_with("base 4")))
        );
    }
}
