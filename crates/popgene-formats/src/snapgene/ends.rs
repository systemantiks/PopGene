//! Native sequence blocks span the union of both strands. The shared model stores
//! the actual top strand and the actual opposite strand separately.
use super::{Document, Molecule, Node, Result, bad, xml};
use popgene_core::{
    restriction::{EndKind, FragmentEnd},
    reverse_complement,
    workflow::DuplexEnds,
};
fn signed(root: &Node, key: &str) -> Result<isize> {
    root.child(key).map_or(Ok(0), |node| {
        node.text.parse().map_err(|_| bad("Invalid end stickiness"))
    })
}
fn shift_features(doc: &mut Document, delta: isize, lower: usize, upper: usize) -> Result<()> {
    for feature in &mut doc.features {
        if feature
            .parts
            .iter()
            .any(|part| part.start < lower || part.end > upper)
        {
            return Err(bad(
                "Native annotation lies outside its coordinate projection",
            ));
        }
        popgene_core::workflow::shift_feature(feature, delta)?;
    }
    Ok(())
}
/// The native union-axis reader partitions annotations using this exact rule.
pub(super) fn outside_top(feature: &popgene_core::Feature, ts: usize, te: usize) -> bool {
    feature
        .parts
        .iter()
        .any(|part| part.start < ts || part.end > te)
        || feature
            .display_gaps()
            .iter()
            .any(|part| part.range.start < ts || part.range.end > te)
        || feature
            .display_cleavage_sites()
            .iter()
            .any(|site| *site < ts || *site > te)
        || feature.location.as_ref().is_some_and(|location| {
            location
                .local_boundaries()
                .iter()
                .any(|site| *site <= ts || *site >= te)
        })
}
pub(super) fn read(doc: &mut Document, blocks: &[(u8, &[u8])]) -> Result<()> {
    let properties = blocks
        .iter()
        .filter(|(kind, _)| *kind == 8)
        .collect::<Vec<_>>();
    if properties.len() > 1 {
        return Err(bad("Duplicate end properties"));
    }
    let Some((_, body)) = properties.first() else {
        return Ok(());
    };
    let root = xml(body)?;
    if root.name != "AdditionalSequenceProperties" {
        return Err(bad("Invalid end properties root"));
    }
    let left = signed(&root, "UpstreamStickiness")?;
    let right = signed(&root, "DownstreamStickiness")?;
    let modification = |key| root.child(key).map(|node| node.text.clone());
    let lm = modification("UpstreamModification");
    let rm = modification("DownstreamModification");
    if left == 0
        && right == 0
        && [&lm, &rm]
            .into_iter()
            .all(|m| m.as_deref().is_none_or(|v| v == "Unmodified"))
    {
        return Ok(());
    }
    if matches!(doc.molecule, Molecule::Dna | Molecule::Rna)
        && doc.circular
        && left == 0
        && right == 0
    {
        // Observed native DNA and RNA circles retain terminal phosphorylation
        // fields. They cannot describe exposed ends on a closed circle;
        // keep the source record without inventing a linear duplex geometry.
        doc.warnings.push("Original terminal-modification fields are retained as source metadata, but this circular molecule has no exposed ends; those fields are not applied as current end chemistry.".into());
        return Ok(());
    }
    if doc.molecule != Molecule::Dna || doc.circular || doc.single_stranded {
        return Err(bad(
            "Nonblunt end geometry requires linear double-stranded DNA",
        ));
    }
    let l = left.unsigned_abs();
    let r = right.unsigned_abs();
    let n = doc.sequence.len();
    if l.checked_add(r).is_none_or(|v| v >= n) {
        return Err(bad("Cohesive ends leave no paired sequence"));
    }
    let ts = if left < 0 { l } else { 0 };
    let te = n - if right > 0 { r } else { 0 };
    let bs = if left > 0 { l } else { 0 };
    let be = n - if right < 0 { r } else { 0 };
    let end = |signed: isize, overhang: String| FragmentEnd {
        kind: if signed > 0 {
            EndKind::FivePrime
        } else if signed < 0 {
            EndKind::ThreePrime
        } else {
            EndKind::Blunt
        },
        overhang,
        enzymes: vec![],
    };
    let left_overhang = if left < 0 {
        reverse_complement(&doc.sequence[..l])
    } else {
        doc.sequence[..l].to_owned()
    };
    let right_overhang = if right > 0 {
        reverse_complement(&doc.sequence[n - r..])
    } else {
        doc.sequence[n - r..].to_owned()
    };
    let ends = DuplexEnds {
        left: end(left, left_overhang),
        right: end(right, right_overhang),
        opposite_sequence: reverse_complement(&doc.sequence[bs..be]),
        left_modification: lm,
        right_modification: rm,
    };
    let mut top_features = Vec::with_capacity(doc.features.len());
    for feature in std::mem::take(&mut doc.features) {
        let outside = outside_top(&feature, ts, te);
        if outside {
            doc.duplex_annotations.features.push(feature);
        } else {
            top_features.push(feature);
        }
    }
    doc.features = top_features;
    if ts != 0 {
        shift_features(doc, -(ts as isize), ts, te)?;
    }
    let mut top_colors = Vec::with_capacity(doc.sequence_colors.len());
    for mut color in std::mem::take(&mut doc.sequence_colors) {
        if color.range.start < ts || color.range.end > te {
            doc.duplex_annotations.sequence_colors.push(color);
        } else {
            color.range.start -= ts;
            color.range.end -= ts;
            top_colors.push(color);
        }
    }
    doc.sequence_colors = top_colors;
    doc.sequence = doc.sequence[ts..te].into();
    ends.validate(&doc.sequence)?;
    doc.duplex_ends = Some(ends);
    Ok(())
}
pub(super) fn encode(doc: &Document) -> Result<(Document, Option<Node>)> {
    let Some(ends) = &doc.duplex_ends else {
        return Ok((doc.clone(), None));
    };
    ends.validate(&doc.sequence)?;
    let signed = |end: &FragmentEnd| match end.kind {
        EndKind::FivePrime => end.overhang.len() as isize,
        EndKind::ThreePrime => -(end.overhang.len() as isize),
        _ => 0,
    };
    let left = signed(&ends.left);
    let right = signed(&ends.right);
    let prefix = if left < 0 {
        reverse_complement(&ends.left.overhang)
    } else {
        String::new()
    };
    let suffix = if right > 0 {
        reverse_complement(&ends.right.overhang)
    } else {
        String::new()
    };
    let mut native = doc.clone();
    if !prefix.is_empty() {
        shift_features(&mut native, prefix.len() as isize, 0, doc.sequence.len())?;
    }
    for color in &mut native.sequence_colors {
        color.range.start += prefix.len();
        color.range.end += prefix.len();
    }
    native
        .features
        .append(&mut native.duplex_annotations.features);
    native
        .sequence_colors
        .append(&mut native.duplex_annotations.sequence_colors);
    native.sequence = format!("{prefix}{}{suffix}", doc.sequence).into();
    let mut node = Node {
        name: "AdditionalSequenceProperties".into(),
        ..Default::default()
    };
    for (key, text) in [
        ("UpstreamStickiness", left.to_string()),
        ("DownstreamStickiness", right.to_string()),
        (
            "UpstreamModification",
            ends.left_modification
                .clone()
                .unwrap_or_else(|| "Unmodified".into()),
        ),
        (
            "DownstreamModification",
            ends.right_modification
                .clone()
                .unwrap_or_else(|| "Unmodified".into()),
        ),
    ] {
        node.children.push(Node {
            name: key.into(),
            text,
            ..Default::default()
        });
    }
    Ok((native, Some(node)))
}
#[cfg(test)]
mod tests {
    use super::*;
    use popgene_core::Segment;
    #[test]
    fn circular_zero_stickiness_retains_inactive_terminal_fields() {
        use crate::snapgene::{packet, read, write};
        let properties = b"<AdditionalSequenceProperties><UpstreamStickiness>0</UpstreamStickiness><DownstreamStickiness>0</DownstreamStickiness><UpstreamModification>FivePrimePhosphorylated</UpstreamModification><DownstreamModification>FivePrimePhosphorylated</DownstreamModification></AdditionalSequenceProperties>";
        let mut bytes = vec![];
        packet(&mut bytes, 9, b"SnapGene\0\x01\0\x0f\0\x14").unwrap();
        packet(&mut bytes, 0, b"\x03AAACCCGGGTTT").unwrap();
        packet(&mut bytes, 8, properties).unwrap();
        let mut doc = read("circular.dna", &bytes).unwrap();
        assert!(doc.circular);
        assert!(doc.duplex_ends.is_none());
        assert!(doc.warnings.iter().any(|w| w.contains("no exposed ends")));
        assert_eq!(write(&doc).unwrap(), bytes);
        doc.description = "Metadata edit".into();
        let edited = write(&doc).unwrap();
        let packets = crate::snapgene_preserve::packets(&edited).unwrap();
        assert_eq!(packets.iter().find(|(k, _)| *k == 8).unwrap().1, properties);
        assert!(read("metadata.dna", &edited).unwrap().duplex_ends.is_none());
        doc.replace(0, 1, "C").unwrap();
        let edited = write(&doc).unwrap();
        assert!(
            !crate::snapgene_preserve::packets(&edited)
                .unwrap()
                .iter()
                .any(|(k, _)| *k == 8)
        );
        let source = read("sequence.dna", &edited).unwrap().source.unwrap();
        let historical = crate::native_history_document(&source, 0).unwrap();
        assert!(
            crate::source_metadata(historical.source.as_ref().unwrap())
                .unwrap()
                .iter()
                .any(|(_, value)| value.contains("FivePrimePhosphorylated"))
        );
        let invalid = String::from_utf8(properties.to_vec())
            .unwrap()
            .replace("<UpstreamStickiness>0", "<UpstreamStickiness>1");
        let mut nonzero = vec![];
        packet(&mut nonzero, 9, b"SnapGene\0\x01\0\x0f\0\x14").unwrap();
        packet(&mut nonzero, 0, b"\x03AAACCCGGGTTT").unwrap();
        packet(&mut nonzero, 8, invalid.as_bytes()).unwrap();
        assert!(read("nonzero.dna", &nonzero).is_err());
    }
    #[test]
    fn circular_rna_zero_stickiness_preserves_inactive_fields_without_inventing_duplex() {
        use crate::snapgene::{packet, read, write};
        let properties = "<AdditionalSequenceProperties><UpstreamStickiness>0</UpstreamStickiness><DownstreamStickiness>0</DownstreamStickiness><UpstreamModification>FivePrimePhosphorylated</UpstreamModification><DownstreamModification>FivePrimePhosphorylated</DownstreamModification></AdditionalSequenceProperties>";
        let native = |flags: u8, properties: &str| {
            let mut bytes = vec![];
            packet(&mut bytes, 9, b"SnapGene\0\x07\0\x0f\0\x14").unwrap();
            let mut sequence = vec![flags];
            sequence.extend_from_slice(b"AAACCCGGGUUU");
            packet(&mut bytes, 32, &sequence).unwrap();
            packet(&mut bytes, 8, properties.as_bytes()).unwrap();
            bytes
        };
        for flags in [1, 3] {
            let bytes = native(flags, properties);
            let mut doc = read("circular.rna", &bytes).unwrap();
            assert_eq!(doc.molecule, Molecule::Rna);
            assert!(doc.circular);
            assert_eq!(doc.single_stranded, flags == 1);
            assert_eq!(doc.sequence.as_ref(), "AAACCCGGGUUU");
            assert!(doc.duplex_ends.is_none());
            assert!(doc.warnings.iter().any(|w| w.contains("no exposed ends")));
            assert_eq!(write(&doc).unwrap(), bytes);
            doc.description = "Unrelated metadata edit".into();
            let edited = write(&doc).unwrap();
            let packets = crate::snapgene_preserve::packets(&edited).unwrap();
            assert_eq!(
                packets.iter().find(|(kind, _)| *kind == 8).unwrap().1,
                properties.as_bytes()
            );
            let reopened = read("metadata.rna", &edited).unwrap();
            assert_eq!(reopened.sequence, doc.sequence);
            assert_eq!(reopened.single_stranded, doc.single_stranded);
            assert!(reopened.duplex_ends.is_none());
            for key in ["UpstreamStickiness", "DownstreamStickiness"] {
                for value in [-1, 1] {
                    let nonzero =
                        properties.replace(&format!("<{key}>0"), &format!("<{key}>{value}"));
                    assert!(read("nonzero.rna", &native(flags, &nonzero)).is_err());
                }
            }
        }
        // The change does not interpret terminal chemistry on linear RNA.
        assert!(read("linear.rna", &native(2, properties)).is_err());
    }
    #[test]
    fn bottom_only_and_cross_boundary_annotations_keep_complete_extent_coordinates() {
        let mut d = Document::new("extent", "TAACCCGGGTTT").unwrap();
        for (index, (label, start, end)) in [
            ("left only", 0, 1),
            ("cross boundary", 0, 4),
            ("right only", 10, 12),
            ("paired", 4, 7),
        ]
        .into_iter()
        .enumerate()
        {
            d.features.push(popgene_core::Feature {
                id: index as u64 + 1,
                label: label.into(),
                kind: "misc_feature".into(),
                parts: vec![Segment { start, end }],
                reverse: true,
                location: None,
                qualifiers: vec![],
                edited: false,
            });
        }
        d.sequence_colors
            .push(popgene_core::metadata::SequenceColor {
                range: Segment { start: 0, end: 3 },
                reverse: true,
                color_rgb: 0xff00ff,
            });
        let xml=b"<AdditionalSequenceProperties><UpstreamStickiness>-1</UpstreamStickiness><DownstreamStickiness>2</DownstreamStickiness></AdditionalSequenceProperties>";
        read(&mut d, &[(8, xml)]).unwrap();
        d.validate().unwrap();
        assert_eq!(d.sequence.as_ref(), "AACCCGGGT");
        assert_eq!(d.features.len(), 1);
        assert_eq!(d.features[0].parts[0], Segment { start: 3, end: 6 });
        assert_eq!(d.duplex_annotations.features.len(), 3);
        assert_eq!(
            d.duplex_annotations.features[1].parts[0],
            Segment { start: 0, end: 4 }
        );
        assert_eq!(
            d.duplex_annotations.sequence_colors[0].range,
            Segment { start: 0, end: 3 }
        );
        let geometry = d.duplex_ends.as_ref().unwrap().geometry(d.sequence.len());
        assert_eq!(geometry.top, Segment { start: 1, end: 10 });
        assert_eq!(geometry.bottom, Segment { start: 0, end: 12 });
        assert_eq!(
            d.duplex_ends
                .as_ref()
                .unwrap()
                .reference_sequence(&d.sequence)
                .unwrap(),
            "TAACCCGGGTTT"
        );
        d.replace(4, 4, "AA").unwrap();
        assert_eq!(
            d.duplex_annotations.features[2].parts[0],
            Segment { start: 12, end: 14 }
        );
        let bytes = super::super::write(&d).unwrap();
        let reread = super::super::read("extent.dna", &bytes).unwrap();
        assert_eq!(reread.sequence, d.sequence);
        assert_eq!(
            reread.duplex_annotations.sequence_colors,
            d.duplex_annotations.sequence_colors
        );
        assert_eq!(
            reread
                .duplex_annotations
                .features
                .iter()
                .map(|f| (&f.label, &f.parts))
                .collect::<Vec<_>>(),
            d.duplex_annotations
                .features
                .iter()
                .map(|f| (&f.label, &f.parts))
                .collect::<Vec<_>>()
        );
        assert!(
            crate::save_warnings(&reread, crate::SaveFormat::GenBank)
                .iter()
                .any(|warning| warning.contains("readable PopGene comment extension"))
        );
        let standard = crate::encode(&reread, crate::SaveFormat::GenBank).unwrap();
        let standard = crate::import("extent.gb", &standard).unwrap();
        assert_eq!(standard.sequence, reread.sequence);
        assert_eq!(standard.duplex_ends, reread.duplex_ends);
        assert_eq!(standard.duplex_annotations, reread.duplex_annotations);
        assert_eq!(standard.sequence_colors, reread.sequence_colors);
    }
    #[test]
    fn observed_signed_geometry_roundtrips_actual_strands_and_feature_coordinates() {
        for (left, right) in [(3, 0), (-1, 0), (-1, -1), (0, 2), (2, -2)] {
            let mut d = Document::new("ends", "TAACCCGGGTTT").unwrap();
            d.features.push(popgene_core::Feature {
                id: 1,
                label: "paired".into(),
                kind: "misc_feature".into(),
                parts: vec![Segment { start: 4, end: 7 }],
                reverse: false,
                location: None,
                qualifiers: vec![],
                edited: false,
            });
            let xml = format!(
                "<AdditionalSequenceProperties><UpstreamStickiness>{left}</UpstreamStickiness><DownstreamStickiness>{right}</DownstreamStickiness></AdditionalSequenceProperties>"
            );
            read(&mut d, &[(8, xml.as_bytes())]).unwrap();
            d.validate().unwrap();
            let (mut native, node) = encode(&d).unwrap();
            assert_eq!(native.sequence.as_ref(), "TAACCCGGGTTT");
            assert_eq!(native.features[0].parts, vec![Segment { start: 4, end: 7 }]);
            let mut serialized = String::new();
            super::super::serialize_node(&node.unwrap(), &mut serialized);
            native.duplex_ends = None;
            read(&mut native, &[(8, serialized.as_bytes())]).unwrap();
            assert_eq!(native.sequence, d.sequence);
            assert_eq!(native.features, d.features);
            assert_eq!(
                native.duplex_ends.as_ref().unwrap().opposite_sequence,
                d.duplex_ends.as_ref().unwrap().opposite_sequence
            );
        }
    }
}
