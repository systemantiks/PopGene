use super::*;

fn native(extra: &[(u8, Vec<u8>)]) -> Source {
    let doc = Document::new("stored records", "ACGTACGTACGT").unwrap();
    let mut bytes = super::super::write(&doc).unwrap();
    for (kind, body) in extra {
        super::super::packet(&mut bytes, *kind, body).unwrap();
    }
    super::super::read("records.dna", &bytes)
        .unwrap()
        .source
        .unwrap()
}

#[test]
fn independently_generated_ztr_can_be_opened_without_inventing_alignment_or_reversal() {
    let ztr =
        include_bytes!("../../../../../validation/fixtures/trace/extended/shared-huffman.ztr");
    let expected = crate::ztr::read("original.ztr", ztr).unwrap();
    let mut container = 1_u32.to_be_bytes().to_vec();
    super::super::packet(&mut container, 18, ztr).unwrap();
    super::super::packet(
        &mut container,
        8,
        b"<TraceProperties name=\"read one\"><Instrument>synthetic</Instrument></TraceProperties>",
    )
    .unwrap();
    let source = native(&[(16, container)]);
    let summaries = native_embedded_traces(&source).unwrap();
    assert_eq!(summaries.len(), 1);
    assert!(summaries[0].reverse);
    assert_eq!(summaries[0].byte_len, ztr.len());
    assert!(
        summaries[0]
            .properties
            .iter()
            .any(|(_, value)| value == "synthetic")
    );
    let mut trace = native_embedded_trace(&source, 0).unwrap();
    assert_eq!(trace.sequence, expected.sequence);
    assert_eq!(trace.trace, expected.trace);
    assert_eq!(trace.source.as_ref().unwrap().bytes.as_ref(), ztr);
    assert!(
        trace
            .warnings
            .iter()
            .any(|comment| comment.contains("orientation: reverse"))
    );
    assert_eq!(crate::encode(&trace, crate::SaveFormat::Ztr).unwrap(), ztr);
    let original_signals = trace.trace.as_ref().unwrap().channels.clone();
    let base = if trace.sequence.starts_with('A') {
        "N"
    } else {
        "A"
    };
    popgene_core::trace::replace_call(&mut trace, 0, base).unwrap();
    let rewritten = crate::encode(&trace, crate::SaveFormat::Ztr).unwrap();
    let edited = crate::ztr::read("edited.ztr", &rewritten).unwrap();
    assert!(edited.sequence.starts_with(base));
    assert_eq!(edited.trace.unwrap().channels, original_signals);
    assert!(native_embedded_trace(&source, 1).is_err());
    let fields = crate::source_metadata(&source).unwrap();
    assert!(
        fields
            .iter()
            .any(|(label, value)| label.contains("embedded trace 1") && value.contains("reverse"))
    );
}

#[test]
fn invalid_trace_flags_and_duplicate_records_are_explicit() {
    let source = native(&[(16, 2_u32.to_be_bytes().to_vec())]);
    assert!(native_embedded_traces(&source).is_err());
    let fields = crate::source_metadata(&source).unwrap();
    assert!(
        fields
            .iter()
            .any(|(label, value)| label.contains("decoding unavailable")
                && value.contains("direction flags"))
    );
    let ztr = b"\xaeZTR\r\n\x1a\n\x01\x02";
    let mut container = 0_u32.to_be_bytes().to_vec();
    super::super::packet(&mut container, 18, ztr).unwrap();
    super::super::packet(&mut container, 18, ztr).unwrap();
    assert!(trace_container(&container).is_err());
}

#[test]
fn documented_xml_and_rna_json_display_preserves_source_values_and_scope() {
    let json = br#"{"bondedPairs":[[22,30],[23,29]],"parameters":{"temperature":37},"probabilityMatrix":[[0,5,0.077]],"stats":{"freeEnergyOptimal":-75.7}}"#;
    let source = native(&[
        (14,b"<CustomEnzymeSets><CustomEnzymeSet type=\"1\" name=\"My panel\" enzymeNames=\"EcoRI BamHI\"/></CustomEnzymeSets>".to_vec()),
        (17,b"<AlignableSequences trimStringency=\"Medium\"><Sequence ID=\"7\" name=\"read\" sequence=\"ACGT\"/></AlignableSequences>".to_vec()),
        (28,b"<EnzymeVisibilities vals=\"EcoRI,1\"/>".to_vec()),
        (34,super::super::compression::xencode(json).unwrap()),
    ]);
    let fields = crate::source_metadata(&source).unwrap();
    for (label, value) in [
        (
            "Original CustomEnzymeSets/CustomEnzymeSet[1]/@enzymeNames",
            "EcoRI BamHI",
        ),
        ("Original AlignableSequences/Sequence[1]/@ID", "7"),
        ("Original EnzymeVisibilities/@vals", "EcoRI,1"),
        ("Original RNA structure/bondedPairs[0][0]", "22"),
        ("Original RNA structure/probabilityMatrix[0][2]", "0.077"),
        ("Original RNA structure/stats/freeEnergyOptimal", "-75.7"),
    ] {
        assert!(
            fields.iter().any(|(key, v)| key == label && v == value),
            "missing {label}"
        );
    }
    assert!(
        fields
            .iter()
            .any(|(_, value)| value.contains("Index origin") && value.contains("not inferred"))
    );
    let mut doc = super::super::read("records.dna", &source.bytes).unwrap();
    doc.description = "edited description".into();
    let rewritten = super::super::write(&doc).unwrap();
    let after = packets(&rewritten).unwrap();
    for (kind, body) in packets(&source.bytes)
        .unwrap()
        .into_iter()
        .filter(|(kind, _)| matches!(kind, 14 | 17 | 28 | 34))
    {
        assert!(after.iter().any(|(k, b)| *k == kind && *b == body));
    }
}

#[test]
fn structure_depth_and_source_display_limits_remain_explicit() {
    assert!(json_depth(format!("{}0{}", "[".repeat(65), "]".repeat(65)).as_bytes()).is_err());
    assert!(json_depth(br#"{"quoted":"[[[\\\""}"#).is_ok());
    let large = format!("<EnzymeVisibilities vals=\"{}\"/>", "α".repeat(10_000));
    let source = native(&[(28, large.into_bytes())]);
    let fields = crate::source_metadata(&source).unwrap();
    assert!(fields.iter().any(
        |(label, value)| label == "Original EnzymeVisibilities/@vals"
            && value.len() <= 4096
            && value.ends_with(" …")
    ));
    assert!(
        fields
            .iter()
            .any(|(label, _)| label == "Source display limit")
    );
    assert!(fields.iter().map(|(k, v)| k.len() + v.len()).sum::<usize>() < 1024 * 1024);
}

#[test]
fn failed_compressed_records_consume_the_aggregate_work_budget() {
    let mut records = vec![(34, b"invalid XZ".to_vec()); 4];
    records.push((
        34,
        super::super::compression::xencode(br#"{"revision":1}"#).unwrap(),
    ));
    let fields = crate::source_metadata(&native(&records)).unwrap();
    assert!(
        fields
            .iter()
            .any(|(_, value)| value.contains("Aggregate ancillary decoding budget is exhausted"))
    );
    assert!(
        !fields
            .iter()
            .any(|(key, _)| key == "Original RNA structure/revision")
    );
}
