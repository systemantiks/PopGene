use super::*;
use crate::Segment;
fn cds(sequence: &str, reverse: bool) -> (Document, Feature) {
    let document = Document::new("translation", sequence).unwrap();
    let feature = Feature {
        location: None,
        id: 1,
        label: "CDS".into(),
        kind: "CDS".into(),
        parts: vec![Segment {
            start: 0,
            end: sequence.len(),
        }],
        reverse,
        qualifiers: vec![],
        edited: false,
    };
    (document, feature)
}
fn qualify(feature: &mut Feature, key: &str, value: &str) {
    feature.qualifiers.push((key.into(), Some(value.into())));
}
#[test]
fn additional_codes_and_context_dependent_termination() {
    for (id, codon, expected) in [
        (3, *b"CTG", b'T'),
        (6, *b"TAA", b'Q'),
        (10, *b"TGA", b'C'),
        (12, *b"CTG", b'S'),
        (22, *b"TCA", b'*'),
        (25, *b"TGA", b'G'),
        (26, *b"CTG", b'A'),
        (29, *b"TAG", b'Y'),
        (30, *b"TAA", b'E'),
        (32, *b"TAG", b'W'),
        (33, *b"TAA", b'Y'),
    ] {
        assert_eq!(
            translate_codon(codon, GeneticCode::from_id(id).unwrap()),
            expected,
            "table {id}"
        );
    }
    for id in [27, 28, 31] {
        let (doc, mut feature) = cds(if id == 31 { "ATGTAATAA" } else { "ATGTGATGA" }, false);
        qualify(&mut feature, "transl_table", &id.to_string());
        let expected = if id == 31 { b"ME*" } else { b"MW*" };
        assert_eq!(translate_cds(&doc, &feature).unwrap().amino_acids, expected);
        assert!(
            GeneticCode::from_id(id)
                .unwrap()
                .has_context_dependent_stops()
        );
    }
}
#[test]
fn exact_exceptions_joined_reverse_partial_stops_and_conflicts() {
    let (doc, mut feature) = cds("ATGTGATAGTAA", false);
    qualify(&mut feature, "transl_except", "(pos:4..6,aa:Sec)");
    qualify(&mut feature, "transl_except", "(pos:7..9,aa:Pyl)");
    let t = translate_cds(&doc, &feature).unwrap();
    assert_eq!(t.amino_acids, b"MUO*");
    assert_eq!(t.exception_codons, vec![1, 2]);
    qualify(&mut feature, "transl_except", "(pos:4..6,aa:Trp)");
    assert!(translate_cds(&doc, &feature).is_err());
    let (doc, mut feature) = cds("TCACAT", true);
    qualify(
        &mut feature,
        "transl_except",
        "(pos:complement(1..3),aa:Sec)",
    );
    assert_eq!(translate_cds(&doc, &feature).unwrap().amino_acids, b"MU");
    let (doc, mut feature) = cds("ATGTCCGA", false);
    feature.parts = vec![Segment { start: 0, end: 4 }, Segment { start: 6, end: 8 }];
    qualify(&mut feature, "transl_except", "(pos:join(4,7..8),aa:Sec)");
    qualify(&mut feature, "ribosomal_slippage", "");
    assert_eq!(translate_cds(&doc, &feature).unwrap().amino_acids, b"MU");
    for (sequence, location) in [("ATGT", "4"), ("ATGTA", "4..5")] {
        let (doc, mut feature) = cds(sequence, false);
        qualify(
            &mut feature,
            "transl_except",
            &format!("(pos:{location},aa:TERM)"),
        );
        let t = translate_cds(&doc, &feature).unwrap();
        assert_eq!(t.amino_acids, b"M*");
        assert!(t.terminal_partial_stop);
        assert_eq!(t.trailing_bases, 0);
    }
    let (doc, mut feature) = cds("ATGG", false);
    qualify(&mut feature, "transl_except", "(pos:4,aa:TERM)");
    assert!(translate_cds(&doc, &feature).is_err());
}
#[test]
fn exceptions_follow_coordinate_changes_and_refuse_overlapping_edits() {
    let (mut doc, mut feature) = cds("ATGTGATAA", false);
    qualify(&mut feature, "transl_except", "(pos:4..6,aa:Sec)");
    doc.features.push(feature);
    doc.replace(0, 0, "C").unwrap();
    assert_eq!(
        translate_cds(&doc, &doc.features[0]).unwrap().amino_acids,
        b"MU*"
    );
    assert_eq!(
        doc.features[0].qualifier("transl_except"),
        Some("(pos:5..7,aa:Sec)")
    );
    doc.revcomp().unwrap();
    assert_eq!(
        translate_cds(&doc, &doc.features[0]).unwrap().amino_acids,
        b"MU*"
    );
    let p = doc.features[0].parts[0].start + 3;
    doc.replace(p, p + 1, "A").unwrap();
    assert!(translate_cds(&doc, &doc.features[0]).is_err());
}
#[test]
fn circular_frames_project_origin_codons_in_both_directions() {
    let sequence = "TGAA";
    let linear = frame_codons_with_options(sequence, 1, 0..4, GeneticCode::STANDARD, false);
    assert_eq!(linear.len(), 1);
    let circular = frame_codons_with_options(sequence, 1, 0..4, GeneticCode::STANDARD, true);
    assert_eq!(
        circular
            .iter()
            .map(|c| (c.center, c.amino_acid))
            .collect::<Vec<_>>(),
        vec![(0, b'M'), (1, b'*')]
    );
    assert_eq!(circular[0].positions, [3, 0, 1]);
    assert_eq!(
        circular[0].selection(),
        vec![Segment { start: 0, end: 2 }, Segment { start: 3, end: 4 }]
    );
    let reverse = frame_codons_with_options("TTCA", -1, 0..4, GeneticCode::STANDARD, true);
    assert!(
        reverse
            .iter()
            .any(|c| c.amino_acid == b'M' && c.positions == [0, 3, 2])
    );
    assert!(frame_codons_with_options("AT", 1, 0..2, GeneticCode::STANDARD, true).is_empty());
}
#[test]
fn orfs_detect_both_strands_circular_starts_limits_and_ambiguous_barriers() {
    let options = OrfOptions {
        min_codons: 2,
        ..OrfOptions::default()
    };
    let doc = Document::new("ORFs", "ATGAAATAATTATTTCAT").unwrap();
    let found = find_orfs(&doc, options);
    assert!(
        found
            .orfs
            .iter()
            .any(|orf| !orf.reverse && orf.start == 0 && orf.protein(&doc.sequence) == "MK*")
    );
    assert!(
        found
            .orfs
            .iter()
            .any(|orf| orf.reverse && orf.start == 17 && orf.protein(&doc.sequence) == "MK*")
    );
    let mut doc = Document::new("origin", "AAATAACATG").unwrap();
    doc.circular = true;
    let found = find_orfs(&doc, options);
    let orf = found
        .orfs
        .iter()
        .find(|orf| !orf.reverse && orf.start == 7)
        .unwrap();
    assert_eq!(orf.protein(&doc.sequence), "MK*");
    assert_eq!(
        orf.segments(10),
        vec![Segment { start: 7, end: 10 }, Segment { start: 0, end: 6 }]
    );
    let doc = Document::new("start", "GTGAAATAA").unwrap();
    assert!(find_orfs(&doc, options).orfs.is_empty());
    let alternative = OrfOptions {
        genetic_code: GeneticCode::BACTERIAL,
        alternative_starts: true,
        ..options
    };
    assert!(
        find_orfs(&doc, alternative)
            .orfs
            .iter()
            .any(|orf| orf.protein(&doc.sequence) == "MK*")
    );
    let doc = Document::new("ambiguous", "ATGTANTAA").unwrap();
    assert!(find_orfs(&doc, options).orfs.is_empty());
    let doc = Document::new("partial", "ATGAAAA").unwrap();
    let result = find_orfs(
        &doc,
        OrfOptions {
            include_partial: true,
            ..options
        },
    );
    assert!(
        result
            .orfs
            .iter()
            .any(|orf| !orf.complete && orf.protein(&doc.sequence) == "MK")
    );
    let limited = find_orfs(
        &Document::new("bounded", "ATGAAATAAATGAAATAA").unwrap(),
        OrfOptions {
            max_results: 1,
            ..options
        },
    );
    assert_eq!(limited.orfs.len(), 1);
    assert!(limited.truncated);
}
#[test]
fn flanks_are_bounded_hypotheses_and_preserve_reverse_and_join_phase() {
    let (doc, mut feature) = cds("TAACCCATGAAATAAGGGTAA", false);
    feature.parts = vec![Segment { start: 6, end: 15 }];
    let flanks = cds_flanking_codons(&doc, &feature, 30).unwrap();
    assert_eq!(
        flanks
            .iter()
            .filter(|c| c.upstream)
            .map(|c| c.codon.amino_acid as char)
            .collect::<String>(),
        "P*"
    );
    assert_eq!(
        flanks
            .iter()
            .filter(|c| !c.upstream)
            .map(|c| c.codon.amino_acid as char)
            .collect::<String>(),
        "G*"
    );
    let mut reversed = doc.clone();
    reversed.features.push(feature);
    reversed.revcomp().unwrap();
    let reverse = cds_flanking_codons(&reversed, &reversed.features[0], 30).unwrap();
    assert_eq!(
        reverse
            .iter()
            .map(|c| c.codon.amino_acid)
            .collect::<Vec<_>>(),
        flanks
            .iter()
            .map(|c| c.codon.amino_acid)
            .collect::<Vec<_>>()
    );
    let mut circular = doc;
    circular.circular = true;
    assert!(
        cds_flanking_codons(&circular, &reversed.features[0], 1)
            .unwrap()
            .len()
            <= 2
    );
}

#[test]
fn fuzzy_boundaries_do_not_force_initiation_or_dual_code_stop() {
    use crate::FeatureLocation;
    let (doc, mut feature) = cds("CTGTGATGA", false);
    qualify(&mut feature, "transl_table", "27");
    feature.location = Some(FeatureLocation::Range {
        start: 0,
        end: 9,
        before: true,
        after: true,
    });
    assert_eq!(translate_cds(&doc, &feature).unwrap().amino_acids, b"LWW");
    feature.location = Some(FeatureLocation::Order(vec![FeatureLocation::Range {
        start: 0,
        end: 9,
        before: false,
        after: false,
    }]));
    assert!(
        translate_cds(&doc, &feature)
            .unwrap_err()
            .0
            .contains("not a resolved local join")
    );
}
#[test]
fn intron_edits_shift_joined_exceptions_without_invalidating_unchanged_codons() {
    let (mut doc, mut feature) = cds("ATGTCCGA", false);
    feature.parts = vec![Segment { start: 0, end: 4 }, Segment { start: 6, end: 8 }];
    qualify(&mut feature, "transl_except", "(pos:join(4,7..8),aa:Sec)");
    doc.features.push(feature);
    doc.replace(5, 5, "A").unwrap();
    assert_eq!(
        translate_cds(&doc, &doc.features[0]).unwrap().amino_acids,
        b"MU"
    );
    assert_eq!(
        doc.features[0].qualifier("transl_except"),
        Some("(pos:join(4,8..9),aa:Sec)")
    );
}
