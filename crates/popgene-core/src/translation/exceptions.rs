//! Exact local INSDC transl_except coordinates, including joined and partial stops.
use super::TranslationError;
use crate::Feature;

pub(super) const REVIEW: &str = "popgene_transl_except_review";
#[derive(Clone, Debug)]
pub(super) struct Exception {
    pub positions: Vec<usize>,
    pub amino_acid: u8,
    pub reverse: bool,
}
fn error(message: &str) -> TranslationError {
    TranslationError(format!("/transl_except {message}"))
}
fn position(text: &str, depth: usize) -> Result<(Vec<usize>, bool), TranslationError> {
    if depth > 8 {
        return Err(error("location nesting exceeds the limit"));
    }
    if let Some(inner) = text
        .strip_prefix("complement(")
        .and_then(|s| s.strip_suffix(')'))
    {
        let (mut p, r) = position(inner, depth + 1)?;
        p.reverse();
        return Ok((p, !r));
    }
    if let Some(inner) = text.strip_prefix("join(").and_then(|s| s.strip_suffix(')')) {
        let mut result = vec![];
        let mut first = 0;
        let mut nesting = 0;
        let mut direction = None;
        for (i, c) in inner
            .char_indices()
            .chain(std::iter::once((inner.len(), ',')))
        {
            if c == '(' {
                nesting += 1;
            } else if c == ')' {
                nesting -= 1;
            } else if c == ',' && nesting == 0 {
                let (p, reverse) = position(&inner[first..i], depth + 1)?;
                if direction.is_some_and(|d| d != reverse) {
                    return Err(error("mixed-strand location"));
                }
                direction = Some(reverse);
                result.extend(p);
                first = i + 1;
                if result.len() > 3 {
                    return Err(error("must identify at most one codon"));
                }
            }
        }
        if result.is_empty() {
            return Err(error("empty location"));
        }
        return Ok((result, direction.unwrap_or(false)));
    }
    let number = |s: &str| s.parse::<usize>().ok().filter(|&n| n > 0);
    let (a, b) = if let Some((a, b)) = text.split_once("..") {
        (number(a), number(b))
    } else {
        (number(text), number(text))
    };
    match (a, b) {
        (Some(a), Some(b)) if b >= a && b - a < 3 => Ok(((a - 1..b).collect(), false)),
        _ => Err(error("requires exact local codon coordinates")),
    }
}
pub(super) fn parse(text: &str) -> Result<Exception, TranslationError> {
    if text.len() > 1024 {
        return Err(error("exceeds the value limit"));
    }
    let text: String = text.chars().filter(|c| !c.is_ascii_whitespace()).collect();
    let inner = text
        .strip_prefix("(pos:")
        .and_then(|s| s.strip_suffix(')'))
        .ok_or_else(|| error("must use (pos:location,aa:amino-acid)"))?;
    let (location, aa) = inner
        .rsplit_once(",aa:")
        .ok_or_else(|| error("has no amino-acid value"))?;
    let amino_acid = match aa.to_ascii_uppercase().as_str() {
        "ALA" => b'A',
        "ARG" => b'R',
        "ASN" => b'N',
        "ASP" => b'D',
        "CYS" => b'C',
        "GLN" => b'Q',
        "GLU" => b'E',
        "GLY" => b'G',
        "HIS" => b'H',
        "ILE" => b'I',
        "LEU" => b'L',
        "LYS" => b'K',
        "MET" => b'M',
        "PHE" => b'F',
        "PRO" => b'P',
        "SER" => b'S',
        "THR" => b'T',
        "TRP" => b'W',
        "TYR" => b'Y',
        "VAL" => b'V',
        "ASX" => b'B',
        "GLX" => b'Z',
        "XLE" => b'J',
        "SEC" => b'U',
        "PYL" => b'O',
        "TERM" => b'*',
        "OTHER" | "XAA" => b'X',
        _ => return Err(error("has an unknown amino-acid abbreviation")),
    };
    let (positions, reverse) = position(location, 0)?;
    if positions.len() != 3 && amino_acid != b'*' {
        return Err(error("only TERM can describe a partial codon"));
    }
    Ok(Exception {
        positions,
        amino_acid,
        reverse,
    })
}
fn format(exception: &Exception) -> String {
    let mut positions = exception.positions.clone();
    if exception.reverse {
        positions.reverse();
    }
    let mut spans: Vec<(usize, usize)> = vec![];
    for p in positions {
        if let Some(last) = spans.last_mut()
            && last.1 == p
        {
            last.1 = p + 1;
            continue;
        }
        spans.push((p, p + 1));
    }
    let mut location = spans
        .iter()
        .map(|(a, b)| {
            if a + 1 == *b {
                (a + 1).to_string()
            } else {
                format!("{}..{}", a + 1, b)
            }
        })
        .collect::<Vec<_>>()
        .join(",");
    if spans.len() > 1 {
        location = format!("join({location})");
    }
    if exception.reverse {
        location = format!("complement({location})");
    }
    let aa = match exception.amino_acid {
        b'A' => "Ala",
        b'R' => "Arg",
        b'N' => "Asn",
        b'D' => "Asp",
        b'C' => "Cys",
        b'Q' => "Gln",
        b'E' => "Glu",
        b'G' => "Gly",
        b'H' => "His",
        b'I' => "Ile",
        b'L' => "Leu",
        b'K' => "Lys",
        b'M' => "Met",
        b'F' => "Phe",
        b'P' => "Pro",
        b'S' => "Ser",
        b'T' => "Thr",
        b'W' => "Trp",
        b'Y' => "Tyr",
        b'V' => "Val",
        b'B' => "Asx",
        b'Z' => "Glx",
        b'J' => "Xle",
        b'U' => "Sec",
        b'O' => "Pyl",
        b'*' => "TERM",
        _ => "OTHER",
    };
    format!("(pos:{location},aa:{aa})")
}
/// Coordinates touched by an edit cannot silently retain their old recoding claim.
pub fn remap_translation_exceptions(
    feature: &mut Feature,
    start: usize,
    end: usize,
    replacement_len: usize,
) {
    let mut review = false;
    for (key, value) in &mut feature.qualifiers {
        if key != "transl_except" {
            continue;
        }
        let Some(text) = value else {
            continue;
        };
        let Ok(mut exception) = parse(text) else {
            continue;
        };
        let overlaps = if start == end {
            feature
                .parts
                .iter()
                .any(|part| part.start < start && start < part.end)
                && exception
                    .positions
                    .windows(2)
                    .any(|pair| pair[0].abs_diff(pair[1]) == 1 && pair[0].max(pair[1]) == start)
        } else {
            exception
                .positions
                .iter()
                .any(|position| (start..end).contains(position))
        };
        if overlaps {
            review = true;
            continue;
        }
        for p in &mut exception.positions {
            if *p >= end {
                *p = (*p as isize + replacement_len as isize - (end - start) as isize) as usize;
            }
        }
        *text = format(&exception);
    }
    if review && !feature.qualifiers.iter().any(|(key, _)| key == REVIEW) {
        feature.qualifiers.push((REVIEW.into(),Some("Sequence edit overlaps a recorded translation exception; review its coordinates and residue.".into())));
    }
}
pub fn reverse_translation_exceptions(feature: &mut Feature, length: usize) {
    for (key, value) in &mut feature.qualifiers {
        if key != "transl_except" {
            continue;
        }
        let Some(text) = value else {
            continue;
        };
        let Ok(mut exception) = parse(text) else {
            continue;
        };
        if exception.positions.iter().any(|p| *p >= length) {
            continue;
        }
        for p in &mut exception.positions {
            *p = length - 1 - *p;
        }
        exception.reverse = !exception.reverse;
        *text = format(&exception);
    }
}

/// Project exact exception coordinates into a fragment or assembled molecule.
/// Call after setting the output feature strand. An unmappable partial codon is
/// retained only as review provenance, never as an active false-coordinate exception.
pub fn project_translation_exceptions(feature: &mut Feature, map: impl Fn(usize) -> Option<usize>) {
    let mut qualifiers = Vec::with_capacity(feature.qualifiers.len());
    let mut review = false;
    for (key, value) in std::mem::take(&mut feature.qualifiers) {
        if key != "transl_except" {
            qualifiers.push((key, value));
            continue;
        }
        let parsed = value.as_deref().map(parse).transpose();
        match parsed {
            Ok(Some(mut exception)) => {
                let mapped = exception
                    .positions
                    .iter()
                    .map(|&p| map(p))
                    .collect::<Vec<_>>();
                if mapped.iter().all(Option::is_none) {
                    continue;
                }
                if mapped.iter().all(Option::is_some) {
                    exception.positions = mapped.into_iter().map(Option::unwrap).collect();
                    exception.reverse = feature.reverse;
                    qualifiers.push((key, Some(format(&exception))));
                } else {
                    review = true;
                    qualifiers.push(("popgene_transl_except_original".into(), value));
                }
            }
            _ => {
                review = true;
                qualifiers.push(("popgene_transl_except_original".into(), value));
            }
        }
    }
    if review && !qualifiers.iter().any(|(k, _)| k == REVIEW) {
        qualifiers.push((REVIEW.into(),Some("A fragment or assembly projection could not retain all translation-exception positions; review the preserved original value.".into())));
    }
    feature.qualifiers = qualifiers;
}

#[cfg(test)]
mod projection_tests {
    use super::*;
    #[test]
    fn projection_keeps_wrapped_and_reversed_codons_and_retires_partial_values() {
        let mut feature = Feature {
            id: 1,
            label: "CDS".into(),
            kind: "CDS".into(),
            parts: vec![],
            location: None,
            reverse: false,
            qualifiers: vec![("transl_except".into(), Some("(pos:9..11,aa:Sec)".into()))],
            edited: false,
        };
        project_translation_exceptions(&mut feature, |p| Some((p + 2) % 10));
        assert_eq!(
            parse(feature.qualifier("transl_except").unwrap())
                .unwrap()
                .positions,
            vec![0, 1, 2]
        );
        feature.reverse = true;
        project_translation_exceptions(&mut feature, |p| Some(9 - p));
        let exception = parse(feature.qualifier("transl_except").unwrap()).unwrap();
        assert_eq!(exception.positions, vec![9, 8, 7]);
        assert!(exception.reverse);
        project_translation_exceptions(&mut feature, |p| (p < 9).then_some(p));
        assert!(feature.qualifier("transl_except").is_none());
        assert!(feature.qualifier(REVIEW).is_some());
        assert!(
            feature
                .qualifier("popgene_transl_except_original")
                .is_some()
        );
        feature.qualifiers = vec![("transl_except".into(), Some("(pos:1..3,aa:Sec)".into()))];
        project_translation_exceptions(&mut feature, |_| None);
        assert!(feature.qualifiers.is_empty());
    }
}
