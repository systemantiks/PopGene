//! Native rust-bio 4.2.2 alignments (MIT), with bounded work and exact residue maps.
//! Pairwise uses affine gaps: open + (length - 1) × extend. Multiple rows are
//! reference-guided center-star alignments; POA consensus is a separate analysis
//! with rust-bio's documented linear gap penalty.
//! https://docs.rs/bio/4.2.2/bio/alignment/index.html
use crate::{Document, Molecule, Segment};
use bio::alignment::{AlignmentOperation, pairwise, poa};
use serde::{Deserialize, Serialize};

const MAX_PAIR_CELLS: usize = 25_000_000;
const MAX_MULTIPLE_CELLS: usize = 100_000_000;
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum AlignmentMode {
    Global,
    Local,
    QueryInReference,
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct AlignmentOptions {
    pub mode: AlignmentMode,
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub both_strands: bool,
    pub protein_blosum62: bool,
}
impl Default for AlignmentOptions {
    fn default() -> Self {
        Self {
            mode: AlignmentMode::QueryInReference,
            match_score: 2,
            mismatch_score: -3,
            gap_open: -5,
            gap_extend: -1,
            both_strands: true,
            protein_blosum62: true,
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[error("{0}")]
pub struct AlignmentError(pub String);
fn error(message: &str) -> AlignmentError {
    AlignmentError(message.into())
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AlignedRow {
    pub name: String,
    pub sequence: String,
    /// Original zero-based residue coordinates; gaps have no invented coordinate.
    pub positions: Vec<Option<usize>>,
    pub qualities: Vec<Option<u8>>,
    pub reverse: bool,
}
#[derive(Clone, Debug)]
pub struct PairAlignment {
    pub molecule: Molecule,
    pub reference: AlignedRow,
    pub query: AlignedRow,
    pub reference_span: Segment,
    /// True when the unrolled reference span crosses its stored circular origin.
    pub crosses_origin: bool,
    pub query_span: Segment,
    pub score: i32,
    pub exact_matches: usize,
    pub aligned_residues: usize,
    pub gaps: usize,
    pub options: AlignmentOptions,
}
impl PairAlignment {
    pub fn identity_percent(&self) -> f64 {
        if self.aligned_residues == 0 {
            0.
        } else {
            100. * self.exact_matches as f64 / self.aligned_residues as f64
        }
    }
    pub fn fasta(&self) -> String {
        alignment_fasta(&[self.reference.clone(), self.query.clone()])
    }
    pub fn pretty(&self, width: usize) -> String {
        pretty_rows(&[self.reference.clone(), self.query.clone()], width)
    }
    pub fn differences(&self) -> String {
        let mut lines =
            vec!["Reference base\tQuery base\tReference\tQuery\tQuery quality".to_string()];
        for (index, (a, b)) in self
            .reference
            .sequence
            .bytes()
            .zip(self.query.sequence.bytes())
            .enumerate()
        {
            if a == b
                || (self.molecule != Molecule::Protein
                    && matches!((a, b), (b'T', b'U') | (b'U', b'T')))
            {
                continue;
            }
            let position = |p: Option<usize>| p.map_or("—".into(), |p| (p + 1).to_string());
            lines.push(format!(
                "{}\t{}\t{}\t{}\t{}",
                position(self.reference.positions[index]),
                position(self.query.positions[index]),
                a as char,
                b as char,
                self.query.qualities[index].map_or("—".into(), |q| q.to_string())
            ));
        }
        lines.join("\n")
    }
}
#[derive(Clone, Debug)]
pub struct MultipleAlignment {
    pub rows: Vec<AlignedRow>,
    /// Majority/IUPAC consensus on the displayed reference-guided columns.
    pub consensus: String,
    pub pair_scores: Vec<i32>,
    pub options: AlignmentOptions,
}
impl MultipleAlignment {
    pub fn fasta(&self) -> String {
        alignment_fasta(&self.rows)
    }
    pub fn pretty(&self, width: usize) -> String {
        pretty_rows(&self.rows, width)
    }
}
fn validate(
    reference: &Document,
    query: &Document,
    options: AlignmentOptions,
) -> Result<(), AlignmentError> {
    if reference.sequence.is_empty() || query.sequence.is_empty() {
        return Err(error("Choose two non-empty sequences"));
    }
    if (reference.molecule == Molecule::Protein) != (query.molecule == Molecule::Protein) {
        return Err(error(
            "Protein and nucleotide sequences cannot share a residue alignment",
        ));
    }
    if options.gap_open >= 0
        || options.gap_extend > 0
        || options.match_score <= 0
        || options.mismatch_score >= 0
        || [
            options.gap_open,
            options.gap_extend,
            options.match_score,
            options.mismatch_score,
        ]
        .iter()
        .any(|v| v.unsigned_abs() > 100)
    {
        return Err(error(
            "Use match scores 1–100, mismatch/open scores −100…−1 and extension scores −100…0",
        ));
    }
    if (reference.sequence.len() + 1)
        .checked_mul(query.sequence.len() + 1)
        .is_none_or(|cells| cells > MAX_PAIR_CELLS)
    {
        return Err(error(
            "Alignment exceeds 25 million dynamic-programming cells; select a smaller reference region",
        ));
    }
    reference
        .molecule
        .normalize(&reference.sequence)
        .map_err(|e| error(&e.to_string()))?;
    query
        .molecule
        .normalize(&query.sequence)
        .map_err(|e| error(&e.to_string()))?;
    Ok(())
}
fn mask(base: u8) -> u8 {
    match base {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' | b'U' => 8,
        b'R' => 5,
        b'Y' => 10,
        b'S' => 6,
        b'W' => 9,
        b'K' => 12,
        b'M' => 3,
        b'B' => 14,
        b'D' => 13,
        b'H' => 11,
        b'V' => 7,
        b'N' => 15,
        _ => 0,
    }
}
fn substitution(a: u8, b: u8, protein: bool, options: AlignmentOptions) -> i32 {
    if protein && options.protein_blosum62 {
        let canonical = |b| {
            if b"ACDEFGHIKLMNPQRSTVWYBZX*".contains(&b) {
                b
            } else {
                b'X'
            }
        };
        return bio::scores::blosum62(canonical(a), canonical(b));
    }
    if a == b {
        options.match_score
    } else if !protein && mask(a) & mask(b) != 0 {
        options.match_score / 2
    } else {
        options.mismatch_score
    }
}
fn quality(document: &Document, position: usize) -> Option<u8> {
    document
        .trace
        .as_ref()
        .filter(|trace| trace.qualities_known)
        .and_then(|trace| trace.qualities.get(position))
        .copied()
        .or_else(|| {
            document
                .sequence_qualities
                .as_ref()
                .and_then(|values| values.get(position))
                .copied()
                .flatten()
        })
}

fn row(
    document: &Document,
    sequence: String,
    positions: Vec<Option<usize>>,
    reverse: bool,
) -> AlignedRow {
    let qualities = positions
        .iter()
        .map(|p| p.and_then(|p| quality(document, p)))
        .collect();
    AlignedRow {
        name: document.name.clone(),
        sequence,
        positions,
        qualities,
        reverse,
    }
}
pub fn align_pair(
    reference: &Document,
    query: &Document,
    options: AlignmentOptions,
) -> Result<PairAlignment, AlignmentError> {
    if reference.circular && options.mode != AlignmentMode::Global {
        return align_circular(reference, query, options);
    }
    align_linear(reference, query, options)
}
fn align_circular(
    reference: &Document,
    query: &Document,
    options: AlignmentOptions,
) -> Result<PairAlignment, AlignmentError> {
    validate(reference, query, options)?;
    let n = reference.sequence.len();
    let project = |mut result: PairAlignment, origin: usize| {
        let length = result.reference_span.end - result.reference_span.start;
        let start = (result.reference_span.start + origin) % n;
        result.reference_span = Segment {
            start,
            end: start + length,
        };
        result.crosses_origin = start + length > n;
        for position in result.reference.positions.iter_mut().flatten() {
            *position = (*position + origin) % n;
        }
        result.reference.qualities = result
            .reference
            .positions
            .iter()
            .map(|p| p.and_then(|p| quality(reference, p)))
            .collect();
        result
    };
    // A valid optimum on the doubled reference is also optimal among all paths
    // that consume at most one circumference. Starts in the second copy collapse
    // back to their actual original coordinates.
    let mut doubled = reference.clone();
    doubled.circular = false;
    doubled.sequence = reference.sequence.repeat(2).into();
    let result = align_linear(&doubled, query, options)?;
    if result.reference_span.end - result.reference_span.start <= n {
        return Ok(project(result, 0));
    }
    // If the unrestricted optimum crosses more than one lap, search each single
    // circumference exactly. Work is capped rather than accepting a duplicated
    // multi-lap match or silently returning a lower, unproven heuristic result.
    let cells = n
        .saturating_mul(n + 1)
        .saturating_mul(query.sequence.len() + 1)
        .saturating_mul(if options.both_strands { 2 } else { 1 });
    if cells > MAX_MULTIPLE_CELLS {
        return Err(error(
            "The best doubled-reference path exceeds one circular lap, and exact one-lap fallback exceeds 100 million cells; use a shorter read or reference region",
        ));
    }
    let mut best = None;
    for origin in 0..n {
        doubled.sequence = format!(
            "{}{}",
            &reference.sequence[origin..],
            &reference.sequence[..origin]
        )
        .into();
        let candidate = project(align_linear(&doubled, query, options)?, origin);
        if best
            .as_ref()
            .is_none_or(|b: &PairAlignment| candidate.score > b.score)
        {
            best = Some(candidate);
        }
    }
    best.ok_or_else(|| error("No one-lap circular alignment"))
}
fn align_linear(
    reference: &Document,
    query: &Document,
    options: AlignmentOptions,
) -> Result<PairAlignment, AlignmentError> {
    validate(reference, query, options)?;
    let protein = reference.molecule == Molecule::Protein;
    let ref_bases = if protein {
        reference.sequence.to_string()
    } else {
        reference.sequence.replace('U', "T")
    };
    let run = |query_bases: &[u8]| {
        let mut aligner = pairwise::Aligner::with_capacity(
            query_bases.len(),
            ref_bases.len(),
            options.gap_open,
            options.gap_extend,
            |a, b| substitution(a, b, protein, options),
        );
        match options.mode {
            AlignmentMode::Global => aligner.global(query_bases, ref_bases.as_bytes()),
            AlignmentMode::Local => aligner.local(query_bases, ref_bases.as_bytes()),
            AlignmentMode::QueryInReference => {
                aligner.semiglobal(query_bases, ref_bases.as_bytes())
            }
        }
    };
    let forward = if protein {
        query.sequence.to_string()
    } else {
        query.sequence.replace('U', "T")
    };
    let mut aligned = run(forward.as_bytes());
    let mut reverse = false;
    let mut query_bases = forward;
    if options.both_strands && !protein {
        let reversed = crate::reverse_complement(&query_bases);
        let candidate = run(reversed.as_bytes());
        if candidate.score > aligned.score {
            aligned = candidate;
            reverse = true;
            query_bases = reversed;
        }
    }
    let mut r = aligned.ystart;
    let mut q = aligned.xstart;
    let mut reference_text = String::new();
    let mut query_text = String::new();
    let mut rp = vec![];
    let mut qp = vec![];
    let mut matches = 0;
    let mut residues = 0;
    let mut gaps = 0;
    let original_query_position = |p| {
        if reverse {
            query.sequence.len() - 1 - p
        } else {
            p
        }
    };
    for operation in &aligned.operations {
        match operation {
            AlignmentOperation::Match | AlignmentOperation::Subst => {
                let a = reference.sequence.as_bytes()[r];
                let mut b = query_bases.as_bytes()[q];
                if query.molecule == Molecule::Rna && b == b'T' {
                    b = b'U';
                }
                if (a == b || (!protein && mask(a) == mask(b)))
                    && (protein && a != b'X' || !protein && mask(a).count_ones() == 1)
                {
                    matches += 1;
                }
                residues += 1;
                reference_text.push(a as char);
                query_text.push(b as char);
                rp.push(Some(r));
                qp.push(Some(original_query_position(q)));
                r += 1;
                q += 1;
            }
            AlignmentOperation::Ins => {
                reference_text.push('-');
                let mut b = query_bases.as_bytes()[q];
                if query.molecule == Molecule::Rna && b == b'T' {
                    b = b'U';
                }
                query_text.push(b as char);
                rp.push(None);
                qp.push(Some(original_query_position(q)));
                q += 1;
                gaps += 1;
            }
            AlignmentOperation::Del => {
                reference_text.push(reference.sequence.as_bytes()[r] as char);
                query_text.push('-');
                rp.push(Some(r));
                qp.push(None);
                r += 1;
                gaps += 1;
            }
            AlignmentOperation::Xclip(_) | AlignmentOperation::Yclip(_) => {
                return Err(error(
                    "Unexpected clipping operation in a standard alignment",
                ));
            }
        }
    }
    let query_span = if reverse {
        Segment {
            start: query.sequence.len() - aligned.xend,
            end: query.sequence.len() - aligned.xstart,
        }
    } else {
        Segment {
            start: aligned.xstart,
            end: aligned.xend,
        }
    };
    Ok(PairAlignment {
        molecule: reference.molecule,
        reference: row(reference, reference_text, rp, false),
        query: row(query, query_text, qp, reverse),
        reference_span: Segment {
            start: aligned.ystart,
            end: aligned.yend,
        },
        crosses_origin: false,
        query_span,
        score: aligned.score,
        exact_matches: matches,
        aligned_residues: residues,
        gaps,
        options,
    })
}
/// Merge independently optimal query/reference alignments on shared reference
/// columns. Insertions at a shared boundary are left aligned and padded, making
/// the reference-guided heuristic explicit rather than claiming a global MSA optimum.
pub fn align_multiple(
    documents: &[Document],
    options: AlignmentOptions,
) -> Result<MultipleAlignment, AlignmentError> {
    if !(2..=32).contains(&documents.len()) {
        return Err(error("Choose 2–32 sequences for multiple alignment"));
    }
    if options.mode == AlignmentMode::Local {
        return Err(error(
            "Multiple alignment needs Global or Query in reference mode to preserve complete input sequences",
        ));
    }
    let reference = &documents[0];
    let n = reference.sequence.len();
    // Global and query-in-reference MSA retain every residue. Even with no
    // inserted columns, the longest input is a hard lower bound on output.
    let longest = documents
        .iter()
        .map(|doc| doc.sequence.len())
        .max()
        .unwrap_or(0);
    if longest.saturating_mul(documents.len()) > 5_000_000 {
        return Err(error("Gapped output exceeds five million displayed cells"));
    }
    let mut budget = 0_usize;
    for query in &documents[1..] {
        validate(reference, query, options)?;
        budget = budget.saturating_add((n + 1).saturating_mul(query.sequence.len() + 1));
    }
    if budget > MAX_MULTIPLE_CELLS {
        return Err(error(
            "Multiple alignment exceeds 100 million cells; use fewer/shorter sequences",
        ));
    }
    struct Projection {
        bases: Vec<(u8, Option<usize>)>,
        insertions: Vec<Vec<(u8, Option<usize>)>>,
        reverse: bool,
    }
    let mut projections = vec![];
    let mut widths = vec![0; n + 1];
    let mut columns = n;
    let mut scores = vec![];
    for query in &documents[1..] {
        // A multiple alignment uses one shared linear column ordering. Unlike
        // pairwise circular alignment, all rows retain the stored origin here.
        let aligned = align_linear(reference, query, options)?;
        scores.push(aligned.score);
        let mut projection = Projection {
            bases: vec![(b'-', None); n],
            insertions: vec![vec![]; n + 1],
            reverse: aligned.query.reverse,
        };
        let mut boundary = aligned.reference_span.start;
        for (index, (a, b)) in aligned
            .reference
            .sequence
            .bytes()
            .zip(aligned.query.sequence.bytes())
            .enumerate()
        {
            if a == b'-' {
                let length = projection.insertions[boundary].len() + 1;
                if length > widths[boundary] {
                    columns += length - widths[boundary];
                    if columns.saturating_mul(documents.len()) > 5_000_000 {
                        return Err(error("Gapped output exceeds five million displayed cells"));
                    }
                    widths[boundary] = length;
                }
                projection.insertions[boundary].push((b, aligned.query.positions[index]));
            } else {
                projection.bases[boundary] = (b, aligned.query.positions[index]);
                boundary += 1;
            }
        }
        projections.push(projection);
    }
    debug_assert_eq!(columns, n + widths.iter().sum::<usize>());
    let mut rows = vec![];
    for (index, document) in documents.iter().enumerate() {
        let mut text = String::with_capacity(columns);
        let mut positions = Vec::with_capacity(columns);
        for (boundary, &width) in widths.iter().enumerate() {
            for offset in 0..width {
                let (base, position) = if index == 0 {
                    (b'-', None)
                } else {
                    projections[index - 1].insertions[boundary]
                        .get(offset)
                        .copied()
                        .unwrap_or((b'-', None))
                };
                text.push(base as char);
                positions.push(position);
            }
            if boundary < n {
                let (base, position) = if index == 0 {
                    (reference.sequence.as_bytes()[boundary], Some(boundary))
                } else {
                    projections[index - 1].bases[boundary]
                };
                text.push(base as char);
                positions.push(position);
            }
        }
        rows.push(row(
            document,
            text,
            positions,
            index > 0 && projections[index - 1].reverse,
        ));
    }
    let consensus = column_consensus(&rows, reference.molecule);
    Ok(MultipleAlignment {
        rows,
        consensus,
        pair_scores: scores,
        options,
    })
}
fn column_consensus(rows: &[AlignedRow], molecule: Molecule) -> String {
    let n = rows.first().map_or(0, |row| row.sequence.len());
    let mut result = String::new();
    for column in 0..n {
        let mut counts = [0_usize; 256];
        for row in rows {
            let mut b = row.sequence.as_bytes()[column];
            if molecule != Molecule::Protein && b == b'U' {
                b = b'T';
            }
            if b != b'-' {
                counts[b as usize] += 1;
            }
        }
        let maximum = *counts.iter().max().unwrap();
        if maximum == 0 {
            continue;
        }
        let winners: Vec<_> = counts
            .iter()
            .enumerate()
            .filter_map(|(b, &count)| (count == maximum).then_some(b as u8))
            .collect();
        let base = if winners.len() == 1 {
            winners[0]
        } else if molecule == Molecule::Protein {
            b'X'
        } else {
            let mask = winners.iter().fold(0, |bits, &b| bits | mask(b));
            b"NACMGRSVTWYHKDBN"[mask as usize]
        };
        result.push(if molecule == Molecule::Rna && base == b'T' {
            'U'
        } else {
            base as char
        });
    }
    result
}
/// Progressive partial-order consensus. POA uses a linear gap penalty (open per
/// gap column), not the pairwise affine penalty. Inputs are oriented beforehand.
pub fn poa_consensus(
    documents: &[Document],
    options: AlignmentOptions,
) -> Result<String, AlignmentError> {
    if !(2..=32).contains(&documents.len()) {
        return Err(error("Choose 2–32 sequences for POA consensus"));
    }
    let reference = &documents[0];
    let protein = reference.molecule == Molecule::Protein;
    let mut sequences = vec![if protein {
        reference.sequence.to_string()
    } else {
        reference.sequence.replace('U', "T")
    }];
    let mut budget = 0_usize;
    let mut total = reference.sequence.len();
    for query in &documents[1..] {
        validate(reference, query, options)?;
        budget = budget.saturating_add(total.saturating_mul(query.sequence.len()));
        total += query.sequence.len();
        if budget > MAX_MULTIPLE_CELLS {
            return Err(error("POA graph work exceeds 100 million cells"));
        }
        let aligned = align_pair(reference, query, options)?;
        let sequence = if protein {
            query.sequence.to_string()
        } else {
            query.sequence.replace('U', "T")
        };
        sequences.push(if aligned.query.reverse {
            crate::reverse_complement(&sequence)
        } else {
            sequence
        });
    }
    let scoring = pairwise::Scoring::new(options.gap_open, 0, |a, b| {
        substitution(a, b, protein, options)
    });
    let mut aligner = poa::Aligner::new(scoring, sequences[0].as_bytes());
    for sequence in &sequences[1..] {
        aligner.global(sequence.as_bytes()).add_to_graph();
    }
    let consensus = String::from_utf8(aligner.consensus())
        .map_err(|_| error("POA returned invalid residue encoding"))?;
    Ok(if reference.molecule == Molecule::Rna {
        consensus.replace('T', "U")
    } else {
        consensus
    })
}
pub fn alignment_fasta(rows: &[AlignedRow]) -> String {
    let mut text = String::new();
    for row in rows {
        text.push_str(&format!(
            ">{}{}\n",
            row.name.replace(['\n', '\r'], " "),
            if row.reverse {
                " [reverse complement]"
            } else {
                ""
            }
        ));
        for chunk in row.sequence.as_bytes().chunks(80) {
            text.push_str(std::str::from_utf8(chunk).unwrap());
            text.push('\n');
        }
    }
    text
}
pub fn pretty_rows(rows: &[AlignedRow], width: usize) -> String {
    let width = width.clamp(20, 160);
    let length = rows.first().map_or(0, |r| r.sequence.len());
    let mut text = String::new();
    for start in (0..length).step_by(width) {
        let end = (start + width).min(length);
        for row in rows {
            let first = row.positions[start..end]
                .iter()
                .flatten()
                .next()
                .map_or("—".into(), |p| (p + 1).to_string());
            let last = row.positions[start..end]
                .iter()
                .flatten()
                .next_back()
                .map_or("—".into(), |p| (p + 1).to_string());
            let name: String = row.name.chars().take(18).collect();
            text.push_str(&format!(
                "{name:<18} {first:>7}  {}  {last}\n",
                &row.sequence[start..end]
            ));
        }
        text.push('\n');
    }
    text
}

#[cfg(test)]
mod tests {
    use super::*;
    fn doc(name: &str, bases: &str) -> Document {
        Document::new(name, bases).unwrap()
    }
    #[test]
    fn affine_gap_scores_and_original_coordinates() {
        let reference = doc("reference", "ACGTACGT");
        let query = doc("query", "ACGTTACGT");
        let result = align_pair(
            &reference,
            &query,
            AlignmentOptions {
                mode: AlignmentMode::Global,
                both_strands: false,
                ..AlignmentOptions::default()
            },
        )
        .unwrap();
        assert_eq!(result.score, 11);
        assert_eq!(result.gaps, 1);
        assert_eq!(
            result.reference.sequence.replace('-', ""),
            reference.sequence.as_ref()
        );
        assert_eq!(
            result.query.sequence.replace('-', ""),
            query.sequence.as_ref()
        );
        let short = align_pair(
            &doc("r", "AAAAAA"),
            &doc("q", "AAA"),
            AlignmentOptions {
                mode: AlignmentMode::Global,
                both_strands: false,
                ..AlignmentOptions::default()
            },
        )
        .unwrap();
        assert_eq!(short.score, 6 - 5 - 2); // gap_open + (3-1)*gap_extend
    }
    #[test]
    fn reverse_trace_coordinates_quality_and_semiglobal_clips() {
        let reference = doc("reference", "GGGATGCCCTTT");
        let mut query = doc("read", "GGGCAT");
        query.sequence_qualities = Some(vec![
            Some(10),
            Some(20),
            Some(30),
            Some(40),
            Some(50),
            Some(60),
        ]);
        let result = align_pair(&reference, &query, AlignmentOptions::default()).unwrap();
        assert!(result.query.reverse);
        assert_eq!(result.reference_span, Segment { start: 3, end: 9 });
        assert_eq!(result.query.sequence, "ATGCCC");
        assert_eq!(
            result.query.positions,
            vec![Some(5), Some(4), Some(3), Some(2), Some(1), Some(0)]
        );
        assert_eq!(result.query.qualities[0], Some(60));
        assert_eq!(result.identity_percent(), 100.);
    }
    #[test]
    fn protein_local_alignment_uses_blosum62() {
        let reference = Document::new_typed("protein", "PEEKSAV", Molecule::Protein).unwrap();
        let query = Document::new_typed("query", "LSPADKTNVKAA", Molecule::Protein).unwrap();
        let result = align_pair(
            &reference,
            &query,
            AlignmentOptions {
                mode: AlignmentMode::Local,
                gap_open: -10,
                gap_extend: -1,
                ..AlignmentOptions::default()
            },
        )
        .unwrap();
        assert_eq!(result.score, 16);
        assert_eq!(result.reference_span, Segment { start: 0, end: 7 });
        assert_eq!(result.query_span, Segment { start: 2, end: 9 });
        assert!(!result.query.reverse);
    }
    #[test]
    fn reference_guided_msa_preserves_all_bases_and_column_maps() {
        let documents = vec![
            doc("reference", "ACGTACGT"),
            doc("insertion", "ACGTTACGT"),
            doc("deletion", "ACGACGT"),
        ];
        let result = align_multiple(
            &documents,
            AlignmentOptions {
                mode: AlignmentMode::Global,
                both_strands: false,
                ..AlignmentOptions::default()
            },
        )
        .unwrap();
        let columns = result.rows[0].sequence.len();
        assert_eq!(columns, 9);
        for (row, source) in result.rows.iter().zip(&documents) {
            assert_eq!(row.sequence.len(), columns);
            assert_eq!(row.sequence.replace('-', ""), source.sequence.as_ref());
            assert_eq!(
                row.positions.iter().flatten().copied().collect::<Vec<_>>(),
                (0..source.sequence.len()).collect::<Vec<_>>()
            );
        }
        assert_eq!(result.consensus.len(), columns);
        assert!(result.fasta().contains(">insertion"));
        let consensus = poa_consensus(
            &[doc("a", "ACGT"), doc("b", "ACGT"), doc("c", "ACCT")],
            AlignmentOptions {
                mode: AlignmentMode::Global,
                both_strands: false,
                ..AlignmentOptions::default()
            },
        )
        .unwrap();
        assert_eq!(consensus, "ACGT");
    }
    #[test]
    fn rejects_unbounded_and_mixed_molecule_work() {
        let protein = Document::new_typed("p", "ACGT", Molecule::Protein).unwrap();
        assert!(align_pair(&doc("d", "ACGT"), &protein, AlignmentOptions::default()).is_err());
        assert!(
            align_pair(
                &doc("a", &"A".repeat(6000)),
                &doc("b", &"A".repeat(6000)),
                AlignmentOptions::default()
            )
            .is_err()
        );
        assert!(
            align_multiple(
                &[doc("a", "ACG"), doc("b", "ACG")],
                AlignmentOptions {
                    mode: AlignmentMode::Local,
                    ..AlignmentOptions::default()
                }
            )
            .is_err()
        );
    }
    #[test]
    fn circular_origin_search_maps_both_strands_and_original_quality_coordinates() {
        let mut reference = doc("circle", "AAAACCCC");
        reference.circular = true;
        reference.sequence_qualities = Some((10..18).map(Some).collect());
        for (sequence, reverse) in [("CCAAAA", false), ("TTTTGG", true)] {
            let result = align_pair(
                &reference,
                &doc("read", sequence),
                AlignmentOptions::default(),
            )
            .unwrap();
            assert_eq!(result.score, 12);
            assert!(result.crosses_origin);
            assert_eq!(result.query.reverse, reverse);
            assert_eq!(
                result.reference.positions,
                vec![Some(6), Some(7), Some(0), Some(1), Some(2), Some(3)]
            );
            assert_eq!(
                result.reference.qualities,
                vec![Some(16), Some(17), Some(10), Some(11), Some(12), Some(13)]
            );
            assert_eq!(result.reference_span, Segment { start: 6, end: 12 });
        }
    }
    #[test]
    fn circular_alignment_never_accepts_a_multi_lap_repeat_and_matches_all_rotations() {
        let mut reference = doc("repeat", "ACGT");
        reference.circular = true;
        let result = align_pair(
            &reference,
            &doc("read", "ACGTACGTACGT"),
            AlignmentOptions::default(),
        )
        .unwrap();
        assert_eq!(result.reference_span.end - result.reference_span.start, 4);
        assert_eq!(
            result
                .reference
                .positions
                .iter()
                .flatten()
                .copied()
                .collect::<std::collections::HashSet<_>>()
                .len(),
            4
        );
        for reference_sequence in ["ACGTTAC", "GGAACT", "ATATTA"] {
            let mut reference = doc("circle", reference_sequence);
            reference.circular = true;
            for query_sequence in ["AACGT", "TACT", "ATATTAAT", "CGGAA"] {
                let query = doc("read", query_sequence);
                for mode in [AlignmentMode::Local, AlignmentMode::QueryInReference] {
                    let options = AlignmentOptions {
                        mode,
                        ..AlignmentOptions::default()
                    };
                    let expected = (0..reference_sequence.len())
                        .map(|origin| {
                            let rotated = doc(
                                "rotation",
                                &format!(
                                    "{}{}",
                                    &reference_sequence[origin..],
                                    &reference_sequence[..origin]
                                ),
                            );
                            align_linear(&rotated, &query, options).unwrap().score
                        })
                        .max()
                        .unwrap();
                    assert_eq!(
                        align_pair(&reference, &query, options).unwrap().score,
                        expected
                    );
                }
            }
        }
    }
    #[test]
    fn mixed_dna_rna_consensus_uses_reference_alphabet_without_changing_rows() {
        let options = AlignmentOptions {
            mode: AlignmentMode::Global,
            both_strands: false,
            ..AlignmentOptions::default()
        };
        for molecule in [Molecule::Dna, Molecule::Rna] {
            let reference = Document::new_typed("reference", "AA", molecule).unwrap();
            let inputs = [
                reference,
                Document::new_typed("rna", "AAU", Molecule::Rna).unwrap(),
                Document::new("dna", "AAT").unwrap(),
            ];
            let alignment = align_multiple(&inputs, options).unwrap();
            assert_eq!(
                alignment.consensus,
                if molecule == Molecule::Rna {
                    "AAU"
                } else {
                    "AAT"
                }
            );
            assert_eq!(alignment.rows[1].sequence.replace('-', ""), "AAU");
            assert!(Document::new_typed("consensus", &alignment.consensus, molecule).is_ok());
        }
        let dna = Document::new("dna", "ACT").unwrap();
        let rna = Document::new_typed("rna", "ACU", Molecule::Rna).unwrap();
        let aligned = align_pair(&dna, &rna, options).unwrap();
        assert_eq!(aligned.exact_matches, 3);
        assert_eq!(aligned.differences().lines().count(), 1);
        let protein_t = Document::new_typed("protein", "ACT", Molecule::Protein).unwrap();
        let protein_u = Document::new_typed("protein", "ACU", Molecule::Protein).unwrap();
        assert_eq!(
            align_pair(&protein_t, &protein_u, options)
                .unwrap()
                .differences()
                .lines()
                .count(),
            2
        );
    }
    #[test]
    fn multiple_output_budget_rejects_before_large_projection_allocation() {
        let long = Document::new("long", &"A".repeat(2_000_000)).unwrap();
        let short = Document::new("short", "A").unwrap();
        for inputs in [
            vec![long.clone(), short.clone(), short.clone()],
            vec![short.clone(), long, short],
        ] {
            let message = align_multiple(&inputs, AlignmentOptions::default())
                .unwrap_err()
                .to_string();
            assert!(message.contains("five million displayed cells"));
        }
    }
    #[test]
    fn nullable_quality_follows_original_reverse_coordinates() {
        let mut query = Document::new("query", "AAAC").unwrap();
        query.sequence_qualities = Some(vec![Some(30), None, Some(0), Some(40)]);
        let reference = Document::new("reference", "GTTT").unwrap();
        let aligned = align_pair(&reference, &query, AlignmentOptions::default()).unwrap();
        assert!(aligned.query.reverse);
        assert_eq!(
            aligned.query.qualities,
            vec![Some(40), Some(0), None, Some(30)]
        );
        assert_ne!(aligned.query.qualities[1], aligned.query.qualities[2]);
    }
}
