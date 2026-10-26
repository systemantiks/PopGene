//! Sequence products for explicitly specified PCR and cloning workflows.
//! Chemistry/geometry references: NEB cloning guide, compatible cohesive ends,
//! Gibson Assembly and Golden Gate manuals. These are sequence constructions,
//! not predictions of reaction yield, polymerase errors or ligase fidelity.
//! https://www.neb.com/tools-and-resources/usage-guidelines/cloning-guide
//! https://www.neb.com/en-ca/tools-and-resources/selection-charts/compatible-cohesive-ends-and-generation-of-new-restriction-sites
//! https://intl.neb.com/en-us/-/media/nebus/files/manuals/manuale1602.pdf
use crate::restriction::{self, DigestConditions, ENZYMES, EndKind, FragmentEnd};
use crate::workflow::{self, DuplexEnds};
use crate::{Document, Error, Feature, Molecule, Segment, reverse_complement};

fn error(message: impl Into<String>) -> Error {
    Error::Operation(message.into())
}
fn dna(document: &Document) -> Result<(), Error> {
    document.validate()?;
    if document.molecule != Molecule::Dna || document.single_stranded || document.trace.is_some() {
        return Err(error("Cloning requires double-stranded DNA documents"));
    }
    Ok(())
}
fn plain_end() -> FragmentEnd {
    FragmentEnd {
        kind: EndKind::Native,
        overhang: String::new(),
        enzymes: vec![],
    }
}
fn ends(document: &Document) -> DuplexEnds {
    document.duplex_ends.clone().unwrap_or_else(|| DuplexEnds {
        left: plain_end(),
        right: plain_end(),
        opposite_sequence: reverse_complement(&document.sequence),
        left_modification: None,
        right_modification: None,
    })
}

#[derive(Clone, Debug)]
pub struct AssemblyInput {
    pub document: Document,
    pub reverse: bool,
}
#[derive(Clone, Copy, Debug, Default)]
pub struct LigationOptions {
    pub circular: bool,
    /// Explicitly model treatment of otherwise unblocked 5′ ends with kinase.
    pub phosphorylate: bool,
    /// Permit one sealed strand per junction, with the remaining nick recorded.
    pub allow_nicks: bool,
}
#[derive(Clone, Copy, Debug, Default)]
pub struct PcrOptions {
    pub forward_phosphorylated: bool,
    pub reverse_phosphorylated: bool,
}

/// Common Golden Gate Type IIS enzymes explicitly documented by NEB. External
/// offsets alone do not establish enzyme class (EcoP15I, for example, is Type III).
/// https://www.neb.com/en-us/faqs/which-restriction-enzymes-are-used-in-golden-gate-assembly
/// https://media.neb.com/m/bd28765a734dba9/original/us_neb_catalog_technicalreference_2025_26.pdf
pub const GOLDEN_GATE_ENZYMES: &[&str] = &[
    "BbsI", "BsaI", "BsmBI", "BspQI", "BtgZI", "Esp3I", "PaqCI", "SapI", "AarI",
];

fn shifted_features(product: &mut Document, source: &Document, offset: usize) -> Result<(), Error> {
    product
        .sequence_colors
        .extend(source.sequence_colors.iter().cloned().map(|mut color| {
            color.range.start += offset;
            color.range.end += offset;
            color
        }));
    let mut references = Vec::new();
    for reference in &source.references {
        let i = product
            .references
            .iter()
            .position(|r| r == reference)
            .unwrap_or_else(|| {
                product.references.push(reference.clone());
                product.references.len() - 1
            });
        references.push(i + 1);
    }
    let mut lineage_bindings = Vec::new();
    for source_feature in &source.features {
        let mut feature = source_feature.clone();
        feature.id = (product.features.len() + 1) as u64;
        lineage_bindings.push((source_feature.id, feature.id));
        feature.project_display_segments(|s| {
            vec![Segment {
                start: s.start + offset,
                end: s.end + offset,
            }]
        });
        for part in &mut feature.parts {
            part.start += offset;
            part.end += offset;
        }
        if let Some(location) = &feature.location {
            feature.location = location
                .replaced(0, 0, offset)
                .map_err(|e| error(e.to_string()))?;
        }
        crate::translation::project_translation_exceptions(&mut feature, |p| p.checked_add(offset));
        for (key, value) in &mut feature.qualifiers {
            if key == "citation"
                && let Some(value) = value
            {
                let original = value
                    .trim()
                    .trim_start_matches('[')
                    .trim_end_matches(']')
                    .parse::<usize>()
                    .ok();
                if let Some(number) = original
                    .and_then(|i| i.checked_sub(1))
                    .and_then(|i| references.get(i))
                {
                    *value = format!("[{number}]");
                } else {
                    *key = "source_citation".into();
                }
            }
        }
        feature
            .qualifiers
            .push(("assembly_source".into(), Some(source.name.clone())));
        product.features.push(feature);
    }
    crate::lineage::append(source, product, offset, &lineage_bindings)?;
    for warning in &source.warnings {
        if !product.warnings.contains(warning) {
            product.warnings.push(warning.clone());
        }
    }
    product.comments.extend(
        source
            .comments
            .iter()
            .map(|c| format!("Source {}: {c}", source.name)),
    );
    Ok(())
}

/// Reverse an entire physical duplex, including annotations on either overhang.
pub fn orient(document: &Document, reverse: bool) -> Result<Document, Error> {
    dna(document)?;
    if !reverse {
        return Ok(document.clone());
    }
    let mut copy = document.clone();
    copy.revcomp()?;
    copy.comments
        .push(format!("Physical reverse orientation of {}", document.name));
    Ok(copy)
}

fn primer_binding(
    template: &Document,
    primer: &Feature,
) -> Result<(usize, usize, String, String), Error> {
    if primer.kind != "primer_bind"
        || primer.parts.is_empty()
        || primer
            .location
            .as_ref()
            .is_some_and(|l| l.uniform_projection().is_none())
    {
        return Err(error("Select bound primers with a definite local interval"));
    }
    let n = template.sequence.len();
    if let Some(ends) = &template.duplex_ends {
        let (top_prefix, top_suffix, _, _) = ends.overhang_lengths();
        if primer
            .parts
            .iter()
            .any(|p| p.start < top_prefix || p.end > n - top_suffix)
        {
            return Err(error(
                "PCR binding intervals must lie in the paired template region; annealing to single-stranded terminal extensions is not modeled",
            ));
        }
    }
    let mut length = 0;
    for (i, p) in primer.parts.iter().enumerate() {
        if p.start >= p.end || p.end > n {
            return Err(Error::Coordinates);
        }
        if i > 0 {
            let previous = &primer.parts[i - 1];
            if previous.end != p.start && !(template.circular && previous.end == n && p.start == 0)
            {
                return Err(error(
                    "PCR primer binding intervals must be contiguous, optionally crossing a circular origin",
                ));
            }
        }
        length += p.end - p.start;
    }
    if length < 2 || length > n {
        return Err(error(
            "Primer annealing interval must contain at least two bases and at most one circumference",
        ));
    }
    let oligo = primer
        .primer_sequence()
        .ok_or_else(|| {
            error(format!(
                "Primer {} has no stored oligo sequence",
                primer.label
            ))
        })?
        .to_ascii_uppercase();
    if !oligo.bytes().all(|b| b"ACGT".contains(&b)) || oligo.len() < length || oligo.len() > 10_000
    {
        return Err(error(
            "PCR oligos must contain 2–10,000 canonical DNA bases and cover their entire annealing interval",
        ));
    }
    let bound = primer.extract(&template.sequence).to_ascii_uppercase();
    if !oligo.ends_with(&bound) {
        return Err(error(format!(
            "Primer {} does not have a perfectly matching 3′ annealing suffix; mismatched PCR products require an explicit mismatch model",
            primer.label
        )));
    }
    let tail = oligo[..oligo.len() - length].to_string();
    // All intended positions are explicit, but additional exact binding sites
    // would make the product ambiguous, so do not silently select one.
    let target = if primer.reverse {
        reverse_complement(&bound)
    } else {
        bound.clone()
    };
    let opposite = reverse_complement(&target);
    let sites = exact_binding_starts(&template.sequence, &target, template.circular);
    if sites != [primer.parts[0].start]
        || !exact_binding_starts(&template.sequence, &opposite, template.circular).is_empty()
    {
        return Err(error(format!(
            "Primer {} does not have one unique exact binding site; PCR product is ambiguous",
            primer.label
        )));
    }
    Ok((primer.parts[0].start, length, oligo, tail))
}

/// KMP search retains at most two positions: enough to establish uniqueness,
/// including overlapping matches, without quadratic repetitive-template scans.
fn exact_binding_starts(sequence: &str, query: &str, circular: bool) -> Vec<usize> {
    let query = query.as_bytes();
    let mut prefix = vec![0; query.len()];
    for i in 1..query.len() {
        let mut j = prefix[i - 1];
        while j > 0 && query[i] != query[j] {
            j = prefix[j - 1];
        }
        if query[i] == query[j] {
            j += 1;
        }
        prefix[i] = j;
    }
    let mut matches = vec![];
    let mut matched = 0;
    let extent = sequence.len() + if circular { query.len() - 1 } else { 0 };
    for i in 0..extent {
        let base = sequence.as_bytes()[i % sequence.len()].to_ascii_uppercase();
        while matched > 0 && base != query[matched] {
            matched = prefix[matched - 1];
        }
        if base == query[matched] {
            matched += 1;
        }
        if matched == query.len() {
            matches.push(i + 1 - query.len());
            if matches.len() == 2 {
                break;
            }
            matched = prefix[matched - 1];
        }
    }
    matches
}

/// Mature double-stranded PCR product. Explicit 5′ tails are incorporated;
/// template sequence is projected once, including circular-origin crossings.
pub fn pcr_product(
    template: &Document,
    forward: &Feature,
    reverse: &Feature,
    options: PcrOptions,
    name: impl Into<String>,
) -> Result<Document, Error> {
    dna(template)?;
    if forward.reverse || !reverse.reverse {
        return Err(error(
            "PCR requires a forward primer and an opposing reverse primer",
        ));
    }
    let (fs, fl, fo, ft) = primer_binding(template, forward)?;
    let (rs, rl, ro, rt) = primer_binding(template, reverse)?;
    let n = template.sequence.len();
    let end = rs + rl;
    let length = if template.circular {
        let distance = (end + n - fs) % n;
        if distance == 0 { n } else { distance }
    } else {
        end.checked_sub(fs)
            .ok_or_else(|| error("Primers face away from the selected product"))?
    };
    if length < fl + rl || length > n {
        return Err(error(
            "Primer binding intervals overlap or do not define one inward-facing amplicon",
        ));
    }
    let reference = workflow::duplex_reference(template)?;
    let source_origin = template
        .duplex_ends
        .as_ref()
        .map_or(0, |e| e.geometry(n).top.start);
    let mut product = workflow::extract_region(&reference, fs + source_origin, length, name)?;
    // These are newly synthesized molecules, not additional observations of
    // the template's per-base read quality.
    product.sequence_qualities = None;
    product.opposite_sequence_qualities = None;
    if template.sequence_qualities.is_some() || template.opposite_sequence_qualities.is_some() {
        product.warnings.push("Template read qualities are not measurements of newly synthesized PCR molecules; product quality scores remain unknown.".into());
    }
    product.duplex_ends = None;
    product.circular = false;
    product.source = None;
    product
        .features
        .retain(|f| f.id != forward.id && f.id != reverse.id);
    if !ft.is_empty() {
        product.replace(0, 0, &ft)?;
    }
    if !rt.is_empty() {
        let n = product.sequence.len();
        product.replace(n, n, &reverse_complement(&rt))?;
    }
    crate::lineage::finish(&mut product, false)?;
    let final_length = product.sequence.len();
    for (original, sequence, part) in [
        (
            forward,
            &fo,
            Segment {
                start: 0,
                end: fo.len(),
            },
        ),
        (
            reverse,
            &ro,
            Segment {
                start: final_length - ro.len(),
                end: final_length,
            },
        ),
    ] {
        let mut f = original.clone();
        f.id = 0;
        // Preserve validated display subdivisions in product coordinates, and
        // extend the terminal subdivision over the newly incorporated 5′ tail.
        f.project_display_segments(|s| {
            let start = (s.start + n - fs) % n + ft.len();
            let end = start + s.end - s.start;
            vec![Segment {
                start: if !original.reverse && start == ft.len() {
                    0
                } else {
                    start
                },
                end: if original.reverse && end == ft.len() + length {
                    final_length
                } else {
                    end
                },
            }]
        });
        f.parts = vec![part];
        f.location = None;
        f.edited = true;
        f.qualifiers
            .retain(|(k, _)| !k.starts_with("binding_") && k != "primer_sequence");
        f.qualifiers
            .push(("primer_sequence".into(), Some(sequence.clone())));
        f.qualifiers.push(("PCR_incorporated".into(), None));
        product.features.push(f);
    }
    let lineage_ids = product
        .features
        .iter()
        .enumerate()
        .filter(|(_, f)| f.id != 0)
        .map(|(i, f)| (f.id, i as u64 + 1))
        .collect::<Vec<_>>();
    for piece in &mut product.annotation_lineage.pieces {
        if let Some((_, id)) = lineage_ids.iter().find(|(id, _)| *id == piece.feature_id) {
            piece.feature_id = *id;
        }
    }
    for (i, f) in product.features.iter_mut().enumerate() {
        f.id = i as u64 + 1;
    }
    crate::lineage::finish(&mut product, false)?;
    crate::lineage::new_copy(&mut product);
    let blunt = FragmentEnd {
        kind: EndKind::Blunt,
        overhang: String::new(),
        enzymes: vec![],
    };
    let chemistry = |p| Some(if p { "5′ phosphate" } else { "5′ hydroxyl" }.into());
    product.duplex_ends = Some(DuplexEnds {
        left: blunt.clone(),
        right: blunt,
        opposite_sequence: reverse_complement(&product.sequence),
        left_modification: chemistry(options.forward_phosphorylated),
        right_modification: chemistry(options.reverse_phosphorylated),
    });
    product.comments.push(format!("PCR: template {}; forward {} ({}-nt annealing, {}-nt 5′ tail); reverse {} ({}-nt annealing, {}-nt 5′ tail). Perfect-match suffix model, explicit primer 5′ phosphorylation F={} R={}; mature blunt product. Polymerase errors, A-tailing and reaction yield are not predicted.",template.name,forward.label,fl,ft.len(),reverse.label,rl,rt.len(),options.forward_phosphorylated,options.reverse_phosphorylated));
    product.validate()?;
    Ok(product)
}

fn chemistry(end: &FragmentEnd, modification: Option<&str>, kinase: bool) -> Result<bool, Error> {
    let state = match modification {
        Some("5′ phosphate" | "5' phosphate" | "5'-phosphate" | "phosphorylated") => Some(true),
        Some("5′ hydroxyl" | "5' hydroxyl" | "dephosphorylated") => Some(false),
        Some(value) => {
            return Err(error(format!(
                "Terminal modification '{value}' is not a defined ligatable 5′ end; resolve its chemistry first"
            )));
        }
        None if !end.enzymes.is_empty() => Some(true),
        None => None,
    };
    if kinase {
        return Ok(true);
    }
    state.ok_or_else(||error("5′ phosphorylation is unspecified; explicitly choose kinase treatment or use documented restriction/PCR ends"))
}
fn geometry(right: &FragmentEnd, left: &FragmentEnd) -> bool {
    let blunt = |e: &FragmentEnd| matches!(e.kind, EndKind::Blunt | EndKind::Native);
    if blunt(right) && blunt(left) {
        return true;
    }
    right.kind == left.kind
        && !right.overhang.is_empty()
        && right
            .overhang
            .bytes()
            .chain(left.overhang.bytes())
            .all(|b| b"ACGT".contains(&b))
        && reverse_complement(&right.overhang) == left.overhang
}
fn junction(a: &DuplexEnds, b: &DuplexEnds, options: LigationOptions) -> Result<bool, Error> {
    if !geometry(&a.right, &b.left) {
        return Err(error(
            "Adjacent ends have incompatible polarity or noncomplementary cohesive sequences",
        ));
    }
    let p = chemistry(
        &a.right,
        a.right_modification.as_deref(),
        options.phosphorylate,
    )?;
    let q = chemistry(
        &b.left,
        b.left_modification.as_deref(),
        options.phosphorylate,
    )?;
    if !(p && q) && !(options.allow_nicks && (p || q)) {
        return Err(error(
            "This junction cannot be sealed with the declared 5′ chemistry; phosphorylate the fragments or explicitly allow a one-strand nick",
        ));
    }
    Ok(!(p && q))
}

/// Ordered, oriented ligation. Input order/direction is explicit; cohesive
/// geometry and both 5′ end chemistries are checked at every joined boundary.
pub fn ligate(
    inputs: &[AssemblyInput],
    options: LigationOptions,
    name: impl Into<String>,
) -> Result<Document, Error> {
    if inputs.is_empty() || inputs.len() > 32 {
        return Err(error("Ligation requires 1–32 fragments"));
    }
    let documents = inputs
        .iter()
        .map(|i| orient(&i.document, i.reverse))
        .collect::<Result<Vec<_>, _>>()?;
    if documents.iter().any(|d| d.circular) {
        return Err(error("Linearize circular inputs before ligation"));
    }
    let all_ends = documents.iter().map(ends).collect::<Vec<_>>();
    if options.phosphorylate {
        for e in &all_ends {
            chemistry(&e.left, e.left_modification.as_deref(), true)?;
            chemistry(&e.right, e.right_modification.as_deref(), true)?;
        }
    }
    let mut nicks = 0;
    for pair in all_ends.windows(2) {
        nicks += usize::from(junction(&pair[0], &pair[1], options)?);
    }
    if options.circular {
        nicks += usize::from(junction(all_ends.last().unwrap(), &all_ends[0], options)?);
    }
    let total = documents.iter().map(|d| d.sequence.len()).sum::<usize>();
    if total > crate::MAX_BASES {
        return Err(Error::Length);
    }
    let sequence = documents
        .iter()
        .map(|d| d.sequence.as_ref())
        .collect::<String>();
    let mut product = Document::new(name, &sequence)?;
    if documents.iter().any(|d| d.sequence_qualities.is_some()) {
        product.sequence_qualities = Some(
            documents
                .iter()
                .flat_map(|d| {
                    d.sequence_qualities
                        .clone()
                        .unwrap_or_else(|| vec![None; d.sequence.len()])
                })
                .collect(),
        );
    }
    if documents
        .iter()
        .any(|d| d.opposite_sequence_qualities.is_some())
    {
        product.opposite_sequence_qualities = Some(
            documents
                .iter()
                .zip(&all_ends)
                .rev()
                .flat_map(|(d, e)| {
                    d.opposite_sequence_qualities
                        .clone()
                        .unwrap_or_else(|| vec![None; e.opposite_sequence.len()])
                })
                .collect(),
        );
    }
    let extent_origin = all_ends[0].geometry(documents[0].sequence.len()).top.start;
    let mut offset = 0;
    for document in &documents {
        let reference = workflow::duplex_reference(document)?;
        let source_origin = ends(document).geometry(document.sequence.len()).top.start;
        let placement = (offset + extent_origin)
            .checked_sub(source_origin)
            .ok_or(Error::Coordinates)?;
        shifted_features(&mut product, &reference, placement)?;
        offset += document.sequence.len();
    }
    let first = &all_ends[0];
    let last = all_ends.last().unwrap();
    let joined_ends = DuplexEnds {
        left: first.left.clone(),
        right: last.right.clone(),
        opposite_sequence: all_ends
            .iter()
            .rev()
            .map(|e| e.opposite_sequence.as_str())
            .collect(),
        left_modification: if options.phosphorylate {
            Some("5′ phosphate".into())
        } else {
            first.left_modification.clone()
        },
        right_modification: if options.phosphorylate {
            Some("5′ phosphate".into())
        } else {
            last.right_modification.clone()
        },
    };
    joined_ends.validate(&product.sequence)?;
    product.circular = options.circular;
    if options.circular {
        if let Some(scores) = product.opposite_sequence_qualities.as_mut() {
            let geometry = joined_ends.geometry(product.sequence.len());
            let shift = (geometry.bottom.end as isize - geometry.top.end as isize)
                .rem_euclid(scores.len() as isize) as usize;
            scores.rotate_left(shift);
        }
        wrap_features(&mut product, extent_origin)?;
    } else {
        product.duplex_ends = Some(joined_ends);
        workflow::partition_duplex_annotations(&mut product)?;
    }
    crate::lineage::reconstruct(&mut product, crate::lineage::Merge::Physical)?;
    product.comments.push(format!("Ordered ligation: {}. Circular {}; kinase treatment {}; permitted nicked junctions {}; {} junction(s) retain a nick. Input direction was explicit; no ligase fidelity, yield or host repair efficiency is predicted.",inputs.iter().map(|i|format!("{} ({})",i.document.name,if i.reverse{"reverse"}else{"forward"})).collect::<Vec<_>>().join(" → "),options.circular,options.phosphorylate,options.allow_nicks,nicks));
    if nicks > 0 {
        product.warnings.push("Product contains explicitly permitted nicked junctions; complete covalent closure requires repair.".into());
    }
    product.validate()?;
    Ok(product)
}

fn unique_overlap(a: &str, b: &str, min: usize, max: usize) -> Result<usize, Error> {
    let matches = (min..=max.min(a.len()).min(b.len()))
        .filter(|&n| {
            a[a.len() - n..] == b[..n]
                && a.as_bytes()[a.len() - n..]
                    .iter()
                    .all(|v| b"ACGT".contains(v))
        })
        .collect::<Vec<_>>();
    match matches.as_slice() {
        [n] => Ok(*n),
        [] => Err(error(
            "No exact terminal overlap satisfies the requested range",
        )),
        _ => Err(error(
            "Repeated terminal sequence permits multiple overlap lengths; assembly is ambiguous",
        )),
    }
}

fn wrap_features(product: &mut Document, origin: usize) -> Result<(), Error> {
    crate::lineage::wrap(product, origin);
    use crate::FeatureLocation as L;
    fn ranges(s: &Segment, n: usize, origin: usize) -> Vec<Segment> {
        let mut start = (s.start + n - origin) % n;
        let mut remaining = s.end - s.start;
        let mut result = vec![];
        while remaining > 0 {
            let length = remaining.min(n - start);
            result.push(Segment {
                start,
                end: start + length,
            });
            remaining -= length;
            start = 0;
        }
        result
    }
    fn location(l: &L, n: usize, origin: usize) -> L {
        match l {
            L::Range {
                start,
                end,
                before,
                after,
            } => {
                let mapped = ranges(
                    &Segment {
                        start: *start,
                        end: *end,
                    },
                    n,
                    origin,
                );
                let last = mapped.len() - 1;
                let nodes = mapped
                    .into_iter()
                    .enumerate()
                    .map(|(i, p)| L::Range {
                        start: p.start,
                        end: p.end,
                        before: i == 0 && *before,
                        after: i == last && *after,
                    })
                    .collect::<Vec<_>>();
                if nodes.len() == 1 {
                    nodes.into_iter().next().unwrap()
                } else {
                    L::Join(nodes)
                }
            }
            L::Between { left, right } => L::Between {
                left: (left + n - origin) % n,
                right: (right + n - origin) % n,
            },
            L::Complement(c) => L::Complement(Box::new(location(c, n, origin))),
            L::Join(c) => L::Join(c.iter().map(|c| location(c, n, origin)).collect()),
            L::Order(c) => L::Order(c.iter().map(|c| location(c, n, origin)).collect()),
            L::Bond(c) => L::Bond(c.iter().map(|c| location(c, n, origin)).collect()),
            L::OneOf(c) => L::OneOf(c.iter().map(|c| location(c, n, origin)).collect()),
            L::Remote { .. } | L::Gap { .. } => l.clone(),
        }
    }
    let n = product.sequence.len();
    if n == 0 || origin >= n {
        return Err(Error::Coordinates);
    }
    let map = |s: &Segment| ranges(s, n, origin);
    for f in &mut product.features {
        f.project_display_segments(map);
        if let Some(l) = &f.location {
            let l = location(l, n, origin);
            (f.parts, f.reverse) = l
                .uniform_projection()
                .unwrap_or_else(|| (l.local_parts(), false));
            f.location = Some(l);
        } else {
            f.parts = f.parts.iter().flat_map(map).collect();
        }
        crate::translation::project_translation_exceptions(f, |p| Some((p + n - origin) % n));
    }
    product.sequence_colors = product
        .sequence_colors
        .iter()
        .flat_map(|color| {
            map(&color.range)
                .into_iter()
                .map(|range| crate::metadata::SequenceColor {
                    range,
                    reverse: color.reverse,
                    color_rgb: color.color_rgb,
                })
        })
        .collect();
    Ok(())
}

/// Ordered exact-overlap assembly, analogous to Gibson assembly's intended
/// sequence product. Every chosen junction must have a unique terminal overlap.
pub fn assemble_overlaps(
    inputs: &[AssemblyInput],
    min_overlap: usize,
    max_overlap: usize,
    circular: bool,
    name: impl Into<String>,
) -> Result<Document, Error> {
    if inputs.len() < 2
        || inputs.len() > 32
        || min_overlap < 10
        || max_overlap < min_overlap
        || max_overlap > 1000
    {
        return Err(error(
            "Overlap assembly requires 2–32 fragments and 10–1,000-base overlap bounds",
        ));
    }
    let documents = inputs
        .iter()
        .map(|i| orient(&i.document, i.reverse))
        .collect::<Result<Vec<_>, _>>()?;
    for document in &documents {
        if document.circular {
            return Err(error("Overlap assembly inputs must be linear"));
        }
        if let Some(e) = &document.duplex_ends
            && e.overhang_lengths() != (0, 0, 0, 0)
        {
            return Err(error(
                "Exact-overlap assembly currently requires fully paired inputs; prepare cohesive ends first",
            ));
        }
    }
    let mut sequence = documents[0].sequence.to_string();
    let mut placements = vec![0];
    let mut overlaps = vec![];
    for pair in documents.windows(2) {
        let overlap = unique_overlap(
            &pair[0].sequence,
            &pair[1].sequence,
            min_overlap,
            max_overlap,
        )?;
        if overlap == pair[0].sequence.len() || overlap == pair[1].sequence.len() {
            return Err(error(
                "An input is wholly contained in an overlap; choose explicit distinct fragments",
            ));
        }
        if sequence.len() + pair[1].sequence.len() - overlap > crate::MAX_BASES {
            return Err(Error::Length);
        }
        placements.push(sequence.len() - overlap);
        sequence.push_str(&pair[1].sequence[overlap..]);
        overlaps.push(overlap);
    }
    let closure = if circular {
        unique_overlap(
            &documents.last().unwrap().sequence,
            &documents[0].sequence,
            min_overlap,
            max_overlap,
        )?
    } else {
        0
    };
    if closure > 0 {
        if closure == documents[0].sequence.len()
            || closure == documents.last().unwrap().sequence.len()
        {
            return Err(error(
                "An input is wholly contained in the circular closing overlap; choose explicit distinct fragments",
            ));
        }
        if closure >= sequence.len() {
            return Err(error("Circular closure consumes the entire product"));
        }
        sequence.truncate(sequence.len() - closure);
    }
    if documents.iter().any(|d| d.sequence.len() > sequence.len()) {
        return Err(error(
            "An input would span more than one product circumference; the assembly is ambiguous",
        ));
    }
    let mut product = Document::new(name, &sequence)?;
    for (document, offset) in documents.iter().zip(placements) {
        shifted_features(&mut product, &workflow::duplex_reference(document)?, offset)?;
    }
    if documents
        .iter()
        .any(|d| d.sequence_qualities.is_some() || d.opposite_sequence_qualities.is_some())
    {
        product.warnings.push("Source read qualities are not measurements of the newly synthesized overlap-assembly product; product quality scores remain unknown.".into());
    }
    product.circular = circular;
    if circular {
        wrap_features(&mut product, 0)?;
    }
    crate::lineage::reconstruct(&mut product, crate::lineage::Merge::VerifiedOverlap)?;
    product.comments.push(format!("Ordered exact-overlap assembly (Gibson-style sequence model): {}; overlaps {:?}; circular closure {} bases. Exact homology only; no mismatch repair, exonuclease kinetics, synthesis errors or reaction yield are inferred.",documents.iter().map(|d|d.name.as_str()).collect::<Vec<_>>().join(" → "),overlaps,closure));
    product.validate()?;
    Ok(product)
}

/// Produce selectable digestion fragments for each Golden Gate input.
pub fn golden_gate_fragments(
    inputs: &[AssemblyInput],
    enzyme_index: usize,
    conditions: DigestConditions,
) -> Result<Vec<Vec<Document>>, Error> {
    let e = ENZYMES
        .get(enzyme_index)
        .ok_or_else(|| error("Invalid Golden Gate enzyme"))?;
    if !GOLDEN_GATE_ENZYMES.contains(&e.name)
        || !e.cleavage_known
        || e.second_cut.is_some()
        || e.top_cut <= e.recognition_length() as i32
        || e.bottom_cut <= e.recognition_length() as i32
        || e.top_cut == e.bottom_cut
    {
        return Err(error(
            "Choose a supported Golden Gate Type IIS enzyme: BbsI, BsaI, BsmBI, BspQI, BtgZI, Esp3I, PaqCI, SapI or AarI",
        ));
    }
    if inputs.is_empty() || inputs.len() > 16 {
        return Err(error("Golden Gate requires 1–16 input documents"));
    }
    inputs.iter().map(|input|{
        dna(&input.document)?;let digest=restriction::digest_document(&input.document,&[enzyme_index],conditions).map_err(error)?;
        if !digest.all_sites_resolved {return Err(error("Golden Gate input has unresolved cleavage sites or methylation assumptions"));}
        if digest.fragments.len()>128 {return Err(error("Golden Gate input exceeds 128 digestion fragments; use a smaller defined construct"));}
        workflow::digest_products(&input.document,&digest.fragments,&[e.name.into()],conditions)
    }).collect()
}
/// Explicit fragment choice for each input, in the intended product order.
/// Remaining recognition sites reject the product because it could be recut.
pub fn golden_gate(
    inputs: &[AssemblyInput],
    fragments: &[usize],
    enzyme_index: usize,
    conditions: DigestConditions,
    name: impl Into<String>,
) -> Result<Document, Error> {
    if fragments.len() != inputs.len() {
        return Err(error("Choose exactly one digestion fragment per input"));
    }
    let candidates = golden_gate_fragments(inputs, enzyme_index, conditions)?;
    let selected = inputs
        .iter()
        .zip(candidates)
        .zip(fragments)
        .map(|((input, choices), &i)| {
            Ok(AssemblyInput {
                document: choices
                    .get(i)
                    .ok_or_else(|| error("Golden Gate fragment index is outside its input digest"))?
                    .clone(),
                reverse: input.reverse,
            })
        })
        .collect::<Result<Vec<_>, Error>>()?;
    let oriented = selected
        .iter()
        .map(|i| orient(&i.document, i.reverse).map(|d| ends(&d)))
        .collect::<Result<Vec<_>, _>>()?;
    for (i, e) in oriented.iter().enumerate() {
        if reverse_complement(&e.right.overhang) == e.right.overhang {
            return Err(error(
                "Palindromic Golden Gate overhang permits ambiguous self-ligation",
            ));
        }
        let mut partners = Vec::new();
        for (j, other) in oriented.iter().enumerate() {
            if geometry(&e.right, &other.left) {
                partners.push((j, false));
            }
            if i != j && geometry(&e.right, &other.right) {
                partners.push((j, true));
            }
        }
        if partners != vec![((i + 1) % oriented.len(), false)] {
            return Err(error(
                "Golden Gate ends do not specify one unique order and orientation; choose distinct complementary overhangs",
            ));
        }
    }
    let mut product = ligate(
        &selected,
        LigationOptions {
            circular: true,
            ..LigationOptions::default()
        },
        name,
    )?;
    if restriction::find_sites_for(&product.sequence, true, &[enzyme_index])
        .map_err(error)?
        .counts[enzyme_index]
        > 0
    {
        return Err(error(
            "Assembled product retains the chosen Type IIS recognition site and can be recut",
        ));
    }
    product.comments.push(format!("Golden Gate sequence construction: {}; explicit fragment choices {:?} (zero-based per input digest). All chosen junctions require complementary ends and documented 5′ phosphate chemistry. Mismatch ligation and enzymatic cycling efficiency are not modeled.",ENZYMES[enzyme_index].name,fragments));
    product.validate()?;
    Ok(product)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn input(document: Document, reverse: bool) -> AssemblyInput {
        AssemblyInput { document, reverse }
    }
    fn enzyme(name: &str) -> usize {
        ENZYMES.iter().position(|e| e.name == name).unwrap()
    }
    fn primer(
        id: u64,
        name: &str,
        sequence: &str,
        start: usize,
        length: usize,
        reverse: bool,
    ) -> Feature {
        Feature {
            id,
            label: name.into(),
            kind: "primer_bind".into(),
            parts: vec![Segment {
                start,
                end: start + length,
            }],
            location: None,
            reverse,
            qualifiers: vec![("primer_sequence".into(), Some(sequence.into()))],
            edited: false,
        }
    }
    #[test]
    fn pcr_tails_and_origin_crossing_preserve_the_insert_annotation() {
        let sequence =
            include_str!("../../../validation/fixtures/primer-design-template.txt").trim();
        let mut source = Document::new("PCR template", sequence).unwrap();
        source.circular = true;
        source.sequence_qualities = Some(vec![Some(20); sequence.len()]);
        source.opposite_sequence_qualities = Some(vec![Some(30); sequence.len()]);
        source
            .annotate(
                "insert".into(),
                Segment {
                    start: 220,
                    end: 230,
                },
                false,
            )
            .unwrap();
        source.references.push(crate::metadata::Reference {
            title: "source publication".into(),
            ..Default::default()
        });
        let mut rotated =
            workflow::extract_region(&source, 225, sequence.len(), "rotated template").unwrap();
        rotated.circular = true;
        rotated.methylation = Some(crate::metadata::DnaMethylation {
            dam: true,
            ..Default::default()
        });
        let mut f = primer(100, "Forward", "GGGCGGATACACTGCCGCAAGTTGGA", 414, 20, false);
        f.qualifiers.push(("note".into(), Some("This forward directional feature has 2 segments: 1: 415 .. 424 / #ff0000 / binding A 2: 425 .. 434 / #0000ff / binding B".into())));
        let r = primer(101, "Reverse", "TTAACCGGGCAAGCTGAGATGTCCTT", 169, 20, true);
        let result = pcr_product(&rotated, &f, &r, PcrOptions::default(), "PCR product").unwrap();
        assert_eq!(
            &*result.sequence,
            format!("GGGCGG{}GGTTAA", &sequence[139..414])
        );
        let feature = result
            .features
            .iter()
            .find(|f| f.label == "insert")
            .unwrap();
        assert_eq!(
            feature.parts,
            vec![
                Segment { start: 87, end: 92 },
                Segment { start: 92, end: 97 }
            ]
        );
        assert_eq!(
            result.feature_sequence(feature).unwrap(),
            &sequence[220..230]
        );
        assert_eq!(result.references, source.references);
        assert!(result.sequence_qualities.is_none());
        assert!(result.opposite_sequence_qualities.is_none());
        assert!(result.methylation.is_none());
        assert!(
            result
                .warnings
                .iter()
                .any(|w| w.contains("newly synthesized PCR"))
        );
        let incorporated = result
            .features
            .iter()
            .find(|p| p.label == "Forward")
            .unwrap();
        let displays = incorporated.display_segments();
        assert_eq!(displays.len(), 2);
        assert_eq!(displays[0].range, Segment { start: 0, end: 16 });
        assert_eq!(displays[1].range, Segment { start: 16, end: 26 });
        assert_eq!(displays[0].color_rgb, Some(0xff0000));
        assert_eq!(
            result
                .duplex_ends
                .as_ref()
                .unwrap()
                .left_modification
                .as_deref(),
            Some("5′ hydroxyl")
        );
        let mut wrong = r.clone();
        wrong.qualifiers[0].1 = Some("TTAACCGGGCAAGCTGAGATGTCCTA".into());
        assert!(pcr_product(&rotated, &f, &wrong, PcrOptions::default(), "bad").is_err());
    }
    #[test]
    fn exact_binding_uniqueness_handles_repeats_and_the_circular_origin() {
        assert_eq!(exact_binding_starts("AAAAC", "AAA", false), [0, 1]);
        assert_eq!(exact_binding_starts("AATGCC", "CCAAT", true), [4]);
        assert!(exact_binding_starts("AATGCC", "CCAAT", false).is_empty());
        assert!(
            exact_binding_starts(
                &"A".repeat(300_000),
                &format!("{}C", "A".repeat(9_999)),
                false
            )
            .is_empty()
        );
    }
    #[test]
    fn ligation_conserves_both_strands_and_orientation_and_checks_phosphates() {
        let source = Document::new("substrate", "AACGGAATTCTGGA").unwrap();
        let conditions = DigestConditions::default();
        let fragments =
            restriction::digest_document(&source, &[enzyme("EcoRI")], conditions).unwrap();
        let products =
            workflow::digest_products(&source, &fragments.fragments, &["EcoRI".into()], conditions)
                .unwrap();
        let joined = ligate(
            &[
                input(products[0].clone(), false),
                input(products[1].clone(), false),
            ],
            LigationOptions::default(),
            "joined",
        )
        .unwrap();
        assert_eq!(joined.sequence, source.sequence);
        assert_eq!(
            joined.duplex_ends.as_ref().unwrap().opposite_sequence,
            reverse_complement(&source.sequence)
        );
        let reversed = ligate(
            &[
                input(products[1].clone(), true),
                input(products[0].clone(), true),
            ],
            LigationOptions::default(),
            "reversed",
        )
        .unwrap();
        assert_eq!(&*reversed.sequence, reverse_complement(&source.sequence));
        let plain = [
            input(Document::new("a", "AACCGG").unwrap(), false),
            input(Document::new("b", "TTAGCA").unwrap(), false),
        ];
        assert!(ligate(&plain, LigationOptions::default(), "unknown chemistry").is_err());
        assert_eq!(
            &*ligate(
                &plain,
                LigationOptions {
                    phosphorylate: true,
                    ..Default::default()
                },
                "phosphorylated"
            )
            .unwrap()
            .sequence,
            "AACCGGTTAGCA"
        );
        let mut blocked = plain[0].clone();
        let mut blocked_ends = ends(&blocked.document);
        blocked_ends.right_modification = Some("5′ biotin".into());
        blocked.document.duplex_ends = Some(blocked_ends);
        assert!(
            ligate(
                &[blocked, plain[1].clone()],
                LigationOptions {
                    phosphorylate: true,
                    ..Default::default()
                },
                "blocked"
            )
            .is_err()
        );
        let mut first = plain[0].clone();
        let mut first_ends = ends(&first.document);
        first_ends.right_modification = Some("5′ phosphate".into());
        first.document.duplex_ends = Some(first_ends);
        let mut second = plain[1].clone();
        let mut second_ends = ends(&second.document);
        second_ends.left_modification = Some("5′ hydroxyl".into());
        second.document.duplex_ends = Some(second_ends);
        assert!(
            ligate(
                &[first.clone(), second.clone()],
                LigationOptions::default(),
                "unsealed"
            )
            .is_err()
        );
        let nicked = ligate(
            &[first, second],
            LigationOptions {
                allow_nicks: true,
                ..Default::default()
            },
            "nicked",
        )
        .unwrap();
        assert!(
            nicked
                .warnings
                .iter()
                .any(|w| w.contains("nicked junctions"))
        );
    }
    #[test]
    fn digestion_and_ligation_preserve_full_duplex_annotations_and_colors() {
        use crate::metadata::SequenceColor;
        for left in [EndKind::FivePrime, EndKind::ThreePrime] {
            for right in [EndKind::FivePrime, EndKind::ThreePrime] {
                let axis = "TTAACCGGAATTCTGGAACCGGTT";
                let top_start = usize::from(left == EndKind::ThreePrime) * 2;
                let top_end = axis.len() - usize::from(right == EndKind::FivePrime) * 3;
                let bottom_start = usize::from(left == EndKind::FivePrime) * 2;
                let bottom_end = axis.len() - usize::from(right == EndKind::ThreePrime) * 3;
                let mut source = Document::new("full duplex", &axis[top_start..top_end]).unwrap();
                source.duplex_ends = Some(DuplexEnds {
                    left: FragmentEnd {
                        kind: left,
                        overhang: if left == EndKind::FivePrime {
                            axis[..2].into()
                        } else {
                            reverse_complement(&axis[..2])
                        },
                        enzymes: vec![],
                    },
                    right: FragmentEnd {
                        kind: right,
                        overhang: if right == EndKind::ThreePrime {
                            axis[axis.len() - 3..].into()
                        } else {
                            reverse_complement(&axis[axis.len() - 3..])
                        },
                        enzymes: vec![],
                    },
                    opposite_sequence: reverse_complement(&axis[bottom_start..bottom_end]),
                    left_modification: Some("5′ phosphate".into()),
                    right_modification: Some("5′ phosphate".into()),
                });
                for (id, (start, end)) in [(0, 1), (1, 5), (axis.len() - 4, axis.len())]
                    .into_iter()
                    .enumerate()
                {
                    source.features.push(Feature {
                        id: id as u64 + 1,
                        label: format!("extent{id}"),
                        kind: "misc_feature".into(),
                        parts: vec![Segment { start, end }],
                        reverse: id % 2 == 0,
                        location: None,
                        qualifiers: vec![],
                        edited: false,
                    });
                }
                for (reverse, start, end) in [
                    (false, top_start, top_end),
                    (true, bottom_start, bottom_end),
                ] {
                    source.sequence_colors.push(SequenceColor {
                        range: Segment { start, end },
                        reverse,
                        color_rgb: if reverse { 0x123456 } else { 0x987654 },
                    });
                }
                workflow::partition_duplex_annotations(&mut source).unwrap();
                source.sequence_qualities = Some(
                    (0..source.sequence.len())
                        .map(|i| if i == 1 { None } else { Some(i as u8) })
                        .collect(),
                );
                source.opposite_sequence_qualities = Some(
                    (0..source.duplex_ends.as_ref().unwrap().opposite_sequence.len())
                        .map(|i| if i == 2 { None } else { Some(i as u8 + 30) })
                        .collect(),
                );
                let reference = workflow::duplex_reference(&source).unwrap();
                let digest = restriction::digest_document(
                    &source,
                    &[enzyme("EcoRI")],
                    DigestConditions::default(),
                )
                .unwrap();
                let fragments = workflow::digest_products(
                    &source,
                    &digest.fragments,
                    &["EcoRI".into()],
                    DigestConditions::default(),
                )
                .unwrap();
                assert!(fragments.iter().any(|d| !d.duplex_annotations.is_empty()));
                let joined = ligate(
                    &fragments
                        .into_iter()
                        .map(|document| input(document, false))
                        .collect::<Vec<_>>(),
                    LigationOptions::default(),
                    "rejoined",
                )
                .unwrap();
                assert_eq!(joined.sequence, source.sequence);
                assert_eq!(joined.sequence_qualities, source.sequence_qualities);
                assert_eq!(
                    joined.opposite_sequence_qualities,
                    source.opposite_sequence_qualities
                );
                assert_eq!(
                    joined.duplex_ends.as_ref().unwrap().opposite_sequence,
                    source.duplex_ends.as_ref().unwrap().opposite_sequence
                );
                let actual = workflow::duplex_reference(&joined).unwrap();
                for original in &reference.features {
                    let retained = actual
                        .features
                        .iter()
                        .find(|f| f.label == original.label)
                        .unwrap();
                    assert_eq!(
                        actual.feature_sequence(retained).unwrap(),
                        reference.feature_sequence(original).unwrap()
                    );
                }
                let colors = |d: &Document| {
                    d.sequence_colors
                        .iter()
                        .flat_map(|c| {
                            (c.range.start..c.range.end).map(move |p| (p, c.reverse, c.color_rgb))
                        })
                        .collect::<std::collections::BTreeSet<_>>()
                };
                assert_eq!(colors(&actual), colors(&reference));
                let restored = orient(&orient(&source, true).unwrap(), true).unwrap();
                assert_eq!(restored.sequence_qualities, source.sequence_qualities);
                assert_eq!(
                    restored.opposite_sequence_qualities,
                    source.opposite_sequence_qualities
                );
                let restored = workflow::duplex_reference(&restored).unwrap();
                let feature_map = |d: &Document| {
                    d.features
                        .iter()
                        .map(|f| (f.label.clone(), f.clone()))
                        .collect::<std::collections::BTreeMap<_, _>>()
                };
                assert_eq!(feature_map(&restored), feature_map(&reference));
                assert_eq!(colors(&restored), colors(&reference));
            }
        }
    }
    #[test]
    fn circular_ligation_rotates_each_measured_strand_quality_axis() {
        for (enzyme_name, sequence) in [
            ("EcoRI", "AACCGAATTCTTGGCAT"),
            ("KpnI", "AACCGGTACCTTGGCAT"),
        ] {
            let mut source = Document::new("quality plasmid", sequence).unwrap();
            source.circular = true;
            source.sequence_qualities = Some((0..sequence.len()).map(|i| Some(i as u8)).collect());
            source.opposite_sequence_qualities = Some(
                (0..sequence.len())
                    .map(|i| if i == 3 { None } else { Some(i as u8 + 30) })
                    .collect(),
            );
            let digest = restriction::digest_document(
                &source,
                &[enzyme(enzyme_name)],
                DigestConditions::default(),
            )
            .unwrap();
            assert_eq!(digest.fragments.len(), 1);
            let fragment = &digest.fragments[0];
            let products = workflow::digest_products(
                &source,
                &digest.fragments,
                &[enzyme_name.into()],
                DigestConditions::default(),
            )
            .unwrap();
            let closed = ligate(
                &[input(products[0].clone(), false)],
                LigationOptions {
                    circular: true,
                    ..Default::default()
                },
                "closed",
            )
            .unwrap();
            let mut expected_top = source.sequence_qualities.clone().unwrap();
            expected_top.rotate_left(fragment.start);
            let mut expected_bottom = source.opposite_sequence_qualities.clone().unwrap();
            expected_bottom.rotate_left((sequence.len() - fragment.start) % sequence.len());
            assert_eq!(closed.sequence_qualities, Some(expected_top));
            assert_eq!(closed.opposite_sequence_qualities, Some(expected_bottom));
            assert!(closed.circular && closed.duplex_ends.is_none());
        }
    }
    #[test]
    fn exact_overlaps_preserve_annotations_and_circular_closure() {
        let ab = "ACGTCAGTGCATCGATGCTA";
        let ba = "TGCACTAGGCTAACGTCGAT";
        let mut a = Document::new("a", &format!("{ba}CCGAGTAC{ab}")).unwrap();
        let mut b = Document::new("b", &format!("{ab}GGTCATCA{ba}")).unwrap();
        a.sequence_qualities = Some(vec![Some(0); a.sequence.len()]);
        b.opposite_sequence_qualities = Some(vec![Some(30); b.sequence.len()]);
        a.methylation = Some(crate::metadata::DnaMethylation {
            dam: true,
            ..Default::default()
        });
        let n = b.sequence.len();
        b.annotate(
            "closing feature".into(),
            Segment {
                start: n - 25,
                end: n,
            },
            false,
        )
        .unwrap();
        b.sequence_colors.push(crate::metadata::SequenceColor {
            range: Segment {
                start: n - 25,
                end: n,
            },
            reverse: true,
            color_rgb: 0x334455,
        });
        let expected_feature = b.feature_sequence(&b.features[0]).unwrap();
        let result = assemble_overlaps(
            &[input(a, false), input(b, false)],
            15,
            25,
            true,
            "assembly",
        )
        .unwrap();
        assert_eq!(&*result.sequence, format!("{ba}CCGAGTAC{ab}GGTCATCA"));
        assert!(result.circular);
        assert!(result.sequence_qualities.is_none());
        assert!(result.opposite_sequence_qualities.is_none());
        assert!(result.methylation.is_none());
        assert!(
            result
                .warnings
                .iter()
                .any(|w| w.contains("newly synthesized overlap-assembly"))
        );
        assert_eq!(result.features[0].parts.len(), 2);
        assert_eq!(result.sequence_colors.len(), 2);
        assert_eq!(result.sequence_colors[0].range, result.features[0].parts[0]);
        assert_eq!(result.sequence_colors[1].range, result.features[0].parts[1]);
        assert!(
            result
                .sequence_colors
                .iter()
                .all(|c| c.reverse && c.color_rgb == 0x334455)
        );
        assert_eq!(
            result.feature_sequence(&result.features[0]).unwrap(),
            expected_feature
        );
        let repeated = [
            input(
                Document::new("repeat a", &format!("CGT{}", "A".repeat(20))).unwrap(),
                false,
            ),
            input(
                Document::new("repeat b", &format!("{}CGT", "A".repeat(20))).unwrap(),
                false,
            ),
        ];
        assert!(assemble_overlaps(&repeated, 10, 20, false, "ambiguous").is_err());
    }
    #[test]
    fn golden_gate_uses_selected_type_iis_products_and_removes_recognition_sites() {
        let cassette = |name: &str, left: &str, body: &str, right: &str| {
            Document::new(name, &format!("AAAAGGTCTCA{left}{body}{right}TGAGACCAAAA")).unwrap()
        };
        let inputs = [
            input(cassette("vector", "AGTC", "CCGATATACCGG", "CGTA"), false),
            input(cassette("insert", "CGTA", "TTAAGGCCAT", "AGTC"), false),
        ];
        assert!(
            golden_gate_fragments(&inputs, enzyme("EcoP15I"), DigestConditions::default()).is_err()
        );
        let result = golden_gate(
            &inputs,
            &[1, 1],
            enzyme("BsaI"),
            DigestConditions::default(),
            "Golden Gate",
        )
        .unwrap();
        assert!(result.circular);
        assert_eq!(&*result.sequence, "AGTCCCGATATACCGGCGTATTAAGGCCAT");
        assert_eq!(
            restriction::find_sites_for(&result.sequence, true, &[enzyme("BsaI")])
                .unwrap()
                .total_sites,
            0
        );
        assert!(
            golden_gate(
                &inputs,
                &[0, 1],
                enzyme("BsaI"),
                DigestConditions::default(),
                "wrong fragment"
            )
            .is_err()
        );
    }

    #[test]
    fn golden_gate_transfers_end_only_annotations_citations_colors_and_strand_scores() {
        use crate::metadata::{Reference, SequenceColor};
        let mut a = Document::new("vector", "AAAAGGTCTCAAGTCCCGATATACCGGCGTATGAGACCAAAA").unwrap();
        let mut b = Document::new("insert", "AAAAGGTCTCACGTATTAAGGCCATAGTCTGAGACCAAAA").unwrap();
        assert_eq!((a.sequence.len(), b.sequence.len()), (42, 40));
        // BsaI selects top intervals [11,27) / [11,25), with opposite
        // intervals [15,31) / [15,29). Each annotation below occupies only
        // the right-hand opposite-strand extension after digestion.
        for (d, label, range, color) in [
            (&mut a, "junction", Segment { start: 27, end: 31 }, 0xff0000),
            (&mut b, "closure", Segment { start: 25, end: 29 }, 0x0000ff),
        ] {
            d.annotate(label.into(), range.clone(), true).unwrap();
            d.references.push(Reference {
                title: format!("{label} publication"),
                ..Default::default()
            });
            d.features[0]
                .qualifiers
                .push(("citation".into(), Some("[1]".into())));
            d.sequence_colors.push(SequenceColor {
                range,
                reverse: true,
                color_rgb: color,
            });
        }
        a.sequence_qualities = Some((0..42).map(|i| Some(i as u8)).collect());
        a.sequence_qualities.as_mut().unwrap()[12] = None;
        a.sequence_qualities.as_mut().unwrap()[13] = Some(0);
        b.sequence_qualities = Some((0..40).map(|i| Some(40 + i as u8)).collect());
        a.opposite_sequence_qualities = Some((0..42).map(|i| Some(90 - i as u8)).collect());
        a.opposite_sequence_qualities.as_mut().unwrap()[15] = Some(0);
        b.opposite_sequence_qualities = Some((0..40).map(|i| Some(50 - i as u8)).collect());
        b.opposite_sequence_qualities.as_mut().unwrap()[13] = None;
        let inputs = [input(a, false), input(b, false)];
        let candidates =
            golden_gate_fragments(&inputs, enzyme("BsaI"), DigestConditions::default()).unwrap();
        for candidates in &candidates {
            assert_eq!(candidates[1].duplex_annotations.features.len(), 1);
            assert_eq!(candidates[1].duplex_annotations.sequence_colors.len(), 1);
        }
        let product = golden_gate(
            &inputs,
            &[1, 1],
            enzyme("BsaI"),
            DigestConditions::default(),
            "assembly",
        )
        .unwrap();
        assert_eq!(&*product.sequence, "AGTCCCGATATACCGGCGTATTAAGGCCAT");
        assert!(product.circular && product.duplex_ends.is_none());
        let junction = product
            .features
            .iter()
            .find(|f| f.label == "junction")
            .unwrap();
        let closure = product
            .features
            .iter()
            .find(|f| f.label == "closure")
            .unwrap();
        assert_eq!(junction.parts, vec![Segment { start: 16, end: 20 }]);
        assert_eq!(closure.parts, vec![Segment { start: 0, end: 4 }]);
        assert_eq!(product.feature_sequence(junction).unwrap(), "TACG");
        assert_eq!(product.feature_sequence(closure).unwrap(), "GACT");
        assert_eq!(junction.qualifier("citation"), Some("[1]"));
        assert_eq!(closure.qualifier("citation"), Some("[2]"));
        assert_eq!(
            product
                .references
                .iter()
                .map(|r| r.title.as_str())
                .collect::<Vec<_>>(),
            ["junction publication", "closure publication"]
        );
        assert_eq!(
            product.sequence_colors,
            vec![
                SequenceColor {
                    range: Segment { start: 16, end: 20 },
                    reverse: true,
                    color_rgb: 0xff0000
                },
                SequenceColor {
                    range: Segment { start: 0, end: 4 },
                    reverse: true,
                    color_rgb: 0x0000ff
                },
            ]
        );
        assert_eq!(
            product.sequence_qualities,
            Some(vec![
                Some(11),
                None,
                Some(0),
                Some(14),
                Some(15),
                Some(16),
                Some(17),
                Some(18),
                Some(19),
                Some(20),
                Some(21),
                Some(22),
                Some(23),
                Some(24),
                Some(25),
                Some(26),
                Some(51),
                Some(52),
                Some(53),
                Some(54),
                Some(55),
                Some(56),
                Some(57),
                Some(58),
                Some(59),
                Some(60),
                Some(61),
                Some(62),
                Some(63),
                Some(64),
            ])
        );
        assert_eq!(
            product.opposite_sequence_qualities,
            Some(vec![
                Some(35),
                Some(34),
                Some(33),
                Some(32),
                Some(31),
                Some(30),
                Some(29),
                Some(28),
                Some(27),
                Some(26),
                Some(79),
                Some(78),
                Some(77),
                Some(76),
                Some(0),
                Some(74),
                Some(73),
                Some(72),
                Some(71),
                Some(70),
                Some(69),
                Some(68),
                Some(67),
                Some(66),
                Some(65),
                Some(64),
                Some(39),
                Some(38),
                None,
                Some(36),
            ])
        );
    }
}
