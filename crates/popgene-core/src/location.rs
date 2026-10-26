//! Lossless INSDC location semantics. Local coordinates are zero-based; ranges
//! are half-open and between-base endpoints identify the two flanking residues.
//! https://www.insdc.org/submitting-standards/feature-table/ (v11.4, 2026).
use crate::{MAX_BASES, Molecule, Segment};
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum FeatureLocation {
    Range {
        start: usize,
        end: usize,
        before: bool,
        after: bool,
    },
    Between {
        left: usize,
        right: usize,
    },
    Complement(Box<Self>),
    Join(Vec<Self>),
    Order(Vec<Self>),
    Bond(Vec<Self>),
    OneOf(Vec<Self>),
    Remote {
        accession: String,
        location: Option<Box<Self>>,
    },
    Gap {
        length: Option<usize>,
        unknown100: bool,
    },
}
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[error("{0}")]
pub struct LocationError(pub String);
fn error(message: &str) -> LocationError {
    LocationError(message.into())
}

impl FeatureLocation {
    /// Validate local coordinates without treating a remote accession as this
    /// molecule. Remote coordinates are checked when explicitly resolved.
    pub fn validate(&self, length: usize, circular: bool) -> Result<(), LocationError> {
        fn walk(
            node: &FeatureLocation,
            length: Option<usize>,
            circular: bool,
            depth: usize,
            count: &mut usize,
        ) -> Result<(), LocationError> {
            *count += 1;
            if depth > 64 || *count > 100_000 {
                return Err(error("Location nesting or size exceeds limits"));
            }
            match node {
                FeatureLocation::Range { start, end, .. } => {
                    if start >= end || *end > i64::MAX as usize || length.is_some_and(|n| *end > n)
                    {
                        return Err(error("Location range is outside its sequence"));
                    }
                }
                FeatureLocation::Between { left, right } => {
                    if left == right
                        || *left >= i64::MAX as usize
                        || *right >= i64::MAX as usize
                        || length.is_some_and(|n| *left >= n || *right >= n)
                        || !(left.checked_add(1) == Some(*right)
                            || (circular && *right == 0 && length == left.checked_add(1))
                            || (length.is_none() && *right == 0))
                    {
                        return Err(error(
                            "Between-base endpoints must be adjacent (or the circular origin)",
                        ));
                    }
                }
                FeatureLocation::Complement(child) => {
                    walk(child, length, circular, depth + 1, count)?
                }
                FeatureLocation::Join(children)
                | FeatureLocation::Order(children)
                | FeatureLocation::Bond(children)
                | FeatureLocation::OneOf(children) => {
                    if children.is_empty() {
                        return Err(error("Empty compound location"));
                    }
                    for child in children {
                        walk(child, length, circular, depth + 1, count)?;
                    }
                }
                FeatureLocation::Remote {
                    accession,
                    location,
                } => {
                    if accession.is_empty()
                        || !accession
                            .bytes()
                            .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'.'))
                    {
                        return Err(error("Invalid remote accession"));
                    }
                    if let Some(child) = location {
                        walk(child, None, false, depth + 1, count)?;
                    }
                }
                FeatureLocation::Gap { length, unknown100 } => {
                    if length.is_some_and(|n| n == 0 || n > MAX_BASES)
                        || (*unknown100 && length.is_some())
                    {
                        return Err(error("Invalid gap length"));
                    }
                }
            }
            Ok(())
        }
        walk(self, Some(length), circular, 0, &mut 0)
    }
    /// Exact local interbase sites. Remote coordinates never become local marks.
    pub fn local_boundaries(&self) -> Vec<usize> {
        fn walk(node: &FeatureLocation, out: &mut Vec<usize>) {
            match node {
                FeatureLocation::Between { left, .. } => {
                    if let Some(site) = left.checked_add(1) {
                        out.push(site);
                    }
                }
                FeatureLocation::Complement(child) => walk(child, out),
                FeatureLocation::Join(children)
                | FeatureLocation::Order(children)
                | FeatureLocation::Bond(children)
                | FeatureLocation::OneOf(children) => {
                    for child in children {
                        walk(child, out);
                    }
                }
                _ => {}
            }
        }
        let mut result = vec![];
        walk(self, &mut result);
        result
    }
    /// Local spans for display only. Their order never implies that order/bond/
    /// alternative locations are a joined biological sequence.
    pub fn local_parts(&self) -> Vec<Segment> {
        fn walk(node: &FeatureLocation, out: &mut Vec<Segment>) {
            match node {
                FeatureLocation::Range { start, end, .. } => out.push(Segment {
                    start: *start,
                    end: *end,
                }),
                FeatureLocation::Complement(child) => walk(child, out),
                FeatureLocation::Join(children)
                | FeatureLocation::Order(children)
                | FeatureLocation::Bond(children)
                | FeatureLocation::OneOf(children) => {
                    for child in children {
                        walk(child, out);
                    }
                }
                _ => {}
            }
        }
        let mut out = vec![];
        walk(self, &mut out);
        out
    }
    /// A projection suitable for the legacy single-strand parts model. Mixed
    /// strands and non-join operators have no such projection.
    pub fn uniform_projection(&self) -> Option<(Vec<Segment>, bool)> {
        fn walk(
            node: &FeatureLocation,
            reverse: bool,
            out: &mut Vec<(Segment, bool)>,
        ) -> Option<()> {
            match node {
                FeatureLocation::Range { start, end, .. } => out.push((
                    Segment {
                        start: *start,
                        end: *end,
                    },
                    reverse,
                )),
                FeatureLocation::Complement(child) => walk(child, !reverse, out)?,
                FeatureLocation::Join(children) => {
                    if reverse {
                        for child in children.iter().rev() {
                            walk(child, reverse, out)?;
                        }
                    } else {
                        for child in children {
                            walk(child, reverse, out)?;
                        }
                    }
                }
                _ => return None,
            }
            Some(())
        }
        let mut out = vec![];
        walk(self, false, &mut out)?;
        let reverse = out.first()?.1;
        if out.iter().any(|(_, strand)| *strand != reverse) {
            return None;
        }
        if reverse {
            out.reverse();
        }
        Some((out.into_iter().map(|(range, _)| range).collect(), reverse))
    }
    /// Partial biological 5′ and 3′ boundaries, after complement/join ordering.
    pub fn partial_ends(&self) -> (bool, bool) {
        match self {
            Self::Range { before, after, .. } => (*before, *after),
            Self::Complement(child) => {
                let (a, b) = child.partial_ends();
                (b, a)
            }
            Self::Join(children) => (
                children.first().is_none_or(|c| c.partial_ends().0),
                children.last().is_none_or(|c| c.partial_ends().1),
            ),
            Self::Remote {
                location: Some(child),
                ..
            } => child.partial_ends(),
            _ => (true, true),
        }
    }
    pub fn expression(&self) -> String {
        fn list(name: &str, children: &[FeatureLocation]) -> String {
            format!(
                "{name}({})",
                children
                    .iter()
                    .map(FeatureLocation::expression)
                    .collect::<Vec<_>>()
                    .join(",")
            )
        }
        match self {
            Self::Range {
                start,
                end,
                before,
                after,
            } => {
                if *end == start + 1 && !before && !after {
                    (start + 1).to_string()
                } else {
                    format!(
                        "{}{}..{}{}",
                        if *before { "<" } else { "" },
                        start + 1,
                        if *after { ">" } else { "" },
                        end
                    )
                }
            }
            Self::Between { left, right } => format!("{}^{}", left + 1, right + 1),
            Self::Complement(child) => format!("complement({})", child.expression()),
            Self::Join(children) => list("join", children),
            Self::Order(children) => list("order", children),
            Self::Bond(children) => list("bond", children),
            Self::OneOf(children) => list("one-of", children),
            Self::Remote {
                accession,
                location,
            } => format!(
                "{accession}{}",
                location
                    .as_ref()
                    .map_or(String::new(), |child| format!(":{}", child.expression()))
            ),
            Self::Gap { length, unknown100 } => format!(
                "gap({})",
                if *unknown100 {
                    "unk100".into()
                } else {
                    length.map_or(String::new(), |n| n.to_string())
                }
            ),
        }
    }
    pub fn extract(&self, sequence: &str, molecule: Molecule) -> Result<String, LocationError> {
        self.extract_with(sequence, molecule, |_| None)
    }
    /// Resolve remote accessions only from a caller-supplied record collection.
    /// No network access, accession guessing, or substitution with local bases.
    pub fn extract_with<'a>(
        &self,
        sequence: &'a str,
        molecule: Molecule,
        resolver: impl Fn(&str) -> Option<&'a str>,
    ) -> Result<String, LocationError> {
        fn walk<'a>(
            node: &FeatureLocation,
            sequence: &'a str,
            molecule: Molecule,
            resolver: &impl Fn(&str) -> Option<&'a str>,
            depth: usize,
        ) -> Result<String, LocationError> {
            if depth > 64 {
                return Err(error("Location nesting exceeds limits"));
            }
            let out = match node {
                FeatureLocation::Range { start, end, .. } => sequence
                    .get(*start..*end)
                    .ok_or_else(|| error("Location range is outside its sequence"))?
                    .to_owned(),
                FeatureLocation::Between { .. } => String::new(),
                FeatureLocation::Complement(child) => {
                    if molecule == Molecule::Protein {
                        return Err(error("Protein locations cannot be complemented"));
                    }
                    let value = walk(child, sequence, molecule, resolver, depth + 1)?;
                    if molecule == Molecule::Rna {
                        crate::reverse_complement(&value.replace('U', "T")).replace('T', "U")
                    } else {
                        crate::reverse_complement(&value)
                    }
                }
                FeatureLocation::Join(children) => {
                    let mut out = String::new();
                    for child in children {
                        let value = walk(child, sequence, molecule, resolver, depth + 1)?;
                        if out.len() + value.len() > MAX_BASES {
                            return Err(error("Extracted location exceeds sequence limit"));
                        }
                        out.push_str(&value);
                    }
                    out
                }
                FeatureLocation::Remote {
                    accession,
                    location,
                } => {
                    let remote = resolver(accession).ok_or_else(|| {
                        LocationError(format!(
                            "Remote accession {accession} has not been supplied"
                        ))
                    })?;
                    if let Some(child) = location {
                        child.validate(remote.len(), true)?;
                        walk(child, remote, molecule, resolver, depth + 1)?
                    } else {
                        remote.to_owned()
                    }
                }
                FeatureLocation::Order(_) => {
                    return Err(error(
                        "An order location does not assert a contiguous sequence",
                    ));
                }
                FeatureLocation::Bond(_) => {
                    return Err(error("A bond location does not specify a joined sequence"));
                }
                FeatureLocation::OneOf(_) => {
                    return Err(error(
                        "An alternative location has no uniquely determined sequence",
                    ));
                }
                FeatureLocation::Gap { .. } => return Err(error("Gap sequence is unknown")),
            };
            if out.len() > MAX_BASES {
                return Err(error("Extracted location exceeds sequence limit"));
            }
            Ok(out)
        }
        self.validate(sequence.len(), true)?;
        walk(self, sequence, molecule, &resolver, 0)
    }
    /// Remap finite local coordinates through a replacement. A deleted location
    /// returns None; remote subtrees retain their original coordinate system.
    /// Between-base sites interrupted by an edit require explicit review.
    pub fn replaced(
        &self,
        start: usize,
        end: usize,
        replacement_len: usize,
    ) -> Result<Option<Self>, LocationError> {
        if start > end {
            return Err(error("Reversed replacement coordinates"));
        }
        let delta = replacement_len as i128 - (end - start) as i128;
        let shift = |position: usize| -> Result<usize, LocationError> {
            usize::try_from(position as i128 + delta)
                .map_err(|_| error("Location coordinate overflow"))
        };
        Ok(match self {
            Self::Range {
                start: a,
                end: b,
                before,
                after,
            } => {
                let x = if *a >= end {
                    shift(*a)?
                } else if *a >= start {
                    start
                } else {
                    *a
                };
                let y = if *b <= start {
                    *b
                } else if *b >= end {
                    shift(*b)?
                } else if *b > start {
                    start + replacement_len
                } else {
                    *b
                };
                (x < y).then_some(Self::Range {
                    start: x,
                    end: y,
                    before: *before || (*a >= start && *a < end),
                    after: *after || (*b > start && *b < end),
                })
            }
            Self::Between { left, right } => {
                if (*left < end && *right >= start) || (start == end && *right == start) {
                    return Err(error(
                        "Edit intersects a between-base location; its new site is ambiguous",
                    ));
                }
                Some(Self::Between {
                    left: if *left >= end { shift(*left)? } else { *left },
                    right: if *right >= end {
                        shift(*right)?
                    } else {
                        *right
                    },
                })
            }
            Self::Complement(child) => child
                .replaced(start, end, replacement_len)?
                .map(|child| Self::Complement(Box::new(child))),
            Self::Join(children)
            | Self::Order(children)
            | Self::Bond(children)
            | Self::OneOf(children) => {
                let mut mapped = vec![];
                for child in children {
                    if let Some(child) = child.replaced(start, end, replacement_len)? {
                        mapped.push(child);
                    }
                }
                if mapped.is_empty() {
                    None
                } else {
                    Some(match self {
                        Self::Join(_) => Self::Join(mapped),
                        Self::Order(_) => Self::Order(mapped),
                        Self::Bond(_) => Self::Bond(mapped),
                        _ => Self::OneOf(mapped),
                    })
                }
            }
            Self::Remote { .. } | Self::Gap { .. } => Some(self.clone()),
        })
    }
    /// Re-express each local residue on the reverse-complemented document while
    /// retaining its biological strand and order. Remote records stay untouched.
    pub fn reverse_complemented(&self, length: usize) -> Result<Self, LocationError> {
        fn complement(child: FeatureLocation) -> FeatureLocation {
            if let FeatureLocation::Complement(inner) = child {
                *inner
            } else {
                FeatureLocation::Complement(Box::new(child))
            }
        }
        Ok(match self {
            Self::Range {
                start,
                end,
                before,
                after,
            } => {
                if start >= end || *end > length {
                    return Err(error("Location outside reverse-complemented sequence"));
                }
                complement(Self::Range {
                    start: length - end,
                    end: length - start,
                    before: *after,
                    after: *before,
                })
            }
            Self::Between { left, right } => {
                if *left >= length || *right >= length {
                    return Err(error("Between-base location outside sequence"));
                }
                complement(Self::Between {
                    left: length - 1 - right,
                    right: length - 1 - left,
                })
            }
            Self::Complement(child) => complement(child.reverse_complemented(length)?),
            Self::Join(children) => Self::Join(
                children
                    .iter()
                    .map(|c| c.reverse_complemented(length))
                    .collect::<Result<_, _>>()?,
            ),
            Self::Order(children) => Self::Order(
                children
                    .iter()
                    .map(|c| c.reverse_complemented(length))
                    .collect::<Result<_, _>>()?,
            ),
            Self::Bond(children) => Self::Bond(
                children
                    .iter()
                    .map(|c| c.reverse_complemented(length))
                    .collect::<Result<_, _>>()?,
            ),
            Self::OneOf(children) => Self::OneOf(
                children
                    .iter()
                    .map(|c| c.reverse_complemented(length))
                    .collect::<Result<_, _>>()?,
            ),
            Self::Remote { .. } | Self::Gap { .. } => self.clone(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn range(start: usize, end: usize) -> FeatureLocation {
        FeatureLocation::Range {
            start,
            end,
            before: false,
            after: false,
        }
    }
    #[test]
    fn complement_compounds_keep_biological_order_and_partial_ends() {
        let location = FeatureLocation::Complement(Box::new(FeatureLocation::Join(vec![
            FeatureLocation::Range {
                start: 0,
                end: 3,
                before: true,
                after: false,
            },
            range(8, 11),
        ])));
        assert_eq!(location.expression(), "complement(join(<1..3,9..11))");
        assert_eq!(location.partial_ends(), (false, true));
        assert_eq!(
            location.uniform_projection(),
            Some((
                vec![Segment { start: 0, end: 3 }, Segment { start: 8, end: 11 }],
                true
            ))
        );
        assert_eq!(
            location.extract("AAACCCGGGTTT", Molecule::Dna).unwrap(),
            "AACTTT"
        );
        let mixed = FeatureLocation::Join(vec![
            range(0, 3),
            FeatureLocation::Complement(Box::new(range(3, 6))),
        ]);
        assert!(mixed.uniform_projection().is_none());
        assert_eq!(
            mixed.extract("AAACCCGGGTTT", Molecule::Dna).unwrap(),
            "AAAGGG"
        );
        let reversed = mixed.reverse_complemented(12).unwrap();
        assert_eq!(
            reversed.extract("AAACCCGGGTTT", Molecule::Dna).unwrap(),
            "AAAGGG"
        );
        assert_eq!(reversed.reverse_complemented(12).unwrap(), mixed);
    }
    #[test]
    fn uncertain_operators_and_remote_accessions_never_become_local_joins() {
        for node in [
            FeatureLocation::Order(vec![range(0, 3), range(6, 9)]),
            FeatureLocation::Bond(vec![range(0, 1), range(8, 9)]),
            FeatureLocation::OneOf(vec![range(0, 3), range(6, 9)]),
            FeatureLocation::Gap {
                length: Some(3),
                unknown100: false,
            },
        ] {
            assert!(node.uniform_projection().is_none());
            assert!(node.extract("AAACCCGGGTTT", Molecule::Dna).is_err());
        }
        let remote = FeatureLocation::Remote {
            accession: "SYNTHETIC.1".into(),
            location: Some(Box::new(range(2, 5))),
        };
        let joined = FeatureLocation::Join(vec![range(0, 3), remote.clone()]);
        assert!(joined.extract("AAACCC", Molecule::Dna).is_err());
        assert_eq!(joined.local_parts(), vec![Segment { start: 0, end: 3 }]);
        assert_eq!(
            joined
                .extract_with("AAACCC", Molecule::Dna, |accession| (accession
                    == "SYNTHETIC.1")
                    .then_some("GGAATT"))
                .unwrap(),
            "AAAAAT"
        );
        assert_eq!(remote.replaced(0, 1, 4).unwrap(), Some(remote.clone()));
        assert_eq!(remote.reverse_complemented(6).unwrap(), remote);
    }
    #[test]
    fn insertion_boundaries_are_outside_and_internal_edits_preserve_fuzziness() {
        let location = FeatureLocation::Range {
            start: 2,
            end: 8,
            before: true,
            after: false,
        };
        assert_eq!(location.replaced(8, 8, 3).unwrap(), Some(location.clone()));
        assert_eq!(
            location.replaced(2, 2, 3).unwrap(),
            Some(FeatureLocation::Range {
                start: 5,
                end: 11,
                before: true,
                after: false
            })
        );
        assert_eq!(
            location.replaced(4, 4, 3).unwrap(),
            Some(FeatureLocation::Range {
                start: 2,
                end: 11,
                before: true,
                after: false
            })
        );
        assert_eq!(
            location.replaced(6, 10, 0).unwrap(),
            Some(FeatureLocation::Range {
                start: 2,
                end: 6,
                before: true,
                after: true
            })
        );
        assert_eq!(location.replaced(0, 10, 0).unwrap(), None);
        let origin = FeatureLocation::Between { left: 9, right: 0 };
        assert!(origin.validate(10, true).is_ok());
        assert!(origin.validate(10, false).is_err());
        assert_eq!(
            FeatureLocation::Between { left: 4, right: 5 }
                .extract("AAAAAA", Molecule::Dna)
                .unwrap(),
            ""
        );
        assert!(
            FeatureLocation::Between { left: 4, right: 5 }
                .replaced(5, 5, 2)
                .is_err()
        );
        assert!(
            FeatureLocation::Between { left: 2, right: 7 }
                .validate(10, true)
                .is_err()
        );
    }
}
