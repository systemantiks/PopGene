//! Structured source metadata shared by interoperable file formats.
use serde::{Deserialize, Serialize};
mod duplex;
pub(crate) use duplex::replaced_features;
pub use duplex::{DuplexAnnotations, DuplexGeometry};

/// Recorded whole-template methylation flags used as analysis assumptions.
/// A clear flag does not establish measured absence, and a set flag does not
/// describe site-specific occupancy or distinguish partial/hemimethylation.
/// `Document::methylation == None` means no current whole-template state is known.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct DnaMethylation {
    pub dam: bool,
    pub dcm: bool,
    pub eco_ki: bool,
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Reference {
    pub title: String,
    pub authors: Option<String>,
    pub consortium: Option<String>,
    pub journal: Option<String>,
    pub pubmed: Option<String>,
    pub remark: Option<String>,
    /// Original reference scope/location text. This is bibliographic provenance,
    /// not an annotation location; sequence edits may require human review.
    pub scope: Option<String>,
    /// Additional bibliographic attributes (DOI, volume, pages, dates, etc.).
    #[serde(default)]
    pub extra: Vec<(String, String)>,
}

/// Custom nucleotide color on one strand. Coordinates use the document's forward
/// reference axis even for the opposite strand; rendering does not change bases.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct SequenceColor {
    pub range: crate::Segment,
    pub reverse: bool,
    pub color_rgb: u32,
}
impl SequenceColor {
    pub(crate) fn replaced(&self, start: usize, end: usize, replacement: usize) -> Option<Self> {
        let delta = replacement as isize - (end - start) as isize;
        let a = if self.range.start >= end {
            self.range.start.checked_add_signed(delta)?
        } else if self.range.start >= start {
            start.checked_add(replacement)?
        } else {
            self.range.start
        };
        let b = if self.range.end <= start {
            self.range.end
        } else if self.range.end >= end {
            self.range.end.checked_add_signed(delta)?
        } else {
            start
        };
        (a < b).then(|| Self {
            range: crate::Segment { start: a, end: b },
            ..self.clone()
        })
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn boundary_insertions_are_outside_and_inner_insertions_inherit_color() {
        let color = SequenceColor {
            range: crate::Segment { start: 2, end: 6 },
            reverse: false,
            color_rgb: 0xff0000,
        };
        assert_eq!(
            color.replaced(2, 2, 2).unwrap().range,
            crate::Segment { start: 4, end: 8 }
        );
        assert_eq!(color.replaced(6, 6, 2).unwrap().range, color.range);
        assert_eq!(
            color.replaced(4, 4, 2).unwrap().range,
            crate::Segment { start: 2, end: 8 }
        );
        assert!(color.replaced(2, 6, 3).is_none());
        assert_eq!(
            color.replaced(1, 4, 2).unwrap().range,
            crate::Segment { start: 3, end: 5 }
        );
    }
    #[test]
    fn document_colors_roundtrip_complement_and_validate_bounds() {
        let mut doc = crate::Document::new("color", "ACGTACGT").unwrap();
        doc.sequence_colors.push(SequenceColor {
            range: crate::Segment { start: 1, end: 4 },
            reverse: false,
            color_rgb: 0xff00ff,
        });
        let original = doc.sequence_colors.clone();
        doc.revcomp().unwrap();
        assert_eq!(
            doc.sequence_colors[0].range,
            crate::Segment { start: 4, end: 7 }
        );
        assert!(doc.sequence_colors[0].reverse);
        doc.revcomp().unwrap();
        assert_eq!(doc.sequence_colors, original);
        doc.sequence_colors[0].color_rgb = 0x1000000;
        assert!(doc.validate().is_err());
    }
}
