//! Per-strand native ColorRange highlighting, using zero-based inclusive ranges.
use super::{Document, Molecule, Node, Result, bad, xml};
use popgene_core::{Segment, metadata::SequenceColor};
fn color(value: &str) -> Option<u32> {
    if value.len() == 7 && value.starts_with('#') {
        return u32::from_str_radix(&value[1..], 16).ok();
    }
    Some(match value {
        "black" => 0x000000,
        "white" => 0xffffff,
        "red" => 0xff0000,
        "green" => 0x008000,
        "blue" => 0x0000ff,
        "cyan" => 0x00ffff,
        "magenta" => 0xff00ff,
        "yellow" => 0xffff00,
        "orange" => 0xffa500,
        "purple" => 0x800080,
        "gray" | "grey" => 0x808080,
        _ => return None,
    })
}
pub(super) fn read(doc: &mut Document, blocks: &[(u8, &[u8])]) -> Result<()> {
    let mut omitted = false;
    for (kind, body) in blocks {
        if *kind != 20 {
            continue;
        }
        let root = xml(body)?;
        if root.name != "StrandColors" {
            return Err(bad("Invalid strand color root"));
        }
        for strand in &root.children {
            let reverse = match strand.name.as_str() {
                "TopStrand" => false,
                "BottomStrand" => true,
                _ => {
                    omitted = true;
                    continue;
                }
            };
            for range in &strand.children {
                if range.name != "ColorRange" {
                    omitted = true;
                    continue;
                }
                let (start, end) = range
                    .attr("range")
                    .split_once("..")
                    .ok_or_else(|| bad("Invalid native color range"))?;
                let start: usize = start
                    .parse()
                    .map_err(|_| bad("Invalid native color coordinate"))?;
                let end: usize = end
                    .parse()
                    .map_err(|_| bad("Invalid native color coordinate"))?;
                let end = end
                    .checked_add(1)
                    .ok_or_else(|| bad("Native color range overflow"))?;
                if start >= end || end > doc.sequence.len() {
                    return Err(bad("Native color range outside sequence"));
                }
                if let Some(rgb) = color(range.attr("colors"))
                    && doc.molecule != Molecule::Protein
                {
                    doc.sequence_colors.push(SequenceColor {
                        range: Segment { start, end },
                        reverse,
                        color_rgb: rgb,
                    });
                } else {
                    omitted = true;
                }
            }
        }
    }
    if omitted {
        doc.warnings.push("Some native color palettes or non-nucleotide colors are not interpreted; their original packets remain retained and metadata-only saves preserve them.".into());
    }
    Ok(())
}
pub(super) fn write(doc: &Document) -> Node {
    let mut root = Node {
        name: "StrandColors".into(),
        ..Default::default()
    };
    for (reverse, key) in [(false, "TopStrand"), (true, "BottomStrand")] {
        let mut strand = Node {
            name: key.into(),
            ..Default::default()
        };
        for color in doc
            .sequence_colors
            .iter()
            .filter(|color| color.reverse == reverse)
        {
            let mut range = Node {
                name: "ColorRange".into(),
                ..Default::default()
            };
            range.attrs.insert(
                "range".into(),
                format!("{}..{}", color.range.start, color.range.end - 1),
            );
            range
                .attrs
                .insert("colors".into(), format!("#{:06x}", color.color_rgb));
            strand.children.push(range);
        }
        if !strand.children.is_empty() {
            root.children.push(strand);
        }
    }
    root
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn custom_color_ranges_roundtrip_and_remap_independently_of_bases() {
        let mut doc = Document::new("colors", "ACGTACGTACGT").unwrap();
        let xml=b"<StrandColors><TopStrand><ColorRange range=\"0..3\" colors=\"magenta\"/></TopStrand><BottomStrand><ColorRange range=\"5..8\" colors=\"#123456\"/></BottomStrand></StrandColors>";
        read(&mut doc, &[(20, xml)]).unwrap();
        let expected = doc.sequence_colors.clone();
        assert_eq!(expected[0].range, Segment { start: 0, end: 4 });
        assert!(expected[1].reverse);
        let native = super::super::write(&doc).unwrap();
        let mut restored = super::super::read("colors.dna", &native).unwrap();
        assert_eq!(restored.sequence_colors, expected);
        restored.replace(2, 2, "AA").unwrap();
        let exported = super::super::write(&restored).unwrap();
        let reread = super::super::read("edited.dna", &exported).unwrap();
        assert_eq!(reread.sequence_colors, restored.sequence_colors);
        assert!(
            crate::save_warnings(&restored, crate::SaveFormat::GenBank)
                .iter()
                .any(|w| w.contains("nucleotide colors"))
        );
    }
}
