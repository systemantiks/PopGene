//! Explicit named enzyme panels from original type-14 XML records.
//! Type 0 and space-delimited enzyme names are present in AutoSnapGene's
//! original FLAG DNA/protein fixtures. Other types and active/visibility state
//! are deliberately not inferred. See the pinned reference provenance notice.
use super::{Result, Source, bad, xml};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NativeEnzymePanel {
    pub name: String,
    /// Exact source spelling; callers must validate every member against their
    /// enzyme catalogue before applying this panel, without dropping unknowns.
    pub enzyme_names: Vec<String>,
    /// Zero-based source packet and child indices, for provenance only.
    pub packet_index: usize,
    pub panel_index: usize,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct NativeEnzymePanels {
    pub panels: Vec<NativeEnzymePanel>,
    pub warnings: Vec<String>,
}

/// Read the verified explicit-panel subset, never changing the current view or
/// source bytes. A caller may offer panels as drafts for an explicit user choice.
pub fn native_enzyme_panels(source: &Source) -> Result<NativeEnzymePanels> {
    let mut result = NativeEnzymePanels::default();
    if source.format != "SnapGene" {
        return Ok(result);
    }
    if source.bytes.len() as u64 > crate::MAX_FILE {
        return Err(bad("Source exceeds file limit"));
    }
    let mut xml_bytes = 0usize;
    let mut panel_count = 0usize;
    let mut member_count = 0usize;
    let mut unsupported = 0usize;
    let mut visibility_records = false;
    for (packet_index, (kind, body)) in crate::snapgene_preserve::packets(&source.bytes)?
        .into_iter()
        .enumerate()
    {
        if kind == 28 {
            visibility_records = true;
        }
        if kind != 14 {
            continue;
        }
        xml_bytes += body.len();
        if xml_bytes > 16 * 1024 * 1024 {
            return Err(bad("Stored enzyme panels exceed 16 MiB"));
        }
        let root = xml(body)?;
        if root.name != "CustomEnzymeSets" {
            return Err(bad("Invalid stored enzyme-panel root"));
        }
        for (panel_index, panel) in root.children.iter().enumerate() {
            panel_count += 1;
            if panel_count > 512 {
                return Err(bad("Too many stored enzyme panels"));
            }
            if panel.name != "CustomEnzymeSet" || panel.attr("type") != "0" {
                unsupported += 1;
                continue;
            }
            let name = panel.attr("name");
            if name.trim().is_empty() || name.len() > 4096 {
                return Err(bad("Invalid stored enzyme-panel name"));
            }
            let members = panel
                .attrs
                .get("enzymeNames")
                .ok_or_else(|| bad("Stored enzyme panel has no member list"))?;
            let mut enzyme_names = vec![];
            for enzyme in members.split_ascii_whitespace() {
                member_count += 1;
                if enzyme.len() > 128 || member_count > 16_384 {
                    return Err(bad("Stored enzyme-panel members exceed limits"));
                }
                enzyme_names.push(enzyme.to_owned());
            }
            result.panels.push(NativeEnzymePanel {
                name: name.to_owned(),
                enzyme_names,
                packet_index,
                panel_index,
            });
        }
    }
    if unsupported != 0 {
        result.warnings.push(format!(
            "{unsupported} stored enzyme panel record(s) have an unsupported type or element; they remain in the original source and are not offered for use."
        ));
    }
    if visibility_records {
        result.warnings.push("Stored enzyme visibility overrides are retained in the original source but are not interpreted. A named panel does not establish the originally active display selection.".into());
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn source(xml: &str) -> Source {
        let mut bytes = vec![];
        super::super::packet(&mut bytes, 14, xml.as_bytes()).unwrap();
        Source {
            file_name: "panels.dna".into(),
            format: "SnapGene".into(),
            bytes: bytes.into(),
            retained_blocks: vec![],
        }
    }

    #[test]
    fn verified_explicit_panel_keeps_unknown_names_and_original_bytes() {
        let mut source = source(
            "<CustomEnzymeSets><CustomEnzymeSet type=\"0\" name=\"Golden Gate Enzymes\" enzymeNames=\"BbsI BsaI BsmBI BtgZI SapI\"/><CustomEnzymeSet type=\"0\" name=\"User choice\" enzymeNames=\"EcoRI UnknownEnzyme EcoRI\"/><CustomEnzymeSet type=\"1\" name=\"Other\" enzymeNames=\"BamHI\"/></CustomEnzymeSets>",
        );
        let mut bytes = source.bytes.to_vec();
        super::super::packet(
            &mut bytes,
            28,
            b"<EnzymeVisibilities vals=\"unknown encoding\"/>",
        )
        .unwrap();
        source.bytes = bytes.clone().into();
        let result = native_enzyme_panels(&source).unwrap();
        assert_eq!(result.panels.len(), 2);
        assert_eq!(result.panels[0].name, "Golden Gate Enzymes");
        assert_eq!(
            result.panels[0].enzyme_names,
            ["BbsI", "BsaI", "BsmBI", "BtgZI", "SapI"]
        );
        assert_eq!(
            result.panels[1].enzyme_names,
            ["EcoRI", "UnknownEnzyme", "EcoRI"]
        );
        assert_eq!(result.panels[1].packet_index, 0);
        assert_eq!(result.panels[1].panel_index, 1);
        assert_eq!(result.warnings.len(), 2);
        assert_eq!(source.bytes.as_ref(), bytes.as_slice());
    }

    #[test]
    fn explicit_empty_panel_is_distinct_from_absent_or_unsupported() {
        let empty = native_enzyme_panels(&source("<CustomEnzymeSets><CustomEnzymeSet type=\"0\" name=\"None\" enzymeNames=\"\"/></CustomEnzymeSets>")).unwrap();
        assert_eq!(empty.panels.len(), 1);
        assert!(empty.panels[0].enzyme_names.is_empty());
        let absent = native_enzyme_panels(&source("<CustomEnzymeSets/>")).unwrap();
        assert!(absent.panels.is_empty());
        let unknown = native_enzyme_panels(&source("<CustomEnzymeSets><CustomEnzymeSet type=\"9\" name=\"Future\" enzymeNames=\"EcoRI\"/></CustomEnzymeSets>")).unwrap();
        assert!(unknown.panels.is_empty());
        assert_eq!(unknown.warnings.len(), 1);
    }

    #[test]
    fn malformed_and_excessive_verified_panels_fail_explicitly() {
        for record in [
            "<CustomEnzymeSet type=\"0\" name=\"\" enzymeNames=\"EcoRI\"/>".into(),
            "<CustomEnzymeSet type=\"0\" name=\"Missing\"/>".into(),
            format!(
                "<CustomEnzymeSet type=\"0\" name=\"Long\" enzymeNames=\"{}\"/>",
                "A".repeat(129)
            ),
            "<CustomEnzymeSet type=\"0\" name=\"Many\" enzymeNames=\"\"/>".repeat(513),
            format!(
                "<CustomEnzymeSet type=\"0\" name=\"Many\" enzymeNames=\"{}\"/>",
                "EcoRI ".repeat(16_385)
            ),
        ] {
            assert!(
                native_enzyme_panels(&source(&format!(
                    "<CustomEnzymeSets>{record}</CustomEnzymeSets>"
                )))
                .is_err()
            );
        }
    }
}
