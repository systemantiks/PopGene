//! Native embedded files: XML manifest and ID-prefixed qCompress data packets.
//! Names are display metadata; callers must choose an explicit output path.
use super::{Result, Source, bad, compression, xml};
use crate::snapgene_preserve::packets;
use std::collections::{BTreeMap, BTreeSet};
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NativeAttachment {
    pub id: u32,
    pub name: String,
    pub size: u64,
    pub modified: Option<String>,
    pub compressed: bool,
}
type AttachmentRecords<'a> = (Vec<NativeAttachment>, BTreeMap<u32, &'a [u8]>);
fn records(source: &Source) -> Result<AttachmentRecords<'_>> {
    if source.format != "SnapGene" {
        return Ok((vec![], BTreeMap::new()));
    }
    if source.bytes.len() as u64 > crate::MAX_FILE {
        return Err(bad("Source exceeds file limit"));
    }
    let mut entries = vec![];
    let mut data = BTreeMap::new();
    let mut ids = BTreeSet::new();
    for (kind, body) in packets(&source.bytes)? {
        if kind != 23 {
            continue;
        }
        let manifest = if body.starts_with(b"<?xml") || body.starts_with(b"<Files") {
            Some(body.to_vec())
        } else if body.get(..4) == Some(&[0, 0, 0, 0]) {
            Some(compression::qdecode(&body[4..], 16 * 1024 * 1024)?)
        } else {
            None
        };
        if let Some(manifest) = manifest {
            let root = xml(&manifest)?;
            if root.name != "Files" {
                return Err(bad("Invalid embedded-file manifest"));
            }
            for file in root.children.iter().filter(|node| node.name == "File") {
                let id = file
                    .attr("id")
                    .parse()
                    .map_err(|_| bad("Invalid embedded-file ID"))?;
                let size = file
                    .attr("size")
                    .parse()
                    .map_err(|_| bad("Invalid embedded-file size"))?;
                if size > crate::MAX_FILE || !ids.insert(id) || entries.len() >= 10_000 {
                    return Err(bad("Invalid duplicate or oversized embedded file"));
                }
                entries.push(NativeAttachment {
                    id,
                    name: file.attr("name").into(),
                    size,
                    modified: file.attrs.get("mtime").cloned(),
                    compressed: match file.attr("compressible") {
                        "1" => true,
                        "0" | "" => false,
                        _ => return Err(bad("Invalid embedded-file compression flag")),
                    },
                });
            }
        } else {
            let id = u32::from_be_bytes(
                body.get(..4)
                    .ok_or_else(|| bad("Truncated embedded-file ID"))?
                    .try_into()
                    .unwrap(),
            );
            if data.insert(id, &body[4..]).is_some() {
                return Err(bad("Duplicate embedded-file payload"));
            }
        }
    }
    for entry in &entries {
        let payload = data
            .get(&entry.id)
            .ok_or_else(|| bad("Embedded-file manifest has no matching payload"))?;
        let size = if entry.compressed {
            u32::from_be_bytes(
                payload
                    .get(..4)
                    .ok_or_else(|| bad("Truncated embedded-file length"))?
                    .try_into()
                    .unwrap(),
            ) as u64
        } else {
            payload.len() as u64
        };
        if size != entry.size {
            return Err(bad("Embedded-file size disagrees with manifest"));
        }
    }
    Ok((entries, data))
}
pub fn native_attachments(source: &Source) -> Result<Vec<NativeAttachment>> {
    Ok(records(source)?.0)
}
pub fn native_attachment_bytes(source: &Source, id: u32) -> Result<Vec<u8>> {
    let (entries, data) = records(source)?;
    let entry = entries
        .iter()
        .find(|entry| entry.id == id)
        .ok_or_else(|| bad("Unknown embedded-file ID"))?;
    if entry.compressed {
        compression::qdecode(data[&id], crate::MAX_FILE as usize)
    } else {
        Ok(data[&id].to_vec())
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    #[test]
    fn observed_attachment_layout_preserves_exact_bytes_and_untrusted_name() {
        let raw = b"Synthetic attachment.\n";
        let mut encoder =
            flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(raw).unwrap();
        let mut payload = 7_u32.to_be_bytes().to_vec();
        payload.extend((raw.len() as u32).to_be_bytes());
        payload.extend(encoder.finish().unwrap());
        let mut bytes = vec![];
        super::super::packet(&mut bytes, 23, &payload).unwrap();
        let manifest = format!(
            "<Files><File id=\"7\" name=\"../../sample.txt\" size=\"{}\" compressible=\"1\"/></Files>",
            raw.len()
        );
        super::super::packet(&mut bytes, 23, manifest.as_bytes()).unwrap();
        let source = Source {
            file_name: "sample.dna".into(),
            format: "SnapGene".into(),
            bytes: bytes.into(),
            retained_blocks: vec![],
        };
        let attachments = native_attachments(&source).unwrap();
        assert_eq!(attachments[0].name, "../../sample.txt");
        assert_eq!(native_attachment_bytes(&source, 7).unwrap(), raw);
        assert!(native_attachment_bytes(&source, 6).is_err());
    }
}

#[cfg(test)]
mod compressed_manifest_tests {
    use super::*;
    use std::io::Write;
    #[test]
    fn observed_id_zero_manifest_and_raw_noncompressible_data() {
        let manifest =
            b"<Files><File id=\"9\" name=\"raw.bin\" size=\"3\" compressible=\"0\"/></Files>";
        let mut encoder =
            flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(manifest).unwrap();
        let mut body = 0_u32.to_be_bytes().to_vec();
        body.extend((manifest.len() as u32).to_be_bytes());
        body.extend(encoder.finish().unwrap());
        let mut bytes = vec![];
        super::super::packet(&mut bytes, 23, &body).unwrap();
        let mut body = 9_u32.to_be_bytes().to_vec();
        body.extend([0, 255, 1]);
        super::super::packet(&mut bytes, 23, &body).unwrap();
        let source = Source {
            file_name: "test.dna".into(),
            format: "SnapGene".into(),
            bytes: bytes.into(),
            retained_blocks: vec![],
        };
        assert_eq!(
            native_attachment_bytes(&source, 9).unwrap(),
            vec![0, 255, 1]
        );
        assert!(!native_attachments(&source).unwrap()[0].compressed);
    }
}
