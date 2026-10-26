//! Reconstruction helpers for native records outside the editable sequence model.
//! Keep sequence-dependent records only while the molecule is identical. Notes
//! and file attachments remain their own source records rather than annotations.
use super::{Error, Result};
use popgene_core::Document;

pub(super) fn packets(bytes: &[u8]) -> Result<Vec<(u8, &[u8])>> {
    let mut result = vec![];
    let mut offset = 0usize;
    while offset < bytes.len() {
        if bytes.len() - offset < 5 {
            return Err(Error::Message("Truncated SnapGene packet".into()));
        }
        let kind = bytes[offset];
        let length = u32::from_be_bytes(bytes[offset + 1..offset + 5].try_into().unwrap()) as usize;
        offset += 5;
        let end = offset
            .checked_add(length)
            .filter(|end| *end <= bytes.len())
            .ok_or_else(|| Error::Message("Truncated SnapGene packet payload".into()))?;
        result.push((kind, &bytes[offset..end]));
        offset = end;
        if result.len() > 10_000 {
            return Err(Error::Message("Too many SnapGene packets".into()));
        }
    }
    Ok(result)
}
pub(super) fn same_molecule(current: &Document, original: &Document) -> bool {
    current.molecule == original.molecule
        && current.sequence == original.sequence
        && current.circular == original.circular
        && current.single_stranded == original.single_stranded
        && current.methylation == original.methylation
        && current.duplex_ends == original.duplex_ends
}
