//! Share the writer's record classification and duplicate-site suppression with
//! envelope preservation and provenance transport. Color defaults do not depend
//! on this packet order.
use super::{Document, native_primer, record_preservation::Records};
use std::collections::BTreeSet;

#[derive(Clone, Copy)]
pub(super) enum Payload<'a> {
    Preserved(usize),
    Primer(&'a str),
    Feature,
}
pub(super) struct Emission<'a> {
    pub feature_index: usize,
    pub payload: Payload<'a>,
}
impl Emission<'_> {
    pub(super) fn is_primer(&self, records: &Records) -> bool {
        match self.payload {
            Payload::Preserved(index) => records.record(index).kind == 5,
            Payload::Primer(_) => true,
            Payload::Feature => false,
        }
    }
    pub(super) fn model_ids<'a>(&self, native: &'a Document, records: &'a Records) -> &'a [u64] {
        match self.payload {
            Payload::Preserved(index) => &records.record(index).model_ids,
            _ => std::slice::from_ref(&native.features[self.feature_index].id),
        }
    }
}

/// Preserve traversal order and emit a retained multi-site record only once,
/// exactly as the writer did before appearance-loss reporting was introduced.
pub(super) fn plan<'a>(native: &'a Document, records: &Records) -> Vec<Emission<'a>> {
    let mut seen = BTreeSet::new();
    native
        .features
        .iter()
        .enumerate()
        .filter_map(|(feature_index, feature)| {
            let payload = if let Some(index) = records.record_index(feature.id) {
                if !seen.insert(index) {
                    return None;
                }
                Payload::Preserved(index)
            } else if let Some(sequence) = native_primer(feature, native) {
                Payload::Primer(sequence)
            } else {
                Payload::Feature
            };
            Some(Emission {
                feature_index,
                payload,
            })
        })
        .collect()
}

#[cfg(test)]
mod tests;
