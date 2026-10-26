//! Retain complete unchanged annotation records when a neighboring annotation
//! changes. Matching requires the original biological model and molecule.
use super::{Document, Node, Result, bad, ranges, xml};
use std::collections::BTreeMap;

pub(super) struct Record {
    pub kind: u8,
    pub node: Node,
    /// Imported binding-site models, in the reader's original distinct-site order.
    pub model_ids: Vec<u64>,
}
#[derive(Default)]
pub(super) struct Records {
    records: Vec<Record>,
    by_feature: BTreeMap<u64, usize>,
    source_primers: Vec<Option<usize>>,
    feature_ids: BTreeMap<u64, u64>,
    primer_ids: BTreeMap<u64, u64>,
    pub next_feature_id: u64,
    pub next_primer_id: u64,
}
impl Records {
    pub fn feature_id(&mut self, model_id: u64) -> Result<u64> {
        allocate(&mut self.feature_ids, &mut self.next_feature_id, model_id)
    }
    pub fn primer_id(&mut self, model_id: u64) -> Result<u64> {
        allocate(&mut self.primer_ids, &mut self.next_primer_id, model_id)
    }
    pub fn contains(&self, feature_id: u64) -> bool {
        self.by_feature.contains_key(&feature_id)
    }
    pub fn record_index(&self, feature_id: u64) -> Option<usize> {
        self.by_feature.get(&feature_id).copied()
    }
    pub fn source_primer_records(&self) -> &[Option<usize>] {
        &self.source_primers
    }
    pub fn record(&self, index: usize) -> &Record {
        &self.records[index]
    }
}
fn allocate(ids: &mut BTreeMap<u64, u64>, next: &mut u64, model_id: u64) -> Result<u64> {
    if let Some(id) = ids.get(&model_id) {
        return Ok(*id);
    }
    let id = *next;
    *next = next
        .checked_add(1)
        .ok_or_else(|| bad("Native annotation identifier exhausted"))?;
    ids.insert(model_id, id);
    Ok(id)
}
pub(super) fn same_annotations(a: &Document, b: &Document) -> bool {
    let same = |a: &[popgene_core::Feature], b: &[popgene_core::Feature]| {
        a.len() == b.len() && a.iter().zip(b).all(|(a, b)| same_feature(a, b))
    };
    same(&a.features, &b.features)
        && same(
            &a.duplex_annotations.features,
            &b.duplex_annotations.features,
        )
}
fn same_feature(a: &popgene_core::Feature, b: &popgene_core::Feature) -> bool {
    let mut a = a.clone();
    let mut b = b.clone();
    a.id = 0;
    b.id = 0;
    a == b
}
pub(super) fn unchanged_records(
    current: &Document,
    original: &Document,
    packets: &[(u8, &[u8])],
) -> Result<Records> {
    if !crate::snapgene_preserve::same_molecule(current, original) {
        return Ok(Records::default());
    }
    let old: BTreeMap<_, _> = original
        .features
        .iter()
        .chain(&original.duplex_annotations.features)
        .map(|feature| (feature.id, feature))
        .collect();
    let now: BTreeMap<_, _> = current
        .features
        .iter()
        .chain(&current.duplex_annotations.features)
        .map(|feature| (feature.id, feature))
        .collect();
    let native_len = original
        .duplex_ends
        .as_ref()
        .map_or(original.sequence.len(), |ends| {
            ends.geometry(original.sequence.len()).extent
        });
    let mut id = 0_u64;
    let mut result = Records::default();
    for (kind, body) in packets {
        if !matches!(kind, 5 | 10) {
            continue;
        }
        let root = xml(body)?;
        let name = if *kind == 10 { "Feature" } else { "Primer" };
        for node in root.children.iter().filter(|node| node.name == name) {
            let count = if *kind == 10 {
                1
            } else {
                let mut seen = Vec::new();
                for site in node
                    .children
                    .iter()
                    .filter(|node| node.name == "BindingSite")
                {
                    let value = (
                        ranges(site.attr("location"), native_len, true, original.circular)?,
                        site.attr("boundStrand") == "1",
                    );
                    if !seen.contains(&value) {
                        seen.push(value);
                    }
                }
                seen.len().max(1)
            };
            let ids = (0..count)
                .map(|_| {
                    id += 1;
                    id
                })
                .collect::<Vec<_>>();
            let unchanged = ids.iter().all(|id| match (old.get(id), now.get(id)) {
                (Some(old), Some(now)) => same_feature(old, now),
                _ => false,
            });
            if let Ok(native_id) = node.attr("recentID").parse::<u64>() {
                let next = native_id
                    .checked_add(1)
                    .ok_or_else(|| bad("Native annotation identifier exhausted"))?;
                if *kind == 10 {
                    result.next_feature_id = result.next_feature_id.max(next);
                    result.feature_ids.insert(ids[0], native_id);
                } else {
                    result.next_primer_id = result.next_primer_id.max(next);
                    if unchanged {
                        for id in &ids {
                            result.primer_ids.insert(*id, native_id);
                        }
                    }
                }
            }
            if *kind == 5 {
                result
                    .source_primers
                    .push(unchanged.then_some(result.records.len()));
            }
            if unchanged {
                let index = result.records.len();
                for id in &ids {
                    result.by_feature.insert(*id, index);
                }
                result.records.push(Record {
                    kind: *kind,
                    node: node.clone(),
                    model_ids: ids,
                });
            }
        }
    }
    if id as usize != old.len() {
        return Err(bad(
            "Native annotation source-record association is inconsistent",
        ));
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn changing_one_feature_keeps_other_feature_and_complete_primer_records() {
        let mut bytes =
            super::super::write(&Document::new("records", "AACCGGTTAACC").unwrap()).unwrap();
        super::super::packet(&mut bytes,10,b"<Features><Feature recentID=\"7\" name=\"preserved\" type=\"misc_feature\" directionality=\"1\" displayPreference=\"custom\"><Segment range=\"1-2\" color=\"#aabbcc\" name=\"domain\" auxiliaryStyle=\"retained\"/></Feature><Feature recentID=\"9\" name=\"change\" type=\"misc_feature\" directionality=\"1\"><Segment range=\"3-4\" color=\"#aabbcc\"/></Feature></Features>").unwrap();
        super::super::packet(&mut bytes,5,b"<Primers><Primer recentID=\"4\" name=\"oligo\" sequence=\"TTAACC\" chemistry=\"source value\"><BindingSite location=\"6-11\" boundStrand=\"0\" meltingTemperature=\"52.3\"><Component hybridizedRange=\"6-11\" bases=\"TTAACC\" customBinding=\"retained\"/></BindingSite></Primer></Primers>").unwrap();
        let mut doc = super::super::read("records.dna", &bytes).unwrap();
        let metadata = super::super::source_metadata(&bytes).unwrap();
        assert!(
            metadata
                .iter()
                .any(|(_, value)| value.contains("displayPreference: custom")
                    && value.contains("auxiliaryStyle: retained"))
        );
        assert!(
            metadata
                .iter()
                .any(|(_, value)| value.contains("chemistry: source value")
                    && value.contains("customBinding: retained"))
        );
        doc.features
            .iter_mut()
            .find(|f| f.label == "change")
            .unwrap()
            .label = "changed name".into();
        doc.annotate(
            "new annotation".into(),
            popgene_core::Segment { start: 4, end: 6 },
            false,
        )
        .unwrap();
        let written = super::super::write(&doc).unwrap();
        let packets = crate::snapgene_preserve::packets(&written).unwrap();
        let features = xml(packets.iter().find(|(kind, _)| *kind == 10).unwrap().1).unwrap();
        let first = features
            .children
            .iter()
            .find(|node| node.attr("name") == "preserved")
            .unwrap();
        assert_eq!(first.attr("displayPreference"), "custom");
        assert_eq!(first.attr("recentID"), "7");
        assert_eq!(features.attr("nextValidID"), "11");
        assert_eq!(
            features
                .children
                .iter()
                .find(|node| node.attr("name") == "new annotation")
                .unwrap()
                .attr("recentID"),
            "10"
        );
        assert_eq!(first.children[0].attr("auxiliaryStyle"), "retained");
        let primers = xml(packets.iter().find(|(kind, _)| *kind == 5).unwrap().1).unwrap();
        assert_eq!(primers.children.len(), 1);
        assert_eq!(primers.children[0].attr("recentID"), "4");
        // No primer record changed, so its original packet is retained exactly;
        // an absent source counter is not invented during a neighboring edit.
        let original_packets = crate::snapgene_preserve::packets(&bytes).unwrap();
        assert_eq!(
            packets.iter().find(|(kind, _)| *kind == 5).unwrap().1,
            original_packets
                .iter()
                .find(|(kind, _)| *kind == 5)
                .unwrap()
                .1
        );
        assert_eq!(primers.children[0].attr("chemistry"), "source value");
        assert_eq!(
            primers.children[0].children[0].children[0].attr("customBinding"),
            "retained"
        );
        let reread = super::super::read("updated.dna", &written).unwrap();
        assert!(
            reread
                .features
                .iter()
                .any(|feature| feature.label == "changed name")
        );
        assert_eq!(
            reread
                .features
                .iter()
                .find(|feature| feature.label == "oligo")
                .unwrap()
                .primer_melting_temperature(),
            Some("52.3")
        );
    }
}
