//! Cleavage marks use interbase coordinates: “after base N” means boundary N.
use super::Feature;
use crate::Segment;
pub(super) fn split(note: &str) -> Option<(&str, Vec<usize>)> {
    let mut found = None;
    for marker in ["Cleavage site after base ", "Cleavage sites after bases "] {
        if let Some(index) = note.rfind(marker) {
            if index > 0 && !note[..index].chars().next_back()?.is_whitespace() {
                return None;
            }
            if found.is_some() {
                return None;
            }
            let tail = &note[index + marker.len()..];
            let mut sites = vec![];
            for value in tail.split(',') {
                let value = value.trim();
                if value.is_empty() || !value.bytes().all(|b| b.is_ascii_digit()) {
                    return None;
                }
                let value: usize = value.parse().ok()?;
                if value > crate::MAX_BASES || sites.contains(&value) {
                    return None;
                }
                sites.push(value);
            }
            if sites.is_empty()
                || sites.len() > 20_000
                || marker == "Cleavage site after base " && sites.len() != 1
            {
                return None;
            }
            found = Some((note[..index].trim_end(), sites));
        }
    }
    found.or(Some((note, vec![])))
}
pub(super) fn suffix(note: &mut String, sites: &[usize]) {
    if sites.is_empty() {
        return;
    }
    if !note.is_empty() {
        note.push(' ');
    }
    note.push_str(if sites.len() == 1 {
        "Cleavage site after base "
    } else {
        "Cleavage sites after bases "
    });
    note.push_str(
        &sites
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join(", "),
    );
}
impl Feature {
    pub fn display_cleavage_sites(&self) -> Vec<usize> {
        let Some(note) = self
            .qualifiers
            .iter()
            .rev()
            .find(|(key, _)| key == "note")
            .and_then(|(_, value)| value.as_deref())
        else {
            return vec![];
        };
        let Some((prefix, sites)) = split(note) else {
            return vec![];
        };
        if !prefix.is_empty()
            && self.segment_metadata().is_none()
            && self.display_note("color").is_none()
            && self.display_note("direction").is_none()
        {
            return vec![];
        }
        if sites.iter().all(|site| {
            self.parts
                .iter()
                .any(|p| p.start <= *site && *site <= p.end)
        }) {
            sites
        } else {
            vec![]
        }
    }
    pub(super) fn remap_cleavage(&mut self, map: impl Fn(usize) -> Option<usize>) {
        let sites = self.display_cleavage_sites();
        if sites.is_empty() {
            return;
        }
        let mapped = sites
            .iter()
            .filter_map(|site| map(*site))
            .collect::<Vec<_>>();
        if let Some((_, Some(note))) = self
            .qualifiers
            .iter_mut()
            .rev()
            .find(|(key, _)| key == "note")
            && let Some((prefix, _)) = split(note)
        {
            let mut text = prefix.to_owned();
            suffix(&mut text, &mapped);
            *note = text;
        }
        if mapped.len() != sites.len() {
            self.qualifiers.push((
                "popgene_cleavage_original".into(),
                Some(
                    sites
                        .iter()
                        .map(usize::to_string)
                        .collect::<Vec<_>>()
                        .join(", "),
                ),
            ));
        }
    }
    pub(super) fn project_cleavage(&mut self, map: &impl Fn(&Segment) -> Vec<Segment>) {
        self.remap_cleavage(|site| {
            let range = if site == 0 {
                Segment { start: 0, end: 1 }
            } else {
                Segment {
                    start: site - 1,
                    end: site,
                }
            };
            let mapped = map(&range);
            (mapped.len() == 1 && mapped[0].end - mapped[0].start == 1).then(|| {
                if site == 0 {
                    mapped[0].start
                } else {
                    mapped[0].end
                }
            })
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Document, FeatureSegmentStyle};
    fn sample() -> Feature {
        Feature{id:1,label:"styles".into(),kind:"misc_feature".into(),parts:vec![Segment{start:10,end:30},Segment{start:50,end:80}],reverse:false,location:None,qualifiers:vec![("note".into(),Some("This forward directional feature has 3 segments: 1: 11 .. 30 / #------ / Line 2: 51 .. 65 / #0000ff / Blue 3: 66 .. 80 / #00cc66 / Green Cleavage sites after bases 30, 65, 80".into())),("snapgene_gap".into(),Some("31 .. 50 / #778899 / connector".into()))],edited:false}
    }
    #[test]
    fn styles_and_cleavage_preserve_biological_bases_and_remap_with_edits() {
        let mut doc = Document::new("styles", &"ACGT".repeat(30)).unwrap();
        doc.features.push(sample());
        let sequence = doc.features[0].extract(&doc.sequence);
        assert_eq!(
            doc.features[0].display_segments()[0].style,
            FeatureSegmentStyle::Line
        );
        assert_eq!(doc.features[0].display_cleavage_sites(), vec![30, 65, 80]);
        assert_eq!(
            doc.features[0].display_gaps()[0].range,
            Segment { start: 30, end: 50 }
        );
        assert_eq!(sequence.len(), 50);
        doc.revcomp().unwrap();
        assert_eq!(doc.features[0].extract(&doc.sequence), sequence);
        assert_eq!(doc.features[0].display_cleavage_sites(), vec![90, 55, 40]);
        assert_eq!(
            doc.features[0].display_gaps()[0].range,
            Segment { start: 70, end: 90 }
        );
        doc.revcomp().unwrap();
        assert_eq!(doc.features[0].display_cleavage_sites(), vec![30, 65, 80]);
        doc.replace(40, 40, "AAA").unwrap();
        assert_eq!(doc.features[0].display_cleavage_sites(), vec![30, 68, 83]);
        assert_eq!(
            doc.features[0].display_gaps()[0].range,
            Segment { start: 30, end: 53 }
        );
        assert_eq!(doc.features[0].extract(&doc.sequence), sequence);
    }
    #[test]
    fn arbitrary_prose_and_invalid_suffixes_are_not_cleavage_metadata() {
        let mut feature = sample();
        for note in [
            "Discuss Cleavage site after base 30",
            "color: #ff0000 Cleavage sites after bases 30, nonsense",
            "color: #ff0000 Cleavage sites after bases 30, 30",
        ] {
            feature.qualifiers[0].1 = Some(note.into());
            assert!(feature.display_cleavage_sites().is_empty());
        }
    }
}
