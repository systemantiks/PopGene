//! Read-only views of imported annotation display metadata.
//!
//! The terminal-note grammar follows the public GenBank–SnapGene conventions.
//! Original qualifiers are never rewritten, and display direction is separate
//! from the biological strand used to extract and translate a feature.
use crate::Feature;
use std::borrow::Cow;
mod cleavage;
mod segments;
pub use segments::{FeatureDisplaySegment, FeatureSegmentStyle};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FeatureDirection {
    Nondirectional,
    Forward,
    Reverse,
    Bidirectional,
}

fn rgb_color(value: &str) -> Option<u32> {
    if value.len() == 7 && value.starts_with('#') {
        return u32::from_str_radix(&value[1..], 16).ok();
    }
    // Named colours belong to the published primer-note vocabulary.
    Some(match value {
        "black" => 0x000000,
        "red" => 0xff0000,
        "orange" => 0xffa500,
        "green" => 0x008000,
        "blue" => 0x0000ff,
        "purple" => 0x800080,
        "gray" => 0x808080,
        _ => return None,
    })
}

fn oligo(value: &str) -> Option<Cow<'_, str>> {
    if value.is_empty()
        || !value.bytes().any(|b| b.is_ascii_alphabetic())
        || !value.bytes().all(|b| {
            b.is_ascii_whitespace() || b"ACGTRYSWKMBDHVN".contains(&b.to_ascii_uppercase())
        })
    {
        return None;
    }
    if value.bytes().any(|b| b.is_ascii_whitespace()) {
        Some(Cow::Owned(
            value.chars().filter(|c| !c.is_ascii_whitespace()).collect(),
        ))
    } else {
        Some(Cow::Borrowed(value))
    }
}

impl Feature {
    pub(crate) fn reverse_display_direction(&mut self, length: usize) {
        // A whole-document reverse complement reverses the displayed arrows too.
        // Preserve all unrelated text, formatting, and unknown metadata.
        let note_direction = self.display_note("direction").map(str::to_owned);
        self.remap_cleavage(|site| length.checked_sub(site));
        self.reverse_display_segments(length);
        for (key, value) in &mut self.qualifiers {
            if key == "snapgene_readingFrame"
                && let Some(frame) = value.as_deref().and_then(|value| value.parse::<i8>().ok())
                && matches!(frame, -3..=-1 | 1..=3)
            {
                *value = Some((-frame).to_string());
            }
            if key == "snapgene_directionality" {
                match value.as_deref() {
                    Some("1") => *value = Some("2".into()),
                    Some("2") => *value = Some("1".into()),
                    _ => {}
                }
            }
        }
        if let Some(old @ ("RIGHT" | "LEFT")) = note_direction.as_deref()
            && let Some((_, Some(note))) =
                self.qualifiers.iter_mut().rev().find(|(k, _)| k == "note")
        {
            *note = note.replacen(old, if old == "RIGHT" { "LEFT" } else { "RIGHT" }, 1);
        }
    }

    pub fn qualifier(&self, key: &str) -> Option<&str> {
        self.qualifiers
            .iter()
            .find_map(|(k, v)| (k == key).then_some(v.as_deref()).flatten())
    }

    /// Interpret only a complete, recognized terminal note, never arbitrary
    /// prose containing words such as "sequence:" or "color:". Unknown fields
    /// leave the entire note available as ordinary annotation text.
    fn parsed_display_note(&self) -> Option<Vec<(&str, &str)>> {
        let note = self
            .qualifiers
            .iter()
            .rev()
            .find(|(k, _)| k == "note")?
            .1
            .as_deref()?;
        let (note, _) = cleavage::split(note)?;
        let mut fields = Vec::new();
        for item in note.split(';').map(str::trim).filter(|s| !s.is_empty()) {
            if item == "5' phosphorylated" && self.kind == "primer_bind" {
                fields.push(("5' phosphorylated", ""));
                continue;
            }
            let (field, value) = item.split_once(':')?;
            let field = field.trim();
            let value = value.trim();
            let valid = match field {
                "color" => rgb_color(value).is_some() || value == "#------",
                "direction" => matches!(value, "RIGHT" | "LEFT" | "BOTH"),
                "sequence" if self.kind == "primer_bind" => oligo(value).is_some(),
                "added" if self.kind == "primer_bind" => {
                    value.len() == 10
                        && value.bytes().enumerate().all(|(i, b)| {
                            if i == 4 || i == 7 {
                                b == b'-'
                            } else {
                                b.is_ascii_digit()
                            }
                        })
                }
                _ => false,
            };
            if !valid {
                return None;
            }
            // Conflicting duplicates invalidate the whole note, including its
            // other fields; a malformed sequence must not authorize its color.
            if fields
                .iter()
                .any(|(key, previous)| *key == field && *previous != value)
            {
                return None;
            }
            fields.push((field, value));
        }
        (!fields.is_empty()).then_some(fields)
    }
    fn display_note(&self, key: &str) -> Option<&str> {
        self.parsed_display_note()?
            .into_iter()
            .find_map(|(field, value)| (field == key).then_some(value))
    }

    /// The recognized terminal display note, including single-span styles,
    /// primers and cleavage-only notes. Writers keep it after prose/review notes.
    pub fn terminal_display_note(&self) -> Option<&str> {
        if self.segment_metadata().is_none()
            && self.parsed_display_note().is_none()
            && self.display_cleavage_sites().is_empty()
        {
            return None;
        }
        self.qualifiers
            .iter()
            .rev()
            .find(|(key, _)| key == "note")?
            .1
            .as_deref()
    }

    /// Published GenBank display vocabulary. Source qualifiers are untouched.
    /// Independent primer sequences are never inferred from binding coordinates.
    pub fn interchange_display_note(&self) -> Option<String> {
        self.interchange_display_note_at(0)
    }

    /// As above, with the document feature index used when a named segment
    /// requires a concrete color even though no source color was recorded.
    pub fn interchange_display_note_at(&self, index: usize) -> Option<String> {
        if let Some(note) = self.segment_display_note() {
            return Some(note.to_owned());
        }
        if self.kind == "primer_bind" {
            if self.qualifier("primer_sequence").is_none()
                && self.qualifier("snapgene_color").is_none()
                && self.qualifier("snapgene_dateAdded").is_none()
                && let Some(note) = self.terminal_display_note()
            {
                return Some(note.to_owned());
            }
            let mut fields = Vec::new();
            let color = self.display_color_rgb().and_then(|rgb| match rgb {
                0x000000 => Some("black"),
                0xff0000 => Some("red"),
                0xffa500 => Some("orange"),
                0x008000 => Some("green"),
                0x0000ff => Some("blue"),
                0x800080 => Some("purple"),
                0x808080 => Some("gray"),
                _ => None,
            });
            if let Some(color) = color {
                fields.push(format!("color: {color}"));
            }
            if let Some(sequence) = self.primer_sequence() {
                fields.push(format!("sequence: {sequence}"));
            }
            if let Some(date) = self.primer_added_date().filter(|date| {
                date.len() == 10
                    && date.bytes().enumerate().all(|(i, b)| {
                        if i == 4 || i == 7 {
                            b == b'-'
                        } else {
                            b.is_ascii_digit()
                        }
                    })
            }) {
                fields.push(format!("added: {date}"));
            }
            if self.display_note("5' phosphorylated").is_some() {
                fields.push("5' phosphorylated".into());
            }
            return (!fields.is_empty()).then(|| fields.join("; "));
        }
        if let Some(note) = self.terminal_display_note() {
            return Some(note.to_owned());
        }
        if let Some(note) = self.named_segment_export_note(index) {
            return Some(note);
        }
        let mut fields = Vec::new();
        if let Some(color) = self.display_color_rgb() {
            fields.push(format!("color: #{color:06x}"));
        }
        if let Some(direction) = self
            .display_direction()
            .or_else(|| {
                (!fields.is_empty()).then_some(if self.reverse {
                    FeatureDirection::Reverse
                } else {
                    FeatureDirection::Forward
                })
            })
            .and_then(|direction| match direction {
                FeatureDirection::Forward => Some("RIGHT"),
                FeatureDirection::Reverse => Some("LEFT"),
                FeatureDirection::Bidirectional => Some("BOTH"),
                _ => None,
            })
        {
            fields.push(format!("direction: {direction}"));
        }
        (!fields.is_empty()).then(|| fields.join("; "))
    }

    pub fn display_color_rgb(&self) -> Option<u32> {
        self.segment_metadata()
            .and_then(|metadata| metadata.spans.first().and_then(|span| span.color_rgb))
            .or_else(|| {
                self.qualifier("snapgene_color")
                    .and_then(rgb_color)
                    .or_else(|| self.display_note("color").and_then(rgb_color))
            })
    }

    /// Effective annotation color. Uncolored primers use PopGene's fixed purple
    /// display default, so native packet regrouping cannot change their appearance.
    /// This is a view fallback, not a recorded primer base-color value.
    /// Other features use their document index, never a filtered-list index.
    /// Recorded feature/segment colors take precedence and remain unmodified.
    pub fn effective_color_rgb(&self, index: usize) -> u32 {
        self.display_color_rgb().unwrap_or_else(|| {
            if self.kind == "primer_bind" {
                0xb985d0
            } else {
                [0x4cc7aa, 0x65abe1, 0xe2ac64, 0xb985d0, 0x5bc7d2, 0xe38299][index % 6]
            }
        })
    }

    pub fn display_direction(&self) -> Option<FeatureDirection> {
        if self
            .location
            .as_ref()
            .is_some_and(|location| location.uniform_projection().is_none())
        {
            return Some(FeatureDirection::Nondirectional);
        }
        if let Some(value) = self.qualifier("snapgene_directionality") {
            return match value {
                "" | "0" => Some(FeatureDirection::Nondirectional),
                "1" => Some(FeatureDirection::Forward),
                "2" => Some(FeatureDirection::Reverse),
                "3" => Some(FeatureDirection::Bidirectional),
                _ => None,
            };
        }
        if let Some(metadata) = self.segment_metadata() {
            return metadata.direction.or_else(|| {
                (self.qualifier("translation").is_none() && self.qualifier("direction").is_none())
                    .then_some(FeatureDirection::Nondirectional)
            });
        }
        match self.display_note("direction") {
            Some("RIGHT") => Some(FeatureDirection::Forward),
            Some("LEFT") => Some(FeatureDirection::Reverse),
            Some("BOTH") => Some(FeatureDirection::Bidirectional),
            _ if self.kind != "primer_bind"
                && self.display_note("color").is_some()
                && self.qualifier("translation").is_none()
                && self.qualifier("direction").is_none() =>
            {
                Some(FeatureDirection::Nondirectional)
            }
            _ => None,
        }
    }

    /// Full imported oligo, preserving case and removing formatting whitespace.
    /// Never inferred from a binding range, which can omit tails or mismatches.
    pub fn primer_sequence(&self) -> Option<Cow<'_, str>> {
        (self.kind == "primer_bind").then_some(())?;
        self.qualifier("primer_sequence")
            .or_else(|| self.display_note("sequence"))
            .and_then(oligo)
    }

    /// An imported value, with no implied calculation or reaction conditions.
    pub fn primer_melting_temperature(&self) -> Option<&str> {
        (self.kind == "primer_bind" && !self.edited).then_some(())?;
        let value = self.qualifier("binding_meltingTemperature")?;
        value.parse::<f64>().ok()?.is_finite().then_some(value)
    }

    pub fn primer_added_date(&self) -> Option<&str> {
        (self.kind == "primer_bind").then_some(())?;
        self.qualifier("snapgene_dateAdded")
            .or_else(|| self.display_note("added"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Segment;

    fn feature(kind: &str, notes: &[&str]) -> Feature {
        Feature {
            location: None,
            id: 1,
            label: "test".into(),
            kind: kind.into(),
            parts: vec![Segment { start: 0, end: 4 }],
            reverse: false,
            qualifiers: notes
                .iter()
                .map(|v| ("note".into(), Some((*v).into())))
                .collect(),
            edited: false,
        }
    }

    #[test]
    fn terminal_metadata_displays_without_rewriting_original_qualifiers() {
        let f = feature(
            "primer_bind",
            &[
                "A primer description",
                "sequence: aaAC GTNN; color: blue; added: 2026-10-21; 5' phosphorylated",
            ],
        );
        let original = f.clone();
        assert_eq!(f.primer_sequence().as_deref(), Some("aaACGTNN"));
        assert_eq!(f.display_color_rgb(), Some(0x0000ff));
        assert_eq!(f.primer_added_date(), Some("2026-10-21"));
        assert_eq!(f.primer_melting_temperature(), None);
        assert_eq!(f, original);
    }

    #[test]
    fn unrecorded_primer_color_is_index_independent_without_creating_metadata() {
        let primer = feature("primer_bind", &["sequence: ACGT"]);
        let before = primer.clone();
        for index in [0, 1, 2, 3, 5, 6, 99, usize::MAX] {
            assert_eq!(primer.display_color_rgb(), None);
            assert_eq!(primer.effective_color_rgb(index), 0xb985d0);
        }
        assert_eq!(
            primer.interchange_display_note().as_deref(),
            Some("sequence: ACGT")
        );
        assert_eq!(primer, before);
        let annotation = feature("misc_feature", &[]);
        assert_eq!(annotation.effective_color_rgb(0), 0x4cc7aa);
        assert_eq!(annotation.effective_color_rgb(1), 0x65abe1);
    }

    #[test]
    fn recorded_primer_colors_override_the_display_default_without_reinterpretation() {
        for (name, color) in [
            ("black", 0x000000),
            ("red", 0xff0000),
            ("orange", 0xffa500),
            ("green", 0x008000),
            ("blue", 0x0000ff),
            ("purple", 0x800080),
            ("gray", 0x808080),
            ("#123456", 0x123456),
        ] {
            let mut primer = feature("primer_bind", &[]);
            primer
                .qualifiers
                .push(("snapgene_color".into(), Some(name.into())));
            let before = primer.clone();
            assert_eq!(primer.display_color_rgb(), Some(color));
            assert_eq!(primer.effective_color_rgb(0), color);
            assert_eq!(primer.effective_color_rgb(7), color);
            assert_eq!(primer, before);
        }
        for value in ["", "unknown-native-color"] {
            let mut primer = feature("primer_bind", &[]);
            primer
                .qualifiers
                .push(("snapgene_color".into(), Some(value.into())));
            assert_eq!(primer.display_color_rgb(), None);
            assert_eq!(primer.effective_color_rgb(2), 0xb985d0);
            assert_eq!(primer.qualifier("snapgene_color"), Some(value));
        }
    }

    #[test]
    fn ordinary_prose_invalid_and_nonterminal_notes_are_not_metadata() {
        for notes in [
            vec!["Contains sequence: ACGT; color: blue"],
            vec!["sequence: ACGT; experimental: true"],
            vec!["sequence: ACGT; color: invalid"],
            vec!["sequence: ACGT; sequence: AAAA"],
            vec!["sequence: ACGT; color: blue", "An ordinary last note"],
        ] {
            let f = feature("primer_bind", &notes);
            assert_eq!(f.primer_sequence(), None);
        }
        assert_eq!(
            feature("misc_feature", &["sequence: ACGT"]).primer_sequence(),
            None
        );
    }

    #[test]
    fn native_metadata_takes_precedence_and_does_not_change_strand() {
        let mut f = feature("misc_feature", &["color: #00FF00; direction: LEFT"]);
        assert_eq!(f.display_color_rgb(), Some(0x00ff00));
        assert_eq!(f.display_direction(), Some(FeatureDirection::Reverse));
        assert!(!f.reverse);
        f.qualifiers
            .push(("snapgene_directionality".into(), Some("".into())));
        assert_eq!(
            f.display_direction(),
            Some(FeatureDirection::Nondirectional)
        );
        f.qualifiers
            .push(("snapgene_color".into(), Some("#993366".into())));
        assert_eq!(f.display_color_rgb(), Some(0x993366));
        assert_eq!(
            feature("misc_feature", &["color: #a6acb3"]).display_direction(),
            Some(FeatureDirection::Nondirectional)
        );
    }

    #[test]
    fn melting_temperature_is_finite_imported_and_hidden_after_edits() {
        let mut f = feature("primer_bind", &[]);
        f.qualifiers
            .push(("binding_meltingTemperature".into(), Some("56".into())));
        assert_eq!(f.primer_melting_temperature(), Some("56"));
        f.edited = true;
        assert_eq!(f.primer_melting_temperature(), None);
        f.edited = false;
        f.qualifiers[0].1 = Some("NaN".into());
        assert_eq!(f.primer_melting_temperature(), None);
    }

    #[test]
    fn reverse_complement_updates_arrows_and_is_reversible() {
        let mut doc = crate::Document::new("test", "ACGTACGT").unwrap();
        doc.features = vec![feature(
            "misc_feature",
            &["color: #00FF00; direction: RIGHT"],
        )];
        doc.features[0]
            .qualifiers
            .push(("snapgene_directionality".into(), Some("1".into())));
        let original = doc.clone();
        doc.revcomp().unwrap();
        assert_eq!(
            doc.features[0].display_direction(),
            Some(FeatureDirection::Reverse)
        );
        assert_eq!(doc.features[0].display_note("direction"), Some("LEFT"));
        doc.revcomp().unwrap();
        assert_eq!(doc, original);
    }
}
