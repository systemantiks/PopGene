//! Typed INSDC location conversion; never replace remote or fuzzy coordinates by
//! exact local bases. The editable model stores unsupported legacy projections.
use super::{Error, Result};
use gb_io::seq::{After, Before, GapLength, Location};
use popgene_core::FeatureLocation;

pub(super) fn from_genbank(value: &Location) -> Result<FeatureLocation> {
    let coordinate = |n: i64| {
        usize::try_from(n).map_err(|_| Error::Message("Negative location coordinate".into()))
    };
    Ok(match value {
        Location::Range((start, before), (end, after)) => FeatureLocation::Range {
            start: coordinate(*start)?,
            end: coordinate(*end)?,
            before: before.0,
            after: after.0,
        },
        Location::Between(left, right) => FeatureLocation::Between {
            left: coordinate(*left)?,
            right: coordinate(*right)?,
        },
        Location::Complement(child) => FeatureLocation::Complement(Box::new(from_genbank(child)?)),
        Location::Join(children) => {
            FeatureLocation::Join(children.iter().map(from_genbank).collect::<Result<_>>()?)
        }
        Location::Order(children) => {
            FeatureLocation::Order(children.iter().map(from_genbank).collect::<Result<_>>()?)
        }
        Location::Bond(children) => {
            FeatureLocation::Bond(children.iter().map(from_genbank).collect::<Result<_>>()?)
        }
        Location::OneOf(children) => {
            FeatureLocation::OneOf(children.iter().map(from_genbank).collect::<Result<_>>()?)
        }
        Location::External(accession, child) => FeatureLocation::Remote {
            accession: accession.clone(),
            location: child
                .as_ref()
                .map(|c| from_genbank(c).map(Box::new))
                .transpose()?,
        },
        Location::Gap(length) => FeatureLocation::Gap {
            length: if let GapLength::Known(n) = length {
                Some(coordinate(*n)?)
            } else {
                None
            },
            unknown100: matches!(length, GapLength::Unk100),
        },
    })
}
pub(super) fn to_genbank(value: &FeatureLocation) -> Location {
    match value {
        FeatureLocation::Range {
            start,
            end,
            before,
            after,
        } => Location::Range(
            (*start as i64, Before(*before)),
            (*end as i64, After(*after)),
        ),
        FeatureLocation::Between { left, right } => Location::Between(*left as i64, *right as i64),
        FeatureLocation::Complement(child) => Location::Complement(Box::new(to_genbank(child))),
        FeatureLocation::Join(children) => {
            Location::Join(children.iter().map(to_genbank).collect())
        }
        FeatureLocation::Order(children) => {
            Location::Order(children.iter().map(to_genbank).collect())
        }
        FeatureLocation::Bond(children) => {
            Location::Bond(children.iter().map(to_genbank).collect())
        }
        FeatureLocation::OneOf(children) => {
            Location::OneOf(children.iter().map(to_genbank).collect())
        }
        FeatureLocation::Remote {
            accession,
            location,
        } => Location::External(
            accession.clone(),
            location.as_ref().map(|c| Box::new(to_genbank(c))),
        ),
        FeatureLocation::Gap { length, unknown100 } => Location::Gap(if *unknown100 {
            GapLength::Unk100
        } else if let Some(length) = length {
            GapLength::Known(*length as i64)
        } else {
            GapLength::Unknown
        }),
    }
}

#[cfg(test)]
mod tests {
    use crate::{SaveFormat, encode, import};
    fn document(locations: &[&str]) -> Vec<u8> {
        let mut value = String::from(
            "LOCUS       LOCATIONS 24 bp DNA circular UNK 25-OCT-2026\nDEFINITION  Synthetic location tests.\nFEATURES             Location/Qualifiers\n",
        );
        for (index, location) in locations.iter().enumerate() {
            value.push_str(&format!("     misc_feature    {location}\n                     /label=\"location_{index}\"\n"));
        }
        value.push_str("ORIGIN\n        1 aaacccgggtttaaacccgggttt\n//\n");
        value.into_bytes()
    }
    #[test]
    fn fuzzy_remote_order_mixed_strand_and_between_survive_standard_roundtrips() {
        let expressions = [
            "<1..>6",
            "complement(<7..12)",
            "order(1..3,10..12)",
            "join(1..3,complement(10..12))",
            "REMOTE.1:1..3",
            "6^7",
            "24^1",
            "one-of(1..3,10..12)",
            "bond(1,12)",
            "join(1..3,gap(3),10..12)",
        ];
        let original = document(&expressions);
        let mut doc = import("locations.gb", &original).unwrap();
        assert_eq!(encode(&doc, SaveFormat::GenBank).unwrap(), original);
        for (feature, expression) in doc.features.iter().zip(expressions) {
            assert_eq!(feature.location.as_ref().unwrap().expression(), expression);
        }
        assert!(doc.features[4].parts.is_empty());
        assert!(doc.features[5].parts.is_empty());
        doc.description = "Edited metadata".into();
        for format in [SaveFormat::GenBank, SaveFormat::Embl] {
            let bytes = encode(&doc, format).unwrap();
            let readback = import("locations", &bytes).unwrap();
            for (feature, expression) in readback.features.iter().zip(expressions) {
                assert_eq!(feature.location.as_ref().unwrap().expression(), expression);
            }
        }
        assert!(encode(&doc, SaveFormat::SnapGene).is_err());
    }
    #[test]
    fn native_between_base_mapping_does_not_claim_two_residues() {
        let mut doc = import("between.gb", &document(&["6^7", "complement(24^1)"])).unwrap();
        let native = encode(&doc, SaveFormat::SnapGene).unwrap();
        let restored = import("between.dna", &native).unwrap();
        assert!(
            restored
                .features
                .iter()
                .all(|feature| feature.parts.is_empty())
        );
        assert_eq!(restored.features[0].location_label(), "6^7");
        assert_eq!(restored.features[1].location_label(), "complement(24^1)");
        assert!(doc.replace(0, 0, "A").is_err());
        assert_eq!(doc.features[0].location_label(), "6^7");
        doc.features.pop();
        doc.replace(0, 0, "A").unwrap();
        assert_eq!(doc.features[0].location_label(), "7^8");
    }
}
