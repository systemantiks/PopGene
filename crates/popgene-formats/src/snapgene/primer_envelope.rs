//! Keep packet-scoped primer settings separate from individual Primer records.
//! Unknown envelopes can be copied intact while their original record set/order
//! remains intact; rebuilding their scope or ID references is not inferred.
use super::{
    Node, Result, annotation_order, bad, packet, record_preservation::Records, serialize_node, xml,
};
use popgene_core::Document;
use quick_xml::{Reader, events::Event};

pub(super) const REFUSAL: &str = "Native primer-packet settings cannot be safely reassigned after this primer or molecule change. The source has an unknown, repeated, interspersed or invalid envelope. Keep the original primer records and molecule unchanged to retain its complete packet(s); changing this primer set requires a verified envelope encoding. No primer settings have been discarded.";

pub(super) enum Envelope<'a> {
    Original(Vec<&'a [u8]>),
    Rebuild {
        settings: Option<Node>,
        minimum_next_id: u64,
        keep_empty: bool,
    },
}
impl Envelope<'_> {
    pub fn minimum_next_id(&self) -> u64 {
        match self {
            Self::Original(_) => 0,
            Self::Rebuild {
                minimum_next_id, ..
            } => *minimum_next_id,
        }
    }
    pub fn emit(&self, out: &mut Vec<u8>, primers: &str, count: usize, next_id: u64) -> Result<()> {
        match self {
            Self::Original(bodies) => {
                for body in bodies {
                    packet(out, 5, body)?;
                }
            }
            Self::Rebuild {
                settings,
                keep_empty,
                ..
            } => {
                if count != 0 || *keep_empty {
                    let mut body = format!(
                        "<?xml version=\"1.0\" encoding=\"UTF-8\"?><Primers nextValidID=\"{next_id}\">"
                    );
                    if let Some(settings) = settings {
                        serialize_node(settings, &mut body);
                    }
                    body.push_str(primers);
                    body.push_str("</Primers>");
                    packet(out, 5, body.as_bytes())?;
                }
            }
        }
        Ok(())
    }
}

pub(super) fn plan<'a>(
    native: &Document,
    records: &Records,
    source_packets: &[(u8, &'a [u8])],
    same_molecule: bool,
) -> Result<Envelope<'a>> {
    let bodies = source_packets
        .iter()
        .filter_map(|(kind, body)| (*kind == 5).then_some(*body))
        .collect::<Vec<_>>();
    if bodies.is_empty() {
        return Ok(Envelope::Rebuild {
            settings: None,
            minimum_next_id: 0,
            keep_empty: false,
        });
    }
    let roots = bodies
        .iter()
        .map(|body| xml(body))
        .collect::<Result<Vec<_>>>()?;
    let source_count = roots
        .iter()
        .map(|root| {
            root.children
                .iter()
                .filter(|node| node.name == "Primer")
                .count()
        })
        .sum::<usize>();
    let emitted = annotation_order::plan(native, records)
        .iter()
        .filter(|emission| emission.is_primer(records))
        .map(|emission| match emission.payload {
            annotation_order::Payload::Preserved(index) => Some(index),
            _ => None,
        })
        .collect::<Vec<_>>();
    let original = records.source_primer_records();
    if same_molecule
        && original.len() == source_count
        && original.iter().all(Option::is_some)
        && original == emitted
    {
        return Ok(Envelope::Original(bodies));
    }
    // These are the complete observed packet-scoped settings in the actual
    // annotated.dna / synthetic-pcr.dna fixtures. Copy values, do not reinterpret
    // them or infer which new binding sites a reference application will show.
    if roots.len() != 1 || opaque_envelope_xml(bodies[0])? {
        return Err(bad(REFUSAL));
    }
    let root = &roots[0];
    if root.name != "Primers"
        || !root.text.trim().is_empty()
        || root.attrs.keys().any(|key| key != "nextValidID")
    {
        return Err(bad(REFUSAL));
    }
    let minimum_next_id = root.attrs.get("nextValidID").map_or(Ok(0), |value| {
        value.parse::<u64>().map_err(|_| bad(REFUSAL))
    })?;
    let mut settings = None;
    let mut saw_primer = false;
    for child in &root.children {
        if child.name == "Primer" {
            saw_primer = true;
            continue;
        }
        if child.name != "HybridizationParams"
            || saw_primer
            || settings.is_some()
            || !child.children.is_empty()
            || !child.text.trim().is_empty()
            || child.attrs.keys().any(|key| {
                !matches!(
                    key.as_str(),
                    "minContinuousMatchLen"
                        | "allowMismatch"
                        | "minMeltingTemperature"
                        | "showAdditionalFivePrimeMatches"
                        | "minimumFivePrimeAnnealing"
                )
            })
        {
            return Err(bad(REFUSAL));
        }
        settings = Some(child.clone());
    }
    Ok(Envelope::Rebuild {
        settings,
        minimum_next_id,
        keep_empty: true,
    })
}

// Node intentionally omits XML comments and processing instructions. Refuse
// rebuilding an envelope containing those rather than silently erasing them.
// Content inside an individual Primer is governed by record preservation.
fn opaque_envelope_xml(body: &[u8]) -> Result<bool> {
    let mut reader = Reader::from_reader(body);
    reader.config_mut().expand_empty_elements = true;
    let mut depth = 0_usize;
    let mut in_primer = false;
    loop {
        match reader
            .read_event()
            .map_err(|error| bad(error.to_string()))?
        {
            Event::Start(node) => {
                if depth == 1 && node.name().as_ref() == b"Primer" {
                    in_primer = true;
                }
                depth += 1;
            }
            Event::End(_) => {
                depth = depth.saturating_sub(1);
                if depth == 1 {
                    in_primer = false;
                }
            }
            Event::Comment(_) | Event::PI(_) | Event::DocType(_) if !in_primer => return Ok(true),
            Event::Eof => return Ok(false),
            _ => {}
        }
    }
}

#[cfg(test)]
mod tests;
