//! Original ancillary records. Layout evidence: sgffp's primary format spec.
//! https://merv1n34k.github.io/sgffp/format-spec.html
//! Decoding is bounded and does not reinterpret these records as edited state.
use super::{Document, Node, Result, Source, bad, compression, xml};
use crate::snapgene_preserve::packets;
use std::collections::BTreeMap;

mod bam;
#[cfg(test)]
mod tests;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NativeEmbeddedTrace {
    /// Zero-based occurrence of the source's type-16 container.
    pub index: usize,
    pub name: String,
    /// Recorded container orientation, not an instruction to alter stored calls.
    pub reverse: bool,
    pub byte_len: usize,
    pub properties: Vec<(String, String)>,
}

fn source_packets(source: &Source) -> Result<Vec<(u8, &[u8])>> {
    if source.format != "SnapGene" {
        return Ok(Vec::new());
    }
    if source.bytes.len() as u64 > crate::MAX_FILE {
        return Err(bad("Stored-record source exceeds the file limit"));
    }
    packets(&source.bytes)
}

type TraceContainer<'a> = (bool, &'a [u8], Vec<(String, String)>);
fn trace_container(body: &[u8]) -> Result<TraceContainer<'_>> {
    let flags = body
        .get(..4)
        .ok_or_else(|| bad("Truncated trace container flags"))?;
    let reverse = match u32::from_be_bytes(flags.try_into().unwrap()) {
        0 => false,
        1 => true,
        _ => return Err(bad("Unrecognized trace-container direction flags")),
    };
    let mut trace = None;
    let mut properties = Vec::new();
    let mut fields = Fields::new(&mut properties);
    for (kind, data) in packets(&body[4..])? {
        match kind {
            18 if trace.is_none() => {
                if !data.starts_with(b"\xaeZTR\r\n\x1a\n") {
                    return Err(bad("Embedded trace is not a ZTR record"));
                }
                trace = Some(data);
            }
            18 => return Err(bad("Duplicate ZTR records in one trace container")),
            8 => fields.xml(&xml(data)?, "Trace properties"),
            _ => {
                fields.push(
                    "Uninterpreted container record",
                    &format!("Packet {kind}: {} bytes retained", data.len()),
                );
            }
        }
    }
    fields.finish();
    Ok((
        reverse,
        trace.ok_or_else(|| bad("Trace container has no ZTR record"))?,
        properties,
    ))
}

/// List stored trace containers without allocating decoded chromatogram signals.
pub fn native_embedded_traces(source: &Source) -> Result<Vec<NativeEmbeddedTrace>> {
    let mut result = Vec::new();
    let mut property_bytes = 0usize;
    for (_, body) in source_packets(source)?
        .into_iter()
        .filter(|(kind, _)| *kind == 16)
    {
        if result.len() >= 1000 {
            return Err(bad("More than 1,000 embedded traces"));
        }
        let index = result.len();
        let (reverse, bytes, properties) = trace_container(body)?;
        property_bytes += properties
            .iter()
            .map(|(key, value)| key.len() + value.len())
            .sum::<usize>();
        if property_bytes > 1024 * 1024 {
            return Err(bad(
                "Embedded trace properties exceed the aggregate display limit",
            ));
        }
        result.push(NativeEmbeddedTrace {
            index,
            name: format!("{} — stored trace {}", source.file_name, index + 1),
            reverse,
            byte_len: bytes.len(),
            properties,
        });
    }
    Ok(result)
}

/// Decode one original embedded ZTR as an independent trace document. The caller
/// must require Save As: it has no external trace path associated with the parent.
pub fn native_embedded_trace(source: &Source, index: usize) -> Result<Document> {
    let (_, body) = source_packets(source)?
        .into_iter()
        .filter(|(kind, _)| *kind == 16)
        .nth(index)
        .ok_or_else(|| bad("Unknown embedded trace index"))?;
    let (reverse, bytes, properties) = trace_container(body)?;
    let name = format!("{}-trace-{}.ztr", source.file_name, index + 1);
    let mut document = crate::ztr::read(&name, bytes)?;
    document.warnings.push(format!("Original embedded trace {} from {}. Container orientation: {}. Stored calls, peaks and signal orientation are retained; no reference alignment or reverse complement has been inferred.", index + 1, source.file_name, if reverse { "reverse" } else { "forward" }));
    for (key, value) in properties {
        document.warnings.push(format!("Original {key}: {value}"));
    }
    document.warnings.push("Opened original embedded chromatogram as an independent document. Choose a new save path; the parent document is not changed.".into());
    Ok(document)
}

/// Shared bounded source display. Counts label bytes as well as field values.
pub(super) struct Fields<'a> {
    output: &'a mut Vec<(String, String)>,
    bytes: usize,
    omitted: bool,
}
impl<'a> Fields<'a> {
    pub(super) fn new(output: &'a mut Vec<(String, String)>) -> Self {
        Self {
            bytes: output.iter().map(|(a, b)| a.len() + b.len()).sum(),
            output,
            omitted: false,
        }
    }
    pub(super) fn full(&self) -> bool {
        self.bytes >= 1024 * 1024 || self.output.len() >= 10_000
    }
    pub(super) fn push(&mut self, label: &str, value: &str) {
        let label = shortened(label, 1024, &mut self.omitted);
        let value = shortened(value, 4096, &mut self.omitted);
        if self.full() || self.bytes + label.len() + value.len() > 1024 * 1024 {
            self.omitted = true;
            return;
        }
        self.bytes += label.len() + value.len();
        self.output.push((label, value));
    }
    fn xml(&mut self, node: &Node, path: &str) {
        if self.full() {
            self.omitted = true;
            return;
        }
        for (key, value) in &node.attrs {
            let key = shortened(key, 256, &mut self.omitted);
            self.push(&format!("{path}/@{key}"), value);
        }
        if !node.text.is_empty() {
            self.push(&format!("{path}/text()"), &node.text);
        }
        if node.attrs.is_empty() && node.children.is_empty() && node.text.is_empty() {
            self.push(path, "(empty)");
        }
        for (index, child) in node.children.iter().enumerate() {
            if self.full() {
                self.omitted = true;
                break;
            }
            let name = shortened(&child.name, 256, &mut self.omitted);
            let path = shortened(
                &format!("{path}/{name}[{}]", index + 1),
                1024,
                &mut self.omitted,
            );
            self.xml(child, &path);
        }
    }
    fn json(&mut self, value: &serde_json::Value, path: &str, depth: usize) -> Result<()> {
        if depth > 64 {
            return Err(bad("RNA structure JSON nesting exceeds 64 levels"));
        }
        if self.full() {
            self.omitted = true;
            return Ok(());
        }
        match value {
            serde_json::Value::Object(object) => {
                if object.is_empty() {
                    self.push(path, "{}");
                }
                for (key, value) in object {
                    if self.full() {
                        self.omitted = true;
                        break;
                    }
                    let key = shortened(key, 256, &mut self.omitted);
                    let path = shortened(&format!("{path}/{key}"), 1024, &mut self.omitted);
                    self.json(value, &path, depth + 1)?;
                }
            }
            serde_json::Value::Array(values) => {
                if values.is_empty() {
                    self.push(path, "[]");
                }
                for (index, value) in values.iter().enumerate() {
                    if self.full() {
                        self.omitted = true;
                        break;
                    }
                    self.json(value, &format!("{path}[{index}]"), depth + 1)?;
                }
            }
            serde_json::Value::String(value) => self.push(path, value),
            _ => self.push(path, &value.to_string()),
        }
        Ok(())
    }
    fn finish(&mut self) {
        if self.omitted {
            self.output.push(("Source display limit".into(),"Additional fields or long values remain in the original source. Display is bounded to 1 MiB / 10,000 fields and 4 KiB per value.".into()));
        }
    }
}
fn shortened(value: &str, limit: usize, omitted: &mut bool) -> String {
    if value.len() <= limit {
        return value.to_owned();
    }
    *omitted = true;
    let mut end = limit - " …".len();
    while !value.is_char_boundary(end) {
        end -= 1;
    }
    format!("{} …", &value[..end])
}

fn json_depth(bytes: &[u8]) -> Result<()> {
    let mut depth = 0usize;
    let mut quoted = false;
    let mut escaped = false;
    for byte in bytes {
        if quoted {
            if escaped {
                escaped = false;
            } else if *byte == b'\\' {
                escaped = true;
            } else if *byte == b'"' {
                quoted = false;
            }
        } else {
            match byte {
                b'"' => quoted = true,
                b'[' | b'{' => {
                    depth += 1;
                    if depth > 64 {
                        return Err(bad("RNA structure JSON nesting exceeds 64 levels"));
                    }
                }
                b']' | b'}' => {
                    depth = depth
                        .checked_sub(1)
                        .ok_or_else(|| bad("Unbalanced JSON delimiters"))?;
                }
                _ => {}
            }
        }
    }
    Ok(()) // serde_json validates syntax and delimiter types after this preflight.
}

pub(super) fn append_metadata(bytes: &[u8], output: &mut Vec<(String, String)>) -> Result<()> {
    let packets = packets(bytes)?;
    let mut fields = Fields::new(output);
    let mut expansion: usize = 32 * 1024 * 1024;
    let mut bgzf_blocks = 10_000;
    let mut trace_index = 0;
    let mut links = BTreeMap::new();
    let mut link_error = None;
    if fields.full() {
        fields.omitted = true;
        fields.finish();
        return Ok(());
    }
    let mut link_input_budget = 16 * 1024 * 1024;
    let mut link_text_budget = 1024 * 1024;
    let has_bam = packets.iter().any(|(kind, _)| *kind == 27);
    for (_, body) in packets.iter().filter(|(kind, _)| *kind == 17 && has_bam) {
        if body.len() > link_input_budget {
            link_error = Some("Alignable-sequence link XML exceeds aggregate decoding limit");
            break;
        }
        link_input_budget -= body.len();
        match xml(body) {
            Ok(root) if root.name == "AlignableSequences" => {
                for sequence in root.children.iter().filter(|node| node.name == "Sequence") {
                    let id = sequence.attr("ID");
                    let name = sequence.attr("name");
                    if id.len() + name.len() > link_text_budget {
                        link_error = Some("Alignable-sequence links exceed aggregate text limit");
                        break;
                    }
                    link_text_budget -= id.len() + name.len();
                    if id.is_empty()
                        || id.len() > 254
                        || sequence.attr("name").len() > 4096
                        || links.len() >= 10_000
                        || links
                            .insert(id.to_owned(), sequence.attr("name").to_owned())
                            .is_some()
                    {
                        link_error = Some("Missing or duplicate alignable-sequence ID");
                    }
                }
            }
            _ => link_error = Some("Invalid alignable-sequence XML"),
        }
    }
    for (kind, body) in packets {
        if !matches!(kind, 14 | 16 | 17 | 27 | 28 | 34) {
            continue;
        }
        if fields.full() {
            fields.omitted = true;
            break;
        }
        let result: Result<()> = (|| {
            match kind {
                14 | 17 | 28 => {
                    expansion = expansion.checked_sub(body.len()).ok_or_else(|| {
                        bad("Ancillary source records exceed aggregate decoding limit")
                    })?;
                    let node = xml(body)?;
                    let expected = match kind {
                        14 => "CustomEnzymeSets",
                        17 => "AlignableSequences",
                        _ => "EnzymeVisibilities",
                    };
                    if node.name != expected {
                        return Err(bad("Unexpected ancillary XML root"));
                    }
                    fields.xml(&node, &format!("Original {expected}"));
                }
                16 => {
                    trace_index += 1;
                    let (reverse, data, properties) = trace_container(body)?;
                    fields.push(&format!("Original embedded trace {trace_index}"),&format!("Stored ZTR: {} bytes; container orientation: {}. Open this trace separately to inspect calls, confidence, peaks and signal.",data.len(),if reverse {"reverse"} else {"forward"}));
                    for (key, value) in properties {
                        fields.push(&format!("Trace {trace_index}/{key}"), &value);
                    }
                }
                34 => {
                    if expansion == 0 {
                        return Err(bad("Aggregate ancillary decoding budget is exhausted"));
                    }
                    // A malformed stream may expand to the limit before failing.
                    // Reserve that work up front and refund only verified bytes.
                    let allowance = expansion.min(8 * 1024 * 1024);
                    expansion -= allowance;
                    let data = compression::xdecode(body, allowance)?;
                    expansion += allowance - data.len();
                    json_depth(&data)?;
                    let value: serde_json::Value = serde_json::from_slice(&data)?;
                    if !value.is_object() {
                        return Err(bad("RNA structure JSON is not an object"));
                    }
                    fields.push("Original RNA structure coordinate scope","Source JSON indices and probabilities are displayed verbatim. Index origin, physical pairing and edited-sequence applicability are not inferred.");
                    fields.json(&value, "Original RNA structure", 0)?;
                }
                27 => {
                    if let Some(error) = link_error {
                        return Err(bad(error));
                    }
                    if expansion == 0 {
                        return Err(bad("Aggregate ancillary decoding budget is exhausted"));
                    }
                    bam::metadata(body, &links, &mut fields, &mut expansion, &mut bgzf_blocks)?;
                }
                _ => unreachable!(),
            }
            Ok(())
        })();
        if let Err(error) = result {
            fields.push(
                &format!("Original packet {kind}: decoding unavailable"),
                &format!("{error}. {} original bytes remain retained.", body.len()),
            );
        }
    }
    fields.finish();
    Ok(())
}
