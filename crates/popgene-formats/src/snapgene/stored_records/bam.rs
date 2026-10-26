//! Read-only BAM display under the SAM/BAM 1.6 specification.
//! https://samtools.github.io/hts-specs/SAMv1.pdf
//! No alignment or reference/read association is inferred from record order.
use super::{Fields, Result, bad};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Read,
};

fn bgzf(bytes: &[u8], budget: &mut usize, blocks: &mut usize) -> Result<(Vec<u8>, bool)> {
    let mut output = Vec::new();
    let mut offset = 0;
    let mut last_empty = false;
    while offset < bytes.len() {
        *blocks = blocks
            .checked_sub(1)
            .ok_or_else(|| bad("More than 10,000 BGZF blocks"))?;
        let remaining = &bytes[offset..];
        if remaining.len() < 18 || remaining[..4] != [31, 139, 8, 4] {
            return Err(bad("Invalid BGZF header"));
        }
        let extra_len = u16::from_le_bytes(remaining[10..12].try_into().unwrap()) as usize;
        let extra = remaining
            .get(12..12 + extra_len)
            .ok_or_else(|| bad("Truncated BGZF extra fields"))?;
        let mut p = 0;
        let mut block_len = None;
        while p < extra.len() {
            let header = extra
                .get(p..p + 4)
                .ok_or_else(|| bad("Truncated BGZF subfield"))?;
            let n = u16::from_le_bytes(header[2..4].try_into().unwrap()) as usize;
            p += 4;
            let value = extra
                .get(p..p + n)
                .ok_or_else(|| bad("Truncated BGZF subfield value"))?;
            if &header[..2] == b"BC" {
                if n != 2 || block_len.is_some() {
                    return Err(bad("Invalid or duplicate BGZF BC field"));
                }
                block_len = Some(u16::from_le_bytes(value.try_into().unwrap()) as usize + 1);
            }
            p += n;
        }
        let n = block_len.ok_or_else(|| bad("Missing BGZF BC field"))?;
        if n < 12 + extra_len + 8 {
            return Err(bad("Invalid BGZF block length"));
        }
        let block = remaining
            .get(..n)
            .ok_or_else(|| bad("Truncated BGZF block"))?;
        let size = u32::from_le_bytes(block[n - 4..].try_into().unwrap()) as usize;
        if size > 65536 || size > *budget {
            return Err(bad("BGZF expansion exceeds limits"));
        }
        let allowance = (*budget).min(65536);
        *budget -= allowance;
        let mut decoded = Vec::with_capacity(size);
        let mut decoder = flate2::bufread::GzDecoder::new(block);
        (&mut decoder)
            .take(allowance as u64 + 1)
            .read_to_end(&mut decoded)?;
        if decoded.len() != size || !decoder.into_inner().is_empty() {
            return Err(bad("Invalid BGZF size, checksum or block boundary"));
        }
        *budget += allowance - size;
        last_empty = decoded.is_empty();
        output.extend(decoded);
        offset += n;
    }
    Ok((output, last_empty))
}

struct Cursor<'a> {
    bytes: &'a [u8],
    offset: usize,
}
impl<'a> Cursor<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, offset: 0 }
    }
    fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        let end = self
            .offset
            .checked_add(n)
            .ok_or_else(|| bad("BAM length overflow"))?;
        let result = self
            .bytes
            .get(self.offset..end)
            .ok_or_else(|| bad("Truncated BAM record"))?;
        self.offset = end;
        Ok(result)
    }
    fn u8(&mut self) -> Result<u8> {
        Ok(self.take(1)?[0])
    }
    fn u16(&mut self) -> Result<u16> {
        Ok(u16::from_le_bytes(self.take(2)?.try_into().unwrap()))
    }
    fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    fn i32(&mut self) -> Result<i32> {
        Ok(i32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    fn cstring(&mut self) -> Result<&'a str> {
        let n = self.bytes[self.offset..]
            .iter()
            .position(|b| *b == 0)
            .ok_or_else(|| bad("Unterminated BAM string"))?;
        let value = text(self.take(n)?)?;
        self.take(1)?;
        Ok(value)
    }
    fn done(&self) -> bool {
        self.offset == self.bytes.len()
    }
}
fn text(bytes: &[u8]) -> Result<&str> {
    std::str::from_utf8(bytes).map_err(|_| bad("Invalid BAM text encoding"))
}
struct Reference {
    name: String,
    length: usize,
    circular: bool,
}
fn lookup_reference(id: i32, refs: &[Reference]) -> Result<Option<&Reference>> {
    if id == -1 {
        Ok(None)
    } else if id >= 0 {
        refs.get(id as usize)
            .map(Some)
            .ok_or_else(|| bad("BAM reference ID is out of bounds"))
    } else {
        Err(bad("Invalid negative BAM reference ID"))
    }
}
fn position(pos: i32, reference: Option<&Reference>) -> Result<()> {
    if pos < -1 || (pos >= 0 && reference.is_none_or(|r| pos as usize >= r.length)) {
        return Err(bad("BAM position is out of reference bounds"));
    }
    Ok(())
}
fn scalar(cursor: &mut Cursor<'_>, kind: u8) -> Result<String> {
    Ok(match kind {
        b'c' => (cursor.u8()? as i8).to_string(),
        b'C' => cursor.u8()?.to_string(),
        b's' => (cursor.u16()? as i16).to_string(),
        b'S' => cursor.u16()?.to_string(),
        b'i' => cursor.i32()?.to_string(),
        b'I' => cursor.u32()?.to_string(),
        b'f' => f32::from_bits(cursor.u32()?).to_string(),
        _ => return Err(bad("Unsupported BAM numeric type")),
    })
}
type AuxFields = (Vec<(String, String)>, Option<Vec<u32>>);
fn auxiliary(cursor: &mut Cursor<'_>) -> Result<AuxFields> {
    let mut result = Vec::new();
    let mut tags = BTreeSet::new();
    let mut long_cigar = None;
    while !cursor.done() {
        if tags.len() >= 10_000 {
            return Err(bad("Too many BAM auxiliary fields"));
        }
        let tag = cursor.take(2)?;
        if !tag[0].is_ascii_alphabetic() || !tag[1].is_ascii_alphanumeric() {
            return Err(bad("Invalid BAM auxiliary tag"));
        }
        let tag = text(tag)?.to_owned();
        if !tags.insert(tag.clone()) {
            return Err(bad("Duplicate BAM auxiliary tag"));
        }
        let kind = cursor.u8()?;
        let value = match kind {
            b'A' => {
                let value = cursor.u8()?;
                if !(33..=126).contains(&value) {
                    return Err(bad("Invalid BAM character value"));
                }
                char::from(value).to_string()
            }
            b'Z' | b'H' => {
                let value = cursor.cstring()?;
                if kind == b'H'
                    && (!value.len().is_multiple_of(2)
                        || !value.bytes().all(|b| b.is_ascii_hexdigit()))
                {
                    return Err(bad("Invalid BAM hexadecimal value"));
                }
                let mut clipped = false;
                super::shortened(value, 4096, &mut clipped)
            }
            b'B' => {
                let subtype = cursor.u8()?;
                let n = cursor.u32()? as usize;
                let size = match subtype {
                    b'c' | b'C' => 1,
                    b's' | b'S' => 2,
                    b'i' | b'I' | b'f' => 4,
                    _ => return Err(bad("Unsupported BAM array type")),
                };
                let data = cursor.take(
                    n.checked_mul(size)
                        .ok_or_else(|| bad("BAM array length overflow"))?,
                )?;
                if tag == "CG" {
                    if subtype != b'I' || n > 1_000_000 {
                        return Err(bad("Invalid or excessive long BAM CIGAR"));
                    }
                    long_cigar = Some(
                        data.as_chunks::<4>()
                            .0
                            .iter()
                            .map(|n| u32::from_le_bytes(*n))
                            .collect(),
                    );
                }
                let mut array = Cursor::new(data);
                let mut text = format!("{}[{n}]:", char::from(subtype));
                for _ in 0..n.min(256) {
                    text.push(' ');
                    text.push_str(&scalar(&mut array, subtype)?);
                }
                if n > 256 {
                    text.push_str(" … (remaining array values retained in source)");
                }
                text
            }
            b'c' | b'C' | b's' | b'S' | b'i' | b'I' | b'f' => scalar(cursor, kind)?,
            _ => return Err(bad("Unsupported BAM auxiliary type")),
        };
        result.push((format!("{tag}:{}", char::from(kind)), value));
    }
    Ok((result, long_cigar))
}

pub(super) fn metadata(
    bytes: &[u8],
    links: &BTreeMap<String, String>,
    fields: &mut Fields<'_>,
    budget: &mut usize,
    blocks: &mut usize,
) -> Result<()> {
    let (decoded, eof) = bgzf(bytes, budget, blocks)?;
    let mut cursor = Cursor::new(&decoded);
    if cursor.take(4)? != b"BAM\x01" {
        return Err(bad("Invalid BAM magic"));
    }
    let header_len = cursor.u32()? as usize;
    if header_len > 8 * 1024 * 1024 {
        return Err(bad("BAM header exceeds 8 MiB"));
    }
    let header = text(cursor.take(header_len)?)?;
    let circular_names: BTreeSet<_> = header
        .lines()
        .filter_map(|line| {
            (line.starts_with("@SQ\t") && line.split('\t').any(|field| field == "TP:circular"))
                .then(|| line.split('\t').find_map(|field| field.strip_prefix("SN:")))
                .flatten()
        })
        .collect();
    let n_ref = cursor.u32()? as usize;
    if n_ref > 10_000 {
        return Err(bad("More than 10,000 BAM references"));
    }
    let mut refs = Vec::new();
    let mut names = BTreeSet::new();
    for _ in 0..n_ref {
        let n = cursor.u32()? as usize;
        if !(2..=4096).contains(&n) {
            return Err(bad("Invalid BAM reference-name length"));
        }
        let raw = cursor.take(n)?;
        if raw[n - 1] != 0 || raw[..n - 1].contains(&0) {
            return Err(bad("Invalid BAM reference name"));
        }
        let name = text(&raw[..n - 1])?.to_owned();
        if !names.insert(name.clone()) {
            return Err(bad("Duplicate BAM reference name"));
        }
        let length = cursor.u32()? as usize;
        if length == 0 || length > i32::MAX as usize {
            return Err(bad("Invalid BAM reference length"));
        }
        refs.push(Reference {
            circular: circular_names.contains(name.as_str()),
            name,
            length,
        });
    }
    // Keep partial display local until every alignment's coordinates and link is
    // validated; malformed records never create apparently valid associations.
    let mut output = Vec::new();
    let mut display = Fields::new(&mut output);
    display.push("Original BAM/SAM header", header.trim_end_matches('\0'));
    if !eof {
        display.push(
            "Original BAM warning",
            "BGZF final empty EOF marker is absent; complete decoded records were still validated.",
        );
    }
    for (index, reference) in refs.iter().enumerate() {
        display.push(
            &format!("Original BAM reference {index}"),
            &format!(
                "{}; declared length {}; circular header {}",
                reference.name, reference.length, reference.circular
            ),
        );
    }
    let mut index = 0;
    while !cursor.done() {
        if index >= 100_000 {
            return Err(bad("More than 100,000 BAM alignments"));
        }
        index += 1;
        let block_len = cursor.u32()? as usize;
        if block_len < 32 {
            return Err(bad("BAM alignment shorter than fixed fields"));
        }
        let mut record = Cursor::new(cursor.take(block_len)?);
        let ref_id = record.i32()?;
        let pos = record.i32()?;
        let name_len = record.u8()? as usize;
        let mapq = record.u8()?;
        let bin = record.u16()?;
        let n_cigar = record.u16()? as usize;
        let flag = record.u16()?;
        let n_seq = record.u32()? as usize;
        if n_seq > popgene_core::MAX_BASES {
            return Err(bad("BAM sequence exceeds supported call limit"));
        }
        let mate_ref_id = record.i32()?;
        let mate_pos = record.i32()?;
        let template_len = record.i32()?;
        let reference = lookup_reference(ref_id, &refs)?;
        let mate_reference = lookup_reference(mate_ref_id, &refs)?;
        position(pos, reference)?;
        position(mate_pos, mate_reference)?;
        if flag & 4 == 0 && (reference.is_none() || pos < 0) {
            return Err(bad("Mapped BAM read lacks a valid position"));
        }
        if name_len < 2 {
            return Err(bad("Invalid BAM read-name length"));
        }
        let raw_name = record.take(name_len)?;
        if raw_name[name_len - 1] != 0 || raw_name[..name_len - 1].contains(&0) {
            return Err(bad("Invalid BAM read name"));
        }
        let name = text(&raw_name[..name_len - 1])?;
        let source_name = links.get(name).ok_or_else(|| {
            bad(format!(
                "Stored BAM read {name:?} has no matching alignable Sequence/@ID"
            ))
        })?;
        let mut cigar = Vec::with_capacity(n_cigar);
        for _ in 0..n_cigar {
            cigar.push(record.u32()?);
        }
        let packed = record.take(n_seq.div_ceil(2))?;
        let qualities = record.take(n_seq)?;
        if !qualities.iter().all(|q| *q <= 93) && !qualities.iter().all(|q| *q == 255) {
            return Err(bad("Invalid or mixed missing BAM base qualities"));
        }
        let (aux, long_cigar) = auxiliary(&mut record)?;
        let mut placeholder_span = None;
        if let Some(long) = long_cigar {
            if cigar.len() != 2 || cigar[0] != ((n_seq as u32) << 4 | 4) || cigar[1] & 15 != 3 {
                return Err(bad("Long BAM CIGAR lacks its required placeholder"));
            }
            placeholder_span = Some((cigar[1] >> 4) as usize);
            cigar = long;
        }
        let mut query_len = 0usize;
        let mut ref_len = 0usize;
        let mut cigar_text = String::new();
        for (i, word) in cigar.iter().enumerate() {
            let op = (word & 15) as usize;
            let len = (word >> 4) as usize;
            if op > 8 || len == 0 {
                return Err(bad("Invalid BAM CIGAR operator/length"));
            }
            if op == 5 && i != 0 && i + 1 != cigar.len() {
                return Err(bad("BAM hard clipping is not terminal"));
            }
            if op == 4
                && !(i == 0
                    || i + 1 == cigar.len()
                    || (i == 1 && cigar[0] & 15 == 5)
                    || (i + 2 == cigar.len() && cigar[cigar.len() - 1] & 15 == 5))
            {
                return Err(bad("BAM soft clipping is not terminal"));
            }
            if matches!(op, 0 | 1 | 4 | 7 | 8) {
                query_len = query_len
                    .checked_add(len)
                    .ok_or_else(|| bad("BAM query length overflow"))?;
            }
            if matches!(op, 0 | 2 | 3 | 7 | 8) {
                ref_len = ref_len
                    .checked_add(len)
                    .ok_or_else(|| bad("BAM reference length overflow"))?;
            }
            if i < 256 {
                cigar_text.push_str(&format!("{len}{}", char::from(b"MIDNSHP=X"[op])));
            }
        }
        if cigar.len() > 256 {
            cigar_text.push_str(" … (remaining CIGAR retained in source)");
        }
        if !cigar.is_empty() && n_seq != 0 && query_len != n_seq {
            return Err(bad("BAM CIGAR query length disagrees with SEQ"));
        }
        if placeholder_span.is_some_and(|span| span != ref_len) {
            return Err(bad(
                "Long BAM CIGAR reference span disagrees with placeholder",
            ));
        }
        if let Some(r) = reference
            && flag & 4 == 0
            && pos >= 0
            && !r.circular
            && (pos as usize)
                .checked_add(ref_len)
                .is_none_or(|end| end > r.length)
        {
            return Err(bad("BAM CIGAR exceeds reference bounds"));
        }
        let label = format!("Original BAM alignment {index} (Sequence/@ID={name})");
        display.push(&label,&format!("Source name: {source_name}\nReference: {} (ID {ref_id}); leftmost 1-based position: {}\nFlags: 0x{flag:04x}; reverse-stored sequence: {}; unmapped: {}\nMAPQ: {}; bin: {bin}; stored bases: {n_seq}; template length: {template_len}\nMate reference ID: {mate_ref_id}; mate 1-based position: {}",reference.map_or("unavailable",|r|r.name.as_str()),if pos<0 {"unavailable".into()} else {(pos as usize+1).to_string()},flag&16!=0,flag&4!=0,if mapq==255 {"unavailable".into()} else {mapq.to_string()},if mate_pos<0 {"unavailable".into()} else {(mate_pos as usize+1).to_string()}));
        display.push(
            &format!("{label}/CIGAR"),
            if cigar_text.is_empty() {
                "unavailable"
            } else {
                &cigar_text
            },
        );
        if flag & 4 != 0 {
            display.push(&format!("{label}/Coordinate scope"),"Read is marked unmapped. Stored position and CIGAR are provenance, not an asserted reference alignment.");
        }
        let mut sequence = String::with_capacity(n_seq.min(4096));
        for i in 0..n_seq.min(4096) {
            let code = if i % 2 == 0 {
                packed[i / 2] >> 4
            } else {
                packed[i / 2] & 15
            };
            sequence.push(char::from(b"=ACMGRSVTWYHKDBN"[code as usize]));
        }
        if n_seq > 4096 {
            sequence.push_str(" … (remaining bases retained in source)");
        }
        display.push(&format!("{label}/Stored sequence"), &sequence);
        let quality_text = if qualities.is_empty() {
            "no stored sequence".into()
        } else if qualities[0] == 255 {
            "unavailable".into()
        } else {
            let mut q = qualities
                .iter()
                .take(256)
                .map(u8::to_string)
                .collect::<Vec<_>>()
                .join(" ");
            if n_seq > 256 {
                q.push_str(" … (remaining scores retained in source)");
            }
            q
        };
        display.push(&format!("{label}/Stored Phred scores"), &quality_text);
        for (tag, value) in aux {
            display.push(&format!("{label}/{tag}"), &value);
        }
    }
    display.finish();
    for (label, value) in output {
        fields.push(&label, &value);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn gzip(data: &[u8]) -> Vec<u8> {
        let mut encoder = flate2::GzBuilder::new()
            .extra(vec![b'B', b'C', 2, 0, 0, 0])
            .write(Vec::new(), flate2::Compression::default());
        encoder.write_all(data).unwrap();
        let mut bytes = encoder.finish().unwrap();
        let bsize = (bytes.len() - 1) as u16;
        bytes[16..18].copy_from_slice(&bsize.to_le_bytes());
        bytes
    }
    fn compressed(raw: &[u8]) -> Vec<u8> {
        // Deliberately split fixed fields, CIGAR and sequence across BGZF blocks.
        let mut bytes = Vec::new();
        for chunk in raw.chunks(23) {
            bytes.extend(gzip(chunk));
        }
        bytes.extend(gzip(&[]));
        bytes
    }
    fn bam(pos: i32, ref_id: i32, long: bool) -> Vec<u8> {
        let header = b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:10\n";
        let mut bytes = b"BAM\x01".to_vec();
        bytes.extend((header.len() as u32).to_le_bytes());
        bytes.extend(header);
        bytes.extend(1_u32.to_le_bytes());
        bytes.extend(4_u32.to_le_bytes());
        bytes.extend(b"ref\0");
        bytes.extend(10_u32.to_le_bytes());
        let mut record = Vec::new();
        record.extend(ref_id.to_le_bytes());
        record.extend(pos.to_le_bytes());
        record.extend([2, 255]);
        record.extend(4681_u16.to_le_bytes());
        record.extend((if long { 2_u16 } else { 3_u16 }).to_le_bytes());
        record.extend(16_u16.to_le_bytes());
        record.extend(4_u32.to_le_bytes());
        record.extend((-1_i32).to_le_bytes());
        record.extend((-1_i32).to_le_bytes());
        record.extend(0_i32.to_le_bytes());
        record.extend(b"7\0");
        for word in if long {
            vec![(4_u32 << 4) | 4, (2 << 4) | 3]
        } else {
            vec![(1 << 4) | 4, 2 << 4, (1 << 4) | 1]
        } {
            record.extend(word.to_le_bytes());
        }
        record.extend([0x12, 0x48]);
        record.extend([40, 30, 20, 10]);
        record.extend(b"NMi");
        record.extend(1_i32.to_le_bytes());
        record.extend(b"RGZgroup\0");
        record.extend(b"XABc");
        record.extend(2_u32.to_le_bytes());
        record.extend([1, 254]);
        if long {
            record.extend(b"CGBI");
            record.extend(3_u32.to_le_bytes());
            for word in [(1_u32 << 4) | 4, 2 << 4, (1 << 4) | 1] {
                record.extend(word.to_le_bytes());
            }
        }
        bytes.extend((record.len() as u32).to_le_bytes());
        bytes.extend(record);
        bytes
    }
    fn read(raw: &[u8], links: &BTreeMap<String, String>) -> Result<Vec<(String, String)>> {
        let mut output = Vec::new();
        metadata(
            &compressed(raw),
            links,
            &mut Fields::new(&mut output),
            &mut (1024 * 1024),
            &mut 10_000,
        )?;
        Ok(output)
    }
    #[test]
    fn independent_bam_layout_retains_cigar_coordinates_strand_and_quality() {
        let links = BTreeMap::from([("7".into(), "Synthetic read".into())]);
        for long in [false, true] {
            let fields = read(&bam(1, 0, long), &links).unwrap();
            assert!(
                fields
                    .iter()
                    .any(|(_, v)| v.contains("leftmost 1-based position: 2")
                        && v.contains("reverse-stored sequence: true")
                        && v.contains("MAPQ: unavailable"))
            );
            for (suffix, value) in [
                ("/CIGAR", "1S2M1I"),
                ("/Stored sequence", "ACGT"),
                ("/Stored Phred scores", "40 30 20 10"),
                ("/NM:i", "1"),
                ("/RG:Z", "group"),
                ("/XA:B", "c[2]: 1 -2"),
            ] {
                assert!(
                    fields
                        .iter()
                        .any(|(k, v)| k.ends_with(suffix) && v == value),
                    "missing {suffix}"
                );
            }
            assert!(
                !fields
                    .iter()
                    .any(|(_, v)| v.contains("EOF marker is absent"))
            );
        }
    }
    #[test]
    fn native_metadata_requires_an_exact_unique_alignment_id_link() {
        for (xml, valid) in [
            (
                "<AlignableSequences><Sequence ID=\"7\" name=\"Synthetic read\"/></AlignableSequences>",
                true,
            ),
            (
                "<AlignableSequences><Sequence ID=\"8\" name=\"Wrong read\"/></AlignableSequences>",
                false,
            ),
            (
                "<AlignableSequences><Sequence ID=\"7\"/><Sequence ID=\"7\"/></AlignableSequences>",
                false,
            ),
        ] {
            let doc = popgene_core::Document::new("native alignment", "ACGTACGTAC").unwrap();
            let mut bytes = crate::encode(&doc, crate::SaveFormat::SnapGene).unwrap();
            crate::snapgene::packet(&mut bytes, 17, xml.as_bytes()).unwrap();
            crate::snapgene::packet(&mut bytes, 27, &compressed(&bam(1, 0, false))).unwrap();
            let doc = crate::import("alignment.dna", &bytes).unwrap();
            let fields = crate::source_metadata(doc.source.as_ref().unwrap()).unwrap();
            assert_eq!(
                fields
                    .iter()
                    .any(|(key, value)| key.ends_with("/CIGAR") && value == "1S2M1I"),
                valid
            );
            assert_eq!(
                fields
                    .iter()
                    .any(|(key, _)| key.contains("packet 27: decoding unavailable")),
                !valid
            );
        }
    }
    #[test]
    fn missing_links_coordinates_and_bgzf_corruption_are_rejected() {
        let links = BTreeMap::from([("7".into(), "Synthetic read".into())]);
        assert!(read(&bam(1, 0, false), &BTreeMap::new()).is_err());
        for (pos, reference) in [(10, 0), (9, 0), (-2, 0), (1, 1), (1, -2)] {
            assert!(read(&bam(pos, reference, false), &links).is_err());
        }
        let mut unmapped = bam(9, 0, false);
        let header_len = u32::from_le_bytes(unmapped[4..8].try_into().unwrap()) as usize;
        let record_offset = 8 + header_len + 4 + 4 + 4 + 4 + 4;
        unmapped[record_offset + 14..record_offset + 16].copy_from_slice(&20_u16.to_le_bytes());
        let fields = read(&unmapped, &links).unwrap();
        assert!(
            fields
                .iter()
                .any(|(_, value)| value.contains("unmapped: true"))
        );
        assert!(
            fields
                .iter()
                .any(|(key, value)| key.ends_with("/CIGAR") && value == "1S2M1I")
        );
        assert!(
            fields
                .iter()
                .any(|(_, value)| value.contains("not an asserted reference alignment"))
        );
        let raw = bam(1, 0, false);
        let mut encoded = compressed(&raw);
        encoded[18] ^= 0xff;
        let mut failed_budget = 1024 * 1024;
        assert!(bgzf(&encoded, &mut failed_budget, &mut 10_000).is_err());
        assert!(failed_budget < 1024 * 1024);
        assert!(bgzf(&compressed(&raw), &mut 10, &mut 10_000).is_err());
        let mut encoded = compressed(&raw);
        encoded.pop();
        assert!(bgzf(&encoded, &mut (1024 * 1024), &mut 10_000).is_err());
        let empty = gzip(&[]).repeat(10_001);
        assert!(
            bgzf(&empty, &mut (1024 * 1024), &mut 10_000)
                .unwrap_err()
                .to_string()
                .contains("10,000 BGZF blocks")
        );
        let mut shared_blocks = 1;
        assert!(bgzf(&gzip(&[]), &mut 1024, &mut shared_blocks).is_ok());
        assert!(bgzf(&gzip(&[]), &mut 1024, &mut shared_blocks).is_err());
    }
    #[test]
    fn malformed_bam_variable_fields_do_not_allocate_or_create_partial_rows() {
        let links = BTreeMap::from([("7".into(), "Synthetic read".into())]);
        let original = bam(1, 0, false);
        let header_len = u32::from_le_bytes(original[4..8].try_into().unwrap()) as usize;
        let record_offset = 8 + header_len + 4 + 4 + 4 + 4 + 4;
        for (relative, replace) in [(16, u32::MAX.to_le_bytes()), (12, 0_u32.to_le_bytes())] {
            let mut invalid = original.clone();
            invalid[record_offset + relative..record_offset + relative + 4]
                .copy_from_slice(&replace);
            assert!(read(&invalid, &links).is_err());
        }
        for n in [0, 3, 12, original.len() - 1] {
            assert!(read(&original[..n], &links).is_err());
        }
        // Query length is still four, but interior soft clipping is not a valid
        // SAM/BAM CIGAR. Reject it rather than presenting a valid association.
        let mut interior = original.clone();
        for (index, word) in [1_u32 << 4, (1 << 4) | 4, 2 << 4].iter().enumerate() {
            let start = record_offset + 32 + 2 + index * 4;
            interior[start..start + 4].copy_from_slice(&word.to_le_bytes());
        }
        assert!(read(&interior, &links).is_err());
        let mut long = bam(1, 0, true);
        let placeholder = record_offset + 32 + 2 + 4;
        long[placeholder..placeholder + 4].copy_from_slice(&((3_u32 << 4) | 3).to_le_bytes());
        assert!(read(&long, &links).is_err());
    }
}
