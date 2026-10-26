//! Public ZTR 1.x chromatograms with bounded nested compression decoding.
//! https://staden.sourceforge.net/ztr.html and the BSD Staden reference codec.
use super::{Error, Result};
use popgene_core::{Document, Source, Trace, TraceCallSet};
mod filters;
const MAGIC: &[u8] = b"\xaeZTR\r\n\x1a\n";
fn bad(message: &str) -> Error {
    Error::Message(format!("ZTR: {message}"))
}
fn integer(bytes: &[u8], offset: usize) -> Result<usize> {
    Ok(u32::from_be_bytes(
        bytes
            .get(offset..offset + 4)
            .ok_or_else(|| bad("Truncated chunk integer"))?
            .try_into()
            .unwrap(),
    ) as usize)
}
struct Chunk<'a> {
    offset: usize,
    kind: [u8; 4],
    metadata: &'a [u8],
    data: &'a [u8],
}
fn chunks(bytes: &[u8]) -> Result<Vec<Chunk<'_>>> {
    if bytes.len() < 10 || !bytes.starts_with(MAGIC) || bytes[8] != 1 {
        return Err(bad("Unsupported header or major version"));
    }
    let mut offset = 10;
    let mut result = vec![];
    while offset < bytes.len() {
        let chunk_start = offset;
        let kind = bytes
            .get(offset..offset + 4)
            .ok_or_else(|| bad("Truncated chunk tag"))?
            .try_into()
            .unwrap();
        let length = integer(bytes, offset + 4)?;
        offset += 8;
        let metadata = bytes
            .get(
                offset
                    ..offset
                        .checked_add(length)
                        .ok_or_else(|| bad("Metadata overflow"))?,
            )
            .ok_or_else(|| bad("Truncated metadata"))?;
        offset += length;
        let length = integer(bytes, offset)?;
        offset += 4;
        let data = bytes
            .get(
                offset
                    ..offset
                        .checked_add(length)
                        .ok_or_else(|| bad("Data overflow"))?,
            )
            .ok_or_else(|| bad("Truncated chunk data"))?;
        offset += length;
        result.push(Chunk {
            offset: chunk_start,
            kind,
            metadata,
            data,
        });
        if result.len() > 10000 {
            return Err(bad("More than 10,000 chunks"));
        }
    }
    Ok(result)
}
fn words(data: &[u8], offset: usize) -> Result<Vec<u16>> {
    if data.len() < offset || !(data.len() - offset).is_multiple_of(2) {
        return Err(bad("Unaligned 16-bit signal"));
    }
    Ok(data[offset..]
        .as_chunks::<2>()
        .0
        .iter()
        .map(|p| u16::from_be_bytes(*p))
        .collect())
}
pub fn read(name: &str, bytes: &[u8]) -> Result<Document> {
    let chunks = chunks(bytes)?;
    let mut tables = filters::HuffmanTables::default();
    let mut decoded = 0;
    let mut warnings = filters::DecodeWarnings::default();
    // Shared tables are file-scoped and may occur after their consumers. Their
    // own compression cannot depend on shared tables (no recursive graph).
    for chunk in chunks.iter().filter(|chunk| chunk.kind == *b"HUFF") {
        let (data, flags) = filters::decode_with_warnings(chunk.data, None)?;
        decoded += data.len();
        if decoded > 16 * 1024 * 1024 {
            return Err(bad("HUFF decoded headers exceed 16 MiB"));
        }
        tables.insert(&data)?;
        warnings.logarithmic |= flags.logarithmic;
        warnings.legacy_float |= flags.legacy_float;
    }
    let mut checksum_start = 0;
    for chunk in chunks.iter().filter(|chunk| chunk.kind == *b"CR32") {
        let (data, _) = filters::decode_with_tables(chunk.data, None, &tables)?;
        if !chunk.metadata.is_empty() || data.len() != 5 {
            return Err(bad("Malformed CR32 checksum chunk"));
        }
        decoded += data.len();
        let expected = u32::from_be_bytes(data[1..5].try_into().unwrap());
        let mut checksum = flate2::Crc::new();
        checksum.update(&bytes[checksum_start..chunk.offset]);
        if checksum.sum() != expected {
            return Err(bad("CR32 checksum mismatch"));
        }
        // The next checksum includes this checksum chunk, per the public spec.
        checksum_start = chunk.offset;
    }
    let base_chunks = chunks
        .iter()
        .filter(|c| c.kind == *b"BASE")
        .collect::<Vec<_>>();
    if base_chunks.len() != 1 {
        return Err(bad(
            "Exactly one BASE chunk is required for the active call set",
        ));
    }
    // TSHIFT signal decoding depends on the original, case-sensitive BASE bytes,
    // irrespective of whether the BASE chunk occurs before or after the signal.
    let (base_data, base_warnings) =
        filters::decode_with_tables(base_chunks[0].data, None, &tables)?;
    if base_data.len() > popgene_core::MAX_BASES + 1 {
        return Err(bad("BASE exceeds ten million calls"));
    }
    warnings.logarithmic |= base_warnings.logarithmic;
    warnings.legacy_float |= base_warnings.legacy_float;
    let mut channels = vec![vec![]; 4];
    let mut sequence = None;
    let mut peaks = None;
    let mut confidence = None;
    let mut comments = vec![];
    let mut text_bytes = 0usize;
    let mut text_omitted = false;
    let mut unknown = 0;
    for chunk in &chunks {
        if !matches!(
            &chunk.kind,
            b"SMP4"
                | b"SAMP"
                | b"BASE"
                | b"BPOS"
                | b"CNF4"
                | b"CNF1"
                | b"TEXT"
                | b"COMM"
                | b"CLIP"
                | b"CR32"
                | b"HUFF"
        ) {
            unknown += 1;
            continue;
        }
        if matches!(&chunk.kind, b"CR32" | b"HUFF") {
            continue;
        } // Verified over stored bytes above.
        let (data, chunk_warnings) = if chunk.kind == *b"BASE" {
            (base_data.clone(), base_warnings)
        } else {
            filters::decode_with_tables(chunk.data, Some(&base_data[1..]), &tables)?
        };
        warnings.logarithmic |= chunk_warnings.logarithmic;
        warnings.legacy_float |= chunk_warnings.legacy_float;
        decoded += data.len();
        if decoded > 128 * 1024 * 1024 {
            return Err(bad("Decoded data exceeds 128 MiB"));
        }
        match &chunk.kind {
            b"SMP4" => {
                let signal = words(&data, 2)?;
                if signal.is_empty()
                    || !signal.len().is_multiple_of(4)
                    || signal.len() / 4 > 2_000_000
                {
                    return Err(bad("Invalid four-channel signal length"));
                }
                let length = signal.len() / 4;
                channels = signal.chunks_exact(length).map(|s| s.to_vec()).collect();
            }
            b"SAMP" => {
                let mut fields = chunk.metadata.split(|b| *b == 0);
                let label = std::iter::from_fn(|| Some((fields.next()?, fields.next()?)))
                    .find(|(key, _)| *key == b"TYPE")
                    .and_then(|(_, value)| (value.len() == 1).then(|| value[0]))
                    .or_else(|| {
                        (chunk.metadata.len() <= 4)
                            .then(|| chunk.metadata.first().copied())
                            .flatten()
                    });
                if let Some(channel) = label.and_then(|b| b"ACGT".iter().position(|v| *v == b)) {
                    let signal = words(&data, 2)?;
                    if signal.len() > 2_000_000 {
                        return Err(bad("Signal exceeds two million samples"));
                    }
                    channels[channel] = signal;
                } else {
                    unknown += 1;
                }
            }
            b"BASE" => {
                if sequence.is_some() {
                    return Err(bad(
                        "Multiple BASE chunks require an explicit alternate-call model",
                    ));
                }
                sequence = Some(
                    std::str::from_utf8(&data[1..])
                        .map_err(|_| bad("Invalid base encoding"))?
                        .to_string(),
                );
            }
            b"BPOS" => {
                if peaks.is_some() || data.len() < 4 || !(data.len() - 4).is_multiple_of(4) {
                    return Err(bad("Duplicate or malformed BPOS"));
                }
                peaks = Some(
                    data[4..]
                        .as_chunks::<4>()
                        .0
                        .iter()
                        .map(|p| u32::from_be_bytes(*p) as usize)
                        .collect::<Vec<_>>(),
                );
            }
            b"CNF4" | b"CNF1" => {
                if confidence.is_some() {
                    return Err(bad("Multiple confidence chunks are ambiguous"));
                }
                confidence = Some((chunk.kind, data[1..].to_vec()));
            }
            b"TEXT" => {
                let mut fields = data[1..].split(|b| *b == 0);
                while let Some(key) = fields.next() {
                    if key.is_empty() {
                        break;
                    }
                    let value = fields
                        .next()
                        .ok_or_else(|| bad("TEXT identifier has no value"))?;
                    if comments.len() >= 10_000
                        || key.len() + value.len() > 1024 * 1024
                        || text_bytes.saturating_add(key.len() + value.len()) > 1024 * 1024
                    {
                        text_omitted = true;
                        break;
                    }
                    text_bytes += key.len() + value.len();
                    comments.push(format!(
                        "{}={}",
                        String::from_utf8_lossy(key),
                        String::from_utf8_lossy(value)
                    ));
                }
            }
            b"COMM" => {
                if comments.len() >= 10_000
                    || text_bytes.saturating_add(data.len() - 1) > 1024 * 1024
                {
                    text_omitted = true;
                } else {
                    text_bytes += data.len() - 1;
                    comments.push(String::from_utf8_lossy(&data[1..]).into_owned());
                }
            }
            b"CLIP" => {
                if data.len() != 9 {
                    return Err(bad("Malformed CLIP coordinates"));
                }
                comments.push(format!("Suggested ZTR clipping (stored call coordinates): left {}, right {}. Not applied to active calls.", integer(&data, 1)?, integer(&data, 5)?));
            }
            _ => unreachable!(),
        }
    }
    let sequence = sequence.ok_or_else(|| bad("No base calls"))?;
    let peaks = peaks.ok_or_else(|| bad("No sample positions for base calls"))?;
    let qualities_known = confidence.is_some();
    let qualities = if let Some((kind, confidence)) = confidence {
        let expected = sequence.len() * if kind == *b"CNF4" { 4 } else { 1 };
        if confidence.len() != expected {
            return Err(bad("Confidence length does not match BASE"));
        }
        confidence[..sequence.len()].to_vec()
    } else {
        vec![0; sequence.len()]
    };
    let mut document = Document::new(
        std::path::Path::new(name)
            .file_stem()
            .unwrap_or_default()
            .to_string_lossy(),
        &sequence,
    )?;
    document.trace = Some(Trace {
        calls_version: 1,
        qualities: qualities.clone(),
        qualities_known,
        peaks: peaks.clone(),
        channel_order: "ACGT".into(),
        channels,
        original_calls: vec![TraceCallSet {
            label: "Original ZTR calls".into(),
            version: 1,
            sequence: document.sequence.to_string(),
            qualities,
            qualities_known,
            peaks,
        }],
        edited: false,
        sample_range: None,
        manual_calls: vec![],
    });
    if text_omitted {
        comments.push("Additional ZTR text exceeds the display limit and remains in the original source bytes.".into());
    }
    document.comments = comments;
    document.description = "Imported ZTR chromatogram".into();
    document.source = Some(Source {
        file_name: name.into(),
        format: "ZTR".into(),
        bytes: bytes.into(),
        retained_blocks: chunks.iter().map(|c| (0, c.data.len())).collect(),
    });
    if unknown > 0 {
        document.warnings.push(format!("{unknown} additional ZTR chunks are retained in the original source; only base calls, confidence, peak positions, A/C/G/T signals and text are displayed."));
    }
    if warnings.logarithmic {
        document.warnings.push("ZTR LOG2 compression is lossy: decoded signal samples are quantized approximations of the original amplitudes.".into());
    }
    if warnings.legacy_float {
        document.warnings.push("Legacy ZTR floating Chebyshev73 decoded with the little-endian reference first-word convention, fixed IEEE754 binary32 intermediates and unfused arithmetic. Historical encoders using different floating-point rounding can produce different signal amplitudes; original bytes are retained.".into());
    }
    document.validate()?;
    Ok(document)
}
fn chunk(output: &mut Vec<u8>, kind: &[u8; 4], metadata: &[u8], data: &[u8]) {
    output.extend(kind);
    output.extend((metadata.len() as u32).to_be_bytes());
    output.extend(metadata);
    output.extend((data.len() as u32).to_be_bytes());
    output.extend(data);
}
pub fn write(document: &Document) -> Result<Vec<u8>> {
    document.validate()?;
    let trace = document
        .trace
        .as_ref()
        .ok_or_else(|| bad("No chromatogram"))?;
    let mut output = MAGIC.to_vec();
    output.extend([1, 2]);
    let mut signal = vec![0, 0];
    for base in b"ACGT" {
        let channel = trace
            .channel_order
            .as_bytes()
            .iter()
            .position(|b| b == base)
            .unwrap();
        for &value in &trace.channels[channel] {
            signal.extend(value.to_be_bytes());
        }
    }
    chunk(&mut output, b"SMP4", &[], &signal);
    let mut calls = vec![0];
    calls.extend(document.sequence.as_bytes());
    chunk(&mut output, b"BASE", &[], &calls);
    let mut positions = vec![0; 4];
    for &peak in &trace.peaks {
        positions.extend((peak as u32).to_be_bytes());
    }
    chunk(&mut output, b"BPOS", &[], &positions);
    let mut quality = vec![0];
    quality.extend(&trace.qualities);
    quality.resize(1 + document.sequence.len() * 4, 0);
    if trace.qualities_known {
        chunk(&mut output, b"CNF4", &[], &quality);
    }
    if !document.comments.is_empty() {
        let mut comments = vec![0];
        comments.extend(document.comments.join("\n").as_bytes());
        chunk(&mut output, b"COMM", &[], &comments);
    }
    if output.len() as u64 > super::MAX_FILE {
        return Err(bad("Output exceeds 128 MiB"));
    }
    Ok(output)
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn portable_trace_roundtrip_retains_signal_and_calls() {
        let document = crate::import(
            "x.ab1",
            include_bytes!("../../../validation/fixtures/synthetic-complete.ab1"),
        )
        .unwrap();
        let bytes = write(&document).unwrap();
        let read = read("x.ztr", &bytes).unwrap();
        assert_eq!(read.sequence, document.sequence);
        let a = document.trace.unwrap();
        let b = read.trace.unwrap();
        assert_eq!(a.peaks, b.peaks);
        assert_eq!(a.qualities, b.qualities);
        for base in b"ACGT" {
            let ai = a
                .channel_order
                .as_bytes()
                .iter()
                .position(|b| b == base)
                .unwrap();
            let bi = b
                .channel_order
                .as_bytes()
                .iter()
                .position(|b| b == base)
                .unwrap();
            assert_eq!(a.channels[ai], b.channels[bi]);
        }
    }
    #[test]
    fn chunk_bounds_fail_before_allocations() {
        for bytes in [
            MAGIC.to_vec(),
            [MAGIC, &[1, 2], b"BASE", &u32::MAX.to_be_bytes()].concat(),
        ] {
            assert!(read("bad.ztr", &bytes).is_err());
        }
    }
    #[test]
    fn absent_confidence_stays_unknown_through_import_export_and_alignment() {
        let document = crate::import(
            "x.ab1",
            include_bytes!("../../../validation/fixtures/synthetic-complete.ab1"),
        )
        .unwrap();
        let bytes = write(&document).unwrap();
        let mut missing = bytes[..10].to_vec();
        for part in chunks(&bytes).unwrap() {
            if part.kind != *b"CNF4" {
                chunk(&mut missing, &part.kind, part.metadata, part.data);
            }
        }
        let mut document = read("unknown.ztr", &missing).unwrap();
        assert!(!document.trace.as_ref().unwrap().qualities_known);
        assert!(!crate::SaveFormat::Fastq.supports(&document));
        assert!(!crate::SaveFormat::Scf.supports(&document));
        assert!(popgene_core::trace::quality_trim_range(&document, 20, 2).is_err());
        popgene_core::trace::replace_call(&mut document, 0, "N").unwrap();
        let reopened = read("unknown.ztr", &write(&document).unwrap()).unwrap();
        assert!(!reopened.trace.as_ref().unwrap().qualities_known);
        let aligned = popgene_core::alignment::align_pair(
            &reopened,
            &reopened,
            popgene_core::alignment::AlignmentOptions::default(),
        )
        .unwrap();
        assert!(aligned.query.qualities.iter().all(Option::is_none));
    }
    #[test]
    fn checksums_cover_header_and_previous_checksum_chunk() {
        // Independent slow ANSI CRC-32 oracle (production uses flate2/CRC32fast).
        fn reference_crc(bytes: &[u8]) -> u32 {
            let mut crc = !0_u32;
            for byte in bytes {
                crc ^= u32::from(*byte);
                for _ in 0..8 {
                    crc = (crc >> 1) ^ if crc & 1 == 1 { 0xedb88320 } else { 0 };
                }
            }
            !crc
        }
        assert_eq!(reference_crc(b"123456789"), 0xcbf43926);
        let doc = crate::import(
            "x.ab1",
            include_bytes!("../../../validation/fixtures/synthetic-complete.ab1"),
        )
        .unwrap();
        let mut data = write(&doc).unwrap();
        let first = data.len();
        let mut crc_data = vec![0];
        crc_data.extend(reference_crc(&data).to_be_bytes());
        chunk(&mut data, b"CR32", &[], &crc_data);
        chunk(&mut data, b"COMM", &[], b"\0Checksum test");
        let mut second = vec![0];
        second.extend(reference_crc(&data[first..]).to_be_bytes());
        chunk(&mut data, b"CR32", &[], &second);
        assert!(read("verified.ztr", &data).is_ok());
        let mut damaged = data.clone();
        damaged[9] ^= 1;
        assert!(
            read("bad.ztr", &damaged)
                .unwrap_err()
                .to_string()
                .contains("checksum mismatch")
        );
        let mut damaged = data.clone();
        damaged[first + 20] ^= 1;
        assert!(read("bad.ztr", &damaged).is_err());
        data.pop();
        assert!(read("truncated.ztr", &data).is_err());
    }
    #[test]
    fn signal_type_metadata_is_order_independent_and_empty_values_do_not_panic() {
        let source = crate::import(
            "trace.ab1",
            include_bytes!("../../../validation/fixtures/synthetic-complete.ab1"),
        )
        .unwrap();
        let bytes = write(&source).unwrap();
        let parts = chunks(&bytes).unwrap();
        let make = |empty: bool| {
            let mut result = bytes[..10].to_vec();
            for part in &parts {
                if part.kind != *b"SMP4" {
                    chunk(&mut result, &part.kind, part.metadata, part.data);
                    continue;
                }
                let size = (part.data.len() - 2) / 4;
                for (i, base) in b"ACGT".iter().enumerate() {
                    let mut metadata = b"TEST\0value\0TYPE\0".to_vec();
                    if !empty {
                        metadata.push(*base);
                    }
                    metadata.push(0);
                    let mut data = vec![0, 0];
                    data.extend_from_slice(&part.data[2 + i * size..2 + (i + 1) * size]);
                    chunk(&mut result, b"SAMP", &metadata, &data);
                }
            }
            result
        };
        let ordinary = read("original.ztr", &bytes).unwrap();
        let separated = read("separate.ztr", &make(false)).unwrap();
        assert_eq!(
            separated.trace.as_ref().unwrap().channels,
            ordinary.trace.unwrap().channels
        );
        assert!(read("empty-type.ztr", &make(true)).is_err());
    }
    #[test]
    fn large_text_is_bounded_without_losing_original_source_bytes() {
        let source = crate::import(
            "trace.ab1",
            include_bytes!("../../../validation/fixtures/synthetic-complete.ab1"),
        )
        .unwrap();
        let mut bytes = write(&source).unwrap();
        let mut text = vec![0];
        for _ in 0..20_000 {
            text.extend_from_slice(b"KEY\0VALUE\0");
        }
        text.push(0);
        chunk(&mut bytes, b"TEXT", &[], &text);
        let document = read("many-fields.ztr", &bytes).unwrap();
        assert!(document.comments.len() <= 10_001);
        assert!(
            document
                .comments
                .last()
                .unwrap()
                .contains("original source bytes")
        );
        assert_eq!(&*document.source.unwrap().bytes, bytes);
    }
}
