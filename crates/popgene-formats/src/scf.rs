//! SCF 2.x/3.x per the public Staden specification. V3 sample channels use
//! two rounds of modular differencing; V2 samples/base records are interleaved.
//! https://staden.sourceforge.net/manual/formats_unix_2.html
use super::{Error, Result};
use popgene_core::{Document, Source, Trace, TraceCallSet};
fn bad(message: &str) -> Error {
    Error::Message(format!("SCF: {message}"))
}
fn integer(bytes: &[u8], offset: usize) -> Result<usize> {
    Ok(u32::from_be_bytes(
        bytes
            .get(offset..offset + 4)
            .ok_or_else(|| bad("Truncated header"))?
            .try_into()
            .unwrap(),
    ) as usize)
}
fn section(bytes: &[u8], offset: usize, length: usize) -> Result<&[u8]> {
    bytes
        .get(
            offset
                ..offset
                    .checked_add(length)
                    .ok_or_else(|| bad("Section overflow"))?,
        )
        .ok_or_else(|| bad("Section extends beyond file"))
}
pub fn read(name: &str, bytes: &[u8]) -> Result<Document> {
    if bytes.len() < 128 || !bytes.starts_with(b".scf") {
        return Err(bad("Invalid header"));
    }
    let samples = integer(bytes, 4)?;
    let bases = integer(bytes, 12)?;
    let version = bytes[36];
    let sample_size = if version < b'2' {
        1
    } else {
        integer(bytes, 40)?
    };
    if samples == 0
        || samples > 2_000_000
        || bases == 0
        || bases > popgene_core::MAX_BASES
        || !matches!(sample_size, 1 | 2)
        || !matches!(version, b'1' | b'2' | b'3')
    {
        return Err(bad("Unsupported version, sample precision, or count"));
    }
    let signal = section(bytes, integer(bytes, 8)?, samples * 4 * sample_size)?;
    let calls = section(bytes, integer(bytes, 24)?, bases * 12)?;
    let mut channels = vec![vec![0u16; samples]; 4];
    for (channel, values) in channels.iter_mut().enumerate() {
        for (i, value) in values.iter_mut().enumerate() {
            let index = if version == b'3' {
                channel * samples + i
            } else {
                i * 4 + channel
            };
            *value = if sample_size == 1 {
                signal[index] as u16
            } else {
                u16::from_be_bytes(signal[index * 2..index * 2 + 2].try_into().unwrap())
            };
        }
        if version == b'3' {
            for _ in 0..2 {
                let mut previous = 0u16;
                for value in values.iter_mut() {
                    *value = value.wrapping_add(previous);
                    if sample_size == 1 {
                        *value &= 255;
                    }
                    previous = *value;
                }
            }
        }
    }
    let mut sequence = String::with_capacity(bases);
    let mut peaks = Vec::with_capacity(bases);
    let mut qualities = Vec::with_capacity(bases);
    for i in 0..bases {
        let peak = if version == b'3' {
            integer(calls, i * 4)?
        } else {
            integer(calls, i * 12)?
        };
        let base = if version == b'3' {
            calls[8 * bases + i]
        } else {
            calls[i * 12 + 8]
        }
        .to_ascii_uppercase();
        let channel = b"ACGT".iter().position(|&b| b == base);
        // Ambiguous SCF calls have no uniquely designated confidence channel.
        // The largest stored estimate is displayed, and all four original values
        // remain in the source. This display convention does not recalibrate them.
        let probabilities = std::array::from_fn::<_, 4, _>(|c| {
            if version == b'3' {
                calls[(4 + c) * bases + i]
            } else {
                calls[i * 12 + 4 + c]
            }
        });
        qualities.push(channel.map_or_else(
            || *probabilities.iter().max().unwrap(),
            |c| probabilities[c],
        ));
        sequence.push(if base == b'-' { 'N' } else { base as char });
        peaks.push(peak);
    }
    let mut document = Document::new(
        std::path::Path::new(name)
            .file_stem()
            .unwrap_or_default()
            .to_string_lossy(),
        &sequence,
    )?;
    let comments = section(bytes, integer(bytes, 32)?, integer(bytes, 28)?)?;
    section(bytes, integer(bytes, 52)?, integer(bytes, 48)?)?;
    document.comments = String::from_utf8_lossy(comments)
        .trim_end_matches('\0')
        .lines()
        .map(str::to_string)
        .collect();
    document.trace = Some(Trace {
        calls_version: version as u32 - b'0' as u32,
        qualities: qualities.clone(),
        qualities_known: true,
        peaks: peaks.clone(),
        channel_order: "ACGT".into(),
        channels,
        original_calls: vec![TraceCallSet {
            label: "Original SCF calls".into(),
            version: version as u32 - b'0' as u32,
            sequence,
            qualities,
            qualities_known: true,
            peaks,
        }],
        edited: false,
        sample_range: None,
        manual_calls: vec![],
    });
    document.source = Some(Source {
        file_name: name.into(),
        format: "SCF".into(),
        bytes: bytes.into(),
        retained_blocks: vec![],
    });
    document.description = "Imported SCF chromatogram".into();
    document.warnings.push("SCF base confidence values are displayed as stored; older instruments may use a different confidence scale. Ambiguous calls show the largest of the four stored estimates. Original four-channel confidence arrays, comments and private records remain in the retained source.".into());
    document.validate()?;
    Ok(document)
}
pub fn write(document: &Document) -> Result<Vec<u8>> {
    document.validate()?;
    let trace = document
        .trace
        .as_ref()
        .ok_or_else(|| bad("No chromatogram"))?;
    let samples = trace.channels[0].len();
    if !trace.qualities_known {
        return Err(bad(
            "No recorded confidence values; use ZTR to preserve an unknown-quality chromatogram",
        ));
    }
    let bases = document.sequence.len();
    let mut output = vec![0; 128];
    output[..4].copy_from_slice(b".scf");
    output[36..40].copy_from_slice(b"3.00");
    let set = |output: &mut Vec<u8>, offset: usize, value: usize| {
        output[offset..offset + 4].copy_from_slice(&(value as u32).to_be_bytes())
    };
    set(&mut output, 4, samples);
    set(&mut output, 8, 128);
    set(&mut output, 12, bases);
    set(&mut output, 40, 2);
    set(&mut output, 44, 2);
    for base in b"ACGT" {
        let channel = trace
            .channel_order
            .as_bytes()
            .iter()
            .position(|b| b == base)
            .unwrap();
        let mut delta = trace.channels[channel].clone();
        for _ in 0..2 {
            let mut previous = 0;
            for value in &mut delta {
                let sample = *value;
                *value = value.wrapping_sub(previous);
                previous = sample;
            }
        }
        for value in delta {
            output.extend(value.to_be_bytes());
        }
    }
    let base_offset = output.len();
    set(&mut output, 24, base_offset);
    for &peak in &trace.peaks {
        output.extend((peak as u32).to_be_bytes());
    }
    // Preserve other per-base confidence channels when source calls have an exact
    // peak/base/quality match. Manually changed calls get only their current score.
    let original = document
        .source
        .as_ref()
        .filter(|s| s.format == "SCF")
        .and_then(|source| {
            read(&source.file_name, &source.bytes)
                .ok()
                .map(|d| (source, d))
        });
    let original_indices = original.as_ref().map(|(_, old)| {
        let old_trace = old.trace.as_ref().unwrap();
        let mut candidates =
            std::collections::HashMap::<(usize, u8, u8), std::collections::VecDeque<usize>>::new();
        for (index, &peak) in old_trace.peaks.iter().enumerate() {
            candidates
                .entry((
                    peak,
                    old.sequence.as_bytes()[index],
                    old_trace.qualities[index],
                ))
                .or_default()
                .push_back(index);
        }
        trace
            .peaks
            .iter()
            .enumerate()
            .map(|(index, &peak)| {
                candidates
                    .get_mut(&(
                        peak,
                        document.sequence.as_bytes()[index],
                        trace.qualities[index],
                    ))
                    .and_then(|indices| indices.pop_front())
            })
            .collect::<Vec<_>>()
    });
    for base in b"ACGT" {
        let channel = b"ACGT".iter().position(|b| b == base).unwrap();
        for (i, &called) in document.sequence.as_bytes().iter().enumerate() {
            let mut score = if called == *base || !b"ACGT".contains(&called) {
                trace.qualities[i]
            } else {
                0
            };
            if let Some((source, old)) = &original
                && let Some(index) = original_indices.as_ref().and_then(|indices| indices[i])
            {
                let offset = integer(&source.bytes, 24)?;
                score = if source.bytes[36] == b'3' {
                    source.bytes[offset + (4 + channel) * old.sequence.len() + index]
                } else {
                    source.bytes[offset + index * 12 + 4 + channel]
                };
            }
            output.push(score);
        }
    }
    output.extend(document.sequence.as_bytes());
    output.resize(output.len() + 3 * bases, 0);
    let comments = format!("{}\0", document.comments.join("\n"));
    let comments_offset = output.len();
    set(&mut output, 28, comments.len());
    set(&mut output, 32, comments_offset);
    output.extend(comments.as_bytes());
    if let Some(source) = document.source.as_ref().filter(|s| s.format == "SCF") {
        let private = section(
            &source.bytes,
            integer(&source.bytes, 52)?,
            integer(&source.bytes, 48)?,
        )?;
        let offset = output.len();
        set(&mut output, 48, private.len());
        set(&mut output, 52, offset);
        output.extend(private);
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
    fn version3_double_delta_roundtrip_and_truncated_sections() {
        let doc = crate::import(
            "trace.ab1",
            include_bytes!("../../../validation/fixtures/synthetic-complete.ab1"),
        )
        .unwrap();
        let bytes = write(&doc).unwrap();
        let reopened = read("trace.scf", &bytes).unwrap();
        let a = doc.trace.unwrap();
        let b = reopened.trace.unwrap();
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
        for end in [0, 127, bytes.len() - 1] {
            assert!(read("bad.scf", &bytes[..end]).is_err());
        }
    }
    #[test]
    fn independent_v2_interleaved_fixture() {
        let mut bytes = vec![0; 128];
        bytes[..4].copy_from_slice(b".scf");
        bytes[36..40].copy_from_slice(b"2.00");
        for (offset, value) in [(4, 3u32), (8, 128), (12, 1), (24, 140), (40, 1)] {
            bytes[offset..offset + 4].copy_from_slice(&value.to_be_bytes());
        }
        bytes.extend([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]);
        bytes.extend([0, 0, 0, 1, 11, 22, 33, 44, b'G', 0, 0, 0]);
        let d = read("v2.scf", &bytes).unwrap();
        assert_eq!(d.sequence.as_ref(), "G");
        let t = d.trace.unwrap();
        assert_eq!(t.peaks, vec![1]);
        assert_eq!(t.qualities, vec![33]);
        assert_eq!(t.channels[0], vec![1, 5, 9]);
        assert_eq!(t.channels[3], vec![4, 8, 12]);
    }
}
