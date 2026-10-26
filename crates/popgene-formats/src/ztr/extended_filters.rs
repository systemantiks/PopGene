//! Bounded XRLE, XRLE2, LOG2, QSHIFT and TSHIFT ZTR layers.
//! Adapted from James Bonfield's Staden io_lib compression.c (BSD-3-Clause).
//! Notices: vendor-notices/staden-io-lib/LICENSE.txt.
//! https://github.com/jkbonfield/io_lib/blob/master/io_lib/compression.c
use super::{LIMIT, Result, bad};

fn append(output: &mut Vec<u8>, word: &[u8], count: usize) -> Result<()> {
    let size = word
        .len()
        .checked_mul(count)
        .and_then(|n| output.len().checked_add(n))
        .filter(|n| *n <= LIMIT)
        .ok_or_else(|| bad("Extended RLE expansion exceeds 64 MiB"))?;
    if size > output.capacity() {
        let capacity = size.max(output.capacity().saturating_mul(2)).min(LIMIT);
        output.reserve_exact(capacity - output.len());
    }
    for _ in 0..count {
        output.extend_from_slice(word);
    }
    Ok(())
}
pub(super) fn xrle(data: &[u8]) -> Result<Vec<u8>> {
    if data.len() < 3 || data[1] == 0 {
        return Err(bad("Invalid XRLE block size/header"));
    }
    let size = data[1] as usize;
    let guard = data[2];
    let mut output = vec![];
    let mut i = 3;
    while i < data.len() {
        let value = data[i];
        i += 1;
        if value != guard {
            append(&mut output, &[value], 1)?;
            continue;
        }
        let count = *data.get(i).ok_or_else(|| bad("Truncated XRLE guard"))? as usize;
        i += 1;
        if count == 0 {
            append(&mut output, &[guard], 1)?;
            continue;
        }
        let word = data
            .get(i..i + size)
            .ok_or_else(|| bad("Truncated XRLE repeated word"))?;
        append(&mut output, word, count)?;
        i += size;
    }
    Ok(output)
}
pub(super) fn xrle2(data: &[u8]) -> Result<Vec<u8>> {
    let size = *data.get(1).ok_or_else(|| bad("Truncated XRLE2 header"))? as usize;
    let header = size.max(2);
    if size == 0 || data.len() < header || !(data.len() - header).is_multiple_of(size) {
        return Err(bad("Invalid XRLE2 word alignment"));
    }
    let mut output = vec![];
    let mut i = header;
    let mut previous: Option<&[u8]> = None;
    let mut run = 0;
    while i < data.len() {
        let word = &data[i..i + size];
        i += size;
        append(&mut output, word, 1)?;
        run = if previous == Some(word) { run + 1 } else { 1 };
        if run == 2 {
            let marker = data
                .get(i..i + size)
                .ok_or_else(|| bad("Truncated XRLE2 repeat count"))?;
            append(&mut output, word, marker[0] as usize)?;
            i += size;
            run = 0;
        }
        previous = Some(word);
    }
    Ok(output)
}
pub(super) fn log2(data: &[u8]) -> Result<Vec<u8>> {
    if data.len() < 2 || !data.len().is_multiple_of(2) {
        return Err(bad("Invalid LOG2 sample alignment"));
    }
    let mut output = Vec::with_capacity(data.len() - 2);
    for word in data[2..].as_chunks::<2>().0 {
        let value = u16::from_be_bytes(*word);
        // The forward transform is floor(10*log2(sample+1)) for u16 samples.
        // Larger codes cannot come from a valid sample and would overflow C's
        // float-to-integer conversion; reject rather than wrap or saturate.
        if value > 160 {
            return Err(bad("LOG2 code exceeds the 16-bit signal range"));
        }
        let sample = (2_f64.powf(f64::from(value) / 10.) as u32 - 1) as u16;
        output.extend(sample.to_be_bytes());
    }
    Ok(output)
}
pub(super) fn qshift(data: &[u8]) -> Result<Vec<u8>> {
    if data.len() < 4 || !data.len().is_multiple_of(4) || data[3] != 0 {
        return Err(bad("Invalid QSHIFT header or alignment"));
    }
    let n = data.len() / 4 - 1;
    let mut output = vec![0; 1 + 4 * n];
    for (i, values) in data[4..].as_chunks::<4>().0.iter().enumerate() {
        output[1 + i] = values[0];
        output[1 + n + 3 * i..1 + n + 3 * i + 3].copy_from_slice(&values[1..]);
    }
    Ok(output)
}
pub(super) fn tshift(data: &[u8], bases: Option<&[u8]>) -> Result<Vec<u8>> {
    let bases = bases.ok_or_else(|| bad("TSHIFT requires decoded BASE calls"))?;
    if data.len() < 8 || !(data.len() - 8).is_multiple_of(8) || (data.len() - 8) / 8 != bases.len()
    {
        return Err(bad("TSHIFT signal length differs from BASE calls"));
    }
    let mut output = vec![0; 2 + bases.len() * 8];
    for (i, (&base, values)) in bases.iter().zip(data[8..].as_chunks::<8>().0).enumerate() {
        let order = match base {
            b'C' => [1, 0, 2, 3],
            b'G' => [2, 0, 1, 3],
            b'T' => [3, 0, 1, 2],
            _ => [0, 1, 2, 3],
        };
        for (column, channel) in order.into_iter().enumerate() {
            let offset = 2 + 2 * (channel * bases.len() + i);
            output[offset..offset + 2].copy_from_slice(&values[column * 2..column * 2 + 2]);
        }
    }
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn independent_staden_extended_vectors_and_complete_files() {
        fn hex(text: &str) -> Vec<u8> {
            text.as_bytes()
                .as_chunks::<2>()
                .0
                .iter()
                .map(|v| u8::from_str_radix(std::str::from_utf8(v).unwrap(), 16).unwrap())
                .collect()
        }
        let facts: serde_json::Value = serde_json::from_str(include_str!(
            "../../../../validation/fixtures/trace/extended/vectors.json"
        ))
        .unwrap();
        for vector in facts["vectors"].as_array().unwrap() {
            let data = hex(vector["input"].as_str().unwrap());
            let expected = hex(vector["expected"].as_str().unwrap());
            let bases = vector["bases"].as_str().map(hex);
            let result = super::super::decode_detailed(&data, bases.as_deref()).unwrap();
            assert_eq!(result.0, expected, "{}", vector["name"]);
            assert_eq!(result.1, vector["lossy"].as_bool().unwrap());
        }
        for fixture in facts["fixtures"].as_array().unwrap() {
            let name = fixture["name"].as_str().unwrap();
            let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("../../validation/fixtures/trace/extended")
                .join(name);
            let data = std::fs::read(path).unwrap();
            let document = crate::import(name, &data).unwrap();
            assert_eq!(
                document.sequence.as_ref(),
                fixture["sequence"].as_str().unwrap()
            );
            let trace = document.trace.as_ref().unwrap();
            let channels: Vec<Vec<u16>> =
                serde_json::from_value(fixture["channels"].clone()).unwrap();
            let qualities: Vec<u8> = serde_json::from_value(fixture["qualities"].clone()).unwrap();
            assert_eq!(trace.channels, channels);
            assert_eq!(trace.qualities, qualities);
            assert_eq!(
                document.warnings.iter().any(|w| w.contains("LOG2")),
                fixture["lossy"].as_bool().unwrap()
            );
        }
    }
    #[test]
    fn rejects_truncation_alignment_and_expansion_overflow() {
        for data in [&[3, 0, 42][..], &[3, 4, 42, 42], &[3, 4, 42, 42, 2, 1]] {
            assert!(xrle(data).is_err());
        }
        for data in [&[4, 0][..], &[4, 2, 0], &[4, 2, 7, 8, 7, 8]] {
            assert!(xrle2(data).is_err());
        }
        assert!(log2(&[75, 0, 0]).is_err());
        assert!(log2(&[75, 0, 0, 161]).is_err());
        assert!(qshift(&[79, 0, 0]).is_err());
        assert!(qshift(&[79, 0, 0, 1]).is_err());
        assert!(tshift(&[80; 8], None).is_err());
        assert!(tshift(&[80; 8], Some(b"A")).is_err());
        let mut output = vec![0; LIMIT - 1];
        assert!(append(&mut output, &[1, 2], 1).is_err());
    }
}
