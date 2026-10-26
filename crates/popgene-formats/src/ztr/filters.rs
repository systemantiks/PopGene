//! Bounded standard ZTR filtering. Header expansion lengths for RLE/zlib are
//! little endian in the Staden reference codec (not chunk integer byte order).
use super::{Result, bad};
use std::io::Read;
#[path = "huffman.rs"]
mod huffman;
pub(super) use huffman::HuffmanTables;
#[path = "extended_filters.rs"]
mod extended_filters;
#[path = "floating_filter.rs"]
mod floating_filter;
#[derive(Clone, Copy, Default)]
pub(super) struct DecodeWarnings {
    pub logarithmic: bool,
    pub legacy_float: bool,
}
const LIMIT: usize = 64 * 1024 * 1024;
fn byte(bytes: &[u8], offset: usize) -> Result<u8> {
    bytes
        .get(offset)
        .copied()
        .ok_or_else(|| bad("Truncated compression layer"))
}
fn expand_length(bytes: &[u8]) -> Result<usize> {
    let n = u32::from_le_bytes(
        bytes
            .get(1..5)
            .ok_or_else(|| bad("Truncated compression header"))?
            .try_into()
            .unwrap(),
    ) as usize;
    if n > LIMIT {
        return Err(bad("Expanded chunk exceeds 64 MiB"));
    }
    Ok(n)
}
#[cfg(test)]
pub(super) fn decode(input: &[u8]) -> Result<Vec<u8>> {
    decode_with_warnings(input, None).map(|result| result.0)
}
/// The flag records a lossy LOG2 layer, including beneath other compression.
#[cfg(test)]
pub(super) fn decode_detailed(input: &[u8], bases: Option<&[u8]>) -> Result<(Vec<u8>, bool)> {
    decode_with_warnings(input, bases).map(|(data, warnings)| (data, warnings.logarithmic))
}
pub(super) fn decode_with_warnings(
    input: &[u8],
    bases: Option<&[u8]>,
) -> Result<(Vec<u8>, DecodeWarnings)> {
    decode_with_tables(input, bases, &HuffmanTables::default())
}
pub(super) fn decode_with_tables(
    input: &[u8],
    bases: Option<&[u8]>,
    tables: &HuffmanTables,
) -> Result<(Vec<u8>, DecodeWarnings)> {
    if input.len() > LIMIT {
        return Err(bad("Compressed chunk exceeds 64 MiB"));
    }
    let mut data = input.to_vec();
    let mut work = 0usize;
    let mut warnings = DecodeWarnings::default();
    for _ in 0..32 {
        work = work.saturating_add(data.len());
        if data.len() > LIMIT || work > 256 * 1024 * 1024 {
            return Err(bad("Compression expansion budget exceeded"));
        }
        let format = byte(&data, 0)?;
        data = match format {
            0 => return Ok((data, warnings)),
            3 => extended_filters::xrle(&data)?,
            4 => extended_filters::xrle2(&data)?,
            1 => {
                let expected = expand_length(&data)?;
                let guard = byte(&data, 5)?;
                let mut result = Vec::with_capacity(expected);
                let mut i = 6;
                while i < data.len() {
                    let value = byte(&data, i)?;
                    i += 1;
                    if value == guard {
                        let count = byte(&data, i)? as usize;
                        i += 1;
                        if count == 0 {
                            result.push(guard);
                        } else {
                            let value = byte(&data, i)?;
                            i += 1;
                            if result.len() + count > expected {
                                return Err(bad("RLE exceeds declared length"));
                            }
                            result.extend(std::iter::repeat_n(value, count));
                        }
                    } else {
                        result.push(value);
                    }
                    if result.len() > expected {
                        return Err(bad("RLE exceeds declared length"));
                    }
                }
                if result.len() != expected {
                    return Err(bad("RLE length mismatch"));
                }
                result
            }
            2 => {
                let expected = expand_length(&data)?;
                let mut result = Vec::with_capacity(expected);
                let mut decoder = flate2::read::ZlibDecoder::new(&data[5..]);
                decoder
                    .by_ref()
                    .take(expected as u64 + 1)
                    .read_to_end(&mut result)
                    .map_err(|_| bad("Invalid zlib stream"))?;
                if result.len() != expected || decoder.total_in() as usize != data.len() - 5 {
                    return Err(bad("Zlib length or trailing data mismatch"));
                }
                result
            }
            64..=66 => {
                let size = 1usize << (format - 64);
                let offset = if size == 4 { 4 } else { 2 };
                let levels = byte(&data, 1)?;
                if !(1..=3).contains(&levels)
                    || data.len() < offset
                    || !(data.len() - offset).is_multiple_of(size)
                {
                    return Err(bad("Invalid delta layer"));
                }
                let mut values = data[offset..]
                    .chunks_exact(size)
                    .map(|chunk| chunk.iter().fold(0u32, |v, b| (v << 8) | *b as u32))
                    .collect::<Vec<_>>();
                let mask = if size == 4 {
                    u32::MAX
                } else {
                    (1u32 << (size * 8)) - 1
                };
                for _ in 0..levels {
                    let mut previous = 0u32;
                    for value in &mut values {
                        *value = value.wrapping_add(previous) & mask;
                        previous = *value;
                    }
                }
                values
                    .into_iter()
                    .flat_map(|v| v.to_be_bytes()[4 - size..].to_vec())
                    .collect()
            }
            70 | 71 => {
                let size = if format == 70 { 2 } else { 4 };
                let mut result = vec![];
                let mut i = 1;
                while i < data.len() {
                    let value = data[i] as i8;
                    i += 1;
                    if value == i8::MIN {
                        let bytes = data
                            .get(i..i + size)
                            .ok_or_else(|| bad("Truncated integer packing escape"))?;
                        result.extend(bytes);
                        i += size;
                    } else {
                        result.extend(&(value as i32).to_be_bytes()[4 - size..]);
                    }
                    if result.len() > LIMIT {
                        return Err(bad("Packed integer expansion limit"));
                    }
                }
                result
            }
            72 => {
                if data.len() < 258 {
                    return Err(bad("Truncated follow predictor"));
                }
                let mut result = Vec::with_capacity(data.len() - 257);
                result.push(data[257]);
                for &residual in &data[258..] {
                    let previous = *result.last().unwrap() as usize;
                    result.push(data[1 + previous].wrapping_sub(residual));
                }
                result
            }
            73 => {
                warnings.legacy_float = true;
                floating_filter::decode(&data)?
            }
            74 => integer_chebyshev(&data)?,
            77 => huffman::decode(&data, tables)?,
            75 => {
                warnings.logarithmic = true;
                extended_filters::log2(&data)?
            }
            79 => extended_filters::qshift(&data)?,
            80 => extended_filters::tshift(&data, bases)?,
            value => {
                return Err(bad(&format!(
                    "Unsupported ZTR compression method {value}; source bytes were not modified"
                )));
            }
        };
    }
    Err(bad("More than 32 compression layers"))
}
// Adapted from Staden io_lib ichebuncomp (James Bonfield), BSD-3-Clause;
// copyright notices and terms: vendor-notices/staden-io-lib/LICENSE.txt.
fn integer_chebyshev(data: &[u8]) -> Result<Vec<u8>> {
    if data.len() < 2 || !data.len().is_multiple_of(2) {
        return Err(bad("Invalid integer Chebyshev layer"));
    }
    let residuals = data[2..]
        .as_chunks::<2>()
        .0
        .iter()
        .map(|p| u16::from_be_bytes(*p))
        .collect::<Vec<_>>();
    let mut decoded = Vec::<u16>::with_capacity(residuals.len());
    for (index, &residual) in residuals.iter().enumerate() {
        let predicted = if index < 4 {
            decoded.last().copied().unwrap_or(0) as i64
        } else {
            let d = &decoded[index - 4..];
            let d = [d[0] as i64, d[1] as i64, d[2] as i64, d[3] as i64];
            let f = [
                d[2] * 11 + d[3] * 139,
                d[2] * 93 + d[3] * 57,
                d[1] * 75 + d[2] * 75,
                d[0] * 57 + d[1] * 93,
                d[0] * 139 + d[1] * 11,
            ];
            let weights = [
                [42, 42, 42, 42, 42],
                [39, 24, 0, -24, -39],
                [33, -12, -42, -12, 33],
                [24, -39, 0, 39, -24],
            ];
            let mut coefficients =
                weights.map(|weights| f.iter().zip(weights).map(|(v, w)| v * w).sum::<i64>());
            let maximum = coefficients.iter().map(|v| v.abs()).max().unwrap();
            let factor = if maximum > 1 << 26 {
                maximum / (1 << 26) + 1
            } else {
                1
            };
            for value in &mut coefficients {
                *value /= factor;
            }
            let [a, b, c, d] = coefficients;
            let dd = (d / 3) * 10 + c;
            (((((dd / 3) * 10 - d + b) / 3) * 5 - dd + a / 2) / (150 * 105) * factor).max(0)
        };
        decoded.push((predicted + residual as i16 as i64) as u16);
    }
    Ok(decoded.into_iter().flat_map(u16::to_be_bytes).collect())
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn independent_delta_and_rle_vectors() {
        assert_eq!(
            decode(&[64, 2, 0, 10, 0, 236, 200, 56, 81]).unwrap(),
            [0, 10, 20, 10, 200, 190, 5]
        );
        assert_eq!(
            decode(&[1, 5, 0, 0, 0, 8, 0, 8, 3, 7, 8, 0]).unwrap(),
            [0, 7, 7, 7, 8]
        );
        assert!(decode(&[1, 1, 0, 0, 0, 8, 8, 255, 7]).is_err());
        assert!(decode(&[65, 1, 0]).is_err());
    }
    #[test]
    fn packed_signed_values_and_predictor() {
        assert_eq!(
            decode(&[70, 0, 0, 5, -5i8 as u8, 128, 0, 200]).unwrap(),
            [0, 0, 0, 0, 0, 5, 255, 251, 0, 200]
        );
        let mut data = vec![72];
        data.extend(0u8..=255);
        data.extend([0, 255, 255, 254]);
        assert_eq!(decode(&data).unwrap(), [0, 1, 2, 4]);
        assert_eq!(
            decode(&[74, 0, 0, 0, 0, 10, 255, 251, 0, 20]).unwrap(),
            [0, 0, 0, 10, 0, 5, 0, 25]
        );
    }
}
