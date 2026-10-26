//! Bounded section-based DNA decoding. Layout/alphabet reference: sgffp 0.22.1;
//! attribution and MIT license: vendor-notices/sgffp-reference/.
use super::{Result, bad};
struct Input<'a> {
    bytes: &'a [u8],
    offset: usize,
}
impl<'a> Input<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        let end = self
            .offset
            .checked_add(n)
            .filter(|end| *end <= self.bytes.len())
            .ok_or_else(|| bad("Truncated packed DNA"))?;
        let out = &self.bytes[self.offset..end];
        self.offset = end;
        Ok(out)
    }
    fn number(&mut self) -> Result<usize> {
        Ok(u32::from_be_bytes(self.take(4)?.try_into().unwrap()) as usize)
    }
    fn byte(&mut self) -> Result<u8> {
        Ok(self.take(1)?[0])
    }
}
pub(super) fn decode(bytes: &[u8]) -> Result<String> {
    let mut input = Input { bytes, offset: 0 };
    let length = input.number()?;
    if length == 0 || length > popgene_core::MAX_BASES {
        return Err(bad("Packed DNA length exceeds limits"));
    }
    let _writer_stamp = input.byte()?;
    let sections = input.number()?;
    let lowercase = input.number()?;
    if sections == 0 || sections > length || lowercase > length {
        return Err(bad("Invalid packed DNA section count"));
    }
    let mut out = Vec::with_capacity(length);
    for _ in 0..sections {
        let kind = input.byte()?;
        let count = input.number()?;
        if count == 0 || out.len().checked_add(count).is_none_or(|n| n > length) {
            return Err(bad("Packed DNA section exceeds declared sequence"));
        }
        match kind {
            1 => {
                let data = input.take(count.div_ceil(4))?;
                for (i, byte) in data.iter().enumerate() {
                    let symbols = (count - i * 4).min(4);
                    for j in 0..symbols {
                        out.push(b"GATC"[((byte >> ((symbols - 1 - j) * 2)) & 3) as usize]);
                    }
                }
            }
            2 => {
                let data = input.take(count.div_ceil(2))?;
                for (i, byte) in data.iter().enumerate() {
                    let symbols = (count - i * 2).min(2);
                    for j in 0..symbols {
                        let symbol = (byte >> ((symbols - 1 - j) * 4)) & 15;
                        let base = match symbol {
                            4 => b'N',
                            5 => b'B',
                            6 => b'D',
                            7 => b'H',
                            8 => b'K',
                            9 => b'M',
                            10 => b'R',
                            11 => b'S',
                            12 => b'V',
                            13 => b'W',
                            14 => b'Y',
                            _ => return Err(bad("Unknown packed IUPAC symbol")),
                        };
                        out.push(base);
                    }
                }
            }
            3 => out.resize(out.len() + count, b'N'),
            _ => return Err(bad("Unknown packed DNA section kind")),
        }
    }
    if out.len() != length {
        return Err(bad("Packed DNA length mismatch"));
    }
    for _ in 0..lowercase {
        let start = input.number()?;
        let end = input.number()?;
        if start > end || end >= length {
            return Err(bad("Invalid packed lowercase range"));
        }
        out[start..=end].make_ascii_lowercase();
    }
    if input.offset != bytes.len() {
        return Err(bad("Trailing packed DNA data"));
    }
    String::from_utf8(out).map_err(|_| bad("Invalid packed sequence"))
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn native_observed_snapshot_and_independent_partial_sections() {
        // Observed native 12-base input (case retained until shared normalization).
        let b = [
            0, 0, 0, 12, 2, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 12, 0x97, 0xf0, 0x2a, 0, 0, 0, 1,
            0, 0, 0, 11,
        ];
        assert_eq!(decode(&b).unwrap(), "Taacccgggttt");
        let mut b = vec![];
        b.extend(8_u32.to_be_bytes());
        b.push(2);
        b.extend(3_u32.to_be_bytes());
        b.extend(0_u32.to_be_bytes());
        b.push(1);
        b.extend(3_u32.to_be_bytes());
        b.push(0x1b);
        b.push(2);
        b.extend(3_u32.to_be_bytes());
        b.extend([0x4e, 0x0a]);
        b.push(3);
        b.extend(2_u32.to_be_bytes());
        assert_eq!(decode(&b).unwrap(), "ATCNYRNN");
        assert!(decode(&b[..b.len() - 1]).is_err());
        b.push(0);
        assert!(decode(&b).is_err());
    }
}
