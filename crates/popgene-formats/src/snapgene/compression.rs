//! Bounded native compression. Limits apply before allocating and after decoding.
use super::{Result, bad};
use std::io::{Read, Write};
pub(super) fn qdecode(bytes: &[u8], limit: usize) -> Result<Vec<u8>> {
    let header = bytes
        .get(..4)
        .ok_or_else(|| bad("Truncated compressed size"))?;
    let size = u32::from_be_bytes(header.try_into().unwrap()) as usize;
    if size > limit {
        return Err(bad("Compressed record expands beyond the limit"));
    }
    let mut decoded = vec![
        0;
        size.checked_add(1)
            .ok_or_else(|| bad("Compressed size overflow"))?
    ];
    let mut decoder = flate2::Decompress::new(true);
    let status = decoder
        .decompress(&bytes[4..], &mut decoded, flate2::FlushDecompress::Finish)
        .map_err(|e| bad(format!("Invalid compressed record: {e}")))?;
    if status != flate2::Status::StreamEnd
        || decoder.total_out() != size as u64
        || decoder.total_in() != (bytes.len() - 4) as u64
    {
        return Err(bad(
            "Compressed record size, checksum or stream boundary is invalid",
        ));
    }
    decoded.truncate(size);
    Ok(decoded)
}
pub(super) fn xdecode(bytes: &[u8], limit: usize) -> Result<Vec<u8>> {
    let mut reader =
        lzma_rust2::XzReader::new_mem_limit(std::io::Cursor::new(bytes), false, 64 * 1024);
    let mut decoded = vec![];
    (&mut reader)
        .take(limit as u64 + 1)
        .read_to_end(&mut decoded)?;
    if decoded.len() > limit || reader.into_inner().position() != bytes.len() as u64 {
        return Err(bad("Native XZ record exceeds limits or has trailing data"));
    }
    Ok(decoded)
}
pub(super) fn xencode(bytes: &[u8]) -> Result<Vec<u8>> {
    let mut writer = lzma_rust2::XzWriter::new(Vec::new(), lzma_rust2::XzOptions::with_preset(6))?;
    writer.write_all(bytes)?;
    Ok(writer.finish()?)
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn compressed_records_verify_lengths_checksums_and_limits() {
        let raw = b"Synthetic embedded data";
        let mut z = flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        z.write_all(raw).unwrap();
        let mut q = (raw.len() as u32).to_be_bytes().to_vec();
        q.extend(z.finish().unwrap());
        assert_eq!(qdecode(&q, raw.len()).unwrap(), raw);
        assert!(qdecode(&q, raw.len() - 1).is_err());
        assert!(qdecode(&q[..q.len() - 1], raw.len()).is_err());
        q.push(0);
        assert!(qdecode(&q, raw.len()).is_err());
        let x = xencode(raw).unwrap();
        assert_eq!(xdecode(&x, raw.len()).unwrap(), raw);
        assert!(xdecode(&x, raw.len() - 1).is_err());
        assert!(xdecode(&x[..x.len() - 1], raw.len()).is_err());
    }
}
