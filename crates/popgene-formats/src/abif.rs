//! Bounded, read-only Sanger ABIF import (edited calls preferred; analysed DATA9..12).
use super::{Error, Result};
use popgene_core::{Document, Source, Trace};
use std::{collections::BTreeMap, path::Path};
fn bad(s: &str) -> Error {
    Error::Message(format!("ABIF import: {s}"))
}
fn u32_at(b: &[u8], p: usize) -> Result<usize> {
    Ok(u32::from_be_bytes(
        b.get(p..p + 4)
            .ok_or_else(|| bad("Truncated integer"))?
            .try_into()
            .unwrap(),
    ) as usize)
}
struct Entry<'a> {
    data_type: u16,
    size: usize,
    count: usize,
    bytes: &'a [u8],
}
fn directory(b: &[u8]) -> Result<BTreeMap<([u8; 4], usize), Entry<'_>>> {
    if b.len() < 34 || &b[..4] != b"ABIF" {
        return Err(bad("Invalid header"));
    }
    let count = u32_at(b, 18)?;
    let offset = u32_at(b, 26)?;
    if count > 10000 || offset.checked_add(count * 28).is_none_or(|n| n > b.len()) {
        return Err(bad("Invalid directory bounds"));
    }
    let mut entries = BTreeMap::new();
    for i in 0..count {
        let p = offset + i * 28;
        let tag: [u8; 4] = b[p..p + 4].try_into().unwrap();
        let number = u32_at(b, p + 4)?;
        let size = u16::from_be_bytes(b[p + 10..p + 12].try_into().unwrap()) as usize;
        let elements = u32_at(b, p + 12)?;
        let length = u32_at(b, p + 16)?;
        if size.checked_mul(elements) != Some(length) {
            return Err(bad("Inconsistent element sizes"));
        }
        let start = if length <= 4 {
            p + 20
        } else {
            u32_at(b, p + 20)?
        };
        let end = start
            .checked_add(length)
            .filter(|&n| n <= b.len())
            .ok_or_else(|| bad("Truncated tag"))?;
        if entries
            .insert(
                (tag, number),
                Entry {
                    data_type: u16::from_be_bytes(b[p + 8..p + 10].try_into().unwrap()),
                    size,
                    count: elements,
                    bytes: &b[start..end],
                },
            )
            .is_some()
        {
            return Err(bad("Duplicate directory tag"));
        }
    }
    Ok(entries)
}
pub fn read(name: &str, b: &[u8]) -> Result<Document> {
    let entries = directory(b)?;
    let get = |tag: &[u8; 4], v| {
        entries
            .get(&(*tag, v))
            .ok_or_else(|| bad("Required analysed trace tag is absent"))
    };
    let version = if entries.contains_key(&(*b"PBAS", 1)) {
        1
    } else {
        2
    };
    let calls = get(b"PBAS", version)?;
    if calls.size != 1 {
        return Err(bad("Unsupported call element size"));
    }
    let sequence = std::str::from_utf8(calls.bytes)
        .map_err(|_| bad("Invalid base-call encoding"))?
        .trim_end_matches('\0');
    let mut d = Document::new(
        Path::new(name)
            .file_stem()
            .unwrap_or_default()
            .to_string_lossy(),
        sequence,
    )?;
    let quality = get(b"PCON", version)?;
    let peaks = get(b"PLOC", version)?;
    if quality.size != 1
        || quality.count < sequence.len()
        || peaks.count < sequence.len()
        || !matches!(peaks.size, 2 | 4)
    {
        return Err(bad("Inconsistent call/quality/peak arrays"));
    }
    let peaks: Vec<usize> = peaks
        .bytes
        .chunks_exact(peaks.size)
        .take(sequence.len())
        .map(|p| {
            if p.len() == 2 {
                u16::from_be_bytes(p.try_into().unwrap()) as usize
            } else {
                u32::from_be_bytes(p.try_into().unwrap()) as usize
            }
        })
        .collect();
    let order = std::str::from_utf8(get(b"FWO_", 1)?.bytes)
        .map_err(|_| bad("Invalid channel order"))?
        .trim_end_matches('\0')
        .to_string();
    let mut channels = vec![];
    for v in 9..=12 {
        let e = get(b"DATA", v)?;
        if e.size != 2 || e.count > 2_000_000 {
            return Err(bad(
                "Unsupported signal element size or more than 2 million samples",
            ));
        }
        channels.push(
            e.bytes
                .as_chunks::<2>()
                .0
                .iter()
                .map(|p| i16::from_be_bytes(*p).max(0) as u16)
                .collect(),
        );
    }
    d.trace = Some(Trace {
        calls_version: version as u32,
        qualities: quality.bytes[..sequence.len()].to_vec(),
        peaks,
        channel_order: order,
        channels,
    });
    d.source = Some(Source {
        file_name: name.into(),
        format: "ABIF".into(),
        bytes: b.into(),
        retained_blocks: vec![],
    });
    d.description = "Imported Sanger sequencing trace".into();
    d.warnings.push(format!("Read-only ABIF trace: PBAS/PCON/PLOC{version} and analysed DATA9–12 are displayed. Original bytes are retained; supported instrument parameters are available in Parameters. Alternate calls, base editing, trimming and alignment are not implemented. Negative signal samples are clipped to zero for display."));
    d.validate()?;
    Ok(d)
}
/// Selected ABIF parameters decoded using the public Applied Biosystems schema.
/// Arrays and proprietary records remain in the original source bytes.
/// Tag-number meanings: Applied Biosystems ABIF File Format, July 2006,
/// pp. 24–28, 50–53 (including SPAC2 and SVER2 analysis-software records).
/// https://archive.gfjc.fiu.edu/workshops/resources/literature/ABIF_File_Format.pdf
pub(crate) fn source_metadata(bytes: &[u8]) -> Result<Vec<(String, String)>> {
    const DISPLAY_LIMIT: usize = 1024 * 1024;
    const LIMIT_LABEL: &str = "Display limit";
    const LIMIT_MESSAGE: &str = "Additional parameters remain in the retained original file (parameter display limited to 1 MiB).";
    let entries = directory(bytes)?;
    let channel_order = entries.get(&(*b"FWO_", 1)).and_then(parameter_value);
    let mut fields = vec![];
    // Include labels and the UI's line separators, reserving room for the notice.
    let mut displayed_bytes = LIMIT_LABEL.len() + LIMIT_MESSAGE.len() + 4;
    for ((tag, number), entry) in entries {
        let label = match &tag {
            b"SMPL" => "Sample name",
            b"MCHN" => "Instrument name / serial number",
            b"MODL" => "Instrument model",
            b"TUBE" => "Well",
            b"LANE" => "Lane",
            b"RunN" => "Run name",
            b"RUND" => match number {
                1 => "Run start date",
                2 => "Run stop date",
                3 => "Data collection start date",
                4 => "Data collection stop date",
                _ => "Date parameter",
            },
            b"RUNT" => match number {
                1 => "Run start time",
                2 => "Run stop time",
                3 => "Data collection start time",
                4 => "Data collection stop time",
                _ => "Time parameter",
            },
            b"CMNT" => "Comment",
            b"CTID" => "Container identifier",
            b"CTNM" => "Container name",
            b"LIMS" => "Sample tracking identifier",
            b"User" => "Plate user",
            b"DySN" => "Dye set name",
            b"DyeN" => "Dye name",
            b"DyeW" => "Dye wavelength",
            b"PDMF" => "Mobility file",
            b"MODF" | b"RMdN" => "Run module",
            b"RPrN" => "Run protocol",
            b"SVER" => match number {
                1 => "Data collection software version",
                2 => "Basecaller version",
                3 => "Data collection firmware version",
                4 => "Sample file format version",
                _ => "Version parameter",
            },
            b"SPAC" => match number {
                1 => "Average peak spacing used in last analysis",
                2 => "Basecaller name / BCP / DLL",
                3 => "Basecaller peak spacing",
                _ => "Spacing / basecaller parameter",
            },
            b"S/N%" => "Signal levels (channel order)",
            b"SCAN" => "Scan count",
            b"Tmpr" => "Run temperature setting (°C)",
            b"FWO_" => "Channel order",
            _ => continue,
        };
        let mut value = parameter_value(&entry)
            .unwrap_or_else(|| "Unavailable: unsupported or malformed parameter encoding".into());
        if tag == *b"S/N%"
            && entry.data_type == 4
            && entry.count == 4
            && entry.size == 2
            && let Some(order) = channel_order.as_ref().filter(|order| {
                order.len() == 4 && b"ACGT".iter().all(|base| order.as_bytes().contains(base))
            })
        {
            value = order
                .chars()
                .zip(entry.bytes.as_chunks::<2>().0.iter())
                .map(|(base, bytes)| format!("{base}={}", i16::from_be_bytes(*bytes)))
                .collect::<Vec<_>>()
                .join(", ");
        }
        let label = format!("{label} ({}{})", String::from_utf8_lossy(&tag), number);
        displayed_bytes += label.len() + value.len() + 4;
        if displayed_bytes > DISPLAY_LIMIT {
            fields.push((LIMIT_LABEL.into(), LIMIT_MESSAGE.into()));
            break;
        }
        fields.push((label, value));
    }
    Ok(fields)
}

fn parameter_value(entry: &Entry<'_>) -> Option<String> {
    let bytes = entry.bytes;
    if bytes.len() > 8192 {
        return Some("Value exceeds the 8 KiB parameter display limit; retained in source".into());
    }
    let text = |b: &[u8]| {
        std::str::from_utf8(b)
            .ok()
            .map(|s| s.trim_end_matches('\0').to_string())
    };
    match entry.data_type {
        2 if entry.size == 1 => text(bytes),
        18 if entry.size == 1 => {
            let (&length, rest) = bytes.split_first()?;
            if length as usize != rest.len() {
                return None;
            }
            text(rest)
        }
        19 if entry.size == 1 && bytes.last() == Some(&0) => text(&bytes[..bytes.len() - 1]),
        10 if entry.size == 4 && entry.count == 1 => {
            let year = u16::from_be_bytes(bytes[..2].try_into().ok()?);
            let month = bytes[2];
            let day = bytes[3];
            let leap =
                year.is_multiple_of(4) && (!year.is_multiple_of(100) || year.is_multiple_of(400));
            let days = match month {
                2 if leap => 29,
                2 => 28,
                4 | 6 | 9 | 11 => 30,
                1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
                _ => return None,
            };
            (year > 0 && day > 0 && day <= days).then(|| format!("{year:04}-{month:02}-{day:02}"))
        }
        11 if entry.size == 4 && entry.count == 1 => {
            (bytes[0] < 24 && bytes[1] < 60 && bytes[2] < 60 && bytes[3] < 100).then(|| {
                format!(
                    "{:02}:{:02}:{:02}.{:02} (instrument local time)",
                    bytes[0], bytes[1], bytes[2], bytes[3]
                )
            })
        }
        _ if entry.count <= 16 => {
            let values: Option<Vec<String>> = bytes
                .chunks(entry.size.max(1))
                .map(|chunk| match (entry.data_type, entry.size) {
                    (1, 1) => Some(chunk[0].to_string()),
                    (3, 2) => Some(u16::from_be_bytes(chunk.try_into().ok()?).to_string()),
                    (4, 2) => Some(i16::from_be_bytes(chunk.try_into().ok()?).to_string()),
                    (5, 4) => Some(i32::from_be_bytes(chunk.try_into().ok()?).to_string()),
                    (7, 4) => {
                        let v = f32::from_be_bytes(chunk.try_into().ok()?);
                        v.is_finite().then(|| v.to_string())
                    }
                    (8, 8) => {
                        let v = f64::from_be_bytes(chunk.try_into().ok()?);
                        v.is_finite().then(|| v.to_string())
                    }
                    _ => None,
                })
                .collect();
            values.map(|v| v.join(", "))
        }
        _ => None,
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    type Parameter<'a> = ([u8; 4], u32, u16, u16, &'a [u8]);

    fn parameter_file(parameters: &[Parameter<'_>]) -> Vec<u8> {
        let mut bytes = vec![0; 34 + parameters.len() * 28];
        bytes[..4].copy_from_slice(b"ABIF");
        bytes[18..22].copy_from_slice(&(parameters.len() as u32).to_be_bytes());
        bytes[26..30].copy_from_slice(&34_u32.to_be_bytes());
        for (i, (tag, number, data_type, size, value)) in parameters.iter().enumerate() {
            let p = 34 + i * 28;
            bytes[p..p + 4].copy_from_slice(tag);
            bytes[p + 4..p + 8].copy_from_slice(&number.to_be_bytes());
            bytes[p + 8..p + 10].copy_from_slice(&data_type.to_be_bytes());
            bytes[p + 10..p + 12].copy_from_slice(&size.to_be_bytes());
            bytes[p + 12..p + 16]
                .copy_from_slice(&((value.len() / *size as usize) as u32).to_be_bytes());
            bytes[p + 16..p + 20].copy_from_slice(&(value.len() as u32).to_be_bytes());
            if value.len() <= 4 {
                bytes[p + 20..p + 20 + value.len()].copy_from_slice(value);
            } else {
                let offset = bytes.len() as u32;
                bytes[p + 20..p + 24].copy_from_slice(&offset.to_be_bytes());
                bytes.extend_from_slice(value);
            }
        }
        bytes
    }

    #[test]
    fn numbered_parameter_meanings_do_not_confuse_basecaller_and_instrument() {
        let date = [7, 234, 10, 22];
        let time = [9, 30, 15, 0];
        let bytes = parameter_file(&[
            (*b"SPAC", 2, 19, 1, b"KB_3730_POP7_BDTv3\0"),
            (*b"SVER", 1, 19, 1, b"Collection v3\0"),
            (*b"SVER", 2, 19, 1, b"KB v1.4\0"),
            (*b"SVER", 3, 19, 1, b"Firmware v2\0"),
            (*b"SVER", 4, 19, 1, b"Format v1\0"),
            (*b"RUND", 2, 10, 4, &date),
            (*b"RUND", 3, 10, 4, &date),
            (*b"RUND", 4, 10, 4, &date),
            (*b"RUNT", 2, 11, 4, &time),
            (*b"RUNT", 3, 11, 4, &time),
            (*b"RUNT", 4, 11, 4, &time),
        ]);
        let fields = source_metadata(&bytes).unwrap();
        for (label, value) in [
            ("Basecaller name / BCP / DLL (SPAC2)", "KB_3730_POP7_BDTv3"),
            ("Data collection software version (SVER1)", "Collection v3"),
            ("Basecaller version (SVER2)", "KB v1.4"),
            ("Data collection firmware version (SVER3)", "Firmware v2"),
            ("Sample file format version (SVER4)", "Format v1"),
            ("Run stop date (RUND2)", "2026-10-22"),
            ("Data collection start date (RUND3)", "2026-10-22"),
            ("Data collection stop date (RUND4)", "2026-10-22"),
            (
                "Run stop time (RUNT2)",
                "09:30:15.00 (instrument local time)",
            ),
            (
                "Data collection start time (RUNT3)",
                "09:30:15.00 (instrument local time)",
            ),
            (
                "Data collection stop time (RUNT4)",
                "09:30:15.00 (instrument local time)",
            ),
        ] {
            assert!(fields.contains(&(label.into(), value.into())), "{label}");
        }
    }

    #[test]
    fn aggregate_parameter_display_is_bounded_and_reports_omission() {
        let value = vec![b'x'; 8192];
        let parameters: Vec<_> = (1..=150)
            .map(|number| (*b"CMNT", number, 2, 1, value.as_slice()))
            .collect();
        let bytes = parameter_file(&parameters);
        let fields = source_metadata(&bytes).unwrap();
        assert!(fields.len() < parameters.len());
        assert_eq!(fields[0], ("Comment (CMNT1)".into(), "x".repeat(8192)));
        assert_eq!(fields.last().unwrap().0, "Display limit");
        assert!(fields.last().unwrap().1.contains("retained original file"));
        let displayed = fields
            .into_iter()
            .map(|(label, value)| format!("{label}\n{value}"))
            .collect::<Vec<_>>()
            .join("\n\n");
        assert!(displayed.len() <= 1024 * 1024);
    }

    #[test]
    fn parameters_decode_strings_numbers_dates_and_channel_associations() {
        let bytes = include_bytes!("../../../validation/fixtures/regression/trace-parameters.ab1");
        let document = read("parameters.ab1", bytes).unwrap();
        let fields = crate::source_metadata(document.source.as_ref().unwrap()).unwrap();
        for (label, value) in [
            ("Sample name (SMPL1)", "PopGene synthetic sample"),
            ("Instrument model (MODL1)", "3730"),
            ("Instrument name / serial number (MCHN1)", "SYNTHETIC-0001"),
            ("Run start date (RUND1)", "2026-10-22"),
            (
                "Run start time (RUNT1)",
                "09:30:15.00 (instrument local time)",
            ),
            ("Lane (LANE1)", "1"),
            ("Average peak spacing used in last analysis (SPAC1)", "3"),
            (
                "Signal levels (channel order) (S/N%1)",
                "G=10, A=20, T=30, C=40",
            ),
        ] {
            assert!(
                fields.contains(&(label.into(), value.into())),
                "missing {label} = {value}"
            );
        }
        assert_eq!(&*document.sequence, "ACGTN");
        assert_eq!(document.trace.unwrap().qualities, vec![40, 35, 30, 20, 5]);
    }
    #[test]
    fn malformed_optional_parameters_are_not_fabricated() {
        let invalid_pascal = Entry {
            data_type: 18,
            size: 1,
            count: 3,
            bytes: &[9, b'A', b'B'],
        };
        assert!(parameter_value(&invalid_pascal).is_none());
        let invalid_date = Entry {
            data_type: 10,
            size: 4,
            count: 1,
            bytes: &[7, 234, 2, 30],
        };
        assert!(parameter_value(&invalid_date).is_none());
        let invalid_time = Entry {
            data_type: 11,
            size: 4,
            count: 1,
            bytes: &[25, 0, 0, 0],
        };
        assert!(parameter_value(&invalid_time).is_none());
        let infinite = Entry {
            data_type: 7,
            size: 4,
            count: 1,
            bytes: &f32::INFINITY.to_be_bytes(),
        };
        assert!(parameter_value(&infinite).is_none());
    }
    #[test]
    fn synthetic_trace_matches_known_arrays() {
        let d = read(
            "trace.ab1",
            include_bytes!("../../../validation/fixtures/synthetic-complete.ab1"),
        )
        .unwrap();
        assert_eq!(&*d.sequence, "ACGTN");
        let t = d.trace.unwrap();
        assert_eq!(t.calls_version, 1);
        assert_eq!(t.peaks, vec![2, 5, 8, 11, 14]);
        assert_eq!(t.qualities, vec![40, 35, 30, 20, 5]);
        assert_eq!(t.channel_order, "GATC");
        assert_eq!(t.channels[0][8], 1000);
    }
    #[test]
    fn original_call_fallback_and_bad_offsets() {
        let b = include_bytes!("../../../validation/fixtures/basecaller-only.ab1");
        assert_eq!(read("x", b).unwrap().trace.unwrap().calls_version, 2);
        let mut b = b.to_vec();
        b[26..30].copy_from_slice(&u32::MAX.to_be_bytes());
        assert!(read("bad", &b).is_err());
    }
}
