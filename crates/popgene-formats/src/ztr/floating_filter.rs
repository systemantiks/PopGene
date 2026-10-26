//! Legacy method73, using explicit IEEE754 binary32 intermediates and unfused
//! arithmetic and first-word byte order from Staden's disabled little-endian
//! reference encoder. Coefficients are fixed bits,
//! eliminating host libm variation. Original encoders with excess precision may
//! disagree; the importer always reports this legacy precision assumption.
//! BSD notices: vendor-notices/staden-io-lib/LICENSE.txt.
use super::{Result, bad};

pub(super) fn decode(data: &[u8]) -> Result<Vec<u8>> {
    if data.len() < 2 || !data.len().is_multiple_of(2) || data[1] != 0 {
        return Err(bad("Invalid floating Chebyshev header or word alignment"));
    }
    let frac = [0x3f6d34a8, 0x3ec36b48, 0x3f000000, 0x3f1e4a5c, 0x3d965ac0].map(f32::from_bits);
    let weights = [
        [0x3ecccccd, 0x3ecccccd, 0x3ecccccd, 0x3ecccccd, 0x3ecccccd],
        [0x3ec2c6c0, 0x3e70c1c0, 0x23e1e84f, 0xbe70c1c0, 0xbec2c6c0],
        [0x3ea5afca, 0xbdfd2590, 0xbecccccd, 0xbdfd2590, 0x3ea5afca],
        [0x3e70c1c0, 0xbec2c6c0, 0xa4a96e3c, 0x3ec2c6c0, 0xbe70c1c0],
    ]
    .map(|row| row.map(f32::from_bits));
    let mut decoded = Vec::<u16>::with_capacity(data.len() / 2 - 1);
    for (index, bytes) in data[2..].as_chunks::<2>().0.iter().enumerate() {
        // The disabled little-endian Staden encoder stores its first
        // word in host order; subsequent residuals are big endian. Signal
        // chunks normally start with the invariant zero/padding word.
        let residual = if index == 0 {
            u16::from_le_bytes(*bytes)
        } else {
            u16::from_be_bytes(*bytes)
        };
        let predicted = if index < 4 {
            i32::from(decoded.last().copied().unwrap_or(0))
        } else {
            let d = &decoded[index - 4..];
            let d = [d[0] as f32, d[1] as f32, d[2] as f32, d[3] as f32];
            let f = [
                d[2] * frac[4] + d[3] * frac[0],
                d[2] * frac[3] + d[3] * frac[1],
                d[1] * frac[2] + d[2] * frac[2],
                d[0] * frac[1] + d[1] * frac[3],
                d[0] * frac[0] + d[1] * frac[4],
            ];
            let c = weights
                .map(|w| f[0] * w[0] + f[1] * w[1] + f[2] * w[2] + f[3] * w[3] + f[4] * w[4]);
            let y = (10.0_f64 / 3.0) as f32;
            let dd = y * c[3] + c[2];
            // C promotes only the outer expression to double; preserve its
            // single-precision inner expression and truncation toward zero.
            let prediction = 0.5 + (5.0_f64 / 3.0) * f64::from(y * dd - c[3] + c[1])
                - f64::from(dd)
                + f64::from(c[0]) / 2.;
            if !prediction.is_finite()
                || prediction > i32::MAX as f64
                || prediction < i32::MIN as f64
            {
                return Err(bad("Floating Chebyshev predictor overflow"));
            }
            (prediction as i32).max(0)
        };
        decoded.push((predicted + i32::from(residual as i16)) as u16);
    }
    Ok(decoded.into_iter().flat_map(u16::to_be_bytes).collect())
}

#[cfg(test)]
mod tests {
    use super::decode;
    #[test]
    fn original_c_encoder_vectors_and_legacy_precision_warning() {
        fn hex(value: &str) -> Vec<u8> {
            value
                .as_bytes()
                .as_chunks::<2>()
                .0
                .iter()
                .map(|pair| u8::from_str_radix(std::str::from_utf8(pair).unwrap(), 16).unwrap())
                .collect()
        }
        let facts: serde_json::Value = serde_json::from_str(include_str!(
            "../../../../validation/fixtures/trace/floating73/vectors.json"
        ))
        .unwrap();
        for item in facts["vectors"].as_array().unwrap() {
            assert_eq!(
                decode(&hex(item["input"].as_str().unwrap())).unwrap(),
                hex(item["expected"].as_str().unwrap()),
                "{}",
                item["name"]
            );
        }
        assert!(decode(&[73]).is_err());
        assert!(decode(&[73, 0, 0]).is_err());
        assert!(decode(&[73, 1]).is_err());
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../validation/fixtures/trace/floating73/legacy.ztr");
        let document = crate::import("legacy.ztr", &std::fs::read(path).unwrap()).unwrap();
        assert!(document.warnings.iter().any(|w| w.contains("IEEE754")));
        assert_eq!(
            document.trace.unwrap().channels,
            vec![
                vec![0, 10, 100, 200, 300, 400],
                vec![1, 2, 3, 4, 5, 6],
                vec![600, 500, 400, 300, 200, 100],
                vec![20, 40, 60, 80, 100, 120]
            ]
        );
    }
}
