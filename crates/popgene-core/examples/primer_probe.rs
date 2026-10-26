//! Production Primer3 integration probe for validation/check_primer3.py.
use popgene_core::primer::*;
use std::io::{self, BufRead};
fn main() {
    for line in io::stdin().lock().lines() {
        let line = line.unwrap();
        let f: Vec<_> = line.split('\t').collect();
        if f[0] == "D" {
            let result = design_pairs(
                f[2],
                &DesignOptions {
                    target: Some(200..250),
                    ..DesignOptions::default()
                },
            )
            .unwrap();
            for (i, p) in result.pairs.iter().enumerate() {
                println!(
                    "D\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    f[1],
                    i,
                    p.forward.sequence,
                    p.forward.binding.start,
                    p.forward.binding.end,
                    p.reverse.sequence,
                    p.reverse.binding.start,
                    p.reverse.binding.end,
                    p.forward.tm_c,
                    p.reverse.tm_c,
                    p.product.len(),
                    p.penalty
                );
            }
            continue;
        }
        let v = |i: usize| f[i].parse::<f64>().unwrap();
        let reaction = Reaction {
            dna_nm: v(3),
            monovalent_mm: v(4),
            magnesium_mm: v(5),
            dntp_mm: v(6),
            dmso_percent: v(7),
            formamide_m: v(8),
            temperature_c: v(9),
            melting_model: if f[10] == "0" {
                MeltingModel::Breslauer
            } else {
                MeltingModel::SantaLucia
            },
            salt_correction: match f[11] {
                "0" => SaltCorrection::Schildkraut,
                "1" => SaltCorrection::SantaLucia,
                _ => SaltCorrection::Owczarzy,
            },
        };
        let a = analyze(f[2], reaction).unwrap();
        println!(
            "T\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            f[1],
            a.tm_c,
            a.hairpin.tm_c,
            a.hairpin.enthalpy_cal_mol,
            a.hairpin.entropy_cal_mol_k,
            a.hairpin.free_energy_cal_mol,
            a.homodimer.tm_c,
            a.homodimer.enthalpy_cal_mol,
            a.homodimer.entropy_cal_mol_k,
            a.homodimer.free_energy_cal_mol
        );
    }
}
