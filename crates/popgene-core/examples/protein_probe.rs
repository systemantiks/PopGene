//! Production protein-properties probe for the independent validation script.
use popgene_core::protein;

use std::io::{self, BufRead};
fn n(v: Option<f64>) -> String {
    v.map_or("-".into(), |v| format!("{v:.12}"))
}
fn main() {
    for line in io::stdin().lock().lines() {
        let sequence = line.unwrap();
        if std::env::args().any(|arg| arg == "--detailed") {
            use protein::{Modification, PkaModel, ProteinOptions};
            let fields = sequence.split('\t').collect::<Vec<_>>();
            let modifications = fields[5]
                .split(',')
                .filter(|s| !s.is_empty())
                .map(|m| {
                    let f = m.split(':').collect::<Vec<_>>();
                    match f[0] {
                        "n" => Modification::NAcetyl,
                        "c" => Modification::CAmide,
                        "k" => Modification::LysAcetyl(f[1].parse().unwrap()),
                        "m" => Modification::MetOxidation(f[1].parse().unwrap()),
                        "p" => Modification::Phosphorylation(f[1].parse().unwrap()),
                        "d" => {
                            Modification::Disulfide(f[1].parse().unwrap(), f[2].parse().unwrap())
                        }
                        _ => panic!("invalid modification"),
                    }
                })
                .collect();
            let p = protein::analyze_with_options(
                fields[0],
                &ProteinOptions {
                    range: Some(fields[1].parse::<usize>().unwrap()..fields[2].parse().unwrap()),
                    model: match fields[3] {
                        "bjellqvist" => PkaModel::Bjellqvist,
                        "emboss" => PkaModel::Emboss,
                        "grimsley2009" => PkaModel::Grimsley2009,
                        _ => panic!("invalid charge model"),
                    },
                    ph: fields[4].parse().unwrap(),
                    modifications,
                    secondary_structure: fields[6] == "1",
                },
            )
            .unwrap();
            println!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                n(p.mass.map(|m| m.average_da)),
                n(p.mass.map(|m| m.monoisotopic_da)),
                n(p.theoretical_pi),
                n(p.charge),
                n(p.extinction_280.map(|e| e as f64)),
                n(p.absorbance_1mg_ml),
                p.secondary.map_or("-".into(), |s| s.assignments)
            );
            continue;
        }
        let p = protein::analyze(&sequence);
        let c = p
            .composition
            .iter()
            .map(|r| format!("{}:{}:{:.12}", r.symbol, r.count, r.percent))
            .collect::<Vec<_>>()
            .join(",");
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            n(p.mass.map(|m| m.average_da)),
            n(p.mass.map(|m| m.monoisotopic_da)),
            n(p.theoretical_pi),
            n(p.charge_at_ph7),
            n(p.absorbance.map(|a| a.reduced_extinction as f64)),
            n(p.absorbance.map(|a| a.paired_cysteine_extinction as f64)),
            p.residue_count,
            p.stop_count,
            c
        );
    }
}
