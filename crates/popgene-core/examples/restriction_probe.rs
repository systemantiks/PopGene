//! Line-oriented production-engine probe used by validation/check_restriction.py.
use popgene_core::restriction::{self, ENZYMES};
use std::io::{self, BufRead};
fn cut(value: Option<usize>) -> String {
    value.map_or("-".into(), |n| n.to_string())
}
fn main() {
    for e in ENZYMES {
        println!(
            "E\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            e.name,
            e.recognition,
            e.cleavage_known,
            e.top_cut,
            e.bottom_cut,
            e.second_cut.map_or("-".into(), |p| p.0.to_string()),
            e.second_cut.map_or("-".into(), |p| p.1.to_string())
        );
    }
    for line in io::stdin().lock().lines() {
        let line = line.unwrap();
        let parts: Vec<_> = line.splitn(4, '\t').collect();
        let indices = if parts[2] == "*" {
            (0..ENZYMES.len()).collect::<Vec<_>>()
        } else {
            vec![ENZYMES.iter().position(|e| e.name == parts[2]).unwrap()]
        };
        let start = std::time::Instant::now();
        let result =
            restriction::find_sites_for(parts[3], parts[1] == "circular", &indices).unwrap();
        println!(
            "T\t{}\t{}\t{}",
            parts[0],
            start.elapsed().as_micros(),
            result.truncated
        );
        for i in indices {
            println!("C\t{}\t{}\t{}", parts[0], ENZYMES[i].name, result.counts[i]);
        }
        for s in result.sites {
            println!(
                "S\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                parts[0],
                ENZYMES[s.enzyme_index].name,
                s.start,
                s.reverse,
                s.recognition_length,
                cut(s.top_cut),
                cut(s.bottom_cut),
                cut(s.second_top_cut),
                cut(s.second_bottom_cut)
            );
        }
    }
}
