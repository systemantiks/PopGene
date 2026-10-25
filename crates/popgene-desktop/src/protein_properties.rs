//! Read-only protein composition and explicitly scoped sequence calculations.
use crate::feature_card::number;
use gpui_kit::component::{Sizable, StyledExt, button::*};
use gpui_kit::{prelude::*, *};
use popgene_core::protein::ProteinProperties;
use std::fmt::Write;

fn stat(label: &'static str, value: String) -> Div {
    div()
        .v_flex()
        .gap_1()
        .min_w(px(160.))
        .child(div().text_xs().text_color(rgb(0x96a6bc)).child(label))
        .child(div().text_lg().font_weight(FontWeight::MEDIUM).child(value))
}

fn report(properties: &ProteinProperties) -> String {
    let mut out = format!(
        "Protein properties\n{} sequence symbols; {} non-stop residues; {} stop symbols\n",
        properties.symbol_count, properties.residue_count, properties.stop_count
    );
    if let Some(mass) = properties.mass {
        let _ = writeln!(
            out,
            "Average neutral mass: {:.4} Da\nMonoisotopic neutral mass: {:.6} Da",
            mass.average_da, mass.monoisotopic_da
        );
    }
    if let Some(pi) = properties.theoretical_pi {
        let _ = writeln!(out, "Theoretical pI: {pi:.2} (Bjellqvist pKa set)");
    }
    if let Some(charge) = properties.charge_at_ph7 {
        let _ = writeln!(out, "Net charge at pH 7: {charge:.2} (Bjellqvist pKa set)");
    }
    if let Some(absorbance) = properties.absorbance {
        let _ = writeln!(
            out,
            "Extinction at 280 nm (reduced Cys): {} M^-1 cm^-1\nExtinction at 280 nm (all possible Cys pairs form cystines): {} M^-1 cm^-1\nA280 at 1 mg/mL, 1 cm (reduced Cys): {:.2}",
            absorbance.reduced_extinction,
            absorbance.paired_cysteine_extinction,
            absorbance.reduced_absorbance
        );
    }
    if let Some(reason) = &properties.unavailable_reason {
        let _ = writeln!(out, "{reason}");
    }
    out.push_str("\nComposition (percent of all stored sequence symbols, including stops)\nSymbol\tResidue\tCount\tPercent\n");
    for residue in &properties.composition {
        let _ = writeln!(
            out,
            "{}\t{}\t{}\t{:.2}%",
            residue.symbol, residue.name, residue.count, residue.percent
        );
    }
    out.push_str("\nModel: one linear unmodified chain; free N/C termini, reduced cysteines. Mass uses average/monoisotopic amino-acid masses minus one water per peptide bond. pI/charge uses Bjellqvist pKa values with residue-specific terminal corrections and charge balance over pH 0–14. Extinction uses Pace/Edelhoch Trp 5500, Tyr 1490, cystine 125; paired Cys is a separate hypothetical optical scenario. No disulfide bonds, modifications or processing are inferred for mass or pI.\nMethods: https://web.expasy.org/compute_pi/pi_tool-doc.html\nExtinction: https://web.expasy.org/protparam/protparam-doc.html\nMass tables: https://github.com/biopython/biopython/blob/biopython-188/Bio/Data/IUPACData.py\n");
    out
}

/// Pass a cached analysis for the current protein document revision.
pub(crate) fn render(properties: &ProteinProperties) -> AnyElement {
    let copy = report(properties);
    let mut summary = div()
        .h_flex()
        .flex_wrap()
        .gap_6()
        .child(stat("SEQUENCE SYMBOLS", number(properties.symbol_count)))
        .child(stat("RESIDUES", number(properties.residue_count)))
        .child(stat("STOP SYMBOLS", number(properties.stop_count)));
    if let Some(mass) = properties.mass {
        summary = summary
            .child(stat("AVERAGE MASS", format!("{:.2} Da", mass.average_da)))
            .child(stat(
                "MONOISOTOPIC MASS",
                format!("{:.4} Da", mass.monoisotopic_da),
            ));
    } else {
        summary = summary.child(stat("MOLECULAR MASS", "Unavailable".into()));
    }
    summary = summary.child(stat(
        "THEORETICAL pI · BJELLQVIST",
        properties
            .theoretical_pi
            .map_or_else(|| "Unavailable".into(), |pi| format!("{pi:.2}")),
    ));
    summary = summary.child(stat(
        "CHARGE AT pH 7 · BJELLQVIST",
        properties
            .charge_at_ph7
            .map_or_else(|| "Unavailable".into(), |charge| format!("{charge:+.2}")),
    ));
    if let Some(absorbance) = properties.absorbance {
        summary = summary
            .child(stat(
                "EXTINCTION 280 nm · REDUCED",
                format!("{} M⁻¹cm⁻¹", number(absorbance.reduced_extinction)),
            ))
            .child(stat(
                "EXTINCTION · CYS PAIRS",
                format!("{} M⁻¹cm⁻¹", number(absorbance.paired_cysteine_extinction)),
            ))
            .child(stat(
                "A280 · 1 mg/mL · 1 cm",
                format!("{:.2}", absorbance.reduced_absorbance),
            ));
    } else {
        summary = summary
            .child(stat("EXTINCTION 280 nm", "Unavailable".into()))
            .child(stat("A280 · 1 mg/mL · 1 cm", "Unavailable".into()));
    }
    div().id("protein-properties-scroll").flex_1().min_h_0().overflow_y_scroll().px_6().pb_6()
        .child(div().v_flex().gap_4().p_5().rounded_lg().border_1().border_color(rgb(0x293a4a)).bg(rgb(0x172231))
            .child(div().h_flex().gap_3()
                .child(div().text_lg().font_weight(FontWeight::SEMIBOLD).child("Protein properties"))
                .child(div().flex_1())
                .child(Button::new("copy-protein-properties").small().label("Copy properties")
                    .on_click(move |_, _, cx| cx.write_to_clipboard(ClipboardItem::new_string(copy.clone())))))
            .child(summary)
            .when_some(properties.unavailable_reason.as_ref(), |d, reason| d.child(div().text_sm().text_color(rgb(0xe7bd7b)).child(reason.clone())))
            .child(div().text_sm().text_color(rgb(0x96a6bc)).child("One linear, unmodified chain with free N/C termini and reduced cysteines. Mass and pI require the 20 standard amino acids; ambiguous residues and stop symbols remain in the composition below."))
            .child(div().text_sm().text_color(rgb(0x96a6bc)).child("Theoretical pI and charge use the Bjellqvist pKa set and residue-specific terminal corrections. These estimates can differ from other software and experiment. Modifications, disulfide bonds and protein processing are not inferred for mass or pI."))
            .child(div().text_sm().text_color(rgb(0x96a6bc)).child("Extinction uses Pace/Edelhoch coefficients. The Cys-pairs value assumes every possible cysteine pair forms a cystine; A280 uses reduced cysteines. No other chromophores are modeled."))
            .child(div().h_flex().flex_wrap().gap_2()
                .child(Button::new("protein-methods").small().ghost().label("Method reference")
                    .on_click(|_, _, cx| cx.open_url("https://web.expasy.org/compute_pi/pi_tool-doc.html")))
                .child(Button::new("protein-mass-tables").small().ghost().label("Mass table reference")
                    .on_click(|_, _, cx| cx.open_url("https://github.com/biopython/biopython/blob/biopython-188/Bio/Data/IUPACData.py")))
                .child(Button::new("protein-optics-reference").small().ghost().label("Extinction reference")
                    .on_click(|_, _, cx| cx.open_url("https://web.expasy.org/protparam/protparam-doc.html")))))
        .child(div().pt_5().pb_3().v_flex().gap_1()
            .child(div().text_lg().font_weight(FontWeight::MEDIUM).child("Amino acid composition"))
            .child(div().text_sm().text_color(rgb(0x96a6bc)).child("Percentages use the full sequence length, including stop symbols. Stop symbols are counted separately above.")))
        .child(div().h_flex().gap_4().px_4().py_2().text_xs().text_color(rgb(0x96a6bc))
            .child(div().w(px(54.)).child("CODE"))
            .child(div().flex_1().child("RESIDUE"))
            .child(div().w(px(100.)).text_right().child("COUNT"))
            .child(div().w(px(110.)).text_right().child("PERCENT")))
        .children(properties.composition.iter().map(|residue| div().h_flex().gap_4().px_4().py_2().border_b_1().border_color(rgb(0x293a4a))
            .child(div().w(px(54.)).font_family("monospace").child(residue.symbol.to_string()))
            .child(div().flex_1().child(residue.name))
            .child(div().w(px(100.)).text_right().child(number(residue.count)))
            .child(div().w(px(110.)).text_right().child(format!("{:.2}%", residue.percent)))))
        .into_any_element()
}
