//! Portable, self-contained vector figures and bounded rasterization.
use std::collections::BTreeSet;
use std::sync::{Arc, Mutex, OnceLock};
mod fonts;
mod pdf;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FigureFormat {
    Svg,
    Pdf,
    Png,
}
impl FigureFormat {
    pub fn label(self) -> &'static str {
        match self {
            Self::Svg => "SVG vector",
            Self::Pdf => "PDF vector",
            Self::Png => "PNG image",
        }
    }
    pub fn extension(self) -> &'static str {
        match self {
            Self::Svg => "svg",
            Self::Pdf => "pdf",
            Self::Png => "png",
        }
    }
}
static FONTS: OnceLock<Arc<svg2pdf::usvg::fontdb::Database>> = OnceLock::new();
/// Generated SVG contains only inline vector paths and escaped text. External
/// images/styles are never resolved. Text is outlined in exported SVG files.
pub fn encode(svg: &str, format: FigureFormat, scale: f32) -> Result<Vec<u8>, String> {
    encode_with_pages(svg, format, scale, None)
}
pub fn encode_with_pages(
    svg: &str,
    format: FigureFormat,
    scale: f32,
    page_height: Option<u32>,
) -> Result<Vec<u8>, String> {
    let fonts = FONTS.get_or_init(|| Arc::new(fonts::database(true)));
    encode_with_fonts(svg, format, scale, page_height, fonts.clone())
}

fn encode_with_fonts(
    svg: &str,
    format: FigureFormat,
    scale: f32,
    page_height: Option<u32>,
    fonts: Arc<svg2pdf::usvg::fontdb::Database>,
) -> Result<Vec<u8>, String> {
    if svg.len() > 32 * 1024 * 1024 {
        return Err("The figure exceeds the 32 MiB vector limit".into());
    }
    let missing = Mutex::new(BTreeSet::new());
    let fallback = svg2pdf::usvg::FontResolver::default_fallback_selector();
    let options = svg2pdf::usvg::Options {
        fontdb: fonts,
        font_family: fonts::DEFAULT_FAMILY.into(),
        font_resolver: svg2pdf::usvg::FontResolver {
            select_font: svg2pdf::usvg::FontResolver::default_font_selector(),
            select_fallback: Box::new(|character, excluded, database| {
                let selected = fallback(character, excluded, database);
                if selected.is_none() {
                    missing.lock().unwrap().insert(character);
                }
                selected
            }),
        },
        ..Default::default()
    };
    let tree = svg2pdf::usvg::Tree::from_str(svg, &options)
        .map_err(|e| format!("Figure rendering failed: {e}"))?;
    let missing = missing.lock().unwrap();
    if !missing.is_empty() {
        let codepoints = missing
            .iter()
            .take(8)
            .map(|c| format!("U+{:04X}", *c as u32))
            .collect::<Vec<_>>()
            .join(", ");
        return Err(format!(
            "Figure labels contain characters unavailable in bundled or installed fonts: {codepoints}. Install a font covering these characters and restart PopGene."
        ));
    }
    drop(missing);
    // A fallback face may contain the requested character yet fail to combine
    // with the original shaping (for example, a ligature changes glyph count).
    // Inspect the finished glyphs too, instead of accepting .notdef as a label.
    let mut unresolved = BTreeSet::new();
    let mut groups = vec![tree.root()];
    while let Some(group) = groups.pop() {
        for node in group.children() {
            match node {
                svg2pdf::usvg::Node::Group(group) => groups.push(group),
                svg2pdf::usvg::Node::Text(text) => {
                    for span in text.layouted().iter().filter(|span| span.visible) {
                        for glyph in &span.positioned_glyphs {
                            if glyph.id.0 == 0 {
                                unresolved.extend(glyph.text.chars());
                            }
                        }
                    }
                }
                _ => {}
            }
        }
    }
    if !unresolved.is_empty() {
        let codepoints = unresolved
            .iter()
            .take(8)
            .map(|c| format!("U+{:04X}", *c as u32))
            .collect::<Vec<_>>()
            .join(", ");
        return Err(format!(
            "Figure labels could not be shaped completely: {codepoints}."
        ));
    }
    match format {
        FigureFormat::Svg => {
            let mut output = tree.to_string(&Default::default());
            // The renderer outlines text and drops non-rendering metadata. Retain
            // the source title/description for accessible, self-described figures.
            let source =
                svg2pdf::usvg::roxmltree::Document::parse(svg).map_err(|e| e.to_string())?;
            let mut metadata = String::new();
            for node in source
                .root_element()
                .children()
                .filter(|n| n.is_element() && matches!(n.tag_name().name(), "title" | "desc"))
            {
                let tag = node.tag_name().name();
                let value = node
                    .descendants()
                    .filter(|n| n.is_text())
                    .filter_map(|n| n.text())
                    .collect::<String>();
                let value = value
                    .replace('&', "&amp;")
                    .replace('<', "&lt;")
                    .replace('>', "&gt;");
                metadata.push_str(&format!("<{tag}>{value}</{tag}>"));
            }
            if let Some(start) = output.find("<svg")
                && let Some(end) = output[start..].find('>')
            {
                output.insert_str(start + end + 1, &metadata);
            }
            Ok(output.into_bytes())
        }
        FigureFormat::Pdf => {
            if let Some(height) = page_height {
                return pdf::encode(&tree, svg, height, &options);
            }
            // Acrobat default-unit page bounds: Adobe PDF1.7 implementation note177.
            // https://opensource.adobe.com/dc-acrobat-sdk-docs/pdfstandards/pdfreference1.7old.pdf
            if tree.size().width() > 19_200. || tree.size().height() > 19_200. {
                return Err(
                    "This figure exceeds the PDF page limit; enable sequence pagination".into(),
                );
            }
            svg2pdf::to_pdf(&tree, Default::default(), svg2pdf::PageOptions { dpi: 96. })
                .map_err(|e| format!("PDF export failed: {e}"))
        }
        FigureFormat::Png => {
            if !scale.is_finite() || !(0.25..=8.).contains(&scale) {
                return Err("Raster scale must be between 0.25 and 8".into());
            }
            let width = (tree.size().width() * scale).ceil() as u32;
            let height = (tree.size().height() * scale).ceil() as u32;
            if u64::from(width) * u64::from(height) > 64_000_000 {
                return Err(
                    "PNG exceeds 64 million pixels. Reduce dimensions or raster scale.".into(),
                );
            }
            let mut pixmap = resvg::tiny_skia::Pixmap::new(width, height)
                .ok_or("Could not allocate the figure raster")?;
            resvg::render(
                &tree,
                resvg::tiny_skia::Transform::from_scale(scale, scale),
                &mut pixmap.as_mut(),
            );
            let mut png = pixmap
                .encode_png()
                .map_err(|e| format!("PNG export failed: {e}"))?;
            // Standard PNG physical-pixel density, matching the scaled 96 dpi layout.
            let ppm = (96. * scale / 0.0254).round() as u32;
            let mut chunk = vec![];
            chunk.extend_from_slice(&9_u32.to_be_bytes());
            chunk.extend_from_slice(b"pHYs");
            chunk.extend_from_slice(&ppm.to_be_bytes());
            chunk.extend_from_slice(&ppm.to_be_bytes());
            chunk.push(1);
            let crc = crc32fast::hash(&chunk[4..]);
            chunk.extend_from_slice(&crc.to_be_bytes());
            png.splice(33..33, chunk);
            Ok(png)
        }
    }
}
#[cfg(test)]
mod tests {
    use super::{FigureFormat, encode, encode_with_fonts, fonts};
    use std::sync::Arc;

    #[test]
    fn bundled_fonts_render_labels_without_any_system_fonts() {
        let database = Arc::new(fonts::database(false));
        assert_eq!(database.faces().count(), 2);
        let mono = database
            .query(&svg2pdf::usvg::fontdb::Query {
                families: &[svg2pdf::usvg::fontdb::Family::Monospace],
                ..Default::default()
            })
            .unwrap();
        assert!(database.face(mono).unwrap().monospaced);

        // Text is the only ink in these inputs: unlike the original regression,
        // geometric paths cannot conceal missing labels. Named-family fallbacks,
        // generic monospace and each scientific symbol are checked separately.
        let examples = [
            (
                "Arial, sans-serif",
                "normal",
                "04-methylation-clear EcoRI 11",
            ),
            ("monospace", "normal", "ACGTNN 5′ → 3′"),
            ("No such installed family", "bold", "α DNA 5′ → 3′"),
            ("sans-serif", "normal", "α"),
            ("sans-serif", "normal", "′"),
            ("sans-serif", "normal", "→"),
        ];
        for (index, (family, weight, label)) in examples.iter().enumerate() {
            let source = format!(
                "<svg xmlns='http://www.w3.org/2000/svg' width='440' height='60'><rect width='440' height='60' fill='white'/><text x='10' y='38' font-family='{family}' font-weight='{weight}' font-size='24'>{label}</text></svg>"
            );
            let render =
                |format| encode_with_fonts(&source, format, 1., None, database.clone()).unwrap();
            let svg = String::from_utf8(render(FigureFormat::Svg)).unwrap();
            assert!(svg.contains("<path"), "No glyph outlines for {label}");
            assert!(!svg.contains("<text"));
            let png = render(FigureFormat::Png);
            let pixels = resvg::tiny_skia::Pixmap::decode_png(&png).unwrap();
            let ink = pixels
                .pixels()
                .iter()
                .filter(|p| p.red() < 200 && p.green() < 200 && p.blue() < 200)
                .count();
            assert!(ink > 0, "No rasterized glyphs for {label}");
            let pdf = render(FigureFormat::Pdf);
            assert!(pdf.windows(11).any(|s| s == b"/Type /Font"));
            assert!(pdf.windows(10).any(|s| s == b"/ToUnicode"));
            if let Ok(dir) = std::env::var("POPGENE_FIGURE_QA_DIR") {
                let dir = std::path::PathBuf::from(dir);
                std::fs::create_dir_all(&dir).unwrap();
                std::fs::write(dir.join(format!("bundled-fonts-{index}.pdf")), pdf).unwrap();
                std::fs::write(dir.join(format!("bundled-fonts-{index}.png")), png).unwrap();
            }
        }
    }

    #[test]
    fn bundled_fonts_survive_paginated_reparsing_and_report_missing_glyphs() {
        let database = Arc::new(fonts::database(false));
        let source = "<svg xmlns='http://www.w3.org/2000/svg' width='440' height='600'><g font-family='monospace'><g data-popgene-row-start='0' data-popgene-row-end='300'><text x='10' y='30'>ACGT α 5′ → 3′</text></g><g data-popgene-row-start='300' data-popgene-row-end='600'><text x='10' y='330'>TTAA EcoRI 11</text></g></g></svg>";
        let pdf =
            encode_with_fonts(source, FigureFormat::Pdf, 1., Some(400), database.clone()).unwrap();
        // Page numbering uses a standard PDF font; a ToUnicode map proves the
        // SVG labels survived the separate page parsing and font embedding too.
        assert!(pdf.windows(10).any(|s| s == b"/ToUnicode"));
        if let Ok(dir) = std::env::var("POPGENE_FIGURE_QA_DIR") {
            let dir = std::path::PathBuf::from(dir);
            std::fs::create_dir_all(&dir).unwrap();
            std::fs::write(dir.join("bundled-fonts-paginated.pdf"), pdf).unwrap();
        }
        let unsupported = "<svg xmlns='http://www.w3.org/2000/svg' width='80' height='60'><text x='10' y='38'>\u{2a6df}</text></svg>";
        for format in [FigureFormat::Svg, FigureFormat::Pdf, FigureFormat::Png] {
            let error =
                encode_with_fonts(unsupported, format, 1., None, database.clone()).unwrap_err();
            assert!(error.contains("U+2A6DF"), "{error}");
        }
    }

    #[test]
    fn incomplete_mixed_font_shaping_does_not_export_missing_glyphs() {
        let database = Arc::new(fonts::database(false));
        // Noto Sans forms the fi ligature and supplies long-s / barred-b. Fira
        // Mono supplies the arrow, but the differing glyph count makes usvg
        // abandon this fallback, despite select_fallback finding a valid face.
        for label in ["fi ſ →", "fi ƀ →"] {
            let source = format!(
                "<svg xmlns='http://www.w3.org/2000/svg' width='200' height='60'><text x='10' y='38' font-family='Noto Sans'>{label}</text></svg>"
            );
            for format in [FigureFormat::Svg, FigureFormat::Pdf, FigureFormat::Png] {
                let error =
                    encode_with_fonts(&source, format, 1., None, database.clone()).unwrap_err();
                assert!(error.contains("could not be shaped completely"), "{error}");
                assert!(error.contains("U+2192"), "{error}");
            }
        }
    }
    #[test]
    fn vector_and_raster_outputs_are_portable_and_bounded() {
        let svg = "<svg xmlns='http://www.w3.org/2000/svg' width='400' height='300'><title>DNA &amp; labels</title><desc>Original positions retained.</desc><rect width='400' height='300' fill='white'/><path d='M20 40L300 40' stroke='black'/><text x='20' y='30' font-family='Arial' font-size='18'>α DNA 5′ → 3′</text></svg>";
        let vector = String::from_utf8(encode(svg, FigureFormat::Svg, 1.).unwrap()).unwrap();
        assert!(!vector.contains("<text"));
        assert!(vector.contains("<path"));
        assert!(vector.contains("<title>DNA &amp; labels</title>"));
        assert!(vector.contains("<desc>Original positions retained.</desc>"));
        let pdf = encode(svg, FigureFormat::Pdf, 1.).unwrap();
        assert!(pdf.starts_with(b"%PDF-"));
        let png = encode(svg, FigureFormat::Png, 2.).unwrap();
        assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n");
        assert_eq!(u32::from_be_bytes(png[16..20].try_into().unwrap()), 800);
        assert_eq!(u32::from_be_bytes(png[20..24].try_into().unwrap()), 600);
        assert_eq!(&png[37..41], b"pHYs");
        assert_eq!(u32::from_be_bytes(png[41..45].try_into().unwrap()), 7559);
        assert_eq!(u32::from_be_bytes(png[45..49].try_into().unwrap()), 7559);
        assert_eq!(png[49], 1);
        assert_eq!(
            u32::from_be_bytes(png[50..54].try_into().unwrap()),
            crc32fast::hash(&png[37..50])
        );
        assert!(encode(svg, FigureFormat::Png, f32::NAN).is_err());
        let huge = svg.replace("width='400' height='300'", "width='8000' height='8000'");
        assert!(encode(&huge, FigureFormat::Png, 2.).is_err());
        if let Ok(dir) = std::env::var("POPGENE_FIGURE_QA_DIR") {
            let dir = std::path::PathBuf::from(dir);
            std::fs::create_dir_all(&dir).unwrap();
            std::fs::write(dir.join("typography.svg"), vector).unwrap();
            std::fs::write(dir.join("typography.pdf"), pdf).unwrap();
            std::fs::write(dir.join("typography.png"), png).unwrap();
        }
    }
}
