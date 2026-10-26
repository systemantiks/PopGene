//! SVG serialization of the same vector geometry used by the native map.
use super::{ACCENT, Draw, INK, MUTED, MapData, MapOptions, MapRestrictionSite, Rect, build_scene};
use gpui_kit::{Hsla, Rgba, rgb};
use popgene_core::Document;
use std::fmt::Write;

pub(crate) fn xml(text: &str) -> String {
    text.chars()
        .filter(|c| !c.is_control() || matches!(c, '\n' | '\r' | '\t'))
        .map(|c| match c {
            '&' => "&amp;".into(),
            '<' => "&lt;".into(),
            '>' => "&gt;".into(),
            '"' => "&quot;".into(),
            '\'' => "&apos;".into(),
            _ => c.to_string(),
        })
        .collect()
}
fn color(mut color: Hsla, light: bool, text: bool) -> (String, f32) {
    if light {
        for (dark, light) in [(INK, 0x172b3b), (MUTED, 0x536579), (ACCENT, 0x087466)] {
            let reference: Hsla = rgb(dark).into();
            if color == reference {
                color = rgb(light).into();
                break;
            }
        }
    }
    if light && text {
        color.l = color.l.min(0.40);
    }
    let color: Rgba = color.into();
    (
        format!(
            "#{:02x}{:02x}{:02x}",
            (color.r * 255.).round() as u8,
            (color.g * 255.).round() as u8,
            (color.b * 255.).round() as u8
        ),
        color.a,
    )
}
pub(crate) fn map_svg(
    document: &Document,
    mut options: MapOptions,
    sites: &[MapRestrictionSite],
    width: u32,
    height: u32,
    light: bool,
    measure: impl Fn(&str, f32) -> f32,
) -> Result<String, String> {
    document.validate().map_err(|e| e.to_string())?;
    if !(400..=8000).contains(&width) || !(300..=8000).contains(&height) {
        return Err("Figure dimensions must be 400–8,000 by 300–8,000 pixels".into());
    }
    options.zoom = 1.;
    let scene = build_scene(
        Rect {
            x: 0.,
            y: 0.,
            w: width as f32,
            h: height as f32,
        },
        &MapData {
            doc: document,
            selected: &[],
            options,
            sites,
            hovered: None,
        },
        measure,
    );
    let mut svg = format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{width}\" height=\"{height}\" viewBox=\"0 0 {width} {height}\"><title>{}</title><desc>{} {}. Coordinates are 1-based inclusive; restriction markers indicate cut boundaries. {}</desc><rect width=\"100%\" height=\"100%\" fill=\"{}\"/><g font-family=\"Arial, sans-serif\">",
        xml(&document.name),
        document.sequence.len(),
        document.unit(),
        xml(&scene.notices.join("; ")),
        if light { "#ffffff" } else { "#111a27" }
    );
    for draw in scene.draws {
        match draw {
            Draw::Stroke(points, width, ink) => {
                let (ink, opacity) = color(ink, light, false);
                write!(svg,"<polyline fill=\"none\" stroke=\"{ink}\" stroke-opacity=\"{opacity:.4}\" stroke-width=\"{width:.3}\" stroke-linejoin=\"round\" stroke-linecap=\"round\" points=\"").unwrap();
                for p in points {
                    write!(svg, "{:.3},{:.3} ", p.x, p.y).unwrap();
                }
                svg.push_str("\"/>");
            }
            Draw::Fill(points, ink) => {
                let (ink, opacity) = color(ink, light, false);
                write!(
                    svg,
                    "<polygon fill=\"{ink}\" fill-opacity=\"{opacity:.4}\" points=\""
                )
                .unwrap();
                for p in points {
                    write!(svg, "{:.3},{:.3} ", p.x, p.y).unwrap();
                }
                svg.push_str("\"/>");
            }
            Draw::Text {
                point,
                text,
                size,
                color: ink,
            } => {
                let (ink, opacity) = color(ink, light, true);
                write!(svg,"<text x=\"{:.3}\" y=\"{:.3}\" font-size=\"{size:.3}\" fill=\"{ink}\" fill-opacity=\"{opacity:.4}\">{}</text>",point.x,point.y+size,xml(&text)).unwrap();
            }
        }
    }
    svg.push_str("</g></svg>");
    Ok(svg)
}

#[cfg(test)]
mod tests {
    use super::*;
    use popgene_core::Segment;
    #[test]
    fn export_retains_geometry_and_escapes_source_labels() {
        let mut doc = Document::new("DNA <test>&", &"ACGT".repeat(100)).unwrap();
        doc.circular = true;
        doc.annotate(
            "marker <script>&".into(),
            Segment {
                start: 20,
                end: 180,
            },
            true,
        )
        .unwrap();
        let original = doc.clone();
        let svg = map_svg(
            &doc,
            MapOptions::default(),
            &[],
            1000,
            800,
            true,
            |s, size| s.chars().count() as f32 * size * 0.58,
        )
        .unwrap();
        let parsed = svg2pdf::usvg::Tree::from_str(&svg, &Default::default()).unwrap();
        assert_eq!(parsed.size().width(), 1000.);
        assert!(svg.contains("DNA &lt;test&gt;&amp;"));
        assert!(svg.contains("marker &lt;script&gt;&amp;"));
        assert!(svg.contains("<polygon"));
        assert!(!svg.contains("<script>"));
        assert_eq!(doc, original);
    }
    #[test]
    fn styled_fixture_figures_render_both_topologies() {
        let mut doc = popgene_formats::import(
            "annotation-styles.gb",
            include_bytes!("../../../../validation/fixtures/regression/annotation-styles.gb"),
        )
        .unwrap();
        for circular in [false, true] {
            doc.circular = circular;
            let svg = map_svg(
                &doc,
                MapOptions::default(),
                &[],
                1000,
                750,
                true,
                |s, size| s.chars().count() as f32 * size * 0.55,
            )
            .unwrap();
            assert!(svg.contains("Between bases"));
            assert!(svg.contains("#0000ff"));
            if let Ok(dir) = std::env::var("POPGENE_FIGURE_QA_DIR") {
                let dir = std::path::PathBuf::from(dir);
                std::fs::create_dir_all(&dir).unwrap();
                for format in [
                    crate::figure::FigureFormat::Svg,
                    crate::figure::FigureFormat::Pdf,
                    crate::figure::FigureFormat::Png,
                ] {
                    std::fs::write(
                        dir.join(format!(
                            "styles-{}.{}",
                            if circular { "circular" } else { "linear" },
                            format.extension()
                        )),
                        crate::figure::encode(&svg, format, 1.).unwrap(),
                    )
                    .unwrap();
                }
            }
        }
    }
}
