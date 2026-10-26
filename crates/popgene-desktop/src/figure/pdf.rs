//! Vector pagination at complete sequence-row boundaries.
use pdf_writer::{Content, Finish, Name, Pdf, Rect, Ref, Str, TextStr};
use std::collections::HashMap;

fn slices(svg: &str, height: f32, page_height: u32) -> Result<Vec<(f32, f32)>, String> {
    if !(400..=8000).contains(&page_height) {
        return Err("PDF page height must be 400–8,000 layout pixels".into());
    }
    let xml = svg2pdf::usvg::roxmltree::Document::parse(svg).map_err(|e| e.to_string())?;
    let mut boundaries = Vec::new();
    let mut previous = 0.;
    for row in xml
        .descendants()
        .filter(|node| node.attribute("data-popgene-row-end").is_some())
    {
        let a = row
            .attribute("data-popgene-row-start")
            .ok_or("Missing sequence row start")?
            .parse::<f32>()
            .map_err(|_| "Invalid sequence row start")?;
        let b = row
            .attribute("data-popgene-row-end")
            .unwrap()
            .parse::<f32>()
            .map_err(|_| "Invalid sequence row end")?;
        if !a.is_finite() || !b.is_finite() || a < previous || b <= a || b > height {
            return Err("Invalid sequence row bounds".into());
        }
        boundaries.push(b);
        previous = b;
    }
    if boundaries.is_empty() {
        return Err("This figure has no sequence-row pagination information".into());
    }
    if previous < height {
        // Keep the final row and its trailing notes/margin together. Otherwise
        // an exactly full page could leave a numbered, blank final page.
        *boundaries.last_mut().unwrap() = height;
    }
    // Reserve 32 layout pixels (24 points) for page numbering.
    let capacity = page_height as f32 - 32.;
    let mut pages = Vec::new();
    let mut start = 0.;
    while start < height {
        let end = boundaries
            .iter()
            .copied()
            .rfind(|&b| b > start && b <= start + capacity)
            .ok_or("The chosen page height cannot fit a complete sequence row; increase it")?;
        pages.push((start, end));
        start = end;
    }
    Ok(pages)
}

fn page_svg(
    svg: &str,
    width: f32,
    full_height: f32,
    start: f32,
    end: f32,
) -> Result<String, String> {
    let xml = svg2pdf::usvg::roxmltree::Document::parse(svg).map_err(|e| e.to_string())?;
    let rows = xml
        .descendants()
        .filter(|n| n.attribute("data-popgene-row-end").is_some())
        .collect::<Vec<_>>();
    let first = rows.first().ok_or("Missing sequence rows")?;
    let last = rows.last().unwrap();
    let body = first.parent().ok_or("Missing sequence body")?;
    if rows.iter().any(|n| n.parent() != Some(body)) || body.parent() != Some(xml.root_element()) {
        return Err("Unsupported sequence figure grouping".into());
    }
    let mut output = format!(
        "<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{}' viewBox='0 0 {width} {}'>",
        end - start,
        end - start
    );
    for node in xml
        .root_element()
        .children()
        .filter(|n| n.is_element() && *n != body)
    {
        output.push_str(&svg[node.range()]);
    }
    output.push_str(&format!("<g transform='translate(0,-{start})'>"));
    let original = &svg[body.range()];
    output.push_str(&original[..=original.find('>').ok_or("Missing sequence group")?]);
    for node in body.children() {
        let include = if let Some(a) = node.attribute("data-popgene-row-start") {
            let a = a.parse::<f32>().map_err(|_| "Invalid row start")?;
            let b = node
                .attribute("data-popgene-row-end")
                .unwrap()
                .parse::<f32>()
                .map_err(|_| "Invalid row end")?;
            a >= start && b <= end
        } else {
            start == 0. && node.range().end <= first.range().start
                || end == full_height && node.range().start >= last.range().end
        };
        if include {
            output.push_str(&svg[node.range()]);
        }
    }
    output.push_str("</g></g></svg>");
    Ok(output)
}

pub(super) fn encode(
    tree: &svg2pdf::usvg::Tree,
    svg: &str,
    page_height: u32,
    options: &svg2pdf::usvg::Options,
) -> Result<Vec<u8>, String> {
    let pages = slices(svg, tree.size().height(), page_height)?;
    let mut alloc = Ref::new(1);
    let catalogue = alloc.bump();
    let page_tree = alloc.bump();
    let font = alloc.bump();
    let info = alloc.bump();
    let page_ids = pages.iter().map(|_| alloc.bump()).collect::<Vec<_>>();
    let mut pdf = Pdf::new();
    pdf.catalog(catalogue).pages(page_tree);
    pdf.pages(page_tree)
        .kids(page_ids.iter().copied())
        .count(pages.len() as i32);
    pdf.type1_font(font).base_font(Name(b"Helvetica"));
    pdf.document_info(info)
        .producer(TextStr("PopGene vector figures"));
    let width = tree.size().width() * 0.75;
    for (i, (&(start, end), page_id)) in pages.iter().zip(page_ids).enumerate() {
        // Remove off-page rows before conversion. Merely clipping the whole
        // figure leaves invisible duplicate text in PDF search/copy extraction.
        let svg = page_svg(svg, tree.size().width(), tree.size().height(), start, end)?;
        let rendered_page =
            svg2pdf::usvg::Tree::from_str(&svg, options).map_err(|e| e.to_string())?;
        let (chunk, graphic) =
            svg2pdf::to_chunk(&rendered_page, Default::default()).map_err(|e| e.to_string())?;
        let mut ids = HashMap::new();
        let chunk = chunk.renumber(|old| *ids.entry(old).or_insert_with(|| alloc.bump()));
        let graphic = ids[&graphic];
        let height = (end - start) * 0.75 + 24.;
        let contents = alloc.bump();
        let mut page = pdf.page(page_id);
        page.media_box(Rect::new(0., 0., width, height))
            .parent(page_tree)
            .contents(contents);
        let mut resources = page.resources();
        resources.x_objects().pair(Name(b"Figure"), graphic);
        resources.fonts().pair(Name(b"PageFont"), font);
        resources.finish();
        page.finish();
        let mut content = Content::new();
        content
            .save_state()
            .rect(0., 24., width, height - 24.)
            .clip_nonzero()
            .end_path()
            .transform([width, 0., 0., height - 24., 0., 24.])
            .x_object(Name(b"Figure"))
            .restore_state();
        let label = format!("Page {} of {}", i + 1, pages.len());
        content
            .begin_text()
            .set_font(Name(b"PageFont"), 9.)
            .next_line(12., 9.)
            .show(Str(label.as_bytes()))
            .end_text();
        pdf.stream(contents, &content.finish());
        pdf.extend(&chunk);
    }
    Ok(pdf.finish())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn pages_break_only_after_complete_rows_and_keep_tail_notes() {
        let svg = "<svg xmlns='http://www.w3.org/2000/svg'><g data-popgene-row-start='80' data-popgene-row-end='250'/><g data-popgene-row-start='250' data-popgene-row-end='500'/><g data-popgene-row-start='500' data-popgene-row-end='800'/></svg>";
        assert_eq!(
            slices(svg, 840., 600).unwrap(),
            vec![(0., 500.), (500., 840.)]
        );
        assert!(slices(svg, 840., 399).is_err());
        assert_eq!(
            slices(svg, 840., 400).unwrap(),
            vec![(0., 250.), (250., 500.), (500., 840.)]
        );
        let edge = "<svg xmlns='http://www.w3.org/2000/svg'><g data-popgene-row-start='80' data-popgene-row-end='400'/><g data-popgene-row-start='400' data-popgene-row-end='807'/></svg>";
        assert_eq!(
            slices(edge, 815., 839).unwrap(),
            vec![(0., 400.), (400., 815.)]
        );
    }
    #[test]
    fn sequence_rows_produce_a_complete_multipage_vector_document() {
        let mut document =
            popgene_core::Document::new("Pagination validation", &"ACGT".repeat(250)).unwrap();
        document
            .annotate(
                "Long annotated region".into(),
                popgene_core::Segment { start: 3, end: 997 },
                false,
            )
            .unwrap();
        let options = crate::sequence_view::SequenceOptions {
            reading_frames: 63,
            ..Default::default()
        };
        let svg =
            crate::sequence_figure::svg(&document, &options, &[], &[], 700., 100, false, None)
                .unwrap();
        let pdf =
            crate::figure::encode_with_pages(&svg, crate::figure::FigureFormat::Pdf, 1., Some(900))
                .unwrap();
        assert!(pdf.starts_with(b"%PDF-"));
        if let Ok(dir) = std::env::var("POPGENE_FIGURE_QA_DIR") {
            let dir = std::path::PathBuf::from(dir);
            std::fs::create_dir_all(&dir).unwrap();
            std::fs::write(dir.join("paginated-sequence.svg"), svg).unwrap();
            std::fs::write(dir.join("paginated-sequence.pdf"), pdf).unwrap();
        }
    }
}
