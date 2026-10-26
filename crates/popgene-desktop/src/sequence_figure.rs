//! Bounded vector sequence figures using the interactive view's coordinate layout.
use crate::{
    feature_card,
    sequence_layout::{self, SequenceColorIndex, SequenceLayout},
    sequence_view::{SequenceOptions, SequenceRestrictionSite, annotation_label_color, base_color},
};
use popgene_core::{
    Document, FeatureSegmentStyle, Molecule, Segment,
    translation::{
        GeneticCode, OrfOptions, cds_flanking_codons, find_orfs, frame_codons_with_options,
    },
};
use std::fmt::Write;

const LEFT: f32 = 64.;
const RIGHT: f32 = 70.;
const MAX_LANES: usize = 4;

fn escape(value: &str) -> String {
    value.chars().fold(String::new(), |mut output, c| {
        match c {
            '&' => output.push_str("&amp;"),
            '<' => output.push_str("&lt;"),
            '>' => output.push_str("&gt;"),
            '"' => output.push_str("&quot;"),
            '\'' => output.push_str("&apos;"),
            '\t'
            | '\n'
            | '\r'
            | '\u{20}'..='\u{d7ff}'
            | '\u{e000}'..='\u{fffd}'
            | '\u{10000}'..='\u{10ffff}' => output.push(c),
            _ => output.push('\u{fffd}'),
        }
        output
    })
}
fn color(value: gpui_kit::Rgba) -> u32 {
    ((value.r * 255.).round() as u32) << 16
        | ((value.g * 255.).round() as u32) << 8
        | (value.b * 255.).round() as u32
}
fn text(output: &mut String, x: f32, y: f32, size: f32, fill: u32, value: &str) {
    let _ = write!(
        output,
        "<text x='{x:.2}' y='{y:.2}' font-size='{size:.1}' fill='#{fill:06x}'>{}</text>",
        escape(value)
    );
}
fn line(output: &mut String, x1: f32, y1: f32, x2: f32, y2: f32, fill: u32, width: f32) {
    let _ = write!(
        output,
        "<path d='M{x1:.2} {y1:.2}L{x2:.2} {y2:.2}' fill='none' stroke='#{fill:06x}' stroke-width='{width}'/>"
    );
}
fn clipped_text(output: &mut String, id: &str, x: f32, y: f32, width: f32, fill: u32, value: &str) {
    let _ = write!(
        output,
        "<clipPath id='{id}'><rect x='{x:.2}' y='{y:.2}' width='{width:.2}' height='20'/></clipPath><g clip-path='url(#{id})'><title>{}</title>",
        escape(value)
    );
    text(output, x + 3., y + 13., 10., fill, value);
    output.push_str("</g>");
}

/// Empty `regions` exports the whole sequence. Nonempty regions are rendered in
/// their supplied biological/selection order, without creating false joins.
/// `flank_feature_id` identifies the selected CDS independently of the crop.
/// Its hypothetical flanks are included only where they lie in the export regions.
/// Width and row limits are enforced before allocating glyphs. Omitted rows,
/// crowded annotations, and unsupported CDS translations remain visibly named.
#[allow(clippy::too_many_arguments)]
pub fn svg(
    document: &Document,
    options: &SequenceOptions,
    sites: &[SequenceRestrictionSite],
    regions: &[Segment],
    width: f32,
    max_rows: usize,
    dark: bool,
    flank_feature_id: Option<u64>,
) -> Result<String, String> {
    document.validate().map_err(|e| e.to_string())?;
    if !width.is_finite() || !(360. ..=4000.).contains(&width) {
        return Err("Sequence figure width must be 360–4000 pixels".into());
    }
    if !(1..=1000).contains(&max_rows) {
        return Err("Sequence figures support 1–1000 rows per export".into());
    }
    let n = document.sequence.len();
    if n == 0 {
        return Err("An empty sequence has no symbols to export".into());
    }
    if regions.len() > 10_000 || regions.iter().any(|r| r.start >= r.end || r.end > n) {
        return Err("Figure selection contains an invalid or excessive number of ranges".into());
    }
    let code = GeneticCode::from_id(options.genetic_code)
        .ok_or_else(|| format!("Unknown genetic code {}", options.genetic_code))?;
    let whole = [Segment { start: 0, end: n }];
    let regions = if regions.is_empty() {
        &whole[..]
    } else {
        regions
    };
    let cw = if options.compact { 9. } else { 11. };
    let per_row = sequence_layout::bases_per_row(width, cw, LEFT + RIGHT);
    let total_rows: usize = regions
        .iter()
        .map(|r| (r.end - r.start).div_ceil(per_row))
        .sum();
    let mut rows = Vec::with_capacity(total_rows.min(max_rows));
    'regions: for region in regions {
        for start in (region.start..region.end).step_by(per_row) {
            if rows.len() == max_rows {
                break 'regions;
            }
            rows.push(Segment {
                start,
                end: (start + per_row).min(region.end),
            });
        }
    }
    let ink = if dark { 0xdce5ed } else { 0x172b3b };
    let muted = if dark { 0x8799ab } else { 0x52697a };
    let accent = if dark { 0x48c9bb } else { 0x087568 };
    let nucleotide = document.molecule != Molecule::Protein;
    let paired = options.show_complement && nucleotide && !document.single_stranded;
    let show_enzymes =
        options.show_enzymes && document.molecule == Molecule::Dna && !document.single_stranded;
    let layout = SequenceLayout::new(document);
    let colors = SequenceColorIndex::new(document);
    let orfs = if options.show_orfs {
        find_orfs(
            document,
            OrfOptions {
                genetic_code: code,
                min_codons: options.orf_min_codons,
                alternative_starts: options.orf_alternative_starts,
                include_partial: options.orf_include_partial,
                ..OrfOptions::default()
            },
        )
    } else {
        Default::default()
    };
    let mut notes = sequence_layout::duplex_display_notice(document)
        .into_iter()
        .collect::<Vec<_>>();
    if total_rows > rows.len() {
        notes.push(format!(
            "TRUNCATED: {} of {total_rows} rows exported; {} rows omitted.",
            rows.len(),
            total_rows - rows.len()
        ));
    }
    if options.show_translations && layout.omitted_translations > 0 {
        notes.push(format!(
            "{} CDS translations omitted:",
            layout.omitted_translations
        ));
        notes.extend(layout.translation_notices.iter().cloned());
    }
    if nucleotide
        && (options.show_translations || options.reading_frames != 0 || options.show_flanking)
    {
        notes.push(format!(
            "Calculated translations · reading-frame code {} · * stop · X unresolved codon",
            code.id()
        ));
    }
    if options.show_orfs {
        notes.push(format!(
            "ORF predictions · ≥{} aa · {} · {} · longest start per stop interval{}",
            options.orf_min_codons,
            if options.orf_alternative_starts {
                "code-specific starts"
            } else {
                "ATG starts"
            },
            if options.orf_include_partial {
                "open ends included"
            } else {
                "complete only"
            },
            if orfs.truncated {
                " · result limit reached"
            } else {
                ""
            }
        ));
    }
    let flanks = if options.show_flanking {
        if let Some(feature) = document
            .features
            .iter()
            .find(|f| f.kind.eq_ignore_ascii_case("CDS") && Some(f.id) == flank_feature_id)
        {
            notes.push(format!(
                "Hypothetical flanks of {} · up to {} codons each side or nearest stop",
                feature.label,
                options.flanking_codons.min(300)
            ));
            cds_flanking_codons(document, feature, options.flanking_codons).unwrap_or_else(
                |error| {
                    notes.push(format!("Flanking translation omitted: {error}"));
                    vec![]
                },
            )
        } else {
            notes.push("Flanking translation requires selection of one CDS annotation".into());
            vec![]
        }
    } else {
        vec![]
    };
    let mut body = String::new();
    text(
        &mut body,
        20.,
        27.,
        17.,
        ink,
        &feature_card::truncate(&document.name, ((width - 40.) / 10.) as usize),
    );
    text(
        &mut body,
        20.,
        46.,
        11.,
        muted,
        &format!(
            "{} · {} {} · {} · original sequence coordinates",
            document.molecule.label(),
            n,
            document.molecule.unit(),
            if document.circular {
                "circular"
            } else {
                "linear"
            }
        ),
    );
    let mut y = 65.;
    for note in notes {
        // Wrap the warning itself: omitting/clipping the warning would conceal data loss.
        let limit = ((width - 40.) / 6.).floor() as usize;
        let chars = note.chars().take(4000).collect::<Vec<_>>();
        for chunk in chars.chunks(limit.max(1)) {
            text(
                &mut body,
                20.,
                y,
                10.,
                muted,
                &chunk.iter().collect::<String>(),
            );
            y += 14.;
        }
        if note.chars().count() > 4000 {
            text(&mut body, 20., y, 10., muted, "[notice text truncated]");
            y += 14.;
        }
    }
    y += 12.;
    let header_end = y - 8.;
    for (row_index, range) in rows.iter().enumerate() {
        let start = range.start;
        let end = range.end;
        let row_width = (end - start) as f32 * cw;
        let row_layout = layout.row(
            start,
            end,
            MAX_LANES,
            options.show_features,
            options.show_primers,
            options.show_translations,
        );
        let row_top = y - 8.;
        let tag_start = body.len();
        let _ = write!(
            body,
            "<g data-row='{row_index}' data-start='{}' data-end='{end}'>",
            start + 1
        );
        let tag_end = body.len();
        line(&mut body, 20., y - 7., width - 20., y - 7., muted, 0.4);
        let mut omitted_enzymes = 0;
        let mut labels = vec![];
        let mut ends = [0_f32; 2];
        if show_enzymes {
            let mut row_sites = sites
                .iter()
                .filter(|site| {
                    let at = site.top_cut.unwrap_or(site.start).min(n - 1);
                    at >= start && at < end
                })
                .collect::<Vec<_>>();
            row_sites.sort_by_key(|site| site.top_cut.unwrap_or(site.start));
            for site in row_sites {
                let label_width = (site.label.chars().count() as f32 * 7. + 12.)
                    .min(155.)
                    .min(row_width);
                let x = ((site.top_cut.unwrap_or(site.start).saturating_sub(start)) as f32 * cw)
                    .min((row_width - label_width).max(0.));
                if let Some(lane) = ends.iter().position(|end| *end <= x) {
                    ends[lane] = x + label_width + 8.;
                    labels.push((site, x, label_width, lane));
                } else {
                    omitted_enzymes += 1;
                }
            }
        }
        for (index, (site, x, w, lane)) in labels.iter().enumerate() {
            clipped_text(
                &mut body,
                &format!("enzyme-{row_index}-{index}"),
                LEFT + x,
                y + *lane as f32 * 21.,
                *w,
                accent,
                &site.label,
            );
        }
        y += labels
            .iter()
            .map(|(_, _, _, lane)| lane + 1)
            .max()
            .unwrap_or(0) as f32
            * 21.;
        text(
            &mut body,
            20.,
            y + 14.,
            10.,
            muted,
            &(start + 1).to_string(),
        );
        for position in start..end {
            if (position + 1).is_multiple_of(10) {
                text(
                    &mut body,
                    LEFT + (position - start) as f32 * cw - 4.,
                    y + 10.,
                    9.,
                    muted,
                    &(position + 1).to_string(),
                );
            }
        }
        let top = y + 18.;
        for (reverse, row_y) in [(false, top), (true, top + 30.)] {
            if reverse && !paired {
                continue;
            }
            text(
                &mut body,
                20.,
                row_y + 15.,
                10.,
                muted,
                if !nucleotide {
                    "N"
                } else if reverse {
                    "3′"
                } else {
                    "5′"
                },
            );
            text(
                &mut body,
                LEFT + row_width + 9.,
                row_y + 15.,
                10.,
                muted,
                if !nucleotide {
                    "C"
                } else if reverse {
                    "5′"
                } else {
                    "3′"
                },
            );
            for position in start..end {
                let raw = document.sequence.as_bytes()[position];
                let symbol = if reverse {
                    let Some(opposite) = sequence_layout::opposite_base(document, position) else {
                        continue;
                    };
                    opposite
                } else {
                    raw
                };
                let fill = colors.color(position, reverse).unwrap_or_else(|| {
                    if options.colored_bases {
                        color(base_color(symbol))
                    } else {
                        ink
                    }
                });
                let _ = write!(
                    body,
                    "<g data-base='{}' data-strand='{}'>",
                    position + 1,
                    if reverse { "reverse" } else { "forward" }
                );
                text(
                    &mut body,
                    LEFT + (position - start) as f32 * cw + 1.,
                    row_y + 16.,
                    if options.compact { 13. } else { 15. },
                    fill,
                    &(symbol as char).to_string(),
                );
                body.push_str("</g>");
            }
        }
        if show_enzymes {
            let opposite_span = sequence_layout::opposite_span(document);
            for site in sites {
                for (boundary, reverse, pair) in [
                    (site.top_cut, false, 1),
                    (site.bottom_cut, true, 1),
                    (site.second_top_cut, false, 2),
                    (site.second_bottom_cut, true, 2),
                ] {
                    let Some(boundary) = boundary else {
                        continue;
                    };
                    if boundary < start
                        || boundary > end
                        || (boundary == end && end != n)
                        || (reverse
                            && (!paired
                                || boundary < opposite_span.start
                                || boundary > opposite_span.end))
                    {
                        continue;
                    }
                    let x = LEFT + (boundary - start) as f32 * cw;
                    let cut_y = top + if reverse { 30. } else { 0. };
                    let _ = write!(
                        body,
                        "<g data-cut='{boundary}' data-pair='{pair}' data-strand='{}'><title>{}</title>",
                        if reverse { "reverse" } else { "forward" },
                        escape(&site.label)
                    );
                    line(&mut body, x, cut_y, x, cut_y + 21., accent, 1.5);
                    line(
                        &mut body,
                        x,
                        cut_y + 21.,
                        x + if reverse { -4. } else { 4. },
                        cut_y + 21.,
                        accent,
                        1.5,
                    );
                    body.push_str("</g>");
                }
            }
        }
        y = top + if paired { 61. } else { 31. };
        if let Some(qualities) = &document.sequence_qualities {
            text(&mut body, 20., y + 10., 9., muted, "Q");
            for (position, quality) in qualities.iter().enumerate().take(end).skip(start) {
                text(
                    &mut body,
                    LEFT + (position - start) as f32 * cw,
                    y + 10.,
                    8.,
                    if quality.is_some_and(|q| q < 20) {
                        0xb26a0e
                    } else {
                        muted
                    },
                    &quality.map_or("—".into(), |q| q.to_string()),
                );
            }
            y += 20.;
        }
        if paired && document.opposite_sequence_qualities.is_some() {
            text(&mut body, 20., y + 10., 9., muted, "Q−");
            for position in start..end {
                if sequence_layout::opposite_base(document, position).is_some() {
                    let quality = sequence_layout::opposite_quality(document, position);
                    text(
                        &mut body,
                        LEFT + (position - start) as f32 * cw,
                        y + 10.,
                        8.,
                        if quality.is_some_and(|q| q < 20) {
                            0xb26a0e
                        } else {
                            muted
                        },
                        &quality.map_or("—".into(), |q| q.to_string()),
                    );
                }
            }
            y += 20.;
        }
        if nucleotide {
            for (bit, frame) in [1, 2, 3, -1, -2, -3].into_iter().enumerate() {
                if options.reading_frames & (1 << bit) == 0 {
                    continue;
                }
                text(&mut body, 20., y + 14., 10., muted, &format!("{frame:+}"));
                let _ = write!(body, "<g data-frame='{frame}'>");
                for codon in frame_codons_with_options(
                    &document.sequence,
                    frame,
                    start..end,
                    code,
                    document.circular,
                ) {
                    text(
                        &mut body,
                        LEFT + (codon.center - start) as f32 * cw,
                        y + 14.,
                        12.,
                        ink,
                        &(codon.amino_acid as char).to_string(),
                    );
                }
                body.push_str("</g>");
                y += 24.;
            }
        }
        for frame in [1, 2, 3, -1, -2, -3] {
            let matching = orfs
                .orfs
                .iter()
                .enumerate()
                .filter(|(_, orf)| orf.frame(n) == frame)
                .flat_map(|(i, orf)| {
                    orf.segments(n).into_iter().filter_map(move |part| {
                        let left = part.start.max(start);
                        let right = part.end.min(end);
                        (left < right).then_some((i, orf, left, right))
                    })
                })
                .collect::<Vec<_>>();
            if matching.is_empty() {
                continue;
            }
            text(&mut body, 6., y + 13., 9., muted, &format!("ORF{frame:+}"));
            for (i, orf, left, right) in matching {
                let x = LEFT + (left - start) as f32 * cw;
                let w = (right - left) as f32 * cw;
                let fill = if orf.reverse { 0x6f9ed8 } else { 0xddaa50 };
                let _ = write!(
                    body,
                    "<rect data-orf='{i}' x='{x}' y='{y}' width='{w}' height='18' fill='#{fill:06x}' fill-opacity='.4'/>"
                );
                clipped_text(
                    &mut body,
                    &format!("orf-{row_index}-{frame}-{i}-{left}"),
                    x,
                    y,
                    w,
                    ink,
                    &format!(
                        "ORF {} · {} aa{}",
                        i + 1,
                        orf.amino_acid_length,
                        if orf.complete { "" } else { " · open end" }
                    ),
                );
            }
            y += 22.;
        }
        for upstream in [true, false] {
            let matching = flanks
                .iter()
                .filter(|f| {
                    f.upstream == upstream && f.codon.center >= start && f.codon.center < end
                })
                .collect::<Vec<_>>();
            if matching.is_empty() {
                continue;
            }
            text(
                &mut body,
                4.,
                y + 13.,
                9.,
                muted,
                if upstream { "Flank5′" } else { "Flank3′" },
            );
            for flank in matching {
                text(
                    &mut body,
                    LEFT + (flank.codon.center - start) as f32 * cw,
                    y + 13.,
                    11.,
                    ink,
                    &(flank.codon.amino_acid as char).to_string(),
                );
            }
            y += 22.;
        }
        for lane in 0..row_layout.lanes {
            for (span_index, span) in row_layout
                .spans
                .iter()
                .enumerate()
                .filter(|(_, span)| span.lane == lane)
            {
                let feature = &document.features[span.feature_index];
                let fill = span.color_rgb.unwrap_or_else(|| {
                    color(feature_card::feature_color(feature, span.feature_index).to_rgb())
                });
                let x = LEFT + span.start as f32 * cw;
                let w = ((span.end - span.start) as f32 * cw).max(3.);
                let direction = feature_card::direction(feature, document.molecule);
                let directional = direction != "Nondirectional" && nucleotide;
                let left_arrow = directional
                    && (direction == "Reverse strand" || direction == "Bidirectional")
                    && !span.continued_left;
                let right_arrow =
                    directional && (direction != "Reverse strand") && !span.continued_right;
                let tip = 7_f32.min(w * 0.5);
                let left_inset = if left_arrow { tip } else { 0. };
                let right_inset = if right_arrow { tip } else { 0. };
                let _ = write!(
                    body,
                    "<g data-feature='{}' data-style='{:?}'><title>{}</title>",
                    feature.id,
                    span.style,
                    escape(&feature.label)
                );
                if span.style == FeatureSegmentStyle::Filled {
                    let _ = write!(
                        body,
                        "<path d='M{} {y}L{} {y}L{} {}L{} {}L{} {}L{x} {}Z' fill='#{fill:06x}'/>",
                        x + left_inset,
                        x + w - right_inset,
                        x + w,
                        y + 9.5,
                        x + w - right_inset,
                        y + 19.,
                        x + left_inset,
                        y + 19.,
                        y + 9.5
                    );
                } else {
                    line(&mut body, x, y + 17., x + w, y + 17., fill, 1.5);
                    if left_arrow {
                        line(&mut body, x + tip, y + 13., x, y + 17., fill, 1.5);
                    }
                    if right_arrow {
                        line(&mut body, x + w - tip, y + 13., x + w, y + 17., fill, 1.5);
                    }
                }
                let label_color = if span.style == FeatureSegmentStyle::Filled {
                    color(annotation_label_color(gpui_kit::rgb(fill).into()))
                } else {
                    ink
                };
                for boundary in &span.cleavage_sites {
                    let x = LEFT + *boundary as f32 * cw;
                    let _ = write!(body, "<g data-cleavage='{}'>", start + boundary);
                    line(&mut body, x, y - 2., x, y + 21., label_color, 2.);
                    body.push_str("</g>");
                }
                if w >= 32.
                    && (span.style != FeatureSegmentStyle::Gap || span.segment_name.is_some())
                {
                    clipped_text(
                        &mut body,
                        &format!("feature-{row_index}-{span_index}"),
                        x,
                        y,
                        w,
                        label_color,
                        &feature_card::truncate(
                            span.segment_name.as_deref().unwrap_or(&feature.label),
                            100,
                        ),
                    );
                }
                body.push_str("</g>");
            }
            for glyph in row_layout
                .translations
                .iter()
                .filter(|glyph| glyph.lane == lane)
            {
                let _ = write!(
                    body,
                    "<g data-cds='{}' data-code='{}'><title>{}</title>",
                    glyph.feature_id,
                    glyph.genetic_code,
                    escape(&format!(
                        "Calculated CDS · bases {}{}",
                        glyph
                            .codon_positions
                            .iter()
                            .take(glyph.nucleotide_count)
                            .map(|p| (p + 1).to_string())
                            .collect::<Vec<_>>()
                            .join(", "),
                        if glyph.exception {
                            " · translation exception"
                        } else {
                            ""
                        }
                    ))
                );
                text(
                    &mut body,
                    LEFT + glyph.base as f32 * cw,
                    y + 34.,
                    12.,
                    ink,
                    &glyph.amino_acid.to_string(),
                );
                body.push_str("</g>");
            }
            y += if row_layout.translation_lanes.contains(&lane) {
                43.
            } else {
                24.
            };
        }
        if row_layout.omitted > 0 || omitted_enzymes > 0 {
            text(
                &mut body,
                LEFT,
                y + 12.,
                10.,
                muted,
                &format!(
                    "{} annotation parts · {} enzyme labels omitted from crowded lanes",
                    row_layout.omitted, omitted_enzymes
                ),
            );
            y += 20.;
        }
        y += 22.;
        body.replace_range(tag_start..tag_end, &format!("<g data-row='{row_index}' data-start='{}' data-end='{end}' data-popgene-row-start='{row_top:.2}' data-popgene-row-end='{:.2}'>", start + 1, y - 8.));
        body.push_str("</g>");
        if body.len() > 30 * 1024 * 1024 {
            return Err(
                "Sequence figure exceeds the 30 MiB content limit; select fewer rows".into(),
            );
        }
    }
    Ok(format!(
        "<svg xmlns='http://www.w3.org/2000/svg' data-popgene-header-end='{header_end:.2}' width='{width:.0}' height='{y:.0}' viewBox='0 0 {width:.0} {y:.0}'><title>{}</title><rect width='100%' height='100%' fill='#{:06x}'/><g font-family='monospace'>{body}</g></svg>",
        escape(&feature_card::truncate(&document.name, 1000)),
        if dark { 0x17212d } else { 0xffffff }
    ))
}

#[cfg(test)]
mod tests {
    use super::svg;
    use crate::sequence_view::{SequenceOptions, SequenceRestrictionSite};
    use popgene_core::{Document, Feature, Segment, metadata::SequenceColor};

    #[test]
    fn vector_layers_use_original_coordinates_and_escape_user_text() {
        let mut doc = Document::new("map & <label>", "ATGAAATAACCCGGGTTTAAAATGA").unwrap();
        doc.features.push(Feature {
            id: 1,
            kind: "CDS".into(),
            label: "gene <a> & \"b\"".into(),
            parts: vec![Segment { start: 0, end: 9 }],
            reverse: false,
            qualifiers: vec![],
            edited: false,
            location: None,
        });
        doc.sequence_colors.push(SequenceColor {
            range: Segment { start: 1, end: 3 },
            reverse: true,
            color_rgb: 0x123456,
        });
        doc.sequence_qualities = Some(vec![Some(30); doc.sequence.len()]);
        let sites = [SequenceRestrictionSite {
            enzyme_index: 0,
            start: 2,
            top_cut: Some(2),
            bottom_cut: Some(4),
            second_top_cut: Some(8),
            second_bottom_cut: Some(10),
            label: "enzyme <2>".into(),
        }];
        let options = SequenceOptions {
            reading_frames: 63,
            show_orfs: true,
            orf_min_codons: 2,
            ..SequenceOptions::default()
        };
        let image = svg(&doc, &options, &sites, &[], 900., 100, false, None).unwrap();
        assert!(image.contains("map &amp; &lt;label&gt;"));
        assert!(image.contains("gene &lt;a&gt; &amp; &quot;b&quot;"));
        assert!(!image.contains("<a>"));
        assert!(image.contains("fill='#123456'"));
        assert_eq!(image.matches("data-cut='").count(), 4);
        assert_eq!(image.matches("data-frame='").count(), 6);
        assert_eq!(image.matches("data-cds='1'").count(), 3);
        assert!(image.contains("data-orf='"));
        svg2pdf::usvg::Tree::from_str(&image, &Default::default()).unwrap();
        if let Ok(directory) = std::env::var("POPGENE_SEQUENCE_FIGURE_QA_DIR") {
            let directory = std::path::Path::new(&directory);
            std::fs::create_dir_all(directory).unwrap();
            std::fs::write(directory.join("sequence.svg"), &image).unwrap();
            std::fs::write(
                directory.join("sequence.png"),
                crate::figure::encode(&image, crate::figure::FigureFormat::Png, 1.).unwrap(),
            )
            .unwrap();
        }
    }
    #[test]
    fn selection_order_and_row_limits_are_explicit() {
        let mut doc = Document::new("circular", &"ACGT".repeat(100)).unwrap();
        doc.circular = true;
        let selected = [
            Segment {
                start: 390,
                end: 400,
            },
            Segment { start: 0, end: 100 },
        ];
        let full = svg(
            &doc,
            &SequenceOptions::default(),
            &[],
            &selected,
            360.,
            100,
            true,
            None,
        )
        .unwrap();
        let first = full.find("data-start='391'").unwrap();
        let next = full.find("data-start='1'").unwrap();
        assert!(first < next);
        let read = |tag: &str, name: &str| -> f32 {
            tag.split(&format!("{name}='"))
                .nth(1)
                .unwrap()
                .split("'")
                .next()
                .unwrap()
                .parse()
                .unwrap()
        };
        let rows = full
            .split("<g data-row=")
            .skip(1)
            .map(|tag| {
                (
                    read(tag, "data-popgene-row-start"),
                    read(tag, "data-popgene-row-end"),
                )
            })
            .collect::<Vec<_>>();
        assert!(!rows.is_empty());
        assert_eq!(rows[0].0, read(&full, "data-popgene-header-end"));
        assert!(rows.iter().all(|(start, end)| start < end));
        assert!(rows.windows(2).all(|pair| pair[0].1 == pair[1].0));

        let bounded = svg(
            &doc,
            &SequenceOptions::default(),
            &[],
            &selected,
            360.,
            1,
            true,
            None,
        )
        .unwrap();
        assert_eq!(bounded.matches("data-row='").count(), 1);
        assert!(bounded.contains("TRUNCATED:"));
        assert!(bounded.contains("5 rows omitted"));
        assert!(
            svg(
                &doc,
                &SequenceOptions::default(),
                &[],
                &[],
                f32::NAN,
                1,
                true,
                None
            )
            .is_err()
        );
        assert!(
            svg(
                &doc,
                &SequenceOptions::default(),
                &[],
                &[],
                900.,
                1001,
                true,
                None
            )
            .is_err()
        );
        assert!(
            svg(
                &doc,
                &SequenceOptions::default(),
                &[],
                &[Segment { start: 0, end: 401 }],
                900.,
                1,
                true,
                None
            )
            .is_err()
        );
    }
    #[test]
    fn physical_duplex_figure_omits_nonexistent_opposite_residues() {
        use popgene_core::{
            restriction::{EndKind, FragmentEnd},
            workflow::DuplexEnds,
        };
        let mut doc = Document::new("physical", "AAACGTCC").unwrap();
        doc.duplex_ends = Some(DuplexEnds {
            left: FragmentEnd {
                kind: EndKind::FivePrime,
                overhang: "AA".into(),
                enzymes: vec![],
            },
            right: FragmentEnd {
                kind: EndKind::ThreePrime,
                overhang: "CC".into(),
                enzymes: vec![],
            },
            opposite_sequence: "ACGT".into(),
            left_modification: None,
            right_modification: None,
        });
        let image = svg(
            &doc,
            &SequenceOptions::default(),
            &[],
            &[],
            900.,
            100,
            false,
            None,
        )
        .unwrap();
        assert_eq!(image.matches("data-strand='forward'").count(), 8);
        assert_eq!(image.matches("data-strand='reverse'").count(), 4);
        assert!(!image.contains("data-base='1' data-strand='reverse'"));
        assert!(!image.contains("data-base='8' data-strand='reverse'"));
        assert!(image.contains("data-base='3' data-strand='reverse'"));
        assert!(image.contains("Top-strand coordinates shown"));
    }
}
