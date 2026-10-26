//! Embedded restriction-enzyme inspector with searchable catalogue navigation.
use crate::feature_card::{number, truncate};
use gpui_kit::component::{
    button::*,
    input::{Input, InputEvent, InputState},
    *,
};
use gpui_kit::{prelude::*, *};
use popgene_core::{
    Document,
    restriction::{ENZYMES, Enzyme, RestrictionAnalysis, RestrictionSite},
};

const BG: u32 = 0x101722;
const PANEL: u32 = 0x17212e;
const BORDER: u32 = 0x2a3646;
const INK: u32 = 0xe4edf4;
const MUTED: u32 = 0x8999ad;
const ACCENT: u32 = 0x48c9bb;
const PAGE_SIZE: usize = 24;

actions!(enzyme_details, [CloseDetails]);

pub enum EnzymeSidebarEvent {
    Close,
    Inspect(usize),
    /// Zero-based recognition start in the current document.
    SelectSite(usize),
}

#[derive(Clone, Debug)]
struct SiteDetails {
    start: usize,
    recognition: String,
    cuts: String,
    strand: &'static str,
    bases: String,
}

/// A bounded snapshot: no large sequence allocation or analysis occurs on opening.
#[derive(Clone, Debug)]
pub struct Details {
    enzyme: &'static Enzyme,
    enzyme_index: usize,
    counts: Vec<usize>,
    document_name: String,
    document_length: usize,
    circular: bool,
    count: usize,
    ambiguous: bool,
    selected: Option<usize>,
    sites: Vec<SiteDetails>,
}

impl Details {
    pub fn from_document(
        document: &Document,
        analysis: &RestrictionAnalysis,
        enzyme_index: usize,
        selected_start: Option<usize>,
        ambiguous: bool,
    ) -> Option<Self> {
        let enzyme = ENZYMES.get(enzyme_index)?;
        let sites = analysis
            .sites
            .iter()
            .filter(|site| site.enzyme_index == enzyme_index)
            .map(|site| site_details(document, enzyme, site))
            .collect::<Vec<_>>();
        let selected = selected_start.and_then(|start| sites.iter().position(|s| s.start == start));
        Some(Self {
            enzyme,
            enzyme_index,
            counts: analysis.counts.clone(),
            document_name: truncate(&document.name, 180),
            document_length: document.sequence.len(),
            circular: document.circular,
            count: *analysis.counts.get(enzyme_index)?,
            ambiguous,
            selected,
            sites,
        })
    }
}

fn cut_label(cut: Option<usize>, length: usize, circular: bool) -> String {
    match cut {
        Some(0) if circular => format!("after {}", number(length)),
        Some(0) => "before base 1".into(),
        Some(value) => format!("after {}", number(value)),
        None => "outside molecule".into(),
    }
}

fn site_details(document: &Document, enzyme: &Enzyme, site: &RestrictionSite) -> SiteDetails {
    let n = document.sequence.len();
    let length = site.recognition_length;
    let end = site.start.saturating_add(length);
    let recognition = if end > n && n > 0 {
        format!(
            "{}–{} / 1–{} · crosses origin",
            number(site.start + 1),
            number(n),
            number(end - n)
        )
    } else {
        format!("{}–{}", number(site.start + 1), number(end))
    };
    let bases = if n > 0 {
        (0..length)
            .map(|i| document.sequence.as_bytes()[(site.start + i) % n] as char)
            .collect()
    } else {
        String::new()
    };
    SiteDetails {
        start: site.start,
        recognition,
        cuts: if !enzyme.cleavage_known {
            "Cleavage positions unknown".into()
        } else {
            let mut cuts = format!(
                "Forward {}; opposite {}",
                cut_label(site.top_cut, n, document.circular),
                cut_label(site.bottom_cut, n, document.circular)
            );
            if enzyme.second_cut.is_some() {
                cuts.push_str(&format!(
                    " · second pair: forward {}; opposite {}",
                    cut_label(site.second_top_cut, n, document.circular),
                    cut_label(site.second_bottom_cut, n, document.circular)
                ));
            }
            cuts
        },
        strand: if site.reverse {
            "Reverse match"
        } else {
            "Forward match"
        },
        bases,
    }
}

fn complement(base: u8) -> char {
    match base {
        b'A' => 'T',
        b'T' => 'A',
        b'G' => 'C',
        b'C' => 'G',
        b'R' => 'Y',
        b'Y' => 'R',
        b'K' => 'M',
        b'M' => 'K',
        b'B' => 'V',
        b'V' => 'B',
        b'D' => 'H',
        b'H' => 'D',
        b'S' => 'S',
        b'W' => 'W',
        _ => 'N',
    }
}

fn palindrome(enzyme: &Enzyme) -> bool {
    enzyme
        .recognition
        .split('|')
        .all(|motif| motif.bytes().rev().map(complement).eq(motif.chars()))
}

fn strand_diagram(sequence: &str, cut: i32, arrow: char) -> String {
    if cut < 0 {
        format!("{arrow} ({} bases before) … {sequence}", -cut)
    } else if cut as usize > sequence.len() {
        format!("{sequence} … {arrow} (+{})", cut as usize - sequence.len())
    } else {
        let (left, right) = sequence.split_at(cut as usize);
        format!("{left}{arrow}{right}")
    }
}

fn overhang(enzyme: &Enzyme) -> String {
    if !enzyme.cleavage_known {
        return "Cleavage positions unknown".into();
    }
    let mut label = end_label(enzyme.top_cut, enzyme.bottom_cut);
    if let Some((top, bottom)) = enzyme.second_cut {
        label.push_str(&format!(" · second pair: {}", end_label(top, bottom)));
    }
    label
}
fn end_label(top: i32, bottom: i32) -> String {
    let delta = bottom - top;
    match delta.cmp(&0) {
        std::cmp::Ordering::Equal => "Blunt ends".into(),
        std::cmp::Ordering::Greater => format!("{}-nt 5′ overhang", delta),
        std::cmp::Ordering::Less => format!("{}-nt 3′ overhang", -delta),
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
enum EndFilter {
    #[default]
    All,
    FivePrime,
    ThreePrime,
    Blunt,
}

fn catalogue_indices(query: &str, ends: EndFilter) -> Vec<usize> {
    let query = query.trim().to_ascii_lowercase();
    let mut indices = ENZYMES
        .iter()
        .enumerate()
        .filter(|(_, enzyme)| {
            let delta = enzyme.bottom_cut - enzyme.top_cut;
            let end_matches = match ends {
                EndFilter::All => true,
                EndFilter::FivePrime => enzyme.cleavage_known && delta > 0,
                EndFilter::ThreePrime => enzyme.cleavage_known && delta < 0,
                EndFilter::Blunt => enzyme.cleavage_known && delta == 0,
            };
            end_matches
                && (enzyme.name.to_ascii_lowercase().contains(&query)
                    || enzyme.recognition.to_ascii_lowercase().contains(&query))
        })
        .map(|(index, _)| index)
        .collect::<Vec<_>>();
    indices.sort_by_key(|index| ENZYMES[*index].name.to_ascii_lowercase());
    indices
}

fn adjacent(indices: &[usize], current: usize, forward: bool) -> Option<usize> {
    if indices.is_empty() {
        return None;
    }
    match indices.iter().position(|index| *index == current) {
        Some(position) if forward => indices.get(position + 1).copied(),
        Some(position) => position
            .checked_sub(1)
            .and_then(|p| indices.get(p).copied()),
        None if forward => indices.first().copied(),
        None => indices.last().copied(),
    }
}

pub struct EnzymeSidebar {
    details: Details,
    focus: FocusHandle,
    page: usize,
    scroll: ScrollHandle,
    catalogue_scroll: ScrollHandle,
    search: Entity<InputState>,
    ends: EndFilter,
}

impl EventEmitter<EnzymeSidebarEvent> for EnzymeSidebar {}

impl EnzymeSidebar {
    pub fn new(details: Details, window: &mut Window, cx: &mut Context<Self>) -> Self {
        let focus = cx.focus_handle();
        focus.focus(window, cx);
        let search = cx
            .new(|cx| InputState::new(window, cx).placeholder("Find enzyme or recognition motif…"));
        cx.subscribe(&search, |this, _, event, cx| {
            if matches!(event, InputEvent::Change) {
                this.catalogue_scroll.set_offset(point(px(0.), px(0.)));
                cx.notify();
            }
        })
        .detach();
        let page = details.selected.unwrap_or(0) / PAGE_SIZE;
        Self {
            details,
            focus,
            page,
            scroll: ScrollHandle::new(),
            catalogue_scroll: ScrollHandle::new(),
            search,
            ends: EndFilter::All,
        }
    }

    pub fn replace(&mut self, details: Details, cx: &mut Context<Self>) {
        let changed_selection = self.details.enzyme_index != details.enzyme_index
            || self.details.selected != details.selected
            || self.details.document_name != details.document_name;
        if changed_selection {
            self.page = details.selected.unwrap_or(0) / PAGE_SIZE;
            self.scroll.set_offset(point(px(0.), px(0.)));
        } else {
            self.page = self
                .page
                .min(details.sites.len().saturating_sub(1) / PAGE_SIZE);
        }
        self.details = details;
        cx.notify();
    }
}

/// Only Escape is local to the inspector. Document and app shortcuts bubble up.
pub fn init(cx: &mut App) {
    cx.bind_keys([KeyBinding::new(
        "escape",
        CloseDetails,
        Some("EnzymeDetails"),
    )]);
}

fn label(text: impl Into<SharedString>) -> Div {
    div().text_xs().text_color(rgb(MUTED)).child(text.into())
}

impl Render for EnzymeSidebar {
    fn render(&mut self, _: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        let d = &self.details;
        let e = d.enzyme;
        let diagrams = e
            .recognition
            .split('|')
            .flat_map(|motif| {
                let opposite = motif.bytes().map(complement).collect::<String>();
                if !e.cleavage_known {
                    return vec![format!("5′ {motif} 3′"), format!("3′ {opposite} 5′")];
                }
                let mut lines = vec![
                    format!("5′ {} 3′", strand_diagram(motif, e.top_cut, '↓')),
                    format!("3′ {} 5′", strand_diagram(&opposite, e.bottom_cut, '↑')),
                ];
                if let Some((top, bottom)) = e.second_cut {
                    lines.extend([
                        format!("Second pair: 5′ {} 3′", strand_diagram(motif, top, '↓')),
                        format!(
                            "Second pair: 3′ {} 5′",
                            strand_diagram(&opposite, bottom, '↑')
                        ),
                    ]);
                }
                lines
            })
            .map(|text| {
                div()
                    .font_family("monospace")
                    .text_size(px(14.))
                    .child(text)
            })
            .collect::<Vec<_>>();
        let indices = catalogue_indices(&self.search.read(cx).value(), self.ends);
        let previous = adjacent(&indices, d.enzyme_index, false);
        let next = adjacent(&indices, d.enzyme_index, true);
        let catalogue_rows = indices
            .iter()
            .map(|index| {
                let index = *index;
                let enzyme = &ENZYMES[index];
                let count = d.counts.get(index).copied().unwrap_or(0);
                Button::new(("catalogue-enzyme", index))
                    .small()
                    .ghost()
                    .w_full()
                    .selected(index == d.enzyme_index)
                    .flex_shrink_0()
                    .child(
                        div()
                            .flex()
                            .w_full()
                            .items_center()
                            .gap_2()
                            .child(div().flex_1().text_left().child(enzyme.name))
                            .child(label(format!(
                                "{} {}",
                                number(count),
                                if count == 1 { "site" } else { "sites" }
                            ))),
                    )
                    .accessibility_label(format!(
                        "Inspect {} restriction enzyme, {} sites",
                        enzyme.name, count
                    ))
                    .on_click(
                        cx.listener(move |_, _, _, cx| cx.emit(EnzymeSidebarEvent::Inspect(index))),
                    )
            })
            .collect::<Vec<_>>();
        let end_filters = [
            ("all", "All", EndFilter::All),
            ("five", "5′ overhang", EndFilter::FivePrime),
            ("three", "3′ overhang", EndFilter::ThreePrime),
            ("blunt", "Blunt", EndFilter::Blunt),
        ]
        .into_iter()
        .map(|(id, name, ends)| {
            Button::new(format!("enzyme-end-filter-{id}"))
                .xsmall()
                .ghost()
                .label(name)
                .accessibility_label(format!("Filter catalogue: {name}"))
                .selected(self.ends == ends)
                .on_click(cx.listener(move |this, _, _, cx| {
                    this.ends = ends;
                    this.catalogue_scroll.set_offset(point(px(0.), px(0.)));
                    cx.notify();
                }))
        })
        .collect::<Vec<_>>();
        let catalogue = div()
            .flex()
            .flex_col()
            .gap_2()
            .p_3()
            .border_b_1()
            .border_color(rgb(BORDER))
            .child(Input::new(&self.search).aria_label("Search enzyme catalogue"))
            .child(div().flex().flex_wrap().gap_1().children(end_filters))
            .child(
                div()
                    .flex()
                    .items_center()
                    .gap_2()
                    .child(label(format!(
                        "{} of {} enzymes",
                        indices.len(),
                        ENZYMES.len()
                    )))
                    .child(div().flex_1())
                    .child(
                        Button::new("previous-enzyme")
                            .xsmall()
                            .ghost()
                            .label("Previous")
                            .accessibility_label("Inspect previous enzyme")
                            .disabled(previous.is_none())
                            .on_click(cx.listener(move |_, _, _, cx| {
                                if let Some(index) = previous {
                                    cx.emit(EnzymeSidebarEvent::Inspect(index));
                                }
                            })),
                    )
                    .child(
                        Button::new("next-enzyme")
                            .xsmall()
                            .ghost()
                            .label("Next")
                            .accessibility_label("Inspect next enzyme")
                            .disabled(next.is_none())
                            .on_click(cx.listener(move |_, _, _, cx| {
                                if let Some(index) = next {
                                    cx.emit(EnzymeSidebarEvent::Inspect(index));
                                }
                            })),
                    ),
            )
            .child(
                div()
                    .id("enzyme-catalogue-scroll")
                    .h(px(112.))
                    .overflow_y_scroll()
                    .track_scroll(&self.catalogue_scroll)
                    .flex()
                    .flex_col()
                    .gap_1()
                    .children(catalogue_rows)
                    .when(indices.is_empty(), |el| {
                        el.child(label("No enzymes match these filters."))
                    }),
            );
        let page_start = self.page * PAGE_SIZE;
        let page_end = (page_start + PAGE_SIZE).min(d.sites.len());
        let mut site_rows = Vec::new();
        for (index, site) in d.sites.iter().enumerate().skip(page_start).take(PAGE_SIZE) {
            let start = site.start;
            site_rows.push(
                div()
                    .p_3()
                    .rounded_md()
                    .bg(rgb(if d.selected == Some(index) {
                        0x203d42
                    } else {
                        PANEL
                    }))
                    .flex()
                    .flex_col()
                    .gap_1()
                    .child(
                        Button::new(("enzyme-site", index))
                            .small()
                            .ghost()
                            .label(format!(
                                "{}{}",
                                if d.selected == Some(index) {
                                    "Selected · "
                                } else {
                                    ""
                                },
                                site.recognition
                            ))
                            .accessibility_label(format!(
                                "Show {} site {} in sequence",
                                e.name, site.recognition
                            ))
                            .on_click(cx.listener(move |_, _, _, cx| {
                                cx.emit(EnzymeSidebarEvent::SelectSite(start))
                            })),
                    )
                    .child(label(site.strand))
                    .child(div().text_sm().child(site.cuts.clone()))
                    .child(
                        div()
                            .font_family("monospace")
                            .text_xs()
                            .text_color(rgb(MUTED))
                            .child(format!("Reference 5′ {} 3′", site.bases)),
                    ),
            );
        }
        let same = ENZYMES
            .iter()
            .enumerate()
            .filter(|(_, candidate)| {
                candidate.name != e.name && candidate.recognition == e.recognition
            })
            .map(|(index, candidate)| {
                Button::new(("same-recognition-enzyme", index))
                    .small()
                    .ghost()
                    .label(candidate.name)
                    .accessibility_label(format!("Inspect {} restriction enzyme", candidate.name))
                    .on_click(
                        cx.listener(move |_, _, _, cx| cx.emit(EnzymeSidebarEvent::Inspect(index))),
                    )
            })
            .collect::<Vec<_>>();
        let mut body = div().flex().flex_col().gap_3().p_4()
            .child(div().flex().items_center().gap_2()
                .child(div().text_2xl().font_weight(FontWeight::BOLD).text_color(rgb(ACCENT)).child(e.name))
                .child(div().flex_1())
                .child(Button::new("copy-enzyme-motif").small().ghost().label("Copy motif")
                    .accessibility_label("Copy enzyme recognition motif")
                    .on_click(move |_, _, cx| cx.write_to_clipboard(ClipboardItem::new_string(e.recognition.into())))))
            .child(label(format!("{} bp{} · {} recognition", e.recognition_length(), if e.recognition.contains('|') { " minimum; alternative motifs" } else { "" }, if palindrome(e) { "Palindromic" } else { "Non-palindromic" })))
            .child(div().font_weight(FontWeight::SEMIBOLD).child(overhang(e)))
            .child(div().p_3().bg(rgb(PANEL)).rounded_lg().flex().flex_col().gap_2()
                .children(diagrams)
                .child(label(if e.cleavage_known { "Arrows mark nominal strand cuts." } else { "This catalogue reports recognition only; cut positions are unspecified." })))
            .when(!same.is_empty(), |el| el.child(div().flex().flex_col().gap_1()
                .child(label("Same recognition in this catalogue"))
                .child(div().flex().flex_wrap().gap_1().children(same))
                .child(label("Cleavage positions may differ."))))
            .when(e.name == "NmeAIII", |el| el.child(div().text_sm().text_color(rgb(0xe7bd7b)).child("NmeAIII nominal cuts can shift by one base with sequence context.")))
            .child(div().border_t_1().border_color(rgb(BORDER)).pt_3().flex().flex_col().gap_1()
                .child(div().font_weight(FontWeight::SEMIBOLD).child(d.document_name.clone()))
                .child(label(format!("{} bp · {}", number(d.document_length), if d.circular { "circular" } else { "linear" })))
                .child(div().mt_2().child(format!("{} recognition {}", number(d.count), if d.count == 1 { "site" } else { "sites" }))))
            .when(d.ambiguous, |el| el.child(div().text_sm().text_color(rgb(0xe7bd7b)).child("Ambiguous bases are excluded from matching; site counts may be incomplete.")))
            .when(d.count > d.sites.len(), |el| el.child(div().text_sm().text_color(rgb(0xe7bd7b)).child(format!("Showing {} retained positions of {} matches. Other positions were not retained.", number(d.sites.len()), number(d.count)))))
            .child(label("Click a site to show its sequence. Ranges are 1-based, inclusive; cuts are after the indicated base on the forward document axis."))
            .child(div().flex().flex_col().gap_2().children(site_rows));
        if d.sites.is_empty() {
            body = body.child(label("No definite recognition positions in this document."));
        } else if d.sites.len() > PAGE_SIZE {
            body = body.child(
                div()
                    .flex()
                    .flex_wrap()
                    .items_center()
                    .gap_2()
                    .child(label(format!(
                        "{}–{} of {} retained positions",
                        page_start + 1,
                        page_end,
                        number(d.sites.len())
                    )))
                    .child(
                        Button::new("enzyme-sites-previous")
                            .small()
                            .label("Previous sites")
                            .disabled(self.page == 0)
                            .on_click(cx.listener(|this, _, _, cx| {
                                this.page = this.page.saturating_sub(1);
                                cx.notify();
                            })),
                    )
                    .child(
                        Button::new("enzyme-sites-next")
                            .small()
                            .label("Next sites")
                            .disabled(page_end >= d.sites.len())
                            .on_click(cx.listener(|this, _, _, cx| {
                                if (this.page + 1) * PAGE_SIZE < this.details.sites.len() {
                                    this.page += 1;
                                    cx.notify();
                                }
                            })),
                    ),
            );
        }
        body = body.child(label("Reference catalogue: Biopython 1.88 / REBASE 404 (2024), plus EcoP15I vendor data. Recognition and nominal cleavage shown here; experimental activity depends on methylation, buffer and enzyme mechanism."));
        div()
            .id("enzyme-details-sidebar")
            .size_full()
            .flex()
            .flex_col()
            .bg(rgb(BG))
            .border_l_1()
            .border_color(rgb(BORDER))
            .text_color(rgb(INK))
            .font_family(cx.theme().font_family.clone())
            .text_sm()
            .key_context("EnzymeDetails")
            .track_focus(&self.focus)
            .on_action(cx.listener(|_, _: &CloseDetails, _, cx| cx.emit(EnzymeSidebarEvent::Close)))
            .child(
                div()
                    .flex()
                    .items_center()
                    .gap_2()
                    .p_3()
                    .border_b_1()
                    .border_color(rgb(BORDER))
                    .child(
                        div()
                            .font_weight(FontWeight::SEMIBOLD)
                            .child("Restriction enzymes"),
                    )
                    .child(div().flex_1())
                    .child(
                        Button::new("close-enzyme-sidebar")
                            .small()
                            .ghost()
                            .label("Close")
                            .accessibility_label("Close enzyme sidebar")
                            .on_click(
                                cx.listener(|_, _, _, cx| cx.emit(EnzymeSidebarEvent::Close)),
                            ),
                    ),
            )
            .child(catalogue)
            .child(
                div()
                    .id("enzyme-details-scroll")
                    .flex_1()
                    .min_h_0()
                    .overflow_y_scroll()
                    .track_scroll(&self.scroll)
                    .child(body),
            )
    }
}

#[cfg(test)]
mod tests {
    use super::{
        Details, EndFilter, adjacent, catalogue_indices, cut_label, overhang, palindrome,
        strand_diagram,
    };
    use popgene_core::{
        Document,
        restriction::{ENZYMES, Enzyme, find_sites},
    };
    fn enzyme(name: &str) -> &'static Enzyme {
        ENZYMES.iter().find(|enzyme| enzyme.name == name).unwrap()
    }
    #[test]
    fn catalogue_search_and_end_filters_cover_all_enzymes() {
        let all = catalogue_indices("", EndFilter::All);
        assert_eq!(all.len(), ENZYMES.len());
        let names = |indices: Vec<usize>| {
            indices
                .into_iter()
                .map(|i| ENZYMES[i].name)
                .collect::<Vec<_>>()
        };
        assert_eq!(
            names(catalogue_indices("  eCoRi  ", EndFilter::All)),
            vec!["EcoRI", "EcoRII"]
        );
        assert!(names(catalogue_indices("GAATTC", EndFilter::All)).contains(&"EcoRI"));
        assert_eq!(
            names(catalogue_indices("kpn", EndFilter::ThreePrime)),
            vec!["KpnI"]
        );
        assert!(catalogue_indices("KpnI", EndFilter::FivePrime).is_empty());
        assert_eq!(
            names(catalogue_indices("EcoRV", EndFilter::Blunt)),
            vec!["EcoRV"]
        );
        assert!(catalogue_indices("does-not-exist", EndFilter::All).is_empty());
        let filtered_count = [
            EndFilter::FivePrime,
            EndFilter::ThreePrime,
            EndFilter::Blunt,
        ]
        .into_iter()
        .map(|f| catalogue_indices("", f).len())
        .sum::<usize>();
        assert_eq!(
            filtered_count,
            ENZYMES.iter().filter(|e| e.cleavage_known).count()
        );
    }
    #[test]
    fn catalogue_navigation_respects_filters_and_boundaries() {
        let indices = vec![4, 9, 1];
        assert_eq!(adjacent(&indices, 9, true), Some(1));
        assert_eq!(adjacent(&indices, 9, false), Some(4));
        assert_eq!(adjacent(&indices, 4, false), None);
        assert_eq!(adjacent(&indices, 1, true), None);
        assert_eq!(adjacent(&indices, 7, true), Some(4));
        assert_eq!(adjacent(&indices, 7, false), Some(1));
        assert_eq!(adjacent(&[], 7, true), None);
    }
    #[test]
    fn complementary_diagrams_and_overhang_polarity() {
        let eco = enzyme("EcoRI");
        assert!(palindrome(eco));
        assert_eq!(strand_diagram(eco.recognition, eco.top_cut, '↓'), "G↓AATTC");
        assert_eq!(strand_diagram("CTTAAG", eco.bottom_cut, '↑'), "CTTAA↑G");
        assert_eq!(overhang(eco), "4-nt 5′ overhang");
        assert_eq!(overhang(enzyme("KpnI")), "4-nt 3′ overhang");
        assert_eq!(overhang(enzyme("EcoRV")), "Blunt ends");
        assert!(!palindrome(enzyme("NmeAIII")));
        assert_eq!(strand_diagram("GCCGAG", 27, '↓'), "GCCGAG … ↓ (+21)");
        assert!(palindrome(enzyme("XcmI")));
    }
    #[test]
    fn snapshot_retains_circular_ranges_and_bounded_positions() {
        let mut document = Document::new("Origin", "TTCGAA").unwrap();
        document.circular = true;
        let analysis = find_sites(&document.sequence, true);
        let index = ENZYMES
            .iter()
            .position(|enzyme| enzyme.name == "EcoRI")
            .unwrap();
        let details = Details::from_document(&document, &analysis, index, Some(3), false).unwrap();
        assert_eq!(details.count, 1);
        assert_eq!(details.selected, Some(0));
        assert_eq!(details.sites[0].recognition, "4–6 / 1–3 · crosses origin");
        assert_eq!(details.sites[0].bases, "GAATTC");
        assert_eq!(cut_label(Some(0), 6, true), "after 6");
        assert_eq!(cut_label(Some(0), 6, false), "before base 1");
        assert_eq!(cut_label(None, 6, false), "outside molecule");
        assert!(Details::from_document(&document, &analysis, ENZYMES.len(), None, false).is_none());
        let repeated = Document::new("Repeated", &"GAATTC".repeat(600)).unwrap();
        let analysis = find_sites(&repeated.sequence, false);
        let details = Details::from_document(&repeated, &analysis, index, None, false).unwrap();
        assert_eq!(details.count, 600);
        assert_eq!(details.sites.len(), 512);
    }
    #[test]
    fn reverse_site_uses_document_coordinates_and_reference_bases() {
        let document = Document::new(
            "Reverse",
            &format!("{}CTCGGC{}", "A".repeat(30), "A".repeat(30)),
        )
        .unwrap();
        let analysis = find_sites(&document.sequence, false);
        let index = ENZYMES
            .iter()
            .position(|enzyme| enzyme.name == "NmeAIII")
            .unwrap();
        let details = Details::from_document(&document, &analysis, index, Some(30), false).unwrap();
        assert_eq!(details.sites[0].recognition, "31–36");
        assert_eq!(details.sites[0].bases, "CTCGGC");
        assert_eq!(details.sites[0].strand, "Reverse match");
        assert_eq!(details.sites[0].cuts, "Forward after 11; opposite after 9");
    }
}
