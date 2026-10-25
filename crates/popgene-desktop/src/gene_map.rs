//! Native GPUI vector maps. Geometry, hit testing and scientific coordinates share one scene.
use crate::feature_card::{direction, feature_color, number, render_feature_card, truncate};
use gpui_kit::{prelude::*, *};
use popgene_core::{Document, Feature, Segment};
use std::{
    cell::RefCell,
    f32::consts::{FRAC_PI_2, TAU},
    rc::Rc,
    time::Duration,
};

const INK: u32 = 0xe2eaf4;
const MUTED: u32 = 0x8394ad;
const ACCENT: u32 = 0x58d4bb;
#[derive(Clone, Copy, Debug, PartialEq)]
struct P {
    x: f32,
    y: f32,
}
impl P {
    fn new(x: f32, y: f32) -> Self {
        Self { x, y }
    }
    fn distance(self, b: Self) -> f32 {
        (self.x - b.x).hypot(self.y - b.y)
    }
    fn offset(self, x: f32, y: f32) -> Self {
        Self::new(self.x + x, self.y + y)
    }
    fn gpui(self) -> Point<Pixels> {
        point(px(self.x), px(self.y))
    }
}
#[derive(Clone, Copy, Debug)]
struct Rect {
    x: f32,
    y: f32,
    w: f32,
    h: f32,
}
impl Rect {
    fn contains(self, p: P, padding: f32) -> bool {
        p.x >= self.x - padding
            && p.x <= self.x + self.w + padding
            && p.y >= self.y - padding
            && p.y <= self.y + self.h + padding
    }
    fn center(self) -> P {
        P::new(self.x + self.w * 0.5, self.y + self.h * 0.5)
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct FeatureHit {
    index: usize,
    segment: Option<usize>,
}
#[derive(Clone, Debug)]
enum HitShape {
    Arc {
        center: P,
        radius: f32,
        start: f32,
        end: f32,
    },
    Box(Rect),
    Line(P, P),
}
impl HitShape {
    fn contains(&self, p: P) -> bool {
        match *self {
            Self::Box(rect) => rect.contains(p, 4.),
            Self::Line(a, b) => {
                let dx = b.x - a.x;
                let dy = b.y - a.y;
                let scale = ((p.x - a.x) * dx + (p.y - a.y) * dy) / (dx * dx + dy * dy).max(0.001);
                p.distance(a.offset(dx * scale.clamp(0., 1.), dy * scale.clamp(0., 1.))) <= 3.
            }
            Self::Arc {
                center,
                radius,
                start,
                end,
            } => {
                let fraction =
                    ((p.y - center.y).atan2(p.x - center.x) + FRAC_PI_2).rem_euclid(TAU) / TAU;
                ((p.distance(center) - radius).abs() <= 6. && fraction >= start && fraction <= end)
                    || p.distance(on_circle(center, radius, start)) <= 6.
                    || p.distance(on_circle(center, radius, end)) <= 6.
            }
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Hit {
    Feature(FeatureHit),
    Site(usize),
}
impl Hit {
    fn same_object(self, other: Self) -> bool {
        match (self, other) {
            (Self::Feature(a), Self::Feature(b)) => a.index == b.index,
            (Self::Site(a), Self::Site(b)) => a == b,
            _ => false,
        }
    }
}
#[derive(Clone)]
struct Region {
    hit: Hit,
    shape: HitShape,
}
fn hit_test(regions: &[Region], point: P) -> Option<Hit> {
    let is_label = |region: &&Region| match (region.hit, &region.shape) {
        (Hit::Feature(_), HitShape::Box(rect)) => rect.h > 12.,
        (Hit::Site(_), HitShape::Box(rect)) => rect.h > 12.,
        _ => false,
    };
    regions
        .iter()
        .rev()
        .filter(is_label)
        .chain(regions.iter().rev())
        .find(|r| r.shape.contains(point))
        .map(|r| r.hit)
}
fn on_circle(center: P, radius: f32, fraction: f32) -> P {
    let angle = fraction * TAU - FRAC_PI_2;
    center.offset(angle.cos() * radius, angle.sin() * radius)
}
enum Draw {
    Stroke(Vec<P>, f32, Hsla),
    Fill(Vec<P>, Hsla),
    Text {
        point: P,
        text: String,
        size: f32,
        color: Hsla,
    },
}
#[derive(Default)]
struct Scene {
    draws: Vec<Draw>,
    regions: Vec<Region>,
    connectors: Vec<(Hit, P, P)>,
    notices: Vec<String>,
}
impl Scene {
    fn line(&mut self, a: P, b: P, width: f32, color: Hsla) {
        self.draws.push(Draw::Stroke(vec![a, b], width, color));
    }
    fn arc(&mut self, center: P, radius: f32, start: f32, end: f32, width: f32, color: Hsla) {
        let steps = (((end - start).abs() * TAU * radius / 2.).ceil() as usize).clamp(2, 720);
        let points = (0..=steps)
            .map(|i| {
                on_circle(
                    center,
                    radius,
                    start + (end - start) * i as f32 / steps as f32,
                )
            })
            .collect();
        self.draws.push(Draw::Stroke(points, width, color));
    }
    fn text(
        &mut self,
        text: String,
        at: P,
        size: f32,
        color: Hsla,
        align: f32,
        measure: &impl Fn(&str, f32) -> f32,
    ) -> Rect {
        let width = measure(&text, size);
        let rect = Rect {
            x: at.x - width * align,
            y: at.y - size * 0.7,
            w: width,
            h: size * 1.4,
        };
        self.draws.push(Draw::Text {
            point: P::new(rect.x, rect.y),
            text,
            size,
            color,
        });
        rect
    }
    fn highlight(&mut self, shape: &HitShape, color: Hsla) {
        match *shape {
            HitShape::Arc {
                center,
                radius,
                start,
                end,
            } => self.arc(center, radius, start, end, 2., color),
            HitShape::Line(a, b) => self.line(a, b, 1.5, color),
            HitShape::Box(rect) => {
                let a = P::new(rect.x - 3., rect.y - 3.);
                let b = P::new(rect.x + rect.w + 3., rect.y + rect.h + 3.);
                self.draws.push(Draw::Stroke(
                    vec![a, P::new(b.x, a.y), b, P::new(a.x, b.y), a],
                    1.2,
                    color,
                ));
            }
        }
    }
}

pub const MAX_FEATURES: usize = 200;
pub const MAX_PRIMERS: usize = 100;
pub const MAX_FEATURE_LABELS: usize = 100;
pub const MAX_SITE_MARKS: usize = 2_000;
pub const MAX_SITE_LABELS: usize = 80;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MapOptions {
    pub linear: bool,
    pub labels: bool,
    pub features: bool,
    pub primers: bool,
    pub positions: bool,
    pub sites: bool,
    /// 1.0 fits the available viewport. Larger values enlarge the scrollable canvas.
    pub zoom: f32,
}
impl Default for MapOptions {
    fn default() -> Self {
        Self {
            linear: false,
            labels: true,
            features: true,
            primers: true,
            positions: true,
            sites: true,
            zoom: 1.,
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MapRestrictionSite {
    pub label: String,
    /// Zero-based cut boundary: p bases precede the cut (not a residue index).
    pub position: usize,
    pub start: usize,
    pub length: usize,
    pub reverse: bool,
    pub enzyme_index: usize,
}
fn site_parts(site: &MapRestrictionSite, n: usize, circular: bool) -> Vec<Segment> {
    if n == 0 || site.start >= n || site.length == 0 || site.length > n {
        return vec![];
    }
    let Some(end) = site.start.checked_add(site.length) else {
        return vec![];
    };
    if end <= n {
        vec![Segment {
            start: site.start,
            end,
        }]
    } else if circular {
        vec![
            Segment {
                start: site.start,
                end: n,
            },
            Segment {
                start: 0,
                end: end - n,
            },
        ]
    } else {
        vec![]
    }
}
fn cut_label(site: &MapRestrictionSite, n: usize, circular: bool) -> String {
    if circular && n > 0 {
        number(if site.position.is_multiple_of(n) {
            n
        } else {
            site.position % n
        })
    } else if site.position == 0 {
        "before 1".into()
    } else {
        number(site.position)
    }
}
fn major_ticks(n: usize) -> Vec<usize> {
    if n == 0 {
        return vec![];
    }
    let target = (n as f64 / 5.).max(1.);
    let magnitude = 10_f64.powf(target.log10().floor());
    let step = [1., 2., 5., 10.]
        .into_iter()
        .map(|v| v * magnitude)
        .find(|v| *v >= target)
        .unwrap_or(target)
        .round()
        .max(1.) as usize;
    let mut ticks: Vec<_> = (step..=n).step_by(step).collect();
    if ticks.last().copied() != Some(n) {
        if ticks.last().is_some_and(|last| n - *last < step / 3) {
            ticks.pop();
        }
        ticks.push(n);
    }
    ticks
}
fn features_overlap(a: &Feature, b: &Feature) -> bool {
    a.parts
        .iter()
        .any(|a| b.parts.iter().any(|b| a.start < b.end && b.start < a.end))
}
/// One stable lane for every ordered, possibly origin-spanning feature.
fn assign_tracks(doc: &Document, indices: &[usize]) -> Vec<usize> {
    let mut lanes: Vec<Vec<usize>> = vec![];
    indices
        .iter()
        .map(|&index| {
            let lane = lanes
                .iter()
                .position(|members| {
                    members
                        .iter()
                        .all(|&other| !features_overlap(&doc.features[index], &doc.features[other]))
                })
                .unwrap_or(lanes.len());
            if lane == lanes.len() {
                lanes.push(vec![]);
            }
            lanes[lane].push(index);
            lane
        })
        .collect()
}
fn feature_indices(doc: &Document, primer: bool, enabled: bool, cap: usize) -> Vec<usize> {
    if !enabled {
        return vec![];
    }
    doc.features
        .iter()
        .enumerate()
        .filter(|(_, f)| (f.kind == "primer_bind") == primer && !f.parts.is_empty())
        .take(cap)
        .map(|(i, _)| i)
        .collect()
}
fn feature_anchor(feature: &Feature, n: f32) -> f32 {
    let part = feature
        .parts
        .iter()
        .max_by_key(|p| p.end - p.start)
        .unwrap();
    (part.start as f32 + (part.end - part.start) as f32 * 0.5) / n
}
struct FeatureMapLabel {
    text: String,
    segment: Option<usize>,
    fraction: f32,
    color: Hsla,
}
/// Named display segments get labels at their own coordinates. The parent name
/// remains in the feature card and table; features without names keep one label.
fn feature_map_labels(feature: &Feature, n: f32, color: Hsla) -> Vec<FeatureMapLabel> {
    let segments = feature.display_segments();
    if segments.iter().any(|segment| segment.name.is_some()) {
        segments
            .iter()
            .enumerate()
            .map(|(index, segment)| FeatureMapLabel {
                text: truncate(segment.name.as_deref().unwrap_or(&feature.label), 28),
                segment: Some(index),
                fraction: (segment.range.start as f32 + segment.range.end as f32) * 0.5 / n,
                color: segment
                    .color_rgb
                    .map(|value| rgb(value).into())
                    .unwrap_or(color),
            })
            .collect()
    } else {
        vec![FeatureMapLabel {
            text: truncate(&feature.label, 28),
            segment: None,
            fraction: feature_anchor(feature, n),
            color,
        }]
    }
}
fn primer_label(feature: &Feature) -> String {
    let arrow = if feature.reverse { "←" } else { "→" };
    let parts = feature
        .parts
        .iter()
        .take(3)
        .map(|p| format!("{}–{}", number(p.start + 1), number(p.end)))
        .collect::<Vec<_>>()
        .join(", ");
    format!("{parts} {arrow}")
}
#[derive(Clone)]
struct Callout {
    hit: Hit,
    text: String,
    detail: Option<String>,
    anchor: P,
    color: Hsla,
    priority: usize,
}
impl Callout {
    fn height(&self) -> f32 {
        if self.detail.is_some() { 32. } else { 19. }
    }
    fn width(&self, measure: &impl Fn(&str, f32) -> f32) -> f32 {
        measure(&self.text, 12.).max(self.detail.as_ref().map_or(0., |d| measure(d, 10.)))
    }
}
struct MapData<'a> {
    doc: &'a Document,
    selected: &'a [Segment],
    options: MapOptions,
    sites: &'a [MapRestrictionSite],
    hovered: Option<Hit>,
}
fn rects_overlap(a: Rect, b: Rect, pad: f32) -> bool {
    a.x < b.x + b.w + pad && b.x < a.x + a.w + pad && a.y < b.y + b.h + pad && b.y < a.y + a.h + pad
}
fn circular_arrow(
    scene: &mut Scene,
    feature: &Feature,
    doc: &Document,
    center: P,
    r: f32,
    color: Hsla,
) {
    if feature.parts.is_empty() {
        return;
    }
    for reverse in [false, true] {
        let dir = direction(feature, doc.molecule);
        if dir != "Bidirectional"
            && dir
                != if reverse {
                    "Reverse strand"
                } else {
                    "Forward strand"
                }
        {
            continue;
        }
        let end = if reverse {
            feature.parts[0].start
        } else {
            feature.parts.last().unwrap().end
        };
        let display_segments = feature.display_segments();
        // Ordered endpoints identify the terminal occurrence even when two
        // biological parts share the same genomic coordinates.
        let color = if reverse {
            display_segments.first()
        } else {
            display_segments.last()
        }
        .and_then(|segment| segment.color_rgb)
        .map(|value| rgb(value).into())
        .unwrap_or(color);
        let fraction = end as f32 / doc.sequence.len() as f32;
        let angle = fraction * TAU - FRAC_PI_2;
        let tip = on_circle(center, r, fraction);
        let sign = if reverse { -1. } else { 1. };
        let t = P::new(-angle.sin() * sign, angle.cos() * sign);
        let radial = P::new(angle.cos(), angle.sin());
        scene.draws.push(Draw::Fill(
            vec![
                tip.offset(t.x * 7., t.y * 7.),
                tip.offset(-t.x * 5. + radial.x * 7., -t.y * 5. + radial.y * 7.),
                tip.offset(-t.x * 5. - radial.x * 7., -t.y * 5. - radial.y * 7.),
            ],
            color,
        ));
    }
}
/// Labels sit on the dark canvas, so retain their hue at a readable lightness.
/// Geometry and imported color values continue to use the exact source color.
fn label_color(mut color: Hsla) -> Hsla {
    color.l = color.l.max(0.72);
    color
}
fn fit_text(text: &str, size: f32, width: f32, measure: &impl Fn(&str, f32) -> f32) -> String {
    if measure(text, size) <= width {
        return text.to_owned();
    }
    let mut chars = text.chars().take(96).collect::<Vec<_>>();
    while !chars.is_empty() {
        chars.pop();
        let text = format!("{}…", chars.iter().collect::<String>());
        if measure(&text, size) <= width {
            return text;
        }
    }
    "…".into()
}
fn pack_callouts(
    scene: &mut Scene,
    rect: Rect,
    center: P,
    radius: f32,
    callouts: Vec<Callout>,
    measure: &impl Fn(&str, f32) -> f32,
) -> (usize, usize, usize) {
    let mut counts = (0, 0, 0);
    for right in [false, true] {
        let mut side: Vec<_> = callouts
            .iter()
            .filter(|c| (c.anchor.x >= center.x) == right)
            .cloned()
            .collect();
        // Primer bindings and short feature labels survive dense enzyme lists.
        side.sort_by_key(|c| c.priority);
        let available = (rect.h - 68.).max(30.);
        let mut used = 0.;
        side.retain(|c| {
            if used + c.height() + 3. <= available {
                used += c.height() + 3.;
                true
            } else {
                false
            }
        });
        side.sort_by(|a, b| a.anchor.y.total_cmp(&b.anchor.y));
        let top = rect.y + 24.;
        let bottom = rect.y + rect.h - 44.;
        let mut positions: Vec<f32> = vec![];
        for (i, c) in side.iter().enumerate() {
            let minimum = if i == 0 {
                top + c.height() * 0.5
            } else {
                positions[i - 1] + (side[i - 1].height() + c.height()) * 0.5 + 3.
            };
            positions.push(c.anchor.y.max(minimum));
        }
        for i in (0..side.len()).rev() {
            let maximum = if i + 1 == side.len() {
                bottom - side[i].height() * 0.5
            } else {
                positions[i + 1] - (side[i + 1].height() + side[i].height()) * 0.5 - 3.
            };
            positions[i] = positions[i].min(maximum);
        }
        for (mut c, y) in side.into_iter().zip(positions) {
            let gutter = (rect.w * 0.5 - radius - 36.).max(56.);
            c.text = fit_text(&c.text, 12., gutter, measure);
            if let Some(detail) = &mut c.detail {
                *detail = fit_text(detail, 10., gutter, measure);
            }
            let width = c.width(measure);
            let x = if right {
                (center.x + radius + 32.).min(rect.x + rect.w - width - 12.)
            } else {
                (center.x - radius - 32.).max(rect.x + width + 12.)
            };
            let elbow = P::new(x, y);
            let end = elbow.offset(if right { 6. } else { -6. }, 0.);
            let text_y = y - if c.detail.is_some() { 6. } else { 0. };
            scene.line(c.anchor, elbow, 1., c.color.opacity(0.5));
            scene.line(elbow, end, 1., c.color.opacity(0.5));
            let mut bounds = scene.text(
                c.text,
                end.offset(0., text_y - y),
                12.,
                label_color(c.color),
                if right { 0. } else { 1. },
                measure,
            );
            if let Some(detail) = c.detail {
                let lower = scene.text(
                    detail,
                    end.offset(0., 10.),
                    10.,
                    label_color(c.color).opacity(0.85),
                    if right { 0. } else { 1. },
                    measure,
                );
                bounds.x = bounds.x.min(lower.x);
                bounds.w = width;
                bounds.h = lower.y + lower.h - bounds.y;
            }
            scene.regions.push(Region {
                hit: c.hit,
                shape: HitShape::Box(bounds),
            });
            scene.connectors.push((c.hit, c.anchor, elbow));
            scene.regions.push(Region {
                hit: c.hit,
                shape: HitShape::Line(c.anchor, elbow),
            });
            match c.hit {
                Hit::Site(_) => counts.0 += 1,
                Hit::Feature(_) => {
                    if c.priority == 0 {
                        counts.1 += 1
                    } else {
                        counts.2 += 1
                    }
                }
            }
        }
    }
    counts
}
fn circular_scene(rect: Rect, data: &MapData<'_>, measure: &impl Fn(&str, f32) -> f32) -> Scene {
    let mut scene = Scene::default();
    let doc = data.doc;
    let n = doc.sequence.len() as f32;
    let options = data.options;
    let center = rect.center();
    let features = feature_indices(doc, false, options.features, MAX_FEATURES);
    let primers = feature_indices(doc, true, options.primers, MAX_PRIMERS);
    let sites: Vec<_> = if options.sites {
        data.sites
            .iter()
            .enumerate()
            .filter(|(_, s)| !site_parts(s, n as usize, true).is_empty())
            .take(MAX_SITE_MARKS)
            .collect()
    } else {
        vec![]
    };
    let site_text = |s: &MapRestrictionSite| {
        format!(
            "{}  ({})",
            truncate(&s.label, 20),
            cut_label(s, n as usize, true)
        )
    };
    let label_width = if options.labels {
        sites
            .iter()
            .take(MAX_SITE_LABELS)
            .map(|(_, s)| measure(&site_text(s), 12.))
            .chain(primers.iter().map(|&i| {
                measure(&truncate(&doc.features[i].label, 24), 12.)
                    .max(measure(&primer_label(&doc.features[i]), 10.))
            }))
            .fold(80_f32, f32::max)
            .min(rect.w * 0.3)
    } else {
        16.
    };
    let minimum_radius = (rect.w * 0.24).min(rect.h * 0.35).max(48.);
    let radius = (rect.w * 0.5 - label_width - 42.)
        .max(minimum_radius)
        .min(rect.h * 0.5 - 48.)
        .max(48.);
    scene.arc(center, radius, 0., 1., 2., rgb(0xb4c0cf).into());
    scene.arc(center, radius - 5., 0., 1., 1., rgb(0x71839a).into());
    if options.positions {
        let ticks = major_ticks(n as usize);
        let minor = ticks.first().copied().unwrap_or(1).div_ceil(5).max(1);
        for p in (0..n as usize).step_by(minor) {
            let f = p as f32 / n;
            scene.line(
                on_circle(center, radius - 6., f),
                on_circle(center, radius + 1., f),
                1.,
                rgb(0x63748d).into(),
            );
        }
        for p in ticks {
            let f = p as f32 / n;
            scene.line(
                on_circle(center, radius - 12., f),
                on_circle(center, radius + 3., f),
                1.5,
                rgb(MUTED).into(),
            );
            scene.text(
                number(p),
                on_circle(center, radius - 23., f),
                11.,
                rgb(INK).into(),
                0.5,
                measure,
            );
        }
    }
    let name = truncate(&doc.name, ((radius * 1.5 / 8.) as usize).clamp(8, 42));
    let center_width = measure(&name, 17.);
    let center_bounds = Rect {
        x: center.x - center_width * 0.5 - 8.,
        y: center.y - 36.,
        w: center_width + 16.,
        h: 84.,
    };
    scene.text(
        name,
        center.offset(0., -13.),
        17.,
        rgb(INK).into(),
        0.5,
        measure,
    );
    scene.text(
        format!("{} {}", number(n as usize), doc.unit()),
        center.offset(0., 13.),
        13.,
        rgb(MUTED).into(),
        0.5,
        measure,
    );
    scene.text(
        format!("CIRCULAR {}", doc.molecule.label().to_uppercase()),
        center.offset(0., 35.),
        10.,
        rgb(ACCENT).into(),
        0.5,
        measure,
    );
    for part in data.selected {
        scene.arc(
            center,
            radius + 5.,
            part.start as f32 / n,
            part.end as f32 / n,
            4.,
            rgb(ACCENT).into(),
        );
    }
    let tracks = assign_tracks(doc, &features);
    let track_capacity = ((radius - 74.) / 21.).floor().max(1.) as usize;
    let mut callouts = vec![];
    let mut occupied = vec![center_bounds];
    let mut drawn_features = 0;
    let mut feature_labels = 0;
    let mut label_candidates = 0;
    let mut expected_feature_labels = 0;
    for (&index, &track) in features.iter().zip(&tracks) {
        if track >= track_capacity {
            continue;
        }
        let f = &doc.features[index];
        let r = radius - 43. - track as f32 * 21.;
        let color = feature_color(f, index);
        drawn_features += 1;
        for (segment, display) in f.display_segments().iter().enumerate() {
            let p = &display.range;
            let segment_color = display
                .color_rgb
                .map(|value| rgb(value).into())
                .unwrap_or(color);
            scene.arc(
                center,
                r,
                p.start as f32 / n,
                p.end as f32 / n,
                10.,
                segment_color,
            );
            scene.regions.push(Region {
                hit: Hit::Feature(FeatureHit {
                    index,
                    segment: Some(segment),
                }),
                shape: HitShape::Arc {
                    center,
                    radius: r,
                    start: p.start as f32 / n,
                    end: p.end as f32 / n,
                },
            });
        }
        circular_arrow(&mut scene, f, doc, center, r, color);
        let map_labels = feature_map_labels(f, n, color);
        expected_feature_labels += map_labels.len();
        for label in map_labels {
            if !options.labels || label_candidates >= MAX_FEATURE_LABELS {
                break;
            }
            label_candidates += 1;
            let FeatureMapLabel {
                text: label,
                segment,
                fraction,
                color,
            } = label;
            let width = measure(&label, 12.);
            let anchor = on_circle(center, r, fraction);
            let mut placed = false;
            for offset in [15., 32., 49.] {
                for angle in [0., -0.025, 0.025] {
                    let at = on_circle(center, (r - offset).max(24.), fraction + angle);
                    let bounds = Rect {
                        x: at.x - width * 0.5,
                        y: at.y - 8.4,
                        w: width,
                        h: 16.8,
                    };
                    let corners = [
                        P::new(bounds.x, bounds.y),
                        P::new(bounds.x + bounds.w, bounds.y),
                        P::new(bounds.x, bounds.y + bounds.h),
                        P::new(bounds.x + bounds.w, bounds.y + bounds.h),
                    ];
                    if corners.iter().all(|p| p.distance(center) < radius - 8.)
                        && !occupied.iter().any(|b| rects_overlap(*b, bounds, 3.))
                    {
                        scene.text(label.clone(), at, 12., label_color(color), 0.5, measure);
                        scene.regions.push(Region {
                            hit: Hit::Feature(FeatureHit { index, segment }),
                            shape: HitShape::Box(bounds),
                        });
                        occupied.push(bounds);
                        feature_labels += 1;
                        placed = true;
                        break;
                    }
                }
                if placed {
                    break;
                }
            }
            if !placed {
                callouts.push(Callout {
                    hit: Hit::Feature(FeatureHit { index, segment }),
                    text: label,
                    detail: None,
                    anchor,
                    color,
                    priority: 1,
                });
            }
        }
    }
    let primer_tracks = assign_tracks(doc, &primers);
    let mut drawn_primers = 0;
    for (&index, &track) in primers.iter().zip(&primer_tracks) {
        if track >= 4 {
            continue;
        }
        drawn_primers += 1;
        let f = &doc.features[index];
        let r = radius + 13. + track as f32 * 10.;
        let color = f
            .display_color_rgb()
            .map(|value| rgb(value).into())
            .unwrap_or_else(|| rgb(if f.reverse { 0xd895f0 } else { 0xb985f0 }).into());
        for (segment, display) in f.display_segments().iter().enumerate() {
            let p = &display.range;
            let segment_color = display
                .color_rgb
                .map(|value| rgb(value).into())
                .unwrap_or(color);
            scene.arc(
                center,
                r,
                p.start as f32 / n,
                p.end as f32 / n,
                3.,
                segment_color,
            );
            scene.regions.push(Region {
                hit: Hit::Feature(FeatureHit {
                    index,
                    segment: Some(segment),
                }),
                shape: HitShape::Arc {
                    center,
                    radius: r,
                    start: p.start as f32 / n,
                    end: p.end as f32 / n,
                },
            });
        }
        circular_arrow(&mut scene, f, doc, center, r, color);
        if options.labels {
            callouts.push(Callout {
                hit: Hit::Feature(FeatureHit {
                    index,
                    segment: None,
                }),
                text: truncate(&f.label, 24),
                detail: Some(primer_label(f)),
                anchor: on_circle(center, r, feature_anchor(f, n)),
                color,
                priority: 0,
            });
        }
    }
    for (ordinal, (index, s)) in sites.iter().enumerate() {
        let fraction = (s.position % (n as usize)) as f32 / n;
        let anchor = on_circle(center, radius + 2., fraction);
        let tip = on_circle(center, radius + 15., fraction);
        scene.line(anchor, tip, 1.2, rgb(0xa9b5c5).into());
        scene.regions.push(Region {
            hit: Hit::Site(*index),
            shape: HitShape::Box(Rect {
                x: tip.x - 3.,
                y: tip.y - 3.,
                w: 6.,
                h: 6.,
            }),
        });
        if options.labels && ordinal < MAX_SITE_LABELS {
            callouts.push(Callout {
                hit: Hit::Site(*index),
                text: site_text(s),
                detail: None,
                anchor: tip,
                color: rgb(0xd9e2ee).into(),
                priority: 2,
            });
        }
    }
    let (site_labels, primer_labels, external_feature_labels) =
        pack_callouts(&mut scene, rect, center, radius, callouts, measure);
    feature_labels += external_feature_labels;
    add_notices(
        &mut scene,
        data,
        drawn_features,
        drawn_primers,
        feature_labels,
        expected_feature_labels,
        primer_labels,
        sites.len(),
        site_labels,
    );
    scene
}

fn linear_arrow(
    scene: &mut Scene,
    f: &Feature,
    doc: &Document,
    left: f32,
    width: f32,
    y: f32,
    color: Hsla,
) {
    let dir = direction(f, doc.molecule);
    for reverse in [false, true] {
        if dir != "Bidirectional"
            && dir
                != if reverse {
                    "Reverse strand"
                } else {
                    "Forward strand"
                }
        {
            continue;
        }
        let end = if reverse {
            f.parts[0].start
        } else {
            f.parts.last().unwrap().end
        };
        let display_segments = f.display_segments();
        // Ordered endpoints identify the terminal occurrence even when two
        // biological parts share the same genomic coordinates.
        let color = if reverse {
            display_segments.first()
        } else {
            display_segments.last()
        }
        .and_then(|segment| segment.color_rgb)
        .map(|value| rgb(value).into())
        .unwrap_or(color);
        let x = left + width * end as f32 / doc.sequence.len() as f32;
        let sign = if reverse { -1. } else { 1. };
        scene.draws.push(Draw::Fill(
            vec![
                P::new(x + sign * 5., y),
                P::new(x - sign * 5., y - 7.),
                P::new(x - sign * 5., y + 7.),
            ],
            color,
        ));
    }
}
fn linear_scene(rect: Rect, data: &MapData<'_>, measure: &impl Fn(&str, f32) -> f32) -> Scene {
    let mut scene = Scene::default();
    let doc = data.doc;
    let options = data.options;
    let n = doc.sequence.len();
    let left = rect.x + 40.;
    let width = (rect.w - 80.).max(30.);
    let y = rect.y + rect.h * 0.45;
    let features = feature_indices(doc, false, options.features, MAX_FEATURES);
    let primers = feature_indices(doc, true, options.primers, MAX_PRIMERS);
    let sites: Vec<_> = if options.sites {
        data.sites
            .iter()
            .enumerate()
            .filter(|(_, s)| !site_parts(s, n, doc.circular).is_empty())
            .take(MAX_SITE_MARKS)
            .collect()
    } else {
        vec![]
    };
    scene.text(
        truncate(&doc.name, 64),
        P::new(left, rect.y + 26.),
        17.,
        rgb(INK).into(),
        0.,
        measure,
    );
    scene.text(
        format!(
            "{} {} · {}",
            number(n),
            doc.unit(),
            if doc.circular {
                "CIRCULAR · LINEAR VIEW"
            } else {
                "LINEAR"
            }
        ),
        P::new(left, rect.y + 49.),
        11.,
        rgb(MUTED).into(),
        0.,
        measure,
    );
    scene.line(
        P::new(left, y),
        P::new(left + width, y),
        2.,
        rgb(0xb4c0cf).into(),
    );
    if options.positions {
        scene.text(
            "1".into(),
            P::new(left, y - 17.),
            11.,
            rgb(MUTED).into(),
            0.,
            measure,
        );
        for p in major_ticks(n) {
            let x = left + width * p as f32 / n as f32;
            scene.line(P::new(x, y - 5.), P::new(x, y + 5.), 1., rgb(MUTED).into());
            scene.text(
                number(p),
                P::new(x, y - 17.),
                11.,
                rgb(MUTED).into(),
                0.5,
                measure,
            );
        }
    }
    for part in data.selected {
        scene.line(
            P::new(left + width * part.start as f32 / n as f32, y + 8.),
            P::new(left + width * part.end as f32 / n as f32, y + 8.),
            4.,
            rgb(ACCENT).into(),
        );
    }
    let primer_tracks = assign_tracks(doc, &primers);
    let mut drawn_primers = 0;
    let mut primer_labels = 0;
    let mut labels = vec![];
    let primer_capacity = ((rect.h * 0.2 - 25.) / 38.).floor().clamp(1., 4.) as usize;
    for (&index, &track) in primers.iter().zip(&primer_tracks) {
        if track >= primer_capacity {
            continue;
        }
        drawn_primers += 1;
        let f = &doc.features[index];
        let at = y + 24. + track as f32 * 38.;
        let color = f
            .display_color_rgb()
            .map(|value| rgb(value).into())
            .unwrap_or_else(|| rgb(if f.reverse { 0xd895f0 } else { 0xb985f0 }).into());
        for (segment, display) in f.display_segments().iter().enumerate() {
            let p = &display.range;
            let segment_color = display
                .color_rgb
                .map(|value| rgb(value).into())
                .unwrap_or(color);
            let a = left + width * p.start as f32 / n as f32;
            let b = left + width * p.end as f32 / n as f32;
            scene.line(P::new(a, at), P::new(b.max(a + 2.), at), 3., segment_color);
            scene.regions.push(Region {
                hit: Hit::Feature(FeatureHit {
                    index,
                    segment: Some(segment),
                }),
                shape: HitShape::Box(Rect {
                    x: a,
                    y: at - 4.,
                    w: (b - a).max(2.),
                    h: 8.,
                }),
            });
        }
        linear_arrow(&mut scene, f, doc, left, width, at, color);
        if options.labels {
            let text = format!("{}  {}", truncate(&f.label, 20), primer_label(f));
            let tw = measure(&text, 11.);
            let x = (left + width * feature_anchor(f, n as f32) - tw * 0.5)
                .clamp(rect.x + 8., (rect.x + rect.w - tw - 8.).max(rect.x + 8.));
            let bounds = Rect {
                x,
                y: at + 7.,
                w: tw,
                h: 15.4,
            };
            if !labels.iter().any(|b| rects_overlap(*b, bounds, 4.)) {
                scene.text(text, P::new(x, at + 14.7), 11., color, 0., measure);
                scene.regions.push(Region {
                    hit: Hit::Feature(FeatureHit {
                        index,
                        segment: None,
                    }),
                    shape: HitShape::Box(bounds),
                });
                labels.push(bounds);
                primer_labels += 1;
            }
        }
    }
    let primer_height = primer_tracks
        .iter()
        .copied()
        .max()
        .map_or(0., |v| ((v + 1).min(primer_capacity)) as f32 * 38.);
    let feature_start = y + 30. + primer_height;
    let feature_tracks = assign_tracks(doc, &features);
    let capacity = ((rect.y + rect.h - 46. - feature_start) / 42.)
        .floor()
        .max(1.) as usize;
    let mut drawn_features = 0;
    let mut feature_labels = 0;
    let mut expected_feature_labels = 0;
    for (&index, &track) in features.iter().zip(&feature_tracks) {
        if track >= capacity {
            continue;
        }
        drawn_features += 1;
        let f = &doc.features[index];
        let at = feature_start + track as f32 * 42.;
        let color = feature_color(f, index);
        for (segment, display) in f.display_segments().iter().enumerate() {
            let p = &display.range;
            let segment_color = display
                .color_rgb
                .map(|value| rgb(value).into())
                .unwrap_or(color);
            let a = left + width * p.start as f32 / n as f32;
            let b = left + width * p.end as f32 / n as f32;
            scene.line(P::new(a, at), P::new(b.max(a + 2.), at), 10., segment_color);
            scene.regions.push(Region {
                hit: Hit::Feature(FeatureHit {
                    index,
                    segment: Some(segment),
                }),
                shape: HitShape::Box(Rect {
                    x: a,
                    y: at - 6.,
                    w: (b - a).max(2.),
                    h: 12.,
                }),
            });
        }
        linear_arrow(&mut scene, f, doc, left, width, at, color);
        let map_labels = feature_map_labels(f, n as f32, color);
        expected_feature_labels += map_labels.len();
        for label in map_labels {
            if !options.labels || feature_labels >= MAX_FEATURE_LABELS {
                break;
            }
            let FeatureMapLabel {
                text,
                segment,
                fraction,
                color,
            } = label;
            let tw = measure(&text, 12.);
            let x = (left + width * fraction - tw * 0.5)
                .clamp(rect.x + 8., (rect.x + rect.w - tw - 8.).max(rect.x + 8.));
            let bounds = Rect {
                x,
                y: at + 10.,
                w: tw,
                h: 16.8,
            };
            if !labels.iter().any(|b| rects_overlap(*b, bounds, 4.)) {
                scene.text(
                    text,
                    P::new(x, at + 18.4),
                    12.,
                    label_color(color),
                    0.,
                    measure,
                );
                scene.regions.push(Region {
                    hit: Hit::Feature(FeatureHit { index, segment }),
                    shape: HitShape::Box(bounds),
                });
                labels.push(bounds);
                feature_labels += 1;
            }
        }
    }
    let mut site_labels = 0;
    let mut site_boxes = vec![];
    let label_rows = ((y - rect.y - 88.) / 23.).floor().max(1.) as usize;
    for (ordinal, (index, s)) in sites.iter().enumerate() {
        let position = if doc.circular {
            if s.position.is_multiple_of(n) {
                n
            } else {
                s.position % n
            }
        } else {
            s.position.min(n)
        };
        let x = left + width * position as f32 / n as f32;
        let anchor = P::new(x, y - 7.);
        scene.line(P::new(x, y + 3.), anchor, 1.2, rgb(0xa9b5c5).into());
        scene.regions.push(Region {
            hit: Hit::Site(*index),
            shape: HitShape::Box(Rect {
                x: x - 3.,
                y: y - 9.,
                w: 6.,
                h: 10.,
            }),
        });
        if options.labels && ordinal < MAX_SITE_LABELS {
            let text = format!(
                "{} ({})",
                truncate(&s.label, 20),
                cut_label(s, n, doc.circular)
            );
            let tw = measure(&text, 12.);
            let lx =
                (x - tw * 0.5).clamp(rect.x + 8., (rect.x + rect.w - tw - 8.).max(rect.x + 8.));
            for row in 0..label_rows {
                let ty = y - 42. - row as f32 * 23.;
                let bounds = Rect {
                    x: lx,
                    y: ty - 8.4,
                    w: tw,
                    h: 16.8,
                };
                if !site_boxes.iter().any(|b| rects_overlap(*b, bounds, 8.)) {
                    let end = P::new(lx + tw * 0.5, ty + 10.);
                    let hit = Hit::Site(*index);
                    scene.line(anchor, end, 1., rgb(0x63748d).into());
                    scene.regions.push(Region {
                        hit,
                        shape: HitShape::Line(anchor, end),
                    });
                    scene.connectors.push((hit, anchor, end));
                    scene.text(text, P::new(lx, ty), 12., rgb(INK).into(), 0., measure);
                    scene.regions.push(Region {
                        hit,
                        shape: HitShape::Box(bounds),
                    });
                    site_boxes.push(bounds);
                    site_labels += 1;
                    break;
                }
            }
        }
    }
    add_notices(
        &mut scene,
        data,
        drawn_features,
        drawn_primers,
        feature_labels,
        expected_feature_labels,
        primer_labels,
        sites.len(),
        site_labels,
    );
    scene
}
#[allow(clippy::too_many_arguments)]
fn add_notices(
    scene: &mut Scene,
    data: &MapData<'_>,
    features: usize,
    primers: usize,
    feature_labels: usize,
    expected_feature_labels: usize,
    primer_labels: usize,
    sites: usize,
    site_labels: usize,
) {
    let available_features = data
        .doc
        .features
        .iter()
        .filter(|f| !f.parts.is_empty() && f.kind != "primer_bind")
        .count();
    let available_primers = data
        .doc
        .features
        .iter()
        .filter(|f| !f.parts.is_empty() && f.kind == "primer_bind")
        .count();
    if data.options.features && features < available_features {
        scene
            .notices
            .push(format!("Feature tracks {features}/{available_features}"));
    }
    if data.options.primers && primers < available_primers {
        scene
            .notices
            .push(format!("Primer tracks {primers}/{available_primers}"));
    }
    if data.options.sites && sites < data.sites.len() {
        scene
            .notices
            .push(format!("Enzyme marks {sites}/{}", data.sites.len()));
    }
    if data.options.labels {
        if data.options.features && feature_labels < expected_feature_labels {
            scene.notices.push(format!(
                "Feature labels {feature_labels}/{expected_feature_labels}"
            ));
        }
        if data.options.primers && primer_labels < primers {
            scene
                .notices
                .push(format!("Primer labels {primer_labels}/{primers}"));
        }
        if data.options.sites && site_labels < data.sites.len() {
            scene
                .notices
                .push(format!("Enzyme labels {site_labels}/{}", data.sites.len()));
        }
    }
}
fn build_scene(rect: Rect, data: &MapData<'_>, measure: impl Fn(&str, f32) -> f32) -> Scene {
    if data.doc.sequence.is_empty() {
        return Scene::default();
    }
    let mut scene = if data.doc.circular && !data.options.linear {
        circular_scene(rect, data, &measure)
    } else {
        linear_scene(rect, data, &measure)
    };
    for region in scene.regions.clone() {
        let selected = !data.selected.is_empty()
            && match region.hit {
                Hit::Feature(hit) => data.doc.features[hit.index].parts == data.selected,
                Hit::Site(index) => {
                    site_parts(
                        &data.sites[index],
                        data.doc.sequence.len(),
                        data.doc.circular,
                    ) == data.selected
                }
            };
        if selected {
            scene.highlight(&region.shape, rgb(ACCENT).into());
        }
        if data.hovered.is_some_and(|hit| hit.same_object(region.hit)) {
            scene.highlight(&region.shape, rgb(INK).into());
        }
    }
    for (hit, a, b) in scene.connectors.clone() {
        if data.hovered.is_some_and(|hovered| hovered.same_object(hit)) {
            scene.line(a, b, 1.5, rgb(INK).into());
        }
    }
    // Explicitly disclose all geometry and label limits. Tables keep every item.
    if !scene.notices.is_empty() {
        let text = format!("{} · Zoom in or use the tables", scene.notices.join(" · "));
        let max_chars = ((rect.w - 32.) / 5.7).max(12.) as usize;
        for (row, chunk) in text
            .chars()
            .collect::<Vec<_>>()
            .chunks(max_chars)
            .take(3)
            .enumerate()
        {
            scene.text(
                chunk.iter().collect(),
                P::new(rect.x + 16., rect.y + rect.h - 16. - (2 - row) as f32 * 13.),
                10.,
                rgb(MUTED).into(),
                0.,
                &measure,
            );
        }
    }
    scene
}

pub(crate) enum MapEvent {
    Select {
        parts: Vec<Segment>,
        feature: u64,
        show_sequence: bool,
    },
    SelectRestriction {
        parts: Vec<Segment>,
        enzyme_index: usize,
    },
}
pub(crate) struct GeneMap {
    document: Rc<Document>,
    selected: Vec<Segment>,
    options: MapOptions,
    sites: Rc<Vec<MapRestrictionSite>>,
    scene: Rc<RefCell<(Rect, Vec<Region>)>>,
    viewport: Rc<RefCell<Rect>>,
    scroll: ScrollHandle,
    hovered: Option<Hit>,
    hover_generation: u64,
    hover_ready: bool,
    pointer: P,
}
impl EventEmitter<MapEvent> for GeneMap {}
impl GeneMap {
    pub(crate) fn new(document: Document, _window: &mut Window, _cx: &mut Context<Self>) -> Self {
        let rect = Rect {
            x: 0.,
            y: 0.,
            w: 0.,
            h: 0.,
        };
        Self {
            document: Rc::new(document),
            selected: vec![],
            options: MapOptions::default(),
            sites: Rc::new(vec![]),
            scene: Rc::new(RefCell::new((rect, vec![]))),
            viewport: Rc::new(RefCell::new(rect)),
            scroll: ScrollHandle::new(),
            hovered: None,
            hover_generation: 0,
            hover_ready: false,
            pointer: P::new(0., 0.),
        }
    }
    pub(crate) fn set_document(
        &mut self,
        document: Document,
        selected: Vec<Segment>,
        mut options: MapOptions,
        restriction_sites: Vec<MapRestrictionSite>,
        cx: &mut Context<Self>,
    ) {
        options.zoom = if options.zoom.is_finite() {
            options.zoom.clamp(1., 3.)
        } else {
            1.
        };
        if *self.document == document
            && self.selected == selected
            && self.options == options
            && *self.sites == restriction_sites
        {
            return;
        }
        if *self.document != document || self.options != options || *self.sites != restriction_sites
        {
            self.hovered = None;
            self.hover_ready = false;
            self.hover_generation += 1;
            self.scene.borrow_mut().1.clear();
        }
        if options.zoom == 1. {
            self.scroll.set_offset(point(px(0.), px(0.)));
        }
        self.document = Rc::new(document);
        self.selected = selected;
        self.options = options;
        self.sites = Rc::new(restriction_sites);
        cx.notify();
    }
    fn clear_hover(&mut self, cx: &mut Context<Self>) {
        if self.hovered.take().is_some() {
            self.hover_ready = false;
            self.hover_generation += 1;
            cx.notify();
        }
    }
}
fn restriction_card(doc: &Document, site: &MapRestrictionSite) -> AnyElement {
    let parts = site_parts(site, doc.sequence.len(), doc.circular);
    let location = parts
        .iter()
        .map(|p| format!("{}–{}", number(p.start + 1), number(p.end)))
        .collect::<Vec<_>>()
        .join(", ");
    let sequence: String = parts
        .iter()
        .flat_map(|p| doc.sequence[p.start..p.end].bytes())
        .take(60)
        .map(char::from)
        .collect();
    let cut = if !doc.circular && site.position == 0 {
        "Cut before base 1".into()
    } else {
        format!(
            "Cut after base {}",
            cut_label(site, doc.sequence.len(), doc.circular)
        )
    };
    let enzyme = popgene_core::restriction::ENZYMES.get(site.enzyme_index);
    div()
        .w(px(342.))
        .p_4()
        .flex()
        .flex_col()
        .gap_2()
        .bg(rgb(0x202a3a))
        .border_1()
        .border_color(rgb(0x3c4b61))
        .rounded_lg()
        .shadow_lg()
        .text_color(rgb(INK))
        .text_size(px(12.))
        .child(
            div()
                .font_weight(FontWeight::SEMIBOLD)
                .text_size(px(16.))
                .child(truncate(&site.label, 120)),
        )
        .when_some(enzyme, |card, enzyme| {
            card.child(format!("Recognition motif: 5′–{}–3′", enzyme.recognition))
                .child(div().text_color(rgb(MUTED)).child(format!(
                    "Motif cut offsets: {}/{}",
                    enzyme.top_cut, enzyme.bottom_cut
                )))
                .when(enzyme.name == "NmeAIII", |card| {
                    card.child(
                        div().text_color(rgb(MUTED)).child(
                            "Nominal cut: NmeAIII can shift by one base with sequence context.",
                        ),
                    )
                })
        })
        .child(cut)
        .child(format!("Recognition: {location} · 1-based, inclusive"))
        .child(if site.reverse {
            "Reverse-strand match"
        } else {
            "Forward-strand match"
        })
        .child(
            div()
                .text_color(rgb(ACCENT))
                .child(format!("Reference bases: {sequence}")),
        )
        .child(
            div()
                .text_color(rgb(MUTED))
                .child("Click to select recognition bases · Double-click to show sequence"),
        )
        .into_any_element()
}

impl Render for GeneMap {
    fn render(&mut self, _window: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        let document = self.document.clone();
        let selected = self.selected.clone();
        let options = self.options;
        let sites = self.sites.clone();
        let hovered = self.hovered;
        let geometry = self.scene.clone();
        let viewport = self.viewport.clone();
        let view = cx.entity().downgrade();
        let surface = canvas(
            move |bounds, window, cx| {
                let rect = Rect {
                    x: bounds.origin.x.as_f32(),
                    y: bounds.origin.y.as_f32(),
                    w: bounds.size.width.as_f32(),
                    h: bounds.size.height.as_f32(),
                };
                let font = window.text_style().font();
                let text_system = window.text_system();
                let scene = build_scene(
                    rect,
                    &MapData {
                        doc: &document,
                        selected: &selected,
                        options,
                        sites: &sites,
                        hovered,
                    },
                    |text, size| {
                        text_system
                            .shape_line(
                                text.to_owned().into(),
                                px(size),
                                &[TextRun {
                                    len: text.len(),
                                    font: font.clone(),
                                    color: rgb(INK).into(),
                                    background_color: None,
                                    underline: None,
                                    strikethrough: None,
                                }],
                                None,
                            )
                            .width
                            .as_f32()
                    },
                );
                let pointer = window.mouse_position();
                if hovered.is_some()
                    && hit_test(
                        &scene.regions,
                        P::new(pointer.x.as_f32(), pointer.y.as_f32()),
                    ) != hovered
                {
                    cx.defer(move |cx| {
                        let _ = view.update(cx, |view, cx| view.clear_hover(cx));
                    });
                }
                *geometry.borrow_mut() = (rect, scene.regions.clone());
                scene
            },
            move |_, scene, window, cx| {
                for draw in scene.draws {
                    match draw {
                        Draw::Stroke(points, width, color) => {
                            let mut path = PathBuilder::stroke(px(width));
                            for (index, p) in points.iter().enumerate() {
                                if index == 0 {
                                    path.move_to(p.gpui());
                                } else {
                                    path.line_to(p.gpui());
                                }
                            }
                            if let Ok(path) = path.build() {
                                window.paint_path(path, color);
                            }
                        }
                        Draw::Fill(points, color) => {
                            let mut path = PathBuilder::fill();
                            for (index, p) in points.iter().enumerate() {
                                if index == 0 {
                                    path.move_to(p.gpui());
                                } else {
                                    path.line_to(p.gpui());
                                }
                            }
                            path.close();
                            if let Ok(path) = path.build() {
                                window.paint_path(path, color);
                            }
                        }
                        Draw::Text {
                            point,
                            text,
                            size,
                            color,
                        } => {
                            let font = window.text_style().font();
                            let line = window.text_system().shape_line(
                                text.clone().into(),
                                px(size),
                                &[TextRun {
                                    len: text.len(),
                                    font,
                                    color,
                                    background_color: None,
                                    underline: None,
                                    strikethrough: None,
                                }],
                                None,
                            );
                            let _ = line.paint(
                                point.gpui(),
                                px(size * 1.4),
                                TextAlign::Left,
                                None,
                                window,
                                cx,
                            );
                        }
                    }
                }
            },
        )
        .size_full();
        let mut root = div()
            .id("gene-map")
            .relative()
            .size_full()
            .min_h(px(260.))
            .bg(rgb(0x111a27))
            .when(hovered.is_some(), |view| view.cursor_pointer())
            .on_hover(cx.listener(|this, inside, _, cx| {
                if !inside {
                    this.clear_hover(cx);
                }
            }))
            .on_mouse_move(cx.listener(|this, event: &MouseMoveEvent, window, cx| {
                let pointer = P::new(event.position.x.as_f32(), event.position.y.as_f32());
                let hit = hit_test(&this.scene.borrow().1, pointer);
                if hit != this.hovered {
                    this.hovered = hit;
                    this.hover_ready = false;
                    this.hover_generation += 1;
                    this.pointer = pointer;
                    cx.notify();
                    if hit.is_some() {
                        let generation = this.hover_generation;
                        cx.spawn_in(window, async move |view, cx| {
                            cx.background_executor()
                                .timer(Duration::from_millis(300))
                                .await;
                            let _ = view.update(cx, |view, cx| {
                                if view.hover_generation == generation && view.hovered.is_some() {
                                    view.hover_ready = true;
                                    cx.notify();
                                }
                            });
                        })
                        .detach();
                    }
                }
            }))
            .on_mouse_down(
                MouseButton::Left,
                cx.listener(|this, event: &MouseDownEvent, _, cx| {
                    let pointer = P::new(event.position.x.as_f32(), event.position.y.as_f32());
                    let hit = hit_test(&this.scene.borrow().1, pointer);
                    if let Some(hit) = hit {
                        let (parts, event) = match hit {
                            Hit::Feature(hit) => {
                                let f = &this.document.features[hit.index];
                                let parts = f.parts.clone();
                                (
                                    parts.clone(),
                                    MapEvent::Select {
                                        parts,
                                        feature: f.id,
                                        show_sequence: event.click_count >= 2,
                                    },
                                )
                            }
                            Hit::Site(index) => {
                                let site = &this.sites[index];
                                let parts = site_parts(
                                    site,
                                    this.document.sequence.len(),
                                    this.document.circular,
                                );
                                (
                                    parts.clone(),
                                    MapEvent::SelectRestriction {
                                        parts,
                                        enzyme_index: site.enzyme_index,
                                    },
                                )
                            }
                        };
                        this.selected = parts;
                        this.hover_ready = false;
                        this.hover_generation += 1;
                        cx.emit(event);
                        cx.notify();
                    }
                }),
            )
            .child(
                canvas(
                    move |bounds, _, _| {
                        *viewport.borrow_mut() = Rect {
                            x: bounds.origin.x.as_f32(),
                            y: bounds.origin.y.as_f32(),
                            w: bounds.size.width.as_f32(),
                            h: bounds.size.height.as_f32(),
                        };
                    },
                    |_, _, _, _| {},
                )
                .absolute()
                .size_full(),
            )
            .child(
                div()
                    .id("map-scroll")
                    .size_full()
                    .overflow_scroll()
                    .track_scroll(&self.scroll)
                    .child(
                        div()
                            .w(relative(options.zoom))
                            .h(relative(options.zoom))
                            .flex_shrink_0()
                            .child(surface),
                    ),
            );
        if self.document.sequence.is_empty() {
            root = root.child(
                div()
                    .absolute()
                    .inset_0()
                    .flex()
                    .flex_col()
                    .items_center()
                    .justify_center()
                    .gap_2()
                    .text_color(rgb(MUTED))
                    .child(
                        div()
                            .text_size(px(22.))
                            .child("A place for your next sequence"),
                    )
                    .child("Add sequence to start exploring this document."),
            );
        }
        if self.hover_ready
            && let Some(hit) = self.hovered
        {
            let bounds = *self.viewport.borrow();
            let x = (self.pointer.x - bounds.x + 18.).clamp(8., (bounds.w - 350.).max(8.));
            let y = (self.pointer.y - bounds.y + 18.).clamp(8., (bounds.h - 340.).max(8.));
            let card = match hit {
                Hit::Feature(hit) => render_feature_card(
                    &self.document,
                    &self.document.features[hit.index],
                    hit.segment,
                ),
                Hit::Site(index) => restriction_card(&self.document, &self.sites[index]),
            };
            root = root.child(
                div()
                    .id("map-feature-card")
                    .absolute()
                    .left(px(x))
                    .top(px(y))
                    .max_h(px((bounds.h - y - 8.).max(160.)))
                    .overflow_y_scroll()
                    .child(card),
            );
        }
        root
    }
}

#[cfg(test)]
mod tests {
    use super::{
        Draw, FeatureHit, Hit, HitShape, MapData, MapOptions, MapRestrictionSite, P, Rect, Region,
        assign_tracks, build_scene, cut_label, hit_test, major_ticks, on_circle, site_parts,
    };
    use popgene_core::Segment;
    fn sample() -> popgene_core::Document {
        popgene_formats::import(
            "annotated.gb",
            include_bytes!("../../../validation/fixtures/annotated.gb"),
        )
        .unwrap()
    }
    fn measure(s: &str, size: f32) -> f32 {
        s.chars().count() as f32 * size * 0.55
    }
    fn rect() -> Rect {
        Rect {
            x: 0.,
            y: 0.,
            w: 950.,
            h: 750.,
        }
    }
    fn site() -> MapRestrictionSite {
        MapRestrictionSite {
            label: "EcoRI".into(),
            position: 2345,
            start: 2344,
            length: 6,
            reverse: false,
            enzyme_index: 7,
        }
    }
    #[test]
    fn circular_hit_testing_handles_origin_gaps_and_tiny_features() {
        let center = P::new(300., 300.);
        let regions = vec![
            Region {
                hit: Hit::Feature(FeatureHit {
                    index: 0,
                    segment: Some(0),
                }),
                shape: HitShape::Arc {
                    center,
                    radius: 100.,
                    start: 0.98,
                    end: 1.,
                },
            },
            Region {
                hit: Hit::Feature(FeatureHit {
                    index: 0,
                    segment: Some(1),
                }),
                shape: HitShape::Arc {
                    center,
                    radius: 100.,
                    start: 0.,
                    end: 0.024,
                },
            },
            Region {
                hit: Hit::Feature(FeatureHit {
                    index: 1,
                    segment: Some(0),
                }),
                shape: HitShape::Arc {
                    center,
                    radius: 100.,
                    start: 0.4,
                    end: 0.40001,
                },
            },
        ];
        assert_eq!(
            hit_test(&regions, on_circle(center, 100., 0.99)),
            Some(Hit::Feature(FeatureHit {
                index: 0,
                segment: Some(0)
            }))
        );
        assert_eq!(
            hit_test(&regions, on_circle(center, 100., 0.02)),
            Some(Hit::Feature(FeatureHit {
                index: 0,
                segment: Some(1)
            }))
        );
        assert_eq!(
            hit_test(&regions, on_circle(center, 100., 0.40001).offset(2., 0.)),
            Some(Hit::Feature(FeatureHit {
                index: 1,
                segment: Some(0)
            }))
        );
        assert_eq!(hit_test(&regions, on_circle(center, 100., 0.5)), None);
        assert_eq!(hit_test(&regions, center), None);
    }
    #[test]
    fn labels_win_over_crossing_leaders() {
        let label = Hit::Site(0);
        let regions = vec![
            Region {
                hit: label,
                shape: HitShape::Box(Rect {
                    x: 0.,
                    y: 0.,
                    w: 30.,
                    h: 17.,
                }),
            },
            Region {
                hit: Hit::Site(1),
                shape: HitShape::Line(P::new(0., 8.), P::new(50., 8.)),
            },
        ];
        assert_eq!(hit_test(&regions, P::new(10., 8.)), Some(label));
        assert_eq!(hit_test(&regions, P::new(40., 8.)), Some(Hit::Site(1)));
    }
    #[test]
    fn natural_ticks_and_restriction_boundaries_are_not_residue_indices() {
        assert_eq!(major_ticks(5000), vec![1000, 2000, 3000, 4000, 5000]);
        assert_eq!(major_ticks(1), vec![1]);
        assert!(major_ticks(0).is_empty());
        assert_eq!(cut_label(&site(), 5000, true), "2,345");
        let mut origin = site();
        origin.position = 0;
        origin.start = 4997;
        assert_eq!(cut_label(&origin, 5000, true), "5,000");
        assert_eq!(cut_label(&origin, 5000, false), "before 1");
        assert_eq!(
            site_parts(&origin, 5000, true),
            vec![
                Segment {
                    start: 4997,
                    end: 5000
                },
                Segment { start: 0, end: 3 }
            ]
        );
        assert!(site_parts(&origin, 5000, false).is_empty());
        origin.length = usize::MAX;
        assert!(site_parts(&origin, 5000, true).is_empty());
    }
    #[test]
    fn overlapping_compound_features_keep_consistent_separate_tracks() {
        let doc = sample();
        let tracks = assign_tracks(&doc, &[0, 1, 2, 3]);
        assert_ne!(tracks[0], tracks[1]);
        assert_ne!(tracks[1], tracks[2]);
        assert_eq!(tracks[0], tracks[2]);
        assert_eq!(tracks[0], tracks[3]);
    }
    #[test]
    fn display_subdivisions_have_distinct_map_colors_hits_and_terminal_arrow() {
        let doc = popgene_formats::import(
            "segments.gb",
            include_bytes!("../../../validation/fixtures/regression/segment-colours.gb"),
        )
        .unwrap();
        let original = doc.clone();
        for linear in [false, true] {
            let scene = build_scene(
                rect(),
                &MapData {
                    doc: &doc,
                    selected: &[],
                    options: MapOptions {
                        linear,
                        ..Default::default()
                    },
                    sites: &[],
                    hovered: None,
                },
                measure,
            );
            for value in [0xff0000, 0x0000ff] {
                assert!(scene.draws.iter().any(|draw| matches!(draw, Draw::Stroke(_, width, color) if *width == 10. && *color == super::rgb(value).into())));
            }
            assert!(scene.draws.iter().any(|draw| matches!(draw, Draw::Fill(points, color) if points.len() == 3 && *color == super::rgb(0x0000ff).into())));
            assert!(scene.regions.iter().any(|region| matches!(
                region.hit,
                Hit::Feature(FeatureHit {
                    index: 0,
                    segment: Some(1)
                })
            )));
            for (segment, name) in ["Red domain", "Blue domain"].iter().enumerate() {
                assert!(scene.draws.iter().any(|draw| matches!(draw, Draw::Text { text, color, .. } if text == name && color.l >= 0.72)), "{name} missing from map mode linear={linear}");
                let hit = Hit::Feature(FeatureHit {
                    index: 0,
                    segment: Some(segment),
                });
                let bounds = scene
                    .regions
                    .iter()
                    .find_map(|region| match &region.shape {
                        HitShape::Box(bounds) if region.hit == hit && bounds.h > 12. => {
                            Some(*bounds)
                        }
                        _ => None,
                    })
                    .expect("each segment label has its own feature-card target");
                assert_eq!(hit_test(&scene.regions, bounds.center()), Some(hit));
            }
            assert!(
                !scene
                    .notices
                    .iter()
                    .any(|notice| notice.starts_with("Feature labels"))
            );
        }
        assert_eq!(doc, original);
        assert_eq!(doc.features[0].parts.len(), 1);
    }
    #[test]
    fn repeated_locations_use_the_terminal_occurrences_arrow_color() {
        let mut doc = popgene_core::Document::new("repeated", "ACGTAC").unwrap();
        doc.circular = true;
        doc.features.push(popgene_core::Feature {
            id: 1, label: "repeated".into(), kind: "misc_feature".into(),
            parts: vec![Segment { start: 0, end: 6 }, Segment { start: 0, end: 6 }],
            reverse: false, edited: false,
            qualifiers: vec![("note".into(), Some("This bidirectional feature has 2 segments: 1: 1 .. 6 / #ff0000 / First occurrence 2: 1 .. 6 / #0000ff / Second occurrence".into()))],
        });
        for linear in [false, true] {
            let scene = build_scene(
                rect(),
                &MapData {
                    doc: &doc,
                    selected: &[],
                    options: MapOptions {
                        linear,
                        ..Default::default()
                    },
                    sites: &[],
                    hovered: None,
                },
                measure,
            );
            let arrows: Vec<_> = scene
                .draws
                .iter()
                .filter_map(|draw| match draw {
                    Draw::Fill(points, color) if points.len() == 3 => Some(*color),
                    _ => None,
                })
                .collect();
            assert_eq!(
                arrows,
                [super::rgb(0x0000ff).into(), super::rgb(0xff0000).into()]
            );
        }
    }

    #[test]
    fn map_modes_preserve_parts_and_read_only_document_with_compact_labels() {
        let doc = sample();
        let before = doc.clone();
        let sites = [site()];
        for linear in [false, true] {
            let options = MapOptions {
                linear,
                ..Default::default()
            };
            let r = Rect {
                w: 530.,
                h: 550.,
                ..rect()
            };
            let scene = build_scene(
                r,
                &MapData {
                    doc: &doc,
                    selected: &[],
                    options,
                    sites: &sites,
                    hovered: None,
                },
                measure,
            );
            assert_eq!(
                scene
                    .regions
                    .iter()
                    .filter(|r| matches!(
                        r.hit,
                        Hit::Feature(FeatureHit {
                            index: 0,
                            segment: Some(_)
                        })
                    ) && !matches!(r.shape, HitShape::Line(..)))
                    .count(),
                2
            );
            for r in scene.regions.iter().filter(|r| {
                matches!(
                    r.hit,
                    Hit::Feature(FeatureHit { segment: None, .. }) | Hit::Site(_)
                )
            }) {
                if let HitShape::Box(bounds) = r.shape {
                    assert!(
                        bounds.x >= 0. && bounds.x + bounds.w <= 530.1,
                        "label out of canvas: {bounds:?}"
                    );
                }
            }
        }
        assert_eq!(doc, before);
    }
    #[test]
    fn display_options_actually_remove_layers_without_mutating_selection() {
        let doc = sample();
        let sites = [site()];
        let selected = doc.features[0].parts.clone();
        for linear in [false, true] {
            let options = MapOptions {
                linear,
                features: false,
                primers: false,
                sites: false,
                positions: false,
                labels: false,
                ..Default::default()
            };
            let scene = build_scene(
                rect(),
                &MapData {
                    doc: &doc,
                    selected: &selected,
                    options,
                    sites: &sites,
                    hovered: None,
                },
                measure,
            );
            assert!(scene.regions.is_empty());
            assert!(!scene.draws.iter().any(|draw|matches!(draw,Draw::Text{text,..} if text=="1,000"||text.contains("EcoRI")||text.contains("audit_forward"))));
            assert_eq!(selected, doc.features[0].parts);
        }
    }
    #[test]
    fn restriction_labels_and_leaders_are_interactive_and_hover_all_geometry() {
        let doc = sample();
        let sites = [site()];
        let options = MapOptions::default();
        let plain = build_scene(
            rect(),
            &MapData {
                doc: &doc,
                selected: &[],
                options,
                sites: &sites,
                hovered: None,
            },
            measure,
        );
        let hovered = build_scene(
            rect(),
            &MapData {
                doc: &doc,
                selected: &[],
                options,
                sites: &sites,
                hovered: Some(Hit::Site(0)),
            },
            measure,
        );
        assert!(hovered.draws.len() > plain.draws.len());
        let label = plain
            .regions
            .iter()
            .rev()
            .find(|r| r.hit == Hit::Site(0) && matches!(r.shape,HitShape::Box(Rect{h,..}) if h>10.))
            .unwrap();
        let HitShape::Box(bounds) = label.shape else {
            unreachable!()
        };
        assert_eq!(
            hit_test(&plain.regions, bounds.center()),
            Some(Hit::Site(0))
        );
        let leader = plain
            .regions
            .iter()
            .find(|r| r.hit == Hit::Site(0) && matches!(r.shape, HitShape::Line(..)))
            .unwrap();
        let HitShape::Line(a, b) = leader.shape else {
            unreachable!()
        };
        assert_eq!(
            hit_test(&plain.regions, P::new((a.x + b.x) * 0.5, (a.y + b.y) * 0.5)),
            Some(Hit::Site(0))
        );
    }
    #[test]
    fn dense_labels_report_omissions_and_zoom_increases_available_tracks() {
        let mut doc = sample();
        let template = doc.features[0].clone();
        doc.features.clear();
        for index in 0..24 {
            let mut f = template.clone();
            f.id = index;
            f.label = format!("overlap {index}");
            doc.features.push(f);
        }
        let sites: Vec<_> = (0..100)
            .map(|i| MapRestrictionSite {
                label: format!("Enzyme{i}"),
                start: i * 40,
                position: i * 40 + 1,
                ..site()
            })
            .collect();
        let data = MapData {
            doc: &doc,
            selected: &[],
            options: MapOptions::default(),
            sites: &sites,
            hovered: None,
        };
        let small = build_scene(rect(), &data, measure);
        let large = build_scene(
            Rect {
                w: 1900.,
                h: 1500.,
                ..rect()
            },
            &data,
            measure,
        );
        assert!(small.notices.iter().any(|n| n.starts_with("Enzyme labels")));
        assert!(
            small
                .notices
                .iter()
                .any(|n| n.starts_with("Feature tracks"))
        );
        let glyphs = |scene: &super::Scene| {
            scene
                .regions
                .iter()
                .filter(|r| {
                    matches!(
                        r.hit,
                        Hit::Feature(FeatureHit {
                            segment: Some(_),
                            ..
                        })
                    )
                })
                .count()
        };
        assert!(glyphs(&large) > glyphs(&small));
    }
}
