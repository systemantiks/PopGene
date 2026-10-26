//! Embedded figure fonts keep exports independent of the host's font packages.
use std::sync::Arc;
use svg2pdf::usvg::fontdb::{Database, Source};

pub(super) const DEFAULT_FAMILY: &str = "Noto Sans";

pub(super) fn database(include_system: bool) -> Database {
    let mut db = Database::new();
    // The original font files and their OFL notices/provenance are retained in
    // assets/fonts and vendor-notices/figure-fonts. Embed bytes in every binary,
    // including Linux builds without a font installation or asset directory.
    for data in [
        include_bytes!("../../../../assets/fonts/NotoSans-Regular.ttf").as_slice(),
        include_bytes!("../../../../assets/fonts/FiraMono-Medium.ttf").as_slice(),
    ] {
        db.load_font_source(Source::Binary(Arc::new(data)));
    }
    if include_system {
        // Extra installed faces can supply scripts missing from the two defaults.
        db.load_system_fonts();
    }
    // Fontconfig may name absent preferred families. A nonempty database alone
    // does not make generic queries resolvable, so bind all generic fallbacks
    // after loading system aliases. usvg also appends Serif as its last resort.
    db.set_sans_serif_family(DEFAULT_FAMILY);
    db.set_serif_family(DEFAULT_FAMILY);
    db.set_cursive_family(DEFAULT_FAMILY);
    db.set_fantasy_family(DEFAULT_FAMILY);
    db.set_monospace_family("Fira Mono");
    db
}
