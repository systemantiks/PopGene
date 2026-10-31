//! Application menus and shortcuts. A command is enabled only while the workspace
//! registers a handler for it, so AppKit greys out commands that do not apply.
use crate::{
    app::*,
    model::{EnzymeFilter, View},
};
use gpui_kit::component::{GlobalState, input};
use gpui_kit::*;

/// Checkmarks and labels that follow the active document.
#[derive(Clone, Copy, Default, PartialEq)]
pub struct MenuState {
    pub view: Option<View>,
    pub features: bool,
    pub primers: bool,
    pub sites: bool,
    pub circular: bool,
    pub navigator: bool,
    /// The site-count filter and whether a custom enzyme list replaces the catalogue.
    pub enzymes: Option<(EnzymeFilter, bool)>,
}

/// Command on macOS, Control elsewhere.
fn key(keys: &str, action: impl Action) -> KeyBinding {
    let cmd = if cfg!(target_os = "macos") {
        "cmd"
    } else {
        "ctrl"
    };
    KeyBinding::new(&format!("{cmd}-{keys}"), action, Some("PopGene"))
}

pub fn bind_keys(cx: &mut App) {
    cx.bind_keys([
        key("n", NewSequence),
        key("o", OpenFiles),
        key("shift-o", OpenFolder),
        key("s", Save),
        key("shift-s", SaveAs),
        key("e", Export),
        key("p", Print),
        key("f", FindSequence),
        key("g", FindNext),
        key("shift-g", FindPrevious),
        key("l", SelectRange),
        key("w", CloseTab),
        key("q", Quit),
        key("z", Undo),
        key("shift-z", Redo),
    ]);
    // Views follow tab order.
    cx.bind_keys([
        key("1", ViewMap),
        key("2", ViewSequence),
        key("3", ViewEnzymes),
        key("4", ViewFeatures),
        key("5", ViewPrimers),
        key("6", ViewTrace),
        key("7", ViewParameters),
        key("8", ViewProperties),
        key("9", ViewHistory),
        key("0", ToggleNavigator),
        key("t", AddFeature),
        key("r", AddPrimer),
    ]);
    cx.bind_keys([
        KeyBinding::new("ctrl-tab", NextTab, Some("PopGene")),
        KeyBinding::new("ctrl-shift-tab", PreviousTab, Some("PopGene")),
        KeyBinding::new("escape", Escape, Some("PopGene")),
    ]);
}

/// Replaces the menu bar. Windows and Linux draw it from the shared menu state.
pub fn install(state: MenuState, cx: &mut App) {
    cx.set_menus(menus(state));
    if !cfg!(target_os = "macos")
        && let Some(menus) = cx.get_menus()
    {
        GlobalState::global_mut(cx).set_app_menus(menus);
    }
}

fn menus(state: MenuState) -> Vec<Menu> {
    let view = |view| state.view == Some(view);
    let filter = |filter| state.enzymes.is_some_and(|(f, _)| f == filter);
    vec![
        Menu::new("PopGene").items([
            MenuItem::action("About PopGene", About),
            MenuItem::action("Check for Updates…", CheckForUpdates),
            MenuItem::separator(),
            MenuItem::action("Quit PopGene", Quit),
        ]),
        Menu::new("File").items([
            MenuItem::action("New Sequence", NewSequence),
            MenuItem::action("Open…", OpenFiles),
            MenuItem::action("Open Folder…", OpenFolder),
            MenuItem::separator(),
            MenuItem::action("Close Tab", CloseTab),
            MenuItem::action("Save", Save),
            MenuItem::action("Save As…", SaveAs),
            MenuItem::separator(),
            MenuItem::action("Export Sequence…", Export),
            MenuItem::action("Export Figure…", ExportFigure),
            MenuItem::separator(),
            MenuItem::action("Print…", Print),
        ]),
        Menu::new("Edit").items([
            MenuItem::action("Undo", Undo),
            MenuItem::action("Redo", Redo),
            MenuItem::separator(),
            // OS selectors keep standard editing shortcuts available
            // in native file panels as well as PopGene text controls.
            MenuItem::os_action("Cut", input::Cut, OsAction::Cut),
            MenuItem::os_action("Copy", input::Copy, OsAction::Copy),
            MenuItem::os_action("Paste", input::Paste, OsAction::Paste),
            MenuItem::os_action("Select All", input::SelectAll, OsAction::SelectAll),
            MenuItem::separator(),
            MenuItem::action("Find Sequence…", FindSequence),
            MenuItem::action("Find Next", FindNext),
            MenuItem::action("Find Previous", FindPrevious),
            MenuItem::separator(),
            MenuItem::action("Select Range…", SelectRange),
        ]),
        Menu::new("View").items([
            MenuItem::action("Map", ViewMap).checked(view(View::Map)),
            MenuItem::action("Sequence", ViewSequence).checked(view(View::Sequence)),
            MenuItem::action("Map and Sequence", ViewSplit).checked(view(View::Split)),
            MenuItem::action("Enzymes", ViewEnzymes).checked(view(View::Enzymes)),
            MenuItem::action("Features", ViewFeatures).checked(view(View::Features)),
            MenuItem::action("Primers", ViewPrimers).checked(view(View::Primers)),
            MenuItem::action("Trace", ViewTrace).checked(view(View::Trace)),
            MenuItem::action("Parameters", ViewParameters).checked(view(View::Parameters)),
            MenuItem::action("Properties", ViewProperties).checked(view(View::Properties)),
            MenuItem::action("History", ViewHistory).checked(view(View::History)),
            MenuItem::separator(),
            MenuItem::action("Project Navigator", ToggleNavigator).checked(state.navigator),
            MenuItem::separator(),
            MenuItem::action("Next Tab", NextTab),
            MenuItem::action("Previous Tab", PreviousTab),
        ]),
        Menu::new("Actions").items([
            MenuItem::action("Edit Sequence…", EditSequence),
            MenuItem::action("Reverse Complement", ReverseComplement),
            MenuItem::action(
                if state.circular {
                    "Make Linear"
                } else {
                    "Make Circular"
                },
                ToggleTopology,
            ),
            MenuItem::separator(),
            MenuItem::action("PCR & Cloning…", CloningTools),
            MenuItem::action("Restriction Digestion…", DigestTools),
            MenuItem::action("Align Open Sequences…", AlignmentTools),
            MenuItem::separator(),
            MenuItem::action("Full Duplex Strands & Annotations…", DuplexView),
            MenuItem::action("Edit / Trim Chromatogram…", TraceTools),
            MenuItem::action("Protein Analysis…", ProteinTools),
            MenuItem::separator(),
            MenuItem::action("References, Attachments & History…", MetadataTools),
        ]),
        Menu::new("Features").items([
            MenuItem::action("Add Feature…", AddFeature),
            MenuItem::separator(),
            MenuItem::action("Show Features", ToggleFeatureLayer).checked(state.features),
        ]),
        Menu::new("Primers").items([
            MenuItem::action("Add or Design Primers…", AddPrimer),
            MenuItem::action("Edit Selected Primer…", EditPrimer),
            MenuItem::separator(),
            MenuItem::action("Show Primers", TogglePrimerLayer).checked(state.primers),
        ]),
        Menu::new("Enzymes").items([
            MenuItem::action("Use All Catalogue Enzymes", UseAllEnzymes)
                .checked(state.enzymes.is_some_and(|(_, custom)| !custom)),
            MenuItem::action("Choose Enzymes…", ChooseEnzymes)
                .checked(state.enzymes.is_some_and(|(_, custom)| custom)),
            MenuItem::separator(),
            MenuItem::action(EnzymeFilter::UniqueSix.label(), FilterUniqueSix)
                .checked(filter(EnzymeFilter::UniqueSix)),
            MenuItem::action(EnzymeFilter::Unique.label(), FilterUnique)
                .checked(filter(EnzymeFilter::Unique)),
            MenuItem::action(EnzymeFilter::OneOrTwo.label(), FilterOneOrTwo)
                .checked(filter(EnzymeFilter::OneOrTwo)),
            MenuItem::action(EnzymeFilter::AllSix.label(), FilterAllSix)
                .checked(filter(EnzymeFilter::AllSix)),
            MenuItem::action(EnzymeFilter::All.label(), FilterAll)
                .checked(filter(EnzymeFilter::All)),
            MenuItem::separator(),
            MenuItem::action("Show Restriction Sites", ToggleSiteLayer).checked(state.sites),
        ]),
    ]
}
