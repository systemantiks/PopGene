#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]
extern crate gpui_kit as gpui;
use gpui_kit::component::{Theme, ThemeMode, TitleBar, input};
use gpui_kit::*;
mod app;
mod enzyme_sidebar;
mod feature_card;
mod forms;
mod gene_map;
mod model;
mod open_dialog;
mod protein_properties;
mod sequence_layout;
mod sequence_view;
mod updater;
use app::*;

fn main() {
    let paths = std::env::args()
        .skip(1)
        .map(std::path::PathBuf::from)
        .collect();
    gpui_kit::application()
        .with_assets(WorkspaceAssets)
        .run(move |cx| {
            gpui_kit::init(cx);
            enzyme_sidebar::init(cx);
            Theme::change(ThemeMode::Dark, None, cx);
            Theme::update(cx, |theme| {
                theme.font_size = px(14.);
                theme.colors.background = rgb(BG).into();
                theme.colors.border = rgb(BORDER).into();
                theme.colors.foreground = rgb(INK).into();
                theme.colors.primary = rgb(ACCENT).into();
                theme.colors.primary_foreground = rgb(BG).into();
                theme.colors.button_primary = rgb(ACCENT).into();
                theme.colors.button_primary_foreground = rgb(BG).into();
                theme.colors.button_primary_hover = rgb(0x68d9cb).into();
                theme.colors.button_primary_active = rgb(0x30ad9f).into();
                theme.colors.input = rgb(BG).into();
            });
            let cmd = if cfg!(target_os = "macos") {
                "cmd"
            } else {
                "ctrl"
            };
            cx.bind_keys([
                KeyBinding::new(&format!("{cmd}-n"), NewSequence, Some("PopGene")),
                KeyBinding::new(&format!("{cmd}-o"), OpenFiles, Some("PopGene")),
                KeyBinding::new(&format!("{cmd}-shift-o"), OpenFolder, Some("PopGene")),
                KeyBinding::new(&format!("{cmd}-s"), Save, Some("PopGene")),
                KeyBinding::new(&format!("{cmd}-shift-s"), SaveAs, Some("PopGene")),
                KeyBinding::new(&format!("{cmd}-e"), Export, Some("PopGene")),
                KeyBinding::new(&format!("{cmd}-f"), FindSequence, Some("PopGene")),
                KeyBinding::new(&format!("{cmd}-w"), CloseTab, Some("PopGene")),
                KeyBinding::new(&format!("{cmd}-q"), Quit, Some("PopGene")),
                KeyBinding::new(&format!("{cmd}-z"), Undo, Some("PopGene")),
                KeyBinding::new(&format!("{cmd}-shift-z"), Redo, Some("PopGene")),
            ]);
            cx.bind_keys([
                KeyBinding::new("ctrl-tab", NextTab, Some("PopGene")),
                KeyBinding::new("ctrl-shift-tab", PreviousTab, Some("PopGene")),
                KeyBinding::new("escape", Escape, Some("PopGene")),
            ]);
            cx.set_menus([
                Menu::new("PopGene").items([
                    MenuItem::action("Check for Updates…", CheckForUpdates),
                    MenuItem::separator(),
                    MenuItem::action("Quit PopGene", Quit),
                ]),
                Menu::new("File").items([
                    MenuItem::action("New Sequence", NewSequence),
                    MenuItem::action("Open…", OpenFiles),
                    MenuItem::action("Open Folder…", OpenFolder),
                    MenuItem::separator(),
                    MenuItem::action("Save", Save),
                    MenuItem::action("Save As…", SaveAs),
                    MenuItem::action("Export…", Export),
                    MenuItem::separator(),
                    MenuItem::action("Close Tab", CloseTab),
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
                    MenuItem::action("Find Sequence", FindSequence),
                ]),
            ]);
            let options = WindowOptions {
                window_bounds: Some(WindowBounds::Windowed(Bounds::centered(
                    None,
                    size(px(1480.), px(920.)),
                    cx,
                ))),
                window_min_size: Some(size(px(1080.), px(700.))),
                app_id: Some("org.popgene.desktop".into()),
                ..TitleBar::window_options()
            };
            let (_, view) = gpui_kit::open_window(options, cx, move |window, cx| {
                window.set_window_title("PopGene");
                let view = cx.new(|cx| PopGene::new(paths, window, cx));
                let weak = view.downgrade();
                window.on_window_should_close(cx, move |_, cx| {
                    weak.update(cx, |view, cx| {
                        view.request_quit(cx);
                        view.model.allow_close
                    })
                    .unwrap_or(true)
                });
                view
            })
            .expect("Open PopGene window");
            // Keep the model alive through the window's Root entity.
            let _ = view;
            cx.activate(true);
        });
}

// Embed just the application icons in addition to the controls' default bundle.
gpui_kit::assets::icon_assets!(
    WorkspaceIcons,
    [
        Save,
        Download,
        Circle,
        MoveHorizontal,
        Scissors,
        ArrowBigRight,
        ArrowRightLeft,
        Ruler,
        Tags,
        Rows2,
        ListFilter
    ]
);
struct WorkspaceAssets;
impl AssetSource for WorkspaceAssets {
    fn load(&self, path: &str) -> Result<Option<std::borrow::Cow<'static, [u8]>>> {
        match WorkspaceIcons.load(path)? {
            Some(bytes) => Ok(Some(bytes)),
            None => gpui_kit::assets::Assets.load(path),
        }
    }
    fn list(&self, path: &str) -> Result<Vec<SharedString>> {
        let mut items = gpui_kit::assets::Assets.list(path)?;
        items.extend(WorkspaceIcons.list(path)?);
        Ok(items)
    }
}
