//! Open a sequence document or a project folder from one entry point.
use gpui_kit::App;
use std::{future::Future, path::PathBuf};

type OpenResult = Result<Option<Vec<PathBuf>>, String>;

/// macOS supports files and folders in the same asynchronous native panel.
/// Start it before spawning the result handler; the platform defers presentation
/// until the caller's GPUI entity lease has been released.
#[cfg(target_os = "macos")]
pub fn prompt(cx: &mut App) -> impl Future<Output = OpenResult> + use<> {
    let response = cx.prompt_for_paths(gpui_kit::PathPromptOptions {
        files: true,
        directories: true,
        multiple: true,
        prompt: Some("Open".into()),
    });
    async move {
        response
            .await
            .map_err(|error| format!("The Open dialog closed unexpectedly: {error}"))?
            .map_err(|error| format!("Could not open the file chooser: {error}"))
    }
}

/// Windows and Linux native pickers cannot select files and folders together.
/// Their Open menu provides a separate folder command.
#[cfg(not(target_os = "macos"))]
pub fn prompt(_cx: &mut App) -> impl Future<Output = OpenResult> + use<> {
    async move {
        Ok(rfd::AsyncFileDialog::new()
            .add_filter(
                "Sequence documents",
                &[
                    "dna", "gb", "gbk", "genbank", "gbff", "gbank", "gbs", "gen", "ape", "gp",
                    "genpep", "genpept", "gpff", "embl", "emb", "ebl", "dat", "fasta", "fa", "fas",
                    "fna", "ffn", "faa", "frna", "ab1", "abi", "abd", "abf", "popgene",
                ],
            )
            .pick_files()
            .await
            .map(|files| {
                files
                    .into_iter()
                    .map(|file| file.path().to_owned())
                    .collect()
            }))
    }
}
