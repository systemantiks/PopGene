//! Open a sequence document or a project folder from one entry point.
use gpui_kit::{App, Window};
use std::{future::Future, path::PathBuf};

type OpenResult = Result<Option<Vec<PathBuf>>, String>;

/// macOS supports files and folders in the same asynchronous native panel.
/// Start it before spawning the result handler; the platform defers presentation
/// until the caller's GPUI entity lease has been released.
#[cfg(target_os = "macos")]
pub fn prompt(_window: &Window, cx: &mut App) -> impl Future<Output = OpenResult> + use<> {
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
pub fn prompt(window: &Window, _cx: &mut App) -> impl Future<Output = OpenResult> + use<> {
    // Capture the owning window before yielding. The owned builder copies its
    // native handles; the chooser is presented only when the future is polled.
    let dialog = rfd::AsyncFileDialog::new().set_parent(window);
    async move {
        Ok(dialog
            .add_filter(
                "Sequence documents",
                &[
                    "dna", "rna", "prot", "gb", "gbk", "genbank", "gbff", "gbank", "gbs", "gen",
                    "ape", "gp", "genpep", "genpept", "gpff", "embl", "emb", "ebl", "dat", "fasta",
                    "fa", "fas", "fna", "ffn", "faa", "frna", "ab1", "abi", "abd", "abf", "scf",
                    "ztr", "fastq", "fq", "popgene",
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
