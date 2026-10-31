//! Sends rendered PDF pages to the platform's print dialog.
#[cfg(target_os = "macos")]
mod native {
    unsafe extern "C" {
        pub fn popgene_print_pdf(
            bytes: *const u8,
            length: usize,
            title: *const std::ffi::c_char,
        ) -> *const std::ffi::c_char;
    }
}

/// Called on GPUI's main thread. macOS shows its print sheet once the current
/// update finishes; Ok means a dialog or viewer is on its way, not that pages printed.
pub fn print_pdf(pdf: &[u8], title: &str) -> Result<String, String> {
    #[cfg(target_os = "macos")]
    {
        let title = std::ffi::CString::new(title.replace('\0', "")).unwrap_or_default();
        let error = unsafe { native::popgene_print_pdf(pdf.as_ptr(), pdf.len(), title.as_ptr()) };
        if error.is_null() {
            Ok("Print dialog opened".into())
        } else {
            Err(unsafe { std::ffi::CStr::from_ptr(error) }
                .to_string_lossy()
                .into_owned())
        }
    }
    #[cfg(not(target_os = "macos"))]
    open_in_viewer(pdf, title)
}

/// GPUI has no print API on Windows or Linux, so the system PDF viewer prints.
#[cfg_attr(target_os = "macos", allow(dead_code))]
fn open_in_viewer(pdf: &[u8], title: &str) -> Result<String, String> {
    use std::io::Write;
    let name: String = title
        .chars()
        .map(|c| if c.is_alphanumeric() { c } else { '-' })
        .take(60)
        .collect();
    let stamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| d.as_nanos());
    let path = std::env::temp_dir().join(format!(
        "popgene-print-{name}-{}-{stamp}.pdf",
        std::process::id()
    ));
    let mut options = std::fs::OpenOptions::new();
    options.write(true).create_new(true);
    // Sequence figures can be confidential; shared temporary folders must not expose them.
    #[cfg(unix)]
    std::os::unix::fs::OpenOptionsExt::mode(&mut options, 0o600);
    options
        .open(&path)
        .and_then(|mut file| file.write_all(pdf))
        .map_err(|e| format!("Could not prepare print pages: {e}"))?;
    let viewer = if cfg!(windows) {
        "explorer"
    } else {
        "xdg-open"
    };
    let mut child = std::process::Command::new(viewer)
        .arg(&path)
        .spawn()
        .map_err(|e| format!("Could not open the PDF viewer for printing: {e}"))?;
    // Reap the launcher; the viewer itself outlives it.
    std::thread::spawn(move || child.wait());
    Ok("Opened print pages in the PDF viewer · print from there".into())
}
