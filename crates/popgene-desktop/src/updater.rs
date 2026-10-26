//! Main-thread bridge to Sparkle and AppKit's termination request.
#[cfg(target_os = "macos")]
mod native {
    unsafe extern "C" {
        pub fn popgene_install_quit_guard() -> std::ffi::c_int;
        pub fn popgene_take_quit_request() -> std::ffi::c_int;
        pub fn popgene_reply_to_quit(approved: std::ffi::c_int) -> std::ffi::c_int;
        pub fn popgene_check_for_updates() -> *const std::ffi::c_char;
    }
}

/// All entry points are called on GPUI's main thread. No native callback borrows
/// GPUI state: the foreground task drains termination requests instead.
#[cfg(target_os = "macos")]
pub fn install_quit_guard() -> bool {
    unsafe { native::popgene_install_quit_guard() != 0 }
}

#[cfg(target_os = "macos")]
pub fn take_quit_request() -> Option<bool> {
    match unsafe { native::popgene_take_quit_request() } {
        1 => Some(true),
        2 => Some(false),
        _ => None,
    }
}

/// True means a postponed Sparkle installation will resume asynchronously.
pub fn reply_to_quit(approved: bool) -> bool {
    #[cfg(target_os = "macos")]
    unsafe {
        native::popgene_reply_to_quit(i32::from(approved)) != 0
    }
    #[cfg(not(target_os = "macos"))]
    {
        let _ = approved;
        false
    }
}

pub fn check_for_updates() -> Result<(), String> {
    #[cfg(target_os = "macos")]
    unsafe {
        let error = native::popgene_check_for_updates();
        if error.is_null() {
            Ok(())
        } else {
            Err(std::ffi::CStr::from_ptr(error)
                .to_string_lossy()
                .into_owned())
        }
    }
    #[cfg(not(target_os = "macos"))]
    Err("Automatic installation is currently available only on macOS. Download releases from https://github.com/systemantiks/PopGene/releases.".into())
}
