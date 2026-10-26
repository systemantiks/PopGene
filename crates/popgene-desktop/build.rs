fn main() {
    println!("cargo:rerun-if-changed=src/updater_macos.m");
    #[cfg(target_os = "macos")]
    if std::env::var("CARGO_CFG_TARGET_OS").as_deref() == Ok("macos") {
        cc::Build::new()
            .file("src/updater_macos.m")
            .flag("-fobjc-arc")
            .flag("-mmacosx-version-min=12.0")
            .warnings_into_errors(true)
            .compile("popgene_updater");
        println!("cargo:rustc-link-lib=framework=AppKit");
    }
}
