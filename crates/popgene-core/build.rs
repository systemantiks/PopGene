fn main() {
    let root = "vendor/primer3";
    let mut build = cc::Build::new();
    build.include(root).include(format!("{root}/klib"));
    build.flag_if_supported("-std=c99").warnings(false);
    // Primer3's masker uses POSIX getline; strict C99 otherwise hides its
    // declaration on glibc, where current GCC rejects implicit declarations.
    if std::env::var_os("CARGO_CFG_UNIX").is_some() {
        build.define("_POSIX_C_SOURCE", "200809L");
    }
    for file in [
        "libprimer3",
        "p3_seq_lib",
        "dpal",
        "thal",
        "thal_parameters",
        "oligotm",
        "amplicontm",
    ] {
        build.file(format!("{root}/{file}.c"));
    }
    // Upstream libprimer3 disables its optional genome masker on Windows.
    // Its POSIX mmap/getline implementation is not needed by our design API.
    if std::env::var("CARGO_CFG_TARGET_OS").as_deref() != Ok("windows") {
        build.file(format!("{root}/masker.c"));
    }
    build.file("src/primer_bridge.c").compile("popgene_primer3");
    println!("cargo:rerun-if-changed=vendor/primer3");
    println!("cargo:rerun-if-changed=src/primer_bridge.c");
    if std::env::var("CARGO_CFG_TARGET_ENV").as_deref() != Ok("msvc") {
        println!("cargo:rustc-link-lib=m");
    }
}
