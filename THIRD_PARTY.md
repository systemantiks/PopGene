# Third-party software

macOS bundles additionally embed **Sparkle 2.10.0**, outside Cargo's dependency graph. Its upstream licence and bundled third-party notices are copied to `Contents/Resources/third-party/Sparkle-LICENSE.txt`, with version, source URL and archive SHA-256 in `Sparkle-version.json`. `scripts/sparkle.json` pins the official distribution. Sparkle supplies native update UI, download verification, installation and relaunch.

PopGene's source is GPL-3.0-or-later. Copyright 2026 PopGene contributors.

The native frontend uses **GPUI and GPUI Kit 0.8.0** (Apache-2.0). GPUI Kit's bundled Lucide icons carry ISC notices and, for Feather-derived icons, MIT notices; both are preserved alongside the Kit Apache licence. The sequence core and adapters continue to use gb-io (MIT), quick-xml (MIT), rfd (MIT), serde/serde_json (MIT OR Apache-2.0), thiserror (MIT OR Apache-2.0), rusqlite (MIT), SQLite (public domain), sha2 (MIT OR Apache-2.0) and tempfile (MIT OR Apache-2.0), with their dependencies. The exact dependency versions are in `Cargo.lock`. There is no JavaScript runtime, webview or egui frontend in the GPUI application.

`scripts/collect-notices.py` runs offline against the resolved Cargo graph for the current host, or `--target <triple>`. It collects upstream licence, copyright, notice and font licence text, then adds reviewed supplements from `vendor-notices/manifest.json`. It writes a fresh `dist/third-party` through a staged directory replacement: removed dependencies do not leave stale notices in a new build. A missing notice/supplement fails collection and preserves the previous complete inventory. `inventory.json` records authors, declared licences, notice paths, supplement files, the target and the Cargo lockfile hash. This is an inventory of resolved dependencies, including build dependencies, rather than a precise linked-code SBOM. The macOS inventory does not establish Windows or Linux notice coverage.

The current GPUI dependency notices include **seven explicit provenance gaps**, recorded in `vendor-notices/provenance-gaps.json` and the generated inventory:

| Packages | Preserved evidence and remaining gap |
| --- | --- |
| `block 0.1.6`, `malloc_buf 0.0.6`, `mac 0.1.1`, `objc_exception 0.1.2` | Their exact upstream release trees declare MIT or MIT/Apache-2.0 but contain no full licence notice. Original declarations, author attribution and source revisions are preserved. No full upstream notice has been invented. |
| `simd_helpers 0.1.0`, `seahash 4.1.0` | The published revisions lack a licence file. Later official upstream licence text is retained with an immutable revision and a clear provenance note. Its applicability to the published package remains a distribution-review item. |
| `pathfinder_geometry 0.5.1` | Its published revision is unavailable from the official repository. Shared project notices from the available revision used by `pathfinder_simd 0.5.6` are retained, together with the geometry crate's original copyright header and an explicit fallback note. |

`python3 scripts/collect-notices.py --strict` rejects these gaps. The local app bundle is not a cleared public distribution: resolve the listed provenance items and check native/generated code and shipped data before public binary redistribution. Historical egui and dispatch supplements remain in the source tree for earlier builds but are excluded from the new inventory when their packages are absent from the target's dependency graph.

Synthetic fixtures are original PopGene synthetic data. No SnapGene binaries, assets, curated sequence databases or source code are bundled. The SnapGene parser uses public format descriptions and Biopython SnapGeneIO as references. Attribution and upstream licence text are retained in [vendor-notices/biopython-reference](vendor-notices/biopython-reference/README.md) and included in the app's notices. Biopython and Python validation/packaging scripts are development tools, not application runtime dependencies.
