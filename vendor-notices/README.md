# Dependency notice supplements

`manifest.json` maps exact crate versions to local supplement directories. The offline collector only includes mappings for packages present in the target's resolved Cargo graph. Historical mappings are retained to reproduce earlier PopGene builds.

New GPUI migration supplements include `SOURCE.json` with the official repository, immutable revision, affected packages and SHA-256 hashes of unmodified downloaded files. `SOURCE.txt` provides a readable provenance record. Existing supplements predate the JSON format and retain their original source-revision notes.

`provenance-gaps.json` distinguishes missing full upstream notices from fallback or later-revision notices. A manifest declaration and an author attribution are evidence, not a reconstructed licence. `collect-notices.py --strict` fails while any such gap is present in the selected target graph. Normal collection includes the evidence and reports each gap for local builds; it does not assert public redistribution clearance.

The GPUI Kit supplement supplies its Apache-2.0 licence to both `gpui-kit` and `gpui-kit-assets`. The assets crate's own Lucide/Feather notice is also collected. Source files whose names happen to contain “copying” are not treated as licence files.
