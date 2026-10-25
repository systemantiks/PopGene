# PopGene

A native, offline sequence workbench built with Rust, GPUI and GPUI Kit. Licensed under [GPL-3.0-or-later](LICENSE).

PopGene opens synthetic example data on launch and works with local files without an account or network connection. Python is used only for developer tooling.

## Features

- Create DNA, RNA and protein documents; open multiple files or browse a project folder.
- Linked circular/linear maps and wrapped sequence views, with an adjustable split view. Thin left and right viewer sidebars hold filters, display options, zoom, search and range tools.
- Inspect features, primers, restriction sites, document metadata, Sanger traces and protein properties.
- Select bases and annotations, search sequences, edit DNA, reverse-complement and change topology, with per-document undo/redo.
- Read and write supported `.dna`, GenBank / GenPept, EMBL and FASTA data. Read ABIF traces and legacy `.popgene` documents. Save and export show format-specific loss warnings.

Displayed coordinates are 1-based inclusive. The restriction catalogue contains 56 enzymes. Calculated CDS translations support genetic codes 1, 2 and 11; six-frame tracks use the standard code.

PopGene is under active development. Unchanged same-format saves preserve source bytes; edited or converted saves preserve only the target codec's supported data. See [file formats and limitations](docs/file-formats.md).

## Implementation status and remaining work

Status updated 23 October 2026. Feature coverage and format support remain incomplete.

| Area | Built | Remaining |
| --- | --- | --- |
| Translations | Calculated CDS translations, six-frame tracks, codon selection and native CDS frame offsets. | ORF discovery, CDS-flanking translations, additional genetic codes, translation exceptions and circular-boundary handling. |
| Protein properties | Composition, mass, named Bjellqvist pI/charge, extinction and A280 estimates. | Additional validated mass and pI/charge calculation models, properties for selected residues, secondary structure and modifications. |
| Annotations and display | Feature segments and qualifiers; segment names and colours in maps, sequence and details. | Remaining segment/cleavage styles, complex locations, extended GenBank annotations and native display settings. |
| Primers | Recorded oligo sequence, binding details, search, copy and labelled imported Tm values. | Validated Tm calculations and reaction parameters, unbound/complex bindings, primer design and editing. |
| Traces | ABIF traces and read-only instrument parameters. | Complete instrument records, more trace formats, alternate base calls, trimming, editing and alignment. |
| Restriction analysis | Site inspection and a 56-enzyme catalogue. | Full enzyme catalogue, methylation/reaction conditions and digestion simulation. |
| File compatibility and metadata | Supported sequence codecs, original headers/references and native Notes fields. | Complete `.dna` record reconstruction, `.rna`/`.prot` containers, multi-record input, fuzzy/remote locations, structured reference editing/conversion, attachments and interpreted procedure history. |

Broader work remains for cloning, sequence alignment, interactive procedure history, publication exports, autosave and session restoration. Windows and Linux still need verification.

**Recommended next priorities:** close the translation gaps and improve file round-trip fidelity. Validate each change with synthetic fixtures and independent expectations, and update this table as gaps are implemented and verified.

## Build

Install Rust **1.99.1** and your platform's C compiler and SDK, then run:

```sh
cargo run --locked --release -p popgene-desktop
cargo test --locked --workspace
cargo clippy --locked --workspace --all-targets -- -D warnings
cargo fmt --all --check
```

The native UI has been exercised on macOS ARM64. Windows and Linux are not yet verified; the [CI workflow](.github/workflows/rust.yml) defines checks and Linux package prerequisites.

For local macOS packaging, run `scripts/package-macos.sh`. The bundle embeds Sparkle for **PopGene → Check for Updates…**; the first packaging run downloads a pinned, checksum-verified Sparkle distribution. Local builds are ad-hoc signed, not notarized. See [releasing and updates](docs/releases.md) for Developer ID signing, notarization, GitHub release drafts and the update feed. [Third-party notices](THIRD_PARTY.md) document unresolved dependency-notice provenance items that need review before public binary redistribution.

## Development

| Crate | Responsibility |
| --- | --- |
| `popgene-core` | Sequence model, coordinates, search, edits and undo |
| `popgene-formats` | File codecs, legacy import, folder inventory and atomic saves |
| `popgene-desktop` | Native workspace, linked views and background file I/O |

See [development and validation](docs/development.md) and [synthetic fixture provenance](validation/fixtures/README.md). Contributions should include focused tests for changed sequence or file semantics and preserve the licence notices in `vendor-notices/`.
