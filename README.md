# PopGene

**Free, open-source software for molecular biology.** Map and annotate plasmids, find restriction sites, design primers, plan cloning, analyze proteins and read Sanger traces in a native desktop app. PopGene works entirely offline: no account, no subscription and no network connection.

![PopGene displaying a circular DNA map with annotations, primers and restriction enzyme details](docs/images/popgene-circular-dna.png)

## Download

**[PopGene 0.0.1 preview](https://github.com/systemantiks/PopGene/releases/tag/v0.0.1)** is an early release. Keep copies of important files while you try it.

| Platform | Download | Requirements |
| --- | --- | --- |
| macOS, Apple silicon | [PopGene-0.0.1-macos-arm64.zip](https://github.com/systemantiks/PopGene/releases/download/v0.0.1/PopGene-0.0.1-macos-arm64.zip) | macOS 12 or later |
| Linux, ARM64 | [PopGene-0.0.1-linux-arm64.tar.gz](https://github.com/systemantiks/PopGene/releases/download/v0.0.1/PopGene-0.0.1-linux-arm64.tar.gz) | glibc 2.39 or later, X11 or Wayland |
| Windows, x64 | [PopGene-0.0.1-windows-x64.zip](https://github.com/systemantiks/PopGene/releases/download/v0.0.1/PopGene-0.0.1-windows-x64.zip) | Experimental: not yet tested on Windows |

Intel Macs, Linux x64 and Windows ARM64 are not packaged yet.

- **macOS:** unzip and move `PopGene.app` to Applications. The preview is not notarized, so macOS blocks the first launch; open **System Settings → Privacy & Security** and choose **Open Anyway**.
- **Linux:** extract the archive and run `./popgene`. Install the system libraries listed in `runtime-dependencies.txt` and a desktop portal for file dialogs.
- **Windows:** extract the archive and run `popgene.exe`. The build is unsigned; if SmartScreen appears, choose **More info → Run anyway**.

Keep the archive's licence and notice files with the app if you share it. PopGene opens with synthetic example data, so you can explore before loading your own files.

## Features

- **Maps and sequences.** Create or open DNA, RNA and protein documents. Linked circular and linear maps sit beside a wrapped sequence view, with features, primers and restriction sites on both, and a full-duplex inspector shows both strands and their overhangs. Search, select and edit bases, reverse-complement and change topology, with undo and redo for each document.
- **Restriction analysis.** 1,089 enzymes, custom enzyme lists, site-count filters, digests and fragment ends. Unknown cut positions are flagged rather than guessed.
- **Primers.** Primer3 design with explicit reaction settings, plus manually added primers and unbound oligos.
- **Cloning.** PCR, ligation, overlap assembly and Golden Gate, with annotations carried into the products.
- **Translation.** CDS and six-frame translation, ORF finding and all 27 current NCBI genetic codes.
- **Proteins.** Mass, charge and isoelectric point (Bjellqvist, EMBOSS or Grimsley), extinction coefficients, modifications and disulfides.
- **Sanger traces.** ABIF, SCF, ZTR and FASTQ, with quality tracks, trimming, base editing and alignment.
- **Alignment.** Pairwise and multiple alignment with consensus.
- **Figures.** Export maps and sequences as SVG, PDF or PNG for papers and posters.
- **Careful with your work.** Unsaved changes survive a crash, saves warn before dropping data a format can't hold, and PopGene won't silently overwrite a file that changed on disk.

## File formats

PopGene opens and saves native `.dna`, `.rna` and `.prot` documents, GenBank/GenPept, EMBL and FASTA, and Sanger traces in ABIF, SCF, ZTR and FASTQ. It also reads legacy `.popgene` documents and can browse a whole project folder.

Format support is partial. Unchanged files keep their original bytes when saved in the same format, and PopGene tells you what will be lost before an edited or converted save. Read [file formats and limitations](docs/file-formats.md) before relying on PopGene for important data.

## Project status

PopGene is under active development, and 0.0.1 is a preview:

- Some format details and workflows are incomplete or unsupported.
- Calculations follow documented models and assumptions; they do not predict experimental yield or measured structure.
- The macOS app is not notarized, and the Windows build has not yet run on Windows.
- An intermittent stall of the macOS Open dialog has been seen; its cause is unknown.

[Status and verification](docs/status.md) records what has been tested on each platform and the remaining work.

## Build from source

PopGene is written in Rust with GPUI and GPUI Kit. Install Rust **1.99.1** and your platform's C compiler and SDK, then run:

```sh
cargo run --locked --release -p popgene-desktop
```

To build a macOS app bundle, run `scripts/package-macos.sh`. The bundle embeds Sparkle for **PopGene → Check for Updates…**; the first packaging run downloads a pinned, checksum-verified Sparkle distribution. Local builds are ad-hoc signed, not notarized. See [releasing and updates](docs/releases.md) for Developer ID signing, notarization, GitHub release drafts and the update feed.

## Contributing

| Crate | Responsibility |
| --- | --- |
| `popgene-core` | Sequence model, coordinates, search, edits and undo |
| `popgene-formats` | File codecs, legacy import, folder inventory and atomic saves |
| `popgene-desktop` | Native workspace, linked views and background file I/O |

Run the checks before sending changes:

```sh
cargo test --locked --workspace
cargo clippy --locked --workspace --all-targets -- -D warnings
cargo fmt --all --check
```

Python is used only for developer tooling. See [development and validation](docs/development.md) for the scientific reference checks and [synthetic fixture provenance](validation/fixtures/README.md). Contributions should include focused tests for changed sequence or file semantics and preserve the licence notices in `vendor-notices/`.

## License

PopGene is free software, licensed under [GPL-3.0-or-later](LICENSE). [Third-party notices](THIRD_PARTY.md) cover bundled dependencies and record the reviewed notice provenance gaps.
