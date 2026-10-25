# Development and validation

Use the toolchain and commands in the [README](../README.md). `scripts/cargo-local` optionally selects an ignored `.tools/` installation; ordinary checkouts can use `cargo` directly. Dependencies must be downloaded once before offline builds work.

## Regression data

Checked-in fixtures are synthetic and self-contained. Native DNA fixtures are generated from inspectable JSON; no third-party application is needed to regenerate them:

```sh
python3 validation/generate_native_fixtures.py --check
python3 validation/generate_fixtures.py
python3 validation/generate_native_fixtures.py
python3 validation/prepare_trace_parameters.py
```

The general generator also makes ignored larger FASTA inputs for local performance checks. `python3 scripts/make-demo-project.py` creates an ignored nested project folder, including a deliberately truncated file for import-error handling.

Independent biological expectation checks use pinned developer dependencies:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r validation/requirements.lock
.venv/bin/python validation/check_restriction.py
.venv/bin/python validation/check_translation_fixtures.py
.venv/bin/python validation/check_protein_properties.py
```

These scripts verify fixture expectations independently; Rust unit tests verify application calculations against those expectations. Reports are written to ignored `validation/results/`. Do not commit screenshots of other applications, local logs, credentials, private sequences or generated project files. See [fixture provenance](../validation/fixtures/README.md) before adding test data.

## UI checks

Exercise Map, Sequence, Enzymes, Features, Primers, History, Trace, Parameters and Properties with the relevant document types. Check selection linking, split-pane resizing, search, filter changes and scrolling. The map and sequence filters belong in thin sidebars beside the viewer; sidebars should scroll independently and controls should not unexpectedly move the viewer boundary. Check narrow windows as well as the default size.

Test new/open/save/export and dirty-document close prompts with disposable synthetic documents. Platform checks and manual observations should be reported with their actual scope; configured CI alone is not evidence of a successful run.

## Packaging and notices

`scripts/package-macos.sh` builds and signs the macOS bundle. `python3 scripts/collect-notices.py` collects licence text from the resolved dependency graph and reviewed supplements. `--strict` also rejects known provenance gaps; see [THIRD_PARTY.md](../THIRD_PARTY.md). Preserve upstream copyright and licence notices when changing dependencies.
