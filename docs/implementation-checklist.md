# Consolidated release acceptance

Status: 26 October 2026, `source0017` baseline plus a final circular-RNA import correction. Built functionality, verified behavior and remaining work are listed separately. An optimized local bundle has been packaged and exercised; this is not a claim of complete feature equivalence.

## Built

- [x] Sequence creation/editing, maps and sequence views, CDS/six-frame translation, ORFs and flanks.
- [x] Protein properties/models, primer design and thermodynamics, trace formats/editing/alignment, restriction analysis and physical digestion products.
- [x] Supported native DNA/RNA/protein and standard formats, multi-record input, references, attachments, original metadata and supported history restoration.
- [x] PCR, ligation, overlap, Golden Gate, alignments/consensus, bounded annotation ancestry and documented physical-fragment persistence.
- [x] Figure exports, transactional session/draft recovery, Undo/Redo and external-file conflict handling.

These checked items mean the supported implementations exist. The limitations below remain part of their contract.

## Verified

- [x] `source0017`: **443 tests on macOS and 443 on Linux ARM64**, strict workspace Clippy, formatting and debug builds; matching frozen-source manifests retained.
- [x] Windows GNU test compilation, strict lint and app cross-linking. This does not establish Windows runtime or MSVC behavior.
- [x] Independent biological/format checks and focused Linux GUI workflows covering translations, protein selections, primers, traces, restriction lists, cloning/alignment, exports, chooser behavior, recovery and conflict refusal.
- [x] macOS crash/relaunch recovered the verified session, unsaved products and original source bytes.
- [x] macOS actual GenBank/EMBL fragment Save/close/reopen/rejoin retained physical strands and reconstructed one phase-2 CDS. Its native saved product closed/reopened and displayed **WNSS** in both PopGene and the reference application using the same file.
- [x] Optimized macOS packaging, code-signature verification and launch of `dist/PopGene.app`; Grimsley and stable uncolored-primer fallback native UI acceptance passed.
- [x] Optimized macOS Open cancellation/retry, multiple-file opening and saved-copy reopen passed after restarting the isolated test app. The earlier invisible-panel stall remains an intermittent issue with undetermined cause.
- [x] Final circular-RNA correction: 172 formats tests, strict scoped Clippy, formatting and exact reference-created RNA import/byte-identical re-encode passed. The full macOS/Linux baseline and Windows cross-checks were not rerun for this narrow correction.
- [x] README separates verified behavior from pending checks and unsupported capabilities without competitor naming.

The evidence proves these cases, not every format field or workflow. Existing detailed checks need not be repeated unless a subsequent change affects them.

## Local delivery

- [x] Rebuild the optimized bundle with the final RNA correction, verify its signature/executable and reopen the exact RNA fixture in the rebuilt application: 16 nt circular RNA, its annotation and zero import errors.
- [x] Deliver the final `dist/PopGene.app` and record the consolidated changes in a local commit without pushing. No separate `/Applications` installation was found; the existing local bundle was replaced.

Preserve pre-existing user changes, including `validation/fixtures/synthetic-pcr.dna`, and retained recovery evidence. Final delivery status must be recorded from the actual bundle and running executable.

## Explicit limits

- Native export refuses fuzzy/remote/mixed-strand and other unsupported location operators. Native line encoding, automatic active-panel/visibility interpretation and unknown history action reconstruction remain unsupported.
- Unknown native segment semantics retain provenance and block unsafe rewrites. Stored RNA structures and original alignments are metadata inspection, not full structure editing or reconstruction.
- Recorded methylation flags are whole-template assumptions. EcoKI sensitivity and partial/site-specific methylation are not modeled.
- Annotation reconstruction is bounded and requires valid identity, mappings and metadata. Canonical Save or reopen of an unpersisted seed may create a new identity. Physical-fragment comment extensions do not make other applications understand that supplemental data.
- Trace precision/conversion and font-shaping limitations are reported explicitly. Experimental yield and measured secondary structure are not predicted.
- **Windows runtime and MSVC checks are unavailable:** the user confirmed there is no Windows test machine or VM. GNU cross-checks are not substitutes.
- Public distribution signing/notarization and third-party notice provenance review are separate from the local ad-hoc build/install task; see [releases](releases.md) and [third-party notices](../THIRD_PARTY.md).
