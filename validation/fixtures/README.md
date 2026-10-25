# Synthetic fixture provenance

All biological sequence content was created for PopGene. No patient, collaborator or experimental data is included. These fixtures exercise parsing and calculations; they are not laboratory constructs.

`../generate_fixtures.py` uses deterministic seed **823** to create a 5,000-base sequence with a small codon example and IUPAC ambiguities. `manifest.json` records expected primers, zero-based feature coordinates and the sequence hash. The GenBank journal line identifies the SnapGene GenBank dialect, not the creator of the original sequence.

| Files | Origin and purpose |
| --- | --- |
| `annotated.gb`, `annotated.fasta` | Generated original sequence; GenBank adds circularity, annotations and notes |
| `annotated.dna` | Minimal native representation of the original annotated sequence |
| `synthetic-pcr.dna` | Synthetic primer product; expected template slice `[1000:2024]` |
| `multi.fasta` | Deliberate multi-record parser edge case |
| `synthetic.ab1`, `synthetic-complete.ab1`, `basecaller-only.ab1` | Generated ABIF calls and signals, with edited/basecaller tag variants |
| `invalid.ab1`, `truncated.dna` | Deliberately malformed input |
| `regression/` | Original segmented annotations, uncached CDS translations, reading-frame offsets, protein sequence and invented ABIF instrument parameters |
| `native-records.json` | Inspectable source for every valid native DNA fixture |

Native syntax originated from serialization of these synthetic sequences by SnapGene 8.2.3. The public fixtures retain only header, sequence, primer, notes and feature packets. Enzyme catalogues, opaque history, display settings, UUIDs and creation/modification dates were removed. This provenance is retained honestly; no application code or curated vendor sequence library is included.

`python3 validation/generate_native_fixtures.py` reproduces all five valid `.dna` files directly from JSON using the Python standard library. Add `--check` to verify checked-in bytes. Regeneration does not require SnapGene. The intentionally truncated parser-error fixture is generated separately by `generate_fixtures.py`.

ABIF containers are written directly using Python `struct`. The calls are `ACGTN`, qualities 40/35/30/20/5 and peaks 2/5/8/11/14. Concrete-base signals contain amplitude-1000 peaks; they do not represent actual instrument output. `../prepare_trace_parameters.py` adds invented instrument tags for the regression example.

The regression GenBank and protein files are hand-authored. Segment expectations and current fixture hashes are in `regression/segment-colours-manifest.json`. Independent translation, restriction and protein expectation scripts live in the parent directory. Generated reports, comparison captures and large scaling inputs are ignored by version control.

Contributed fixtures must document creation method, redistribution permission, transformations and expected semantics. Use synthetic data and exclude personal metadata, opaque vendor catalogues and private experimental material.
