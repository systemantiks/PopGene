# Annotation ancestry convention — development version 1

This describes PopGene's readable provenance extension for annotations carried through sequence construction. It is a new PopGene convention in ordinary comments, not a claim that another application implements the same ancestry model. Implementation and acceptance are still in progress; see the [release checklist](implementation-checklist.md).

Ordinary ligation reconstructs an annotation only when its pieces have the same recorded origin and copy identity, unchanged annotation and citation baselines, an unambiguous ordered placement, and matching sequence. Verified exact-overlap synthesis may combine distinct copied occurrences of the same origin into a fresh product copy under its explicit overlap rules. Repeated or overlapping location leaves remain distinct occurrences. Equal labels or sequence alone do not establish shared ancestry. Complementary cohesive strands may cover the same product coordinate; duplicated coverage of the same physical strand does not establish a valid merge.

The registry stores annotation descriptors, coordinate mappings and digests. It does not store an ancestor `Document`, omitted parent DNA, original file bytes or an encoded document snapshot. Current product sequence, measured qualities and physical end chemistry remain the product's own data.

## Text carrier

Supported carrier locations are GenBank `COMMENT`, EMBL `CC`, and the existing native Notes/Comments field. No private packet type or XML attribute is introduced.

```text
[PopGene lineage v1 begin]
bytes=N;sha256=HEX;json=JSON
[PopGene lineage end]
```

`N` is the UTF-8 byte length of `JSON`; `HEX` is its lowercase SHA-256 digest. The JSON is compact, with whitespace inside string values escaped using JSON escapes (`\n`, `\r`, `\t` or Unicode escapes as appropriate). The writer wraps the header and payload at UTF-8 character boundaries into lines of at most 48 bytes. The reader removes ASCII whitespace from the header/payload before checking the length and digest; other wrapping whitespace is rejected. The marker lines are not part of the digest.

Exactly one matching marker pair is accepted. Unknown versions, duplicate or colliding markers, malformed JSON, failed checksums and inconsistent records remain readable comment text with a warning and no active ancestry. A checksum detects corruption; it does not establish authenticity or authorize restoration.

The payload object has exactly four keys:

| Key | Meaning |
| --- | --- |
| `format` | `genbank`, `embl` or `native`, matching the actual file codec. |
| `records_sha256` | Digest of the actual emitted annotation records, including their order. |
| `context_sha256` | Digest of the canonical decoded sequence, geometry, annotations, colors and references. |
| `registry` | The typed registry described below. |

For native files, the record digest includes each Feature/Primer packet in file order: its one-byte type, its payload length as an unsigned 64-bit big-endian integer, then the payload bytes. For GenBank it includes complete lines from `FEATURES` up to, but excluding, `ORIGIN` or `CONTIG`. For EMBL it includes each complete `FT` line. Line endings are included. Thus external record reordering or rewriting can make otherwise readable provenance inactive; the importer does not guess associations by matching identical annotations.

The context object contains `sequence_sha256`, `physical_reference_sha256`, `geometry`, `molecule`, `single_stranded`, `circular`, `features`, `duplex_features`, `sequence_colors`, `duplex_colors` and `references`. Sequence digests cover the current UTF-8 sequence letters. Geometry is `[extent, top_range, opposite_range]` on the shared physical axis. New context digests use compact JSON with object keys recursively sorted, independent of build features; array order and scalar values remain unchanged. For earlier development files, the reader also accepts the exact insertion-order encoding: the context fields in the order listed here, with nested typed structs in their declared serialization order. It reconstructs both hashes from the same decoded model; no alternate biological interpretation or arbitrary object ordering is accepted. This is a semantic binding in addition to the raw-record digest.

## Registry and coordinate records

All coordinates are zero-based and half-open. A range is `{"start": S, "end": E}`. UUIDs identify instances and copies independently of sequence content. Digests are hexadecimal SHA-256 strings.

| Registry field | Meaning |
| --- | --- |
| `version` | Integer `1`. |
| `instance_id` | UUID identifying this document instance, including before its first digestion. |
| `tracking_limited` | Whether some ancestry is unavailable; prevents an untracked clipped descendant from being silently seeded as a new complete ancestor. |
| `feature_epochs` | Pairs of current feature ID and creation UUID, distinguishing reused local IDs. |
| `origins` | Flat list of original annotation descriptors. |
| `pieces` | Associations between current annotations and original ordered occurrences. |

Each origin contains `id`, `source_id`, `occurrence_id`, `top`, `bottom`, `feature`, `references`, `extent`, `molecule`, `circular`, `source_digest` and `sequence_digest`. `feature` is the original annotation descriptor, including its ordered location tree, qualifiers and display metadata. `references` preserves cited objects at their original numeric slots, with empty placeholders for unused earlier slots. Strand ranges describe the original physical extent. The sequence digest covers the original annotation's ordered extracted sequence; the source digest fingerprints the source molecule without retaining its bases.

Each piece contains `mapping_digest`, `origin_id`, `copy_id`, `feature_id`, `axis`, `runs`, `baseline`, `references_digest` and `sequence_digest`. `axis` is `Top` for an ordinary sequence annotation or `Duplex` for a full-extent annotation. `baseline` is the exact current annotation descriptor. Citation and sequence digests bind that baseline to the current document. These are checked before any coordinate transformation or codec normalization refreshes a baseline.

Each run contains:

| Run field | Meaning |
| --- | --- |
| `leaf` | Zero-based depth-first Range occurrence in the original location tree; duplicate leaves have different indices. |
| `offset` | Offset from that original Range's low coordinate. |
| `length` | Positive number of mapped residues. |
| `target_start` | Low coordinate on the current full physical extent. |
| `reverse` | Whether the placement reverses original coordinate order. |
| `strands` | Bitmask: `1` top physical residues, `2` opposite physical residues, `3` both. |

An annotation descriptor has `id`, `label`, `kind`, `parts`, optional `location`, `reverse`, `qualifiers` and `edited`. Qualifiers are ordered pairs of key and optional value. The location tree uses readable tagged JSON such as `{"Range":{"start":0,"end":15,"before":false,"after":false}}`, `{"Complement": ...}` or `{"Join":[...]}`. Current reconstruction supports Range/Complement/Join trees and retains genuine fuzzy endpoints; other operators are not guessed. The full model definitions are in [core locations](../crates/popgene-core/src/location.rs), [annotations](../crates/popgene-core/src/lib.rs) and [references](../crates/popgene-core/src/metadata.rs).

Core tuple digests use compact JSON arrays; nested structs use their declared field order and omit fields explicitly marked for omission in their model definitions. The digest inputs are:

| Digest | Input |
| --- | --- |
| `occurrence_id` | `[instance_id, original_feature_id, epoch_or_null]` |
| `source_digest` | `[top_sequence, opposite_5to3_or_null, molecule, circular]`; only the resulting digest is retained. |
| Origin `id` | `[source_id, occurrence_id, feature, references, extent, molecule, circular, source_digest, sequence_digest, top, bottom]` |
| Piece `mapping_digest` | `[runs, axis, [extent, top, bottom], molecule, circular, single_stranded]` |
| `references_digest` | Ordered citation entries, each `[qualifier_value_or_null, resolved_reference_or_null]`. |
| Both `sequence_digest` fields | Raw UTF-8 biological extraction, including repeated ordered leaves, rather than JSON. |

The exact struct field order, validation and coordinate transforms are defined by the readable [core registry implementation](../crates/popgene-core/src/lineage.rs). File record association, canonical projection and comment parsing are defined by the [carrier implementation](../crates/popgene-formats/src/lineage.rs). Nested JSON fields not recognized by the current registry are rejected rather than silently ignored when activating a carrier.

## Bounds and incomplete ancestry

The comment payload is limited to 4 MiB and 96 JSON nesting levels. Core limits are 2,048 origins, 8,192 pieces, 8,192 feature-creation epochs, 65,536 runs, 8 MiB of serialized registry data, two million original and mapped residue occurrences, and sixteen million occurrence-validation steps. Ordinary sequence operations retain their existing size limits: optional ancestry tracking must not turn a previously valid large sequence operation into a failure. A bounded tracking failure retains ordinary projected annotations and reports unavailable ancestry. Exhausting the epoch limit disables optional new tracking without preventing annotation creation.

Session/recovery serialization retains the document identity even when no pieces exist. Ordinary unchanged public source files are not rewritten solely to add a seed. Canonical Save or reopening such a source without a persisted seed creates a new identity; its newly digested fragments cannot be merged with earlier fragments merely because their sequence matches. Existing valid seed carriers are retained.

Files that cannot carry active ancestry report its loss. Save uses the actual decoded saved representation, and Undo can recover the richer prior document; exporting a separate copy does not change the open document. A limited-tracking marker is retained where the target has capacity; it does not contain enough evidence to reconstruct annotations. If even that small marker cannot fit, ordinary supported data is still saved with an explicit warning that the marker, tracking limitation and source identity will not survive.

The carrier does not repair a target format's other limitations. Native fuzzy-location export remains unsupported. GenBank/EMBL now have a separate [physical-fragment comment convention](physical-fragments.md) that preserves supported current geometry before ancestry is validated; ordinary readers still see the actual top-strand record. Codec tests verify saving and reopening cohesive CDS fragments before rejoining. Desktop, platform and GUI acceptance remain in progress. Fully represented annotations and physical geometry must survive the codec before their ancestry can remain active. Third-party retention and reference-application behavior require separate acceptance evidence.
