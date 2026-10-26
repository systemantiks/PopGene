# Physical fragment comment convention — v1

This is a new, public PopGene convention in GenBank `COMMENT` and EMBL `CC`
text. It is not a native file-format attribute, INSDC physical-end standard, or
claim that another application implements the convention. This convention is implemented in the development tree; final platform and
application acceptance remain in progress.

The ordinary record's `ORIGIN`/`SQ` is **the actual top strand**, 5′→3′. Ordinary
feature coordinates remain on that strand. No complement is inserted into
ordinary sequence positions where a top physical residue is absent. Readers
that ignore this extension see the supported top-strand record; full duplex
geometry and annotations may be visible only as comment text in those readers.

## Purpose and scope

The supplement records a current linear, double-stranded DNA fragment without
trace data. It carries actual opposite-strand residues, end geometry/chemistry,
current qualities and colors, full-extent annotations, and current bibliography.
It does not contain an ancestor sequence, an ancestor `Document`, an encoded
document snapshot, original file bytes, compressed/base64 data, or executable
instructions. Annotation ancestry remains a separate independently validated
carrier. Physical metadata is restored before ancestry is considered.

Native fuzzy-location writing is unchanged. GenBank/EMBL location syntax remains
ordinary syntax for top-strand features; a full-extent feature that cannot fit on
the actual top sequence is carried explicitly as a supplemental descriptor.
Fuzzy endpoints, ordered/repeated location leaves and biological strand are not
flattened to fit an ordinary top feature.

## Envelope

```text
[PopGene physical fragment v1 begin]
bytes=N;sha256=HEX;json=JSON
[PopGene physical fragment end]
```

`N` is the UTF-8 byte count of `JSON`; `HEX` is its lowercase SHA-256 digest.
JSON is compact. The serializer's standard JSON escapes (such as `\n`, `\r` and
`\t`) are retained; any whitespace still literal in serialized string values is
escaped as JSON `\uXXXX` sequences (including surrogate pairs where needed).
No literal whitespace remains inside strings, and the compact serializer
inserts no whitespace outside strings. The writer wraps at UTF-8
character boundaries, with at most 48 bytes per payload/header line. The reader
removes only ASCII wrapping whitespace. Non-ASCII wrapping is invalid. Checksums
provide integrity, not authenticity.

Exactly one complete marker pair is recognized. Unknown versions, duplicate or
colliding markers, malformed/deep/oversized JSON, checksum failures and semantic
inconsistencies remain readable text with a warning and **no activated physical
metadata**. Staging is atomic: a failed block cannot partly change any strand,
annotation, reference, chemistry, color or confidence value.

## JSON fields

The envelope has exactly four keys:

| Key | Meaning |
| --- | --- |
| `format` | `genbank` or `embl`, matching the actual codec. |
| `records_sha256` | Digest of the ordinary emitted annotation records, in order. |
| `context_sha256` | Digest of the canonical ordinary top-record context below. |
| `fragment` | The current physical metadata object described below. |

For GenBank the record digest covers complete lines from `FEATURES` up to but
excluding `ORIGIN`/`CONTIG`. For EMBL it covers complete `FT` lines. Line endings
are included. New context hashes use compact JSON with recursively sorted object keys, independent of build features, containing
`name`, `description`, `accession`, `sequence_sha256`, `molecule`, `circular`,
`single_stranded`, `features`, and `references` from the *ordinary* decoded
record. It excludes comments to avoid recursive cross-binding of comment
carriers. Sequence digests cover the canonical current UTF-8 letters. For earlier
development files, the reader also accepts the exact insertion-order encoding:
the context fields in the order listed above, with nested typed structs in
their declared serialization order. Both hashes use the same decoded model;
array order, scalar values and every existing semantic check remain unchanged.

Therefore external edits to sequence, topology, strand kind, ordinary feature
records/order, names/descriptions/accession or reference objects make the old
supplement inactive. Arbitrary prose comment changes do not change biological
context. The reader never repairs association by matching labels or sequence.

The `fragment` object has exactly these fields:

| Key | Meaning |
| --- | --- |
| `ends` | Current `DuplexEnds`: left/right end kind, actual protruding sequence, source enzyme names, actual `opposite_sequence` in 5′→3′ order, optional left/right modification strings. |
| `duplex_annotations` | Full-extent Feature descriptors and nucleotide colors. Coordinates span the union of actual physical strands. |
| `sequence_colors` | Current top-axis nucleotide colors; each records its physical strand. |
| `sequence_qualities` | Nullable vector in top 5′→3′ order, or `null` when unavailable. |
| `opposite_sequence_qualities` | Nullable vector in opposite 5′→3′ order, or `null`. |
| `methylation` | Recorded whole-template Dam/Dcm/EcoKI assumptions, or `null`; no occupancy or efficiency is inferred. |
| `references` | Exact current structured Reference objects in numeric citation order. This preserves current bibliographic fields that ordinary header syntax may normalize. |

Core descriptors use their readable serde fields as defined by the existing
model. Nested unknown fields are rejected through strict typed reserialization;
none are silently discarded during activation. Quality values are integers
0–93 or `null`. Phred zero remains distinct from unknown. Optional whole vectors
remain distinct from vectors of unknown entries. Unknown terminal modification
strings are retained as stated metadata; cloning still refuses chemistry it
cannot interpret.

Ordinary decoded features receive IDs `1..N` in actual serialized occurrence
order. Supplemental feature IDs must be exactly `N+1..N+M` in payload order.
Duplicate labels/content remain separate occurrences. Writers map each original
top/full-extent feature ID explicitly to this canonical occurrence; ancestry
transport uses this mapping only after physical reimport succeeds. Reference
count and slots must agree with the ordinary reference table; original exact
Reference objects are restored before validating every current citation. When
the GenBank display writer appends its known software reference, the bridge
verifies every field and its terminal ordinal before appending it to the
supplement. Original scientific slots remain unchanged, and their citations
are validated before the append so an invalid citation cannot acquire meaning.

## Biological validation and limits

Before activation, existing `DuplexEnds::validate` verifies paired reverse
complement, actual overhang residues, nonempty paired region, end polarity and
strand lengths. `Document::validate` verifies linear dsDNA intent, local location
trees/projections, coordinate bounds, unique feature IDs, RGB/color bounds,
quality-vector lengths and scores. All current citation numbers must resolve.
The full extent is a coordinate reference, not an additional physical strand.

Payload size is at most 4 MiB, nesting at most 96, and supplemental feature count
at most 8,192. Borrowed field serialization is counted against the budget before
cloning metadata or building a JSON value tree, so an optional carrier cannot
allocate unbounded metadata just to reject it later. Existing sequence/file and
core annotation limits remain in force.

If optional preparation cannot preserve the supplement or it exceeds budget,
the ordinary valid top record is still saved with explicit physical-data loss
warnings. Ancestry is then evaluated against the *actual decoded* result and
becomes limited/inactive if required geometry is absent. A conflicting existing
marker together with active current physical ends instead refuses the save:
returning text containing a stale valid block could otherwise reactivate old
chemistry. The user must preserve or remove that conflicting text explicitly.
When no physical state is active, ordinary serialization is also checked before
saving retained marker-like text: if whitespace or annotation formatting would
make formerly inactive metadata active again, saving is refused. The analogous
guard applies to retained inactive ancestry, including seed-only carriers.
Unknown versions and still-inactive comments remain savable.

Save canonicalization, close/reopen and unchanged-byte reencoding must use the
actual file result. Export does not mutate the open document. Undo/recovery are
separate capabilities. Removing or rewriting the carrier in another application
necessarily loses the supplement; equivalent sequence alone cannot restore it.

## Validation coverage and remaining acceptance

- Actual disk Save/reopen of both EcoRI fragments of `ATGGAATTCAAGTAA`, then
  ligation to one original CDS, phases1/2/3 giving `MEFK*`, `WNSS`, `GIQV`.
- Reverse reassembly and circular single-cut closure; all four terminal
  polarity combinations; staged and mixed-source products; duplicate/repeated
  Join occurrences and genuine original partial endpoints.
- Exact current strand sequences, chemistry, nullable scores, colors,
  citations and unchanged-byte reencode, including independent JSON/raw
  record verification and a standards-only parser seeing the actual top.
- Failed/stale/stripped/unknown/duplicate blocks and forged inconsistent
  payloads: ordinary data stays readable, physical state is never partly
  activated, and ancestry never survives missing geometry.
- Capacity fallback with genuinely actionable input ancestry; explicit Save
  canonicalization and Undo/recovery; final platform and GUI acceptance.

The physical-carrier checkpoint passed 167 formats tests. The subsequent desktop checkpoint passed 137 tests, including actual Save As/reopen/rejoin and SQLite recovery with Undo/Redo for both text formats. These richer cases preserve original citation slots alongside a generated display reference, both strand colors and nullable qualities. An independent Biopython/raw-JSON check passes 18 emitted files. A new platform checkpoint and GUI save/reopen/rejoin remain separate acceptance gates; later fixes require fresh integrated checks.
