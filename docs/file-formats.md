# File formats

PopGene saves interoperable sequence files. New annotated documents default to GenBank / GenPept. Save updates a supported source format, Save As chooses a destination format, and Export writes a copy without changing the open document's path.

| Format | Support and limits |
| --- | --- |
| SnapGene `.dna` | DNA, supported features, primers and notes. An unchanged same-format save preserves original bytes. Edited saves rebuild supported blocks and omit unsupported history, settings and other records after warning. `.rna` and `.prot` containers are unsupported. |
| GenBank / GenPept | Single-record sequence, supported locations, annotations and qualifiers. Some colour, direction and segment-note extensions are interpreted. Fuzzy, remote and complex locations are not fully supported. |
| EMBL | Single-record sequence and supported annotation fields. Conversion can lose source-specific metadata. |
| FASTA | Single sequence with header; annotations are not representable. Protein and RNA extensions preserve molecule intent. Multi-record input is unsupported. |
| ABIF `.ab1` / `.abi` | Read-only called bases, qualities, traces and selected instrument parameters. Unchanged source bytes can be copied. No trace editing or trimming. |
| Legacy `.popgene` | Read-only migration import. Save As to a supported interoperable format; no new private containers or sidecars are written. |

Format conversion warnings describe known losses. Original import metadata remains available in History, and original source bytes are retained in memory for the current session. This is not persistent backup or autosave. Keep original files when editing data with unsupported records.

Editing annotations marks them for review and removes stale cached translations. Calculated CDS translations use current bases, strand, ordered segments, `codon_start` and supported genetic codes; unsupported codes and translation exceptions are reported. Six-frame tracks use the standard code. Circular-seam translation and ORF discovery are not implemented.

Primer inspection shows recorded oligos and imported melting temperatures when available; it does not calculate validated reaction-specific Tm or design primers. Protein physical estimates are withheld for ambiguous or stop-containing input. Restriction analysis uses a limited catalogue and does not model methylation or digestion conditions.

SnapGene is referenced to identify a supported file format. PopGene is an independent project; no third-party application binaries, screenshots, curated sequence libraries or application source code are included.
