# Public format reference

The Rust SnapGene reader was independently implemented using public TLV/XML descriptions and Biopython 1.88's `Bio/SeqIO/SnapGeneIO.py` as format references. That file credits Copyright 2017–2019 Damien Goutte-Gattat and is dual licensed under the Biopython License Agreement or BSD 3-Clause License. The upstream licence text is retained here as attribution. No Biopython code or runtime is linked into PopGene.

Source: https://github.com/biopython/biopython/blob/biopython-188/Bio/SeqIO/SnapGeneIO.py

ABIF field descriptions: https://incenp.org/dvlpt/docs/binary-sequence-formats/binary-sequence-formats.pdf

The compiled restriction catalogue is derived data from Biopython 1.88's
`Bio/Restriction/Restriction_Dictionary.py`, which credits Copyright (C) 2004
Frederic Sohm under the Biopython License Agreement retained in this directory.
That snapshot represents REBASE version 404 (2024), including recognition motifs,
nominal cleavage offsets, database identifiers and reference temperatures.
`validation/generate_restriction_catalogue.py` reproduces the Rust data without
a Python runtime dependency. EcoP15I is a separate factual vendor entry.

Catalogue source: https://github.com/biopython/biopython/blob/biopython-188/Bio/Restriction/Restriction_Dictionary.py

EcoP15I source: https://www.neb.com/en/products/r0646-ecop15i
