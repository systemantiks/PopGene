# Protein calculation models

The Properties tab uses the Bjellqvist pKa model. In **Analyze → Protein analysis**,
choose a model, pH, sequence range and explicit modifications before calculating.
The report identifies the selected assumptions. A selected interval is treated as
an excised linear peptide with new termini; no processing or disulfide bond is
inferred from sequence annotations.

| Model | Parameter source |
| --- | --- |
| Bjellqvist | Residue-specific terminal corrections and the [ExPASy method](https://web.expasy.org/compute_pi/pi_tool-doc.html). This remains the default. |
| Original EMBOSS | The published [Epk.dat values](https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/emboss/data/Epk.dat). |
| Grimsley 2009 (IPC Arg) | [Measured protein-group means](https://doi.org/10.1002/pro.19), completed with Arg pKa 12.04 as specified in [Kozlowski 2016, Table 4](https://doi.org/10.1186/s13062-016-0159-9). Arg was not measured in the former summary. |

The Grimsley option uses N-terminus 7.7, C-terminus 3.3, Cys 6.8, Asp 3.5,
Glu 4.2, His 6.6, Lys 10.5, Arg 12.04 and Tyr 10.3. Charge uses independent
ionizable groups and Henderson–Hasselbalch equations. A pI requires a bracketed
root: charge must be positive at pH 0 and negative at pH 14. A fully blocked
peptide with zero charge throughout that range has no unique reported pI.
These are estimates for the stated model,
not measured properties of a folded protein.

Changing the charge model does not change molecular mass. Average and
monoisotopic mass use the pinned Biopython IUPAC mass tables and explicit
modification deltas. Different mass tables can give small numerical differences.
The 42-residue comparison fixture, for example, gives 5029.7561 Da with this
table; another application's display gave 5029.83 Da. No undocumented mass
table is substituted to force a match.

The same fixture gives pI 7.6030815 and charge at pH 7 of 1.1785952 with the
Grimsley/IPC parameters. The independent checker is
`validation/check_protein_models.py`; GUI acceptance is recorded separately in
[development and validation](development.md).
