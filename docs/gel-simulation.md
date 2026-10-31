# Agarose gel simulation

After **Analyze → Restriction digestion…** calculates products, the dialog draws a
simulated agarose gel: a ladder lane followed by the current digest. Choose
**Keep lane for comparison** to retain a digest, change the enzymes or conditions,
and simulate again; kept lanes stay on the gel until cleared. Up to ten digests can
be kept beside the ladder and the current result. **Copy gel report** copies every
lane's bands, run positions, relative intensities and notes as tab-separated text.

## Model

| Input | Source |
| --- | --- |
| Separation range | Efficient separation of linear dsDNA in standard agarose, 0.5–1.5%: Ausubel et al. (1995) *Current Protocols in Molecular Biology*, as [tabulated by Thermo Fisher](https://www.thermofisher.cn/cn/en/home/references/ambion-tech-support/rna-electrophoresis-markers/general-articles/gel-electrophoresis-tables.html). Concentrations between rows are interpolated geometrically. |
| GeneRuler 1 kb | Fourteen fragments and three reference bands from the [product manual](https://documents.thermofisher.com/TFS-Assets/LSG/manuals/MAN0013004_GeneRuler_1kb_DNALadder_250ug_UG.pdf). |
| λ DNA/HindIII | Eight fragments from the [product manual](https://documents.thermofisher.com/TFS-Assets/LSG/manuals/MAN0012986_Lambda_DNA_HindIII_Marker2_250ug_UG.pdf), including its note that unheated cohesive fragments can form an extra 27,491 bp band. |
| Custom ladder | Sizes entered in bp, drawn with uniform intensity. |

Band position is a logistic in log10(length) that places the largest and smallest
tabulated lengths at 15% and 85% of the run. Lengths outside the range are
compressed towards the well or the end of the gel and reported. Fragments closer
than 1.2% of the run are drawn as one band and reported as possibly unresolved.
Digest products are equimolar, so band intensity is proportional to fragment
length; the lane's strongest band is 100%. Fragment sizes are top-strand lengths,
as in the products list.

## Limits

Positions compare lanes on one gel; they are not measured mobilities. Voltage,
buffer, run time, stain, load and gel thickness are not modeled. Circular DNA
migrates according to its topology (supercoiled, nicked or relaxed), so an uncut
circular product is drawn with a dashed outline at its linear length for reference
only. Ladder band masses are not modeled beyond the manufacturer's reference bands.
Partial digests, star activity and agarose outside 0.5–1.5% are not simulated;
products marked provisional by the digestion model are labeled the same way on the gel.
