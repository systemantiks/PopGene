# Agarose gel simulation

**Actions → Simulate Agarose Gel…** opens a gel with the default ladder and the
current document loaded uncut. Up to twelve lanes can hold ladders or any open
double-stranded DNA document.

- **Lanes.** Add open sequences, catalogue ladders or a custom ladder. Click a lane
  label or band to select that lane, then change its label, enzymes (empty loads
  uncut DNA) or Dam/Dcm/CpG methylation assumptions and choose **Apply**.
  **Apply enzymes to all lanes** gives every sample lane the same digest. Lanes can
  be moved and removed, and a ladder can sit in any lane.
- **Gel.** Agarose percentage, 1X TAE or 1X TBE, relative run length (0.5–2×) and
  whether the calculated position marks the top, middle or bottom of a band.
  Fragments that run past the end of the gel are listed.
- **Files and figures.** **Save gel…** writes a `.popgel` file containing the
  settings and every loaded sequence; **Open gel…** restores it for further editing.
  The gel can be exported as SVG, PDF or PNG, copied as an image, printed, or copied
  as a tab-separated report that lists bands, run-off fragments and enzymes without
  a site. **Dark bands on white** inverts exported images; printing always uses it.
- **Defaults.** **Save as defaults** keeps the first catalogue ladder and the gel
  settings for new gels, in `gel-defaults.json` beside the recovery database.

**Actions → Restriction Digestion…** also draws a gel beneath its products, using
the same controls and defaults. **Keep lane for comparison** holds a digest while
the enzymes or conditions change. Clicking a band in the current digest lane selects
or clears its products for **Create selected products**. In **Actions → PCR & Cloning…**,
**Show on gel** opens the gel tool with the calculated product beside the default
ladder.

## Model

| Input | Source |
| --- | --- |
| Resolution range | Fragment size ranges for standard agarose in 1X TAE (0.6–2.0%) and 1X TBE (0.5–1.75%) from the [UltraPure Agarose protocol](https://www.thermofisher.com/jp/ja/home/references/protocols/nucleic-acid-purification-and-analysis/dna-protocol/ultrapure-agarose.html), also printed in the Invitrogen 1 Kb Plus manual. Concentrations between rows are interpolated geometrically. |
| GeneRuler 1 kb | [MAN0013004](https://documents.thermofisher.com/TFS-Assets/LSG/manuals/MAN0013004_GeneRuler_1kb_DNALadder_250ug_UG.pdf): 14 fragments; 6,000, 3,000 and 1,000 bp reference bands. |
| GeneRuler 1 kb Plus | [MAN0013047](https://documents.thermofisher.com/TFS-Assets/LSG/manuals/MAN0013047_GeneRuler_1kb_Plus_DNALadder_250ug_UG.pdf): 15 fragments; 5,000, 1,500 and 500 bp reference bands. |
| GeneRuler 50 bp | [MAN0013016](https://documents.thermofisher.com/TFS-Assets/LSG/manuals/MAN0013016_GeneRuler_50bp_DNALadder_50ug_UG.pdf): 13 fragments; 500 and 250 bp reference bands. |
| Invitrogen 1 Kb Plus | [MAN0000898 Rev. A.0](https://documents.thermofisher.com/TFS-Assets/LSG/manuals/1Kb_Plus_DNA_ladder_man.pdf): 18 fragments; 1,500 bp reference band. Earlier versions of this product had a different band set. |
| λ DNA/HindIII | [MAN0012986](https://documents.thermofisher.com/TFS-Assets/LSG/manuals/MAN0012986_Lambda_DNA_HindIII_Marker2_250ug_UG.pdf): 8 fragments, including the note that unheated cohesive fragments can form an extra 27,491 bp band. |
| Custom ladder | Sizes entered in bp, drawn with uniform intensity. |

At run length 1, band position is a logistic in log10(length) that places the
largest and smallest tabulated lengths at 15% and 85% of the gel. Lengths outside the
range are compressed towards the well or the end of the gel and reported. Run length
scales every distance; fragments beyond the end of the gel have run off. Fragments
closer than 1.2% of the gel are drawn as one band and reported as possibly
unresolved. Digest products are equimolar, so band intensity is proportional to
fragment length; the lane's strongest band is 100%. Sizes are top-strand lengths.
Sample lanes use the default digestion reaction (37 °C, pH 7.9, 10 mM Mg²⁺, 5%
glycerol) with the lane's methylation assumptions. The digestion dialog uses its own
entered conditions.

## Gel files

A `.popgel` file is UTF-8 JSON with `"format": "popgene-gel"` and `"version": 1`,
the gel settings and an ordered `lanes` list. Each lane is a catalogue `ladder` (by
name), a `custom_ladder` (label and sizes) or a `sample` (label, the loaded sequence
with its topology and any physical duplex ends, enzyme names and methylation
assumptions). Opening validates every lane before anything is shown. Unknown
ladders, enzymes or newer versions are refused rather than approximated. Files are
limited to 256 MiB.

## Limits

Positions compare lanes on one gel; they are not measured mobilities. Voltage, run
time in minutes, stain, load and gel thickness are not modeled. Circular DNA
migrates according to its topology (supercoiled, nicked or relaxed), so an uncut
circular molecule is drawn with a dashed outline at its linear length for reference
only. Ladder band masses are not modeled beyond the manufacturers' reference bands.
Partial digests and star activity are not simulated.

## Data still needed

These are not implemented because no citable data has been found yet. Each needs a
primary source, ideally the manufacturer's manual or a peer-reviewed paper, before
it is added.

| Item | What is needed | Where to look |
| --- | --- | --- |
| Supercoiled migration | Apparent linear size of supercoiled plasmids by agarose percentage and buffer (TAE, TBE, SB), so uncut plasmids can be placed. | Electrophoresis literature on supercoiled versus linear mobility; supercoiled ladder manuals (NEB N0472, ELPIS). |
| Sieving agarose, 2.5–4% | Resolution ranges by percentage for high-resolution agarose. | Lonza NuSieve 3:1 and MetaPhor datasheets; Thermo/Invitrogen E-Gel EX documentation. |
| NEB ladders | Band sizes and reference bands for 1 kb (N3232), 1 kb Plus (N3200), 100 bp (N3231), 2-Log (N3270/N0550), Low Molecular Weight (N3233), Quick-Load variants and the Supercoiled DNA Ladder (N0472). | Product pages and manuals at neb.com. The site refused automated retrieval, so these need to be copied from a browser. |
| Further Thermo/Invitrogen ladders | GeneRuler 100 bp, 100 bp Plus, Low Range, High Range and DNA Ladder Mix; O'GeneRuler and ZipRuler; Invitrogen TrackIt, 100 bp and Low DNA Mass. | Product information sheets linked from each catalogue page at thermofisher.com. |
| Other suppliers | DyneBio, PHENIX Research, ELPIS Biotech, Promega, Bioline HyperLadder, Takara, Bio-Rad and Qiagen ladders. | Each manufacturer's product manual or certificate of analysis. |

Adding a ladder means adding one entry to `LADDERS` in
`crates/popgene-core/src/gel.rs` with its sizes (largest first), reference bands,
whether it is an equimolar digest, a note and the source URL. The catalogue test
checks ordering, reference bands and unique names.
