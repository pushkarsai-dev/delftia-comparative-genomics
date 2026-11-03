# Figure 1 caption, alternative text, and methods

## Caption

```
Fig. 1. Phylogeny and directional average nucleotide identity of the 40-genome Delftia panel. (a) Species tree of the 40 panel genomes, obtained by removing the single reference genome from the inferred tree. Branch lengths are substitutions per site; the scale bar is 0.004 substitutions per site. Zero-length branches and multifurcations are shown as inferred. Filled circles, local support ≥ 0.95; open circles, 0.70 ≤ support < 0.95. Internal nodes without a circle either have support < 0.70 or lack a reported support value; the root also has no support value. (b) Directional ANI. Rows are query genomes, columns are reference genomes; ANI is not symmetric and the matrix is not averaged. The color scale runs from 95 to 100%, and no value falls outside it.

Rows and columns follow the same genome order. Numbers beside the tree identify the corresponding heatmap rows and columns; a tick marks every column and selected column numbers are labeled for readability. These are figure-local numbers, not strain or accession identifiers. Squares and blue type, D. tsuruhatensis; diamonds and brown type, D. acidovorans; red star, genome sequenced or reassembled in this study. Symbols mark documented culture-collection deposit groups, which are not an exhaustive audit of deposit equivalences across the panel: † rows 38, 39, and 40 are three deposits of one strain; ‡ rows 30 and 31 are two deposits of one strain; § rows 33 and 34 are two deposits of one strain. Each assembly is plotted as a separate tip.
```

## Alternative text

```
Phylogenetic tree of 40 Delftia genomes beside a 40 by 40 directional ANI heatmap sharing the same numbered rows, shaded 95 to 100 percent. Stars mark the four genomes sequenced or reassembled in this study; circles mark branch support.
```

## Phylogenetic methods

```
A species tree was inferred with the KBase Build Microbial SpeciesTree v1.6.0 app (Arkin et al., 2018) from the 40 panel genomes together with one nearest reference genome, Delftia acidovorans SPH-1 (GCF_000018665.1), and without a reference skeleton set. The app maps proteomes to 49 curated universal COG profile alignments with RPS-BLAST (Camacho et al., 2009), retains alignment columns with at least 95% occupancy, and infers a midpoint-rooted tree with FastTree 2 under its default JTT+CAT model using -fastest (Price et al., 2010). For Fig. 1a, SPH-1 was removed from the exported 41-tip tree using distance-preserving pruning, retaining the original root and the support labels of retained clades, with unreported support kept missing. No branch lengths were re-estimated. SH-like local support values are shown in Fig. 1a; ten internal nodes carry no support value in the untrimmed tree and are shown unmarked.
```

