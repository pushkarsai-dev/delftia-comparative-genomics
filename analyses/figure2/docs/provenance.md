# Inputs, provenance, and calculation definitions

> **Status update — 7 October 2026:** Native pangenome and species-tree objects, and available service provenance, are now preserved in the [service-evidence release](https://github.com/pushkarsai-dev/delftia-comparative-genomics/releases/tag/service-evidence-2026-10-04). The accumulation report `266589/260/1` remains represented by recovered screen-derived provenance rather than a native report export. Approved component licenses are described in the [licensing scope](../../../docs/licensing-scope.md). The dated record below retains the original calculation and provenance boundaries.

## Source mapping

| Packaged input | Archived source | Use |
|---|---|---|
| `pangenome.json.gz` | A5 Pangenome Analysis: `Delftia_panel_40_pangenome.json` | Complete ortholog-family membership export; deterministic gzip compression only |
| `genome_name_map.tsv` | A5: same basename | KBase genome-reference to panel identifier mapping |
| `families.tsv` | A5: same basename | Saved family summaries, independently checked against the export |
| `presence_absence.tsv` | A5: same basename | Saved copy-number matrix, independently checked cell by cell |
| `node_table.tsv` | A8 Pangenome Accumulation: same basename | Historical clade membership and counts |
| `accumulation_curve.tsv` | A8: same basename | Historical rounded aggregates, checked against recomputed means |
| `tree_40tips_historical.json` | A8/raw: `speciestree_v2_full.json` | Historical tree used for accumulation, preserved without modification |
| `pan_phylo_table.tsv` | A8/raw: same basename | Exported KBase clade counts |
| `tree_40tips_corrected.newick` | Figure 1 module: same basename | Confirm identical clade membership, not a replacement accumulation input |

`input_checksums.json` pins all nine input files. The uncompressed pangenome JSON
SHA-256 is `6061bd966a8a0465bf884f2a18b0fbedfd58b71edb332c8adb2bf04c8b1fb738`.
No assignments were filtered, renamed, or reclustered during packaging.

## Record crosswalk

The BAA-554 Comparative Genomics Methods & Data Record §38 records the
5 August 2026 accumulation run, its strict requirement that all tree genomes
be in the pangenome, and report object `266589/260/1`. It names Phylogenetic
Pangenome Accumulation v1.7.2. Section 109 records the subsequent accumulation
paragraph verification. Section 161 records the Figure 1 tree correction;
Figure 2 uses clade membership only. The supplied archived objects and tables
are not a claim that complete raw KBase execution provenance has been deposited.

## Definitions

For each family, genome occupancy is the number of distinct genome references;
gene count includes multiple member genes in the same genome. A repeated
gene/genome identifier is rejected. The saved matrix records copy numbers,
despite its historical `presence_absence.tsv` filename.

For each clade, total families are the union, perfect core families occur in
every member genome, singletons occur in one member genome, and partial families
occur in more than one but fewer than all members. The 36 clades contain at
least two genomes. Means are arithmetic means among clades of equal size;
the plot uses exact means, while the historical aggregate table is rounded.

For this 40-genome panel, categories are core = 40, soft core = 38–39,
shell = 6–37, and cloud = 1–5 genomes. All original families remain included,
including possible clustering artifacts discussed in the manuscript.

## Reproduction boundary

This module rebuilds counts and figure outputs from the archived family object.
It does not execute KBase, infer orthology, retrain gene prediction, or infer a
tree. Correct postprocessing does not validate the biological identity of every
family or establish that sampling is unbiased. Licensing, deposition approval,
and full upstream provenance remain separate release tasks.
