# Figure 4 provenance and correction

> **Status update, 7 October 2026:** Figure 4 is integrated into `main`.
> Original-code MIT and original-study-data CC BY 4.0 terms are approved, with
> [third-party exclusions](../../../docs/licensing-scope.md) preserved. Earlier
> licensing holds below describe their dated preparation stages. Public release
> and final archival publication remain on hold; see the
> [current analysis index](../../../docs/analysis-index.md).

> Current reporting, 29 September 2026: the complete saved execution in
> `reconstruction/full_panel_repeat/results/` is the primary pair-level evidence
> for the manuscript. See the [current-analysis record](../CURRENT_ANALYSIS.md).
> The development and comparison history below is retained as an audit record.

## Historical inputs

- Methods and Data Record §50 (19 August 2026): reference region identification,
  45 proteins, and panel-wide reciprocal-hit summaries. The figure's reference
  FASTA, conservation CSV, and per-genome CSV are byte-preserved from A11.
- §54 (21 August): completion of the two antiSMASH reruns supplied with
  Pyrodigal-generated GFF3 features. All 40 genomes have a BGC0000984 top match.
- A14 Region Gene Order (24 August): region-restricted reciprocal searches,
  stored as `region_gene_order.tsv`, distinct from the whole-proteome search.
- A12 `fig4.py`: historical 19-genome selection and three-panel display.

GenBank region comments identify antiSMASH 8.0.4 and the source run. The supplied
manifest identifies the exact region selected for each genome. Some fragmented
assemblies have multiple regions with the same top ClusterBlast match; we retain
the record actually used in A14. A check of alternative candidates found none
with enough CDS to add another at-least-43 track beyond the corrected set.

## Independently checked on 21 September 2026

1. All 45 reference FASTA proteins match the GenBank translations and coordinates.
   The two historical locus tags also match their archived InterProScan query
   sequences exactly after removing terminal stop characters.
2. All 1,293 historical gene-coordinate rows agree with the corresponding CDS
   features. A fresh historical-parser search reproduced every row exactly.
3. The corrected feature-aware search retained all 1,293 rows and added 51.
   Forward and reverse best-hit identities, E-values, and bit scores are saved.
4. Panel-wide summary marginals both total 1,774 matches. Per-genome counts are
   43 in four genomes, 44 in 18, and 45 in 18; 41 genes have counts of 40, with
   four others at 25, 33, 37, and 39. These are consistency checks, not a fresh
   whole-proteome RBH run.
5. The first BAA-554 ClusterBlast block contains 20 unique matching query
   proteins in 22 alignments, summing to 38,447. Other cluster blocks are excluded.

## Extraction defect and scientific effect

The historical parser searched a block extending from one CDS to the next for
`locus_tag`. In the two GFF3-supplied reruns, CDS features used `gene`/`ID`
instead; a locus tag was sometimes found only in a following domain feature.
Consequently, CDS without such domain features were omitted from the target set.
The new parser restricts extraction to each CDS and accepts all three identifier
fields. This corrects parsing, not gene prediction, annotation, or the threshold.

GCA_029835695.1 has 30 translated CDS in the selected region, rather than the
eight extracted historically, and its reciprocal matches rise from seven to 22.
GCA_037198095.1 has 44, rather than eight, and its matches rise from eight to 44.
The latter becomes the twentieth plotted assembly. No previously plotted track
is removed or changes its match set. Both reference-to-target and reverse searches
use the complete corrected target set; results were not manually patched.

The corrected 20 tracks retain collinear gene order and relative strands.
The same three records are reversed for display: ATCC13751_hybrid,
GCA_029590315.1, and GCA_965137545.1. This is a display normalization, not
evidence of a biological inversion. Block spans are 77,927–84,738 bp,
or 77.9–84.7 kb at the manuscript precision. BAA-554 is 84,632 bp (84.6 kb).

## Coordinates and graphics

The reference region is zero-based [2159931, 2245340), equivalent to 1-based
inclusive 2,159,932–2,245,340 (85,409 bp). Archived CDS coordinates are 1-based
inclusive. The new plot uses interval edges start-1 and end, and computes span
as end-start+1. The historical plot omitted that final base from lengths;
correcting it changes none of the displayed one-decimal spans.

The old arrowhead length mixed vertical row units with genomic base-pair units,
producing effectively invisible heads. New polygon heads use genomic x units,
bounded by the actual CDS length. The whole CDS footprint, strand, and spacing
are preserved. The 10 kb scale shares the track x-axis. Short CDS necessarily
remain short; no gene is artificially lengthened. The page is taller and labels
are at least 8 pt. Gray ClusterBlast ticks avoid overlapping triangle markers.

## Limits and release holds

The original panel-wide per-protein reciprocal-hit files and executable were
not saved. The subsequent whole-proteome reconstruction independently
corroborates the counts, but cannot establish historical cell/pair identity.
All 45 per-gene and 40 per-genome totals match; these are matrix marginals,
not proof of the lost historical matrix. The earlier bounded repeat covered six
genomes; a separate [full-panel repeat](../reconstruction/full_panel_repeat/README.md)
now records all 40, with identical selected pairs and decisions but some
lower-level numerical and alignment differences.
All 40 input assemblies, saved search evidence, and reviewed checking code
are now packaged under `reconstruction/`. Four search choices are explicitly
documented as reconstruction settings, not recovered historical facts.
This does not reproduce antiSMASH prediction or establish metabolite production.

The full internal Methods Record and correspondence are not included. Public
release, licensing, ownership approval, and DOI deposition remain separate holds.

## 22 September packaging addendum

The 21 September summary checks above describe the earlier audit. The later
whole-proteome reconstruction is now integrated without changing any plotting
input, figure output, or caption. `reconstruction/README.md` identifies the
input identity limits, chosen settings, raw-hit checks, and six-genome repeat.
The package was tested from a clean, relocated directory. That packaging test
did not repeat the 40 sequence searches.

## 22 September full-panel execution addendum

Following the packaging test, the unchanged deposited runner was executed from
commit `2eeb444ea6dbc0dbe7d75477b9063bc40e154f78` in a fresh directory. All 40
proteomes were predicted anew, with no reused cache or checkpoints. The full
repeat yielded 1,774 reciprocal best hits and identical selected protein pairs,
decisions, and matrix bytes. `reconstruction/full_panel_repeat/` preserves the
new results and reports differences in scores, lower-ranked ordering, and
domain coordinates. No selected best-hit coordinate changed. The cause of
cross-platform numerical differences was not isolated experimentally.
This extends repeat evidence, not historical provenance, and changes no
plotting inputs, reference figure, caption, or search calculation.
