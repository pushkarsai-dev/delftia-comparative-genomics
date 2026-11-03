# Figure 4c: current analysis

The complete saved execution in `reconstruction/full_panel_repeat/results/`
is the primary pair-level evidence for the current manuscript's Figure 4c.
This selection changes the reporting and verification path. It does not change
any sequence, search parameter, stored call, plotted count, or reference figure.

The plotting CSVs in `data/region002_*.csv` are retained from the 19 August
analysis, whose executable and pair-level tables were not located; see the
[preserved provenance record](docs/provenance.md). The current derived tables
are computed independently from the selected complete run and agree with
the plotted counts. They do not establish recovery of the earlier execution.

## Execution and inputs

- Executed on 22 September 2026 with fresh protein predictions and no reused checkpoints.
- Search source commit: `2eeb444ea6dbc0dbe7d75477b9063bc40e154f78`.
- Saved evidence introduced in commit `59d19550bc1a170ae892cd6b2f4eaa6035ee75f5`.
- Update prepared from `figure4-review` commit `27e6c454135d5cc93ab9fc930cf839b663e995a0`.
- Forty assemblies, 45 verified queries, and 239,999 predicted proteins.
- Pyrodigal 3.7.1, with per-genome training; PyHMMER 0.12.2, with single-sequence profiles and `hmmsearch` in both directions.
- Sequence E < 1e-5 in both directions, reciprocal agreement, lowest E then highest bit score, exact ties excluded, and no additional coverage cutoff.
- Seed 42; actual target-protein counts define database sizes. Reverse searches use the complete BAA-554 proteome.

The saved `compatibility_manifest.json`, `proteome_manifest.json`,
`query_id_check.json`, and `environment.json` record input hashes, effective
settings, and software. `reconstruction/full_panel_repeat/execution.json`
records the source commit and execution. Scientific evidence and comparison
records remain unchanged; the public-copy provenance metadata were corrected on 4 October 2026.

## Verify before using the summaries

From `analyses/figure4`, in the pinned environment:

```sh
python reconstruction/bin/validate_evidence.py --results reconstruction/full_panel_repeat/results
python scripts/verify_current_analysis.py --module . --output build/current-analysis-check
```

Use a new or empty output directory. The second command runs the raw-hit
validator itself and stops on failure. It requires the exact 40 by 45 manifest
grid and exact agreement between every TSV field and its validated JSON call,
including rejected pairs. Accepted contigs must match the contig prefix in the validated target-protein ID, using the predictor's recorded `{contig}_{1-based ordinal}` convention. It then checks all 45 per-gene counts and all 40
per-genome recovery and contig counts against the plotting inputs, and writes
two aggregate CSVs plus a JSON report. The first command remains useful as a
standalone diagnostic; it is no longer an unchecked prerequisite.

Expected: 1,800 unique decisions, 1,774 accepted reciprocal matches, 41 genes
present in all 40 genomes, and variable counts of 25, 33, 37, and 39. Four
genomes have 43 matches, 18 have 44, and 18 have 45; the median is 44.
Twenty-three genomes have all recovered matches on one contig.

The saved report and derived tables are in `validation/current-analysis/`.
The derived per-gene CSV is byte-identical to its plotting input. The derived
per-genome CSV contains the three verified columns only; it does not copy
unverified columns from the original summary.

## Columns outside the current verification scope

`data/region002_rbh_matrix.csv` is preserved intact, including the older
`largest_block_genes` and `block_span_kb` columns. Those two columns are not
validated by this verifier and do not support a current manuscript claim.
Do not use them for Figure 4b spans: those come from the corrected region
coordinates in `data/region_gene_order_corrected.tsv`. The three supported
whole-proteome columns are explicit in the derived per-genome CSV.

The public copy of `reconstruction/evidence/profiles/queries.hmm` had its
45 COM command-history lines normalized on 1 October 2026 to remove local
paths. `reconstruction/evidence/source_files.json` records distinct original
and public-copy hashes and marks it as not byte-preserved. The selected
complete run's profile copy also omits local COM paths, as documented in
`reconstruction/full_panel_repeat/execution.json`. These metadata changes
do not alter profile parameters, saved calls, or raw-hit results. DATE/COM
fields are command-history metadata, not a substitute for checking results.

## Separate panels and supporting evidence

Figure 4b uses the separate region-restricted `phmmer` search, its corrected
20-track table, and the selected GenBank regions. It is not generated from
Figure 4c's whole-proteome calls. Figure 4a and the region selection are also
supported by the [antiSMASH evidence package](../../supporting-analysis/antismash/README.md).
The native-output archive is a separate companion file, with its digest and
verification command in that package; no public archive identifier is claimed.

The [existing comparison record](reconstruction/full_panel_repeat/README.md)
reports agreement in selected pairs and decisions and the observed numerical
differences. The causes of those differences were not isolated.

## Software references (APA)

Larralde, M. (2022). Pyrodigal: Python bindings and interface to Prodigal, an efficient method for gene prediction in prokaryotes. *Journal of Open Source Software, 7*(72), Article 4296. https://doi.org/10.21105/joss.04296

Larralde, M., & Zeller, G. (2023). PyHMMER: A Python library binding to HMMER for efficient sequence analysis. *Bioinformatics, 39*(5), Article btad214. https://doi.org/10.1093/bioinformatics/btad214

The dated `docs/manuscript-crosswalk-20260929.tsv` pins the 29 September manuscript. Its printed line numbers are historical locators. In the current manuscript, use the phrases “Cluster conservation was assessed” and “For the gene-order comparison” to locate the corresponding methods.
