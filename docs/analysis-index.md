# Analysis index

The five figure modules in this release were integrated on 6 October 2026. The original integration identifier `1fdfbb38e0ab638415d87d9783ced2d8fddef47b` is preserved in historical records. Reproduce the released work from the [v1.0.0 source](https://github.com/pushkarsai-dev/delftia-comparative-genomics/releases/tag/v1.0.0).

- [Figure 1: Phylogeny and directional ANI](../analyses/figure1/README.md).
- [Figure 2: Pangenome structure](../analyses/figure2/README.md).
- [Figure 3: Detected protein families](../analyses/figure3/README.md).
- [Figure 4: Delftibactin-associated locus](../analyses/figure4/README.md).
- [Figure 5: Protein-domain assignments](../analyses/figure5/README.md).
- [antiSMASH native-output support](../supporting-analysis/antismash/README.md).

The [v1.0.0 release notes](release-status-2026-11-03.md) link the archive and public sequence records. See the [30 October provenance record](release-status-2026-10-30.md#species-tree-records) for the KBase response.

## Analysis packages

Start with the readable [evidence download table](evidence-downloads.md). The [release asset index](release-assets.json) gives the release URLs, asset names, byte sizes, and GitHub SHA-256 digests for the eight tagged evidence prereleases. Large data remain release assets rather than repository duplicates. Use named evidence assets; automatic source archives are the tagged source snapshots.

| Register | Package or module | Current preservation scope |
|---|---|---|
| 1 | Sequencing support | Run reports, assignments, read support, and input identities |
| 2 | SeqCoast processing records | EXTERNAL GAP: order-specific QC and processing records unavailable; no further author-account search planned |
| 3, 13 | 18-mer and locus crosswalk | Recovered computation, selected reads, and supporting annotation crosswalk |
| 4–7 | Assembly, graph, A6b, and FCS-GX service packages | Preserved native records and accepted evidence with explicit historical limitations |
| 8–10 | Panel selection and ANI | Selection records, type comparisons, and complete saved 40-by-40 comparison |
| 11 | Prediction evidence | Shared 239,999-protein FASTA and exact normalized Figure 4 sequence crosswalk |
| 12, 15, 17–18 | KBase service evidence and Figure 1 | Native objects, filtered assemblies, and available execution provenance |
| 14 | Figure 5 module | Two-protein InterProScan output and figure reconstruction |
| 16 | Keyword-screen evidence release | Recovered August execution snippets/output, native family input, and September four-combination records; historical runtime pins remain limited |
| 19 | antiSMASH release and support module | Selected native results and archive verification |
| 20–21 | Family-search evidence | Seventeen saved models, family results, Table S1 source, and pooled control |
| 22–23 | Direct-search evidence | Queries, saved comparisons, Table S3 source correction, and targeted UbiD/UbiX evidence |
| 24 | Naming survey and dated negative-query supplement | September 9 survey plus recovered September 8 original query code, returned header values, NCBI ESearch bodies, and accession summary; verified against downloaded evidence |
| 25–26 | Figure 4 module | Corrected reciprocal-hit verification and region-limited gene order |
| 27 | Figures and supplementary sources | Figures 1–5, Table 1 and the [current Table S2 and sequence identifiers](../supporting-analysis/manuscript-tables/README.md), and deposited S1/S3 sources; 8 October consistency review and 30 October technical checks completed; final submission templates, author approval, and verification of the final edited export remain open |
| 28 | Optional diagnostic history | Local; no dedicated deposit or current required claim identified |

See [register-to-deposit mapping](register-deposits.json) for machine-readable package locations. Package preservation does not recover every historical executable, cache, or service environment. Existing READMEs state those limits; no full study rerun is implied.

## Coverage and limits

The release covers register items #1 and #3–26. The original order-specific SeqCoast processing records (#2) remain unavailable. Item #27 includes final manuscript proof, author approval, and public access to the BAA-554 and ATCC 13751 assemblies. Optional historical diagnostics (#28) have no dedicated deposit. See the [release checklist](release-checklist.md) for the distinction between the code and evidence release and the remaining manuscript checks.
