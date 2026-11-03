# Delftia comparative genomics

Code and supporting evidence for comparative analyses of 40 *Delftia* genomes, including phylogeny, directional average nucleotide identity (ANI), pangenome structure, protein-family searches, and a delftibactin-associated locus. The accepted outputs for manuscript Figures 1–5 are included below.

**Release v1.0.0 (3 November 2026).** Download the [complete code and evidence archive](https://doi.org/10.5281/zenodo.23219657) (DOI `10.5281/zenodo.23219657`) or browse the [tagged source](https://github.com/pushkarsai-dev/delftia-comparative-genomics/releases/tag/v1.0.0). The Zenodo archive includes the exact tagged source and all 44 named evidence files. To download one analysis package, use the [evidence table](docs/evidence-downloads.md).

The [release notes](docs/release-status-2026-11-03.md) give sequence accessions and describe the limits of the archived evidence.

## Start here

| I want to… | Start with |
|---|---|
| View a manuscript figure | The [figure gallery](#figure-gallery) below; open the PDF for the vector version |
| Reproduce a figure | [Reproduction guide](docs/reproduce.md), then that figure's module README |
| Download supporting data or analysis evidence | [Evidence downloads](docs/evidence-downloads.md), with named archives and checksum instructions |
| Check which analyses are preserved and what is missing | [Analysis index](docs/analysis-index.md) and [release checklist](docs/release-checklist.md) |
| Cite or reuse the work | [Citation metadata](CITATION.cff), [code license](LICENSE), [data license](DATA_LICENSE.md), and [third-party terms](docs/licensing-scope.md) |

## Figure gallery

| Figure | Reference output | Reproduction instructions |
|---|---|---|
| Figure 1: Phylogeny and directional ANI | [PNG](analyses/figure1/reference/figure1.png) · [PDF](analyses/figure1/reference/figure1.pdf) | [Module](analyses/figure1/README.md) |
| Figure 2: Pangenome structure | [PNG](analyses/figure2/reference/figure2.png) · [PDF](analyses/figure2/reference/figure2.pdf) | [Module](analyses/figure2/README.md) |
| Figure 3: Detected protein families | [PNG](analyses/figure3/reference/figure3.png) · [PDF](analyses/figure3/reference/figure3.pdf) | [Module](analyses/figure3/README.md) |
| Figure 4: Delftibactin-associated locus | [PNG](analyses/figure4/reference/figure4.png) · [PDF](analyses/figure4/reference/figure4.pdf) | [Module](analyses/figure4/README.md) |
| Figure 5: Protein-domain assignments | [PNG](analyses/figure5/reference/figure5.png) · [PDF](analyses/figure5/reference/figure5.pdf) | [Module](analyses/figure5/README.md) |

## Reproduction scope

Each figure module supplies its own inputs, dependency instructions, commands, expected outputs, and tests. Use a separate environment for each module. Figure reproduction from saved results is distinct from repeating upstream annotation, clustering, tree inference, or sequence searches; each module describes which operations it implements.

Additional analysis code and evidence are preserved in versioned release assets. A source-code ZIP contains the repository files, but not these separately attached assets. Use the [evidence download table](docs/evidence-downloads.md) to retrieve the packages needed for a particular analysis. Raw read and instrument-signal deliveries are not bundled here; their accession and assignment records are identified in the supporting packages and manuscript Data Summary.

The [manuscript table files](supporting-analysis/manuscript-tables/README.md) include the current Table S2, its source attribution, and the assembly and BioSample identifiers of the 40-genome panel.

## Evidence status

The release preserves 25 register items (#1 and #3–26) within their documented scope. The original batch-specific SeqCoast QC and processing records (#2), original species-tree alignment, and exact historical FastTree point version were unavailable.

The [analysis index](docs/analysis-index.md) maps the preserved evidence to each item. Dated records inside the evidence packages describe their status when prepared.

## Citation and licensing

Sai, P., & Goller, C. C. (2026). *Delftia comparative genomics: analysis code and supporting evidence* (Version 1.0.0) [Software and supporting data]. Zenodo. [10.5281/zenodo.23219657](https://doi.org/10.5281/zenodo.23219657)

Use [CITATION.cff](CITATION.cff) for citation metadata. [MIT](LICENSE) applies to original project code and associated documentation. [CC BY 4.0](DATA_LICENSE.md) applies to original study data where the creators hold the rights. [Third-party terms](docs/licensing-scope.md) apply to material from other sources.

## Repository history and contributions

The public Git history starts with the v1.0.0 source snapshot. The eight evidence tags group the named downloads by analysis and preparation date. Their automatic source downloads point to the v1.0.0 snapshot. The [source map](docs/integration-sources.json) and [changelog](CHANGELOG.md) preserve the original review identifiers and dates as provenance. Use the public v1.0.0 tag to reproduce analyses. See [CONTRIBUTING.md](CONTRIBUTING.md) for future changes.
