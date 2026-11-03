# Release v1.0.0

**3 November 2026.** [Zenodo DOI 10.5281/zenodo.23219657](https://doi.org/10.5281/zenodo.23219657) identifies the source for [GitHub tag v1.0.0](https://github.com/pushkarsai-dev/delftia-comparative-genomics/releases/tag/v1.0.0) and the supporting evidence. Use the [evidence table](evidence-downloads.md) to find an individual package or the [reproduction guide](reproduce.md) to render a figure.

## Contents

The complete archive includes the tagged source, 44 named evidence assets, an inventory, and an offline checksum verifier. The eight evidence tags group the 17 analysis archives and their companion files by preparation date. All eight tags identify the same source snapshot as v1.0.0. Use the named asset and its digest to identify an evidence package. Original branch names and commit hashes in dated records refer to the earlier review history.

The 18-mer/crosswalk archive has two revised README files and an updated internal manifest. Its other 27 members match the original package byte for byte. The archive, its external checksum file, and its download instructions have new digests in the [release asset index](release-assets.json). The 41 other named evidence assets retain their original bytes.

The [manuscript tables](../supporting-analysis/manuscript-tables/README.md) supply the current Supplementary Table S2 and the 40-genome panel crosswalk. The dated panel-selection package preserves the original worksheet and native GTDB-Tk evidence. The current table identifies the comparator and measurement source for GCA_903684635.1; numerical values and exclusions are unchanged.

## Sequence records

The [study BioProject PRJNA1512546](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1512546) links the sequencing records. On 3 November 2026, [SRA study SRP726995](https://www.ncbi.nlm.nih.gov/sra/?term=SRP726995) returned 11 public runs across the three laboratory strains. The public assembly search under that BioProject returned [ATCC 15668, GCA_060658745.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_060658745.1/), with chromosome [CM192903.1](https://www.ncbi.nlm.nih.gov/nuccore/CM192903.1). Public assembly records for BAA-554 and ATCC 13751 had not been verified under the study BioProject on that date.

The [ALG19 reassembly](https://doi.org/10.5281/zenodo.21651127) is a separate deposit of the assembly made from [previously published reads, SRR31881480](https://www.ebi.ac.uk/ena/browser/view/SRR31881480). The sequencing-support package contains mappings and records. The consolidated archive excludes raw instrument signals and full read deliveries.

## Provenance and rights

The [analysis index](analysis-index.md) describes the preserved scope. Original batch-specific SeqCoast QC records and the historical 49-COG species-tree alignment remain unavailable. The [KBase provenance record](release-status-2026-10-30.md#species-tree-records) explains the limits on identifying the exact historical FastTree build.

Pushkar Sai and Carlos C. Goller are the creators. Original project code and usage documentation use MIT. Original study data use CC BY 4.0 where the creators hold the rights. [Component-specific rights and third-party attribution](licensing-scope.md) apply throughout the archive.
