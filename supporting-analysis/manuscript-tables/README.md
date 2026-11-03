# Genome-panel accessions and Supplementary Table S2

The tab-separated table contains the 16 excluded assemblies and six columns. The accompanying caption and table match the 3 November 2026 manuscript. All numerical values and exclusion decisions are unchanged from the dated 4 October panel-selection evidence package. Use this table and caption in place of the earlier display copies in that package.

## Evidence for the GCA_903684635.1 row

The study BAA-554 comparison (90.38%) is recorded in the saved selection worksheet, preserved as `panel-selection-evidence/source-worksheet/05-exclusions.tsv` and its normalized extract `panel-selection-evidence/exclusions.tsv` (source row 19). The worksheet identifies the result but does not establish a particular historical executable or invocation.

The comparison with the separate public NBRC 16741 assembly (GCF_001571325.1) is 90.22% ANI. The match to GCF_903815225.1 is 100% ANI with alignment fraction 0.999. Both are FastANI results in the retained GTDB-Tk v2.3.2 output against GTDB R214. Within `panel-selection-evidence/native-qc/GTDB-Tk/output/`, `id_to_name.map` maps the excluded assembly to `id31`. The corresponding row of `gtdbtk.bac120.summary.tsv` lists GCF_001571325.1 under other related references and GCF_903815225.1 as the FastANI reference.

The analysis-register item #10 whole-panel results cover the 40 retained assemblies. Their 1,600 directional comparisons and fragment counts are in `whole-panel-ani-evidence/all-1600-comparisons.tsv` in the whole-panel ANI archive. The excluded assembly does not belong to that comparison matrix.

The panel-selection and whole-panel ANI evidence archives are available from the [panel and ANI evidence release](https://github.com/pushkarsai-dev/delftia-comparative-genomics/releases/tag/panel-ani-evidence-2026-10-04).

## Accession crosswalk for the retained panel

`panel-accessions.tsv` lists the 40 retained genomes in Table 1 order. It links the strain label in the manuscript to the BioSample, BioProject, and versioned public assembly accession or assembly DOI. For the 37 rows with a public GenBank assembly accession, the BioSample, BioProject, and submitting organization were retrieved from the [NCBI Datasets genome reports](https://www.ncbi.nlm.nih.gov/datasets/genome/) on 3 November 2026 using the versioned accession in that row. The original submitter is given as recorded by NCBI. “Not designated” is the Table 1 strain label for four records; distinguish them using the assembly and BioSample accessions.

For the three laboratory strains, the table identifies all 11 source read runs from [SRA study SRP726995](https://www.ncbi.nlm.nih.gov/sra/?term=SRP726995). As of 3 November 2026, the public study BioProject linked the ATCC 15668 assembly. No public INSDC assembly was verified for BAA-554 or ATCC 13751. Their GenBank submission identifiers and the assigned ATCC 13751 WGS project identifier are shown in separate columns. The assembly-service evidence archive contains the study FASTAs. The original read BioSample and BioProject are given for the separately published [ALG19 reassembly](https://doi.org/10.5281/zenodo.21651127) made from [SRR31881480](https://www.ebi.ac.uk/ena/browser/view/SRR31881480); the reads were produced by [Wang et al. (2025)](https://doi.org/10.3389/fmicb.2025.1567054).

The two accession columns distinguish an INSDC assembly identifier from a dataset DOI. Blank cells indicate an identifier unavailable on the observation date or a field outside the scope of the row. Run identifiers are given only for this study's three assemblies and the ALG19 reassembly. URLs lead to the public assembly, reassembly, or study evidence file corresponding to each row.
