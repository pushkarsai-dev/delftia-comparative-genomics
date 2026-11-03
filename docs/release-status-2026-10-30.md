# Historical release status and provenance, 30 October 2026

This page describes the repository and draft on 30 October. For version 1.0.0, see the [3 November release notes](release-status-2026-11-03.md).

Checked 30 October 2026. This dated update supersedes earlier current-status statements about held reads, source packaging, and manuscript checks; original scientific evidence and historical records retain their dates.

## Access and archive

The GitHub repository remains private. The unpublished Zenodo draft reserves DOI `10.5281/zenodo.23219657`; it is restricted and is not a public deposit. Its four uploaded files match the local archive, inventory, README, and checksum file. The archive contains the 7 October main commit `1e5912f3fa54f7a0e61300cb7b12167b66e6e253`, 390 source files, and all 44 named evidence assets from eight prereleases. All archived source bytes match that commit, and all evidence hashes and sizes match the live release index.

This documentation update follows that archive. Before publication, refresh the source snapshot and affected manifests, verify the consolidated package, and replace the draft upload as needed. Record the resulting approved commit, version, access arrangements, and publication date. A reserved DOI and uploaded files do not establish public availability.

## Public sequences

NCBI publicly returns all 11 study runs under [SRP726995](https://www.ncbi.nlm.nih.gov/sra/?term=SRP726995). All 11 public normalized SRA data URLs returned HTTP 200 and the expected file size on 30 October. This checked access without downloading and rehashing the complete read datasets.

| Strain and public BioSample | Nanopore runs | Illumina run |
|---|---|---|
| BAA-554 [SAMN62386786](https://www.ncbi.nlm.nih.gov/biosample/SAMN62386786) | SRR40157013, SRR40157014, SRR40157015, SRR40157017 | SRR40157018 |
| ATCC 15668 [SAMN62386788](https://www.ncbi.nlm.nih.gov/biosample/SAMN62386788) | SRR40157008, SRR40157016 | SRR40157009 |
| ATCC 13751 [SAMN62386787](https://www.ncbi.nlm.nih.gov/biosample/SAMN62386787) | SRR40157010, SRR40157011 | SRR40157012 |

The ATCC 15668 chromosome [CM192903.1](https://www.ncbi.nlm.nih.gov/nuccore/CM192903.1) is public and is 6,796,693 bp. BAA-554 submission SUB16410709 and ATCC 13751 submission SUB16410659 remain Processing in the submission portal. The assigned ATCC 13751 master accession JCERJF000000000 does not yet resolve publicly. Public reads and BioSamples do not establish public access to those assemblies. The empty ENA study-file report does not negate verified NCBI read access.

## Species tree records

KBase support confirmed on 12 October 2026 that the intermediate alignment files are not retained. The available original evidence remains the saved app-cell outputs, logs, native tree objects, and pinned source material. No original concatenated 49-COG alignment was recovered, and no rerun was requested.

For the logged image tag `dockerhub-prod.kbase.us/kbase:speciestreebuilder.dce166f6d1673018a001b750c191b9a2deda0c71`, support returned digest `sha256:37d9858430a74b4fe58207e2ce367c67df081b5578e0cd06fcac880770c29ff1` from an inspection of that image. Support could not distinguish its original build from possible later rebuilds, although they believed it was unchanged. This digest is therefore a returned image reference, not verified proof of the original job's execution image. Neither the precise historical FastTree version nor its binary SHA-256 was supplied. Retain the wording FastTree 2 without an unverified point release.

## Validation and remaining work

Freshly downloaded source at the commit above passed all 231 tests, all source manifests, and all 139 relative documentation links. The five figure workflows and their deposited numerical inputs passed the 30 October technical checks. This does not substitute for final author approval or completion of the manuscript's submission templates and declarations. Final wording and references must be checked again after author edits.

Items #1 and #3–26 remain privately preserved within their documented scope. Item #2 remains an unavailable external SeqCoast QC/processing-record gap; item #27 retains author, access, final-proof, and release requirements; item #28 is optional historical material. KBase's response closes the records inquiry with a documented limitation, rather than recovering the missing alignment or historical runtime. No new scientific inference is implied.
