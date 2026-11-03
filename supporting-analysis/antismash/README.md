# antiSMASH supporting evidence for Figure 4

The [version 1.0.0 archive](https://doi.org/10.5281/zenodo.23219657) includes the native-result archive, its verification kit, and this support module. Original project contributions are covered by the repository [code](../../LICENSE) and [data](../../DATA_LICENSE.md) licenses, subject to [third-party terms](../../docs/licensing-scope.md). The native results and metadata were prepared in September 2026.

## Final-run mapping and inputs

`final-runs.tsv` names exactly 40 final outputs: 38 from the original `results/` folders and two from `rerun_two/out_GCA_029835695.1` and `out_GCA_037198095.1`. The two initial `results/` folders are not final runs. All 40 JSONs report antiSMASH 8.0.4, bacteria, and relaxed strictness, and contain KnownClusterBlast, SubClusterBlast, ClusterCompare, and active-site-finder evidence.

All 40 assembly FASTAs are already supplied in `../../analyses/figure4/reconstruction/inputs/assemblies/*.fna.gz`. They decompress byte-identically to the A3 panel inputs. The two rerun FASTAs also match their panel inputs. `packaged_input` paths in the TSV are relative to this folder. Do not substitute newer accessions without documenting an input change.

The two GFF3s are preserved under `historical/rerun_two/`. Their headers identify Pyrodigal 3.7.1, Single mode, and 5,965/5,865 CDS. Their generation script was not located in this audit. This is an explicit provenance limit, not a request to silently regenerate them.

## Companion archive location

Download `final-native-results.tar.gz` and the verification kit from the [antiSMASH evidence release](https://github.com/pushkarsai-dev/delftia-comparative-genomics/releases/tag/antismash-evidence-2026-09-30). The release also contains checksums. This module's `archive-descriptor.json` gives the byte size and SHA-256 of the native archive.

From this repository root, verify the companion without extracting it:

```sh
python3 supporting-analysis/antismash/verify_archive.py --archive /path/to/final-native-results.tar.gz
```

This verifies the archive digest, all 875 native-file hashes, and the exact
member set. All small metadata, saved scripts, GFF3 files, and logs are present
in this folder. Only `packaged_input` paths in `final-runs.tsv` resolve from this folder.
The `*_source` columns preserve paths relative to the original A3 directory,
not this repository. Packaged scripts and logs are under `historical/`
(including `historical/rerun_two/` for the two replacement runs). Map native
result sources to archive members with `native-output-manifest.json`.

## Native results and archive verification

`final-native-results.tar.gz` contains 875 final native output files: per-genome JSON and full GenBank, every region GenBank, and KnownClusterBlast/SubClusterBlast text. `native-output-manifest.json` records original relative paths, archive member paths, byte sizes, and SHA-256. Every archive member was read back and hash-verified. No archive extraction is needed for the already packaged Figure 4.

HTML/JS/CSS, images, duplicate ZIP exports, failed outputs, installation-repair scripts, and strictness trials are excluded. This is the complete selected native evidence set, not a byte-for-byte copy of every antiSMASH output-directory asset. Original directories remain intact.

The selected region and KnownClusterBlast text for every genome match the Figure 4 module byte-for-byte. BGC0000984.4 is the top match in at least one region in each genome; it is not the top match in every biosynthetic region. BAA-554 has 20 distinct matching query proteins in 22 alignments, cumulative score 38,447. All 45 reference translations agree after removing terminal stop markers from the query FASTA.

## Selection rule and panel (b)

`region-selection-audit.tsv` enumerates 54 qualifying regions across 40 genomes. Ten genomes have multiple candidates. The historical extractor returned the first qualifying result from an unsorted glob. This candidate preserves the exact selections recorded in the module's `region_manifest.json`; it does not reinterpret them as longest or highest-scoring choices. Every unselected alternative contains at most 24 translated CDS, below the 43-hit display threshold. Thus the 20 displayed tracks are unaffected by those alternatives. No new region search was run.

## Scripts, logs and environment

`historical/05_run_panel40_conda.sh` and `historical/rerun_two/run_rerun_two.sh` are the actual saved recipes. `04_setup_antismash_conda.sh` is installation provenance. These scripts are preserved, not executed or promoted as a fresh validated installer. Keep the GFF3 replacements distinct from the other 38 genomes. Do not infer a single Prodigal mode for all 40.

Thirty-six selected final logs are empty (34 original and both replacement logs). Other logs are limited, and the CPU count was not captured. Script options and result metadata support the settings; empty logs do not independently record the commands. The initial failure logs are preserved separately under `historical/superseded-initial-logs/`.

`historical/antismash_region_counts.tsv` is a later edited summary, not the runner's untouched stdout. Its `.bak` predecessor is also retained. The original A3 README describes an obsolete Docker route and is not used as the run recipe.

The installed environment still exists. `environment-observed-20260929.json`, its explicit-package text, and `conda-history.txt` record 134 packages and a creation entry dated 2 August 2026. This is a current observation, not independently certified historical state. No fresh conda installation was tested. `database-observed-20260929.json` gives hashes for six installed MIBiG 4.0 database files; contents are not redistributed. Other historical antiSMASH database versions and exact historical bytes remain unverified.

A future full antiSMASH reproduction would require an approved fresh environment, lawful exact database restoration, new output directories, the original 40 inputs, and two GFF3s, followed by comparison with this archive. It is not needed to inspect the deposited Figure 4 evidence and was not performed here.

## Rights and integration at the 30 September deposit

The author approved private GitHub release-asset storage on 30 September 2026. The native archive is an asset, not an ordinary Git object. This tag covers the Figure 4 review branch snapshot and does not integrate other figure branches. Public access, data/code/third-party rights, licensing, and a public archive identifier remain separate decisions. No licence has been assigned.

The archived setup/run scripts received wording-only troubleshooting-message cleanup on 1 October 2026. Executed commands and scientific parameters are unchanged; unmodified originals are preserved in the local audit archive.
