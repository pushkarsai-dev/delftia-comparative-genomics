# Source and reproduction provenance

## Historical records

Methods and Data Record §48 describes the 19 August 2026 13-family scan; §49 describes the four additional family searches; §59.2 describes Figure 3. The archived plot reads `panel40_family_counts.csv` and adds the PF18631 zero row from `pictet_spenglerase_search.txt`.

`data/panel40_family_counts.csv` and `data/pictet_spenglerase_search.txt` are byte-preserved historical inputs. The search report is retained as evidence of the recorded PF18631 negative result, not as current scientific prose. Its broader conclusions and old reference identifications are superseded by the current manuscript; do not copy them into new captions or Results text.

The CSV contains 40 genomes, a protein total column, and 13 family-count columns. The figure selects ten of those families plus the separately searched PF18631 family. Genome order is preserved from the old figure: four study assemblies first, then other assembly identifiers in lexical order. No tree ordering, clustering, or species inference is introduced here.

## New independent execution, 21 September 2026

`data/independent_search.json` contains the input FASTA SHA-256 values, per-genome protein totals, extracted-profile hashes, tool versions, and complete counts under both target-database arrangements. `data/independent_hits.tsv` contains one row per protein-family hit, including the target arrangement, sequence E-value, and bit score. These are newly generated evidence, not recovered original output.

Prediction: Pyrodigal 3.7.1, a newly trained `GeneFinder(meta=False)` for each genome; all contigs in the genome were supplied for training. A new prediction was made from each FASTA; no previous proteome cache was loaded. The 40 FASTA hashes were checked against the 8 September deposition package's `01_proteomes_manifest.json`. All 40 counts and the 239,999 total matched.

Search: PyHMMER 0.12.2 `hmmsearch`, reporting threshold E=1e-3, retaining protein hits with sequence E-value strictly below 1e-3. Profiles were extracted from the local Pfam-A file with MD5 `2cacd81054f0c267718e752ff9322bf1`. The archived record identifies this file as Pfam 37.3. All 17 families were evaluated. Each protein counts once for a matching family, not once per domain.

**Target-database size matters.** Per-genome searches exactly reproduced all 520 original CSV cells and the separate report's PF18631/PF02317/PF10100 zero totals and PF02423 total of 198. Pooled searches altered 92 of the 520 archived cells, all in families omitted from Figure 3. The pooled PF02423 total was 85. None of the 440 displayed values changed between arrangements. This establishes the new reproduction's parameters; exact agreement alone does not recover undocumented original code.

## Optional upstream reproduction

From this module directory, with the original project archive and full Pfam-A database available:

```sh
python scripts/rerun_family_search.py \
  --project '/path/to/July 2026 MGen Paper Work' \
  --pfam '/path/to/Pfam-A.hmm' \
  --outdir '/path/to/new-search-output'
```

Requires `pyrodigal==3.7.1` and `pyhmmer==0.12.2` in addition to the plotting dependencies. It reads `A3 antiSMASH Local (2 Aug 2026)/panel40_fasta`, the A10 family-count CSV, and `A10 InterProScan (14 Aug 2026)/Deposition package (8 Sep 2026)/reference_run/results/01_proteomes_manifest.json` beneath `--project`. It refuses existing output directories and mismatched FASTA hashes or Pfam MD5. It creates new proteomes, extracted profiles, raw hits, counts, and mismatch lists; it never reuses or overwrites historical caches. The external archive layout is an explicit dependency, not a claim that the complete original genome archive is deposited in this module.

The submitted figure is family-level evidence. It does not independently establish the biochemical activities, substrate specificities, or pathway products of individual proteins. Identification of UbiX versus CoaC/DFP-like proteins relies on the separate sequence-level work described in the manuscript.
