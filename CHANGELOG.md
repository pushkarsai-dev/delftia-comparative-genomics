# Changelog

## 1.0.0 (3 November 2026)

Tagged the integrated source for Figures 1–5 and antiSMASH support. The [versioned Zenodo archive](https://doi.org/10.5281/zenodo.23219657) combines this source and all 44 named supporting-evidence files. Updated citation metadata, download instructions, release status, and their manifests. The 18-mer/crosswalk package READMEs and accompanying checksum documentation were updated. All scientific inputs, computation code, and reference outputs retain their original bytes. The public Git history begins with this version. Earlier entries below document work in the original review repository, including identifiers from that history.

## 2026-10-07 - Reader navigation and status reconciliation

Added a short start page, a figure-reproduction guide, and a human-readable download table for all 17 evidence archives across eight private prereleases. Linked the checksum and verification instructions. Corrected outdated KBase recovery wording, added dated current-status notices to the Figure 2, 3, and 4 records and the antiSMASH README, and reconciled the release checklist with the 7 October technical review. Scientific code, input data, tests, accepted figures, and existing release-asset bytes are unchanged. Private access and preserved history remain in force; no public release or final DOI is claimed.

## 2026-10-07 - Approved credits, licenses, and release preparation

Recorded Pushkar Sai and Carlos C. Goller as creators in that order. Applied MIT to original project code and CC BY 4.0 to original study data, preserving third-party terms. Added citation metadata without unconfirmed identifiers or a release date. Public release remains explicitly on hold. Prepared a checksum-verified local delivery candidate, documented current sequence-access gates and the unavailable external SeqCoast record gap, and closed redundant Figure 2/3 PRs without deleting branches. Scientific module files and existing release assets are unchanged.

## 6 October 2026 - Main integration

Merged the reviewed integration into main at `1fdfbb38e0ab638415d87d9783ced2d8fddef47b`, preserving both parents and the review branch. Freshly downloaded merged source is byte-identical to the reviewed head; scientific module files and the applicability of 231 passing tests are unchanged. Original item 24 September 8 query evidence is privately deposited and fresh-download verified in dated-negative-query-evidence-2026-10-06, complementing the September 9 naming survey. Public access, licensing, and archival deposition remain separate.

## 5 October 2026

Prepared a local candidate combining exact Figure 1–5 modules and antiSMASH support. Reconciled shared navigation, all seven current evidence releases, register-to-package links, release boundaries, and the root checksum manifest. Prediction, search, and keyword-screen evidence are now privately deposited as prereleases; item 24 remains partial. The keyword package preserves recovered original August execution snippets/output and September checks, with limited historical runtime pins. No scientific module files changed, and no main merge or public release occurred.

## Preserved source history

Duplicate inherited entries are consolidated below. Distinct source entries and their dates are preserved. The current analysis index governs present deposit and integration status.

## 4 October 2026 - Figure 4 verification corrections

- Require raw-hit validation, an exact manifest decision grid, and keyed TSV/JSON field agreement before producing Figure 4c summaries.
- Add six rejection checks (including accepted contig identity), correct original/public profile hashes and metadata wording, and record original antiSMASH script identities.
- Label the September manuscript crosswalk as historical. Saved scientific calls, summary tables, and reference figures are unchanged.
- Reconcile the Figure 4 review base into the current-analysis branch without rewriting existing history.

## 2026-09-30 - Private antiSMASH evidence deposit

- Linked the versioned private GitHub asset for all 875 selected native files across 40 final runs.
- Added archive retrieval metadata; release assets include a portable verification kit and checksums.
- Kept the large archive outside Git history, preserved its existing bytes, and retained the separate public-deposition and licensing gates.

## 2026-09-29 - Figure 4 current-analysis evidence

- Identified the complete saved 40-genome execution supporting current Figure 4c reporting.
- Added independent aggregation of pair-level results, expected tables, and failure checks.
- Added the antiSMASH final-run map, scripts, GFF3 files, logs, archive manifest, and archive verifier; the native-output archive remains a separate companion file.
- Preserved scientific inputs, stored hits, earlier comparison records, and reference figures unchanged.
- Local candidate preparation does not imply a remote push, merge, public archive, licence selection, or sequence release.

## 2026-09-22 - Figure 4c full-panel local repeat

- Deposited a fresh 40-genome repeat of the reconstruction: 239,999 predicted proteins, 1,800 decisions, and 1,774 reciprocal matches.
- All selected protein pairs and decisions match; the matrix is byte-identical. Scores and some lower-ranked hit orders and domain coordinates differ and are documented.
- Preserved original detailed results, all 40 checkpoints, versions, timings, logs, and a read-only comparison tool.
- Updated current verification wording from a six-genome bounded repeat to a full-panel repeat, while preserving the earlier records and the historical identity limitation.
- Figure 4 inputs, reference PDF/PNG, caption, and search calculations are unchanged. No merge, licence decision, or DOI release is implied.

## 2026-09-22 - Figure 4c reconstruction integration

- Integrated saved forward/reverse evidence, input hashes, six repeat records,
  and the reviewed shared validator and runner. Included 40 compressed genome
  assemblies with a hash-verifying restoration command.
- Panel (c): independently corroborated by reconstruction (1,774 reciprocal
  best hits); exact historical cell/pair identity remains unestablished.
- Figure, caption, plotting data, and all Figure 1 files are unchanged.
- Validation-only clean-directory testing does not constitute a new 40-genome
  search. The previous 21 September pending status below is superseded.
- Prepared for the figure4-review branch; no upload, merge, or release made.

## 2026-09-22 - Figure 5 review module

- Added a self-contained archived-annotation plotting module under `analyses/figure5/`.
- Preserved both sequences, eight reported annotations and E-values, and the motif at residues 211-219.
- Displayed the archived discontinuous SSF102405 match as two segments rather than one continuous span.
- Enlarged the stacked panels and text, used a common residue scale, and embedded actual TrueType fonts.
- Added corruption tests, independent exported-PDF checks, pinned dependencies, and a relocated clean-directory repeat.
- Linked Figures 2-4 as separate review branches; no merge, visibility change, licence decision, or archival release.

## 2026-09-21 - Figure 3 review module

- Added `analyses/figure3/` on a separate, unmerged review branch from main.
- Added a taller 40-genome by 11-family heatmap, explicit zeros, larger labels, and vector PDF and 600-dpi PNG outputs, preserving all 440 historical displayed counts.
- Included independent-search evidence, a bounded upstream rerun script with external-input requirements, 27 tests, export validation, caption, and provenance.
- Verified the downloaded GitHub module against the prepared package and reran its tests and figure validation.
- Figure 1 is unchanged; Figure 2 remains in its separate review branch. No merge, public release, licensing decision, or DOI deposit was performed.

## 2026-09-21 - Figure 4 review contribution

- Added a self-contained Figure 4 module with 40 selected GenBank regions,
  source hashes, two reproducible parser modes, best-hit evidence, and tests.
- Reproduced all 1,293 historical region matches. Correct CDS-feature parsing
  retained every row and added 51, making GCA_037198095.1 eligible for panel (b).
- Corrected plot: 20 tracks, 933 source-checked gene arrows, readable heads,
  shared physical scale, 8-point labels, and vector PDF. Reference gene content
  and whole-proteome conservation summaries are unchanged.
- Fresh dependency installation and relocated execution passed 32 tests and
  export validation; fresh PDF, PNG, and plot-data files match exactly.
- Panel (c) upstream whole-proteome reproduction remains pending; no release,
  licence, publication, or whole-study reproduction is implied.
- This contribution starts from main db9fb53; Figures 2 and 3 are not merged.

## 2026-09-20 - Self-contained analysis modules

- Moved Figure 1 code, tests, data, reference outputs, dependencies, and documentation into `analyses/figure1/`.
- Added a Figure 1 README and prominent homepage links; updated contribution instructions and the module template.
- Updated input paths and checksum-manifest keys only; scientific input bytes, input digest values, and accepted reference outputs are unchanged.
- Removed superseded root-level Figure 1 locations; the previous layout remains recoverable in Git history.
- Retained pending statuses for ANI recomputation, other analyses, upstream KBase provenance, licensing, and release.

## 2026-09-20 - Figure 2 pangenome module

- Added the archived family object, saved tables, historical tree, and corrected-tree topology cross-check.
- Added full-panel count recomputation, plotting, export validation, negative controls, dependencies, and documentation in `analyses/figure2/`.
- Preserved the approved 170 × 105 mm two-panel figure and all Figure 1 module files.
- This contribution reproduces downstream family/clade counts and rendering, not gene prediction, OrthoMCL clustering, or phylogenetic inference.
- Licensing, public release, and an archival DOI remain pending.

## 2026-09-19 — initial Figure 1 contribution

- Added the tested Figure 1 rendering and tree-correction package, saved
  inputs, provenance, accepted reference outputs, and validation records.
- Added an incremental analysis index, module record, reusable documentation
  template, and contribution conventions.
- Preserved Figure 1 scientific code, inputs, tests, caption, and accepted
  outputs without modification during repository organization.
- No upstream ANI calculation or phylogenetic inference was rerun.
- Repository remains private; licensing, public release, and archival DOI
  are not established by this contribution.
