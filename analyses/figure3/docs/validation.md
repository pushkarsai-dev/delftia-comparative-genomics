# Figure 3 validation, 21 September 2026

> **Current status, 7 October 2026:** the record below describes the 21 September
> validation snapshot. Figure 3 is now integrated into `main`, and related
> prediction/search evidence is privately deposited. The upload and merge
> statements below are historical, not pending work. See the
> [analysis index](../../../docs/analysis-index.md) and
> [evidence downloads](../../../docs/evidence-downloads.md). Public release remains on hold.

## Numerical reproduction

- Original panel FASTA hashes: 40/40 match the archived 8 September manifest.
- New gene prediction: all per-genome totals and the total of 239,999 proteins match.
- Per-genome HMM searches: 520/520 original CSV cells match.
- Separate families: PF18631, PF02317, and PF10100 remain zero in every genome; PF02423 totals 198, matching the separate historical report.
- Figure 3: all 440 old exported PDF numerals/dashes agree with the archived inputs. All 440 new exported PDF numerals agree with both the archived inputs and independently counted new hit records.
- Target-specific controls: PF00282 finds no panel hit and exactly one hit after A7B1V0 is added (225.2577 bits, E=1.3177474e-65). PF18631 finds no panel hit and detects all four reference proteins (68.8824–172.6368 bits). These were rerun with the newly predicted proteins, not taken on trust from an earlier narrative.
- Pooled versus per-genome searches: 92 archived CSV values differ under pooling, all in the omitted PF03088, PF20067, and PF00155 families. No displayed value differs. Database arrangement is documented explicitly; it must not be confused with the separate pooled positive-control experiment.

## Figure and package checks

The module includes 27 regression and failure tests. Tests cover changed counts, missing and duplicate genomes, missing PF18631 evidence, altered input hashes, a changed family count that preserves the total, wrong plot labels, shifted rows, incorrect colors, small text, and protected/existing output paths.

The export validator independently counts raw hit records, extracts and checks all 440 numerical PDF glyphs by cell position, checks dimensions and page boundaries, checks embedded fonts with Poppler, and compares the complete PNG and PDF rasterization against a fresh same-environment render. Exact pixels can vary between font installations or library versions; numerical checks are exact.

Design: 170 × 220 mm; 8-point cell numbers and labels; Arial TrueType on the tested Mac, with TrueType fallbacks for other environments. No text extends beyond the exported page. The historical plot had two out-of-page characters and embedded CFF/TrueType mismatch warnings; the new export has neither.

## Limits and release status

This verifies computation and presentation under the recorded parameters. It does not demonstrate biochemical activity, prove sensitivity to unrelated enzyme folds, or establish completeness of all possible biological pathways. The original historical scan code was not recovered; the independent scan script is a new reproduction, with that status retained in its output.

The package is prepared locally. GitHub upload, a review branch, downloaded-branch verification, and a pull request remain pending. Do not describe it as deposited, publicly released, or merged. Keep the later Figure 3 review branch unmerged until explicitly approved. Figures 1 and 2 and all original analysis files are unchanged.
