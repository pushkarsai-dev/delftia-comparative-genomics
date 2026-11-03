# Verification record — 19 September 2026

## Module relocation verification — 20 September 2026

The self-contained `analyses/figure1/` layout was tested from a clean copy and
an unrelated working directory in the same environment recorded below.
All 29 tests passed. Input verification, tree reconstruction, figure generation,
and exported-figure validation passed, including all 780 retained-tip distances
and all 1,600 ANI cells. Nonempty and protected output destinations were refused.

Twelve scientific input/reference files retained their original bytes; input
checksum keys changed only to reflect the new paths. A fresh regeneration of
the previous layout in the same environment produced identical PNG pixels and
identical PDF pixels at 200 dpi. Mapping and support-audit CSVs were identical;
the mapping JSON differed only in its relocated source path. No scientific
values, tree inference, or figure layout changed.

See the [relocation test log](reorganization-tests.txt) and
[relocation figure-validation log](reorganization-validation.txt). This does
not constitute new tree inference or ANI recomputation. The dated checks below
are retained as the original package's verification record.

## Executed checks

The package was copied to a new temporary directory outside the project
archive and run from an unrelated working directory. No original project path
or authenticated service was needed during execution.

| Check | Result |
|---|---|
| Input integrity | 7 input files match pinned SHA-256 values |
| Regression suite | 29 tests passed, none skipped |
| Tree reconstruction | Byte-identical to the accepted corrected Newick |
| Retained-tip distances | 780/780 equal the 41-tip source |
| Support metadata | 3 filled, 17 open, 5 low-support, 10 unreported, and 1 unlabelled root |
| Plotted ANI | 1,600/1,600 cells equal independently indexed source values |
| Genome identities | Forty unique, consistently ordered tips, panel IDs, and matrix names |
| Saved CSV/JSON | Consistent with the actual plotting inputs and support audit |
| PDF | One page, 170 × 220 mm, embedded Arial TrueType faces, no clipping |
| PDF scale bar | Error below the 0.000001-point tolerance |
| PNG | 600-dpi export; measured scale bar agrees with PDF within one pixel |
| Full PNG correspondence | Zero differing pixels against a fresh in-memory export |
| Full PDF correspondence | Zero differing pixels at 200 dpi against a fresh in-memory export |
| Caption and alt text | Deposit rows, support definitions, scale, and sparse labels agree; alt text 236 characters |
| Existing/protected outputs | Refused without overwriting files or creating the protected destination |

See [figure validation output](validation-output.txt) and
[regression output](test-output.txt). The package's code, inputs, and documents
remained byte-identical during the clean-directory execution.

## Independent refactor check

The compact metadata extract matches all 40 retained workbook rows and all
40 source tip mappings. The matrix copy is byte-identical. Comparing the
accepted plotting source with the portable version in the same environment
gave identical panel structures, tree positions, ANI arrays, legend strings,
support audits, mappings, and rendered pixels. The original accepted PDF/PNG
remain under `reference/`.

## Validation improvement

The archived validator checked in-memory plotted data and exported geometry
but could pass an exported image with an obscured heatmap. This package adds
complete saved-image comparison against fresh exports from the checked artists,
plus verification of saved support-count summaries and mapping metadata.

Negative controls now reject a transposed matrix, one changed ANI cell, a
wrongly filled support marker, one changed PNG pixel, a blanked PNG region,
and changed PDF content. Separate full-figure probes rejected a blanked heatmap
in either exported format and a support count altered to 999.

## Environment and limits

Python 3.12.14, macOS arm64, Matplotlib 3.10.9, NumPy 2.2.6, pdfplumber 0.11.10,
Pillow 12.3.0, pypdfium2 5.13.0, and Poppler 26.05.0. Complete Python dependency
versions appear in `requirements-lock.txt`.

These checks verify this package and these archived inputs; they are not a
claim of universal correctness. PDF image comparison is at 200 dpi, not an
identity check of every PDF object. Exact output comparison requires the same
font/software environment for generation and validation. Cross-platform
execution has not been tested. The complete upstream KBase export, new tree
inference, ANI recomputation, other manuscript analyses, public deposition,
and licensing approval are outside this completed check.
