# Figure 2 verification — 20 September 2026

## Executed checks

The module was copied to a new temporary directory and run from outside that
directory with a fresh Python 3.12.14 virtual environment. All pinned packages
were installed from `requirements-lock.txt`; no existing scientific environment
or external project-data path was required. Poppler supplied the PDF font check.

| Check | Observed result |
|---|---|
| Figure 2 regression and failure tests | 35 passed |
| Existing Figure 1 tests | 29 passed; all 32 Figure 1 module files byte-identical |
| Raw family object | 21,853 families, 252,945 unique gene/genome members, 40 genomes |
| Saved copy-number matrix | All 874,120 cells matched the raw object |
| Clade records | All 36 clades, 144 counts, and 15 aggregate rows checked |
| Tree cross-check | Historical and corrected trees had identical rooted clade membership |
| Occupancy partitions | 4,263 core; 307 soft core; 2,576 shell; 14,707 cloud |
| Endpoint comparison | 312 additional families; core decreased by three |
| Plot values | 72 individual-clade points, 30 mean points, and 40 bars checked |
| Exported PDF bars | Positions, log-scaled heights, and category colors matched; maximum geometry error below 0.000001 pt |
| PDF | One 170 × 105 mm page, embedded TrueType fonts, no font warnings, text within page bounds |
| Fresh export versus accepted reference | Complete PNG pixel arrays identical; complete PDF page render arrays identical at scale 2 |
| Output safety | Existing and protected destinations rejected with exit 2; existing results and inputs unchanged |

The comparison used Arial on macOS. PDF byte equality is not required because
creation metadata may differ; rendered page equality was tested separately.
The accepted two-panel reference remains unchanged. The caption contains no
literature citation, and its alt text passes the 255-character check.

## Negative controls

Tests reject changed matrix cells, family categories/counts, clade membership or
counts, duplicate families/genes/genome mappings, unknown genomes, empty families,
malformed numeric values, nonfinite JSON values, duplicate JSON keys, altered
input hashes, changed plotted means, incorrect bar heights/colors, and altered
exported bar geometry. Safety checks also cover symlink destinations and Python
optimization mode. Tests are bounded checks, not a proof that no defect can exist.

## Limits

The calculation starts from saved ortholog-family assignments; it does not
repeat upstream annotation, OrthoMCL clustering, KBase execution, or phylogenetic
inference. The figure retains the original family definitions and sampling.
No fitted model of pangenome openness is provided. A verified local package is
not proof of a completed GitHub upload; a published copy must be downloaded and
checked separately. Complete upstream provenance and release permissions remain
separate tasks.
