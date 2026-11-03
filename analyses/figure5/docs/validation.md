# Validation record - 22 September 2026

## Verified

- Both FASTA sequences match the archived JSON, sequence MD5s, and protein lengths.
- All eight TSV annotations reconcile with the JSON at the correct E-value scope.
- Nine matched fragments, two full-protein bars, and the motif at residues 211-219 match the exported PDF geometry.
- Every accession, residue-range label, and E-value is checked on its own rendered row.
- All text is at least 7.5 pt at the native 170 x 160 mm export size; embedded fonts are TrueType, and no text falls outside the page.
- The PDF was rasterized and visually inspected separately from the PNG export.
- All 35 tests pass. These include deliberately filled gaps, shifted segments, incorrect or swapped E-values, swapped accessions, shifted motif positions, source corruption, and overwrite attempts.
- A fresh pinned environment and a relocated source copy produced byte-identical PDF, PNG, and plotted-data JSON outputs on macOS arm64 with Python 3.12.14. Cross-platform encoded PDF byte identity is not promised.
- The preserved historical plotting script also runs on an isolated copy of its original relative-directory layout. Its continuous SSF102405 bar is the known presentation limitation corrected here.

Machine-readable checks and test output accompany this record. Passing tests
does not establish biochemical activity, reproduce the historical web service,
or exclude every conceivable defect. External source archives were not changed.

## Manuscript relationship

The 22 September manuscript has the correct eight annotation values and scoped
conclusions. Existing Methods references [46-52] cover InterProScan, InterPro,
Pfam, PANTHER, Gene3D/CATH, and SUPERFAMILY; Results reference [68] identifies
Naseem, Bencurova, and Dandekar (2018) for the LOG/decarboxylase distinction.
No new citation is required for displaying already archived fragment boundaries.
Reference numbers are observations of that manuscript version, not permanent IDs.

The internal Methods and Data Record sections 127-129 document the scope of the
protein interpretation and prior Figure 5 revision. This module preserves the
absence of unverified active-site markers and the sequence-prediction caveat.
The new display and its caption supersede the earlier continuous-bar rendering,
not the source annotations. The manuscript itself has not been edited here.
