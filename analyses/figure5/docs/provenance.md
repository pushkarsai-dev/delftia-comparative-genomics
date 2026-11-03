# Provenance

The inputs are the archived InterProScan results dated 14 August 2026 and the
submitted two-protein FASTA. All three are copied byte-for-byte, without rewriting
JSON, TSV, headers, or sequence order. SHA-256 hashes are in `data/input_checksums.json`.

| Component | Archived version |
|---|---|
| InterProScan | 5.78-109.0 |
| Gene3D | 4.3.0 |
| SUPERFAMILY | 1.75 |
| Pfam | 38.2 |
| PANTHER | 19.0 |

The original figure script's SHA-256 was
`32f345eaf54adcd341d60fbb6ca160f055fab5bf8ffc99491f88254706c9aa87`.
Its preserved copy was run separately to verify the historical plotting route.
The replacement is a portable plotting implementation, not the historical script.
Original source files and shared plotting styles were not modified.

The FASTA sequences match the JSON sequences and MD5s. The motif PGGFGTLDE
occurs once at residues 211-219 in HFJBPOEG_02821. Eight TSV records reconcile
with eight JSON matches. The JSON, not the TSV's outer bounds alone, supplies
fragment structure. It lists the SSF102405 fragments out of order; sorting by
start coordinate changes their presentation, not their values.

The absence statement in the caption is limited to the archived signatures
detected for this protein. It is not proof that a protein has no unrecognized
activity. Family predictions are not biochemical validation.

Sources for format and submission conventions:

- [InterProScan v5 output formats](https://interproscan-docs.readthedocs.io/en/v5/OutputFormats.html).
- [Microbiology Society article preparation](https://www.microbiologyresearch.org/prepare-an-article).

The figure's chosen dimensions and type sizes are readability decisions, not
claimed journal-mandated dimensions. The separate PDF is the vector submission
file; the PNG is provided for manuscript insertion and convenient preview.
