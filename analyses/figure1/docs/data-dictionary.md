# Compact input dictionary

All files are UTF-8 text. Trees use Newick with numeric internal support labels.
An absent internal label means support was not reported, not zero or one.

## ANI matrix

`ani40_matrix.tsv` has a header row, a first column of query identifiers,
and 40 reference-genome columns. Cells are directional ANI percentages.
The saved row/column order is explicitly remapped to tree order for display.
The diagonal is recorded data, not a display replacement.

## Panel metadata

`panel_metadata.tsv` retains the source worksheet's seven column names:

- `Panel ID`: project genome identifier.
- `Assembly accession`: source assembly accession or study-assembly descriptor.
- `Strain`: archived strain label; shortened only for display where documented.
- `Species (NCBI)`: archived public taxonomic label.
- `GTDB call (fill)`: archived GTDB label, when provided.
- `Tier`: panel inclusion tier, T1, T2, or T3.
- `Type material`: archived type-material annotation; not a substitute for
  independent deposit-equivalence evidence.

## Tip mapping and output mappings

`tip_mapping.tsv` links KBase `tip_id` to its source `label`; all included
rows have `in_panel=yes`.

Generated `row_mapping.csv` links each figure-local `row_number` to its
`display_label`, `species`, `study` status, `panel_id`,
`assembly_accession`, `ani_matrix_name`, `kbase_tip_id`,
`deposit_symbol`, and `column_number_printed` flag.

Generated `support_audit.csv` records the number of tips in each retained
internal clade, its support class, its value (blank if unreported), and its
space-separated tip IDs. The root has its own `root_no_value` row, with a
blank value, and is also counted separately in `row_mapping.json`.
