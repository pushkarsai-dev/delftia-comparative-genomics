# Template for one analysis

Copy to `analyses/<analysis-id>/README.md`. Replace every bracketed field.
Leave unresolved fields explicitly pending; do not invent provenance.

## Identity and scope

- Analysis ID: [stable lowercase identifier]
- Question or result supported: [specific purpose]
- Manuscript heading, figure, or table: [stable anchor]
- Methods and Data Record: [section IDs, exact headings, and record date]
- Superseded record entries: [IDs and what changed, or none]
- Computational status: [planned / staged—not tested / tested—bounded scope]
- Release status: [private / approval pending / public / archived release]
- Reproduction boundary: [primary-data calculation, saved-result rendering,
  fixture test, or other precisely bounded scope]

## Inputs and provenance

| Repository path or accession | Source/version/date | SHA-256 or archive identifier | Redistribution status |
|---|---|---|---|
| [input] | [origin] | [checksum/identifier] | [permission or restriction] |

State all transformations, filtering, mappings, and any data that must be
downloaded separately. Distinguish original inputs from derived intermediates.

## Software and effective parameters

[Versions, environment file, tools/databases, parameters, seeds if applicable,
and third-party references. Separate configured from execution-verified values.]

## Run instructions

[Exact commands, working-directory assumptions, new output directory,
approximate runtime/memory/disk needs, and any explicit download step.]

## Outputs

[Filenames, schema/units, expected counts, and relationship to the manuscript.]

## Verification actually performed

[Date, environment, commands, independent expected outputs, tolerances, test
results, negative controls, and figure inspection if relevant. State what
was not rerun. Do not use future-tense plans as evidence of completed tests.]

## Limitations and pending work

[Missing source/provenance, unresolved discrepancies, access/rights limitations,
and next action. Distinguish scientific, computational, and release limitations.]

## Change history and references

[Changes from the archived workflow, affected manuscript claims, scientific
references, tool attribution, and actual contributors.]
