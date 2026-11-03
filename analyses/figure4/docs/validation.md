# Validation scope

The automated validator checks 933 drawn gene polygons (45 reference genes
and 888 genes in the 20 comparative tracks) against source GenBank coordinates,
45 conservation bar heights, 20 ClusterBlast ticks, row labels/counts/spans,
the exported PDF polygons, the physical 10 kb scale, page dimensions,
embedded TrueType fonts, and fresh PDF/PNG rendering in the tested environment.
Essential upright text is at least 8 pt. The final PDF and PNG were visually
inspected for overlap and clipping.

Tests include intentionally shifted gene arrows, a changed conservation bar,
an altered scale-axis width, invalid coordinates/strands, a checksum failure,
the two identifier-fallback cases, and preservation of the historical table.
Successful tests mean those checks passed; they are not a universal guarantee.

The historical-parser search was rerun for all 40 selected regions and
reproduced 1,293/1,293 rows. The feature-aware search was also rerun for all 40,
yielding 1,344 rows with no lost historical match. Those comparisons and the
underlying best-hit evidence are in `data/`.

Panel (c) retains the historical summary files. It is independently
corroborated by the packaged reconstruction: 1,800 calls and matrix cells
validate against raw hits, and all per-gene and per-genome counts agree with
the plotting sources. The 54 reconstruction checks and four integration
checks supplement the 32 figure tests. Six archived repeat records also agree
with the primary reconstruction. A subsequent [full-panel repeat](../reconstruction/full_panel_repeat/README.md)
processed all 40 genomes with fresh predictions and no reused checkpoints.
All 1,800 selected pairs and decisions match, and the matrix is byte-identical.
Numerical and lower-ranked alignment differences are explicitly recorded.
No historical cell-level table survives, so historical pair identity remains
unverifiable even though the deposited reconstruction has now been repeated in full.

Five additional full-repeat comparator tests verify the saved comparison and
reject altered target IDs, duplicate calls, changed matrix cells, and altered
protein hashes. The current suite therefore contains 95 tests.

Journal guidance: [Microbiology Society figure preparation](https://www.microbiologyresearch.org/prepare-an-article)
calls for independently understandable figures, defined symbols, and suitable
editable high-resolution files; it does not require this particular page size
or a commercial figure application. The PDF is vector output; the PNG is for
manuscript placement. Alt text is below 255 characters.
