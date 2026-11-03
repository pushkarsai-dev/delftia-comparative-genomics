# Reproduce and verify a figure

Use this guide to choose a figure workflow and check downloaded files. The module README is the authority for that figure's commands, inputs, outputs, and scientific scope.

## 1. Obtain the source

Download the source ZIP from [release v1.0.0](https://github.com/pushkarsai-dev/delftia-comparative-genomics/releases/tag/v1.0.0) and extract it, or clone that tag:

```sh
git clone --branch v1.0.0 --depth 1 https://github.com/pushkarsai-dev/delftia-comparative-genomics.git
cd delftia-comparative-genomics
```

Record the commit or release tag you used. In a Git clone, `git rev-parse HEAD` prints the commit. Use the named [evidence archives](evidence-downloads.md) for the separate supporting data. Historical source identifiers inside packages describe the original review history. The public checkout target is `v1.0.0`.

## 2. Check file integrity

From the extracted repository root on macOS, run:

```sh
shasum -a 256 -c SHA256SUMS.txt
```

Every listed file should report `OK`. A missing file or checksum mismatch must be resolved before using the package. This manifest covers tracked files except itself; it checks file integrity, not scientific correctness.

## 3. Select one module

The documented figure workflows were tested with Python 3.12.14 on macOS. Figures 1–4 use Poppler's `pdffonts` for export validation; install Poppler through your package manager and make it available on `PATH`. Figure 5's validator reads its archived data and exported PDF directly. Each module pins its own dependencies. Installation needs network access; the saved-input figure workflows run offline afterward.

| Figure | Module instructions | Main generated output | Scope |
|---|---|---|---|
| 1 | [Phylogeny and ANI](../analyses/figure1/README.md) | `build/figure1/figure1.pdf` | Saved tree and ANI matrix; distance-preserving pruning; no original tree inference or ANI rerun |
| 2 | [Pangenome structure](../analyses/figure2/README.md) | `build/figure2/figure2.pdf` | Recompute summaries from archived ortholog-family assignments; no annotation or clustering rerun |
| 3 | [Protein-family counts](../analyses/figure3/README.md) | `build/figure3/figure3.pdf` | Plot archived counts and supplied search checks; optional upstream search has additional inputs |
| 4 | [Locus conservation](../analyses/figure4/README.md) | `build/figure4/figure4.pdf` | Plot and validate recorded evidence; separate commands repeat region-restricted searches |
| 5 | [Protein-domain assignments](../analyses/figure5/README.md) | `build/figure5/figure5.pdf` | Reconstruct the figure from archived InterProScan results; no new service run |

Paths in the output column are relative to the selected module directory. Each run also produces a PNG and module-specific records.

## 4. Follow the README in a fresh environment

Enter the selected module before following its setup instructions, for example:

```sh
cd analyses/figure2
```

Create and activate the module's virtual environment, install its specified requirements, and run its documented commands in order. Do not combine dependency sets across modules. Use a new output directory for every run; accepted files in `reference/` and source inputs must remain unchanged. When switching modules, deactivate the first environment before setting up the next.

Figure 2 is a compact first check. Its quick start, 35 tests, and export validator passed again on 7 October 2026 in a fresh Python 3.12.14 environment. All 874,120 matrix cells, 36 clades, and 40 exported bars passed; fresh PNG and PDF-render comparisons passed. This was a check of the documented saved-input workflow, not a rerun of upstream clustering.

## 5. Interpret success

A successful workflow exits without errors, produces the named outputs, passes the module's tests and export validator, and yields a visually correct figure. Read the validation output rather than relying only on the presence of a PDF. Fonts and library versions can alter pixels across platforms; the modules state which numerical and rendering comparisons are expected.

Passing computational checks does not establish biochemical function or recover a missing historical service record. Review the module's interpretation and limitations with the [analysis index](analysis-index.md).

## If a step fails

- Check that you are in the intended module directory and have activated its environment.
- Check the documented Python version, pinned dependencies, and any Poppler requirement.
- Choose a new output directory if the existing directory is refused.
- If an input hash differs, retrieve the correct source version; do not edit the manifest to hide the mismatch.
- Record the commit, module, command, software versions, and error when reporting a problem. Do not include credentials or private correspondence.

Return to the [repository home](../README.md) or [evidence downloads](evidence-downloads.md).
