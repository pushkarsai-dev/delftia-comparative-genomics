# Adding an analysis

Add one self-contained, verified analysis at a time. `main` contains the
accepted working versions for the whole manuscript, not a separate figure.

## 1. Give each analysis its own folder

Choose a stable lowercase identifier such as `figure2`, `panel-ani`,
`pangenome`, or `enzyme-family-searches`:

```text
analyses/<analysis-id>/
  README.md                 purpose, manuscript mapping, and exact commands
  scripts/                  executable code
  data/                     small inputs or accession/download manifests
  reference/                accepted expected outputs, never overwritten
  tests/                    regression tests and negative controls
  docs/                     provenance, parameter records, and validation
  requirements.txt          pinned Python dependencies, if applicable
  environment.yml           alternative environment definition, if needed
  build/                    regenerated outputs, excluded from Git
```

Create only the files and directories the analysis needs. Copy the
[README template](docs/analysis-template.md) into the new module.
Do not add empty placeholders for work that has not been deposited.
Keep dependencies specific to each module; do not silently change Figure 1's
environment to make another analysis run.

## 2. State provenance and computational scope

Record relevant Methods and Data Record section identifiers, record date,
exact manuscript heading/figure/table, and superseded entries. Do not guess
section numbers or use manuscript page/line numbers as permanent anchors.

Distinguish plotting saved results from recalculating them from primary data.
Record each input's source, version/date, checksum, and permission to redistribute.
When analyses share inputs, identify the authoritative source and exact checksum;
do not let duplicate copies diverge silently.

Do not upload internal correspondence, the full internal Methods Record,
private preparation logs, credentials, or unrelated files. Large sequence
datasets normally belong in sequence archives; provide accessions, retrieval
instructions, and checksums instead of bundling raw reads.

## 3. Make execution portable and safe

- Resolve inputs relative to the module or an explicit input argument.
- Accept a fresh output directory; never overwrite inputs or accepted outputs.
- State software versions, effective parameters, and expected resources.
- Separate downloads from offline computation where possible.
- Test the exact documented commands from a clean location.
- Include independent output comparisons and useful failure tests.
- Visually inspect regenerated figures as well as checking their values.

Run tests explicitly for the module, for example:

```sh
python -m unittest discover -s analyses/figure1/tests -v
```

Use that module's environment. Do not assume a test-discovery command at the
repository root finds nested tests automatically. Rerun affected modules when
shared inputs or dependencies change. A fixture test is not a full-panel rerun.

## 4. Review and record the addition

1. Update [the analysis index](docs/analysis-index.md), with an honest scope and status.
2. Describe the change in [CHANGELOG.md](CHANGELOG.md).
3. Review the diff for unintended source changes and private information.
4. Check existing hashes before refreshing `SHA256SUMS.txt`; explain intentional changes.
5. Commit one coherent module. Prefer a feature branch and review before merging scientific changes.

The root manifest covers all tracked files except itself. Exclude generated
outputs, virtual environments, caches, and Git internals. Input manifests are
scientific provenance: do not regenerate hashes merely to silence a mismatch.

## Credit and release

A collaborator invitation is not manuscript authorship or a contribution in
GitHub's commit-based contributor list. Credit actual contributions accurately.
Original project code is covered by MIT, and original study data by CC BY 4.0,
subject to the component scope and third-party exclusions in docs/licensing-scope.md.
Preserve third-party provenance and notices when contributing. For a new release, update the version, archive DOI, citation metadata, and checksums together.
