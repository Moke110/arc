# AGENTS.md — stage2/02_cm_data

Routing index; the human record is `README.md`, the contract ADR 0026.

- `execution/build_cm_datasets.py` — the only entry point: `--check` (source
  counts), `--corpus-only` (verify/freeze `output/corpus/v1`), `--bundle
  datasets-1` (rebuild views), `--rebuild-corpus` (deliberate v2 promotion).
- `execution/cmdata/` — library: `sources`, `genes`, `drugs`, `views`, `corpus`,
  `qa`, `paths` (repo-root walk-up + frozen expectations).
- `execution/tests/test_cmdata.py` — 23 stdlib-unittest tests; `BuiltBundleTests`
  reads `output/` and skips if a bundle is absent.
- `output/corpus/v1/` — frozen gene + drug universes, no matrices (tracked).
- `output/datasets-1/` — five h5ad views + `manifest.csv` (authoritative),
  `compounds.csv`, `qa_report.md`, `gene_crosswalk.csv`. `*.h5ad` and
  `gene_crosswalk.csv` are gitignored; the CSVs/markdown stay tracked.

Run from the repository root; paths resolve by walking up to `.git`, never by
counting `../..`:

```bash
PY=/home/chang/miniconda3/envs/arc/bin/python
$PY stage2/02_cm_data/execution/build_cm_datasets.py --corpus-only
$PY -m unittest discover -s stage2/02_cm_data/execution/tests -v
```

Fetch the identity assets first if missing:
`pert_data/Annotation/script/fetch_annotation.py`. The build makes no network
call. Never hand-edit `output/corpus/v1/` or `qa_report.md` — the build
regenerates them; promoting a corpus version is a deliberate `--rebuild-corpus`
edit.
