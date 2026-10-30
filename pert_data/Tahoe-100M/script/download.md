# Download — Tahoe-100M

**method: dataset** — HuggingFace dataset `tahoebio/Tahoe-100M`, pinned
revision `349e2a55f5538a3e73b32a06c6bdd8e9e87011e8`, license CC0-1.0.
**repo:** <https://huggingface.co/datasets/tahoebio/Tahoe-100M>

## Canonical one-liner

The whole atlas as an unmodified HuggingFace snapshot into `source/`:

```bash
cd pert_data/Tahoe-100M
hf download tahoebio/Tahoe-100M \
  --revision 349e2a55f5538a3e73b32a06c6bdd8e9e87011e8 \
  --local-dir source/
```

## The import actually used (resumable)

`script/fetch_tahoe.py` is the importer that produced the tree here; it writes
the same snapshot plus `MANIFEST.jsonl` (HF `lfs.sha256` per file), `AXES.json`
and `DATASET.json`, and keeps its ledger/status/HF cache under `local/_state/`:

```bash
# from the repository root
python pert_data/Tahoe-100M/script/fetch_tahoe.py --scope full   # 3,388 shards + 1,026 DE shards
python pert_data/Tahoe-100M/script/fetch_tahoe.py --verify-only deep
```

`script/tahoe-dl.sh {start|pause|status|logs}` runs that fetch as a transient
systemd `--user` unit; `script/run_fetch.sh` is its wrapper. The importer accepts
`hf-mirror.com` via `HF_ENDPOINT`.

## What lands in `source/`

4,430 files / 428,811,080,261 bytes: `data/` (3,388 expression parquet shards),
`metadata/` (1,033 files incl. the authors' `pseudobulk_differential_expression`),
`AXES.json`, `DATASET.json`, `MANIFEST.jsonl`, `README.md`, `LICENSE.md`,
`tutorials/`, `dataset_size_comparison.png`, `.gitattributes`.

`local/` holds `Tahoe_summary.md` and `_state/` (ledger, status, probes, HF/xet
cache); the HuggingFace transfer cache `tahoe/.cache/` of the pre-split tree was
deleted. `script/verify_import.py` is the offline import check; `script/verify.sh`
is the manifest verifier used by `script/verify_all.sh`.
