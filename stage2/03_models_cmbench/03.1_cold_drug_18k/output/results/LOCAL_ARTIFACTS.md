# Local-only artifacts — 03.1_cold_drug_18k

The conclusion tables next to this file (`leaderboard.csv`, `per_drug.csv`,
`fidelity.csv`, `run_cards.csv`, and `bootstrap.npz` where present) are tracked:
they are the results of the 60-unit run. The items below are deliberately
**not** in git (large and/or regenerable); their path and digest are recorded
here so a local copy is never an anonymous blob.

| artifact | files | bytes | sha256 | why local |
|---|---:|---:|---|---|
| `metrics.csv` | 1 | 53,415,110 | `3a551a46b9fc04506f6c117982784ab69c9c12630a7bfb487985161b0d5c6fdb` | 53 MB long-format metrics (group/drug/fold/headline, all and atc-moa subsets); regenerate with `execution/script/evaluate.py` |
| `../input/` | 14 | 1,147,667,040 | per-file hashes in the tracked `../input/manifest.json` | frozen benchmark input arrays; the manifest is the durable digest, the arrays are local |
| `../units/` | 1108 | 34,478,207,145 | — | per-unit fit/predict artifacts and the shared `_cases/` bundles; regenerable by `execution/script/run_queue.sh` |

Digests measured 2026-10-30 with `sha256sum` and `stat -c%s`. `03.1_cold_drug_18k` is the
18,164-gene primary substep.
