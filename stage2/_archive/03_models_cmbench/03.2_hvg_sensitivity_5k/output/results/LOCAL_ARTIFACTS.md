# Local-only artifacts — 03.2_hvg_sensitivity_5k

The conclusion tables next to this file (`leaderboard.csv`, `per_drug.csv`,
`fidelity.csv`, `run_cards.csv`, and `bootstrap.npz` where present) are tracked:
they are the results of the 60-unit run. The items below are deliberately
**not** in git (large and/or regenerable); their path and digest are recorded
here so a local copy is never an anonymous blob.

| artifact | files | bytes | sha256 | why local |
|---|---:|---:|---|---|
| `metrics.csv` | 1 | 53,163,051 | `95155eab63e9690b557892b8741ccb7fb7f616b6f973f4e3314b074a80d5571b` | 53 MB long-format metrics (group/drug/fold/headline, all and atc-moa subsets); regenerate with `execution/script/evaluate.py` |
| `../input/` | 17 | 320,392,653 | per-file hashes in the tracked `../input/manifest.json` | frozen benchmark input arrays; the manifest is the durable digest, the arrays are local |
| `../units/` | 1106 | 14,456,630,976 | — | per-unit fit/predict artifacts and the shared `_cases/` bundles; regenerable by `execution/script/run_queue.sh` |

Digests measured 2026-10-30 with `sha256sum` and `stat -c%s`. `03.2_hvg_sensitivity_5k` is the
5,000-gene HVG sensitivity substep.
