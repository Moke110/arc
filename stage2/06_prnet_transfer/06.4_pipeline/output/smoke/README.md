# 06.4 smoke record — local artifacts + tracked numbers

`make gpu-smoke` runs `execution/smoke_real.py` on the frozen 06.3 input and
writes its artifacts here. Only the small record is tracked:

| file | tracked | content |
|---|---|---|
| `smoke_record.json` | yes (`git add -f`) | per-run wall time, peak RSS, per-epoch validation cost split, cold/warm cache epochs, both initialisations, headline AUPRC sanity, checkpoint sha256, git SHA, data hashes, concurrency probe |
| `smoke_summary.csv` | yes (`git add -f`) | one row per run (predictions/evaluations/checkpoints stay local) |
| `pretrain_*/`, `cm_fold*/`, `logs/`, `hash_cache.json`, `*.npy`, `*.pt` | no | run trees, weights, raw predictions, logs |

**A smoke run is pipeline evidence, not a science result**: 2 epochs, a small
Tahoe row subset, one CM fold, no selection claims.
