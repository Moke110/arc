# 06.4 smoke record — Phase-B real-data evidence

`make gpu-smoke` runs `execution/smoke_real.py` on the frozen 06.3 input.  The
**Phase-B run is complete** (2026-11-02 15:04–15:09 Asia/Shanghai, commit
`99650d8`, clean tree, frozen input freshly hash-verified: 8 artifacts / 15.30 GB
/ hash-cache hits 0).  All 17 steps passed: pretrain (subset) → CM finetune fold 0
with **both** initialisations → predict Δ → evaluate direction AUPRC@200 → model
card, plus a warm-cache epoch and a 2-concurrent pretrain probe.

| file | tracked | content |
|---|---|---|
| `smoke_record.json` | yes (`git add -f`) | the authoritative numbers: wall time, peak RSS, per-epoch validation cost split, cold/warm cache epochs, both initialisations, headline AUPRC sanity, checkpoint sha256, git SHA/dirty, data hashes, concurrency probe |
| `smoke_summary.csv` | yes (`git add -f`) | one row per run |
| `pretrain_seed*/model_card.md`, `cm_fold0_random/model_card.md`, `cm_fold0_tahoe/model_card.md` | yes (`git add -f`) | model cards with the evaluation attached (each states it is a pipeline smoke) |
| `pretrain_*/`, `cm_fold*/`, `concurrency_probe_*/`, `logs/`, `hash_cache.json`, `*.npy`, `*.pt` | no | run trees, weights, raw predictions, per-step logs |

Measured peaks for this run: **GPU 2,945 MiB** of 16,376; **host 8.54 GiB** (the
Tahoe-valid evaluation over 8,524 rows / 8,110 conditions — the run's peak;
training peaks at 2.70 GiB pretrain / 1.33 GiB CM).  The `MemoryMax=12G` scope cap
was never tripped, swap stayed unused, and total wall was 4:29 including the fresh
15.30 GB artifact hash.

Recorded checks: both training checkpoints `reload_identical: true`; the
Tahoe-initialised CM run reloads the pretrain checkpoint in a separate process
through the strict axis/structure-checked loader; all five prediction runs report
`repeated_inference_identical: true` (zero noise, `model.eval()`, fixed batch
order); every CM evaluation re-derives the input baseline and would fail unless the
prediction used this fold's fixed c\* (max abs diff ≤ 1e-4) — all passed.  The
commit's 120-test suite including the 20 leakage tests also passes.

**A smoke run is pipeline evidence, not a science result**: 2 epochs, a small Tahoe
row subset (3,000 train / 1,000 valid), one CM fold, no selection claims.  The
headline direction AUPRC@200 values in the record (0.0466 pretrain-valid; 0.0894
CM test for both initialisations) are sanity numbers only.
