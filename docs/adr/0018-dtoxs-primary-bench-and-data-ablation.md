# Step-05 benchmark restructure: dtoxs-primary bench + data-ablation experiment

**Status: accepted**

The official step-05 fold benchmark's test sets all come from the **dtoxs_primary
arm** (byte-identical to v1.1/v2, ADR 0014/0015/0017; gse244740/gse262419 are
train-only), so it measures primary-arm test performance. Its outputs move to
`05_xpert_finetune_bench/dtoxs-primary/` and `benchmark.png` is renamed
`dtoxs-primary-bench.png` to say what it is. To measure how the **finetune-corpus
scale** affects XPert performance on the model's actual training distribution,
step 05 additionally gains a **data-ablation experiment** (outputs under
`05_xpert_finetune_bench/data-ablation/`) that is **decoupled from the official
folds**: 3 seeds × (warm-start | cold-drug) × 5 train-size fractions
{1/16, 1/8, 1/4, 1/2, 1.0}, each seed holding one **unified test set** (~15% of
the 4-dataset union, uniform random, composition mirroring the union) shared by
every fraction of both modes. Train = a dataset-stratified f-subsample of the
seed's train pool (warm-start: union ∖ test; cold-drug: union ∖ test ∖ all
test-drug samples across all four datasets). The training recipe is identical to
the official fold bench (4 epochs, lr 0.002, batch 32, `cm.yaml`, step-04 L1000
pretrain init); the L1000-pretrain weight anchors the curve at 0 rows. Single
deterministic seed per (seed, mode, fraction); per-run train/test well sets are
persisted as a mask manifest so repeats can be added later.

Why not reuse the official folds? Their test sets are primary-only and never
include gse244740/gse262419 rows (ADR 0017 contract), so a train-row subsample
over them could only ever measure primary-arm prediction — it cannot answer a
scale question about the union the model actually trains on. The ablation's own
union-wide test sets are a deliberate deviation from the ADR 0017 test-set
identity; the two benchmarks are independent and answer different questions.

Decisions within the experiment:

- **Warm-start purity degrades with scale by design.** f=1 guarantees every test
  drug is present in the train pool; at f<1 the dataset-stratified subsample may
  drop whole test drugs, and the number of absent test drugs per (seed, mode,
  fraction) is recorded as a diagnostic. cold-drug stays strict at every
  fraction.
- **Train subsamples are dataset-stratified** so the union composition is held
  constant and the horizontal axis is purely sample scale.
- **1.0 is re-run** on the ablation's own split scheme (it does not reuse the
  existing dtoxs-primary bench, whose test sets differ).
- The ablation drives `src.train` through temporary h5ads and **does not pollute
  `cm_union.h5ad`** (the 21 official fold columns stay intact).
- The ablation is a **diagnostic analysis artifact**: it does not change the
  final screen model (step 06) or any downstream step.

Also: `pert_data/cm_drugs.html` coverage table gains a per-dataset **in-HG 样本数**
column (662 / 296 / 2424 / 936 / 6, union 4,318) next to the existing
treated-sample counts.

**Consequences:** the cross-workflow bench page link
(`workflow/xpert-bench-v1.html`) is updated to the renamed png; the official
bench stays authoritative for cross-workflow comparison; step-05 gotchas and
entry points in `workflow/xpert-fda-miri-v3/AGENTS.md` document both benchmark
methods.
