# v2 step-05 restructure: full/ bench + drug-partitioned data-ablation

**Status: accepted**

v2's step `05_xpert_finetune_bench` gains the same two-benchmark structure v3
adopted under ADR 0018. The official fold benchmark — whose fold test sets all
come from the dtoxs_primary arm (byte-identical to v1.1, ADR 0014/0015), the
CONTEXT.md *dtoxs-primary benchmark* — is relocated **verbatim** into
`05_xpert_finetune_bench/full/` (scripts, cached `output/`, `benchmark.png`;
no re-run, the fold columns in `cm_union.h5ad` are untouched). The pretrain-eval
config stays at the step-top-level `config/` (`cm_pretrain_eval.yaml`, shared
with the ablation — the XPert loader resolves a config's relative dataset paths
against its third parent, so configs must sit two levels under the workflow
root). Step 05 additionally gains a **data-ablation experiment** (outputs under
`05_xpert_finetune_bench/data-ablation/`, ADR 0019) that measures how the
finetune-corpus **scale** affects XPert performance on v2's own union.

The ablation mirrors v3's design (ADR 0018) **adapted to v2's 45-drug union**
by partitioning test sets **by drug, not by row**: a row-uniform ~15% test
(v3's recipe on its ~180-drug union) would draw ~40 of v2's 45 drugs into
test and collapse the cold-drug train pool to a handful of rows. Instead the
45 union drugs are split into a **disjoint tripartition of 3×15 test drugs**
(fixed master RNG): seed s tests drug-third s, so every union drug is tested
exactly once across the 3 seeds and the seeds' test sets never overlap. Each
seed's **unified test set** is ~15% of the union (~144 rows) drawn only from
its 15 test drugs (per-drug allocation proportional to row count minus one;
every test drug keeps ≥1 row in the pool, so warm-start f=1 overlap is 1.0).
Train = a dataset-stratified f-subsample of the seed's pool over 5 train-size
fractions **{20%, 40%, 60%, 80%, 100%}** (warm-start: union ∖ test, ~814 rows,
all 45 drugs; cold-drug: the 30 hold-out drugs only, ~610–655 rows, strictly
drug-disjoint at every fraction). Training recipe identical to the official
fold bench (4 epochs, lr 0.002, batch 32, `cm.yaml`, step-04 L1000 pretrain
init); the L1000-pretrain weight anchors the curve at 0 rows (3 zero-shot
evaluations). Outputs: mask manifest `ablation_masks.csv`, `ablation.csv`
(30 finetune + 3 pretrain rows), the 4×4 `ablation.png` (warm/cold ×
deg `RMSE|PCC|Spearman|R2` + VCC `PDS|JAC|FID|REACH`), and
`ablation_report.md`.

Decisions within the experiment:

- **Test sets are drug-partitioned (not row-uniform) because v2's union is
  small.** Rationale as above; the tripartition also makes the 3-seed mean ±
  std average over three complementary drug-thirds rather than overlapping
  resamples. This is a deliberate deviation from v3's row-uniform recipe
  (which works there because ~45 non-test drugs remain after a 15% draw).
- **Fractions are arithmetic {20…100}%, not v3's {1/16…1.0} log curve**, since
  the pools are small; the cold curve's f=100% point is the hold-out-drug
  scale (~30 drugs), not union scale — the warm/cold x-ranges differ by design
  and are documented in the report.
- **1.0 is re-run** on the ablation's own split scheme (never the full/
  bench, whose test sets differ and are primary-only).
- The ablation drives `src.train` through temporary h5ads (fold columns
  `ablation_{20..100}`) and **does not pollute `cm_union.h5ad`** (its 21
  official fold columns stay intact). Per-run YAMLs are generated under
  `05_xpert_finetune_bench/config/` (the loader's third-parent nesting requires
  configs two levels under the workflow root — `full/config/` would misresolve);
  the pretrain-eval config `config/cm_pretrain_eval.yaml` is shared with the
  full/ bench (no copy).
- The ablation is a **diagnostic analysis artifact**: it does not change the
  final screen model (step 06) or any downstream step.

**Consequences:** `full/` is v2's official *dtoxs-primary benchmark* in the
CONTEXT.md sense but keeps the name `full/` (a deliberate asymmetry with v3's
`dtoxs-primary/` directory — v2 has no added train-only datasets, so the
directory name says "the whole prior bench, moved intact"). Entry points and
step-05 gotchas in `stage1/xpert-fda-miri-v2/AGENTS.md`, the screen notebook
(`_lib/build_screen_notebook.py`), root `AGENTS.md`, and CONTEXT.md document
both benchmark methods. The v3 ablation keeps its row-uniform recipe
unchanged.
