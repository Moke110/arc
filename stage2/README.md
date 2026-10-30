# Stage 2 — model setup and selection

Stage 2 chooses the open-source whole-transcriptome drug-perturbation model that
will produce the perturbomes feeding future MIRI drug scoring. Setup ("does the
model run in its own environment?") is deliberately separate from selection
("which model is best?"): `01_setup/` answers the first, `03_models_cmbench/`
the second.

## Current state (2026-10-30)

| step | question | what it produced | record |
|---|---|---|---|
| [`01_setup/`](01_setup/README.md) | does each model run? | 9/9 run cards, one pinned environment per model; **no model emitted a whole-transcriptome CM perturbome** | ADR 0021, ADR 0024 |
| [`02_cm_data/`](02_cm_data/README.md) | what do we evaluate on? | frozen corpus: **18,164 genes × 525 drug keys**; five `datasets-1` h5ad views | ADR 0026 |
| [`03_models_cmbench/`](03_models_cmbench/README.md) | which model is best? | executed **60-unit cold-drug benchmark** in two substeps; primary 18k leaderboard + 5k HVG sensitivity fork | ADR 0027, ADR 0028 |
| [`04_tahoe_data/`](04_tahoe_data/README.md) | what do we pretrain on? | **Tahoe-100M pseudobulk corpus**: 67,018 profiles × 62,710 genes, 25.8 GB, QA gates pass | ADR 0029 |

The stage is now in **selection follow-up**: the benchmark evidence exists, the
pretrain corpus exists, and no model has yet been trained on it inside Stage 2.

## Selection constraints in force

- **Delivery requirement: 15,000+ genes.** The eventual MIRI perturbome must be
  whole-transcriptome, so a 978-landmark model cannot be the delivery model.
- **978 quarantine.** PRnet's native demo is the 978-gene LINCS matrix; it is
  allowed only inside `01_setup/execution/prnet/` and its outputs are labelled
  *setup-only*. No later Stage 2 evaluation, metric, fold or comparison may read
  a 978 artifact (ADR 0024). The whole-transcriptome guard still governs selection.
- **No entrant beats the constant baseline at 18,164 genes.** On the primary
  metric (k=200 direction-aware macro AUPRC) `delta_mean` scores **0.1989**; the
  best entrant, PerturbNet, scores **0.1861**. On the 5,000-gene HVG axis PRnet
  reaches **0.2735** against `delta_mean` **0.2491** (+0.0243 paired), so the
  18k conclusion is **gene-space dependent** — the 5k run is a sensitivity
  result, not a replacement leaderboard (ADR 0028). Absolute AUPRC is not
  comparable across the two axes; compare paired gaps only.
- **Squidiff is excluded from further selection, tuning and debugging**
  (user-approved, 2026-10-26) after both the 18k and 5k recipes failed to
  converge and further remediation offered uncertain benefit. This is an
  opportunity-cost decision, not an architectural gene limit or a demonstrated
  compute-budget failure.
- **Remaining candidates:** PRnet, CPA, chemCPA, biolord, PerturbNet, CMonge;
  keep the CPU baselines; keep every historical seven-model artifact and
  registry for reproduction. New queues exclude Squidiff unless the user reopens
  the decision. The governing record is
  [`audit/model_selection_decisions.md`](audit/model_selection_decisions.md).

## What each step produced

- **`01_setup/`** — one run card per deployed model, each with its own conda env
  or venv overlay, a run command and a saved artifact. "9/9" means the code path
  runs, **not** that the model is usable for selection: chemCPA's public weights
  are L1000 977/2,000-gene; biolord ships no weights; several others had no
  checkpoint at all. This is why selection needed its own benchmark and a
  gene-space/retraining decision rather than a leaderboard.
- **`02_cm_data/`** — the model-agnostic CM corpus and its filtered views:
  18,164 genes, 525 resolvable drug keys (`atc-moa` 177, `atc-moa-mesh` 316,
  `all-smiles` 525), five `datasets-1` views including the benchmark's
  `temp_env_all-smiles.h5ad` (4,478 wells, 460 drug keys, 16 plates, three doses,
  24 h, `ipsc_cm`). Reruns reproduce the frozen `manifest.csv`.
- **`03_models_cmbench/`** — the executed benchmark: five CPU baselines plus
  seven real adapters × five cold-drug folds = 60/60 units, one fixed
  group/fold/metric contract, 12-entry registry, per-unit fidelity ledgers, and
  the conclusion tables. `03.1_cold_drug_18k/` is the primary leaderboard;
  `03.2_hvg_sensitivity_5k/` re-runs the same wells, folds, entrants, seeds and
  budgets on a 5,000-gene HVG axis to test whether the ranking depends on gene
  space.
- **`04_tahoe_data/`** — the pretraining substrate: Tahoe-100M aggregated to
  `(plate, well, cell_line)` pseudobulk profiles, 62,710 native genes, summed
  raw counts plus `log1p_cp10k`, controls as rows, Δ deliberately not
  materialized. QA: median 25,290 detected genes/profile; 3,755 contrasts vs
  the authors' DESeq2 DE give Pearson 0.779 at `baseMean > 500` (gate ≥ 0.60)
  and 0.812 at `baseMean > 800`; all 700 `(plate, cell_line)` pairs carry a
  vehicle.

## Deployed vs reference-only models

Nine models are deployed in Stage 2 setup (PRnet, CPA, chemCPA, biolord,
PerturbNet, scVIDR, Squidiff, CMonge, DECANT); five are reference-only (WAVE —
no dose input, ExPO — 978-only, State, TranSiGen, XPert) and have no run card.
Preprint models (X-Pert, DORA, BaiZe, PerturbLDM) are excluded and not vendored.
scVIDR and DECANT exit at selection (no cold-drug structure path); the
per-item detail, origins and licences are in
[`../pert_model/README.md`](../pert_model/README.md).

## Environments

One environment per model, so a model's pins never leak into another's; the
shared `arc` environment is not modified (ADR 0021). Seven models use the
`arc-stage2-bench-*` conda environments, biolord and CMonge use the venvs under
[`_envs/`](_envs/README.md). `env.lock` / `*.lock` files are `pip freeze` /
`conda list --explicit` **records** of the environments the runs used — the
environments are already built, so the locks are documentation, not
re-resolution instructions. The venv trees themselves are machine-local and
gitignored (multi-GB, interpreter-specific); rebuild recipes are in
`_envs/README.md` and each run card.

## Where the decisions live

- ADR [0021](../docs/adr/0021-per-model-environments-stage2.md) — one conda environment per model.
- ADR [0026](../docs/adr/0026-stage2-cm-data-contract.md) — frozen CM corpus and filtered views.
- ADR [0027](../docs/adr/0027-stage2-benchmark-contract.md) — the 18,164-gene cold-drug benchmark contract.
- ADR [0028](../docs/adr/0028-stage2-5k-hvg-benchmark.md) — the 5,000-gene HVG sensitivity fork.
- ADR [0029](../docs/adr/0029-tahoe-pseudobulk-pretrain-corpus.md) — the Tahoe-100M pseudobulk pretrain corpus.
- [`audit/model_selection_decisions.md`](audit/model_selection_decisions.md) — the live selection override and what follows from it.
- [`_archive/PLAN.md`](_archive/PLAN.md) — the superseded plan, archived as-is; this page carries its live content.

Stage 1 (`stage1/`) produces the MIRI disease inputs and DEG lists Stage 2's
models will serve; it is not a dependency of this stage's setup or benchmark.
