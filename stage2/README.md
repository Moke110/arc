# Stage 2 — model setup and selection

Stage 2 chooses the open-source whole-transcriptome drug-perturbation model that
will produce the perturbomes feeding future MIRI drug scoring. Setup ("does the
model run in its own environment?") is deliberately separate from selection
("which model is best?"): `01_setup/` answers the first, `03_models_cmbench/`
the second.

## Current state (2026-11-04)

| step | question | what it produced | record |
|---|---|---|---|
| [`01_setup/`](01_setup/README.md) | does each model run? | 9/9 run cards, one pinned environment per model; **no model emitted a whole-transcriptome CM perturbome** | ADR 0021, ADR 0024 |
| [`02_dataset_corpus/`](02_dataset_corpus/README.md) | what data substrate do we build on? | **five source-dataset corpora** — four CM native-axis + `tahoe.h5ad`, 67,018 × 62,710, 25.8 GB — raw counts, every raw well (controls + blanks), no pharmacology filter and no Δ; **seven crosswalk tables** + coverage README (18,164-gene CM axis as `tables/genes.csv`); `--verify` **G1–G6 PASS**, with G4 reproducing 5/5 frozen v1 view content hashes | ADR 0033 (amends ADR 0026) |
| [`03_models_cmbench/`](03_models_cmbench/README.md) | which model is best? | executed **60-unit cold-drug benchmark** in two substeps; primary 18k leaderboard + 5k HVG sensitivity fork. **Frozen pending a redo on the merged corpus**: its inputs now live read-only in `stage2/_frozen/cm-data-v1/`, and no output is rebuilt or rewritten | ADR 0027, ADR 0028 |
| [`05_prnet_pretrain/`](05_prnet_pretrain/README.md) | can we train on it at all? | **executed 2026-11-01**: PRnet pretrained on **its own** CM-schema Tahoe view (`output/pretrain_view/tahoe_cm_schema.h5ad`, 67,018 × 18,164, 16.32 GB) — **68 epochs**, early stop at 68, best valid Gaussian NLL **−2.628080** at epoch 47, **956 s (15.9 min)** on an RTX 4070 Ti SUPER; pipeline claim only. Since the merge this step also owns the view (moved byte-identically) and its builder (moved, then repointed at the merged corpus and given a `--verify` replay) | ADR 0031, ADR 0033 |
| [`06_prnet_transfer/`](06_prnet_transfer/README.md) | does Tahoe pretraining improve unseen-drug CM prediction? | **executed 2026-11-02 and concluded — negative first round**: eligibility + QC audits, frozen five-fold input with consensus controls, 3 formal Tahoe pretrains, 90 CM candidates, 30-unit paired evaluation. Tahoe-init −0.0063 AUPRC@200 vs random (95% CI [−0.0074, −0.0052], n=118 drugs) and both inits far below `delta_mean` (−0.0800 / −0.0738); no transfer gain demonstrated under the pre-registered rule. The capacity (#10) and env-compound (#11) ablations are deferred. Paths only were reconnected after the merge; no result, hash or digest key changed | ADR 0032 |

**The two data steps are now one, and `04` is a deliberate gap.** On 2026-11-04
the retired Tahoe data step was folded into `02_dataset_corpus/`: one information-maximal
corpus layer, with every model-facing view owned by its consumer (ADR 0033). The
number `04` is not reused — a missing number is cheaper misdirection than renaming a
concluded study's nine substeps — so `05`/`06` keep their numbers. The five retired v1
CM views and the old identity tables live read-only in `stage2/_frozen/cm-data-v1/`
as a **temporary bridge**, deleted once `03`/`05`/`06` are redone.

The stage is now in **first-round conclusion**: the benchmark evidence exists, the dataset
corpora exist, PRnet has been pretrained on the CM-schema view, and the formal-transfer
study has finetuned that checkpoint and evaluated it against the benchmark. It **answered**
the pretrain→finetune question for this round, and the answer is **negative** — see
[`06_prnet_transfer/README.md`](06_prnet_transfer/README.md) and the ticket #9 results block
in [`06.9_evaluation/README.md`](06_prnet_transfer/06.9_evaluation/README.md). The
post-first-round ablations (#10/#11) are deferred, not planned work in flight, and the
benchmark redo on the merged corpus is a separate project.

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
- **`02_dataset_corpus/`** — the information-maximal data substrate: five raw-count
  h5ad corpora on their sources' **native** gene axes (`dtoxs_primary` 779 × 26,328,
  `dtoxs_ipsc` 312 × 26,328, `sctox` 3,072 × 31,880, `tenv` 4,510 × 19,613, `tahoe`
  67,018 × 62,710), one row per raw well including vehicle controls and technical
  blanks, no pharmacology filter, no integer encoding and no Δ. `output/corpus/v2/`
  also holds the seven crosswalk tables and the coverage README; the frozen
  18,164-gene CM axis survives as the projection target (`tables/genes.csv`), not as
  a corpus universe. The `--verify` suite passes G1–G6, including the reconstruction
  gate: the five retired v1 CM views are rebuilt from the v2 corpora and their
  content hashes match the archived manifest 5/5.
- **`03_models_cmbench/`** — the executed benchmark: five CPU baselines plus
  seven real adapters × five cold-drug folds = 60/60 units, one fixed
  group/fold/metric contract, 12-entry registry, per-unit fidelity ledgers, and
  the conclusion tables. `03.1_cold_drug_18k/` is the primary leaderboard;
  `03.2_hvg_sensitivity_5k/` re-runs the same wells, folds, entrants, seeds and
  budgets on a 5,000-gene HVG axis to test whether the ranking depends on gene
  space. The step is **frozen**: both substeps read `stage2/_frozen/cm-data-v1/`
  and nothing rebuilds or rewrites their `output/`.
- **`05_prnet_pretrain/`** — the first model trained on the Tahoe corpus inside
  Stage 2: PRnet's `PGM` from random init on the CM-schema view that this step now
  owns (`output/pretrain_view/tahoe_cm_schema.h5ad`, 67,018 × 18,164, 16.32 GB,
  moved byte-identically from the retired Tahoe step; 7,185,000 parameters;
  68 epochs, early stop at epoch 68, best valid Gaussian NLL −2.628080 at epoch
  47; 956 s (15.9 min) and 15,150 MiB peak RSS on an RTX 4070 Ti SUPER). A pipeline
  bring-up — no finetune, no cold-drug evaluation, no comparison with the
  benchmark.

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
- ADR [0026](../docs/adr/0026-stage2-cm-data-contract.md) — frozen CM corpus and filtered views (amended by ADR 0033).
- ADR [0027](../docs/adr/0027-stage2-benchmark-contract.md) — the 18,164-gene cold-drug benchmark contract.
- ADR [0028](../docs/adr/0028-stage2-5k-hvg-benchmark.md) — the 5,000-gene HVG sensitivity fork.
- ADR [0029](../docs/adr/0029-tahoe-pseudobulk-pretrain-corpus.md) — the Tahoe-100M pseudobulk corpus.
- ADR [0031](../docs/adr/0031-tahoe-cm-schema-view-and-prnet-pretrain.md) — the CM-schema Tahoe pretrain view and the PRnet pretrain bring-up.
- ADR [0033](../docs/adr/0033-stage2-dataset-corpus.md) — one dataset-corpus layer, views downstream; the Tahoe data step is folded in and `04` is a gap.
- [`audit/model_selection_decisions.md`](audit/model_selection_decisions.md) — the live selection override and what follows from it.
- [`_archive/PLAN.md`](_archive/PLAN.md) — the superseded plan, archived as-is; this page carries its live content.

Stage 1 (`stage1/`) produces the MIRI disease inputs and DEG lists Stage 2's
models will serve; it is not a dependency of this stage's setup or benchmark.
