# Stage 2 — model setup and selection

Stage 2 chooses the open-source whole-transcriptome drug-perturbation model that
will produce the perturbomes feeding future MIRI drug scoring. Setup ("does the
model run in its own environment?") is deliberately separate from selection
("which model is best?"): `01_setup/` answers the first; the active `03_modes_tenv/` comparison answers
the second under ADR 0042.

## Current state (2026-11-08)

| step | question | what it produced | record |
|---|---|---|---|
| [`01_setup/`](01_setup/README.md) | does each model run? | nine historical run cards plus the separately dated PerturbDiff native setup addition; one pinned environment per model; setup remains separate from whole-gene runtime admission | ADR 0021, ADR 0024 |
| [`02_dataset_corpus/`](02_dataset_corpus/README.md) | what data substrate do we build on? | **five source-dataset corpora** — four CM native-axis + `tahoe.h5ad`, 67,018 × 62,710, 25.8 GB — raw counts, every raw well (controls + blanks), no pharmacology filter and no Δ; **seven crosswalk tables** + coverage README (18,164-gene CM axis as `tables/genes.csv`); `--verify` **G1–G6 PASS**, with G4 reproducing 5/5 frozen v1 view content hashes | ADR 0033 (amends ADR 0026) |
| [`03_modes_tenv/`](03_modes_tenv/README.md) | which admissible model predicts tenv with shared Morgan inputs? | final 63-slot comparison across five predictors and two references. PerturbDiff led seven of eight nine-run metric means, including positive/negative AP@100 of 0.23914/0.23463. The training-group mean reference led mean R² | ADR 0042 |
| [`04_pdiff_hyperparameter/`](04_pdiff_hyperparameter/README.md) | how does standalone PerturbDiff train and score on Tahoe and tenv? | complete nine-fit matrices for both datasets, learning curves, eight held-out metrics with matched references, and nine verified tenv EMA checkpoints for handoff | ADR 0046 |
| [`05_pdiff_transfer/`](05_pdiff_transfer/README.md) | does Tahoe training improve tenv prediction? | verified common cold tenv fold with reused CM checkpoints and deterministic references | ADR 0047 |
| [`06_pdiff_transfer_optimize/`](06_pdiff_transfer_optimize/README.md) | what accounts for the Tahoe Transfer deficit? | paired validation replay, independent interventions and source-domain audit | ADR 0048 |

The retired Tahoe data step was folded into `02_dataset_corpus/` (ADR 0033).
The old CM benchmark and PRnet steps 05 and 06 were removed on 2026-11-07. Their tracked records are available in Git at
`f95fe1e`. The retained v1 CM views and identity tables live in
[`_frozen/cm-data-v1/`](_frozen/cm-data-v1/README.md).

## Selection constraints in force

- **Delivery requirement: 15,000+ genes.** The eventual MIRI perturbome must be
  whole-transcriptome, so a 978-landmark model cannot be the delivery model.
- **978 quarantine.** PRnet's native demo is the 978-gene LINCS matrix; it is
  allowed only inside `01_setup/execution/prnet/` and its outputs are labelled
  *setup-only*. No later Stage 2 evaluation, metric, fold or comparison may read
  a 978 artifact (ADR 0024). The whole-transcriptome guard still governs selection.
- The tenv run uses the shared chemical inputs, runtime admission gates and
  frozen cold-drug folds documented in [ADR 0042](../docs/adr/0042-tenv-model-benchmark-redesign.md).
  [The 2026-10-26 audit](audit/model_selection_decisions.md) records the
  earlier CM study's selection decisions.

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
- **`03_modes_tenv/`** — the final held-out tenv tables and figures across
  three cold-drug folds and three training seeds. The step owns the frozen
  18,164-gene tenv input release, model adapters and evaluation artifacts.

- **`04_pdiff_hyperparameter/`** — the two standalone PerturbDiff studies
  selected EMA checkpoints using each fold's validation loss. Both
  datasets have complete held-out scoring. Tenv reused the frozen Step 03
  test memberships. Two Tahoe predictions interrupted by a power cut
  recovered from their sealed fits, and the completed Tahoe evaluation
  replayed byte-for-byte. [Results and recovery](04_pdiff_hyperparameter/README.md#published-results)

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
- ADR [0033](../docs/adr/0033-stage2-dataset-corpus.md) — one dataset-corpus layer, views downstream; the old Tahoe data step is folded in.
- [`audit/model_selection_decisions.md`](audit/model_selection_decisions.md) — the earlier CM study's selection decisions.
- ADR [0042](../docs/adr/0042-tenv-model-benchmark-redesign.md) — the completed tenv comparison.
- ADR [0046](../docs/adr/0046-perturbdiff-standalone-training.md) — the standalone Tahoe and tenv PerturbDiff studies.
- ADR [0047](../docs/adr/0047-perturbdiff-tahoe-transfer-to-tenv.md) — Tahoe training strategies across three tenv split regimes.
- [`_archive/PLAN.md`](_archive/PLAN.md) — the superseded plan, archived as-is; this page carries its live content.

Stage 1 (`stage1/`) produces the MIRI disease inputs and DEG lists Stage 2's
models will serve; it is not a dependency of this stage's setup or benchmark.
