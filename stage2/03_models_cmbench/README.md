# 03_models_cmbench — Stage 2 cold-drug benchmark (two substeps)

The executed benchmark is a **60-unit cold-drug run**: five CPU baselines plus
seven real models × five cold-drug folds, all scored against one fixed
group/fold/metric contract (ADR 0027). The step is named for the work it does,
not for the dataset, and splits into two independently-readable substeps:

| substep | question it answers | headline (k=200 macro AUPRC) | contract |
|---|---|---|---|
| [`03.1_cold_drug_18k/`](03.1_cold_drug_18k/README.md) | Can any model beat the mean baseline at 18,164 genes? | `delta_mean` **0.1989**, best entrant PerturbNet 0.1861 | ADR 0027 |
| [`03.2_hvg_sensitivity_5k/`](03.2_hvg_sensitivity_5k/README.md) | Does that ranking depend on the gene axis? | PRnet **0.2735**, ahead of `delta_mean` 0.2491 | ADR 0028 |

**Conclusion.** At the full transcriptome no entrant beats the constant
training-mean baseline; on a 5,000-gene HVG axis PRnet moves above it (+0.0243
paired). The 18k conclusion is therefore **gene-space dependent**, and the 5k
run is a sensitivity result, not a new winner: absolute AUPRC is not comparable
across the two runs (k=200 prevalence 1.1% vs 4.0%) — compare paired gaps only.
Squidiff failed to converge on both axes and its rows are documented failures.

Selection follow-up (2026-10-26): Squidiff is excluded from further
selection/tuning; the remaining six candidates and the preservation policy are
in [`../audit/model_selection_decisions.md`](../audit/model_selection_decisions.md).

## Layout

```
03_models_cmbench/
├── AGENTS.md  README.md              # this step: routing + the two substeps
├── 03.1_cold_drug_18k/
│   ├── AGENTS.md  README.md  summary.png  summary_table.csv
│   ├── execution/{script/, configs/, models/, tests/, Makefile}
│   └── output/{input/, units/, results/}
└── 03.2_hvg_sensitivity_5k/          # identical shape
```

`execution/` is everything an agent runs (harness, registry, run cards, tests);
`output/input/` is the frozen benchmark input a rerun must read, `output/units/`
the per-unit fit/predict artifacts, `output/results/` the evaluated tables. The
human figure and table sit at substep level, per the ADR 0030 step contract.

The two substeps duplicate the harness **on purpose**: each must be reproducible
standalone, and each resolves the repository root by walking up to `.git`
(never by counting directory levels). Their only shared dependencies are the
CM corpus (`stage2/02_cm_data/output/`) and the per-model environments
(`stage2/_envs/`).

## Reproduce

Each substep is three commands from the repository root:

```bash
# 1. offline contract suite (data hashes, splits, metrics, adapter wiring)
python -m unittest discover -s stage2/03_models_cmbench/03.1_cold_drug_18k/execution/tests -v

# 2. regenerate the human summary from the evaluated results
python stage2/03_models_cmbench/03.1_cold_drug_18k/execution/script/make_summary.py

# 3. independent acceptance checks over the 60 units
python stage2/03_models_cmbench/03.1_cold_drug_18k/execution/script/validate_run.py
```

Swap `03.1_cold_drug_18k` for `03.2_hvg_sensitivity_5k` for the sensitivity
run. A full rerun is `bash <substep>/execution/script/run_queue.sh` (GPU, hours);
the frozen input rebuilds with `<substep>/execution/script/build_data.py`.

## Notes

- Frozen evidence keeps historical path strings: each substep's
  `output/input/manifest.json` (`build_command`) is the 2026-10 pre-restructure
  form, and `03.2`'s `"benchmark": "temp_env_all_5k"` is the historical
  benchmark id. The dataset identifiers (`temp_env_all`, `temp_env_all-smiles`,
  `temp_env_all_5k`) are data values and are never renamed.
- `output/units/` is regenerable and gitignored; `output/results/` keeps the
  conclusion tables (leaderboard, per-drug, fidelity, run cards, bootstrap).
