# 03_modes_tenv — tenv shared-Morgan comparison

This step implements [specification #40](https://github.com/Moke110/arc/issues/40)
under [ADR 0042](../../docs/adr/0042-tenv-model-benchmark-redesign.md).
The requested spelling `03_modes_tenv` is intentional. Historical benchmark
artifacts are sealed under `stage2/_archive/03_models_cmbench/`; their evidence
and the `_frozen/cm-data-v1/` bridge are retained.

The input release projects tenv raw counts to the ordered CM gene axis before
CP10K/log1p normalization. Same-plate DMSO means are allowed baseline inputs.
Training and validation targets are Δ from those controls; test-treated Δ and
directional labels live only under `fold_N/evaluator/`. Predictor manifests
contain chemical, dose, case and approved control paths only. Each group's
well predictions and truth are averaged before equally weighted group scoring.

All entrants share chirality-aware radius-2 binary Morgan-2048 vectors and three
parent-connectivity-block folds (split seed 20261105), crossed with training
seeds 17, 42 and 2026. Model admission requires separate runtime evidence on all
18,164 genes; the discovery ledger does not establish successful admission.

Build and verify inputs from the repository root using the `arc` interpreter:

```bash
python stage2/03_modes_tenv/execution/script/inputs.py build
python stage2/03_modes_tenv/execution/script/inputs.py verify
python -m unittest discover -s stage2/03_modes_tenv/execution/tests -p test_inputs.py -v
```

`output/input/manifest.json` records exact membership counts, source/config hashes,
RDKit version, structure normalization and all release file hashes. Arrays use
float32 `.npy` files for bounded memory mapping. Source corpora and concluded
pretraining/transfer results remain owned by their existing steps.

Frozen input evidence (#42): **460 drugs, 4,206 treated wells, 1,380 drug–dose
groups and 18,164 genes**. Controls separate into 192 DMSO, 80 media and 32
blank rows. Test drug counts are 154/153/153; each validation set has 46 drugs;
training counts are 260/261/261. Two complete builds have identical manifests
and array hashes. Three synthetic CLI tests pass, including normalization
order and malformed-control rejection. Build peaks and commands are in
`execution/input_record/verification.json`; arrays remain local-only.
