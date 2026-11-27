# Step 10 CarPS-Transcriptome dose-input comparison

Read the study report and ADR 0103 before changing the model or protocol.
Use the shared training and evaluation APIs for fits and replay, and follow
ADR 0109 for comparison figures. The Step 07 alpha release and accepted
comparison bundles are immutable inputs.
[Study report](README.md) · [Dose design](../../docs/adr/0103-carps-vae-dose-inputs.md) · [Figure contract](../../docs/adr/0109-cpst-seven-metric-comparison-figures.md) · [Shared APIs](../util/README.md#transigen-training-and-evaluation)

## Prepare the inputs

- Use the frozen ordered Denoise-DEG 5k axis, measured same-plate DMSO and the
  three `cold_dose_drug` folds. Test the intersection of held parent and held
  dose. Train and validate on other parents at the other doses. Exclude both
  single-holdout strips.
  [Training release manifest](../03_data_train/03.7_tenv_denoise-5k/output/manifest.json) · [Step 03 guide](../03_data_train/AGENTS.md) · [Input verifier](execution/dose_inputs.py)
- Use the Step 07 `XΔ-XΔ` backbone and its selected objective with frozen UniMol
  for the nine CarPS-Transcriptome arms. The TranSiGen reference uses the
  author `XX-XX` architecture and objective, frozen KPGT and no dose.
  Materialize the representation bundle and the author source before fitting.
  Match same-shaped initialized parameters and training random streams across
  CarPS-Transcriptome arms within fold and seed.
  [Alpha release](../07_cpst_alpha/output/release.json) · [Representation pins](../util/output/api/cpst_util/representations.py) · [Author source manifest](../../pert_model/TranSiGen/manifest.sha256) · [Protocol](execution/dose_protocol.py) · [Model adapter](execution/dose_model.py)
- Select each checkpoint by minimum validation objective. Require exact
  independent replay of the validation objective and both validation and
  held-out predictions. For scoring, subtract measured same-plate DMSO from
  the treated predictions. Score full-axis Δ PCC and Δ RMSE and calibrate MMD
  on each fold's complete training partition. Share the bandwidth across
  arms and seeds, and pair CarPS-Transcriptome arms against no-dose within
  fold and seed. Summarize defined coordinate values with equal weight and
  sample SD.
  [Coordinate runner](execution/dose_runs.py) · [Scorer](execution/dose_report.py) · [Independent audit](execution/audit_dose_report.py)

## Set up a measured run

Run from the checkout root with the pinned `arc` environment active. Use CUDA
for fit and replay, and keep the pinned Git revisions accessible for sealed
evidence verification. For a linked worktree, set `MATERIALIZED` to the
checkout with the Step 03 release, frozen representations and TranSiGen source.
Check live WSL and Windows host RAM, GPU and disk before launching a matrix.
Run the heavy phases serially. Use a fresh scratch directory for each recipe.
[Execution limits](../../docs/agents/execution-environment.md) · [Migration receipt](../util/log/385-release-migration.json)

```bash
set -e
STEP="$PWD/cpst/10_cpst_dose"
EXECUTION="$STEP/execution"
SAMPLER="$PWD/cpst/06_tsg_shape/06.1_round1/execution/sample_resources.py"
MATERIALIZED="${MATERIALIZED:-$PWD}"
STEP03="$MATERIALIZED/cpst/03_data_train"
REPRESENTATIONS="$MATERIALIZED/cpst/util/output/representations"
SOURCE="$MATERIALIZED/pert_model/TranSiGen/source"
mkdir -p "$PWD/tmp"
WORK=$(mktemp -d "$PWD/tmp/step10.XXXXXX")
mkdir "$WORK/scratch" "$WORK/matplotlib"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib" LC_ALL=C
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2
capped() {
  free -g | awk '$1 == "Mem:" {available=$7} END {exit available<12}' || return
  systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=8G -- "$@"
}
```

For each CUDA matrix, start the coordinator sampler in a second terminal with
the same `WORK` and `SAMPLER` paths. Wait for its first sample before starting
preflight and keep it active until the first sample after the measured matrix
command exits. The gate checks the ten-second cadence and coverage.
[Coordinator sampler](../06_tsg_shape/06.1_round1/execution/sample_resources.py)

```bash
python "$SAMPLER" --output "$WORK/coordinator.jsonl"
```

## Run the bounded CUDA gate

Exercise all ten arms on fold 0 and seed 17 for two epochs using the first
64 training, 24 validation and 24 test conditions. After the post-matrix
coordinator sample, score the sealed fits, reproduce the three figures twice
from their plotted CSVs and require a `PASS` gate receipt.
[Matrix runner](execution/dose_matrix.py) · [Gate assembler](execution/dose_gate.py) · [Accepted gate](log/377-bounded.json)

```bash
capped /usr/bin/time -v -o "$WORK/command.resources" \
  python "$EXECUTION/dose_matrix.py" --bounded --step03 "$STEP03" \
  --representations "$REPRESENTATIONS" --source "$SOURCE" \
  --output "$WORK/matrix" --logs "$WORK/matrix-audit"
# Continue after the coordinator records the first post-matrix sample.
capped /usr/bin/time -v -o "$WORK/report.resources" \
  python "$EXECUTION/dose_report.py" --bounded --step03 "$STEP03" \
  --matrix "$WORK/matrix" --output "$WORK/report" --figure-data "$WORK/figure/data"
capped /usr/bin/time -v -o "$WORK/figure.resources" \
  python "$EXECUTION/dose_figures.py" build --data "$WORK/figure/data" \
  --output "$WORK/figure" --scratch "$WORK/scratch" > "$WORK/figure-replay.json"
capped python "$EXECUTION/dose_gate.py" --bundle "$WORK" --output "$WORK/377-bounded.json"
```

Keep the scratch receipt with the bounded evidence. Formal preflight reads
the accepted `log/377-bounded.json` and checks its configuration digest against
the historical configuration identity reconstructed for the current roster.

## Reproduce and verify the formal study

Set up a fresh `WORK` and start a new sampler. Fit all ten arms across folds
0–2 and seeds 17, 42 and 2026. After the post-matrix sample, build and verify
the candidate. The publisher scores in a fresh process and independently
audits all held-out Δ values, metric eligibility, fold and seed statistics,
and plotted values against the report and sealed learning traces. Numeric
audit tolerances are 10⁻¹² absolute and relative. Figure replay is exact.
[Formal publisher](execution/dose_formal.py) · [Accepted formal receipt](log/378-formal.json)

```bash
capped /usr/bin/time -v -o "$WORK/command.resources" \
  python "$EXECUTION/dose_matrix.py" --step03 "$STEP03" \
  --representations "$REPRESENTATIONS" --source "$SOURCE" \
  --output "$WORK/matrix" --logs "$WORK/matrix-audit"
# Continue after the coordinator records the first post-matrix sample.
capped python "$EXECUTION/dose_formal.py" build --bundle "$WORK" \
  --publication "$WORK/publication" --scratch "$WORK/scratch" --step03 "$STEP03"
capped /usr/bin/time -v -o "$WORK/verify.resources" \
  python "$EXECUTION/dose_formal.py" verify --bundle "$WORK" \
  --publication "$WORK/publication" --scratch "$WORK/scratch" --step03 "$STEP03"
```

Require `PASS` for all 90 exact checkpoint replays, all 900 scores and all three
CSV-only figures before accepting a publication. Install the verified bundle
at `output/formal` when that path is free. Use an unused `--destination` for a
candidate with a different identity while the accepted bundle remains active.
Keep the bundle and install destination on the same filesystem. The installer
hardlinks immutable fit files, copies the phase audit and report tables, and
installs them atomically after checking the passing receipts. Publish the
candidate's `figure/` and `log/` files together under this step.
[Figure receipt](log/figure-render.json)

```bash
capped python "$EXECUTION/dose_formal.py" install --bundle "$WORK" \
  --publication "$WORK/publication" --destination "$STEP/output/formal"
```

## Rebuild a report from sealed fits

For report-only work, use a fresh scratch bundle and the original coordinator
samples. Point the bundle at the accepted matrix and per-phase audits, then
use the same formal `build` and `verify` commands with this `WORK`. The publisher
writes the rebuilt report to scratch.
[Accepted coordinator samples](log/378-coordinator-resources.jsonl)

```bash
SEALED="$MATERIALIZED/cpst/10_cpst_dose/output/formal"
ln -s "$SEALED/matrix" "$WORK/matrix"
ln -s "$SEALED/matrix-audit" "$WORK/matrix-audit"
cp "$STEP/log/378-coordinator-resources.jsonl" "$WORK/coordinator.jsonl"
```

## Verify the installed release

Use the original checkout at `a7613b74fc23f9f48d8ea0d4479542c008d701db`
with the materialized `output/formal` bundle for the full scientific rescore. The pinned guide contains the command
and the sealed report is validated against its original tracked receipts.
The active figure receipt records the projection of those accepted score tables.
Run the active CSV-only figure check from this checkout.
[Historical verification](https://github.com/Moke110/carps/blob/a7613b74fc23f9f48d8ea0d4479542c008d701db/cpst/10_cpst_dose/AGENTS.md#verify-the-installed-release) · [Active figure receipt](log/figure-render.json)

```bash
python "$EXECUTION/dose_figures.py" verify --data "$STEP/figure/data" \
  --output "$STEP/figure" --scratch "$WORK/scratch"
```
