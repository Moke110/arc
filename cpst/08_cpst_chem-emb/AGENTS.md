# Step 08 CarPS-Transcriptome chemical representations

Read the study report and ADR 0099 before changing the roster or protocol.
For model, fit and replay changes, use the shared training and evaluation APIs
through the numbered commands. [Study report](README.md) · [ADR 0099](../../docs/adr/0099-carps-vae-alpha-and-chemical-representations.md) · [Shared APIs](../util/README.md#transigen-training-and-evaluation)

## Prepare frozen inputs

- Use the Step 07 alpha release, the Step 03 tenv Denoise-DEG 5k release and
  the frozen tenv vectors, keyed by desalted parent InChIKey. Fit with the ordered
  gene axis, measured same-plate DMSO, the training-fold dose transform and the
  alpha objective. Calibrate MMD bandwidth separately on the frozen tenv HVG-DEG
  5k training fold. [Alpha guide](../07_cpst_alpha/AGENTS.md) · [Fit release manifest](../03_data_train/03.7_tenv_denoise-5k/output/manifest.json) · [Bandwidth release manifest](../03_data_train/03.1_tenv_hvg-5k/output/manifest.json) · [Representation manifest pins](../util/output/api/cpst_util/representations.py) · [Author source manifest](../../pert_model/TranSiGen/manifest.sha256)
- Vary only the frozen vectors and their first projection. Use each arm's native
  width and a 400-wide projection output. Pair downstream initialization and
  CPU and CUDA RNG state within each fold and seed. The formal roster covers
  seven arms, cold-drug folds 0–2 and training seeds 17, 42 and 2026. Keep the
  validation and held-out sampling seeds at 22617 and 100042. [Coordinate runner](execution/chem_runs.py)
- Run from the materialized repository root with `arc` active, CUDA available
  for training and checkpoint replay, and pinned Git revisions accessible. The
  commands default to the installed inputs. Supply identical relocated inputs
  with `--training`, `--representations`, `--source`, `--bandwidth-training` or
  `--round-four` where the command exposes that flag. [Execution limits](../../docs/agents/execution-environment.md)

Use a fresh scratch directory for each recipe. Run heavy commands within the
12 GiB RAM and 8 GiB swap scope, one CUDA coordinate at a time.

```bash
set -e
STEP="$PWD/cpst/08_cpst_chem-emb"
SAMPLER="$STEP/../06_tsg_shape/06.1_round1/execution/sample_resources.py"
mkdir -p "$PWD/tmp"
WORK=$(mktemp -d "$PWD/tmp/carps-chemical.XXXXXX")
mkdir "$WORK/scratch" "$WORK/matplotlib"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib" PYTHONDONTWRITEBYTECODE=1
capped() {
  LC_ALL=C free -g | awk '$1 == "Mem:" {available=$7} END {exit available<12}' || return
  systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=8G -- \
    env OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 "$@"
}
```

## Run the bounded CUDA gate

For a fit or protocol change, exercise all seven arms on fold 0 and seed 17
for two epochs with 64 train, 24 validation and 24 held-out conditions. Start
`python "$SAMPLER" --output "$WORK/coordinator.jsonl"` in a second terminal,
with the paths from the setup above. Confirm a first ten-second sampler
reading before training and retain measurements through completion.
[Coordinator sampler](../06_tsg_shape/06.1_round1/execution/sample_resources.py) · [Gate receipt](log/363-gate.json)

```bash
capped python "$STEP/../07_cpst_alpha/execution/verify_alpha.py" > "$WORK/alpha-before.json"
capped python "$STEP/execution/chem_matrix.py" --smoke \
  --output "$WORK/matrix" --logs "$WORK/matrix-audit"
for report in report replay-report; do
  capped python "$STEP/execution/chem_report.py" --smoke \
    --matrix "$WORK/matrix" --output "$WORK/$report"
done
diff -r "$WORK/report" "$WORK/replay-report"
python "$STEP/execution/chem_figures.py" build --data "$WORK/report" \
  --output "$WORK/figures" --scratch "$WORK/scratch"
capped python "$STEP/../07_cpst_alpha/execution/verify_alpha.py" > "$WORK/alpha-after.json"
cmp "$WORK/alpha-before.json" "$WORK/alpha-after.json"
```

Require exact independent checkpoint replay for all seven fits, 63 matching
scores across both scoring processes, paired initialization and RNG, two
matching CSV-only renders and identical alpha receipts. Compare the bounded
KPGT files under `matrix/runs/KPGT/fold-0/seed-17` to the Step 06.4 digests in
`bounded_kpgt_matches_step_06_4.identical_files_sha256` in the gate receipt.
Retain the phase and resource records. The formal reuse runner validates the
accepted gate receipt and its alpha and roster identities before fitting.

## Reproduce the formal study

Start the sampler as above. Keep it running until it has recorded at least
one sample after the matrix finishes. Register the full roster before
training. The reference runner tries each KPGT coordinate against the sealed
Step 06.4 ledgers and pinned code, verifies model, input and protocol identity
against the frozen alpha, then independently replays the checkpoint and
rescores all nine metrics and eligibility counts exactly. It saves the reuse
plan and fits any coordinate that fails the gate. For standalone replay with
the original Step 06.4 runner, follow its pinned checkout procedure.
[Reference runner](execution/chem_reference.py) · [Pinned Step 06.4 replay](../06_tsg_shape/06.4_round4/README.md#replay-a-sealed-step-064-fit) · [Formal runner](execution/chem_formal.py)

```bash
mkdir "$WORK/formal"
capped python "$STEP/execution/chem_reference.py" \
  --matrix "$WORK/formal/matrix" --logs "$WORK/formal/matrix-audit" \
  --evidence "$WORK/formal/reuse-preflight.json"
# Wait for the sampler's next reading before building the scientific report.
capped python "$STEP/execution/chem_formal.py" build \
  --bundle "$WORK/formal" --coordinator "$WORK/coordinator.jsonl" \
  --scratch "$WORK/scratch" --publication "$WORK/publication"
```

Require a complete 63-coordinate roster, exact replay for every fresh fit,
exactly matching independent reports for all 567 scores, and a resource audit
that brackets the matrix. `build` uses the historical figure set. Keep the
candidate and its phase logs under `tmp/`.

## Verify and refresh the installed release

Run the formal verifier against the materialized `output/formal` bundle. It
rechecks recorded KPGT reuses, independently scores the complete roster and
audits the paired report, resource logs and canonical figures.
[Formal receipt](log/364-formal.json) · [Fit ledger](log/runs.csv)

```bash
capped python "$STEP/execution/chem_formal.py" verify --scratch "$WORK/scratch"
```

Read `working_default` and compare it to the inherited criterion on the
unrounded paired summary. Change the working default from KPGT only when one
arm uniquely leads both primary metrics, both metrics have nine defined
pairs, and each improvement strictly exceeds its own paired sample SD. Use
mean arm-minus-KPGT Δ PCC and the negative of mean arm-minus-KPGT Δ RMSE as the
two improvements. Follow ADR 0100 for the downstream configuration.
[Paired summary](figure/data/paired-summary.csv) · [ADR 0100](../../docs/adr/0100-carps-vae-unimol.md)

Build the canonical figures from the accepted bundle, matched controls and
registered Step 06.4 reuses. The exporter checks the sealed JAC-top200 scores. For a CSV-only check of installed PNGs, run the
second command by itself. Require two identical renders and matching PNG
hashes. Keep the formal evidence available for Step 09 replay.
[Comparison exporter](execution/build_comparison.py) · [Figure verifier](execution/chem_figures.py) · [Figure receipt](log/figure-render.json) · [ADR 0109](../../docs/adr/0109-cpst-seven-metric-comparison-figures.md)

```bash
capped python "$STEP/execution/build_comparison.py" \
  --output "$WORK/comparison" --scratch "$WORK/scratch"
python "$STEP/execution/chem_figures.py" verify --comparison \
  --data "$STEP/figure/data" --output "$STEP/figure" --scratch "$WORK/scratch"
```
