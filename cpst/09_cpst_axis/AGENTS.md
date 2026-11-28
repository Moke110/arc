# Step 09 CarPS-Transcriptome prediction axes

Read the study report and ADR 0101 before changing the roster or fit protocol.
Apply ADR 0102 for full-axis scoring, figure publication and axis selection.
Use the shared training and evaluation APIs for model, fit and replay changes.
[Study report](README.md) · [ADR 0101](../../docs/adr/0101-carps-vae-prediction-axes.md) · [ADR 0102](../../docs/adr/0102-stage2-axis-selection-and-comparison-figures.md) · [Shared APIs](../util/README.md#transigen-training-and-evaluation)

## Prepare frozen inputs

- Use all eight axis-local, plate-pooled tenv training views with the same
  raw membership, exclusions and ordered cold-drug partitions. Step 03 owns
  building and verifying the 2k views. Preserve the ordered genes and the
  measured same-plate control on each axis.
  [Axis catalog](../02_gene_axis/README.md) · [Step 03 guide](../03_data_train/AGENTS.md)
- Use the fixed Step 07 alpha, the Step 08 UniMol vectors and the sealed Step
  08 fit evidence for `denoise-5k`. Materialize the representation bundle at
  `cpst/util/output/representations` and the author TranSiGen source at
  `pert_model/TranSiGen/source`. For standalone replay of the alpha lineage,
  also materialize the Step 06.4 formal bundle.
  [Alpha guide](https://github.com/Moke110/carps/blob/a7613b74fc23f9f48d8ea0d4479542c008d701db/cpst/07_cpst_alpha/AGENTS.md) · [UniMol decision](../../docs/adr/0100-carps-vae-unimol.md) · [Representation pins](../util/output/api/cpst_util/representations.py) · [Author source manifest](../../pert_model/TranSiGen/manifest.sha256) · [Step 08 receipt](../08_cpst_chem-emb/log/364-formal.json)
- Keep the frozen UniMol projection, internal model widths, loss and fit policy.
  Vary the input widths of the control and measured Δ encoders and the output
  widths of the control and shared Δ decoders. Pair the other initialized
  tensors and the training random streams by cold-drug fold and seed. Keep the
  training-fold dose transform and the fixed validation and held-out sampling.
  [Axis adapter](execution/axis_model.py) · [Coordinate runner](execution/axis_runs.py)
- For scoring, materialize Step 04's `output/report-matched-control/bandwidths.csv`.
  Verify its digest and independently reconstruct its calibration on the
  HVG-DEG 5k training rows. Use one bandwidth per fold across axes and seeds.
  Score held-out conditions on the full prediction axis under the Step 04
  eligibility rules. Pair defined scores with `denoise-5k` within each fold
  and seed, then summarize coordinates with equal weight and sample SD.
  [Step 04 receipt](../04_tsg_axis/log/319-formal.json) · [Scorer](execution/axis_report.py) · [Independent audit](execution/audit_axis_report.py)

Run from the materialized repository root with `arc` active, CUDA available for
fits and replay, and pinned Git revisions accessible. For a linked worktree, use
`--step03`, `--representations` and `--source` for the matrix, reference and
formal fit commands, `--step08` for reference and formal commands, and
`--step03` and `--step04` for scoring, build and verification. Use the original
provenance when reconstructing the ignored payloads.
[Migration receipt](../util/log/385-release-migration.json)

## Set up a measured run

Check live Windows host RAM, WSL RAM, GPU and disk. Run heavy phases serially.
Use a fresh scratch directory for each recipe. For the bounded and formal fit
recipes, start the coordinator sampler in a second terminal with the same
`WORK` path. Confirm its first sample before the matrix starts and keep it
running through the first post-matrix sample. Retain the phase timings and
sampler records with the candidate.
[Execution limits](../../docs/agents/execution-environment.md) · [Coordinator sampler](../util/training/sample_resources.py)

```bash
set -e
STEP="$PWD/cpst/09_cpst_axis"
FORMAL="$STEP/execution/axis_formal.py"
SAMPLER="$STEP/../util/training/sample_resources.py"
SOURCE="$PWD/pert_model/TranSiGen/source"
mkdir -p "$PWD/tmp"
WORK=$(mktemp -d "$PWD/tmp/carps-axis.XXXXXX")
mkdir "$WORK/scratch" "$WORK/matplotlib"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib" LC_ALL=C
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2
capped() {
  free -g | awk '$1 == "Mem:" {available=$7} END {exit available<12}' || return
  systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=8G -- "$@"
}
```

Run this command in the coordinator terminal after setting the same `WORK`
and `SAMPLER` paths.

```bash
python "$SAMPLER" --output "$WORK/coordinator.jsonl"
```

## Run the bounded CUDA gate

Exercise all eight axes on fold 0 and seed 17. Replay the nine Step 07 alpha
coordinates and the nine sealed Step 08 UniMol reference coordinates. The
preflight binds all 72 formal coordinates and checks the release and partition
identities. The gate requires exact independent replay, paired epoch RNG,
identical initialized state at equal widths, full-axis scores, two matching
CSV-only renders and coordinator coverage.
[Matrix runner](execution/axis_matrix.py) · [Reference audit](execution/axis_reference.py) · [Gate verifier](execution/verify_axis_gate.py) · [Accepted gate](log/367-bounded.json)

```bash
capped /usr/bin/time -v -o "$WORK/alpha.time" \
  env PYTHONPATH="$STEP/.." python -m util.training.cpst_reference --source "$SOURCE" \
  > "$WORK/alpha.json"
capped /usr/bin/time -v -o "$WORK/gate.time" \
  python "$STEP/execution/axis_matrix.py" --bounded --source "$SOURCE" \
  --output "$WORK/matrix" --logs "$WORK/matrix-audit"
capped /usr/bin/time -v -o "$WORK/reference.time" \
  python "$STEP/execution/axis_reference.py" --source "$SOURCE" \
  --step08 "$STEP/../08_cpst_chem-emb" > "$WORK/reference.json"
capped /usr/bin/time -v -o "$WORK/score.time" \
  python "$STEP/execution/axis_report.py" --bounded --matrix "$WORK/matrix" \
  --output "$WORK/report" --figure-data "$WORK/figures/data"
capped /usr/bin/time -v -o "$WORK/figures.time" \
  python "$STEP/execution/axis_figures.py" build --data "$WORK/figures/data" \
  --output "$WORK/figures" --scratch "$WORK/scratch" > "$WORK/figures.json"
capped python "$STEP/execution/verify_axis_gate.py" --work "$WORK" \
  --coordinator "$WORK/coordinator.jsonl" --source "$SOURCE" --output "$WORK/gate"
```

Require a `PASS` receipt. The formal runner checks the accepted gate and its
companion digests under `log/` before accepting a formal preflight. Keep the
fresh gate evidence in scratch for review.

## Reproduce and publish the formal study

Register all eight axes across folds 0–2 and seeds 17, 42 and 2026. The formal
runner saves a reuse decision for each Step 08 reference coordinate. Reuse requires
matching input and protocol provenance, initialization and exact validation and
held-out checkpoint replay. Fit any rejected reference coordinate and retain
its rejection reason. Build a candidate after the first post-matrix sampler
reading. `build` independently audits the scores, pairings, resource evidence
and CSV-only figures before writing a publication.
[Formal runner](execution/axis_formal.py) · [Accepted reuse plan](log/368-reference.json)

```bash
capped /usr/bin/time -v -o "$WORK/fit.time" \
  python "$FORMAL" fit --source "$SOURCE" --output "$WORK/bundle"
# Continue once the coordinator records the first post-matrix sample.
capped python "$FORMAL" build --bundle "$WORK/bundle" \
  --coordinator "$WORK/coordinator.jsonl" --publication "$WORK/publication" \
  --scratch "$WORK/scratch"
capped /usr/bin/time -v -o "$WORK/verify.time" \
  python "$FORMAL" verify --bundle "$WORK/bundle" \
  --publication "$WORK/publication" --scratch "$WORK/scratch"
```

Require exact replay evidence for all coordinates, the complete 72-coordinate
full-axis report and three independently reproduced figures. After acceptance,
install the verified bundle at the local `output/formal` path. If that path
already holds a different fit identity, use an unused `--destination` for the
candidate and preserve the accepted bundle through publication. Publish the
candidate's `log/` and `figure/` files together under the step. Keep the local
bundle available for later verification.
[Formal receipt](log/371-formal.json) · [Figure receipt](log/figure-render.json)

```bash
capped python "$FORMAL" install --bundle "$WORK/bundle" \
  --publication "$WORK/publication" --destination "$STEP/output/formal"
```

## Verify the installed release

Use the original checkout at `a7613b74fc23f9f48d8ea0d4479542c008d701db`
with the materialized `output/formal` bundle for the full scientific rescore. The pinned guide contains the command
and the sealed report is validated against its original tracked receipts.
The active figure receipt records the projection of those accepted score tables.
Run the active CSV-only figure check from this checkout.
[Historical verification](https://github.com/Moke110/carps/blob/a7613b74fc23f9f48d8ea0d4479542c008d701db/cpst/09_cpst_axis/AGENTS.md#verify-the-installed-release) · [Active figure receipt](log/figure-render.json)

```bash
python "$STEP/execution/axis_figures.py" verify --data "$STEP/figure/data" \
  --output "$STEP/figure" --scratch "$WORK/scratch"
```

## Rebuild the report from sealed fits

For a report-only change, point a scratch bundle at the accepted matrix,
phase audits and reuse plan. Build with the original coordinator samples. Use
the same candidate verification and installation commands with this `WORK`
path. The installer verifies the fit identity before replacing the report
bound to the recognized predecessor digest.
[Original coordinator samples](log/368-coordinator-resources.jsonl)

```bash
SEALED="$STEP/output/formal"
mkdir "$WORK/bundle"
for item in matrix matrix-audit reuse-preflight.json; do
  ln -s "$SEALED/$item" "$WORK/bundle/$item"
done
capped python "$FORMAL" build --bundle "$WORK/bundle" \
  --coordinator "$STEP/log/368-coordinator-resources.jsonl" \
  --publication "$WORK/publication" --scratch "$WORK/scratch"
```
