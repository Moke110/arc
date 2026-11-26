# Step 08 CarPS-VAE chemical representations

[ADR 0099](../../docs/adr/0099-carps-vae-alpha-and-chemical-representations.md) compares ECFP4C, ECFP6C, FCFP4C, FCFP6C, E3FP, KPGT and UniMol on the fixed [CarPS-VAE alpha](../07_cpsv_alpha/README.md). Each arm uses the same 1,383 plate-pooled tenv conditions from 459 parent compounds, the Denoise-DEG 5k axis, measured matched DMSO, and alpha dose handling and objective. Only the frozen chemical vectors and their first projection vary. Native widths range from 512 to 2,304 and all project to 400. The other initialized weights and CPU and CUDA random streams are paired within each fold and seed.

## Bounded CUDA gate

The gate checks all seven arms with cold-drug fold 0 and training seed 17, two epochs and the first 64 training, 24 validation and 24 held-out conditions from the inherited protocol. The formal preflight roster names all 63 coordinates across folds 0–2 and seeds 17, 42 and 2026. The [gate receipt](log/363-gate.json) records resource limits, independent train/predict/replay phases, inherited scoring and CSV-only figure verification. The KPGT gate also exactly matches the Step 06.4 bounded run’s validation and held-out predictions, condition and gene order, signed Δ arrays and traces. Validation and held-out sampling use seeds 22617 and 100042 respectively.

Run from the materialized checkout with the `arc` environment active and at least 12 GiB of WSL RAM available. Start the [ten-second coordinator sampler](../06_tsg_shape/06.1_round1/execution/sample_resources.py) before training and keep it running through completion. Follow the [execution limits](../../docs/agents/execution-environment.md).

```bash
set -e
STEP="$PWD/stage2_CarPS-VAE/08_cpsv_chem-emb"
WORK=$(mktemp -d "$PWD/tmp/carps-chemical-gate.XXXXXX")
mkdir "$WORK/scratch"
export TMPDIR="$WORK/scratch" PYTHONDONTWRITEBYTECODE=1
capped() {
  LC_ALL=C free -g | awk '$1 == "Mem:" {available=$7} END {exit available<12}' || return
  systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=8G -- \
    env OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 "$@"
}
capped python "$STEP/../07_cpsv_alpha/execution/verify_alpha.py" > "$WORK/alpha-before.json"
capped python "$STEP/execution/chem_matrix.py" --smoke \
  --output "$WORK/matrix" --logs "$WORK/matrix-audit"
for report in report replay-report; do
  capped python "$STEP/execution/chem_report.py" --smoke \
    --matrix "$WORK/matrix" --output "$WORK/$report"
done
diff -r "$WORK/report" "$WORK/replay-report"
python "$STEP/execution/chem_figures.py" build --data "$WORK/report" \
  --output "$WORK/figures" --scratch "$WORK/scratch"
capped python "$STEP/../07_cpsv_alpha/execution/verify_alpha.py" > "$WORK/alpha-after.json"
cmp "$WORK/alpha-before.json" "$WORK/alpha-after.json"
```

The matrix, phase and scoring commands accept alternate paths for the same frozen training release, vector manifest and author source. Scoring also accepts the frozen HVG-DEG 5k bandwidth training release path. The registration captures content identities and ordered partitions so accepted runs can be replayed after moving the checkout.
