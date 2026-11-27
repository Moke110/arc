# Step 08 CarPS-VAE chemical representations

[ADR 0099](../../docs/adr/0099-carps-vae-alpha-and-chemical-representations.md) compares ECFP4C, ECFP6C, FCFP4C, FCFP6C, E3FP, KPGT and UniMol on the fixed [CarPS-VAE alpha](../07_cpsv_alpha/README.md). Each arm uses the same 1,383 plate-pooled tenv conditions from 459 parent compounds, the Denoise-DEG 5k axis, measured matched DMSO, and alpha dose handling and objective. Only the frozen chemical vectors and their first projection vary. Native widths range from 512 to 2,304 and all project to 400. The other initialized weights and CPU and CUDA random streams are paired within each fold and seed.

## Bounded CUDA gate

The gate checks all seven arms with cold-drug fold 0 and training seed 17, two epochs and the first 64 training, 24 validation and 24 held-out conditions from the inherited protocol. The formal preflight roster names all 63 coordinates across folds 0–2 and seeds 17, 42 and 2026. The [gate receipt](log/363-gate.json) records resource limits, independent train/predict/replay phases, inherited scoring and CSV-only figure verification. The KPGT gate also exactly matches the Step 06.4 bounded run’s validation and held-out predictions, condition and gene order, signed Δ arrays and traces. Validation and held-out sampling use seeds 22617 and 100042 respectively.

Run from the materialized checkout with the `arc` environment active and at least 12 GiB of WSL RAM available. Start the [ten-second coordinator sampler](../06_tsg_shape/06.1_round1/execution/sample_resources.py) before training and keep it running through completion. Follow the [execution limits](../../docs/agents/execution-environment.md).

```bash
set -e
STEP="$PWD/stage2_CarPS-VAE/08_cpsv_chem-emb"
mkdir -p "$PWD/tmp"
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

## Formal results

The accepted matrix covers all 63 coordinates. It reuses 9 verified
KPGT fits from Step 06.4 and exactly replays 54 fresh fits.
Each point is one cold-drug fold and training-seed coordinate, and each mean gives equal
weight to the nine coordinates. Uncertainty is sample SD across the repeated coordinates.

| Arm | Δ PCC | Δ RMSE | Paired Δ PCC vs KPGT | Paired Δ RMSE vs KPGT |
| --- | ---: | ---: | ---: | ---: |
| ECFP4C | 0.424806 ± 0.014891 | 0.242734 ± 0.001844 | -0.004024 ± 0.019202 | +0.001982 ± 0.003108 |
| ECFP6C | 0.421715 ± 0.016410 | 0.242380 ± 0.001609 | -0.007115 ± 0.015305 | +0.001628 ± 0.001952 |
| FCFP4C | 0.428691 ± 0.024270 | 0.240448 ± 0.001044 | -0.000139 ± 0.019621 | -0.000304 ± 0.002016 |
| FCFP6C | 0.423095 ± 0.013477 | 0.241313 ± 0.001584 | -0.005735 ± 0.016716 | +0.000560 ± 0.002741 |
| E3FP | 0.435100 ± 0.021071 | 0.238626 ± 0.002190 | +0.006270 ± 0.014760 | -0.002126 ± 0.002559 |
| KPGT | 0.428830 ± 0.018574 | 0.240752 ± 0.002598 | +0.000000 ± 0.000000 | +0.000000 ± 0.000000 |
| UniMol | 0.441262 ± 0.013620 | 0.238599 ± 0.003729 | +0.012432 ± 0.013493 | -0.002153 ± 0.002016 |

UniMol had the highest mean Δ PCC (0.441262).
UniMol had the lowest mean Δ RMSE (0.238599).
UniMol’s mean paired Δ PCC gain was smaller than its sample SD (0.012432 versus 0.013493). The comparison is closed, and [ADR 0100](../../docs/adr/0100-carps-vae-unimol.md) selects UniMol for subsequent CarPS-VAE studies.

- [All nine metrics by fold and seed](figure/data/fold-scores.csv), [foldwise paired summaries](figure/data/paired-by-fold.csv) and [overall paired summaries](figure/data/paired-summary.csv)
- [Nine held-out metrics](figure/delta-metrics.png), [nine paired differences vs KPGT](figure/delta-metrics-paired-difference.png) and [CSV-only figure replay receipt](log/figure-render.json)
- [Fit and checkpoint ledger](log/runs.csv), [epoch objectives](figure/data/epoch-objectives.csv) and [bandwidths](log/bandwidths.csv)
- [Formal acceptance receipt](log/364-formal.json), [matrix resources](log/364-resources.jsonl) and [coordinator resources](log/364-coordinator-resources.jsonl)

## Formal reproduction

Run the formal comparison under the same resource limits and environment as
the bounded gate. The formal runner checks each KPGT reference separately
against the frozen Step 07 release, the Step 06.4 ledgers, a fresh exact
checkpoint replay and independently rescored held-out metrics. It checks
the inherited loss and fit modules against the original code hashes in the
alpha release. It registers a new fit for any rejected coordinate.
Each new fit has separate training,
prediction and exact replay processes.

Use the `capped` function and `STEP` defined above. Set `COORDINATOR` to the
location of the running ten-second sampler’s JSONL file. Start the sampler
before the matrix and keep it running until at least one sample after the
matrix completes.

```bash
mkdir -p "$PWD/tmp"
WORK=$(mktemp -d "$PWD/tmp/carps-chemical-formal.XXXXXX")
mkdir "$WORK/scratch" "$WORK/formal" "$WORK/matplotlib"
export TMPDIR="$WORK/scratch" MPLCONFIGDIR="$WORK/matplotlib" PYTHONDONTWRITEBYTECODE=1
capped python "$STEP/execution/chem_reference.py" \
  --matrix "$WORK/formal/matrix" --logs "$WORK/formal/matrix-audit" \
  --evidence "$WORK/formal/reuse-preflight.json"
sleep 10
capped python "$STEP/execution/chem_formal.py" build \
  --bundle "$WORK/formal" --coordinator "$COORDINATOR" \
  --scratch "$WORK/scratch" --publication "$WORK/publication"
```

The comparison publishes one score per metric, arm, fold and training
seed. The paired differences subtract the matching KPGT score. Higher
Δ PCC and lower Δ RMSE indicate improvement. Each summary weights
its coordinates equally and reports sample SD. The original Step 08
criterion for changing the working configuration required the same
representation to lead both primary metrics and both paired improvements
to exceed their sample SD across the nine pairs. The successor decision
is recorded in [ADR 0100](../../docs/adr/0100-carps-vae-unimol.md).

To verify the installed release from the materialized checkout, use the
canonical bundle in `output/formal` and pass a fresh scratch directory under
`tmp/`. This replays the retained KPGT predictions, independently rescores
the complete matrix, audits the arithmetic and verifies the PNGs using
only their tracked CSV.

```bash
capped python "$STEP/execution/chem_formal.py" verify --scratch "$WORK/scratch"
```
