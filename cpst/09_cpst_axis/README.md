# Step 09 CarPS-Transcriptome prediction axes

**Status** Denoise-DEG 5k selected. See [ADR 0102](../../docs/adr/0102-stage2-axis-selection-and-comparison-figures.md) and the [full-axis acceptance receipt](log/371-formal.json).

Compare CarPS-Transcriptome across eight [Step 02 prediction axes](../02_gene_axis/README.md), using the [Step 08 UniMol representation](../../docs/adr/0100-carps-vae-unimol.md) and pairing each arm with `denoise-5k` by cold-drug fold and training seed. Use folds 0–2 and seeds 17, 42 and 2026, for 72 axis × fold × seed coordinates.

[Step 03](../03_data_train/README.md) owns the two additional 2k training views and the six established views. Pool matched-plate counts and normalize on each selected axis with the same method. Verify common condition keys, raw membership, exclusions and frozen partition order. Shared genes can have different normalized values across axes.

Resize the input projections of both expression encoders and the output projections of the control and shared Δ decoders. Keep the internal widths and architecture, UniMol input and projection, matched-control protocol, fold-local dose transform, loss and fit policy fixed. Use common initialized values for same-shaped parameters and a paired training random stream. Anchor pairing to the Step 08 UniMol reference so equivalent reference fits can be reused. Verify each reference coordinate's provenance, protocol equivalence and exact checkpoint replay. With all nine reference coordinates accepted, 63 coordinates require new training. Freshly fit any rejected reference coordinate.

## Formal results

**Denoise-DEG 5k is the selected CarPS-Transcriptome prediction axis with UniMol.** It had the highest mean full-axis Δ PCC (0.4413 ± 0.0136) and lowest mean full-axis Δ RMSE (0.2386 ± 0.0037). Denoise-DEG 4k came closest in Δ PCC. Its paired changes were −0.0017 ± 0.0083 in Δ PCC and +0.0249 ± 0.0012 in Δ RMSE. Axis-local normalization gives shared genes different numerical targets across axes.

The matrix covers 72 coordinates. All nine Step 08 UniMol Denoise-DEG 5k fits passed the [reuse audit](log/368-reference.json). The accepted evidence records exact replay of validation and held-out predictions for the 63 new and nine reused fits. The table shows the equal-weight mean ± sample SD over nine fold and seed coordinates per axis. [All scores and paired changes](log/summary.csv) · [Fit lineage](log/runs.csv) · [MMD calibrations](log/bandwidths.csv)

| Axis | Full-axis Δ PCC ↑ | Full-axis Δ RMSE ↓ |
| --- | ---: | ---: |
| HVG-DEG 2k | 0.4209 ± 0.0123 | 0.4712 ± 0.0050 |
| HVG-DEG 3k | 0.4160 ± 0.0088 | 0.3371 ± 0.0047 |
| HVG-DEG 4k | 0.4243 ± 0.0110 | 0.3000 ± 0.0051 |
| HVG-DEG 5k | 0.4247 ± 0.0099 | 0.2764 ± 0.0043 |
| Denoise-DEG 2k | 0.4212 ± 0.0183 | 0.3210 ± 0.0050 |
| Denoise-DEG 3k | 0.4310 ± 0.0138 | 0.2865 ± 0.0039 |
| Denoise-DEG 4k | 0.4396 ± 0.0146 | 0.2635 ± 0.0043 |
| Denoise-DEG 5k | 0.4413 ± 0.0136 | 0.2386 ± 0.0037 |

The [original dual-scope publication](https://github.com/Moke110/carps/blob/263240e896444197bf045739f5bdb54e23f3be53/stage2_CarPS-VAE/09_cpsv_axis/README.md#formal-results) preserves the earlier results, figures and replay route.

## Figures

| Figure | Plotted values |
| --- | --- |
| [Learning curves](figure/learning-curve.png) | Training and validation objectives by epoch for each fit, including Step 08 reuses. |
| [Full-axis metrics](figure/metrics.png) | Eight held-out metrics on the full prediction axes. |
| [Paired full-axis metrics](figure/metrics-paired.png) | Each score minus the Denoise-DEG 5k score for the same cold-drug fold and seed. |

Use the Step 04 condition-level eligibility rules and foldwise distribution metrics. Check the Step 04 full-axis MMD bandwidth once per fold and keep it fixed across all axes and seeds. Record the bandwidth, ordered training rows and gene indices. Compute signed wJAC-2000 per held-out condition at cutoffs 200, 500, 1,000 and 2,000 before averaging over defined conditions. Summarize coordinates with equal weight and sample SD. Pair coordinates with two defined scores. Each metric panel includes zero and shows each defined replicate as a dot. The plotted CSVs and [figure replay receipt](log/figure-render.json) support independent reproduction.

## Bounded CUDA verification

[Ticket #367](https://github.com/Moke110/carps/issues/367) verified all eight axes on fold 0 and seed 17. The [bounded receipt](log/367-bounded.json) records independently replayed checkpoints, initialized parameter and per-epoch random-stream hashes, the metrics and figure replay for the original publication. The [preflight](log/367-preflight.json) binds all 72 formal coordinates and identifies the nine [Step 08 UniMol reference candidates](log/367-step08-reference.json). [Step 07 exact replay](log/367-step07-replay.json) also passed for all nine alpha coordinates.

Run `execution/axis_matrix.py --bounded` with `--output` and `--logs` under separate scratch paths, using the project's CUDA and process memory limits. Set `OMP_NUM_THREADS=2`, `MKL_NUM_THREADS=2`, `OPENBLAS_NUM_THREADS=2` and `NUMEXPR_NUM_THREADS=2`. Override `--step03`, `--representations` and `--source` when using a linked worktree with shared materialized inputs. Run `execution/axis_reference.py --step08` against the sealed Step 08 checkout path and the same frozen inputs. Run `execution/axis_report.py --bounded` with `--matrix`, `--output` and `--figure-data` to derive held-out scores and the plotted CSVs. The report command accepts `--step03` and `--step04` path overrides. Run `execution/axis_figures.py build` or `verify` with `--data`, `--output` and `--scratch` to check CSV-only figure reproduction.

## Reproduce and verify

Use the pinned `arc` environment with the materialized Step 02, Step 03, Step 04, Step 08, and UniMol inputs. Apply the [execution limits](../../docs/agents/execution-environment.md) and check live host RAM, WSL RAM, GPU and disk. Run the coordinator sampler concurrently in a separate shell, starting before the fit and finishing after its last sample.

```bash
mkdir -p "$PWD/tmp"
work=$(mktemp -d "$PWD/tmp/step09-formal.XXXXXX")
mkdir -p "$work/scratch" "$work/matplotlib"
export TMPDIR="$work/scratch" MPLCONFIGDIR="$work/matplotlib"
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2
stage=cpst
formal=$stage/09_cpst_axis/execution/axis_formal.py

# Separate shell, with the same work path.
python "$stage/06_tsg_shape/06.1_round1/execution/sample_resources.py" \
  --output "$work/coordinator.jsonl"

# Fit and validate one coordinate at a time in a capped scope.
systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=8G -- \
  /usr/bin/time -v python "$formal" fit --output "$work/bundle"

# Wait for the coordinator's first post-fit sample before building the candidate.
systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=8G -- \
  python "$formal" build --bundle "$work/bundle" --coordinator "$work/coordinator.jsonl" \
  --publication "$work/publication" --scratch "$work/scratch"
```

Before verifying a materialized publication, check `free -g` and confirm at least 12 GiB is available in WSL. Run the independent rescore alone under the same memory and swap cap and record its peak RSS.

```bash
systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=8G -- \
  /usr/bin/time -v -o "$work/verify.resources" \
  python "$formal" verify --scratch "$work/scratch"
```

The verifier reopens the 72 fit records and checks the full-axis scores, equal-weight summaries and three figures against the published bytes. To verify a candidate, pass its `--bundle` and `--publication` paths. Use `install` with the verified bundle to materialize the local fit products at `output/formal` once the tracked figure and report receipts are available. For linked worktrees, specify the materialized roots with `--step03`, `--step08`, `--representations` and `--source` for `fit` and with `--step03`, `--step04` and `--step08` for scoring.

To rebuild only the report from the accepted fit evidence, point the candidate bundle at the existing fit matrix, fit audit and reuse plan. Pass the original coordinator samples to verify the recorded fit resources.

```bash
sealed="$PWD/$stage/09_cpst_axis/output/formal"
mkdir "$work/report-bundle"
for item in matrix matrix-audit reuse-preflight.json; do
  ln -s "$sealed/$item" "$work/report-bundle/$item"
done
python "$formal" build --bundle "$work/report-bundle" \
  --coordinator "$stage/09_cpst_axis/log/368-coordinator-resources.jsonl" \
  --publication "$work/full-axis-publication" --scratch "$work/scratch"
python "$formal" install --bundle "$work/report-bundle" \
  --publication "$work/full-axis-publication" --destination "$sealed"
```

The installer validates the fit identity and replaces the sealed predecessor report with the verified full-axis report. The fit matrix and reuse lineage are kept in place.
