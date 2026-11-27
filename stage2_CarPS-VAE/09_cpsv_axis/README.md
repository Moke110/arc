# Step 09 CarPS-VAE prediction axes

**Status** Formal comparison complete. See [ADR 0101](../../docs/adr/0101-carps-vae-prediction-axes.md) and the [formal acceptance receipt](log/368-formal.json).

Compare CarPS-VAE across eight [Step 02 prediction axes](../02_gene_axis/README.md), using the [Step 08 UniMol representation](../../docs/adr/0100-carps-vae-unimol.md) and pairing each arm with `denoise-5k` by cold-drug fold and training seed. Use folds 0–2 and seeds 17, 42 and 2026, for 72 axis × fold × seed coordinates.

| Axis family | Axis IDs | Disease-union overlap at 2k / 3k / 4k / 5k |
| --- | --- | --- |
| HVG-DEG | `hvg-2k`, `hvg-3k`, `hvg-4k`, `hvg-5k` | 1,258 / 1,889 / 2,525 / 3,149 |
| Denoise-DEG | `denoise-2k`, `denoise-3k`, `denoise-4k`, `denoise-5k` | 906 / 1,225 / 1,375 / 1,513 |

[Step 03](../03_data_train/README.md) owns the two additional 2k training views and the six established views. Pool matched-plate counts and normalize on each selected axis with the same method. Verify common condition keys, raw membership, exclusions and frozen partition order. Shared genes can have different normalized values across axes.

Resize the input projections of both expression encoders and the output projections of the control and shared Δ decoders. Keep the internal widths and architecture, UniMol input and projection, matched-control protocol, fold-local dose transform, loss and fit policy fixed. Use common initialized values for same-shaped parameters and a paired training random stream. Anchor pairing to the Step 08 UniMol reference so equivalent reference fits can be reused. Verify each reference coordinate's provenance, protocol equivalence and exact checkpoint replay. With all nine reference coordinates accepted, 63 coordinates require new training. Freshly fit any rejected reference coordinate.

## Formal results

The matrix covers 72 coordinates. All nine Step 08 UniMol Denoise-DEG 5k fits passed the [reuse audit](log/368-reference.json). All 63 new fits and nine reused checkpoints replayed their validation and held-out predictions exactly. The table shows the equal-weight mean ± sample SD over nine fold–seed coordinates per axis. [All metrics and paired changes](log/summary.csv) · [fit lineage](log/runs.csv) · [MMD calibrations](log/bandwidths.csv)

| Axis | Full Δ PCC ↑ | Full Δ RMSE ↓ | Disease Δ PCC ↑ | Disease Δ RMSE ↓ | Disease direction-correct coverage ↑ |
| --- | ---: | ---: | ---: | ---: | ---: |
| HVG-DEG 2k | 0.4209 ± 0.0123 | 0.4712 ± 0.0050 | 0.4355 ± 0.0132 | 0.4247 ± 0.0068 | 23.32 ± 0.18% |
| HVG-DEG 3k | 0.4160 ± 0.0088 | 0.3371 ± 0.0047 | 0.4322 ± 0.0078 | 0.2895 ± 0.0059 | 34.72 ± 0.30% |
| HVG-DEG 4k | 0.4243 ± 0.0110 | 0.3000 ± 0.0051 | 0.4399 ± 0.0091 | 0.2566 ± 0.0057 | 46.56 ± 0.42% |
| HVG-DEG 5k | 0.4247 ± 0.0099 | 0.2764 ± 0.0043 | 0.4407 ± 0.0083 | 0.2362 ± 0.0049 | 57.98 ± 0.47% |
| Denoise-DEG 2k | 0.4212 ± 0.0183 | 0.3210 ± 0.0050 | 0.4249 ± 0.0172 | 0.3067 ± 0.0062 | 16.91 ± 0.21% |
| Denoise-DEG 3k | 0.4310 ± 0.0138 | 0.2865 ± 0.0039 | 0.4359 ± 0.0126 | 0.2710 ± 0.0051 | 22.92 ± 0.20% |
| Denoise-DEG 4k | 0.4396 ± 0.0146 | 0.2635 ± 0.0043 | 0.4464 ± 0.0124 | 0.2486 ± 0.0053 | 25.84 ± 0.27% |
| Denoise-DEG 5k | 0.4413 ± 0.0136 | 0.2386 ± 0.0037 | 0.4480 ± 0.0135 | 0.2241 ± 0.0047 | 28.40 ± 0.33% |

Denoise-DEG 5k had the highest mean Δ PCC and lowest mean Δ RMSE in both scopes. Denoise-DEG 4k was closest on Δ PCC, with paired changes of −0.0017 ± 0.0083 in each scope. Its Δ RMSE rose by 0.0249 ± 0.0012 on the full axis and 0.0245 ± 0.0012 in the disease scope. Denoise-DEG 2k had the highest JAC_top200 in both scopes. HVG-DEG 5k had the highest disease-restricted Δ R² (0.0923), a paired gain of 0.0282 ± 0.0199.

HVG-DEG 5k had the highest direction-correct disease coverage, gaining 29.57 ± 0.33 percentage points over Denoise-DEG 5k. Their disease-union overlaps were 3,149 and 1,513 genes, respectively. The coverage gain reflects both the broader gene overlap and the number of correct nonzero directions. Axis-local normalization gives shared genes different numerical targets across axes. The full score and disease-restricted score comparisons therefore describe the targets available on each fitted axis, and their variability is the sample SD over the registered folds and seeds.

## Figures

Use the nine Stage 2 metrics in their [Step 04 report order](../04_tsg_axis/README.md#score-and-replay-an-axis-comparison): JAC_top200, JAC_pos100, JAC_neg100, Δ PCC, Δ Spearman, Δ RMSE, Δ E-distance, Δ MMD and Δ R².

Let `D` be the 3,481 distinct frozen gene indices in the eligible DEG1000 union for nHF, ICM, genetic cardiomyopathy and TOX, and let `A` be the prediction axis for a run. Define `Δpred_d` and `Δ_d` by slicing predicted and measured Δ to `A ∩ D` in axis order. The disease references define membership. The two response vectors provide the directions and magnitudes that are compared.

| Figure | Plotted values |
| --- | --- |
| [Learning curves](figure/learning-curves.png) | Training and validation objectives by epoch, with one panel per axis and both traces for each fold and seed. |
| [Full-axis metrics](figure/axis-metrics.png) | Nine held-out metrics on each full prediction axis. |
| [Paired full-axis metrics](figure/axis-metrics-paired.png) | Each full-axis score minus the matching `denoise-5k` score. |
| [Disease-restricted metrics](figure/disease-metrics.png) | Nine held-out metrics between `Δpred_d` and `Δ_d`. |
| [Paired disease-restricted metrics](figure/disease-metrics-paired.png) | Each disease-restricted score minus the matching `denoise-5k` score. |
| [Disease direction-correct coverage](figure/disease-correct-coverage.png) | Per held-out condition, `100 × G / 3481`, where `G` counts genes with both values positive or both values negative in the disease-restricted vectors. A zero in either vector contributes zero to `G`. |
| [Paired direction-correct coverage](figure/disease-correct-coverage-paired.png) | Coverage minus the matching `denoise-5k` coverage, in percentage points. |

Use the Step 04 condition-level metrics, eligibility rules and foldwise distribution metrics for both scopes. Use the same full-axis MMD bandwidth as Step 04 for each fold. For the disease-restricted scope, use the same bandwidth estimator on the HVG-DEG 5k measured training Δ sliced to its disease-union intersection. Freeze one bandwidth per fold and scope across axes and seeds. Record the bandwidth, training rows and gene indices used.

Average defined condition-level scores and coverage within each coordinate. Compute distribution scores from all held-out rows of that coordinate. Summarize defined coordinates with equal weight and sample SD. Subtract the reference only for matching fold and seed coordinates with two defined scores. Retain eligible counts and undefined values. Display axis overlap counts with the disease results. The `denoise-5k` direction-correct coverage ceiling is 1,513 / 3,481, or 43.46%.

Publish all seven figures, their replicate values and plotted CSVs, and CSV-only renderers under the [Stage 2 figure contract](../AGENTS.md#step-contract). Use one dot per defined coordinate, means and sample SDs. Include zero on the paired plots and the Δ R² panels.

## Bounded CUDA verification

[Ticket #367](https://github.com/Moke110/carps/issues/367) verified all eight axes on fold 0 and seed 17. The [bounded receipt](log/367-bounded.json) records independently replayed checkpoints, initialized parameter and per-epoch random-stream hashes, both held-out scoring scopes and two exact CSV-only renders of all seven figures. The [preflight](log/367-preflight.json) binds all 72 formal coordinates and identifies the nine [Step 08 UniMol reference candidates](log/367-step08-reference.json). [Step 07 exact replay](log/367-step07-replay.json) also passed for all nine alpha coordinates.

Run `execution/axis_matrix.py --bounded` with `--output` and `--logs` under separate scratch paths, using the project's CUDA and process memory limits. Set `OMP_NUM_THREADS=2`, `MKL_NUM_THREADS=2`, `OPENBLAS_NUM_THREADS=2` and `NUMEXPR_NUM_THREADS=2`. Override `--step03`, `--representations` and `--source` when using a linked worktree with shared materialized inputs. Run `execution/axis_reference.py --step08` against the sealed Step 08 checkout path and the same frozen inputs. Run `execution/axis_report.py --bounded` with `--matrix`, `--output` and `--figure-data` to derive held-out scores and the plotted CSVs. The report command accepts `--step03`, `--step04` and `--disease-release` path overrides. Run `execution/axis_figures.py build` or `verify` with `--data`, `--output` and `--scratch` to check CSV-only figure reproduction.

## Reproduce and verify

Use the pinned `arc` environment with the materialized Step 02, Step 03, Step 04, Step 08, and UniMol inputs. Apply the [execution limits](../../docs/agents/execution-environment.md) and check live host RAM, WSL RAM, GPU and disk. Run the coordinator sampler concurrently in a separate shell, starting before the fit and finishing after its last sample.

```bash
mkdir -p "$PWD/tmp"
work=$(mktemp -d "$PWD/tmp/step09-formal.XXXXXX")
mkdir -p "$work/scratch" "$work/matplotlib"
export TMPDIR="$work/scratch" MPLCONFIGDIR="$work/matplotlib"
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2
stage=stage2_CarPS-VAE
formal=$stage/09_cpsv_axis/execution/axis_formal.py

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

To validate a materialized publication, run `python "$formal" verify --scratch "$work/scratch"`. The verifier reopens the 72 fit records and checks both scoring scopes, the signed coverage, the equal-weight summaries and all figures against the published bytes. To verify a candidate, pass its `--bundle` and `--publication` paths. Use `install` with the verified bundle to materialize the local fit products at `output/formal` once the tracked figure and report receipts are available. For linked worktrees, specify the materialized roots with `--step03`, `--step08`, `--representations` and `--source` for `fit` and with `--step03`, `--step04`, `--step08` and `--disease-release` for scoring.
