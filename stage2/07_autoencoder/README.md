# 07_autoencoder — whole-transcriptome reconstruction

Step 07 compares methods that reconstruct a measured expression profile from
itself. It uses the corrected [Step 02 shared CM release](../02_dataset_corpus/README.md)
with 18,164 ordered genes and measured FAM153B in both Tahoe and tenv.
[ADR 0052](../../docs/adr/0052-stage2-autoencoder-reconstruction-boundary.md)
defines the benchmark. Step 08 owns perturbation prediction and its integration
checks. The [historical Latent XPert workflow](../_archive/07_lxpert/README.md)
documents the earlier checkpoint, figure and treatment smoke.

## Frozen split

Plates are ranked separately within each source using SHA-256 over the original
JSON serialization of `[20261109, source, plate]`. The first fifth, rounded
down, is validation. The next equally sized block is test, and the remaining
plates are training. The [frozen split](output/split.json) lists every plate
by role, the ranked-plate digest, and each role's ordered profile membership
digest. Membership uses expression row, native row, profile ID and plate,
in release row order. Loading the split verifies its corpus manifest, gene
order, ranks, memberships and full coverage.

| Source | Train plates / profiles | Validation plates / profiles | Test plates / profiles |
| --- | ---: | ---: | ---: |
| Tahoe | 10 / 47,822 | 2 / 9,600 | 2 / 9,596 |
| tenv | 10 / 2,810 | 3 / 843 | 3 / 825 |

## Held-out reconstruction comparison

The [AAE ablation seal](output/test_seal_aae_ablation.json) binds 11 configurations
to their selected validation evidence and fitted payloads. All configurations
reconstructed the same 9,596 Tahoe and 825 tenv test profiles on the 18,164-gene
axis. The [test scores](output/test_scores_aae_ablation.json) give full-precision
metrics, profile counts and provenance. The [comparison record](execution/aae_ablation_record/comparison.json)
records the seal commit, resource measurements and byte-identical independent
test replay. The nine earlier configurations reproduced their published
metrics, counts and provenance exactly.

![Held-out reconstruction metrics for eleven configurations on Tahoe and tenv](output/reconstruction_test_aae_ablation.png)

[ADR 0054](../../docs/adr/0054-step07-separate-reconstruction-and-nuisance-tokens.md)
labels the published cross-attention model `AAE_512_SA` in this comparison. Its artifact
ID is `AAE_512`. The new `AAE_512` (artifact ID `AAE_512_Transformer`) adds
two encoder Transformer blocks to the 512 × 128 representation.
`AAE_512+16` gives those blocks 16 additional encoder tokens for a 528 × 128
representation. Both new decoders read the first 512 tokens, or 65,536 scalars.

| Display name | Tahoe PCC | Tahoe RMSE | Tahoe relative squared error | tenv PCC | tenv RMSE | tenv relative squared error |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `AAE_512_SA` | 0.96793 | 0.09617 | 1.07843 | 0.67100 | 0.30180 | 0.35547 |
| `AAE_512` | 0.94682 | 0.13160 | 1.78404 | 0.55299 | 0.34073 | 0.44497 |
| `AAE_512+16` | 0.94330 | 0.13257 | 1.80668 | 0.53666 | 0.34637 | 0.45901 |

The signed deltas subtract the first model's score from the second model's
score using the full-precision values.

| Contrast | Source | Δ PCC | Δ RMSE | Δ relative squared error |
| --- | --- | ---: | ---: | ---: |
| `AAE_512_SA` → `AAE_512` (two encoder blocks) | Tahoe | −0.02111 | +0.03544 | +0.70561 |
| `AAE_512_SA` → `AAE_512` (two encoder blocks) | tenv | −0.11801 | +0.03893 | +0.08949 |
| `AAE_512` → `AAE_512+16` (16 encoder tokens) | Tahoe | −0.00352 | +0.00096 | +0.02263 |
| `AAE_512` → `AAE_512+16` (16 encoder tokens) | tenv | −0.01633 | +0.00564 | +0.01404 |

Both changes reduced PCC and increased both error metrics on both sources.
The added blocks produced the larger deterioration within each source and
metric. `AAE_512` reached the 10,000-update cap and selected its checkpoint
at 6,332 updates. `AAE_512+16` stopped at 9,498 updates after three validation
checks without improvement and selected its checkpoint at 4,749 updates.
Both selections used tenv validation MSE and reproduced their validation
reports byte-for-byte. The [fit and seal record](execution/aae_ablation_record/fit-and-seal.json)
binds the checkpoints to their validation reports and epoch histories.

These are measured-expression `log1p(CP10K)` reconstruction results on one
frozen plate split with one seeded fit per configuration. Variation across
training seeds and plate splits remains unmeasured. The extra tokens are
exploratory, and assigning them a biological nuisance role would require
separate evidence.

### Historical nine-configuration comparison

Each of the nine configurations reconstructed the same 9,596 Tahoe and 825
tenv test profiles. The [selection seal](output/test_seal.json) was committed
before scoring. It binds the selected validation scores, fitted manifest and
payload hashes to the corrected corpus, ordered gene axis and frozen split.

`MatrixAE_64x8` improved on `AAE_512` across all three metrics on both
sources. The relative squared error fell from 1.0784 to 0.9269 on Tahoe and
from 0.3555 to 0.1653 on tenv. Mean profile PCC rose from 0.9679 to 0.9725 on
Tahoe and from 0.6710 to 0.8571 on tenv. `PCA_512` led all three metrics on
both sources, with relative squared errors of 0.4087 and 0.1208.
The eight prior score reports match the previous comparison at `1654ec0`.

The [test scores](output/test_scores.json) give full-precision metrics,
including the number of defined PCCs and the profile count for every source
and configuration. The dotted line marks the training-gene mean in each panel.
The [original comparison](../_archive/07_lxpert/README.md) remains in the archive
with its checkpoint manifest and perturbation smoke.

### Replay the published scores

The [Step 02 shared CM release instructions](../02_dataset_corpus/README.md#shared-cm-expression-release)
give the build and verification commands for `output/corpus/shared-cm-v1`,
including `genes.csv` in model order. The release's
[verification record](../02_dataset_corpus/execution/build_record/shared_cm_v1.json)
pins its manifest. Step 07 publishes [the frozen split](output/split.json)
and [the fitted artifact manifests](output/models/). The fitted parameters,
including the mean, PCA arrays and neural checkpoints, stay local. To replay
the published scores from a fresh checkout, restore the original payloads to
the paths listed in the manifests and verify their digests against the seal.
A new fit can be evaluated on validation and may produce different bytes.
Fit and validation resource records
are available for [mean](execution/mean_record/run.json),
[PCA](execution/pca_record/run.json), [MLP](execution/mlp_record/run.json),
[AAE](execution/aae_record/run.json), and [MatrixAE](execution/matrixae_record/run.json).

From the repository root in the `arc` environment, run replay with fresh
resource admission as recorded in the [comparison run](execution/comparison_record/run.json).
Replay verifies the seal, loads all nine local artifacts and checks the
published scores with relative tolerance `1e-6` and absolute tolerance `1e-8`.
Mean and PCA score on CPU. MLP, AAE and MatrixAE score on CUDA. Run CUDA jobs serially
and check free GPU memory with `nvidia-smi` before starting. The replay record
uses a 5,120 MiB free GPU memory floor. Use `run_step07` from this block for
the fit and validation examples below in the same shell.

```bash
run_step07() {
  local available
  available=$(free -g | awk '/^Mem:/ {print $7}')
  if [ "${available:-0}" -lt 12 ]; then
    printf 'Step 07 requires at least 12 GiB of available RAM.\n' >&2
    return 1
  fi
  systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=8G -- \
    /usr/bin/time -v env OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 \
    NUMEXPR_NUM_THREADS=2 python "$@"
}

step=stage2/07_autoencoder
cli="$step/execution/script/compare.py"
scores="$step/output/test_scores.json"
run_step07 "$cli" reload --allow-test --expected "$scores"
python "$cli" figure --expected "$scores" --output /tmp/arc-reconstruction-test.png
```

To replay the ablation, restore the fitted payloads for all 11 configurations,
including the checkpoints under `output/models/AAE_512_Transformer/` and
`output/models/AAE_512+16/`. Use the same `run_step07` helper in the same
shell. The checks below verify the committed seal and GPU admission before
test access. Replay checks the published report at the same numerical
tolerances as the historical comparison.

```bash
(
  set -euo pipefail
  step=stage2/07_autoencoder
  cli="$step/execution/script/compare.py"
  seal="$step/output/test_seal_aae_ablation.json"
  scores="$step/output/test_scores_aae_ablation.json"
  seal_commit=657da8d8ba2e4790f692d0bc5c87fcacb282f744
  git merge-base --is-ancestor "$seal_commit" HEAD
  git show "$seal_commit:$seal" | cmp - "$seal"
  git show "HEAD:$seal" | cmp - "$seal"
  gpu_free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits)
  test "$gpu_free" -ge 5120
  run_step07 "$cli" reload --seal "$seal" --allow-test --expected "$scores"
  python "$cli" figure --expected "$scores" --output /tmp/arc-reconstruction-test-aae-ablation.png
)
```

## Fit, score and replay

From the repository root in the `arc` environment, schedule real commands
serially under the [resource limits](../../docs/agents/execution-environment.md).
The CLI checks every release payload against the manifest previously verified
by Step 02. It memory-maps the expression matrices and reads bounded batches.

```bash
cli=stage2/07_autoencoder/execution/script/reconstruct.py
run_step07 "$cli" freeze-split --split /tmp/arc-split.json
run_step07 "$cli" fit --split /tmp/arc-split.json --output /tmp/arc-mean
run_step07 "$cli" score --split /tmp/arc-split.json --artifact /tmp/arc-mean \
  --role validation --output /tmp/arc-mean-validation.json
run_step07 "$cli" reload --split /tmp/arc-split.json --artifact /tmp/arc-mean \
  --role validation --expected /tmp/arc-mean-validation.json
```

The published split and [local mean artifact](output/models/mean/manifest.json) are
the defaults for scoring. The mean reduces all joint training rows in float64,
casts the result to float32 once and broadcasts the same vector to every
profile. Each published artifact carries the corpus, ordered gene axis, split
and file digests. Replay checks this binding and independently recomputes the
score report with relative tolerance `1e-6` and absolute tolerance `1e-8`.

The scorer computes PCC and RMSE across genes per profile, then averages
profiles equally within each source. A profile has undefined PCC if either
vector has zero variance. The report includes the number of defined PCCs and
uses JSON `null` when none are defined. MSE averages all squared residuals
within a source. Relative squared error divides the model's summed squared
residuals by the training-gene mean's summed squared residuals on the same
profiles and is `null` if the denominator is zero.

The library exposes `Corpus.read(source, rows)`, verified row memberships
at `Split.rows[source][role]`, and `score(corpus, split, role, predict, mean,
batch_size)`. Entrants pass a batch prediction function returning an array
in the same row and gene order. `publish_artifact` provides a staging directory
for the fitted payloads and publishes them with a manifest as one directory.
`verify_artifact` checks the binding and digests before reload. Epoch-based
entrants select fitted artifacts using tenv validation MSE. Each corpus read
and configurable batch is limited to 1,024 profiles. For final publication,
`compare.py score` and `compare.py reload` require `--allow-test` and verify
the selection seal before scoring real test profiles. Per-artifact test
scoring is available for bounded fixtures.

Run bounded fixtures with
`python -m unittest discover -s stage2/07_autoencoder/execution/tests`.

The mean fit used 50,632 profiles across the two sources. Independent validation
replay reproduced the [score report](execution/mean_record/validation.json)
byte-for-byte. The [run record](execution/mean_record/run.json) gives the
admission measurements and peak RSS.

| Source | Validation profiles | Mean PCC | Mean RMSE | Relative squared error |
| --- | ---: | ---: | ---: | ---: |
| Tahoe | 9,600 | 0.978971 | 0.074663 | 1.0 |
| tenv | 843 | 0.184086 | 0.512094 | 1.0 |

### PCA reconstruction

One joint training fit estimates the per-gene center and a 512-component basis.
The solver starts with a Gaussian sketch seeded at 20261109, oversamples by 64,
and applies eight streamed covariance products with reduced QR after each pass.
Rayleigh–Ritz orders the resulting directions by variance. The covariance
products use float32 at the highest matmul precision. Each artifact carries
its exact prefix of the basis and the fitted center. Inference projects the
centered profile, maps it back into expression space and adds the center.

```bash
cli=stage2/07_autoencoder/execution/script/reconstruct.py
run_step07 "$cli" fit --model PCA --device cuda --batch-size 512 --output /tmp/arc-pca
for rank in 32 128 512; do
  artifact="stage2/07_autoencoder/output/models/PCA_$rank"
  scores="/tmp/arc-pca-validation-$rank.json"
  run_step07 "$cli" score --artifact "$artifact" --role validation --output "$scores"
  run_step07 "$cli" reload --artifact "$artifact" --role validation --expected "$scores"
done
```

The example scores the original local artifacts. To validate a new fit, use the
rank directories written under `/tmp/arc-pca`. CPU scoring and replay use the
same published training-gene mean reference as the other methods. Use
`--mean-artifact` if the reference was fitted to a fixture or another release.

### MLP reconstruction

Each `MLP_k` fits a separate fully connected autoencoder with widths
18,164 → 512 → *k* → 512 → 18,164. GELU follows the two 512-wide hidden
layers. The bottleneck and output are linear, and *k* counts scalars.
The frozen recipe uses float32, full-axis MSE, AdamW with learning rate
0.001 and weight decay 0.01, and 20 epochs of 512-profile batches. A
PCG64 shuffle seeded by 20261109 and the epoch visits every joint training
profile once at the natural source ratio. Each fit has a 30-minute ceiling.
Validation uses 64-profile batches, records MSE separately by source, and
selects the earliest epoch at the minimum tenv MSE.

```bash
cli=stage2/07_autoencoder/execution/script/reconstruct.py
for width in 32 128 512; do
  artifact="/tmp/arc-MLP_$width"
  scores="/tmp/arc-MLP_$width-validation.json"
  run_step07 "$cli" fit --model "MLP_$width" --device cuda --output "$artifact"
  run_step07 "$cli" score --artifact "$artifact" --device cuda --role validation --output "$scores"
  run_step07 "$cli" reload --artifact "$artifact" --device cuda --role validation --expected "$scores"
done
```

### Attention Autoencoder reconstruction

`AAE_512_SA` (artifact ID `AAE_512`) uses 512 latent queries, 128 channels per
query and eight attention heads. It contains 65,536 latent scalars per profile.
Fixed sinusoidal positions encode the order of the 18,164 genes. The encoder
queries all measured expression values and the decoder uses the same ordered gene positions to
reconstruct the full profile.

The frozen recipe uses full-axis MSE, AdamW with learning rate 0.001,
weight decay 0.01 and batches of 32 profiles. A PCG64 shuffle seeded by
20261109 and the epoch visits every joint training profile once per full
epoch. Validation runs after each epoch and at the 10,000-update cap,
including when the cap falls within an epoch. Training stops after three
consecutive validation checks without a lower tenv MSE and has a 30-minute
wall time ceiling. The earliest checkpoint with the lowest tenv validation
MSE is selected. Tahoe validation MSE is recorded separately.
Real fitting and scoring require native CUDA BF16 and use
FlashAttention exclusively. Validation, scoring and replay use the same
32-profile batch size. The bounded CPU fixture uses float32, a five-update
cap and two-profile batches.

![AAE joint training MSE and source-specific validation MSE by optimizer updates](output/aae_loss.png)

The [epoch history](output/models/AAE_512/epochs.json) records the joint
training MSE from the optimizer's batch losses and the two validation MSEs
at each check. The run reached 10,000 updates after six complete epochs and
502 batches of the seventh. The selected checkpoint was at the cap, with
validation MSE of 0.009339 on Tahoe and 0.093597 on tenv. The updated comparison
evaluates that checkpoint on the same previously reported test plates.

```bash
cli=stage2/07_autoencoder/execution/script/reconstruct.py
artifact=/tmp/arc-AAE_512
scores=/tmp/arc-AAE_512-validation.json
run_step07 "$cli" fit --model AAE_512 --device cuda --output "$artifact"
run_step07 "$cli" score --artifact "$artifact" --device cuda --role validation --output "$scores"
run_step07 "$cli" reload --artifact "$artifact" --device cuda --role validation --expected "$scores"
python stage2/07_autoencoder/execution/script/compare.py aae-loss-figure \
  --history "$artifact/epochs.json" --output /tmp/arc-aae-loss.png
```

### Compact matrix autoencoder reconstruction

`MatrixAE_64x8` uses 64 learned latent queries and eight bottleneck channels
per query. Its 512 latent scalars match the budget of `PCA_512` and
`MLP_512`. Learned gene identities and projected measured expression feed
the encoder, and gene identity queries read the latent tokens to reconstruct
the full axis. Two latent Transformer blocks on either side of the
bottleneck use 128 internal channels. The model adapts the deterministic
autoencoder path of [scLDM](https://arxiv.org/html/2511.02986v2) to the
full measured Step 07 input and MSE objective. [ADR 0053](../../docs/adr/0053-step07-compact-matrix-autoencoder.md)
records the design.

The fit uses the AAE optimizer, batch size, early stopping and update cap.
It reached 10,000 updates and selected the checkpoint at 6,332 updates,
with validation MSE of 0.007422 on Tahoe and 0.044412 on tenv.
The [epoch history](output/models/MatrixAE_64x8/epochs.json) records the
joint training batch losses and both validation losses.

![MatrixAE joint training MSE and source-specific validation MSE by optimizer updates](output/matrix_ae_loss.png)

```bash
cli=stage2/07_autoencoder/execution/script/reconstruct.py
artifact=/tmp/arc-MatrixAE_64x8
scores=/tmp/arc-MatrixAE_64x8-validation.json
run_step07 "$cli" fit --model MatrixAE_64x8 --device cuda --output "$artifact"
run_step07 "$cli" score --artifact "$artifact" --device cuda --role validation --output "$scores"
run_step07 "$cli" reload --artifact "$artifact" --device cuda --role validation --expected "$scores"
python stage2/07_autoencoder/execution/script/compare.py matrix-ae-loss-figure \
  --history "$artifact/epochs.json" --output /tmp/arc-matrix-ae-loss.png
```
