# 07_autoencoder — whole-transcriptome reconstruction

Step 07 compares methods that reconstruct a measured expression profile from
itself. It uses the corrected [Step 02 shared CM release](../02_dataset_corpus/README.md)
with 18,164 ordered genes and measured FAM153B in both Tahoe and tenv.
[ADR 0052](../../docs/adr/0052-stage2-autoencoder-reconstruction-boundary.md)
defines the benchmark. Step 08 owns perturbation prediction and its integration
checks. The [historical Latent XPert workflow](../_archive/07_lxpert/README.md)
contains the earlier checkpoint, figure and treatment smoke.

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

## Fit, score and replay

From the repository root in the `arc` environment, schedule real commands
serially under the [resource limits](../../docs/agents/execution-environment.md).
The CLI checks every release payload against the manifest previously verified
by Step 02. It memory-maps the expression matrices and reads bounded batches.

```bash
cli=stage2/07_autoencoder/execution/script/reconstruct.py
python "$cli" freeze-split --split /tmp/arc-split.json
python "$cli" fit --split /tmp/arc-split.json --output /tmp/arc-mean
python "$cli" score --artifact /tmp/arc-mean --role validation --output /tmp/arc-mean-validation.json
python "$cli" reload --artifact /tmp/arc-mean --role validation --expected /tmp/arc-mean-validation.json
```

The published split and [mean artifact](output/models/mean/manifest.json) are
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
entrants select fitted artifacts using tenv validation MSE. Real test scoring
requires `--allow-test` and is reserved for final publication.

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
python "$cli" fit --model PCA --device cuda --batch-size 512 --output /tmp/arc-pca
for rank in 32 128 512; do
  artifact="stage2/07_autoencoder/output/models/PCA_$rank"
  scores="/tmp/arc-pca-validation-$rank.json"
  python "$cli" score --artifact "$artifact" --role validation --output "$scores"
  python "$cli" reload --artifact "$artifact" --role validation --expected "$scores"
done
```

The example scores the published artifacts. To validate a new fit, use the
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
  python "$cli" fit --model "MLP_$width" --device cuda --output "$artifact"
  python "$cli" score --artifact "$artifact" --device cuda --role validation --output "$scores"
  python "$cli" reload --artifact "$artifact" --device cuda --role validation --expected "$scores"
done
```
