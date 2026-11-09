# 07_lxpert — whole-transcriptome Latent XPert

Step 07 trains a transcriptome autoencoder and uses the refit checkpoint in a
Latent XPert treatment smoke. The development split and full-data refit serve
reconstruction. The treatment branch exercise checks integration.
[ADR 0051](../../docs/adr/0051-latent-xpert.md)

## Expression input release

From the repository root, use the `arc` Python environment. Run each real-data
command serially under the [repository resource limits](../../docs/agents/execution-environment.md).

```bash
python stage2/07_lxpert/execution/script/lxpert.py prepare
python stage2/07_lxpert/execution/script/lxpert.py verify
```

The default source is the Step 02 v2 corpus. `prepare --corpus-dir PATH` selects
another copy of the same corpus and checks it against the frozen source digests.
Both commands take `--release PATH`, which defaults to `output/input/` within
this step. Prepare requires a new destination and publishes the directory after
the build succeeds. Verify checks the release files, expression values, masks,
memberships and declared source hashes. Add `verify --corpus-dir PATH` to rehash
the original files as well.

The release includes 67,018 Tahoe and 4,478 tenv profiles. The 32 tenv technical
blanks are excluded. Raw counts are projected onto the ordered 18,164-gene CM
axis, then normalized by the projected library size with CP10K and `log1p`.
The original release maps Tahoe by exact Ensembl ID and tenv by HGNC ID.
It stores Tahoe FAM153B as zero and masks it from the reconstruction loss.
The Step 02 crosswalk identifies a measured Tahoe column under a different
Ensembl ID. [ADR 0052](../../docs/adr/0052-stage2-autoencoder-reconstruction-boundary.md)
records the corrected mapping for the next release. Both sources in the
original release use the complete ordered decoder axis.

For each source, `prepare` ranks plate IDs by a hash of the seed, source and
plate. The first fifth of the plate list, rounded down with a minimum of one
plate, is held out for development validation. The refit uses every released
profile. The default seed is 20261109.

| Release file | Contract |
| --- | --- |
| `genes.csv` | Frozen ordered CM gene axis |
| `{tahoe,tenv}_expression.npy` | Row-major float32 matrices readable via `np.load(..., mmap_mode="r")` |
| `{tahoe,tenv}_target_mask.npy` | Boolean arrays in gene order. True identifies observed reconstruction targets |
| `{tahoe,tenv}_profiles.csv` | One row per expression row, with the original zero-based HDF5 row, profile ID, plate, control flag and development role |
| `manifest.json` | Source, axis, file and membership SHA-256 digests, plus normalization and population counts |

The [verified input release record](execution/input_record/release.json) captures
the published manifest and resource measurements for the real-data build and
verification.

Run bounded CLI tests with
`python -m unittest discover -s stage2/07_lxpert/execution/tests -p test_inputs.py`.

## Bounded autoencoder fit and checkpoint replay

```bash
python stage2/07_lxpert/execution/script/lxpert.py fit \
  --steps 8 --batch-size 32 --output stage2/07_lxpert/execution/output/short_fit
python stage2/07_lxpert/execution/script/lxpert.py reload \
  --checkpoint stage2/07_lxpert/execution/output/short_fit/checkpoint.pt
```

Both commands verify the input release and use CUDA BF16 fused FlashAttention
for the full gene axis. The short fit draws distinct development training
profiles from both sources for each batch and uses AdamW at a learning rate of
0.001. It averages observed-target MSE per profile, with the Tahoe-specific mask.
This balanced sample exercises both data paths. The development training recipe
in ADR 0051 uses the natural source ratio.

The checkpoint binds its weights to the input manifest, ordered gene axis,
normalization and architecture. It stores a reconstruction probe for one
validation profile per source. A separate `reload` process recomputes the probe
using the fit device and precision and checks every output within `atol=1e-6`
and `rtol=1e-5`. The fit report records source counts, losses, step timings,
peak RSS and CUDA allocator peaks. `reload` prints its probe and resource report
as JSON.

Use `--fixture --device cpu` for the bounded fixture. Run its CLI test with
`python -m unittest discover -s stage2/07_lxpert/execution/tests -p test_autoencoder.py`.

The [verified short fit record](execution/short_fit_record/run.json) covers
128 Tahoe and 128 tenv profiles. The independent reload reproduced both
18,164-gene probe predictions exactly. The scheduled fit took 13.51 seconds
wall time and peaked at 1,971,380 KiB host RSS and 1,707,081,728 bytes of CUDA
allocator reservation. After the first step, batch-of-32 steps took 43.6–46.3 ms.

## Formal autoencoder training and release

Run the four commands as separate processes under the scheduler, starting from
the repository root with the `arc` environment. Both fits use CUDA BF16
FlashAttention and the admitted batch size of 32.

```bash
cli=stage2/07_lxpert/execution/script/lxpert.py
run_dir=stage2/07_lxpert/execution/output/formal_fit_run
autoencoder=stage2/07_lxpert/execution/output/autoencoder
python "$cli" freeze --run-dir "$run_dir" --batch-size 32 --seed 20261109 \
  --max-epochs 4 --max-seconds 1800
python "$cli" develop --run-dir "$run_dir"
python "$cli" refit --run-dir "$run_dir" --output "$autoencoder"
python "$cli" reload --checkpoint "$autoencoder/checkpoint.pt"
```

`freeze` seals the AdamW settings, seed, training limits, plate memberships and
input and code digests. Each development epoch shuffles the combined training
population and visits each profile once at the natural source ratio. Its report
gives equal-profile observed-target reconstruction MSE separately for the
held-out Tahoe and tenv plates. The minimum tenv validation MSE selects the
duration, with the earliest epoch breaking ties.

`refit` initializes a new model and optimizer from the frozen seed and uses all
71,496 released profiles for the selected number of epochs. `reload` checks
the checkpoint and provenance digests and reproduces the checkpoint's two
full-axis reconstruction probes. Both training commands emit one progress
record per completed epoch and report elapsed time, peak host RSS, and peak
CUDA allocated and reserved memory.

The Git-tracked release in `output/autoencoder/` contains the 926 KiB
`checkpoint.pt` and the frozen plan, development scores and refit record with
their SHA-256 files. To replay it, prepare and verify the pinned expression
release as described above and use the `arc` environment on CUDA with native
BF16 support.

```bash
python stage2/07_lxpert/execution/script/lxpert.py reload \
  --checkpoint stage2/07_lxpert/output/autoencoder/checkpoint.pt
```

To reuse the published frozen plan for `develop` or `refit`, use training
revision `9032f291dac7b4800551f54b2a76bb9faa9f2bfc` and the input and runtime
recorded in `plan.json`. Copy `plan.json` and `plan.sha256` to a fresh local
run directory, then run `develop` and `refit` there. The plan verifies the exact
CLI and training source digests. The current CLI can replay the published
checkpoint. A fresh `freeze` records a new timestamp, so its refit can have a
different checkpoint digest.

Development trained on 57,418 Tahoe and 3,635 tenv profiles per epoch. The
held-out plates contributed 9,600 Tahoe and 843 tenv profiles to the validation
scores.

| Development epoch | Tahoe validation MSE | tenv validation MSE |
| --- | ---: | ---: |
| 1 | 0.021118 | 0.176949 |
| 2 | 0.010911 | 0.138069 |
| 3 | 0.015075 | 0.174446 |
| 4 (selected) | 0.009845 | 0.126428 |

The [verified formal run record](execution/formal_fit_record/run.json) includes
the frozen plan, per-epoch scores, sample-order verification, release digests
and scheduler measurements. The refit ran for four epochs and independent
reload reproduced both full-axis probe predictions exactly.

| Fit | Process wall time | Peak host RSS | Peak CUDA reserved |
| --- | ---: | ---: | ---: |
| Development | 399.67 s | 6,573,888 KiB | 1,709,178,880 bytes |
| Full-data refit | 444.56 s | 6,623,460 KiB | 1,715,470,336 bytes |

The plate-held-out MSE describes the development fit. The full-data refit
includes treated tenv expression from the earlier cold-drug test folds. A
future comparable cold-drug experiment must first refit the autoencoder with
its held-out treated profiles excluded.

Use `--fixture --device cpu` for the bounded training workflow. Run its CLI tests
with `python -m unittest discover -s stage2/07_lxpert/execution/tests -p test_formal_fit.py`.

## Latent XPert treatment smoke

Run from the repository root with the `arc` environment under the scheduler.

```bash
python stage2/07_lxpert/execution/script/lxpert.py smoke
```

The smoke defaults to the Git-tracked full-data refit and verifies its
published checkpoint digest. Pass `--checkpoint PATH` to smoke-test a new
full-data refit alongside its own `plan`, `development` and `refit` records
and digest files. The smoke verifies the provenance of either checkpoint.
It requires the pinned expression release and the verified Step 03 input
release. The latter supplies a measured tenv well,
its same-plate DMSO mean, its shared Morgan-2048 vector, and its measured
positive dose. The default case is `temp_env::Plate1-C02`. Use `--case-id` to
select another well. The JSON report identifies the input artifacts and rows.

The frozen encoder converts the measured baseline to 512 latent tokens.
The treatment branch applies cross-attention, two self-attention blocks, and
a final cross-attention block. The two cross-attention blocks read the same
two condition tokens. One is a linear projection of the shared fingerprint,
and the other is an affine projection of `log10(dose / 0.05 µM)`. The frozen
decoder maps the treated latent to 18,164 expression values. The smoke checks
that subtracting the measured baseline yields the reported full-axis delta.

A single backward pass against a synthetic target exercises full-axis delta
MSE and Pearson loss. The report checks finite treatment gradients, including
both adapters, and unchanged autoencoder parameters with no gradients.
The command performs no optimizer updates and loads no evaluator truth.
The synthetic probe loss is an integration diagnostic.

Run the bounded CPU fixture with
`python -m unittest discover -s stage2/07_lxpert/execution/tests -p test_treatment_smoke.py`.

The [verified smoke record](execution/treatment_smoke_record/run.json) uses
Plate1's 12 DMSO wells and a 10 µM dose. It produced two 128-wide condition
tokens and finite post-treatment and delta vectors of length 18,164. Delta
matched the independent subtraction exactly. The Morgan and dose adapter
gradient norms were 0.167 and 0.0755. The scheduled run took 11.94 seconds,
peaked at 1,606,172 KiB host RSS and reserved 94,371,840 bytes in the CUDA
allocator.
