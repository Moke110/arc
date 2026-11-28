# Step 06 TranSiGen module studies

Read the [CarPS-Transcriptome workflow](../AGENTS.md) and the
[architecture report](README.md). For runner or scoring changes, use the
[shared TranSiGen training and evaluation APIs](../util/README.md#transigen-training-and-evaluation)
through the numbered entry points.

## Register, fit and replay

Keep the frozen tenv cold-drug folds, ordered Denoise-DEG 5k genes,
measured same-plate DMSO, frozen KPGT and training-fold standardized log
dose fixed across the architecture studies. Use each round's registered
objective, seed roster and reference. For new fits, use the pinned
TranSiGen source and the HVG-DEG 5k release for fold bandwidth calibration.
[Training manifest](../03_data_train/03.7_tenv_denoise-5k/output/manifest.json) ·
[Calibration manifest](../03_data_train/03.1_tenv_hvg-5k/output/manifest.json) ·
[KPGT reader and pins](../util/README.md#shared-input-and-release-apis) ·
[Author source manifest](../../pert_model/TranSiGen/manifest.sha256)

Follow the round guide to register, qualify a bounded run, fit, score and
verify. Use its pinned revision and paths when replaying sealed evidence.
Keep the accepted releases read-only and use scratch copies for fit replay.
[Installed release inventory](../util/training/verification.md)

| Study guide | Formal verification and replay |
| --- | --- |
| [06.1 Expression widths and depths](06.1_round1/AGENTS.md) | [Sealed #347 source](06.1_round1/AGENTS.md#verify-the-accepted-347-source) and [encoder-width extension](06.1_round1/AGENTS.md#reproduce-the-formal-matrix) |
| [06.2 Prediction paths](06.2_round2/AGENTS.md) | [Formal matrix](06.2_round2/AGENTS.md#reproduce-the-formal-matrix) |
| [06.3 Architecture combinations](06.3_round3/AGENTS.md) | [Revised formal run with normalized FiLM](06.3_round3/AGENTS.md#reproduce-the-revised-formal-run) |
| [06.4 Measured encoder targets](06.4_round4/AGENTS.md) | [Formal study](06.4_round4/AGENTS.md#reproduce-the-formal-study) and [sealed fit replay](06.4_round4/AGENTS.md#replay-a-sealed-step-064-fit) |

## Locate the installed releases

The current checkout locates each round's installed bundle at
`cpst/06_tsg_module/<round>/output/formal`. When migrating a materialized
checkout, stage the four `output/formal` directories outside the active tree
before the Git cutover, then place them under their corresponding round names.
Keep their contents byte-for-byte and verify each against its formal receipt.
The [verification inventory](../util/training/formal-receipts.sha256) records
the current receipt paths. The original replay guides use the paths inside
their pinned Git revisions.

## Rebuild the comparison figures

Run the [comparison builder](execution/build_comparison.py) from the
materialized repository root in `arc`. It uses the Denoise-DEG 5k training
release and the accepted Step 06 matrices and reports, including the
sealed #347 fits and the retained completed fits from the first Round 3
attempt. To read those inputs from a separate canonical checkout, pass
`--evidence-root` with that checkout's root. Pass `--step 06.4_round4`,
for example, to rebuild a single study. Keep its predecessor matrices
available for reused fits.

The builder verifies the formal receipts, fit identities and learning
curves. It checks the signed Δ arrays against the held-out predictions
and matched controls and requires parity with the sealed JAC-top200 scores. It subtracts the study
reference within each fold and seed. Require `PASS` for every requested
round and a candidate `log/figure-render.json` recording two matching
renders from its plotted CSVs. The candidate receipt records current
source digests. Published evidence digests pin the original figure build.
Follow the [execution environment](../../docs/agents/execution-environment.md) and
[comparison figure API](../util/README.md#comparison-figures).

```bash
(
  set -e
  LC_ALL=C free -g | awk '$1 == "Mem:" { available = $7 } END { exit available < 12 }'
  mkdir -p "$PWD/tmp"
  work=$(mktemp -d "$PWD/tmp/tsg-architecture-comparisons.XXXXXX")
  mkdir "$work/scratch" "$work/matplotlib"
  export TMPDIR="$work/scratch" MPLCONFIGDIR="$work/matplotlib"
  systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=8G -- \
    env OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 \
    python cpst/06_tsg_module/execution/build_comparison.py --output "$work/publication"
)
```

For CSV-only replay, use the procedure in the numbered guide. To verify
all four published comparisons together, use the
[current CarPS-Transcriptome figure verifier](../README.md#verify-the-current-figures).
