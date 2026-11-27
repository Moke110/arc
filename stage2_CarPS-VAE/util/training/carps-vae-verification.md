# Steps 07–10 shared-runtime verification

The [receipt inventory](carps-vae-receipts.sha256) binds the installed
Step 07 release and the formal rosters, checkpoint ledgers, report tables,
lineage and figure evidence for Steps 08–10. The matrix manifests bind each
fresh registration, phase manifest and checkpoint through its recorded
path and digest. The reuse plans bind the exact predecessor registrations
and replay evidence. The scientific and figure receipts record the
generated tables, eligible counts, objective traces, CSVs, PNGs and
rendering environment.

| Study | Accepted result | Original registration code | Fresh and reused coordinates |
| --- | --- | --- | --- |
| [07 alpha](../../07_cpsv_alpha/README.md) | `output/release.json`, `log/362-alpha.json`, Step 06.4 formal matrix and report | `825f48c96a0968c11b0bffa9a8b06686b0d16e03` and Step 06.4 `f387f2c370eb87cd27c02502aa98f165851491a7` | 9 Step 06.4 replays |
| [08 chemical representations](../../08_cpsv_chem-emb/README.md) | `output/formal/{matrix,report,reuse-preflight.json}`, `log/364-formal.json`, `figure/` | `ac62a177a54aa99e6ac626728af12394f3e0d940` | 54 fresh, 9 Step 06.4 reuses |
| [09 prediction axes](../../09_cpsv_axis/README.md) | `output/formal/{matrix,report,reuse-preflight.json}`, `log/371-formal.json`, `figure/` | `d3395ca0e0a055dc28e8e431740f5ec640ee4247` | 63 fresh, 9 Step 08 reuses |
| [10 dose inputs](../../10_cpsv_dose/README.md) | `output/formal/{matrix,matrix-audit,report}`, `log/378-formal.json`, `figure/` | `b8be1ad9244c9a3e03eba3e0c77ec64a59ed695d` | 90 fresh across CarPS-VAE and native TranSiGen |

The current registration API exposes architecture, objective coefficients,
frozen chemical representation, dose conditioning and transform, ordered
gene axis, split, fold and seed. The bounded and formal preflights bound all
225 formal and 25 short configurations for Steps 08–10. All 207 fresh
published registrations matched the complete scientific configuration
and the Git blobs at their original revisions. The Step 09 and 10 formal
preflights also matched their pinned bounded gates. The Step 07 release
verifier checked its 30 source hashes and the Step 06.4 pin. The Steps 04–06
[consolidation verification](verification.md) records their original source
revisions and comparison scope.

## History-isolated acceptance

The ticket worktree contained the current tracked code and links to installed active
payloads. The retired Stage 2 archive was absent. The active-input verifier opened
the raw corpus, checked a measured tenv window with its paired same-plate
DMSO controls, checked native-axis projections, verified all eight ordered
tenv training releases and their folds, and read both ends of the active Tahoe
release. The shared figure gate rendered all 45 current PNGs twice from the
published CSVs and matched the tracked receipts, covering 13 comparisons,
five Step 02 figures and the palette.

The installed scientific verifiers passed exactly. Step 08 verified 567
scores and receipt `56b8983678f72adc2deb9a34613561428b3bd77507f44b0ae091b1bc7d40a6f2`.
Step 09 verified 720 scores, 33,192 condition rows and receipt
`b8d20c28479c49348a90ef7d6c2071d5c741451fc8f5bd3ee4568de5f09472eb`.
Step 10 verified 900 scores, 13,830 condition rows and receipt
`c39a2d6a4b08c6b8c22b2d45fa79e7a72fd28f925301e8b23bec68ebe1dbbd1d`.
Each verifier also reproduced its accepted CSV-only figures.

The Step 08 registry and preflight validated the original code with the pinned
Git revision on the first pass. Its score reader still bound the working-tree
runtime hash and the additional `experiment` field when comparing a sealed
registration. That reader now performs the exact same historical code check and
scientific-configuration comparison as the phase runner. Step 09's original
renderer also binds the axis protocol source byte for byte. The shared-score
import is initialized by its caller so the renderer retains that dependency.

Run from the repository root with the materialized active releases and the
`arc` environment. Set `WORK` to a fresh directory under the repository
`tmp/`. Use the [execution limits](../../../docs/agents/execution-environment.md)
and run the commands serially.

```bash
stage=stage2_CarPS-VAE
sha256sum -c "$stage/util/training/carps-vae-receipts.sha256"
python "$stage/04_tsg_axis/execution/verify_active_inputs.py"
python "$stage/07_cpsv_alpha/execution/verify_alpha.py"
python "$stage/08_cpsv_chem-emb/execution/chem_formal.py" verify --scratch "$WORK"
python "$stage/09_cpsv_axis/execution/axis_formal.py" verify --scratch "$WORK"
python "$stage/10_cpsv_dose/execution/dose_formal.py" verify --scratch "$WORK"
python "$stage/util/figure/script/verify_report_figures.py" --work "$WORK/figures"
```

## Bounded CUDA fits and sealed replays

Ran each numbered runner with `--fold 0 --seed 17` and the listed study mode
and arm. The pinned original checkout was
`f0d5bb673435b1a0341301c802ceac5094708a54`. Both implementations ran
`register`, `train`, `predict` and `replay` as separate processes.
Each pair matched all tensors, initialization, objective and learning traces,
validation predictions, held-out predictions and Delta arrays, and condition
and gene ordering. Step 09 and 10 training stream records also matched.
Independent condition and fold scoring produced byte-identical tables.

| Runner | Registration flags | Equal checkpoint tensors | Fold metrics CSV SHA-256 |
| --- | --- | ---: | --- |
| `08_cpsv_chem-emb/execution/chem_runs.py` | `--arm UniMol --smoke` | 51 | `d3b4fd200b7420e28d09a4cb5f6f645b95f545a75af0c66890deb17da3b6eb2a` |
| `09_cpsv_axis/execution/axis_runs.py` | `--axis hvg-2k --bounded` | 51 | `0999940f9168d73ff7e2e42c8071e18567ae37f332374826c911c2f170bdb067` |
| `10_cpsv_dose/execution/dose_runs.py` | `--arm Zmol-FiLM --bounded` | 55 | `9faec134eeddefea23ee3f5f66270921e956d471980ba5e9de5dec7b3b8e89bd` |
| `10_cpsv_dose/execution/dose_runs.py` | `--arm TranSiGen --bounded` | 51 | `ed756a57f26c1a3663ebf8bb450ddf2b42e083c69d300e6491b4620d12921f10` |

The current shared API replayed the corresponding four original formal
checkpoints from isolated scratch copies. Each passed the pinned source,
validation objective and prediction checks and matched the sealed held-out
predictions, signed Delta arrays and row identities. It also reproduced
all nine accepted Step 07 alpha coordinates exactly.

The fit and replay scope used Python 3.13.13, PyTorch 2.12.0+cu130, CUDA 13.0,
NumPy 2.4.6 and two Torch threads on an RTX 4070 Ti SUPER. The cap was
10 GiB RAM and 8 GiB swap. The complete eight-fit and four-replay parity
run took 139.64 seconds with peak process RSS of 1,802,480 KiB.

All 32 focused Python tests, both path fixtures, the tracked-file check,
and the Stage 2 path audit passed. The isolated checks leave materialized
storage retirement to the deployment gate.

## Integration review verification

After moving the KL diagnostics, paired architecture scores and chemical
alignment into shared modules, 114 targeted tests passed. The configuration
gate rebound all 225 formal and 25 bounded coordinates and verified all 207
sealed fresh registrations. New Step 10 identities include the current
early-stopping, RNG and input runtime. The historical Step 10 gate reads its
original four-file identity from Git.

Both receipt checksum ledgers and the original Step 08 UniMol and Step 10
Zmol-FiLM CUDA checkpoint replays passed. A fresh two-epoch Step 08 UniMol
fit matched all 51 tensors, initialization, traces, KL diagnostics, predictions,
Delta arrays and row identities. The replay and fit scope peaked at
1,800,656 KiB process RSS.
