# Formal matrix terminal block — 2026-11-16

The #226/#231 acceptance contract is not met. Frozen formal implementation
`551a3c91c217c83440bae8e5b3b3ca9502a5de10` (integrated as dev `958e668`) completed
18/90 registered runs with exact fresh-process native sampled replay: 15 cold-drug
PRnet runs and cold-drug TranSiGen fold0 seeds17/42/2026. The 19th identity,
`cold-drug / TranSiGen / fold1 / seed17`, failed with
`non-finite native TranSiGen loss`. No incomplete result was published. Old02
remains accepted and intact; tickets/spec remain open and main was not deployed.
Tracker receipts: [#231](https://github.com/Moke110/arc/issues/231#issuecomment-6920746120),
[#226](https://github.com/Moke110/arc/issues/226#issuecomment-6920747417).

## Preserved evidence and numerical boundary

Canonical original candidate (do not modify or relabel):
`/tmp/arc-tenv5k-231/stage2/02_models_tenv-5k/execution/checkpoints/matrix-231`.
`plan.json` SHA256 is
`9b747f83868e305361897dadd76e50afa9f339a87653a3f5a7abd8d98d545858`;
failed `runs/cold-drug/transigen/1/17/terminal.json` SHA256 is
`e91bdf9c8076d3e5dbdd1b3042748edd79c21db0dfb3ce19215f2e9e21b180e2`.
The failed candidate retains config, terminal and its finite best checkpoint
(best epoch130, validation MSE0.06892538476141653). It cannot be declared complete:
epoch155 is only25 stale epochs after the best, below registered patience50.

Independent diagnostic evidence is under `/tmp/arc-tenv5k-execution/`.
`transigen-nonfinite-diagnosis.md` SHA256
`65342c0e9be95e21975ecc5d3571c69b7ac627f80619cf325c33ed271af04bd8`
records four same-identity probes, all failing at epoch155 with identical preceding
losses: baseline, capture-only instrumentation, cuDNN-disabled and full author
backend flags. The inputs and all parameters remain finite. Native predicted
log-variance reaches222.34149; float32 `exp(logvar/2)` becomes Inf before decoding,
so the forward prediction is already NaN. Loss-only float64 arithmetic makes the
KL term finite but cannot repair that non-finite forward prediction. A native
1×1 scalar sampling probe independently reproduces the overflow. No OOM was
observed; available RAM stayed40–41GiB.

Diagnostic SHA256 receipts (filenames relative to that evidence directory):

| Evidence | SHA256 |
|---|---|
| `repro-transigen-1.log` | `0af896aef6a0ec4c343296e81e3b02d9707680e6f56a280e67978ec27df938cc` |
| `repro-transigen-2.log` (capture run) | `fe0f5e40a45fff69197856bfde7aa421ea3da2a25daa12e871c7d6ff47ed8fb2` |
| `repro-transigen-cudnn-off.log` | `8dedecc1d1615209c522b905ed8fd1c0cb5b53b43972d016e45df68520664e58` |
| `repro-transigen-author-flags.log` | `d4827c21b7827c52ed589bea996d217b3df1c221db9d60cba36ca844aa4cd07c` |
| `probe-transigen-loss.log` | `5a86e9b157d30f1e7e1c1364651327c6b8ff486bd654c6a8894866fe44314045` |
| `nonfinite-loss-inputs.pt` | `dc87cca40af4d7c60d8ba1e0dfb2856eb2ba54b696bd5ba50ef07f2652a059f5` |
| `transigen-contract-audit.md` | `ce6c39ac4d8c26b7f44991c669872a7aea55b9a804357b8955e6d0e5600d709c` |

The same directory preserves scripts `repro-transigen.py`,
`repro-transigen-capture.py`, `repro-transigen-cudnn-off.py`,
`repro-transigen-author-flags.py` and `probe-transigen-loss.py`, timed receipts,
formal execution status and10-second resource monitoring. This compact record
does not duplicate full traces/tensors. Preserve those payloads and the original
worktree; a new worktree does not contain them automatically.

## Engineering correction and continuation boundary

Author `utils.setup_seed` CUDA settings were missing: cuDNN enabledFalse,
benchmarkFalse, deterministicTrue. The isolated backend correction applies and
records these in the real train/replay lifecycle before native construction,
leaving CPU backend flags unchanged. Three lightweight tests observe those entry
points, including a fresh subprocess, while stopping before real native model
construction. They do not claim a successful scientific rerun. The full author
three-flag diagnostic also fails at epoch155, so this correction does not solve
the formal block.

Original completed CUDA runs remain valid evidence of their recorded original
execution and can be replayed using frozen code551a3c9. The new corrected replay
entry rejects CUDA checkpoints lacking its backend identity; this prevents
relabeling old evidence as having used corrected settings. Never overwrite the
original plan, weights or terminal evidence to bypass this identity boundary.

Reconstruction starts from frozen code551a3c9, the candidate plan and the explicit
inputs `/tmp/arc-tenv5k-inputs`, Means `/tmp/arc-tenv5k-means`, and immutable native
sources `/home/chang/arc/pert_model/{PRnet,TranSiGen}/source`. Matrix `--plan-only`
is read-only; `--accept-only` must reject this incomplete matrix. No formal retry
or resume is authorized by this engineering correction.

No clamp, gradient clip, learning-rate/dtype/seed change, early stopping shortcut,
batch omission or mean/MC inference change was made. A numerical-method variant
requires a new declared scientific identity and honest affected-run reruns after
user authorization; no demonstrated repair preserving the current frozen path
has been established. Stage3 and all external source bytes remain unchanged.
