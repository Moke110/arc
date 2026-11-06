# CMonge first formal post-training failure diagnosis

No repository/source/config/environment/formal artifact edits. No retraining, test truth or test scoring. All probes used saved fold0/seed17 `best.msgpack` and predictor-safe training inputs, with 10GiB RAM/8GiB swap scope and time-v. Active fold0/seed42 remained independent.

## Tight reproduction

`stage2/_envs/cmonge/bin/python /tmp/tenv_cmonge_gate_probe_gpu.py` reproduces the exact failing four-observable predicate on CUDA in4.22s, peak1559656KiB RSS. Log/JSON adjacent to script. CPU variant `/tmp/tenv_cmonge_gate_probe.py` is a differential control (2.74s,909116KiB).

| observable | exact CUDA | CPU |
|---|---:|---:|
| adapter parameter change |1.0975415707|1.0975416899|
| baseline sensitivity |0.1000179052|0.1000003815|
| dose sensitivity |0|2.384185791e-7|
| reload max absolute error |0|0|

Both produce finite predictions. The CUDA gate fails solely because its selected checkpoint has exactly zero output dose effect on the single first training condition. The CPU effect is only float32 ULP scale, so it is not reliable evidence of useful dose dependence.

## Discrimination

Ranked hypotheses before follow-up: local GPU-rounded dose response; NumPy-vs-JAX reload mismatch; GPU inference nondeterminism/layout; a single unrepresentative condition versus global missing dose wiring.

`/tmp/tenv_cmonge_dose_probes.py` repeats exact CUDA gate and tests offsets without training. First condition remains identical for logdose offsets+2,+10,-1,-10. But three of twelve deterministic evenly spaced training conditions respond to logdose+1 (max changes5.125999e-6,1.072884e-5,1.545995e-5). Other nine respond0. Follow-up4.79s,1637624KiB RSS. JSON/log adjacent.

`/tmp/tenv_cmonge_context_probe.py` uses only NumPy parameter inspection: dose input coefficient in upstream second context Dense is0.1715848148, nonzero; first condition context preactivation changes0.332885→0.504470 for logdose+1. Thus dose reaches its native context layer. Do not attribute this to GELU saturation at that first context layer (preactivations are modest). Downstream dependence is weak/locally numerically invisible. Global missing dose wiring and checkpoint reload defect are falsified by the saved model's other-condition response and exact reload.

## Finding and correction boundary

The post-fit implementation gate at cmonge_adapter.py119 incorrectly treats zero local output sensitivity on one fixed condition as an execution failure. Such a selected checkpoint may be a poor/locally insensitive predictor; the comparison should retain that behavior rather than censor its run. Admission previously established actual input participation. Qualification checks can remain strict; formal checkpoint metadata should retain zero sensitivities explicitly as diagnostics. Keep adapter-change, finite outputs, valid shape/identity and exact reload checks strict. Do not retrain or change native objective/architecture/dose pathway merely to force a sensitivity score.

## Suggested regression seam (proposal only)

A deterministic cheap CLI regression can use the existing lightweight adapter artifact seam with a tiny fixture adapter calling a newly factored real checkpoint-observable validation policy. Supply finite predictions, changed adapter, exact reload and zero dose/baseline sensitivity: formal policy must allow completed artifacts with explicit zero diagnostics; qualification policy must reject. A malformed/nonfinite/unchanged-adapter/reload mismatch case must still reject. That tests the repaired policy through the artifact CLI without relying on a two-hour trained checkpoint or arbitrary synthetic training saturation. Keep existing real CMonge qualification regression intact for model integration. Parent owns design/fix/identity and retry-budget accounting.

## Implemented correction and native regression

Parent subsequently requested the fix. Unstaged changes are limited to `cmonge_adapter.py`, native `test_cmonge.py`, and CMonge README/FIDELITY documentation. No formal identities/queues/environments/source changed here.

A real native synthetic full18,164gene one-iteration/batch2 CLI case with finite dose1e-30 makes the dose-context response exactly zero while adapter parameters change. Old production gate rejects it (RED3.893s,720868KiB); corrected gate completes and validates finite2×18164predictions with dose_sensitivity0 and exact reload (GREEN4.002s,723320KiB). A four-gene native qualification companion still rejects the same locally flat input. Existing real native fidelity regression remains green. All three tests pass15.997s,wall19.11s,721376KiB peakRSS under10G/8Gscope. Logs and machine-readable record: `/tmp/tenv-cmonge-fix-record.json`.

The same two new artifact CLI regressions also pass through the actual pinned CUDA adapter entrypoint:2tests25.132s,wall28.19s,1362648KiB RSS under10G/8Gscope. The driver selects the public CUDA entrypoint without model mocks or production changes. Production adapter SHA256 remains99c71c28a6c47149b493c6c0c7aa5981c6ffd38de5dc4e5945c87e45ed7685eb.
