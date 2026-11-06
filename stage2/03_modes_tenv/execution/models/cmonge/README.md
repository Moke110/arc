# CMonge runtime admission

[Ticket #50](https://github.com/Moke110/arc/issues/50): **admitted on CUDA** against
the frozen tenv full 18,164-gene contract. The final 100-step native batch512 profile
produced finite 1401 × 18164 predictions, exact checkpoint reload, trained Morgan
adapter parameters, and nonzero baseline/dose sensitivity. No test score was
computed. See `admission.json` and `profile_gpu512_*` evidence.
Qualification retains strict input-participation checks. Formal trained
checkpoints report single-condition baseline/dose sensitivities as diagnostics,
including exact zeros; adapter change, exact reload and valid finite predictions
remain mandatory. See `FIDELITY.md` for this distinction.

The native objective is Sinkhorn divergence plus Euclidean Monge gap. Validation
samples one held-out condition every 100 iterations using an independent per-run
seeded RNG; checkpoint history records its group ID. Formal configuration keeps
the author's 10,000 iterations and batch 512 in `formal.json`. Bounded profiling
used 100 iterations, not a completed formal run.

GPU batch 512 profile: 72.350 s training/validation, 0.750 s inference; peak RSS 1,695,352 KiB,
conservative sampled process-tree peak 1.711 GiB; JAX allocator peak 0.649 GiB and
pool peak 1 GiB. Workstation GPU total peaked 3.316 GiB, including background use;
the scheduling profile conservatively uses that total. Linear extrapolation is
2.010 h per 10,000-step run, within the enforced four-hour training allowance;
this remains an estimate. Nine formal runs remain separate.

The first formal run exposed a post-training validation bug: its selected
checkpoint was numerically insensitive to dose on the single fixed probe
condition, although other training conditions responded to dose. Exact reload,
adapter updates and baseline response passed. Qualification still requires
nonzero input sensitivity; formal runs record these local sensitivities as
diagnostics, including zero, rather than rejecting a fitted predictor on that
basis. Finite diagnostics, adapter updates, exact reload and prediction integrity
remain required. Architecture, native objective, checkpoint selection, random
draws, dose conditioning and training settings are unchanged.

The original failed attempt and the gracefully interrupted second attempt are
preserved under the formal repair record. Corrected fits carry their elapsed
training time into the same four-hour cumulative budgets. Only CMonge's nine
expected code identities are revised; the entrant set and all other identities
remain fixed. See `../../formal/cmonge_repair_record/` for diagnosis, old
identities, checkpoints' hashes and migration evidence.

The CPU native 100-step profile took 787.053 s, projecting 21.86 h; its batch 512 OOM
and unique batch-256 retry remain preserved. GPU batch 512 was separately measured with further
OOM retry disabled and succeeded. GPU batch 256 is a measured fallback, not a third
CPU retry. The first CUDA qualifier found the system's old CUDA 12.0 assembler;
one recorded environment repair pins overlay nvcc 12.9.86 through PATH and XLA's
CUDA tools directory. No dependency or upstream source was replaced during that
repair. See `CUDA.md`, `FIDELITY.md` and the initial failure evidence.

Formal adapter command (the benchmark appends `--request`):

```bash
stage2/_envs/cmonge/bin/python \
  stage2/03_modes_tenv/execution/script/cmonge_gpu_adapter.py \
  --config stage2/03_modes_tenv/execution/models/cmonge/formal.json
```

The GPU wrapper loads the isolated CUDA overlay and runs the same native
adapter. Its sidecar verifies live package metadata and 341 actual file hashes,
including explicitly allowed NVIDIA-only symlink targets. CPU environment and
upstream 85-file source manifest remain unchanged. Generated unmanifested native
bytecode was relocated; both entrypoints now suppress bytecode writes.

The accepted subprocess/artifact test exercises native objective, seeded
multi-condition validation replay, joint adapter updates, input participation,
checkpoint reload and explicit device memory evidence. Last run: PASS, 7.438 s,
642,568 KiB peakRSS. Run through the repository memory gate:

```bash
systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=8G -- /usr/bin/time -v \
  python -m unittest discover -s stage2/03_modes_tenv/execution/tests -p test_cmonge.py -v
```

`PAUSE.md` and `pause_record.json` preserve the historical pause state; the
current runtime outcome is `admission.json`. All formal execution must use the
resource scheduler. Qualification output arrays/checkpoints remain local-only.
