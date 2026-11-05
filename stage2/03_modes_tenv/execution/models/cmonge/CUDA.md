# CMonge isolated CUDA runtime

The native 100-step, batch256 CPU profile took 787.053 seconds for training and
one sampled validation event. Linear extrapolation to the unchanged 10000-step
recipe is 21.86 hours, beyond the four-hour training allowance. The upstream
trainer has a fixed iteration range and no early stopping. See
`native_resource_assessment.json` and `profile_native_*` evidence. This resource
profile is qualification evidence, not a formal run or test score.

The CUDA overlay preserves JAX/jaxlib0.4.38 and the original isolated CPU
environment. Matching CUDA12 plugin/PJRT packages and nvcc12.9.86 are installed
under `pert_model/CMonge/local/tenv_cuda/`. Nine NVIDIA-only component directories
and their package metadata are symlinked from the existing chemCPA CUDA12.8
runtime. No chemCPA Python model packages are imported, and its environment is
not modified. Each reused version satisfies the JAX0.4.38 plugin's published
`with_cuda` requirements. The workstation driver617.42 supports CUDA12.

The complete pinned package set is `cuda_overlay.lock`; exact reused paths are
in `cuda_overlay_reused_libraries.json`; `cuda_overlay.sha256` hashes341 files,
including native binaries reached through the symlinks. Verify from repo root:

```bash
sha256sum -c stage2/03_modes_tenv/execution/models/cmonge/cuda_overlay.sha256
```

The initial broad dependency download was stopped before installation because
large CUDA wheels were slow to download. The needed three packages had already
entered the package cache. The successful bounded installation used:

```bash
uv pip install --offline --no-deps --python stage2/_envs/cmonge/bin/python \
  --target pert_model/CMonge/local/tenv_cuda \
  jax-cuda12-plugin==0.4.38 jax-cuda12-pjrt==0.4.38 \
  nvidia-cuda-nvcc-cu12==12.9.86
```

On a fresh machine, the entire overlay can instead be installed from the pinned
`cuda_overlay.lock` using the same `--target` path without `--offline`; the
version pins replace local symlink reuse. Actual installed files must then be
rehashed before establishing new run identities. Logs retain the bounded
installation (44376KiB peakRSS) and hash-generation measurements.

`execution/script/cmonge_gpu_adapter.py` selects the CUDA backend and inserts
only the isolated overlay on its import path, then runs the same native adapter.
Its identity sidecar pins the wrapper, adapter, upstream source, CPU environment
lock, CUDA lock and full overlay file manifest. CUDA setup alone is not runtime
admission: full18164-gene train/save/reload/predict and formal-batch throughput
qualification are separate required gates.

Primary references: [official JAX installation documentation](https://docs.jax.dev/en/latest/installation.html),
[matching CUDA plugin release](https://pypi.org/project/jax-cuda12-plugin/0.4.38/).
