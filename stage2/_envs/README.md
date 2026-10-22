# Stage 2 environments (ticket 10/13/14/15)

Shared Stage 2 infrastructure: the per-model bench environments live here so
later stages can reuse them (originally under `03_bench_models/envs/`).  The
`biolord` and `cmonge` venvs moved with the directory; their recorded paths were
rewritten to `stage2/_envs/`, while the four `arc-stage2-bench-*` conda envs
remain in the conda root and are referenced by name.

The seven bench environments are isolated per model (ADR 0021); the shared
`arc` environment is not modified. `*.lock` is the final freeze (`pip freeze`
for conda envs and venvs, `conda list --explicit` when noted); the sha256 is
recorded here for provenance. The three ticket-13 envs were cloned from the
untouched bringup envs and given the CUDA wheel required by ticket 10.

| model | env | Python | framework | packages | lock sha256 |
|---|---|---|---|---|---|
| PRnet | `arc-stage2-bench-prnet` | 3.9.25 | torch 2.0.1+cu118 | 45 | `e79dc086b8886f578f143ba1f0b1f3f036f5a45225995d1b5032aaae6caaa9bc` |
| CPA | `arc-stage2-bench-cpa` | 3.9.25 | torch 2.0.1+cu118 | 104 | `8a7d8d373652128b0566e490f11a3c49cc0462ff03add4ff8615e8bee7ecaf25` |
| chemCPA | `arc-stage2-bench-chemcpa` | 3.9.25 | torch 2.8.0+cu128 | 97 | `7fa2808ed8c4af81205cb3f7aa5cc9ab175532ef2d76b9fc1487a207e4961dfa` |
| biolord | `stage2/_envs/biolord` (venv) | 3.13.13 | torch 2.12.0+cu130, scvi-tools 1.1.6 | 362 | `b4a5f6b8a787a1c49468ac208ae88ecf8ccd290f33028757517b2359befe0438` |
| PerturbNet | `arc-stage2-bench-perturbnet` | **3.7.16** | **torch 1.13.1+cu117, scvi-tools 0.7.1, numpy 1.19.5** | 78 | `b559d36a251dc9059012c4e5971c1c682245d64e73269dc9274fd157b0df0920` |
| Squidiff | `arc-stage2-bench-squidiff` | 3.13.16 | torch 2.14.1+cu130 | 34 | `74fdd59fe1c5c78155c37554bdfecaa1a5dd30c387b754c6eb9f6d924a98aace` |
| CMonge | `stage2/_envs/cmonge` (venv) | 3.11.15 | jax/jaxlib 0.4.38, **CPU backend** | 84 | `dc7dcd41c0c57ceaeec8f321dd18fef93bc6fb4319ef7fcae5e7b56779e19f7d` |

Rebuild commands (one-off; not counted in the 16 h/model budget):

```bash
conda create -y -n arc-stage2-bench-prnet  --clone arc-stage2-bringup-prnet
conda create -y -n arc-stage2-bench-cpa    --clone arc-stage2-bringup-cpa
conda create -y -n arc-stage2-bench-chemcpa --clone arc-stage2-bringup-chemcpa
conda create -y -n arc-stage2-bench-squidiff --clone arc-stage2-bringup-squidiff
# PRnet / CPA only: replace the CPU wheel with the CUDA build
pip install --force-reinstall --no-deps "torch==2.0.1+cu118" \
  --index-url https://download.pytorch.org/whl/cu118
# Squidiff: CUDA 13 build
pip install --force-reinstall --no-deps "torch==2.14.1+cu130" \
  --index-url https://download.pytorch.org/whl/cu130

# biolord: venv overlay on the arc interpreter
/home/chang/miniconda3/envs/arc/bin/python -m venv --system-site-packages \
  stage2/_envs/biolord
stage2/_envs/biolord/bin/pip install \
  scvi-tools==1.1.6 anndata==0.10.9 zarr==2.18.7 biolord==0.0.3 \
  descriptastorus==2.8.0 pandas_flavor==0.7.0 xarray==2024.3.0

# PerturbNet: the authors' faithful py3.7 pin (ticket 10)
conda create -y -n arc-stage2-bench-perturbnet python=3.7
conda run -n arc-stage2-bench-perturbnet pip install \
  "torch==1.13.1+cu117" --extra-index-url https://download.pytorch.org/whl/cu117
conda run -n arc-stage2-bench-perturbnet pip install \
  -r pert_model/perturbnet/requirements.txt   # TensorFlow 1.15 is not installed; see run card

# CMonge: the bringup uv venv copied to the bench path
cp -a pert_model/cmonge/.venv stage2/_envs/cmonge
```

Driver 610.88 supports cu117/cu118/cu128/cu130. PerturbNet runs the authors'
py3.7 / torch 1.13.1 / scvi-tools 0.7.1 pin; the ticket-10 py3.10 fallback was
not needed. TensorFlow 1.15 is deliberately absent: only the unused
normalised-expression cell-VAE branch imports it, and its CUDA 10 runtime
segfaults beside torch 1.13/cu117. CMonge stayed on CPU JAX after the CUDA 12
wheel install exceeded its timebox; the 10k-iteration recipe is ~32 min/fold on
CPU, inside the budget.
