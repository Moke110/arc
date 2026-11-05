# PerturbDiff native setup

Source is immutable `pert_model/PerturbDiff/source/`, pinned to
`31be89cb8cce5dddef0782bac29d5c748ce55931`. The local isolated overlay
`stage2/_envs/perturbdiff/` inherits arc without modifying it; `env.lock` records
the effective packages. No pretrained checkpoint or upstream dataset is used.

`run.sh` exercises native `PlModel`, CovEncoder, 8×512 Cross_DiT and native
GaussianDiffusion MSE + energy MMD on synthetic2,000-gene populations: two
optimization steps, save/load and native100-step DDIM. Synthetic one-hot conditions
and one cell type are native configuration options. This setup proof is separate
from the tenv jointly trained Morgan adapter and full18,164-gene admission.

Run from the repository root under a granted heavy-job slot:

```bash
systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=8G -- \
  /usr/bin/time -v bash stage2/01_setup/execution/perturbdiff/run.sh
```

Native setup PASS:2×3×2,000 finite output, native losses12.9897/14.4546,
checkpoint reload error0,100-step DDIM. `/usr/bin/time -v` peak RSS2,525,488KiB,
CUDA allocated peak2,021,450,240bytes, wall9.64s. `output/meta.json` records the
actual environment, checkpoint hash, shape and resource evidence. No tenv data
or formal benchmark scores were used.
