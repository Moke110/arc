# DECANT setup smoke

**Status:** synthetic 2,000-gene forward with two dose values

- source: `pert_model/DECANT/source/` (commit `e50642a9`), MIT
- weights/data: none — the repository is a reproduction skeleton
- env: `arc-stage2-bringup-decant` (Python 3.13, CPU torch, numpy, pandas)
- command:

```bash
bash stage2/01_setup/execution/decant/run.sh
```

Expected output: `output/meta.json`, `output/run.log`. `pred_gene_dose{1,100}.npy`
from the 2026-10-19 run were deleted in ticket 06; a rerun recreates them.

Caveat: this proves the model and dose context run at 2,000 genes; it does not
validate mechanism recovery, which needs matched treated/control cell bags.

**Result (2026-10-19):** passed. 162,384-parameter model; forward produced (4, 2,000) at both doses, with dose effect L1 = 0.00495.
