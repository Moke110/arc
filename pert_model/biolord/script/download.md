# Download — biolord

**method: repo + dataset** — upstream clone (BSD-3-Clause); the demo/training h5ad
was fetched separately.

## Repository

```bash
git clone https://github.com/nitzanlab/biolord pert_model/biolord/source
git -C pert_model/biolord/source checkout b7688790e49728d7f8d0906b980a629de484b19b
# vendored read-only: the nested .git directory was removed on 2026-10-18
```

## Dataset

`source/data/GSM_new.h5ad` — 14,811 cells x 4,999 genes, 5 drugs x 7 doses
(A549), 178,640,710 B, sha256
`11b4363d46c00bcb9109c232b4f18411d01889d7a218f459a7d3feaa9b8f23da`.
It was downloaded on 2026-10-18 for the Stage 2 setup smoke; the old
`STAGE2_DATA.md` inventory did not record its upstream URL, so the author log
says so explicitly.

## Notes

Deployed model; no published weights (per-dataset training). Run card
`stage2/01_setup/execution/biolord/`.
