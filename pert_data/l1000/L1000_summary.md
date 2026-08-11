# L1000 pretrain dataset (pert_data/l1000)

The L1000 landmark-gene perturbation pretrain data used by the XPert pretrain
recipe (`workflow/xpert-fda-miri-v1/03_xpert_pretrain/l1000.yaml`). Moved here
from the model library during the workspace refactor (ADR 0008): like
`pert_data/dtoxs/`, this is a static dataset home, not a model artifact.

- `data/l1000_mdmt_268022_pretrain.h5ad` — pretrain split (268,022 samples).
- `data/l1000_mdmt_full_336852.h5ad` — full pretrain set (336,852 samples).
- `data/l1000_mdmt_example.h5ad` — small example subset.

All three are log2 expression space (no z-score conversion in training). The
matching drug embeddings (UniMol array, SMILES mapping, author HG reference)
are workflow-owned in `workflow/xpert-fda-miri-v1/input/`.
