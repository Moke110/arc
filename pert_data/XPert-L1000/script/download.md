# Download — XPert-L1000

**method: unknown** (version unknown, ADR 0030 spec §2). **license: unknown.**

No upstream download record survives for these three h5ads. They were produced
for the XPert L1000 pretrain recipe and moved into the data lake from the XPert
model library during the workspace refactor (ADR 0008: "the three L1000 pretrain
h5ads move to the dataset home"), then known under the dataset-home name `l1000-xpert` (2026-08-12), and now
live at `pert_data/XPert-L1000/source/data/`.

Files (as present in `source/data/`):

| file | bytes |
|---|---:|
| `l1000_mdmt_268022_pretrain.h5ad` | 4,362,329,844 |
| `l1000_mdmt_full_336852.h5ad` | 5,464,216,676 |
| `l1000_mdmt_example.h5ad` | 6,811,568 |

They are log2 expression space and are read by the workflow-owned pretrain YAMLs
(`stage1/xpert-fda-miri-v*/04_xpert_pretrain/config/l1000.yaml`). The matching
drug embeddings are workflow-owned (`stage1/input/`), not part
of this item. `local/L1000_summary.md` records what is known; the upstream
provenance (source URL / producing commit) is still to be filled in here and in
`author_log.json` (`version: "unknown"`).
