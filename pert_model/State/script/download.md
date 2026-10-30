# Download — State

**method: repo + released weights** — upstream clone plus the Arc Institute
released pretrained weights.

## Repository

```bash
git clone https://github.com/ArcInstitute/state pert_model/State/source
git -C pert_model/State/source checkout 9bbfe78a434a55205e4de834e1ea99f85f7a3add
# vendored read-only: the nested .git directory was removed on 2026-10-02
```

The clone shipped a 0-byte `.gitmodules` (not a live submodule); it was removed
as debris in ticket 04.

## Weights

`source/weights/` holds the released pretrained weights, materialized into the
repository on 2026-10-02 (ticket 02) and no longer a symlink into the external
weight archive:

- `SE-600M/model.safetensors` (2,860,276,816 B), `config.yaml`,
  `protein_embeddings.pt` (410,886,729 B);
- `ST-HVG-Tahoe/checkpoints/final.ckpt` (1,071,778,471 B), `pert_onehot_map.pt`,
  `config.yaml`, `var_dims.pkl`, `hvg_gene_names_recovered.json`,
  `data_module.torch`, batch/cell-type one-hot maps.

Licence: **Arc State Model Non-Commercial License** for the weights and
**CC BY-NC-SA 4.0** for the code. Non-commercial use only.

## Notes

Reference-only for Stage 2; the Tahoe one-hot limits cold-drug generalisation.
