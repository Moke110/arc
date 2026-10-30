# Download — CMonge

**method: repo (bundled data)** — upstream clone (MIT); the demo data and demo
embeddings are shipped in the repository and are not a separate download.

## Repository

```bash
git clone https://github.com/AI4SCR/Conditional-Monge pert_model/CMonge/source
git -C pert_model/CMonge/source checkout 195f5cd1822e7c7de7c7ee842a7e4082992ccc04
# vendored read-only: the nested .git directory was removed on 2026-10-19
```

## Bundled content

`source/data/` (`dummy_data.h5ad`, `lincs_trapnell.smiles`,
`trapnell_drugs_smiles.csv`), `source/models/` (demo embeddings/checkpoints) and
`source/assets/` are repository content — "bundled demo data, no download" as the
deleted `download/manual.md` said. No pretrained weights are published by design.

## Notes

Deployed model; per-dataset HVG output, RDKit structure + categorical
`drug-dose` conditioning. Run card `stage2/01_setup/execution/cmonge/`. That card's
`run.sh` used the model-local `.venv`; ticket 04 deleted it as a duplicate of
`stage2/_envs/cmonge` (see `stage2/_envs/README.md`) and repointed the card at
the shared bench env. Upstream bug: `configs/demo_config.yml` embeds the author's
macOS path `/Users/alicedriesen/Projects/conditional-monge`; patch a copy when
running the demo.
