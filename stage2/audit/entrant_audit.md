# Stage 2 entrant audit — candidate models (ticket 02 re-run)

Capability gates come first: a whole-transcriptome output on the 20,570-gene universe, a cold-drug prediction from chemical structure, and a dose input for the cold-dose column. **Tier S** models may be trained locally on the fold's train wells (≤16 h for the full benchmark on the RTX 4070 Ti SUPER); **Tier L** models require a public pretrained checkpoint. Evidence was read from the cited papers and repos at audit time. This survey supersedes the pre-re-audit status lines in `03_*/README.md`; the verdicts below are the model-selection record.

| model | tier | status | cold-drug | dose | gene space | est h | license | checkpoint |
|---|---|---|---|---|---|---|---|---|
| chemCPA | S | admitted | yes | yes | 977 / 2,000 (public checkpoints); 20,570 only after retraining | 10 | MIT | https://f003.backblazeb2.com/file/chemCPA-models/chemCPA_models.zip |
| CPA | S | admitted | yes | yes | 20,570 (VAE gene dimension) | 9 | BSD-3-Clause | Google Drive folder linked from the CPA README ('Datasets and Pre-trained models') |
| biolord | S | admitted | yes | yes (ordered attribute) | 20,570 (VAE gene dimension) | 12 | BSD-3-Clause | — |
| PRnet | S | admitted | yes | yes | 20,570 (x_dimension is configurable; the published config is 978/5,000) | 8 | Apache-2.0 | in-repo checkpoint/lincs_best_epoch_all.pt and sciplex_best_epoch_all.pt; also Zenodo 10.5281/zenodo.13751384 |
| State | L | reference | no | yes (inside the Tahoe one-hot) | 2,000 HVG (1,398 in the common universe) | 16 | Arc Research Institute State Model Non-Commercial License | Arc Institute State pretrained weights (SE-600M + ST-HVG-Tahoe) |
| scGen | S | excluded | no | no | 20,570 (VAE gene dimension) | — | GPL-3.0 | — |
| CellOT | S | excluded | no | no | 20,570 (input dimension) | — | BSD-3-Clause | — |
| scPRAM | S | excluded | no | no | 20,570 (VAE gene dimension) | — | MIT | — |
| CFM-GP | S | excluded | no | no | 20,570 (input dimension) | — | GPL-3.0 | README lists pretrained models |
| TxPert | L | excluded | no | no | full transcriptome (gene nodes) | — | Recursion Non-Commercial End User License (CC BY-NC-SA-like) | in-repo weights for the top public-data variants |
| PertDiff (roadmap: MCDP) | S | excluded | yes | no | 978 landmark genes (then a linear expansion to 10,174) | — | none declared in the repo (Zenodo license not verified) | Zenodo 10.5281/zenodo.18427848 (model.pt) |
| LDM-Drug | S | excluded | yes | yes | 978 landmark genes only | — | Pre-publication Research License (view-and-study only; no use without written permission) | Google Drive (best_ldm.pt), linked from the README |
| scGPT | L | excluded | no | no | 1,200 input HVGs in the perturbation tasks (paper) | — | MIT | public scGPT checkpoints (figshare/Hugging Face) |
| scFoundation | L | excluded | no | no | ~19,000 genes (close to the universe) | — | Apache-2.0 | public scFoundation checkpoint (~19k genes) |

## Admitted (usable entrants)

- **chemCPA** (S tier, ~10 h) — chemical encoder generalises to unseen drugs and the model is dose-aware; the public checkpoints are L1000-only, so whole-transcriptome use requires Tier S retraining inside the 16 h budget. Recipe: public checkpoint proves load/forward on L1000; for a 20,570-gene CM run, reuse the chemical encoder and retrain the VAE decoder on the fold's CM train wells with drug = chemical embedding and condition = per-plate baseline; predict held-out drugs.
- **CPA** (S tier, ~9 h) — CPA's OOD result covers completely unseen drugs via a chemical encoder; per-dataset training is minutes-to-hours on the target GPU. Recipe: train CPA on each fold's CM train wells; reuse chemCPA's chemical encoder for the held-out drug (the CPA paper's unseen-drug result uses that extension); dose is a native covariate.
- **biolord** (S tier, ~12 h) — the paper's OOD split holds out nine unseen drugs and predicts them from chemical features; single VAE, tractable on one 16 GB card. Recipe: train biolord on each fold's CM train wells; drug = RDKit chemically informed embedding, dose = ordered attribute, baseline well = control cell; counterfactual prediction for held-out drugs.
- **PRnet** (S tier, ~8 h) — explicitly designed for novel compounds and dosages; Apache-2.0; tiny model, so a 20,570-gene re-train fits the budget. Recipe: set x_dimension=20,570 in the training config and train on each fold's CM train wells; drug = rFCFP chemical embedding, dose = native dose encoder; predict held-out drugs.

## Reference only (runnable, not a contender)

- **State** — perturbation identity is a fixed 1,138-way Tahoe drug-dose one-hot: 62/79 CM drugs are unseen and fall back to control, and there is no CM cell type. Kept as a reference row, not a contender.

## Excluded

- **scGen** — scGen transfers a latent perturbation shift to unseen cell types; it has no drug encoder and no dose input, so it cannot predict a held-out drug
- **CellOT** — CellOT learns a separate transport map per perturbation and generalises to unseen patients/cells, not unseen drugs; it has no chemical or dose input
- **scPRAM** — scPRAM's stated generalisation is to unseen cell types; the perturbation is a categorical embedding with no chemical structure or dose input
- **CFM-GP** — the vector field is conditioned on cell type only (gene perturbation across cell types); there is no drug identity or dose in the model
- **TxPert** — TxPert predicts genetic perturbations (CRISPR KO/OE); a drug would first have to be mapped to target genes, which changes the task and is not a chemical-structure cold-drug prediction
- **PertDiff (roadmap: MCDP)** — the model is 978-landmark by construction (attention shapes are hard-coded at 978) and the repo declares no license; both violate the Stage 2 gates
- **LDM-Drug** — trained and evaluated on 978 L1000 landmarks only, and the license forbids use without written permission; excluded until a 20,570-gene refactor and a usable license exist
- **scGPT** — the perturbation head is gene-token based with a ~1,200-HVG input cap and no chemical or dose conditioning; adding a drug encoder would be a new model, not a fine-tune
- **scFoundation** — a cell-embedding foundation model with no drug-conditioned perturbation head; using it would mean training a new drug encoder and decoder head, not fine-tuning the released model

## Shared adapter risks

- **Input scale.** The benchmark freezes log2(CPM+1), while chemCPA / CPA / biolord natively model counts. Each adapter must either invert to counts or switch to a Gaussian likelihood; the fold-0 pilot must show a sane reconstruction before the full run.
- **Bulk-as-cell.** DToxS and DToxS-iPSC are bulk 3'-DGE and gse244740 is pseudo-bulk. The single-cell VAEs treat each well as one cell; this is an approximation to document, not a hidden one.
- **Drug identity.** chemCPA / CPA need a chemical encoder, biolord needs RDKit features, PRnet uses rFCFP. All consume SMILES, but the encoder must be pinned and its gene-independence verified.
- **RAM, not VRAM.** The host has 15 GB RAM and 16 GB VRAM. The precomputed bench matrices are ~0.4 GB, so training fits; raw h5ad loading must be streamed/chunked.

## What would change a verdict

- A State adapter that re-initialises the perturbation embedding could improve cold-dose, but cold-drug remains impossible without chemical initialisation.
- LDM-Drug / PertDiff could re-enter after a 20,570-gene refactor and a usable license; that is a new model, not the published checkpoint.
- scGPT / scFoundation could re-enter only with a drug-conditioned head trained on CM data — a research project, not a fine-tune.
