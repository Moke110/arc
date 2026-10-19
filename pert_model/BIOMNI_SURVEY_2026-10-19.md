> Archived from the `biomni` subagent run of 2026-10-19. It is the raw survey
> behind `SURVEY_2026-10-19.md`. Caveat: a few name-level negative claims
> (e.g. `BioInfoVAE`, `MCDP`) were not reconciled with `reference/markdowns/`;
> use the gate analysis, repository URLs, and licence checks, not those
> name-level claims.

# Survey: open-source whole-transcriptome drug-perturbation ("virtual cell") models, 2023–2026

Scope: models that map **(drug structure / SMILES) → (genome-scale expression profile)**, ideally
supporting **dose** as an input and generalising to **cold (held-out) drugs**. Screened against four gates:

- **G1** whole-transcriptome output (ideally ≈20k genes, not just 978 L1000 landmarks)
- **G2** cold-drug generalisation from chemical structure (not a learned per-drug embedding)
- **G3** dose as input / covariate
- **G4** public code, trainable on one 16 GB GPU in ≲16 h

Legend for evidence: **[C]** = confirmed by me directly (paper full text / repo file fetch / GitHub API);
**[U]** = unverified (asserted by a paper or README but I could not independently confirm);
**[I]** = inference.

---

## 1. Ranking table (ordered by relevance to the four gates)

| # | Model | Paper (year, DOI) | Code + license | Predicts | Gene space | Cold-drug mechanism | Dose | Weights | Verdict |
|---|-------|-------------------|----------------|----------|-----------|---------------------|------|---------|---------|
| 1 | **Perturbed-GE-LDM** | "Predicting Condition-Aware Drug-Induced Transcriptional Responses via a Latent Diffusion Model" (Kim & Yoo; year **[U]**, repo 2026) | `github.com/bmil-jnu/Perturbed-GE-LDM` **[C]**; "Pre-publication Research License" **[C]** (restrictive, non-OSI) | bulk L1000 drug | 978 **[I]** (from PRnet preprocessing) | MolFormer pretrained SMILES embeddings **[C]** | **Yes — explicit dose + time + cell-type conditioning** **[C]** | Pretrained LDM checkpoint on Google Drive **[C]** (link in README) | **Strongly relevant** — closest functional match to all four gates; caveats: restrictive licence, 978 genes, no peer-reviewed paper found |
| 2 | **WAVE** | "Learning chemical-induced gene expression perturbations with WAVE", *Nat Commun* 2026, `10.1038/s41467-026-77645-3` **[C]** | `github.com/ZJUFanLab/WAVE` **[C]**; **MIT** **[C]** (article itself CC BY-NC-ND) | bulk + single-cell drug | **978 or 3000** (configurable) **[C]** | Morgan fingerprint → DrugNN MLP (2048→512→128) **[C]** | **No** — trained at fixed 10 µM **[C]** (reviewers raised this) | Pretrained model + datasets on FigShare `29650832` **[U]** (author response in peer-review file) | **Strongly relevant** — cleanest β-VAE baseline with released weights; gene space and dose are the limits |
| 3 | **PerturbNet** | "PerturbNet predicts single-cell responses to unseen chemical and genetic perturbations", *Mol Syst Biol* 2025, `10.1038/s44320-025-00131-3` **[C]** | `github.com/welch-lab/PerturbNet` **[C]**; **GPL-3.0** **[C]** | single-cell **distributions**, chemical + genetic | **978** (LINCS-Drug) / **5087** (sci-Plex) **[C]** | ChemicalVAE over one-hot SMILES, pretrained on 250k ZINC molecules **[C]** | **Yes — dose & cell type as covariates in the cINN mapping** **[C]** | No author weights in repo **[C]**; third-party `huggingface.co/cyclopeta/PerturbNet_reproduce` **[C]** | **Strongly relevant** — only model here that predicts a *distribution* for unseen drugs; 978–5087 genes |
| 4 | **DeepICER** | "DeepICER: A deep learning framework for predicting compound-induced gene expression profiles", *Acta Pharm Sin B* 2026, `10.1016/j.apsb.2026.01.046` **[C]** | `github.com/GreatChenLab/DeepICER` **[C]**; **no LICENSE file** **[C]** (all rights reserved by default) | bulk L1000 drug + cell line | **978** **[C]** | GCN over RDKit 2D molecular graph from isomeric SMILES **[C]** | **Yes — dose (1/10/20 µM) and duration (24/48 h) as explicit bilinear-attention inputs** **[C]** | None released; web server only **[C]** | **Strongly relevant** — exactly your dose+duration formulation; licensing is a blocker for redistribution |
| 5 | **ExPO** | "ExPO: an exposure-conditioned neural operator for L1000 signature prediction", *J Cheminform* 18:94, 2026, `10.1186/s13321-026-01226-1` **[C]** | `github.com/MLBC-lab/ExPO` **[C]**; **MIT** **[C]** | bulk L1000 z-scores | **978** **[C]** | ChemBERTa-2 SMILES encoder + LoRA; **scaffold-held-out** evaluation **[C]** | **Yes — continuous (dose, time) as DeepONet trunk over Fourier features; arbitrary dose–time pairs** **[C]** | None stated (code + benchmark dataset) **[U]** | **Strongly relevant** — best dose–time surface handling and the most honest cold-drug (scaffold) split; 978 genes only |
| 6 | **CMonge** (Conditional Monge Gap) | "Conditional Monge Gap enables generalizable single-cell perturbation modelling", *Nat Mach Intell* 2026, `10.1038/s42256-026-01242-8` (arXiv 2504.08328) **[C]** | `github.com/AI4SCR/Conditional-Monge` **[C]**; **MIT** **[C]**; also PyPI `cmonge`, Zenodo `10.5281/zenodo.17854338` **[C]** | single-cell distributions; drug, dose, combinations; also 40-plex protein (4i) | per-dataset HVG (configurable) **[C]** | RDKit structure embedding **[C]** (MoA embedding works better — honest negative result) | **Yes — dose is a first-class covariate** **[C]** | No pretrained weights — by design (~22.8k params, 10–50 MB GPU, ~10.6 s per new condition) **[C]** | **Strongly relevant methodologically** — cheapest to reproduce (minutes, not 16 h), true cold-drug, but per-dataset retraining and no 20k output |
| 7 | **DECANT** | "DECANT: decoupling mechanism from context in single-cell drug perturbation representation", *Bioinformatics* 2026, `10.1093/bioinformatics/btag662` **[C]** | `github.com/bliulab/DECANT` **[C]**; **MIT** **[C]**; Zenodo `10.5281/zenodo.21216567`; web server `bliulab.net/DECANT` **[C]** | single-cell drug (Sci-Plex3), gene- and program-level | fixed HVG space, configurable **[C]** | SMILES/molecular graph encoder + LINCS teacher distillation **[C]** | Dose + time treated as **context nuisance to be removed**; they benchmark stability *across* dose **[C]** (not a generative dose input) **[I]** | **No** — "large trained weights … are not included" **[C]**; repo is a reproduction skeleton | **Possibly relevant→strong** — right problem (mechanism vs context), but no weights, no data, dose not generative |
| 8 | **Squidiff** | "Squidiff: predicting cellular development and responses to perturbations using a diffusion model", *Nat Methods* 2026, `10.1038/s41592-025-02877-y` **[C]** | `github.com/siyuh/Squidiff` **[C]**; **BSD-3-Clause** **[C]** (+ `Squidiff_reproducibility`) | single-cell; drug, gene perturbation, differentiation | **configurable `--gene_size` / `--output_dim`** (200 used for sci-Plex; 500 default) **[C]** | optional drug-structure conditioning flag `--use_drug_structure True` **[C]**; encoder details **[U]** | No explicit dose input found **[U]/[I]** | Checkpoints on FigShare `10.6084/m9.figshare.27948633` **[C]** | **Possibly relevant** — only model here whose gene space is genuinely open-ended; drug-structure path is secondary |
| 9 | **MiTCP** | "Predicting transcriptional changes induced by molecules with MiTCP", *Brief Bioinform* 2024, `10.1093/bib/bbaf006` **[C]** | `github.com/yangkaiyuan1027/MiTCP` **[C]**; **no LICENSE** **[C]** | bulk L1000 CTPs | **978** **[C]** | GNN over molecular graph + gene co-expression graph **[C]** | No **[C]** | Zenodo `13991478` **[U]** (cited in paper) | **Possibly relevant** — simple, small, but no dose and no licence |
| 10 | **MultiDCP** *(likely what "MCDP" meant)* | "Deep learning prediction of chemical-induced dose-dependent and context-specific multiplex phenotype responses…", *PLoS Comput Biol* 2022, `10.1371/journal.pcbi.1010367` **[C]** | repo **unverified** (I did not confirm an official repository) | bulk L1000 + viability, dose-dependent | 978 **[I]** | chemical structure features **[U]** | **Yes — dose-dependent by design** **[C abstract]** | unknown | **Possibly relevant** — the 2022 dose-aware ancestor; code availability needs a manual check |
| 11 | **CIGER** | "Chemical-induced gene expression ranking and its application to pancreatic cancer drug repurposing", *Patterns* 2022, `10.1016/j.patter.2022.100441` **[C]** | `github.com/pth1993/CIGER` **[C]**; **no LICENSE** **[C]** | bulk L1000 **rankings** | 978 **[C]** | SMILES → structure encoder **[C]** | No | Not released **[U]** | **Possibly relevant** — ranking objective is useful for connectivity scoring, but 2022-era and unlicensed |
| 12 | **DeepCE** | "A deep learning framework for high-throughput mechanism-driven phenotype compound screening…", *Nat Mach Intell* 3:247, 2021, `10.1038/s42256-020-00285-9` **[C]** | `github.com/pth1993/DeepCE` **[C]**; Zenodo `10.5281/zenodo.3978774`; **no LICENSE file** **[C]** | bulk L1000 | 978 **[C]** | GCN over molecular graph + multi-head substructure–gene attention **[C]** | No | No weights; releases 11,179 DrugBank predicted profiles **[C]** | **Possibly relevant** — canonical cold-drug baseline (already used as such by XPert/ExPO) |
| 13 | **BaiZe** | "A multi-view conditional diffusion framework for simulating cellular responses across perturbation contexts" — paper DOI/year **[U]** | `github.com/simonzqw/BaiZe` **[C]**; **no LICENSE** **[C]** | single-cell; genetic + drug + dose + time + species (+ATAC) | configurable (dataset gene space) **[C]** | drug structure from SMILES via RDKit into a shared condition encoder **[C]** | **Yes — dose** **[C]** | None released **[U]** | **Possibly relevant** — architecturally the closest to "everything you want" (Δx formulation, multi-context), but 2★, no paper, no licence |
| 14 | **X-Pert** | "Unified multimodal learning enables generalized cellular response prediction to diverse perturbations", bioRxiv 2025.11.13.688367 **[C]** | `github.com/Chen-Li-17/X-Pert` **[C]**; **MIT** **[C]** | single-cell + bulk; genetic/chemical/combinatorial | scGPT gene vocabulary **[U]** | "Perturbation Perceiver" over heterogeneous perturbation encodings **[C]** | **Yes — quantitative perturbation strength / dose** **[C]** | None **[U]** | **Possibly relevant** — dose-aware unified model; depends on scGPT, so 16 GB is tight |
| 15 | **scDPR** | "Deep generative framework for modeling single-cell drug perturbation response", *Neural Networks* 202:109005, 2026, `10.1016/j.neunet.2026.109005` **[C]** | **No public repo found** — searched GitHub by name and by title; nothing **[C]** | single-cell drug (L1000 + sci-Plex3) | **[U]** | drug molecular features + OT causal graph **[C abstract]** | **Yes — dosage information** **[C abstract]** | None found | **Not usable now** — paper claims SOTA over chemCPA on unseen compounds but ships no code |
| 16 | **CellPert** | "Graph Variational Learning with Biological Priors Enables Cross-Dataset Perturbation Prediction" — venue/year **[U]** | `github.com/QSong-github/CellPert` **[C]**; **MIT** **[C]** | single-cell chemical perturbation, cross-dataset | dataset gene space + STRING PPI graph **[C]** | **Not structure-based** — transfers a reference set's latent displacement **[C]** → cold drug only via nearest reference compound **[I]** | Not modelled **[U]** | HuggingFace checkpoint + LINCS latents + Tahoe mini **[C]** | **Possibly relevant** — good engineering (released weights, tutorials), wrong generalisation mechanism |
| 17 | **GPO-VAE** | "GPO-VAE: modeling explainable gene perturbation responses utilizing GRN-aligned parameter optimization", *Bioinformatics* 2025, `10.1093/bioinformatics/btaf256` **[C]** | `github.com/dmis-lab/GPO-VAE` **[C]**; licence `NOASSERTION` **[C]** | single-cell **genetic** perturbation | **8,563 / 8,749 / 17,035 genes** (full measured set) **[C]** | none (gene identity) | No | No weights | **Possibly relevant only as an architecture reference** — proof that a VAE can decode ~17k genes on one GPU; not a drug model |
| 18 | **scVIDR** | "Single Cell Variational Inference of the Dose Response" — *Patterns*, DOI **[U]** | `github.com/BhattacharyaLab/scVIDR` **[C]**; **MIT** **[C]** | single-cell chemical perturbation, multi-dose, across cell types | dataset gene space **[C]** | **chemical identity as a categorical covariate — no structure encoder** **[C]** | **Yes — single- and multi-dose models** **[C]** | Pretrained models listed in README **[C]** | **Possibly relevant** — best-in-class dose handling, but cannot do a genuinely novel compound |
| 19 | **W1OT** | "Fast and scalable Wasserstein-1 neural optimal transport solver for single-cell perturbation prediction", *Bioinformatics* 2025, `10.1093/bioinformatics/btaf253` (arXiv 2411.00614) **[C]** | `github.com/poseidonchan/w1ot` **[C]**; **GPL-3.0** **[C]** | single-cell, per-condition OT map | configurable (anndata) **[C]** | **None** — one transport map per perturbation; no chemical encoder **[C]** | Not modelled | No | **Not relevant for cold drug** — a fast CellOT replacement; useful only as a solver inside another model |
| 20 | **CellCap** | "Modeling interpretable correspondence between cell state and perturbation response with CellCap", bioRxiv 2024.03.14.585078 **[C]** | `github.com/broadinstitute/CellCap` **[C]**; **BSD-3-Clause** **[C]** | single-cell, genetic **and** compound | dataset gene space | perturbation is a **categorical index** `p` with a relevance matrix `H_pq` — no structure **[C]** | No | No weights **[U]** | **Not relevant for cold drug** — interpretable response-program decomposition; good for MoA readout, not for novel compounds |
| 21 | **scCausalVI** | "scCausalVI disentangles single-cell perturbation responses with causality-aware generative model", *Cell Syst* 2025, `10.1016/j.cels.2025.101443` **[C]** | `github.com/ShaokunAn/scCausalVI` **[C]**; **MIT** **[C]** | single-cell, treatment vs control (IFN-β, COVID) | dataset gene space | none | No | No | **Not relevant** — causal disentanglement, not structure→expression |
| 22 | **scPerb** | "scPerb: Predict single-cell perturbation via style transfer-based variational autoencoder", *J Adv Res* 2025, `10.1016/j.jare.2024.10.035` **[C]** | `github.com/QSong-github/scPerb` **[C]**; **MIT** **[C]** | single-cell (Perturb-seq) | dataset gene space | style-transfer latent, no structure | No | No | **Not relevant** — genetic perturbation |
| 23 | **MultiPert** | "MultiPert: An adversarial alignment and dual attention framework for single-cell multi-omics perturbation prediction", *PLoS Comput Biol* 2026, `10.1371/journal.pcbi.1014054` **[C]** | `github.com/MengyuanZhaoo/MultiPert` **[C]**; **no LICENSE** **[C]** | single-cell RNA + protein (THP-1, kidney) | dataset gene space | perturbation as condition vector, no structure | No | Zenodo `7041849` **[U]** | **Not relevant** — genetic/multi-omics, not chemical |
| 24 | **CRADLE-VAE** | "CRADLE-VAE: Enhancing Single-Cell Gene Perturbation Modeling with Counterfactual Reasoning-based Artifact Disentanglement", AAAI 2025, `10.1609/aaai.v39i15.33695` (arXiv 2409.05484) **[C]** | `github.com/dmis-lab/CRADLE-VAE` **[C]**; **no LICENSE** **[C]** | single-cell genetic | dataset gene space | none | No | No | **Not relevant** — QC-artefact disentanglement |
| 25 | **SAMS-VAE** | "Modelling Cellular Perturbations with the Sparse Additive Mechanism Shift VAE", NeurIPS 2023, arXiv 2311.02794 **[C]** | `github.com/insitro/sams-vae` **[C]**; **CC (non-commercial)** **[C]** | single-cell genetic | dataset gene space | sparse global latent per perturbation — no structure | "dosage" = binary indicator vector only **[C]** | No | **Not relevant** — conceptual ancestor of CPA; no chemical generalisation |
| 26 | **CellBox** | "CellBox: Interpretable Machine Learning for Perturbation Biology…", *Cell Syst* 2021, `10.1016/j.cels.2020.11.013` **[C]** | `github.com/sanderlab/CellBox` (≡ `dfci/CellBox`) **[C]**; **MIT** **[C]**; PyTorch port `sanderlab/cellbox_torch` **[C]** | **RPPA proteomics** (99 nodes), ODE dynamics | n/a — 82 proteins, not transcriptome **[C]** | none | drug concentration as experimental condition only | No | **Not relevant** — proteomic ODE model, wrong readout |
| 27 | **Dr.VAE** | "Dr.VAE: improving drug response prediction via modeling of drug perturbation effects", *Bioinformatics* 2019, `10.1093/bioinformatics/btz158` **[C]** | `github.com/rampasek/DrVAE` **[C]**; **GPL-3.0** **[C]** | drug **sensitivity** (viability), not profile generation | n/a **[C]** | none | no | No | **Not relevant** — predicts response labels, not transcriptomes; repo unmaintained since 2019 |
| 28 | **TT-WOPT ("TensorTrain")** | "Predicting drug-induced transcriptome responses of a wide range of human cell lines by a novel tensor-train decomposition algorithm", *Bioinformatics* 2019, `10.1093/bioinformatics/btz313` **[C]** | No public code found **[C]** | bulk L1000 tensor completion | 978 **[C]** | **none — drug is a tensor index, so cold drugs are impossible by construction** **[C]** | time axis, not dose | n/a | **Not relevant** — cannot generalise to a new compound |
| 29 | **LPM** | "In silico biological discovery with large perturbation models", *Nat Comput Sci* 2025, `10.1038/s43588-025-00870-1` **[C]** | repo **[U]** (not verified) | heterogeneous perturbations, decoder-only P/R/C | readout-defined | **perturbation is a symbolic in-vocabulary token** **[C] → no cold drug** | no | unknown | **Not relevant** — elegant, but explicitly limited to in-vocabulary perturbations |
| 30 | **AROMA** | "AROMA: Augmented Reasoning Over a Multimodal Architecture for Virtual Cell Genetic Perturbation Modeling", ACL Findings 2026, `10.18653/v1/2026.findings-acl.1594` (arXiv 2604.20263) **[C]** | `github.com/blazerye/AROMA` **[C]**; **no LICENSE** **[C]**; HF `blazerye/AROMA` **[C]** | **ternary classification** (up/down/ND) per (perturbation gene, target gene) | n/a — per-gene classification **[C]** | gene embeddings + KG + ESM-2 | n/a | Yes (HF) **[C]** | **Not relevant** — not a transcriptome generator, genetic only |
| 31 | **Stable-Shift** | "Stable-Shift: Predicting Transcriptional Responses of Unseen Gene Perturbations Using GNNs with Biological Priors", ACM-BCB 2026, `10.1145/3807503.3820871` (arXiv 2606.24940) **[C]** | likely `github.com/Sajib-006/PerturbGraph` **[U]**; MIT | pseudo-bulk **genetic** perturbation | 8,563 genes (K562) **[C]** | none — STRING/GO/graph priors | no | unknown | **Not relevant** — the "cold" axis is cold *gene*, not cold drug |
| 32 | **PertFlow** | "Generating Joint Transcriptomic and Morphological Responses to Drug Perturbations via Rectified Flow", bioRxiv 2026.02.02.703189 **[C]** | `github.com/wangmengbo/PertFlow` — "**will be** available" **[C] → unverified/empty** | bulk RNA-seq + Cell Painting | **[U]** | drug metadata conditioning **[C abstract]** | **[U]** | promised | **Watch list** — not yet reproducible |

### Names from your list I could **not** substantiate

| Name | Finding |
|---|---|
| **BioInfoVAE** | No paper, preprint, or repository found in Europe PMC, arXiv, or GitHub. **[C]** — likely spurious. |
| **GeneRegDrug** | Same: zero hits anywhere. **[C]** — likely spurious. |
| **CyclePert** | Zero hits in Europe PMC/arXiv/GitHub; only unrelated CR9000X manual noise. **[C]** — likely spurious. |
| **HeteroPert** | Zero hits. **[C]** — likely spurious. |
| **MCDP** | No perturbation-biology model by this exact acronym. Two possibilities: (a) a typo for **MultiDCP** (row 10); (b) **MCDP** = "Modality-Composable Diffusion Policy" (arXiv 2503.12466) — a **robotics** paper, not biology. **[C]** I found no biological MCDP. **PertDiff**, the other half of your "MCDP/PertDiff" entry, *is* real (*Bioinformatics* 2026, `10.1093/bioinformatics/btag482`, 978 genes, MolGNet structure encoder, LLM cell semantics, conditional diffusion) and you already know it. |
| **TensorTrain** | Real but 2019 and not a neural model — see row 28 (TT-WOPT). |
| **VEGA** | The only single-cell VEGA is Seninge et al., *Nat Commun* 2021 (`10.1038/s41467-021-26017-0`), an interpretable VAE for **network activity inference**, not perturbation prediction. **[C]** Not relevant. |
| **scPILOT** | Real: "Predicting Single-Cell Perturbation Responses Across Biological Contexts With a Deep Generative Model Integrating Optimal Transport", *Adv Sci* 2026, `10.1002/advs.77461`. It **transfers previously observed perturbations across contexts** — no chemical encoder, so no cold drug. Repo not verified. **Not relevant.** |
| **scDPR** | Real paper, **no code** (row 15). |

---

## 2. Short notes on the ones that matter

### The uncomfortable headline: nobody hits all four gates
No drug-structure model in this survey emits ~20,000 genes. The ceiling among *drug* models is
**978** (L1000) or **3000** (WAVE's extended panel). Models that do emit 8k–17k genes
(GPO-VAE, Stable-Shift, PerturbNet-on-sci-Plex at 5087) are **genetic-perturbation** models with no
chemical encoder. The only models with genuinely open gene dimensionality are the configurable-HVG
single-cell ones (Squidiff, BaiZe, CMonge, State), and they top out at whatever the dataset measured
(typically 2k–5k HVG). **[C for each individual claim; the synthesis is [I]]**

Practical consequence: if you need ~20k-gene output, you will be **composing** — e.g. a 978-landmark
drug model (WAVE / ExPO / DeepICER) plus an imputation/latent-projection step, or a configurable-HVG
diffusion model (BaiZe / Squidiff) retrained on your own full-gene matrix.

### WAVE — the best-supported "just run it" option
β-VAE: `GeneVAE` encodes the **basal** profile, `DrugNN` maps a Morgan fingerprint to a 128-d
embedding, `GeneDrugFusion` predicts the **delta**, added back to the reconstructed baseline. MIT
licence, FigShare weights + preprocessed data. Two hard limits: (i) **dose is not an input** — it was
trained on LINCS 2020 at 10 µM only, and a reviewer explicitly challenged this; (ii) 978 or 3000 genes.
Its own peer-review file contains a damaging reviewer point you should read before trusting the
benchmark numbers: DeepCE and MultiDCP scored ≈0 in several panels and the authors conceded they used
*default* hyperparameters with no tuning. Treat WAVE's SOTA claims as "competitive, not established".

### PerturbNet — the only distributional, dose-aware, structure-driven model
Predicts a *distribution* of cell states via a conditional invertible network, so it gives you
heterogeneity rather than a mean vector — directly useful for your efficacy ranking. ChemicalVAE on
one-hot SMILES pretrained on 250k ZINC molecules; dose and cell type enter the mapping network as
covariates. GPL-3.0 (copyleft — check this against your distribution plans). No author checkpoints;
the HF `cyclopeta/PerturbNet_reproduce` repo is third-party. Independently, the 2026 *Sci Adv*
benchmark (PMC13557097) found PerturbNet had the **best DE-signature recovery** of 13 methods — a
meaningful external validation.

### DeepICER and ExPO — the two that actually take dose seriously
- **DeepICER** is the closest match to your formulation: dose (1/10/20 µM) × duration (24/48 h) × cell
  line, bilinear attention over (compound, cell line), GCN on the molecular graph. **But there is no
  LICENSE file**, so it is legally "all rights reserved" — fine for internal evaluation, not for
  redistribution.
- **ExPO** models the response as a *continuous field* over (dose, time) with a DeepONet, uses
  **Bemis–Murcko scaffold splitting** for the cold-drug test (the strictest chemical OOD protocol I
  saw in this literature), and is MIT. If you want a defensible cold-drug number to quote, this is the
  protocol to copy.

### CMonge — the cheap, honest baseline you should not skip
22.8k parameters, 10–50 MB GPU, ~10.6 s per new condition, MIT, on PyPI. It conditions on **dose**
and on **RDKit structure**, and reports the honest result that structure embeddings only close the gap
to MoA embeddings once you have many conditions. For a 16 GB / 16 h budget this is essentially free,
and it is a strong sanity floor for any heavier model.

### DECANT — right idea, wrong artefact state
The framing (responses are entangled with cell line/dose/time; learn a mechanism-aligned drug
representation and *prove* it is stable across context) is exactly the inductive bias a repurposing
screen wants. But the repo is explicitly a **reproduction skeleton**: no weights, no processed data,
empty placeholder directories. Budget for reimplementation.

### Watch list
- **Perturbed-GE-LDM** would be my #1 pick on function alone (SMILES + dose + time + cell type, latent
  diffusion, pretrained checkpoint) — but it carries a bespoke "Pre-publication Research License" and
  I found no peer-reviewed paper. Verify the licence text against your use case before investing.
- **BaiZe** is architecturally the most complete multi-context design (Δx prediction, drug structure +
  dose + time + species + optional ATAC) but is a 2-star repo with no paper and no licence.
- **X-Pert** (MIT, dose-aware, unified) is the most credible "next gen" competitor, but it leans on
  scGPT, so a 16 GB single-GPU run is doubtful.

---

## 3. What is *not* here, and why it matters

Every model in rows 19–32 fails on the **same** axis: it treats the perturbation as a **learned
per-perturbation latent** (CPA-style) or a **categorical index** (CellCap, scVIDR, LPM) rather than
mapping from chemistry. That is fine for interpolation but mathematically cannot extrapolate to a
compound absent from training, no matter how much data you add. When you audit a new candidate paper,
the single fastest discriminator is therefore: *is the drug input a SMILES/graph/fingerprint that
passes through a fixed encoder, or an `nn.Embedding(n_drugs, d)`?* If the latter, cold-drug
generalisation is impossible and the paper's "unseen drug" split is probably a strawman.

---

## 4. Verification log (what I actually did)

- Europe PMC REST `search` with `TITLE:` / `ABSTRACT:` queries for every named candidate; `resultType=core`
  for DOI/journal/year, `fullTextXML` + regex for `github.com|gitlab.com|zenodo.org|huggingface.co` to
  extract the **paper's own** code-availability URL rather than guessing.
- arXiv Atom API for names not indexed in PubMed.
- GitHub: `api.github.com/repos/{owner}/{repo}` for stars/licence/`pushed_at`; `raw.githubusercontent.com`
  for README and `LICENSE` bytes (licence strings below are from the actual file contents, not badges).
- bioRxiv API for CellCap metadata.
- **No repositories or licences were inferred from memory.** Anything I could not fetch is marked
  `[U]` / "unverified" / "no public code found".

### Licence summary (copy-pasteable)
| Licence | Models |
|---|---|
| **MIT** | WAVE, ExPO, CMonge, DECANT, scPerb, CellPert, scVIDR, CellBox, scCausalVI, X-Pert, PerturbGraph **[U]** |
| **GPL-3.0** (copyleft — check before redistribution) | PerturbNet, W1OT, Dr.VAE |
| **BSD-3-Clause** | CellCap, Squidiff |
| **CC non-commercial** | SAMS-VAE |
| **Bespoke / restrictive** | Perturbed-GE-LDM ("Pre-publication Research License") |
| **No licence file → legally all rights reserved** | DeepICER, MiTCP, CIGER, DeepCE, MultiPert, CRADLE-VAE, BaiZe, AROMA, GPO-VAE (`NOASSERTION`) |
