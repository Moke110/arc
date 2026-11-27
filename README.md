# CarPS

CarPS (Cardiomyocyte Perturbation Prediction System) studies how chemical
perturbations change gene expression and phenotypes in cardiomyocytes.
CarPS-Transcriptome predicts signed gene-expression changes. CarPS-Phenotype
has established a structure-and-dose baseline for phenotype direction and
defines a future comparison using predicted gene-expression changes.
[CarPS-Transcriptome report](cpst/README.md) ·
[CarPS-Phenotype report](cpsp/README.md)

## CarPS-Transcriptome (CPST)

The selected Stage 2 variational autoencoder predicts treatment-induced gene
expression change (Δ) from measured same-plate DMSO control expression and a
chemical representation. The selected `XΔ-XΔ` architecture uses one signed Δ
decoder for both reconstruction and prediction, frozen UniMol features and the
ordered Denoise-DEG 5k gene axis. The `no-dose` prediction head had the highest
mean Δ PCC and lowest mean Δ RMSE in the tenv comparison on held-out drug
parents and doses. Dose still defines measured treatment conditions and
evaluation partitions.
[Selected configuration](cpst/README.md#selected-carps-transcriptome-configuration) ·
[Dose comparison](cpst/10_cpst_dose/README.md#formal-results)

The fixed Step 07 alpha release preserves KPGT and training-fold standardized
log dose. Its parity check reproduced the predictions and validation objective
from all nine sealed Step 06.4 checkpoints. The CPST report traces that
reference release through the later representation, gene-axis and dose
studies. [Alpha replay evidence](cpst/07_cpst_alpha/README.md#checkpoint-and-replay-evidence) ·
[Experiment sequence](cpst/README.md#evidence-in-experiment-order)

## CarPS-Phenotype (CPSP)

Stage 3 froze two Heliste cell-viability corpora from HL-1 cardiomyocytes
measured by CellTiter-Glo, one under duration-matched normoxia and one under
hypoxia–reoxygenation. The verified Morgan baseline uses chiral fingerprints
and training-partition standardized log dose. All 60 formal fits were
replayed across cold-drug, cold-dose and warm-start evaluations. Cold drug
had the lowest mean macro-F1 in both conditions, with zero recall for the
increase class in every cold-drug fit. [Corpus report](cpsp/01_corpus/README.md) ·
[Morgan results and verification](cpsp/02_base_morgan/README.md#formal-results-and-release-evidence)

The prospective predictor would take signed CPST Δ on the selected ordered
gene axis plus assay dose. Each condition–phenotype pair has its own model
with three output probabilities for decrease, within-band and increase.
The paired comparison will reuse the frozen split, labels and metrics on
the same covered held-out measurements. Selecting the CPST prediction
release for Heliste, checking coverage by parent and dose, and integrating
the inputs remain pending. [Prospective interface](cpsp/README.md#prospective-carps-transcriptome-interface) ·
[Comparison protocol](docs/adr/0045-stage3-morgan-direction-baseline.md)

## Data and reproduction

| Source | Contents |
| --- | --- |
| [Perturbation data](pert_data/README.md) | Chemical perturbation transcriptomes and annotations |
| [Disease data](disease_data/README.md) | Myocardial disease transcriptomes used to develop CPST gene axes |
| [Phenotype data](pheno_data/README.md) | Assay source archive, including Heliste |
| [Perturbation models](pert_model/README.md) | Upstream model implementations and assets |
| [Drug representation models](drug_emb_model/README.md) | Upstream chemical representation libraries and assets |

The public tree carries code, reports and tracked provenance. Large source
downloads and generated payloads are materialized in the working checkout.
Each external item records its source and declared local additions in
`author_log.json`, with source checksums in `manifest.sha256`. Use the owning
source guide for acquisition and offline verification.
[Provenance contract](docs/adr/0030-external-data-provenance-and-step-layout.md) ·
[Data inventory](docs/data-inventory.md)

- For repository setup and offline checks, see [Start work](AGENTS.md#start-work)
  and [External data and verification](AGENTS.md#external-data-and-verification).
- CPST source preparation, shared APIs, fit replay and figure verification
  follow the [CPST guide](cpst/AGENTS.md).
- CPSP corpus preparation, Morgan replay and release verification follow the
  [CPSP guide](cpsp/AGENTS.md).

Stage 1, the local predecessor, and its workflow documentation are kept
privately outside the public Git tree. Required public checks run independently
of Stage 1. [Public release boundary](docs/adr/0106-cpst-and-cpsp-public-release-boundaries.md)

See the [glossary](GLOSSARY.md) for domain terminology and the
[decision index](docs/adr/README.md) for the design record.
