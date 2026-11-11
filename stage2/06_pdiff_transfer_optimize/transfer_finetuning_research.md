# PerturbDiff Transfer fine-tuning options

This targeted search used Biomni's arXiv and PubMed tools to find literature,
then checked candidate methods against the original papers and the local
training and scoring code. The options below are hypotheses for Tahoe-to-tenv
transfer.

## Local reference

The Step 06 cold-drug fold 0 screen (seed 17) starts from the globally cold
Tahoe EMA at step 8,000. Full-model tenv fine-tuning uses a fresh AdamW
optimizer, scheduler and EMA, stops at step 7,500 and selects the step 6,500
EMA by validation diffusion loss. Each stage has a 30,000-step cap and a
base learning rate of `2e-4`. The 200-step linear warmup is followed by a
20,000-step plateau on a 200,000-step scheduler. Both runs stop on the
plateau, before cosine decay affects their learning rates.
[Checkpoint index](../05_pdiff_transfer/execution/formal_record/checkpoints.json) ·
[Transfer rules](../../docs/adr/0047-perturbdiff-tahoe-transfer-to-tenv.md#training-and-checkpoint-reuse) ·
[Training configuration](../03_modes_tenv/execution/models/perturbdiff/v2.config.json) ·
[Scheduler call](../04_pdiff_hyperparameter/execution/script/native_training.py#L191) ·
[Scheduler implementation](https://github.com/DeepGraphLearning/PerturbDiff/blob/31be89cb8cce5dddef0782bac29d5c748ce55931/src/common/utils.py#L85-L98)

Historical retrieval scores were removed when Stage 2 adopted full-CM
significance matching. The research recommendations below describe training
mechanisms; current comparisons use the shared three-metric tenv evaluator.

## Options

| Path | What the paper tested | How it would differ from local Transfer |
| --- | --- | --- |
| **1. Tune the conditioning pathway first** | scDCA adds drug-conditioned adapters to frozen scGPT blocks for transcriptional response prediction. PertAdapt trains a condition-sensitive adapter over frozen scFoundation or AIDO.Cell for genetic perturbations. LP-FT warms up a task head before full fine-tuning in image classification under distribution shift. [scDCA, MLGenX 2025](https://openreview.net/pdf/08e45b41bba2972564303604f20f42041fa79f12.pdf) · [PertAdapt, Bioinformatics 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC13341120/) · [LP-FT, ICLR 2022](https://openreview.net/forum?id=UYneFzXSJWh) | Freeze the Tahoe denoiser while adapting the positive-dose embeddings and, in a separate expansion, the dose-input part of the covariate projection. Evaluate staged full unfreezing as a separate contrast. Each arm uses the **same** dose initialization as its full-model comparator. This tests adaptation through the existing conditioning pathway under the diffusion loss. |
| **2. Use separate learning rates for transferred and dose parameters** | Sun et al. evaluated decreasing learning rates toward lower layers of BERT for text classification. [BERT fine-tuning, §5.3.4](https://arxiv.org/pdf/1905.05583) | Keep the full model trainable. An initial contrast could use `2e-4` for positive-dose parameters and `1e-4` for the transferred denoiser, with other groups assigned in advance. Test with the same dose initialization and schedule shape as its control. Step 06 lowered every trainable group's base rate, so it does not test this contrast. A block-depth decay is a separate variant if the simpler grouping is informative. |
| **3. Use parameter-efficient denoiser adaptation** | DiffFit tunes selected biases, normalization and class conditioning plus learned scales initialized to 1 in image DiT. Its published Table 1 has lower average FID for DiffFit than full fine-tuning, while the two tested attention-LoRA configurations performed worse. LoRA itself learns a low-rank update to a frozen projection, starting from a zero update. [DiffFit, ICCV 2023](https://openaccess.thecvf.com/content/ICCV2023/papers/Xie_DiffFit_Unlocking_Transferability_of_Large_Diffusion_Models_via_Simple_Parameter-efficient_ICCV_2023_paper.pdf) · [LoRA, ICLR 2022](https://openreview.net/forum?id=nZeVKeeFYf9) | If conditioning-only has too little capacity, compare it with one predeclared addition such as unit-initialized residual scales. Several local block LayerNorms have no affine parameters, so inventory the actual trainable modules for a DiffFit-style test. Low-rank updates on the treated and control adaLN projections are a separate option. Report trainable parameters, VRAM and wall time alongside AUPRC. |
| **4. Replace the categorical dose bridge** | CPA uses a nonlinear dose-dependent weight and chemCPA uses a chemical-feature-conditioned dosage scaler. These are perturbation autoencoders. [CPA, Mol Syst Biol 2023](https://link.springer.com/article/10.15252/msb.202211517) · [chemCPA, NeurIPS 2022](https://proceedings.neurips.cc/paper_files/paper/2022/file/aa933b5abc1be30baece1d230ec575a7-Paper-Conference.pdf) | Compare a small shared continuous log-dose encoder with an explicit zero-dose control branch against the categorical bridge and target dose reset. Fit a matched tenv-only model with the same encoder to separate representation changes from transfer gains. A shared source-to-target encoder requires a matched Tahoe refit. Three source levels constrain the learned curve, and tenv 10 µM is outside Tahoe's 0.05–5 µM range. Report that dose separately. |

PerturbDiff's own paper evaluates a different source objective called
marginal pretraining. The authors learn context-conditioned cell distributions before
introducing perturbation supervision and report data-efficiency gains on
Replogle and downsampled PBMC. Local Tahoe Transfer instead starts from a
drug-conditioned Tahoe fit. An unperturbed or marginal source fit is a
separate source-objective experiment if the same-checkpoint adaptations
underperform. The published results do not estimate gains for Tahoe-to-tenv
transfer or AUPRC. [PerturbDiff, §4.6](https://arxiv.org/html/2602.19685v1#S4.SS6) ·
[§5.3](https://arxiv.org/html/2602.19685v1#S5.SS3)

Custom Diffusion adapts a diffusion U-Net through its text-conditioning
key/value projections. Here the denoiser uses joint attention over treated
and control representations. Drug and dose enter through the covariate
projection and adaLN, making those projections the corresponding local
adaptation sites. [Custom Diffusion, CVPR 2023](https://openaccess.thecvf.com/content/CVPR2023/papers/Kumari_Multi-Concept_Customization_of_Text-to-Image_Diffusion_CVPR_2023_paper.pdf) ·
[Pinned covariate encoder](https://github.com/DeepGraphLearning/PerturbDiff/blob/31be89cb8cce5dddef0782bac29d5c748ce55931/src/models/covariate_encoding.py#L278-L307) ·
[Pinned MM-DiT block](https://github.com/DeepGraphLearning/PerturbDiff/blob/31be89cb8cce5dddef0782bac29d5c748ce55931/src/models/cross_dit/cross_dit_blocks.py#L28-L75)

PerturbLDM is a direct Tahoe precedent for representing the same three
positive concentrations as a scalar log dose in a latent diffusion model.
Its cross-dataset analysis calibrates precomputed Tahoe responses. It
supplies a representation precedent, with a different denoiser and
adaptation procedure. [PerturbLDM preprint](https://www.biorxiv.org/content/10.64898/2026.08.07.743610v1.full)

## Comparison plan

Start with a new paired conditioning-only versus full-model contrast using
the same Tahoe EMA and the same tenv initialization for each pair. Compare
the dose-reset/full-model arm with original Transfer to measure the
initialization effect, and compare the dose-reset/conditioning-only arm
with the dose-reset/full-model arm to measure the effect of constraining
updates. Use an explicit parameter slice or positive-dose residual so the
0 µM row remains fixed in the conditioning-only arm, including under AdamW
weight decay. Add the mapped-dose/conditioning-only cell if the interaction
between the two factors is worth estimating. A staged unfreeze or two-rate
experiment is the next same-checkpoint contrast. Reserve PEFT expansion and
a continuous-dose source refit for separate questions.

Fix the tenv split, drug–dose groups, source exclusion, controls, Morgan
fingerprints, training objective and inference protocol. Use the same
planned update budget across paired arms. Report both a predeclared
reachable matched step and each arm's diffusion-loss-selected EMA, plus
paired overall and per-dose AUPRC against original Transfer and CM. Keep
tuning on validation and predeclare the confirmation rule for any new
screen. The completed Step 06 decision remains governed by ADR 0048.
[Baseline protocol](README.md#validation-baseline) ·
[ADR 0048](../../docs/adr/0048-perturbdiff-transfer-diagnosis.md)
