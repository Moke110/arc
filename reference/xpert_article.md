# Modelling drug-induced cellular perturbation responses with a biologically informed dual-branch transformer

**Received:** 8 May 2025  
**Accepted:** 2 December 2025  
**Published online:** 26 January 2026  

**Yue Guo¹, Hao Zhang², Haitao Hu², Jialu Wu², Ji Cao¹,4,5,6, Chang-Yu Hsieh®, Bo Yang®¹,4,5,7**

---

## Abstract

Systematic mapping of chemical perturbation responses is revolutionizing polypharmacological drug discovery, yet remains constrained by experimental scalability. Here we introduce XPert, a biologically informed dual-branch transformer model designed to model gene-specific perturbation effects and dose‑time dynamics. The dual‑branch architecture separately encodes pre‑perturbation and post‑perturbation cellular states, allowing the model to disentangle intrinsic transcriptional patterns from regulatory shifts triggered by perturbations. By leveraging context‑aware gene network modelling, XPert overcomes the over‑denoising issues inherent in dominant variational‑autoencoder‑based approaches, achieving 36.7% higher Pearson's correlation coefficient and 78.2% lower mean square error in cold‑cell generalization under single‑dose‑single‑time scenarios. Through extension to multi‑dose‑multi‑time prediction, XPert precisely resolves pharmacodynamic trajectories and uncovers key molecular events underlying drug effects. To address real‑world data scarcity, we apply knowledge transfer from large‑scale preclinical screens to clinical contexts, achieving up to 15.04% improvement in patient‑specific response predictions. XPert provides mechanistic interpretability, as evidenced by the identification of clinically validated resistance biomarkers.

---

## Introduction

Drug discovery is shifting from the traditional 'one‑drug‑one‑target' paradigm to a more complex 'one‑drug‑multiple‑targets' framework, which stems from the realization that drugs often interact with multiple molecular targets and pathways, triggering complex signalling cascades that result in diverse phenotypic outcomes. This shift necessitates understanding genome‑wide perturbation effects to elucidate mechanisms of action (MoA) and optimize therapy. To this end, researchers have increasingly turned to large‑scale perturbational profiles across compounds, which facilitate the identification of high‑affinity targets and context‑specific effectors. However, progress remains limited by the scarcity of high‑quality perturbation data—particularly in clinical settings—and by confounding factors within perturbation profiles that obscure mechanistic insights.

Deep learning approaches have emerged as powerful tools in this domain, enabling the mapping of broader perturbation spaces from limited observations. To achieve this, a primary challenge lies in effectively representing context‑specific pre‑ and post‑perturbation states from noisy high‑throughput sequencing data. Early efforts primarily focused on autoencoder‑based methods. Pioneering works, such as CPA and chemCPA, leveraged the autoencoder's denoising capabilities with an adversarial network to decouple and reconstruct perturbation profiles. Building on this, TranSiGen and PRnet introduced variational autoencoders (VAEs), showing improved generalization across cellular contexts. Although effective at eliminating confounders, these methods risk over‑denoising, which can obscure critical biological information. Another line of research, including DeepCE and CIGER, explored attention‑based architectures to model gene‑drug interactions but overlooked the cellular context, yielding suboptimal results. Although the advent of single‑cell large models has advanced the characterization of cellular context, they struggle with zero‑shot adaptation of the post‑treatment space. This motivates our approach, which separately addresses pre‑ and post‑perturbation states.

Another fundamental challenge is how to translate chemical perturbations into biological perturbation signals. State‑of‑the‑art (SOTA) approaches typically concatenate chemical and cellular features, capturing global cell‑state alterations but failing to resolve gene‑specific responses. This limitation necessitates advanced fusion strategies that integrate prior knowledge (for example, drug‑target interactions, DTIs) to bridge chemical and biological spaces. A further gap is the inadequate modelling of the well‑established dose‑ and time‑dependent nature of drug effects. Previous attempts have relied on simplistic encodings (for example, one‑hot encoding), which is insufficient for modelling nonlinear dose‑response relationships (for example, inverted U‑shaped curves), restricting a full understanding of transcriptional pharmacodynamics.

To address these challenges, we introduce XPert, a transformer‑based solution that predicts drug‑induced transcriptional perturbations by jointly modelling cellular contexts, multiscale drug properties and dose‑time dynamics. XPert features a dual‑branch architecture that concurrently captures intrinsic gene‑gene interactions and extrinsic chemical perturbation effects. This architecture is further augmented by a knowledge‑informed heterogeneous graph (HG) to bridge chemical‑biological spaces and by condition tokens that encode nonlinear dose‑time responses. XPert excels in diverse scenarios, including generalization to unseen drugs and cells and multi‑dose‑multi‑time (mdmt) predictions. It also supports a pretrain‑fine‑tune framework that extends large‑scale preclinical data for clinical predictions. By enabling precise in silico modelling of perturbation dynamics across different scales, XPert establishes a transformative tool for preclinical drug discovery and personalized therapeutic development.

---

## Results

### Overview of XPert

XPert is a transformer‑based model designed to predict drug‑induced transcriptional perturbation effects. Given the unperturbed cell gene expression and perturbation attributes, it simultaneously outputs the post‑perturbation cell expression (x_pert) and the difference between post‑ and pre‑perturbation gene expression (x_deg = x_pert – x_base).

XPert leverages a dual‑branch architecture to encode both pre‑ and post‑perturbation cell states (Fig. 1). The base encoder, built with stacked self‑attention layers, models complex gene‑gene interactions under diverse cellular contexts, whereas the perturb encoder branch uses cross‑attention to capture cell‑drug interactions and condition‑dependent perturbation effects. Each cell is represented as a 'sentence' of gene tokens, along with a `<cls>` token representing the global cell state. Each gene token is initialized with its functional representation and binned expression value, and is then dynamically refined based on regulatory interactions with other genes and the constraints of the perturbation.

XPert also models four key perturbation attributes: drug's chemical properties, biological properties, perturbation time and dosage. All features are tokenized and fed into the model. Chemical features are extracted using UniMol, a powerful three‑dimensional (3D) molecular model. Biological tokens are derived from a knowledge‑informed HG built on the drug's MoA, encompassing three key relationships: DTI, protein‑protein interaction (PPI) and drug‑drug structure similarity (DDS). Given the sparsity of known DTIs, this graph infers potential drug‑gene interactions, informed by two biological intuitions: (1) genes close in the PPI network respond similarly to perturbations and (2) structurally similar drugs often yield comparable effects. Through unsupervised HG pretraining, XPert bridges the chemical and biological spaces of drugs, generating embeddings that reflect their biological effects. Furthermore, to account for the dose‑dependent and time‑dependent biological responses of drugs, XPert introduces condition tokens (for example, dose and time) to capture the nonlinear transcriptional effects to varying conditions.

By leveraging intra‑ and cross‑modal attention mechanisms and enhanced by knowledge graph representations, XPert captures the intricate interplay between drugs and genes, leading to more accurate predictions of gene expression changes in response to perturbations.

### Benchmarking drug perturbation prediction in single‑dose‑single‑time scenario

We benchmarked XPert against existing methods on the L1000 dataset, a major resource for studying transcriptomic responses to perturbations. In this paper, we aim to conduct a comprehensive benchmark of existing methods and XPert based on the L1000 dataset. To ensure a fair comparison against models lacking explicit dose‑time modelling, we first focused on a simpler scenario: the single‑dose‑single‑time (sdst) prediction task. For this, we created the L1000_sdst subset by filtering for the most common perturbation time and dose. We then performed a strict fivefold cross‑validation using three split strategies: (1) warm‑start: random splits to test generalization to unseen cell‑drug pairs; (2) cold‑drug: excluding test drugs from training; and (3) cold‑cell: excluding test cell lines from training.

Our benchmark includes four SOTA models: two VAE‑based models (TranSiGen and PRnet) and two attention‑based models (DeepCE and CIGER). In particular, both DeepCE and CIGER do not account for the cellular context. All baselines solely focus on drug chemical features. Additionally, we implemented two simple MLP baselines (MLP_UniMol and MLP_Morgan) that differ by their input drug features. Finally, three mean‑based baselines (Mean, Mean_cell, Mean_drug) were implemented to assess whether models learned beyond simple population‑level averages.

To facilitate a systematic comparison, we adopted a diverse set of metrics inspired by recent single‑cell perturbation benchmarks covering four categories: error, goodness‑of‑fit, correlation and distribution metrics. Metrics were separately computed for two prediction targets: perturbed profile (x_pert) and gene expression changes (x_deg). Predicting x_deg is more challenging as it requires capturing subtle transcriptional shifts. To assess this, we also included precision metrics measuring the proportion of correctly predicted top up‑ and downregulated genes, providing a multifaceted evaluation of model performance.

XPert consistently outperformed all baselines, particularly in the challenging x_deg prediction task. We attribute this success to its dual‑branch design, which effectively learns pre‑ to post‑treatment state changes, a finding validated by ablation studies. For instance, XPert's Pearson's correlation coefficient (PCC) surpassed the next‑best model, TranSiGen, by 8.2% (warm‑start), 15.9% (cold‑drug) and 36.7% (cold‑cell). In particular, the context‑specific mean baselines (Mean_cell and Mean_drug) were highly competitive, outranking some complex models. Additionally, XPert's performance was stable across different random seeds.

The cold‑cell setting proved the most challenging due to cell‑specific drug responses, with an average performance drop of 121% versus the warm‑start scenario. Although the performance varied significantly across cell lines, XPert showed the lowest variance, indicating stable predictions. This variability correlated partially with cell‑line similarity to the training set, as performance dropped from high‑ to low‑similarity groups (mean square error ranged from 0.24 to 0.66). Despite this, XPert still achieves an average gain of 67.54% over the current SOTA, TranSiGen, in the cold‑cell setting, demonstrating a substantial advance in generalization.

Moreover, our results show an unreported limitation of VAE‑based models: a lack of robustness in blind tests relative to attention‑based approaches. For example, the leading VAE model, TranSiGen, performed well in warm‑start tests but its performance deteriorated in cold‑cell settings, scoring negative R² values despite good correlation, suggesting a failure to adapt to unseen cellular contexts. We attribute this failure to two intrinsic VAE properties: (1) the Kullback‑Leibler divergence regularizer forces information compression that can lead to over‑denoising, erasing critical cellular context features; (2) VAEs are constrained by their training data, leading to low‑fidelity outputs when encountering out‑of‑distribution samples like unseen cell lines.

Plotting the predicted expression changes of the top‑ten highly variable genes (HVGs) in cold settings visually confirms these findings. XPert most accurately captured the mean and range of gene expression changes and was the only model to predict correct trends for key genes like AARS and GRN. By contrast, the VAE‑based TranSiGen captured the distributional shape but failed on the magnitude of the effect. These validated the strong advantage of XPert in terms of generalization ability.

### Knowledge‑informed XPert exhibits superior generalization and interpretability

To explore XPert's learned latent features and mechanisms behind its performance, we further analysed its handling of batch effects—a key challenge that hinders generalization in high‑throughput sequencing datasets. Unlike prior VAE‑based models that rely on denoising, XPert explicitly distinguishes true biological signals from noise in a supervised manner, enhancing cell‑specific representation and overall generalization.

We applied uniform manifold approximation and projection (UMAP) to project the raw post‑treatment expression and the `<cls>` token embeddings obtained from XPert in the test dataset. Compared with the raw expression, XPert partially mitigates plate‑related noise, aggregating subclusters of the same cell line (for example, HCCS15 and HAIE) more cohesively. Quantitative scIB benchmark confirmed its strong biological conservation ability. These highlight that XPert captures intrinsic cell identity features, concurrently preserving the perturbations‑induced biological differences. Furthermore, the model incorporates biologically relevant gene embeddings as prior knowledge, which guides it to focus on gene interrelations over sequencing noise.

XPert further benefits from incorporating drug‑gene interaction prior knowledge. Although structurally similar drugs may imply similar properties, biological activity often does not directly correlate with chemical structure. Drugs with the same MoAs are dispersed in the chemical space, indicating a limitation of solely relying on chemical features. XPert overcomes this by using a pretrained‑knowledge HG to create biologically coherent drug representations, where drugs with the same MoA naturally cluster together. Our ablation studies validate that this prior improves predictive performance in the cold‑drug scenario.

Furthermore, what XPert learns is inherently interpretable due to its reliance on attention mechanisms, which explicitly reveal intramolecular or intracellular interactions. Analysis of atom‑level attention for several widely used clinical drugs shows that the model learns chemically meaningful local structures that align with known structure‑activity relationships (SAR). For instance, with epidermal growth factor receptor (EGFR) inhibitors like gefitinib and erlotinib, XPert highlights the quinazoline ring core and its key hydrogen‑bonding N1 and N3 atoms, which are crucial for EGFR binding. For histone deacetylase (HDAC) inhibitors, the zinc‑binding group receives higher attention, as it chelates the catalytic zinc ion at the active site.

By incorporating multidimensional biological prior knowledge, XPert enhances its ability to capture the intricate biological mechanisms driving drug perturbations, offering a more comprehensive and interpretable model for drug‑induced transcriptional responses.

### XPert supports robust transcriptional responses prediction in the mdmt scenario

Understanding dose‑ and time‑dependent responses is fundamental to pharmacodynamic research. Recent advances in experiment‑driven mdmt perturbation studies enable detailed profiling of drug‑induced cellular dynamics, providing critical insights into temporal molecular drivers and potential off‑target effects at high doses. The complexity of such data presents a rigorous test for predictive models. We therefore benchmarked XPert against established baseline methods in this realistic scenario.

For the mdmt scenario, we used the L1000_mdmt subset, which contains drug‑cell pairs with varied dose‑time points, including 40 cell lines and 1,977 drugs. As L1000 includes many pharmacologically equivalent doses from minor experimental variations, we avoided coarse methods like one‑hot encoding. Instead, we aggregated similar doses into ten discrete ranges and encoded them as conditional tokens, applying a parallel strategy to time attributes. By modelling the interplay between these condition tokens and gene networks, XPert's context‑aware framework captures complex dose‑ and time‑dependent response patterns.

We used the three partition strategies—warm‑start, cold‑cell and cold‑drug—to deploy comparative experiments. None of the previous models simultaneously accounted for both dose and time attributes. Some attempts were made by models such as DeepCE and CIGER to encode perturbation dose using one‑hot encoding, whereas PRnet utilized logarithmic doses as the weight of drug features. To adapt baseline models for the mdmt scenario, we concatenated one‑hot encoded dose and time features with their standard inputs.

The results confirmed XPert's dominant performance, as it ranked first in most metrics, followed by TranSiGen and CIGER. For the x_deg task, XPert's PCC improved on the next‑best model by 8.34% (warm‑start), 5.85% (cold‑drug) and 30.54% (cold‑cell). In particular, in the cold‑cell scenario, only XPert and DeepCE avoided negative R² values. A similar phenomenon was observed—although TranSiGen excelled in capturing correlations, it was readily surpassed by XPert and other attention‑based models in terms of fitting ability, error and distribution metrics. CIGER's sharp performance drop in the cold‑cell setting highlights the inadequacy of its simple one‑hot encoding for cell lines. This observation further underscores XPert's notable advance in modelling cell‑specific responses.

To assess whether XPert captures subtle gene expression changes, we conducted a case study on vorinostat, a pan‑HDAC inhibitor with well‑documented dose‑ and time‑dependent biological effects. First, with a fixed 24‑h time point, we presented the transcriptional impact of vorinostat at different doses across ten cell lines. Increasing the dose generally leads to stronger effects on genes. A principal component analysis (PCA) confirmed a clear dose‑response gradient along the first principal component (PC1) that strongly correlated with increasing dose. We also observed that changes in dosage could reverse transcriptional effects. For instance, increasing the dose shifted genes like NRIP1 and ELOVL6 from upregulation to downregulation. Moreover, cell‑type‑specific expression effects were noted, where the drug had opposite effects on the same gene across different cell lines. Crucially, XPert accurately captured these nuanced patterns, consistent with experimental measurements, demonstrating its ability to model complex dose‑response relationships in a cell‑specific context.

Next, we examined x_deg changes in drug biomarkers under various doses using a set of known pharmacodynamic genes for HDAC inhibitors, covering several critical cellular processes, such as proliferation, apoptosis, metastasis and immunogenicity. Transcriptional responses were not uniform across doses; different genes responded at different concentrations. For instance, vorinostat downregulates TP53 and alters cell‑cycle genes at lower doses, preceding changes in its direct targets, suggesting potential combinatorial therapies at lower doses.

To jointly investigate the role of treatment time, we applied radial basis function interpolation to fit dose‑ and time‑dependent effects of specific genes in the A549 cell line. Genes like HDAC6 and TP53 have both dose‑ and time‑dependent responses, whereas others like NRIP1 are mainly dose‑dependent. These results underscore the importance of jointly modelling dose and time to capture transcriptional perturbation dynamics and elucidate drug mechanisms.

### Few‑shot learning enhances prediction in unseen dose‑time condition

A practical challenge in profiling chemical perturbation responses is the high cost associated with measuring multiple time points and doses, resulting in datasets that include measurements taken at only a single dose or time point. Our analysis of the L1000 dataset indicates that only 6.2% of cell‑drug pairs include mdmt measurements, with most containing only a single dose or time point. To address this issue, we propose leveraging large‑scale mdmt datasets like L1000 for pretraining, followed by fine‑tuning with limited target data of specific cell‑drug pairs, which can yield high‑precision predictions for unmeasured dose‑time conditions. This is based on the idea that transferability between dose and time may be easier than between different drugs and cell contexts, ultimately aiding in the construction of dynamic drug perturbation maps and reducing experimental burdens.

As a proof of concept, we split the complete L1000 dataset into two parts: the L1000_mdmt subset was used as the fine‑tuning dataset, whereas the remaining served as the pretraining dataset. Using stratified sampling, we generated five random splits, where for each cell‑drug pair, one data point was assigned to the test set and a proportion of the remaining data was used for fine‑tuning.

We next compared training from scratch with pretraining‑fine‑tuning to assess the gains from pretraining. Under both experimental settings, XPert consistently demonstrated optimal performance across all metrics. For x_deg prediction, XPert gained 5.64%–12.45% in PCC under various fine‑tuning ratios, with improvements diminishing as the ratio increased. The power of this strategy was most evident with one‑shot fine‑tuning, which substantially improved performance over the zero‑shot setting for most models. For top models like XPert and TranSiGen, one‑shot fine‑tuning matched or surpassed training from scratch on 80% of the dataset. This substantiates our hypothesis that the complete dose‑time response landscape for novel cell‑drug pairs can be accurately inferred from a minimal set of experimental measurements. However, models like DeepCE and CIGER gained little from pretraining, highlighting the necessity for rational model design to maximize the benefits of few‑shot learning.

### XPert bridges preclinical datasets to clinical prediction

Given the challenges associated with obtaining clinical perturbation data, we next investigated transferring knowledge from the large‑scale L1000 preclinical dataset to smaller, high‑fidelity clinical datasets via our pretrain‑fine‑tune framework. We hypothesized that despite its technical noise, L1000 could serve as a valuable low‑fidelity pretraining resource to bolster model performance where clinical data are scarce, offering a promising strategy to bridge the preclinical‑to‑clinical gap.

To test this paradigm, we turned to a clinical dataset—CDS_DB—which includes paired pre‑ and post‑treatment clinical transcriptomic data from cancer patients. Similar to the preclinical setting, we evaluated three different partitioning strategies: unseen‑patient, unseen‑drug and unseen‑cancer. Given the imbalanced cancer type distribution in CDS_DB, we specifically focused on two predominant cancer types: breast cancer and leukaemia.

Surprisingly, despite a notable domain shift between preclinical and clinical data, pretraining consistently enhanced prediction for unseen patients. Specifically, XPert achieved performance gains of 2.51% for pan cancer, 15.04% for breast cancer and 12.58% for leukaemia. The limited gain in pan cancer probably reflects that most cancer types had very few samples (<20) for fine‑tuning, constraining overall performance.

Both XPert and XPert (pretrain) accurately predicted the distribution of x_deg in unseen patients; however, XPert (pretrain) demonstrated more precise capture of extreme values, exhibiting lower error on genes with large expression changes, as confirmed by stratified analysis. This demonstrates that deep learning can learn transferable representations from preclinical data that are effectively refined using only a few clinical profiles. This strategy is, therefore, promising for developing specialized models for specific cancer types.

Moreover, XPert is the only model benefiting from pretraining in both the unseen‑drug and unseen‑cancer settings. We attribute this advantage to XPert's mechanistic foundation in learning fine‑grained drug‑gene interaction patterns, which enables seamless transfer of pretrained pharmacological knowledge across the preclinical‑to‑clinical domain by focusing on conserved interaction mechanisms rather than context‑specific patterns.

Further, we explored the link between drug‑induced transcriptomic changes and clinical responses. Our focus was on a subset of CDS_DB, specifically GSE20181, which provides records of patient responses to letrozole treatment. Responders exhibited a stronger transcriptomic response than non‑responders, characterized by a more pronounced long‑tail distribution of x_deg and a greater number of enriched HVGs. This motivated exploration of patients' pretreatment states to identify key drivers of drug resistance.

To explore this, we conducted additional analyses using gene‑level attention scores captured by the base encoder of XPert. Only a subset of genes showed notable intersample variability and attention patterns remained stable across folds. We next performed a differential attention analysis between two response groups, comparing it with conventional differential expression analysis. FGFR2 was enriched in both analyses, consistent with its reported role in enhancing PI3K/AKT pathway and promoting antioestrogen resistance in breast cancer. More importantly, our attention‑based analysis uniquely identified other key resistance biomarkers, such as TIAM1, RCPP, HK1 and CDKN1B, that were invisible to expression‑level analysis. These results underscore the power of attention‑based methods to reveal latent gene‑phenotype associations beyond mere expression changes, thereby providing a new layer of insight into drug resistance mechanisms.

---

## Discussion

In this study, we introduce XPert, a knowledge‑guided, dual‑branch attention framework for predicting drug‑induced transcriptional responses. Through a comprehensive evaluation on both x_pert and x_deg prediction tasks across multidimensional metrics, our results highlight the exceptional capability of attention‑based frameworks in context‑aware cellular modelling. When tasked with inferring responses in unseen cellular states, XPert outperformed the next‑best model by an average of 67.54% across all metrics in the sdst scenario. Moreover, our analysis reveals a previously underappreciated limitation of dominant VAE‑based approaches—excessive overcorrection that obscures cellular context in blind‑test scenarios—resulting in substantial deficits in both error metrics and expression‑change distribution fidelity.

Another pivotal contribution of this work lies in addressing dose‑time dynamics in drug‑induced transcriptional effects. We propose a universal encoding method for perturbation attributes (e.g., dose and time), enabling interpretable modelling of nonlinear pharmacodynamic relationships. XPert represents the most effective framework currently available for mdmt scenarios. We demonstrated this by generating gene‑specific 3D dose‑time response maps (e.g., for vorinostat) that reveal dynamic gene network reorganization induced by chemical perturbations. Furthermore, our proof‑of‑concept experiments establish that with few‑shot learning strategies, deep learning algorithms can assist in constructing more comprehensive dynamic maps by interpolating unmeasured dose‑time conditions. This approach promises to substantially reduce experimental burden and accelerate construction of large‑scale perturbation omics landscapes.

Beyond that, we extend the XPert framework by applying transfer learning to overcome clinical data scarcity and translational roadblocks in perturbation studies. By modelling conserved drug‑gene interactions, it enables reliable knowledge transfer from larger‑scale preclinical screens to patient transcriptomes, supporting personalized response prediction and biomarker identification. To ensure completeness, we extended this strategy to PANACEA, an independent preclinical dataset. Despite discrepancies in measured perturbation signals that can induce negative transfer, XPert demonstrated superior robustness over other models.

Although XPert demonstrates strong performance, we identify key avenues for future development. One limitation is its computational cost, which could be mitigated by memory‑efficient training strategies (e.g., DeepSpeed) and scalable architectures like Hyena designed for long‑range dependencies. Biologically, extending the framework's scope is the primary goal. This includes transitioning from bulk‑level to single‑cell‑level predictions as large‑scale datasets become available, and broadening the model beyond small‑molecule transcriptional effects to encompass biologics, genetic perturbations and multiomics integration, contingent on future data availability.

In summary, XPert represents a substantial step forward in modelling drug‑induced perturbation effects through an interpretable and generalizable deep learning framework. With further development, XPert holds substantial promise as a core component of next‑generation in silico drug discovery pipelines and precision medicine platforms.

---

## Methods

### Dataset preprocessing

To systematically evaluate performance of XPert and SOTA models in drug perturbation prediction, we utilized three benchmark datasets, including two preclinical datasets—LINCS L1000 (referred to as L1000) and PANACEA—as well as one clinical dataset, namely, the cancer‑drug‑induced gene expression signature database (CDS‑DB).

**LINCS L1000 dataset.** The L1000 dataset contains gene expression profiles resulting from various drug treatments across different cell lines. We extracted drug‑induced perturbation and control samples from L1000 level‑3 data. The L1000 platform measures mRNA transcript abundance of 978 'landmark' genes. Data cleaning removed low‑quality data: perturbations with missing/ambiguous information, low‑frequency time points (retained 3 h, 6 h, 24 h), and profiles that did not pass quality control. Each expression profile was matched with a randomly selected DMSO control sample from the same plate to create paired pre‑/post‑treatment profiles. Replicate‑collapsed z‑score vectors were computed. Several subsets were created: L1000_full, L1000_sdst (10 µM, 24 h), L1000_mdmt (multiple dose‑time points per cell‑drug pair), L1000_mdmt_pretrain (L1000_full excluding profiles in L1000_mdmt). Due to thousands of perturbation doses, doses were grouped into ten discrete intervals.

**PANACEA.** A resource from Columbia Cancer Target Discovery and Development Center, including dose‑response and RNA‑seq profiles for 25 cell lines exposed to ~400 clinical oncology drugs. Perturbations at IC20 for 24 h. Processed to log2(TPM+1) and filtered to 978 landmark genes. Doses mapped to ten predefined ranges.

**CDS‑DB.** Contains patient‑derived paired pre‑ and post‑treatment clinical transcriptomic data across 78 treatment‑specific datasets, 85 therapeutic regimens, 39 cancer subtypes, 3,628 patient samples. Retained only microarray data (larger sample size), excluded combination therapies or non‑chemical drugs. Final dataset: 613 paired profiles, 14 cancer subtypes, 14 drugs, restricted to 978 landmark genes. Standardized administration dosage and treatment time into unified intervals.

### Transcript profile embedding

Inspired by transformer architectures in single‑cell large language models, we encode gene expression profiles for pre‑perturbation cells. Each cell is a 'sentence' of genes plus a `<cls>` token capturing global state. Input representation: tensor X ∈ R^{N×(M+1)×d}, where N = number of cells, M = number of genes, d = embedding dimension.

**Input gene embedding:** For gene j, two components: (1) gene token (g_j) using biologically meaningful gene embeddings from CellLM (GraphMAE on PPI network); (2) binned expression value (e_j) using value binning technique (scGPT). Binning converts expression counts into relative values with shared bin edges across dataset. PyTorch embedding layers map to dimension d.

**Cell embedding:** Special `<cls>` token initialized with Gaussian distribution, appended to beginning of gene token sequence.

Final input embedding for cell i: C_i = [C_cls, G_1, G_2, ..., G_M].

### Drug tokenization

**Chemical tokens:** UniMol (universal 3D molecular pretraining framework) used to derive chemical tokens. For each drug, SMILES → canonical SMILES (RDKit), extract atom types, pretrained model outputs mol token and atom tokens (512‑d), projected to d dimensions. Alternative features: Morgan Fingerprints, KPGT fingerprints (evaluation in Supplementary Table 26).

**Biological token based on prior‑knowledge HG:** Construct heterogeneous knowledge graph with three relationships: DTI (drug‑target interactions), DDS (drug‑drug structure similarity, Tanimoto similarity >0.5), PPI (STRING database, score >700). Drug nodes initialized with UniMol `<mol>` token embeddings, protein nodes with PPI‑derived gene embeddings. Heterogeneous graph neural network (HeteroConv layers with SAGEConv) under unsupervised contrastive learning (InfoNCE loss) to learn latent relationships. Outputs d‑dimensional biological token vectors for drugs.

**Condition tokens:** Encode dose and time by discretizing raw values into predefined ranges, enabling cross‑dataset normalization.

**Integration of tokens:** For drug j, D_j = [<ConditionTokens>, <BiologicalTokens>, <ChemicalTokens>]. Learnable positional embeddings added.

### XPert architecture overview

**Base encoder branch:** Captures unperturbed cell state using stacked self‑attention layers. Given initial representation C_base ∈ R^{(k+1)×d}, applies n layers of self‑attention: C_0_base = C_base, C_i_base = SelfAttention(C_{i-1}_base). Final output C_n_base.

**Pert encoder branch:** Integrates drug molecular features with cellular context via cascaded cross‑attention and self‑attention. Cross‑attention: cell representation as query, tokenized drug representation as key/value. After m layers: C_n_pert = Pert_Encoder(C_base, D).

**Multi‑objective learning:** Three tasks jointly optimized:
1. **Perturbation gene expression prediction (x_pert):** MLP on C_n_pert → ŷ_pert; loss L_pert = α × MSE(x_pert, ŷ_pert).
2. **Gene expression delta prediction (x_deg):** MLP on (C_n_pert – C_n_base) → ŷ_deg; loss L_deg = β × MSE(x_deg, ŷ_deg) + γ × (1 – PCC(x_deg, ŷ_deg)).
3. **Cell‑type classification:** Auxiliary task using `<cls>` token representations, cross‑entropy loss L_cls = δ × CrossEntropy(y_true, y_pred).

Total loss: L_total = L_pert + L_deg + L_cls.

**Training and testing:** Implemented in PyTorch, Adam optimizer (lr=4e‑3, weight decay=1e‑3), LambdaLR scheduler (reduce by 0.5 after warm‑up epochs), early stopping (50 epochs). Flash attention for efficiency. Hyperparameter search (Supplementary Table 11). Fivefold cross‑validation with four split strategies: warm‑start, cold‑drug, cold‑cell, cold‑dose‑time.

**Pretraining and fine‑tuning:** Pretrain on large‑scale dataset (L1000_full or L1000_mdmt_pretrain), then full‑parameter fine‑tuning on downstream datasets.

### Mean baseline models

- **Global mean (Mean):** average expression profile across all training data.
- **Cell‑specific mean (Mean_cell):** average of training samples from same cell line.
- **Drug‑specific mean (Mean_drug):** average of training samples treated with same drug.

### Evaluation metrics

Four categories (ten metrics): error (MSE, RMSE, MAE), goodness‑of‑fit (R²), correlation (PCC, Spearman's ρ), distributional similarity (Wasserstein distance, MMD). Additionally, precision metrics for x_deg: Positive Precision(@K) and Negative Precision(@K) measuring overlap of top‑K up‑/downregulated genes with ground truth.

### UMAP and t‑SNE visualizations

PCA to reduce to 40 dimensions, then UMAP or t‑SNE (k=15 neighbours for UMAP).

### Statistics and reproducibility

Paired t‑test for model comparisons; two‑sample t‑test for differential expression. Significance levels: **P ≤ 0.001; 0.001 < **P ≤ 0.01; 0.01 < *P ≤ 0.05; n.s., P > 0.05.

---

## Data availability

L1000 level‑3 data from https://clue.io/data/CMap2020#LINCS2020. PANACEA from GEO GSE186341. CDS‑DB from http://cdsdb.ncpsb.org.cn/. Processed datasets at figshare (https://doi.org/10.6084/m9.figshare.28955141) and Zenodo (https://doi.org/10.5281/zenodo.15357711). Minimal dataset on GitHub (https://github.com/GSanShui/XPert).

## Code availability

Source code at Zenodo (https://doi.org/10.5281/zenodo.15357711), figshare, and GitHub (https://github.com/GSanShui/XPert).

---