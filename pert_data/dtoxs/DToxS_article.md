# Transcriptomic profiling of human cardiac cells predicts protein kinase inhibitor-associated cardiotoxicity

**Authors:** J.G. Coen van Hasselt^1,2,7, Rayees Rahman^1,7, Jens Hansen^1, Alan Stern^1, Jaehee V. Shim^1, Yuguang Xiong^1, Amanda Pickard^1, Gomathi Jayaraman^1, Bin Hu^1, Milind Mahajan^3, James M. Gallo^1,4, Joseph Goldfarb^1, Eric A. Sobie^1, Marc R. Birtwistle^1,5, Avner Schlessinger^1,8, Evren U. Azeloglu^1,6,8 & Ravi Iyengar^1,8

**Journal:** *Nature Communications* (2020) 11:4809 | https://doi.org/10.1038/s41467-020-18396-7

**Affiliations:**

1. Department of Pharmacological Sciences and Systems Biology Center New York, Icahn School of Medicine at Mount Sinai, New York, NY, USA.
2. Division of Systems Biomedicine and Pharmacology, Leiden Academic Centre for Drug Research, Leiden University, Leiden, Netherlands.
3. Department of Genetics and Genomic Sciences, and Icahn Institute for Genomic Sciences and Multiscale Biology, Icahn School of Medicine at Mount Sinai, New York, NY, USA.
4. Department of Pharmaceutical Sciences, School of Pharmacy and Pharmaceutical Sciences, University at Buffalo, Buffalo, NY, USA.
5. Department of Chemical and Biomolecular Engineering, Clemson University, Clemson, SC, USA.
6. Department of Medicine, Division of Nephrology, Icahn School of Medicine at Mount Sinai, New York, NY, USA.
7. These authors contributed equally: J.G. Coen van Hasselt, Rayees Rahman.
8. These authors jointly supervised this work: Avner Schlessinger, Evren U. Azeloglu, Ravi Iyengar.

## Abstract

Kinase inhibitors (KIs) represent an important class of anti-cancer drugs. Although cardiotoxicity is a serious adverse event associated with several KIs, the reasons remain poorly understood, and its prediction remains challenging. We obtain transcriptional profiles of human heart-derived primary cardiomyocyte like cell lines treated with a panel of 26 FDA-approved KIs and classify their effects on subcellular pathways and processes. Individual cardiotoxicity patient reports for these KIs, obtained from the FDA Adverse Event Reporting System, are used to compute relative risk scores. These are then combined with the cell line-derived transcriptomic datasets through elastic net regression analysis to identify a gene signature that can predict risk of cardiotoxicity. We also identify relationships between cardiotoxicity risk and structural/binding profiles of individual KIs. We conclude that acute transcriptomic changes in cell-based assays combined with drug substructures are predictive of KI-induced cardiotoxicity risk, and that they can be informative for future drug discovery.

## Introduction

Protein kinase inhibitors (KIs) are an important class of therapeutics used for the treatment of various forms of cancer and other diseases. There are currently more than 48 KIs approved for clinical use by the U.S. Food and Drug Administration (FDA) and other regulatory agencies, and more than 250 KIs are undergoing clinical trials or are in development. The clinical effectiveness of KIs as cancer drugs has led to a broad effort to develop drugs that are more efficacious and have reduced the propensity for adverse events. Cardiotoxicity (CT) is a clinically important adverse event associated with several KIs. KI-associated CT manifests as loss of cardiomyocyte function, which can lead to heart failure. Given the extensive therapeutic potential of KIs, approaches to identify and subsequently mitigate the risk for CT during early development of novel KIs and during clinical administration are urgently required.

We do not yet sufficiently understand the mechanisms underlying KI-associated CT. The human kinome consists of more than 500 protein kinases. Given that many KIs exhibit multitarget pharmacology, inhibition of multiple protein kinases in cardiomyocytes may lead to adverse drug effects such as CT. For individual KIs, pathways involved in mitochondrial function, endoplasmic reticulum stress response, and AMPK inhibition, have been shown to be associated with KI-induced CT. Overall, however, the general mechanisms of KI-induced CT are still poorly understood.

Obtaining quantitative clinical risk scores for KI-associated CT is also challenging, as the risk for KI-associated CT has not been systematically studied. The FDA adverse event report system (FAERS) database has been previously applied to quantify the risk of adverse drug reactions (ADRs). The FAERS database contains over 9 million individual drug-associated adverse-event reports reported by industry and physicians. Through statistical analyses of the FAERS database, relatively unbiased estimates for the relative risk of specific ADRs can be computed. Such risk scores are clinically relevant as they are based on real-life patient population, and they are not solely based on selected patient cohorts. We previously used such analyses of the FAERS database in combination with systems' pharmacology-based approaches to obtain mechanistic insights into adverse-event mechanisms.

In the current study, generated as part of the NIH-funded Library of Integrated Network-Based Cellular Signatures (LINCS) Drug Toxicity Signature Generation Center (DToxS), we take a top-down global approach to determine if a comprehensive profiling of gene expression changes in human cardiomyocytes can provide insight into pathways associated with KI-induced CT, and to potentially predict the risk of CT. The rationale for this approach is based on the central assumption that CT largely originates from cardiomyocytes where one or more protein kinases contribute to the pathophysiology. Since progression to heart failure takes several months to manifest, it is not immediately obvious if gene expression changes measured after drug treatment for a few days would have predictive value. Thus, a second important assumption is that early changes in gene expression upon drug treatment of cardiomyocytes are indicative of later physiological events. We test the validity of our assumptions by experimentally obtaining gene-expression patterns for the different KIs, and if these patterns could be selectively associated with the clinical risk of CT for each KI, thereby providing gene-expression signatures for KI-associated CT.

We report the generation of transcriptomic profiles from four human primary cardiomyocyte-like cell lines. These profiles are generated using 23 KIs that were FDA-approved and used extensively at the time of experimental design, such that an adequate number of clinical reports have been collected. Drugs are used at their imputed therapeutic concentrations. Through this pan-KI transcriptomic profiling, we obtained insights into the affected pathways that may be related to KI-associated CT. We show that selective patterns of gene expression can be associated with the FAERS-derived clinical risk for KI-associated CT, which may be highly relevant to identify KI drug candidates at risk for showing clinical CT. We also describe the relationships between KI CT risk and structural properties of KIs, highlighting the potential for re-engineering small molecules that exhibit a high risk for CT.

## Results

### Differences in CT risk of kinase inhibitors

In order to obtain unbiased estimates of clinical risk of KI-associated CT, we analyzed individual adverse-event reporting data from FAERS (Fig. 1a). Reporting odds ratios (RORs) were derived based on the relative frequencies of AE occurrence of each KI compared to all KIs. These risk scores provide a relative ranking of KI-associated toxicity. Kinase inhibitors were shown to have pronounced differences in the relative risk of CT (Fig. 1b). When comparing the ranking of risk scores derived from FAERS with adverse drug-reaction (ADR) reporting data from the World Health Organization (WHO) ADR reporting database, we find that the ranking from these databases largely agrees (Fig. 1c), indicating the general consistency of the clinical risk scores across databases.

#### Figure 1 — Cardiotoxicity of protein kinase inhibitors

> **Fig. 1a** — 流程示意图：从 FDA FAERS 数据库量化激酶抑制剂临床心脏毒性风险评分的方法概述。左侧展示了 FAERS 数据库中包含超过 900 万份不良事件报告，通过构建 2×2 列联表（特定 KI 与其他 KI 的心力衰竭相关 ADR 频率），计算报告比值比（Reporting Odds Ratio, ROR），从而得出每种 KI 相对于其他 KI 的心脏毒性相对风险。右侧列出了所分析的 23 种激酶抑制剂名称。
>
> **Fig. 1b** — 水平条形图：展示 23 种激酶抑制剂按心脏毒性风险排序的报告比值比（ROR）。图中每种 KI 旁标注了 ROR 的均值和 95% 置信区间。Ponatinib、Bosutinib 和 Dasatinib 等显示出最高的心脏毒性风险（ROR > 10），而 Tofacitinib、Cabozantinib 和 Afatinib 等风险较低（ROR 接近 1 或更低）。x 轴为 ROR 对数刻度（0.5–50）。
>
> **Fig. 1c** — 散点图：比较 FAERS 与 WHO 药物警戒数据库得出的激酶抑制剂心脏毒性风险排名的一致性。x 轴为 FAERS 的 ROR 值，y 轴为 WHO 的 ROR 值。两种数据库的风险排名总体一致，验证了临床风险评分的可靠性。
>
> **Fig. 1d** — 热图/矩阵：展示文献报道的体外和体内临床前实验对 KI 相关心脏毒性的预测结果与 FAERS 临床风险评分的对比。每行为一项研究（标注了模型系统、处理方式和参考文献），每列为一种激酶抑制剂（按心脏毒性风险排序）。红色/深色表示在临床药物浓度下与对照组相比毒性变化 >20%（预测有心脏毒性），蓝色/浅色表示变化 <20%（无显著毒性）。该图清晰地显示：临床前表型实验（如细胞存活率、搏动频率、细胞凋亡、组织病理学等读数）与 FAERS 衍生的临床风险评分之间没有可识别的关联——高临床风险的药物在实验中可能显示无毒性，反之亦然。
>
> **Fig. 1e** — 热图/散点图：展示作者自己进行的体外剂量-响应实验结果。选用具有不同 CT 风险的 KI 处理原代心肌细胞 48 小时，检测细胞存活率和线粒体应激。左上热图显示各 KI 在两种读数下的毒性（>20% 或 <20% 变化），右上散点图比较存活率读数与 FAERS ROR，右下比较线粒体应激读数与 FAERS ROR。结果显示这些表型实验同样与临床风险评分缺乏相关性（相关系数分别为 1.27 和 0.08），再次强调了需要基于早期分子 signatures 的替代方法来预测心脏毒性。

### Phenotypic assays poorly correlate with CT

We performed a literature review for in vitro and in vivo experimental datasets that aimed to predict CT risk based on phenotypic readouts, such as cell viability or beating rate from in vitro cardiomyocyte or animal models, to determine if such phenotypic experiments can predict the clinical risk scores for CT. Studies in which drugs at the clinical concentration induced more than a 20% change in various phenotypic readouts compared to control experiments were classified as predicting potential CT (Fig. 1d). Across these studies, it was apparent that there was no identifiable relationship between apparent experimental toxicity in comparison to the relative incidence of CT in patients as derived from FAERS.

We conducted dose–response experiments with selected KIs that had varying risks for CT using the cardiomyocyte cell lines that were used in the current study for transcriptomic profiling, quantifying cell viability, and mitochondrial stress after 48 h of exposure to the selected KIs. We again assessed if drugs caused more than a 20% change in cell viability and mitochondrial stress at the typical clinically used concentration (Supplementary Table 1). These studies showed a similar lack of correlation with clinical risk (Fig. 1e, Supplementary Fig. 1). These findings underscore the need for alternative approaches such as early molecular signatures for CT. This identified lack of the predictiveness of preclinical in vitro and in vivo phenotypic assays, as has been noted by others.

### Transcriptomic profiling of human primary cardiomyocyte-like cell lines

To study the transcriptomic response to KIs associated with CT, we obtained four primary cardiomyocyte lines that were isolated from ventricles of healthy adult human heart (two male and two female, PromoCell GmbH, Germany). Culture conditions, detailed phenotypic characterization of each cell line, including gene and protein expression, morphology, and functional assays, can be found on the DToxS Center website (www.dtoxs.org) under the "Cellular Metadata" section. Confluent cardiomyocyte-like cells were treated with drugs for 48 h at concentrations similar to their clinical concentration (Supplementary Table 1) with 3–4 replicates and 3–4 cell lines (Supplementary Table 2), after which RNA was extracted and sequenced using the 3′ digital gene-expression method (Fig. 2a).

#### Figure 2 — Overview of pan-KI transcriptomic profiling in human primary cardiomyocyte-like cells

> **Fig. 2a** — 实验流程示意图：展示转录组数据生成的实验方法。四种人原代心肌样细胞系（来自两名男性和两名女性）分别用 23 种 FDA 批准的蛋白激酶抑制剂处理 48 小时，随后进行 3′-DGE（数字基因表达）mRNA 测序。经过质量控制后，计算差异基因表达值。每种药物按绝对平均倍数变化值排序，保留前 250 个基因，共完成 1309 次实验。
>
> **Fig. 2b** — 箱线图：展示 PromoCell 心肌细胞系（102 个样本）的基因表达特征与 GTEx 数据库中各人体组织（17,382 个总样本）基因表达特征的 Jaccard 相似性比较。箱线图显示各组织类型的两两 Jaccard 系数的上下四分位范围。结果显示：PromoCell 心肌细胞的基因表达与血液（rank 2）、肌肉（rank 4）和心脏组织（rank 10）最为相似，说明这些细胞系能够提供与心脏细胞可比较的基因表达变化。
>
> **Fig. 2c** — 热图：展示所有 KI 两两之间在排名前 250 的差异表达基因上的 Jaccard 指数。热图颜色编码表示基因重叠程度。结果表明各 KI 之间排名靠前的差异基因重叠度有限（Jaccard 指数 <0.25），意味着不同 KI 引起的转录组响应在很大程度上是不同的。
>
> **Fig. 2d** — 三维主成分分析（PCA）散点图：基于所有 KI 的完整平均倍数变化基因表达谱绘制的前三个主成分（PC1、PC2、PC3）。不同 KI 用不同颜色和标签表示。结果显示：9 种 KI（如 Ponatinib、Bosutinib、Vandetanib、Axitinib、Sunitinib、Dasatinib、Tofacitinib、Afatinib、Erlotinib 等）在 PC 空间中表现出较大的基因表达变异，而其余 KI（如 Imatinib、Gefitinib、Crizotinib、Regorafenib、Sorafenib、Vemurafenib、Pazopanib、Lapatinib、Nilotinib、Trametinib 等）基因表达变化较小——尽管其中一些 KI 的心脏毒性已被充分证实。这表明仅依靠排名差异基因表达值无法充分提供与心脏毒性相关的基因表达谱信息。

We investigated if transcriptomic profiles of PromoCell cardiomyocytes are related to human heart tissue and hence a good model to study CT. We compared the gene-expression similarity of untreated PromoCell cardiomyocytes against tissues available in the Genotype-Tissue Expression (GTEx) project, which contains gene-expression data from many human tissues, including the heart. Using the Jaccard distance for the top expressed 250 genes (based on transcript per million counts) for both untreated PromoCell and GTEx tissues, we observe that PromoCell cardiomyocytes' expression exhibits a gene expression similar to blood (rank 2), muscle (rank 4), and heart (rank 10) tissue. Based on these results, we conclude that the PromoCell cardiomyocytes can offer comparable gene-expression changes to that of cardiomyocytes.

### Limited overlap in differentially expressed genes across KIs

Differential gene-expression fold-change values were computed across the four cell lines. Initial analyses showed that the DEGs generally clustered more strongly by drugs than by cells. We calculated median fold-change values for each KI across cell lines, resulting in a single gene-expression profile for each KI. Ranked gene lists for each KI were generated by ranking by differential gene-expression p value and keeping the top 250 genes. To assess the similarity between genes present in the top 250 genes for each KI, the Jaccard index was calculated for each ranked list of KI-specific genes, which indicated a limited overlap (<0.25) between the top 250 genes across KIs (Fig. 2c). Principal component analysis showed variable gene-expression patterns for nine KIs, while for the remaining KIs, little variation in gene expression was seen (Fig. 2d), even though these remaining KIs included drugs for which CT is well established. We concluded that ranked differential gene-expression values would not be sufficient to provide clear insights into gene-expression profiles associated with CT.

### Pathways correlated with KI-associated CT

To identify pathways and subcellular processes across KIs and their potential involvement with CT, we performed enrichment analysis for protein kinases and KEGG terms using the top 250 differentially expressed genes ranked by p value across cell lines and KIs. We then correlated p values of enriched terms with clinical FAERS-derived risk scores to identify potential kinases and pathways associated with CT risk (Fig. 3a). The protein kinase LIMK2, which is involved in actin cytoskeleton reorganization pathways, ranked the highest in its correlation specifically enriched for KIs with a higher risk score (Fig. 3b). Sucrose- and pyruvate-metabolism pathways were the most strongly enriched pathways correlating with high risk scores (Fig. 3c). However, since no directionality in pathways is considered in these enrichment analyses, both the positively and negatively correlated processes may play a role in the development of CT. When considering enriched protein kinases and KEGG processes across all KIs, multiple pathways were identified (Supplementary Fig. 2). These findings indicate that there is likely substantial complexity underlying the action of KI in cardiomyocytes, although currently these analyses remain correlational and do not offer proof of causal relationships.

#### Figure 3 — Analysis of transcriptomic profiling data in relation to cardiotoxicity risk

> **Fig. 3a** — 流程示意图：展示转录组数据分析流程。每种激酶抑制剂在所有细胞系中的前 250 个差异表达基因（按 p 值排序）→ 进行富集分析（KEA 用于激酶富集，KEGG 用于通路富集）→ 将富集 p 值与临床心脏毒性风险评分进行相关性分析。
>
> **Fig. 3b** — 横向条形图：展示与心脏毒性风险评分显著相关的富集激酶。x 轴为与心脏毒性风险评分的相关系数（r），图中列出了 25 个蛋白激酶及其富集 p 值。LIMK2（r ≈ 0.45）和 ERN1（r ≈ 0.40）排名最高，与较高的心脏毒性风险正相关。其他显著激酶包括 TRIO、BCR、PDGFRB、CDK6、CSNK2A1、TESK1、PIK3CG 等。
>
> **Fig. 3c** — 横向条形图：展示与心脏毒性风险评分显著相关的富集 KEGG 通路。x 轴为相关系数（r），图中列出了 25 个 KEGG 通路。淀粉和蔗糖代谢（Starch and sucrose metabolism，r ≈ 0.50）和丙酮酸代谢（Pyruvate metabolism，r ≈ 0.45）与高风险评分呈最强正相关。类固醇激素生物合成、胰高血糖素信号通路、神经营养因子信号通路等也显示正相关。这些结果提示代谢相关通路在 KI 诱导的心脏毒性中可能发挥重要作用。

### Transcriptomic signature to predict CT risk

We tested if our KI-wide fold-change gene-expression profiles correlated with the KI-specific clinical risk scores for CT to identify a predictive transcriptomic signature for CT risk. Given the limited similarity between top-ranking gene-expression profiles across KIs, the entirety of the gene-expression profiles for different KIs were considered as potential predictors for KI-associated CT risk. KI-specific expression profiles of 10,749 genes were available as potential predictors for KI-specific CT risk scores. To identify genes most strongly associated with CT risk, we used an elastic net-penalized regression approach, which aims to select the most predictive variables while avoiding overfitting.

A two-stage regression analysis was performed (Fig. 4a). From the available 23 KIs with the associated clinical CT risk scores, we randomly left out 2 KIs for external validation of the model (test set, 10% of data). The differential gene-expression profiles of 21 remaining KIs were then used to train the model. Given the limited number of available drugs, small changes in expression patterns for drug were expected to affect the identity of the overall set of predictor genes. Therefore, we generated bootstrap datasets by random resampling of KI risk and the associated gene-expression profiles. These bootstrapped datasets were then fit using elastic net models. This first step was performed to identify gene-based predictors that could consistently predict CT risk and contributed significantly to the prediction of this risk. The bootstrap analysis resulted in stable selection of potential predictors. Predictors to be included in the final elastic net regression model were selected based on their minimal root-mean-squared prediction error (RMSE) after cross-validation. Based on this cross-validation, the gene-expression-based predictors in the final elastic net models consisted of 26 genes with the associated variable importance values (Fig. 4b).

Repeated cross-validation analyses indicated good predictive performance of the model for left-out KIs (Fig. 4c). We evaluated our 26-gene signature for predicting CT risk on an independent validation set of six KIs, of which three KIs were previously untested (Fig. 4d). We note that the independent validation set was performed 1 year after the original signatures were generated, using a different experimental protocol for the transcriptomic assay that was based on mRNA detection using random primers. We observed accurate predictive performance for five out of six KIs tested. The outlier, ibrutinib, had the lowest, albeit acceptable, predictive performance, with an error of 0.493 between the predicted and observed risk scores. Taken together, the developed signature can be of relevance to support risk prioritization of newly developed KIs. When we tested which of the 21 KIs drove the prediction strength of the model, we found that excluding any of four low-CT risk drugs (cabozantinib, tofacitinib, pazopanib, and erlotinib) increased the error substantially, indicating that these KIs contribute distinct information to the signature. In contrast, several of the high-ranking CT drugs could be excluded without sacrificing accuracy (Supplementary Fig. 3).

#### Figure 4 — Regression analysis for transcriptomic signatures to predict clinical risk

> **Fig. 4a** — 流程示意图：展示弹性网络回归分析的整体流程。21 种激酶抑制剂的中位数倍数变化基因表达谱（10,749 个基因）→ 基因 × 药物差异表达矩阵 → 自举弹性网络回归 → 生成心脏毒性基因 signature（26 个基因）→ 用 6 个留出的激酶抑制剂验证 signature 的预测能力。
>
> **Fig. 4b** — 横向条形图：展示 26 个被选为预测心脏毒性风险评分的转录组 signature 基因及其相对重要性（变量重要性百分比）。图中列出了所有 26 个基因：MFN2、DNAJB14、UBE2Q2、FAM45BP、URM1、FASTKD2、ZSCAN25、ARHGAP12、SNX14、SCO1、PRKRIP1、ZNF776、ZNF567、PMS2、CSE1L、NFATC2IP、DSTNP2、DDX3Y、GOLPH3L、BCAT2、RPL13P5、DHRS4、CASK、ZNF611、SNORD47。部分基因（如 MFN2、DNAJB14、UBE2Q2 等）显示出较高的相对重要性。
>
> **Fig. 4c** — 散点图：展示弹性网络交叉验证分析中观察到的风险评分与预测风险评分的比较。x 轴为观察到的风险评分（0–2.0），y 轴为预测的风险评分（0–1.5）。图中每种 KI 用不同颜色的点标注，对角虚线表示完美预测。结果显示良好的预测性能（RMSE = 0.137，R² = 0.9）。Ponatinib（观察风险最高）和 Bosutinib 被准确预测，低风险药物如 Tofacitinib、Vemurafenib 等也落在对角线附近。误差线显示了重复交叉验证的标准差。
>
> **Fig. 4d** — 条形图：展示 6 个激酶抑制剂外部验证的结果。x 轴为 KI 名称（REG=Regorafenib, SUN=Sunitinib, IBR=Ibrutinib, LEN=Lenvatinib, NIN=Nintedanib, CER=Ceritinib），y 轴（上方两个图）分别为观察到的风险评分和预测的风险评分，y 轴（下方图）为绝对误差。五种 KI（REG、SUN、LEN、NIN、CER）的预测准确性良好，仅 Ibrutinib 的预测误差较大（0.493），这与 Ibrutinib 作为独特的 VI 型共价结合抑制剂且高度多靶点结合的特性有关。

We then used the 26-gene signature to construct a protein–protein interaction network analysis to identify protein kinases and transcription factors associated with the signature (Supplementary Fig. 4). Several protein kinases were retrieved that are both known targets of the studied KIs, and which may be associated with the occurrence of KI-induced CT.

### Chemical structures of KIs inform CT risk

Off-target binding or polypharmacology is commonly observed in KIs. Since off-target binding is dependent on the structure of the drug, we investigated the relationship between kinase inhibitor chemical structure, binding target profile, and CT risk. To do this, we generated a structure–activity–similarity (SAS) map of the 26 tested inhibitors (in both the training and validation set) and their CT-risk score (Fig. 5a). SAS maps can be divided into four quadrants: the upper-left quadrant shows KI pairs with low chemical similarity and large changes in CT risk. The lower-left quadrant describes largely dissimilar KI pairs with small changes in CT risk. The lower-right quadrant describes KI pairs that exhibit a "smooth" structure–activity relationship, that is, small changes in chemical similarity are associated with small changes in CT risk. Finally, the upper-right quadrant indicates highly chemically similar compounds with large changes in CT risk. KI pairs in the upper-right region represent activity cliffs, that is, that small changes in chemical structure are associated with large changes in CT risk. In this region, we find several KI pairs, in particular, we observe large activity cliffs between afatinib and bosutinib as well as bosutinib and erlotinib. Here, all four compounds have the same chemical core (Fig. 5b); however, both afatinib and erlotinib show respectively lower CT risk scores compared to bosutinib. We hypothesized that harmonization of drug substructure, similarity, and promiscuity in the context of kinase inhibitor type may inform on our ability to predict CT risk (Fig. 5c).

By investigating their KI target profiles, we observe that both afatinib and erlotinib are less promiscuous KIs compared to bosutinib (which is one of the most promiscuous KIs in this set, Fig. 5d), and they both inhibit EGFR at nanomolar concentrations. On the other hand, less promiscuous KIs, such as lapatinib and gefitinib, exhibit a comparably lower CT risk score (Fig. 5e). Indeed, we observe a correlation between kinase inhibitor promiscuity and the observed CT risk score (Supplementary Fig. 5). However, KI promiscuity may not be the sole determinant of CT risk. For example, KIs such as imatinib and nilotinib are not as promiscuous as bosutinib; however, both exhibit relatively high CT risk scores. In this case, both imatinib and nilotinib CT may be explained due to their similar chemical structure and high specificity for protein kinases such as DDR1 and ABL.

Finally, kinase inhibitors have distinct binding modes against their targets. Kinase inhibitors that bind their kinase targets can be classified based on their binding mode, including the kinase conformation they bind and/or type of interactions they make with their kinase targets (e.g., covalent vs. non-covalent). For example, type I inhibitors bind an active kinase conformation, while Type I1/2, II–V bind distinct inactive states; type VI KI binds the kinase target covalently. We do not observe a clear relationship between kinase inhibitor-binding mode and CT. For example, the type II inhibitors imatinib and nilotinib are observed to have a high CT risk, while the type II inhibitors sorafenib and regorafenib have comparatively lower observed CT risk. However, both pairs of inhibitors are highly chemically similar and have similar binding targets. Taken together, the observed CT risk of a KI may be related to both a kinase inhibitor's selectivity and its chemical structure. Furthermore, we observe a relationship between chemical structure and binding target similarity to the predictive performance of our signature (Fig. 5e–g).

#### Figure 5 — Structure–Activity–Similarity (SAS) maps of kinase inhibitor activity and cardiotoxicity

> **Fig. 5a** — SAS 散点图：展示 26 种激酶抑制剂两两之间的化学相似性（x 轴，Tanimoto 系数 Tc，基于 ECFP4、ECFP2、Daylight 和 MACCS 四种化学指纹的加权平均）与心脏毒性评分差异（y 轴，DCS）的关系。图被虚线分为四个象限：左上（化学不相似但毒性差异大）、左下（化学不相似且毒性差异小）、右下（化学相似且毒性差异小，即"平滑"的结构-活性关系）、右上（化学相似但毒性差异大，即"活性悬崖"）。图中展示了多种 KI 的配对点，右上象限中 Afatinib–Bosutinib 和 Bosutinib–Erlotinib 等配对显示出明显的活性悬崖。化学相似性阈值为 Tc = 0.38，DCS 阈值为最大值的一半 = 0.82。
>
> **Fig. 5b** — 化学结构式：高亮展示右上象限（活性悬崖区域）和右下象限中关键 KI 的化学骨架。图中展示了 Afatinib、Bosutinib、Erlotinib、Gefitinib、Vandetanib、Lapatinib 等具有相同或相似化学核心的药物结构，以及 Sorafenib/Regorafenib（化学结构相似）和 Imatinib 的化学结构。这些结构比较说明了具有相同化学核心的药物可能因结合选择性（promiscuity）的不同而导致心脏毒性风险评分的显著差异。
>
> **Fig. 5c** — 热图：展示激酶抑制剂基于化学相似性的层次聚类和其结合模式标注（如 Type I、Type I1/2、Type II、Type III、Type IV、Type VI）。热图颜色编码表示激酶抑制剂对 500 nM 浓度下各激酶的结合选择性分数。各 KI 按化学相似性层次聚类排列，并用罗马数字标注其结合模式。该图直观展示了具有相似化学结构的 KI 往往具有类似的激酶结合谱。
>
> **Fig. 5d** — 条形图：展示各激酶抑制剂在 500 nM 下的选择性分数（kinase selectivity score），反映 KI 的多靶点结合程度（即"promiscuity"）。Bosutinib、Ibrutinib、Dasatinib、Ponatinib 等显示出较高的 promiscuity（选择性分数较低），而 Gefitinib、Erlotinib、Lapatinib、Afatinib 等显示出较高的选择性（选择性分数较高）。
>
> **Fig. 5e** — 聚类热图：展示各激酶抑制剂观察到的（归一化为零的）心脏毒性风险评分，按化学相似性层次聚类排序。高风险药物（如 Bosutinib、Ponatinib、Sorafenib、Nilotinib、Imatinib）聚在一起，低风险药物（如 Gefitinib、Afatinib、Erlotinib、Lapatinib）聚在一起。
>
> **Fig. 5f** — 聚类热图：展示各激酶抑制剂预测的（归一化为零的）心脏毒性风险评分，按与 5e 相同的层次聚类排序。预测风险与观察风险的总体模式一致，验证了 signature 的预测能力。
>
> **Fig. 5g** — 聚类热图：展示观察到的与预测的心脏毒性风险评分之间的绝对误差。大部分 KI 的绝对误差较小，表明 signature 具有良好的预测准确性。化学结构独特（如 Ibrutinib，VI 型共价抑制剂）的药物误差相对较大。

## Discussion

The occurrence of drug treatment-associated CT, leading to decreased cardiac function, follows the therapeutic effects of the drugs, and is only observed in a subset of the patients using the drug. This raises the question of whether it would be possible to obtain early cell-based signatures predictive for drug toxicity. Here we addressed this question by attempting associating drug treatment-induced gene-expression patterns with the clinical risk for the adverse events of interest.

By estimating clinical risk from the FAERS database, our method utilizes a relevant and unbiased approach for the quantification of CT risk. As a result, our CT risk scores lack notable pitfalls such as selection bias associated with tightly controlled clinical trials, which underestimate adverse-event risks due to cohort size, trial duration, and selective inclusion criteria for subjects. Nevertheless, there are limitations to the FAERS database as well, which we have discussed and addressed in previous work. Specifically, use of the FAERS resource may confound demographics information such as age and sex, which was observed not to vary across different KIs. Moreover, CT risk score does not reflect absolute risk for developing CT. Rather, it reflects the relative risk for a subset of patients for which drug-associated adverse events were reported. In addition, there may be some systematic biases based on the sampling frequency of drugs by institution.

It remains unclear if all KIs induce CT through similar mechanisms, and to what extent ultimate clinical pathologies are similar. While the FAERS database allows us to distinguish between different types of CT, the annotation is not uniform and may either refer to distinct pathophysiological descriptions or rather more general clinical presentations of heart failure. To this end, we chose to lump all forms of heart failure, while excluding cardiac AEs that have known and unrelated origin such as coronary artery disease and arrhythmias.

We compared KI-associated transcriptomic response profiles generated from cultured human primary cardiomyocyte-like cells with clinical CT risk scores to obtain a reduced set of genes that may predict the relative risk for KI-associated CT. Using the clinically weighted signatures and the associated regression coefficients identified in the elastic net model, the relative risk for CT can be predicted. The risks predicted by our signatures and the associated regression model can be used in drug development to rank the potential risk of novel KIs with respect to existing KIs with better characterized clinical risks for CT.

The signatures generally showed good prediction of CT risk during cross-validation as well as on an independent set of KIs (Fig. 4), while the only poorly performant KI, ibrutinib, inhibitor of Bruton nonreceptor protein-tyrosine kinase, represents a unique KI in terms of binding mode (i.e., type VI inhibitor) and high promiscuity (Fig. 5). Specifically, it is a member of an emerging class of kinase inhibitor drugs that bind their targets covalently (type VI KIs). These drugs are highly underrepresented in the databases used in this analysis, explaining the misclassification of ibrutinib.

The four cell lines we studied are insufficient to fully capture such human variability to KIs. Therefore, in our analysis, we used median fold-change gene-expression profiles across multiple cell lines. The resulting averaged gene-expression profiles thus reflect relatively consistent changes in gene expression across cell lines, i.e., changes in gene expression that are less likely to be highly variable across cell lines, yet may also reflect a set of predictors that may be more consistent in the population. Given that the FAERS CT risk scores also reflect a population-level CT risk, the use of these median values in fold-change gene-expression values is a reasonable starting point for our analyses.

The experimental underpinning of the transcriptomic profiles generated in this study makes them likely to be of value in selecting drug candidates with a low risk for CT as an adverse event. Our analysis is based on primary human heart-derived cardiomyocyte-like cells. Although these cell lines do have phenotypic limitations due to dedifferentiation, the signatures obtained from the cells could be relevant for prediction of clinical drug effects. These cell lines may be reflective of human cardiac pharmacology, i.e., in comparison with animal-derived cardiomyocytes, even though further characterization and standardization are still needed. Detailed characterization of these cell lines is available as metadata to the RNAseq datasets at www.dtoxs.org. Our analyses used drug exposures similar to clinically reported maximum plasma concentrations of the individual KIs, rather than using the same concentrations for all KIs, even though we did not correct for protein binding. We expect that the duration of 48-h exposure may reflect transcriptomic changes that are likely related to early changes in subcellular processes associated with the adverse event of interest. Unfortunately, in this study, it is not feasible or ethically possible, due to lack of prior informed consent, to compare cardiac gene-expression signatures with gene-expression profiles from patients who receive KI-therapy and/or developed KI-associated CT. We considered whether we could compare our gene-expression signatures to cardiac gene-expression data from patients with heart failure who undergo surgery. Typically these are patients with advanced disease, and the gene expression in tissue from advanced disease is not likely to be of relevance to acute drug-induced CT.

By investigating the chemical structure and binding profile similarity of KIs, we are able to observe that chemical components and scaffolds that lead to promiscuous binding of KIs to multiple binding targets are correlated with higher CT values. This is consistent with the notion that a portion of CT risk of KIs can be attributed to higher levels of off-target interactions. Indeed, when we investigate the binding profile of three chemically similar KIs: afatinib, erlotinib, and gefitinib, we find that their binding profiles are fairly specific compared to other KIs, and they have a lower normalized CT risk score. One limitation we have observed with our approach is that chemically distinct KIs (e.g., in terms of binding profile, substructural similarity, and type), such as the type IV inhibitor ibrutinib, exhibit diminished predictive performance. However, we think that using the guidelines we provide herein, this signature could still assist in the development and prioritization of KIs with lower toxicity risks.

We cautiously anticipate that clinically weighted transcriptomic signatures such as those developed in this study may be of relevance to guide safety assessment in early drug development. Unlike the relatively well-established assessment of electrophysiological safety issues such as QT prolongation, the assessment of non-QT type of CT associated with KI and other novel drugs lacks reliable biomarkers. The transcriptomic signature for CT identified in this study may help fill this gap, especially if its structure and binding profiles are closely represented within the inhibitors in this study. One could anticipate that after initial selection of promising KIs with apparent efficacy in preclinical screens, transcriptomic profiling using the signatures developed here may possibly be used to rank drugs for the expected CT risk and exclude those with high CT risk scores (Supplementary Fig. 6).

While beyond the scope of this study, future extension of our studies could explore the idea of studying individualized risk scores for CT. That is, do baseline gene-expression profiles of larger libraries of patient-derived cardiomyocyte cell lines predict the difference in risk for CT between individual patients? Ideally, such an analysis would be conducted using induced pluripotent stem cell-derived cardiomyocytes from patients, who have received KIs and experienced different levels of CT, such as was recently described for anthracycline chemotherapeutics. This would then further enable the development of precision medicine approaches to KI therapy that could minimize the risk for CT.

## Methods

### Cell culture and drug treatment

Adult human cardiomyocytes (Cat#: C12810) were purchased from PromoCell GmbH (Heidelberg, Germany) and grown in culture as per the manufacturer's instructions. Four different cell line lots (Lot#: 3042901.2, 4031101.3, 2082801.2, and 2120301.2) isolated from two male and two female subjects were cultured under serum-free differentiation conditions for at least 28 days prior to drug treatment. Details regarding metadata information, including cell line metadata and the quality control and assurance metrics, can be found on www.dtoxs.org.

### Dose–response experiments

For two of the four cell lines, dose–response experiments were conducted treating cells for 48 h with eight increasing perturbagen concentrations (5 nM, 50 nM, 100 nM, 500 nM, 1 µM, 5 µM, 10 µM, and 100 µM) and vehicle-treated control, in quadruplicates. We assayed for viability through image-based analysis of nuclear counts with Hoechst 33342 (ThermoFisher, Cat#: H3570) and MitoTracker Red (ThermoFisher, Cat#: M22426) for mitochondrial toxicity. Details of the experimental protocols for cell culture, drug treatment, and transcriptomics have been described as step-by-step standard operating procedures for the various experiments available on www.dtoxs.org.

### Transcriptomics

Cells were treated for 48 h with a single perturbagen concentration around the maximal concentration (Supplementary Table 1). After drug treatment, the cells were lysed, RNA was collected using TRIzol, and gene-expression profiles were measured using the 3′ digital gene-expression method.

### Sequence alignment and processing of gene-expression data

The raw sequence reads were aligned to (99.75%) genes. Subsequently, the resulting alignment files were parsed to identify the fragments with acceptable alignment quality, to remove duplicate fragments, and to assign accepted fragments to the corresponding genes. The resulting read-count (i.e., transcript count) table was then subjected to correlation analysis at each treatment condition, to identify and remove outlier samples, determined by predefined thresholds for robust multi-array average (RMA) normalization, gene filtering, and outlier removal using R (version 3.4.3) and the edgeR R package. These procedures are described elsewhere, and step-by-step protocols are available on www.dtoxs.org. The resulting normalized and log-transformed fold-change gene-expression values for each sample are also deposited for public access to the DToxS data repository (www.dtoxs.org).

### Processing and exploratory analysis of gene-expression data

The median log-transformed gene-expression fold-change value was calculated across all cell lines for each individual KI. The resulting matrix of gene fold-change values by KIs was used for the regression analysis. To obtain insight into the general patterns present in this KI-perturbed transcriptomics dataset, we generated rankings of the top 500 genes for each drug, by their absolute mean fold-change value, i.e., whether positive or negative. For each of these KI-associated rankings, we determined the frequency of these changes being also present in the ranking of other drugs, e.g., the similarity in genes present in the top 250 gene lists for each KI. This was visualized using the Jaccard index, and by plotting the most highly differentially expressed genes for 26 kinase inhibitors. For each pair of kinase inhibitors, we computed Tc based on RNA sequencing data from PromoCell cardiomyocyte exposures to kinase inhibitors. The top 500 genes for a KI were set as 1, while genes that were not in the top 500 were set as 0.

### Calculation of clinical risk RORs

Adverse-event frequencies from the FDA Adverse Event Reporting System (FAERS) were obtained from the AERSmine resource, which contains a curated version of the FAERS database. ADRs in the FAERS database are organized according to MedDRA, which is a hierarchical ontology to classify ADRs from high-level organs associated with the pathology to reported low-level specific pathological conditions. We downloaded the frequencies of the occurrence of ADRs for all protein KIs available in FAERS, together with all other frequencies of ADRs reported for these KIs. A time-stamped record of this download to reproduce this analysis was retained. RORs were then computed for each KI using the frequency f_dt of the ADR of interest, the frequency f_dn of any other ADR occurring, the frequencies f_nt of occurrence of the ADR of interest for any other protein kinase inhibitor, and the frequency f_nn for all other ADRs and KIs. The ROR was calculated using Eq. (1):

ROR = (f_dt / f_dn) / (f_nt / f_nn)

whereas the standard error (SE) of the log ROR was calculated using Eq. (2):

SE(log ROR) = sqrt(1/f_dt + 1/f_dn + 1/f_nt + 1/f_nn)

with the log-transformed confidence interval (CI) being calculated as follows:

CI = log(ROR) ± 1.96 × SE(log ROR)

Adverse events in FAERS are mapped to the MEDDRA dictionary. CT events related to heart failures and cardiomyopathies, excluding arrhythmogenic ADRs and coronary artery disorders, were selected from the main MEDDRA cardiac ADR group. The selected ADRs primarily reflected different stages of heart failure, which were grouped together.

### Elastic net regression analysis

The FAERS-derived risk RORs for CT were regressed against the KI-associated vectors of mean fold-change values across the four cell lines. A two-step regression procedure was then used to select predictor genes reducing the sensitivity to changes in dataset composition. For this, we first generated 1000 bootstrap datasets with replacements for gene expression–KI risk score pairs. Each of these bootstrap datasets was fit using an elastic net regression model (R version 3.4.3, package glmnet, version 2.0-16). The genes that were selected as predictors (i.e., nonzero regression coefficient) and the scaled values of the gene-associated coefficients were saved for each bootstrap dataset. Across all bootstrap datasets, the relative frequency of the selection of gene-based predictors, and the mean-scaled coefficient value was computed. We then calculated the product of the mean frequency and scaled coefficient value, rank predictors by their importance with respect to robustness (selection frequency). A large number of percentiles of these rankings were evaluated using leave-one-out cross-validation. The selection percentiles that yielded minimum RMSE were selected to generate expression signatures. The selected predictor genes were then ranked by their relative importance, and by their median fold-change values, and displayed as clustered heatmaps. We finally evaluated the predictive value of the resulting regression model to predict CT risk scores for the two left-out KIs.

When using this approach to analyze similar datasets of cardiomyocyte transcriptomics together with risk scores, it is possible that potentially different gene signatures are identified than those described in the current report. This difference is associated with the intrinsic property of penalized regression approaches that select predictors from potentially highly correlated sets of predictor candidates. Hence, small changes in either risk scores or gene-expression datasets may affect correlation structures of the data and thereby the list of genes for a signature.

### Enrichment and network analyses

Enrichment analysis was performed based on a one-tailed Fisher's exact test using R (package stats), in order to identify enrichment of specific genes in predefined gene lists. For enrichment of pathways and biological processes, we used the KEGG database (2016), and for enrichment of protein kinases, we used the KEA database (2015). Diseases were excluded from the KEGG list of processes (e.g., diabetes, depression, and cancer), in order to only evaluate general biological processes or pathways. We used the top 250 DEGs ranked by p value for each KI to perform enrichment analysis. Subsequently enriched term p values were correlated with CT risk scores to identify kinases and pathways associated with CT risk.

The gene part of the signature for CT identified in the regression analysis was used as seed note to perform a protein–protein interaction network (PPI) analysis, conducted using the web application X2K, which aims to identify associated kinases and transcription factors based on multiple PPI databases. Computed as a result of this analysis from computational methods, the relationships between kinases and transcription factors were derived through analysis of a computational compendium of binary PPI data.

### Calculation of chemical similarity

RDKit was used to generate chemical fingerprints for 26 kinase inhibitors. For each pair of kinase inhibitors, we computed Tc values using four chemical fingerprints, including Morgan_2 (2048-bit, ECFP4), Morgan_1 (2048-bit, ECFP2), Daylight-like, and MACCS. Because each of these fingerprints capture distinct chemical properties, we computed a weighted Tc average of the three fingerprints: 30% ECFP4, 30% ECFP2, 30% Daylight-like, and 10% MACCS, which exhibited the most optimal spread of the distribution of the pairwise distances. To generate the SAS maps (Fig. 5a), we plotted the pairwise-weighted Tc values with their difference in CT scores (DC_T). Finally, 0.35 was set as the threshold for chemical similarity, while half of the maximum difference was set as the threshold for DC_S. Chemical structures were drawn using Marvin based on SMILES strings obtained from PubChem.

### Calculation of KI-binding target similarity

Kinome-wide kinase inhibitor-binding (K) profiling data were obtained from Klaeger et al., which consisted of kinome-wide binding (Kd) profiling data for all of the tested kinase inhibitors across 242 kinases. A heatmap was generated for selected kinase inhibitors based on the negative log of the Kd values from Klaeger et al. Notably, the Kd values were scaled by 100,000 to avoid negative log values.

### Reporting summary

Further information on research design is available in the Nature Research Reporting Summary linked to this article.

## Data availability

All processed RNAseq data and the curated version-controlled standard operating procedures featured in this study can be downloaded freely at (www.dtoxs.org) or the LINCS Data Portal. Raw transcriptomics data can be accessed through the Gene Expression Omnibus (GEO) repository with accession numbers GSE146096 and GSE146097. Source data for each figure are provided with this paper. All remaining data will be available from the corresponding author upon reasonable request. Source data are provided with this paper.

## Code availability

All scripts are open-source and available from the DToxS GitHub repository (https://github.com/dtoxs).

## Table 1 — Overview of KIs included in this analysis

| Drug | Three-letter code | Approval year^a^ | Therapeutic targets | Concentration (µM)^b^ |
|------|-------------------|-----------------|---------------------|----------------------|
| Afatinib | AFA | 2013 | ErbB2 and EGFR | 0.05 |
| Axitinib | AXI | 2012 | VEGFR1/VEGFR2/VEGFR3/PDGFRB/c-KIT | 0.2 |
| Bosutinib | BOS | 2012 | Bcr-Abl and SRC | 0.1 |
| Cabozantinib | CAB | 2012 | c-Met and VEGFR2 | 2 |
| Ceritinib | CER | 2014 | ALK | 1 |
| Crizotinib | CRI | 2011 | ALK and HGFR | 0.25 |
| Dabrafenib | DAB | 2013 | BRAF | 2.5 |
| Dasatinib | DAS | 2006 | ABL, ARG, KIT, PDGFRα/β, and SRC | 0.1 |
| Erlotinib | ERL | 2004 | ErbB1 | 3 |
| Gefitinib | GEF | 2003 | ErbB1 | 1 |
| Imatinib | IMA | 2001 | Bcr-Abl | 5 |
| Lapatinib | LAP | 2007 | ErbB1 | 2 |
| Nilotinib | NIL | 2007 | Bcr-Abl | 3 |
| Pazopanib | PAZ | 2009 | VEGFR2, PDGFRα/β, and KIT | 10 |
| Ponatinib | PON | 2012 | Bcr-Abl, BEGFR, PDGFR, FGFR, EPH, SRC, c-KIT, RET, TIE2, and FLT3 | 0.1 |
| Regorafenib | REG | 2012 | RET, VEGFR, and PDGFR | 1 |
| Ruxolitinib | RUX | 2011 | JAK | 1 |
| Sorafenib | SOR | 2005 | BRAF, VEGFRs, PDGFRα/β, FLT3, and KIT | 0.5 |
| Sunitinib | SUN | 2006 | VEGFR, PDGFR, CSF1R, FLT3, and KIT | 1 |
| Trametinib | TRA | 2013 | MEK1 and MEK2 | 0.1 |
| Tofacitinib | TOF | 2012 | JAK | 1 |
| Vandetanib | VAN | 2011 | RET, VEGFR, and EGFR | 0.33 |
| Vemurafenib | VEM | 2011 | BRAF | 2 |

^a^ US approval date, first indication.
^b^ Derived from maximum total (bound + free) plasma concentrations in humans as reported in the literature. Table S3 lists the purity and literature references to clinical concentrations.
