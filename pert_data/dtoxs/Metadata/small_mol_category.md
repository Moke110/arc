# DTOX 小分子药物结构分类

DToxS 数据集 36 种小分子按化学母核/药理合并为 **6 大类**，用于 cold-drug split。

分类信息来源：`Small_Molecule_Metadata.txt` 的 SMILES/InChI 数据 + `DToxS_article.md`（Nature Communications 2020）中激酶抑制剂的化学核心/骨架聚类 + 非 KI 药物的已知药理学分类。

## Category 1: 喹唑啉类激酶抑制剂 (Quinazoline-based KIs)

共享 4-苯胺基喹唑啉 (4-anilinoquinazoline) 母核，主要抑制 EGFR/HER2。

| 代码 | 名称 | 靶点 | 样本数 |
|------|------|------|--------|
| AFA | Afatinib | ErbB2, EGFR | 16 |
| ERL | Erlotinib | ErbB1 (EGFR) | 14 |
| GEF | Gefitinib | ErbB1 (EGFR) | 18 |
| LAP | Lapatinib | ErbB1 (EGFR), HER2 | 27 |
| VAN | Vandetanib | VEGFR, EGFR, RET | 16 |

**5 药，91 样本**

## Category 2: 嘧啶母核激酶抑制剂 (Pyrimidine-core KIs)

以嘧啶/氨基嘧啶为母核，靶向 BCR-ABL、JAK、ALK、MEK、B-Raf、VEGFR 等。

| 代码 | 名称 | 靶点 | 样本数 |
|------|------|------|--------|
| CER | Ceritinib | ALK | 16 |
| DAB | Dabrafenib | BRAF | 16 |
| DAS | Dasatinib | ABL, SRC, KIT | 16 |
| IMA | Imatinib | BCR-ABL | 20 |
| NIL | Nilotinib | BCR-ABL | 16 |
| PAZ | Pazopanib | VEGFR, PDGFR | 16 |
| PON | Ponatinib | BCR-ABL, SRC | 16 |
| RUX | Ruxolitinib | JAK | 16 |
| TOF | Tofacitinib | JAK | 18 |
| TRA | Trametinib | MEK1, MEK2 | 16 |

**10 药，166 样本**

## Category 3: 其他杂环激酶抑制剂 (Other heterocyclic KIs)

母核多样（二苯脲、吡啶、吲哚、吲唑、喹啉、噻唑/吡啶并嘧啶、磺酰胺等），靶向 VEGFR/Raf/c-Met/PDGFR 等多靶点。

| 代码 | 名称 | 靶点 | 样本数 |
|------|------|------|--------|
| AXI | Axitinib | VEGFR, PDGFR | 16 |
| BOS | Bosutinib | BCR-ABL, SRC | 16 |
| CAB | Cabozantinib | c-Met, VEGFR2 | 14 |
| CRI | Crizotinib | ALK, HGFR (c-Met) | 19 |
| REG | Regorafenib | VEGFR, PDGFR, RET | 20 |
| SOR | Sorafenib | BRAF, VEGFR, PDGFR | 28 |
| SUN | Sunitinib | VEGFR, PDGFR, KIT | 22 |
| VEM | Vemurafenib | BRAF V600E | 16 |

**8 药，151 样本**

## Category 4: 类固醇及相关 (Steroid-related)

| 代码 | 名称 | 类别 | 样本数 |
|------|------|------|--------|
| EST | Estradiol | 雌激素 (steroid hormone) | 18 |
| URS | Ursodiol | 胆汁酸 (bile acid, steroid-related) | 12 |

**2 药，30 样本**

## Category 5: 神经系统药物 (Neuro/psychoactive drugs)

| 代码 | 名称 | 类别 | 样本数 |
|------|------|------|--------|
| DOM | Domperidone | 多巴胺 D2 拮抗剂 (benzimidazole) | 3 |
| DTP | Diethylpropion | 中枢兴奋剂 (phenethylamine) | 3 |
| FLE | Flecainide | 钠通道阻滞剂 (benzamide) | 9 |
| LOP | Loperamide | 阿片受体激动剂 (piperidine) | 4 |
| PAR | Paroxetine | SSRI (benzodioxole) | 7 |

**5 药，26 样本**

## Category 6: 其他治疗药物 (Other diverse therapeutics)

| 代码 | 名称 | 类别 | 样本数 |
|------|------|------|--------|
| ALE | Alendronic acid | 双膦酸盐 (bisphosphonate) | 3 |
| AMI | Amiodarone | 抗心律失常 (benzofuran) | 8 |
| CYC | Cyclosporine | 免疫抑制 (cyclic peptide) | 17 |
| ENT | Entecavir | 抗病毒 (nucleoside analog) | 3 |
| MIL | Milrinone | PDE3 抑制剂 (pyridine-carbonitrile) | 8 |
| OLM | Olmesartan | 血管紧张素 II 受体拮抗剂 (biphenyl-tetrazole) | 18 |

**6 药，57 样本**

---

**总计: 6 类，36 药，521 样本**
