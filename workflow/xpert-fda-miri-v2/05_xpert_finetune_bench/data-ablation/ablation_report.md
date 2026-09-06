# XPert finetune-corpus data-ablation — summary

Diagnostic experiment (ADR 0019), decoupled from the official full/ folds. v2's union is 958 rows / 45 drugs across the two finetune datasets (dtoxs_primary + dtoxs_ipsc), so the per-seed test sets are **drug-partitioned** (not row-uniform as in v3's 4-dataset ablation): the 45 drugs are split into a disjoint tripartition of 3 x 15 test drugs; each seed tests one third (every union drug tested exactly once across the 3 seeds), with a ~15%-of-union unified test set drawn only from that third. Train = a dataset-stratified f-subsample of the seed's pool over 5 train-size fractions {20%, 40%, 60%, 80%, 100%}; warm-start = union minus test (~814 rows, all 45 drugs); cold-drug = the 30 hold-out drugs only (~610-655 rows, strictly drug-disjoint at every fraction). Recipe identical to the official fold bench (4 epochs, lr 0.002, batch 32, cm.yaml, unimol+HG, step-04 L1000 pretrain init). The L1000-pretrain weight anchors the curve at 0 rows. This report changes no model used downstream.

Runs aggregated: 33 (30 finetune + 3 pretrain anchors).

## Scale curve (mean of 3 seeds per mode)

| mode | metric | f=20% | f=40% | f=60% | f=80% | f=100% | anchor |
|---|---|---|---|---|---|---|---|
| warm-start | RMSE | 0.628 | 0.619 | 0.621 | 0.611 | 0.603 | 0.685 |
| warm-start | PCC | 0.152 | 0.190 | 0.211 | 0.227 | 0.242 | 0.030 |
| warm-start | Spearman | 0.134 | 0.171 | 0.189 | 0.202 | 0.213 | 0.025 |
| warm-start | R2 | -5893753.042 | -2671934.281 | -963273.708 | -1217336.281 | -671136.318 | -735076448.000 |
| warm-start | PDS | 0.559 | 0.622 | 0.678 | 0.683 | 0.740 | 0.510 |
| warm-start | JAC | 0.075 | 0.091 | 0.095 | 0.100 | 0.119 | 0.041 |
| warm-start | FID | 0.694 | 0.714 | 0.736 | 0.739 | 0.758 | 0.531 |
| warm-start | REACH | 0.202 | 0.305 | 0.316 | 0.316 | 0.350 | 0.028 |
| cold-drug | RMSE | 0.629 | 0.632 | 0.636 | 0.631 | 0.638 | 0.685 |
| cold-drug | PCC | 0.136 | 0.157 | 0.171 | 0.180 | 0.184 | 0.030 |
| cold-drug | Spearman | 0.118 | 0.143 | 0.157 | 0.163 | 0.167 | 0.025 |
| cold-drug | R2 | -3342820.958 | -854645.313 | -1666132.021 | -1302432.516 | -4785267.000 | -735076448.000 |
| cold-drug | PDS | 0.537 | 0.562 | 0.589 | 0.651 | 0.646 | 0.510 |
| cold-drug | JAC | 0.074 | 0.076 | 0.082 | 0.083 | 0.089 | 0.041 |
| cold-drug | FID | 0.672 | 0.662 | 0.670 | 0.713 | 0.706 | 0.531 |
| cold-drug | REACH | 0.156 | 0.172 | 0.209 | 0.168 | 0.215 | 0.028 |

PDS / JAC / FID / REACH are the four raw (unscaled) higher-is-better Virtual Cell Challenge 2026 metrics (cell_eval2 vcc2026 profile), computed per drug from the per-cell predicted DEGs and averaged over the test drugs. Because the test set holds ~10 cells per test drug, the DE gene set for JAC/FID/REACH uses the top-50 genes by |log2FC| (not a Wilcoxon p_adj<0.05 set, which is underpowered at n=10); PDS is a cosine-retrieval rank and needs no significance test. No target-gene exclusion applies (drug perturbations). Drugs with very few rows (< 10) contribute proportionally fewer test cells.

## Warm-start drug-coverage diagnostic (n_absent_drugs)

Number of test drugs absent from the train subsample at that fraction (warm-start; cold-drug is strictly drug-disjoint by design, so its diagnostic always equals n_test_drugs).

| seed | mode | fraction | n_train | n_test | n_absent_drugs |
|---|---|---|---|---|---|
| 1 | cold-drug | 0.2 | 131 | 144 | 15 |
| 1 | cold-drug | 0.4 | 261 | 144 | 15 |
| 1 | cold-drug | 0.6 | 392 | 144 | 15 |
| 1 | cold-drug | 0.8 | 522 | 144 | 15 |
| 1 | cold-drug | 1.0 | 653 | 144 | 15 |
| 1 | pretrain | 0.0 | 0 | 144 | 15 |
| 1 | warm-start | 0.2 | 163 | 144 | 3 |
| 1 | warm-start | 0.4 | 326 | 144 | 0 |
| 1 | warm-start | 0.6 | 488 | 144 | 1 |
| 1 | warm-start | 0.8 | 651 | 144 | 0 |
| 1 | warm-start | 1.0 | 814 | 144 | 0 |
| 2 | cold-drug | 0.2 | 121 | 144 | 15 |
| 2 | cold-drug | 0.4 | 241 | 144 | 15 |
| 2 | cold-drug | 0.6 | 362 | 144 | 15 |
| 2 | cold-drug | 0.8 | 482 | 144 | 15 |
| 2 | cold-drug | 1.0 | 603 | 144 | 15 |
| 2 | pretrain | 0.0 | 0 | 144 | 15 |
| 2 | warm-start | 0.2 | 163 | 144 | 0 |
| 2 | warm-start | 0.4 | 326 | 144 | 0 |
| 2 | warm-start | 0.6 | 488 | 144 | 0 |
| 2 | warm-start | 0.8 | 651 | 144 | 0 |
| 2 | warm-start | 1.0 | 814 | 144 | 0 |
| 3 | cold-drug | 0.2 | 132 | 144 | 15 |
| 3 | cold-drug | 0.4 | 264 | 144 | 15 |
| 3 | cold-drug | 0.6 | 396 | 144 | 15 |
| 3 | cold-drug | 0.8 | 528 | 144 | 15 |
| 3 | cold-drug | 1.0 | 660 | 144 | 15 |
| 3 | pretrain | 0.0 | 0 | 144 | 15 |
| 3 | warm-start | 0.2 | 163 | 144 | 4 |
| 3 | warm-start | 0.4 | 326 | 144 | 1 |
| 3 | warm-start | 0.6 | 488 | 144 | 0 |
| 3 | warm-start | 0.8 | 651 | 144 | 0 |
| 3 | warm-start | 1.0 | 814 | 144 | 0 |

## Interpretation guide

- A metric that stops changing between f=80% and f=100% marks where more CM data stops paying off on the union-wide test distribution.
- warm-start vs cold-drug spread shows whether scale helps known-drug (memorised) prediction differently from unseen-drug generalisation.
- The anchor column is the un-finetuned L1000 pretrain weight (0 training rows), so the curve spans no-finetune to full-finetune.
- A warm-start drop at small f may be 'forgot a known drug' — read the n_absent_drugs diagnostic before attributing it to scale.
- **Warm vs cold scale ranges differ by design.** The cold-drug train pool strips every row of the 15 test drugs, so it holds only the 30 hold-out drugs (~610-655 rows vs ~814 warm): the cold curve's f=100% point is the hold-out-drug scale, NOT union scale, and the two rows' x-ranges must not be read as comparable full-scale spans.