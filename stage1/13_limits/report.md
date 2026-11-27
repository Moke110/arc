# v2 limitations analysis (step 13)

Generated: 2026-10-14 (deterministic, step 13)

Read-only diagnostic of the `xpert-fda-miri-v2` screen. Three figures, one structural limitation each.

## 1. 978-landmark pathway coverage

The pathway scorer only sees the genes of the 978 L1000-landmark panel. Resolving each of the 40 scoring pathways against that panel shows how little of each pathway is actually measured:

- Landmark coverage per pathway: **5.5% to 36.4%**, mean **16.0%** (median 15.0%).
- Median pathway: 88 genes total, only 8 scored.
- The 40 pathways contribute 607 scored gene-slots (with cross-pathway duplication) out of 4333 pathway genes; their union is 275 of the 978 landmarks.

**Limitation.** Pathway-level reversal is judged on a biased ~16% slice of each pathway. Pathways with a higher landmark share (e.g., GO-Fc Receptor Mediated Stimulatory Signaling Pathway at 36%) are favoured over those mostly absent from the panel; direction calls on the missing genes are impossible by construction.

## 2. Critical-pathway overlap

The step-10 pathway filter keeps drugs reversing >=6/8 *critical* pathways. Pairwise overlap is measured as shared genes over the ROW pathway's landmark subset (asymmetric; diagonal = 1):

- Off-diagonal overlap: **max 0.38** (Necroptosis -> Apoptosis, 9 genes), mean **0.08**.
- 22 of 56 ordered pairs share no gene at all.

**Limitation.** The 8 critical pathways are largely *independent*, not redundant: they are not restatements of one another. The >=6/8 filter is therefore a genuinely multi-axis bar (a drug must reverse six distinct processes), and the critical count is not inflated by duplicated gene sets. The corollary is that the filter is stringent and heavily trims the candidate list.

## 3. Screened-drug embedding structure

Each drug passing the filter is projected to PC1/PC2 from two related spaces (both per-dimension z-scored across the 64 drugs), with the 4 direct adrenergic agonists and 6 antagonists highlighted:

- **Left — the model's own input** `[HG(256) | mean-pooled projected UniMol tokens(256)]` (the tensor the finetuned XPert transformer consumes): PC1 50.2% + PC2 10.8% of variance.
- **Right — the HG component only**: PC1 13.3% + PC2 8.5% of variance.

**Limitation (a negative result, not a cluster).** In the space the model actually consumes the adrenergic ligands are NOT clustered (adrenergic-vs-rest Cohen's d on PC1 = +0.34); agonists and antagonists are intermixed. The class signal is real but lives in the HG branch alone (HG-only d = +1.37); it is diluted by the 512 UniMol-token dimensions in the joint input (UniMol-only d = +0.33). So the HG embedding that shapes the prediction is class-correlated, while the joint input the transformer sees does not expose the class linearly — the screen can carry class-level structure that is invisible in a joint view yet still biases per-compound calls.

## Conclusion

v2's pathway scoring is bounded by the 978-landmark panel (~16% of each pathway), its critical-pathway filter is strict rather than redundant, and its drug-embedding input hides a real pharmacologic class signal (strong in the HG branch, diluted by the 512 UniMol-token dimensions in the joint input the model consumes). All three are properties of the screen's *inputs and aggregation*, independent of the XPert finetune-data change that defines v2. They bound how much biological claim the ranking can support.
