# Ticket 03 — Author-reproduction gate: result record

**Status: documented blocker → gate A repurposed (see docs/adr/0006).**

## What was tested

Whether a frozen, layer-by-layer forward pass of the author-given HG model over
the original L1000 graph (8,981 drugs, 19,392 genes, author node features,
original PPI/DTI/DDS edges) reproduces the author reference embeddings
(`drug_emb_model/hg/embeddings/article/HG_drug_embeddings.npy`) exactly.

## Evidence

| Forward output (frozen `model/20260518_102121/global_fold_early_stop.pth`) | max\|diff\| vs author ref | stats |
|---|---|---|
| Original-graph forward (this experiment) | 10844.76 | mean 6.70, std 51.37 |
| Local phase1 reference (`model/.../phase1_drug_embeddings.npy`) | 10844.76 | mean 6.70, std 51.37 |
| Author reference | — | mean 22.61, std 226.31, 97.8% zeros |

Every available checkpoint (pretrain + all finetune runs) produces mean ≈6.7–8.2,
std ≈51–56, max ≈1500 — statistically incompatible with the author reference
(mean 22.6, std 226, max 10844). The article's figshare release
(`drug_emb_model/hg/figshare/`) ships only XPert checkpoints plus the HG
*embeddings* and input data — **no HG checkpoint**. The author HG model that
produced the reference was never released; the local checkpoint is a
reproduction of it.

The frozen checkpoint IS deterministic and re-usable: it reproduces its own
phase1 reference within float tolerance (maxdiff 7.32e-4) and two consecutive
runs agree to 2.4e-4.

## Conclusion

Gate A as specified ("frozen forward reproduces the author reference, max|diff|
== 0") **cannot pass**: no checkpoint or released artifact produces the author
reference. Per the feature spec this would block any extension inference. After
consultation, the HG-aligned feature proceeds as a **hybrid** (docs/adr/0006):
rows 0–8980 = author reference byte-identical by construction, FDA rows =
author-pinned inference in the reproduction model's space + full-CDF quantile
alignment to the author-embedding scale. Gate A is repurposed to assert the
checkpoint is deterministic and re-usable for pinning (reproduces phase1), which
passes. All four gates (A–D) pass in `_lib/verify_gates.py`.
