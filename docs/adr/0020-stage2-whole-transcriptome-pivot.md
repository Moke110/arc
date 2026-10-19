# Stage 2 pivots to whole-transcriptome prediction (no 978-landmark subset)

**Status: accepted**

Stage 1 built the MIRI drug-repositioning method entirely on the 978 L1000 landmark genes:
XPert is pretrained and finetuned on 978-dim perturbomes, and every Step 05/08/09 evaluation
and score consumes that 978 space. But the raw CM datasets are full-transcriptome (`~12k` /
38k / 44k genes), and the landmark projection (together with XPert's training) is the root
of the known bias — ~97.2% down-regulated DEGs, low inter-drug variance, poor drug
discrimination.

**Decision.** Stage 2 benchmarks **whole-transcriptome** drug-perturbation models on
cardiomyocytes, and uses **no 978-landmark subset anywhere**; Stage 1's evaluation is a
*method reference only*. The two incumbent models in `pert_model/` (XPert, TranSiGen) are
landmark-only and are therefore excluded as Stage 2 entrants.

**Consequences.** Stage 1 stays the historical landmark deliverable, unchanged. The MIRI
scorer (`08/09`) is 978-based, so producing a whole-transcriptome MIRI drug ranking requires
rebuilding the MIRI reference gene sets and scorer on the whole transcriptome — deferred to
a follow-on stage informed by the Stage 2 winner. The benchmark's ground truth and metric
are computed on the full-transcriptome common gene universe (intersection of the datasets).
See `docs/roadmap.md` (Stage 2) and the archived spec at
`stage2/_archive/cm_bench_v1/`. The bring-up-stage 978 exception for PRnet is
recorded in ADR 0024; it does not change this evaluation rule.
