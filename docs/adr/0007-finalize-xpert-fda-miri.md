# Finalize the MIRI screening workflow as xpert-fda-miri-v1

The three overlapping MIRI screening workflows (`xpert-fda-v1` FDA extension via HG finetune + norm mapping, `xpert-fda-v2` FDA extension via author-pinned inference, `xpert-fda-v3` author-aligned L1000-contained screen) are consolidated into exactly one canonical workflow: `xpert-fda-v3` becomes `stage1/xpert-fda-miri-v1`, `xpert-fda-v1` and `xpert-fda-v2` are removed, and the legacy `xpert-fda-v0` is kept frozen as the sole historical snapshot (`stage1/xpert-fda-miri-v0`).

## Context

Three near-identical workflows coexisted in the old repositioning tree (renamed from `miri_reposition/` to `stage1/` here), so the canonical screening approach was ambiguous. v3 was the author-aligned, L1000-contained screen: it screens only the 1,004 FDA drugs whose exact canonical SMILES matches the L1000 embedding table (drug_idx < 8981), uses the author-reference HG table byte-identical (no HG training / FDA-extension inference), and owns the authoritative Entrez-verified DToxS→L1000 gene map. It is the cleanest author alignment: every screened drug embedding is exactly what the author's XPert was pretrained on.

## Decision

- Finalize on the v3 basis. `xpert-fda-v3` is renamed to `stage1/xpert-fda-miri-v1` and reorganized into the finalized 11-step layout: `01_drug_embed`, `02_cm_baseline`, `03_xpert_pretrain`, `04_xpert_finetune_bench`, `05_xpert_finetune`, `06_xpert_predict`, `07_miri_phenotype`, `08_miri_efficacy_scores`, `09_filter`, `10_integrate`, `11_report` (benchmark-before-finetune; perturbation prediction + z-score merged; new agentic bilingual report step).
- The workspace root is renamed `miri_reposition/` → `stage1/`. `xpert-fda-v1` and `xpert-fda-v2` are removed (git rm); `xpert-fda-v1_history` is deleted; `xpert-fda-v0` → `stage1/xpert-fda-miri-v0` stays gitignored and frozen.
- The finalized workflow uses only XPert-author-given drug embeddings (author HG reference byte-identical + precomputed UniMol array indexed by drug_idx) and the author pretrained XPert weight — no HG training, no FDA-extension inference, no feature regeneration from SMILES. Scope stays 1,004 L1000-contained FDA drugs.
- ADR 0006 (HG-aligned FDA extension) is superseded: the FDA-extension machinery is archived (see ADR 0008) and no canonical workflow consumes it.

## Consequences

- Exactly one canonical MIRI workflow (`stage1/xpert-fda-miri-v1`) and one frozen snapshot (`stage1/xpert-fda-miri-v0`).
- Deterministic workflow steps reproduce v3's archived outputs (parity verified; GPU steps 04/05/06 are not re-run).
- The 11_report step produces the bilingual final report over the `09_filter` survivors with integrated-score percentile > 0.5 (41 drugs), driven by a checked-in `prompt.md` and resumable per-drug JSON cache.
