# 0084. Retire Stage 2 Bases and reserve a loss comparison step

**Status:** accepted (2026-11-23), implementation pending.

ADR 0083 retained Bases to bridge active consumers and support historical replay. Retire that replay and remove the Bases step after active consumers have migrated. Keep accepted study results and their audit records. The corpus owns source counts and readers, the gene-axis step owns axis evidence and projection, the view step owns expression views and their readers, and util owns shared helpers, including derivation kernels still used by external DEG builders.

Name Step 02 `02_gene_axis` and Step 04 `04_model_view`. Step 04 keeps its existing tenv single-dose comparison. Reserve `05_tsg_loss` with routing and a placeholder README. Shift the drug representation, dose, width and dose, and transfer studies to Steps 06–09, respectively. Reuse public util APIs when they match the operation and keep study-specific fold validation, aggregation, and tables in the owning study.

## See

- [ADR 0083](0083-stage2-corpus-axis-view-boundaries.md)
- [ADR 0062](0062-stage2-accepted-products-and-reproducible-figures.md)
