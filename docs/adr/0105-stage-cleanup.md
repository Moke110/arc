# Retain the Stage 1 v2 screen and consolidate Stage 2 model experiments

Accepted 2026-11-27. Supersedes the local archive requirement in ADR 0087.

## Decision

Complete the following blocks in order. Start the next block only after verifying the retained formal products of the preceding block. Implement and validate each ticket in a worktree based on dev. Record inventories and acceptance evidence outside disposable scratch before cleanup.

1. Promote the complete v2 workflow to the Stage 1 root. Retain its two predecessor H5AD inputs under Stage 1 and preserve the former v3 audit CSV used by Annotation under its declared local additions. Verify digests, external weight links, schema checks, official folds, downstream inputs and reproducible outputs with retired siblings unavailable. Remove the other local screens, including the unique v0 reports.
2. Consolidate model training and evaluation for active Steps 04–10 under Stage 2 util/training and util/eval. Expose architecture, objective terms and coefficients, chemical representation, dose transform and conditioning, ordered gene axis, split, fold and seed as configuration. Preserve study-specific evaluation populations, weighting, axis projections, calibration and random state. Keep numbered execution entry points as thin callers. Steps 01–03 own data construction. Remove the retired archive route and delete history only after proving all active products work independently.
3. Following integration into dev, deploy the accepted tree to the materialized checkout through main. Clear repository tmp, then remove current-user ARC/CarPS scratch in WSL temporary directories. Inspect registered worktrees, active processes, ignored payloads and external targets before removal. Preserve unrelated user work. Compare before/after free space and verify the installed workflows with retired and temporary inputs absent. Publish main after verification.

## Reproduction

Inventory every retained formal product, its inputs, command, environment and verification method. Compare checksums for preserved bytes, verify published metrics and figures, and exercise every supported experiment configuration. Compare representative sealed checkpoints and short training runs with their original predictions, objective and metrics. Use pinned Git code for historic runtime digest checks. Preserve original evidence. Record failures and resolve them before moving to the next block. Update affected README and AGENTS routes, inventory, ignores and health checks at each cutover.

## Rationale

The local archive and scratch occupy most of the reclaimable space. Deleting them loses ignored historical bytes. Freezing the small predecessor inputs and verifying the retained studies before deletion preserves their reproducibility within the available disk budget.
