# Name the public stages CPST and CPSP

Accepted 2026-11-27. This updates the naming decision in ADR 0099 and the public
release boundary in ADR 0105.

Use `cpst/` for Stage 2, CarPS-Transcriptome, and `cpsp/` for Stage 3,
CarPS-Phenotype. The selected Stage 2 model and its fixed alpha release also
take the CarPS-Transcriptome name. Use the new names in current paths, public
imports, commands, documentation and machine-readable identities.

Reissue current registrations, receipts and checksum chains with the new
identities. Verify them against originals from pinned Git revisions or sealed
copies using a recorded mapping from original to migrated evidence. Preserve
the scientific protocols and the contents of unchanged data, checkpoints,
predictions and figures. Check the migrated entry points and replay on the
materialized checkout before publishing.

Stage 1 is a local predecessor. Keep its working tree and documentation local,
and exclude it from new public trees and required public checks.

Rewrite the existing README and AGENTS files in CPST and CPSP from the deepest
directories upward. Rewrite the repository AGENTS and create its README last.
Give each document to a separate agent that inspects its owning directory before
writing. Readmes summarize completed work and downstream selections in step
order. Agent instructions own execution and verification guidance, with the
common training, evaluation and plotting APIs documented at their owner.
