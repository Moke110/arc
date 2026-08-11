# HG-aligned FDA extension via author-pinned frozen inference

> **Status: SUPERSEDED** by [ADR 0007](0007-finalize-xpert-fda-miri.md). The
> FDA-extension workflow variants (v1 HG finetune, v2 author-pinned inference)
> are removed; the canonical screen is the L1000-contained
> `workflow/xpert-fda-miri-v1`. The FDA-extension machinery and its data/weights
> are preserved for provenance in the gitignored `reference/hg-fda-extension-archive/`.

The `xpert-fda-v2` workflow extends the HG drug table to the 2,123 FDA-extension drugs (rows 8981–11103) using **author-pinned frozen inference** instead of the v1 low-LR finetune, keeping rows 0–8980 byte-identical to the author reference by construction.

## Context

The author-given HG *embeddings* (now workflow-owned at `workflow/xpert-fda-miri-v1/input/HG_drug_embeddings.npy`, byte-identical to the article's figshare release) are the L1000 rows the author-given XPert was trained on. v1 produced its FDA rows by low-LR finetuning the HG model on the extended graph and patching scale with an undocumented quantile-interpolation norm-mapping plus a cosine-guard — a trust gap on every axis that defines an embedding.

The intent of the HG-aligned feature was: run the *same frozen model* that produced the L1000 embeddings over the extended graph, pinning the L1000/gene states per layer so only FDA drug states evolve. **That model does not exist.** The article's figshare release ships only XPert checkpoints, the HG embeddings, and the HG input data — no HG checkpoint. The available `model/20260518_102121/global_fold_early_stop.pth` is a local reproduction: a frozen forward over the original graph reproduces its own reference (`phase1_drug_embeddings.npy`, maxdiff 7.3e-4) but NOT the author reference (maxdiff 10844, mean/std 6.7/51 vs the author's 22.6/226).

## Decision

Build the HG-aligned FDA extension as a **hybrid** (the only path with real artifacts):

- **L1000 rows (0–8980)** are the author reference, byte-identical **by construction** (copied into the table, verified by gate B) — preserving v1 parity and the DToxS benchmark.
- **FDA rows (8981–11103)** are produced by author-pinned frozen inference using the local reproduction as the frozen model: capture per-layer L1000-drug and gene states on the original graph (the pins), then run the frozen model layer-by-layer over the extended graph overwriting the L1000/gene states with the pins so only the FDA drug states evolve.
- FDA node features use the UniMol mol-token read from the frozen workflow UniMol array (ticket 02 verdict: fresh UniMol generation does not reproduce the L1000 node features, cos 0.16; the frozen array's mol-token is cos 0.97 close). FDA DTI edges are the existing ChEMBL pChEMBL≥6.0 edges; FDA DDS edges are recomputed with the pretraining recipe `Tanimoto > 0.5` (v2-owned artifact).
- The FDA states live in the reproduction model's space; they are **full-CDF quantile-aligned** to the author-embedding distribution (mean/std/density), which is the documented, principled replacement for v1's norm-mapping hack. No v1-style L1000 norm-mapping; the L1000 rows are never rescaled.

**Gate repurposing (recorded deviation):** the spec's gate A ("frozen forward reproduces the author reference, max|diff| == 0") is impossible because the author model was never released. Gate A now asserts the checkpoint is deterministic and re-usable for pinning: a frozen forward over the original graph reproduces the checkpoint's own phase1 reference within float tolerance (maxdiff 7.3e-4) and two runs agree. This is documented in the feature's tickets and in this ADR; gates B (L1000 rows byte-identical), C (FDA/train std ratio ≈ 1.0), and D (FDA DTI signal ≥ 2× random) remain as specified.

## Alternatives considered and rejected

- **Low-LR finetune on the extended graph (v1)** — rejected: weight updates, cosine guard, and norm-mapping that the alignment is meant to remove; not reproducible from a frozen artifact.
- **Plain frozen forward over the extended graph** — rejected: it perturbs the L1000/gene states, so rows 0–8980 would drift from the author reference and invalidate the pretrained XPert.
- **Explicit restore step** — rejected: restoring the L1000 rows post-hoc reintroduces the "silent drift" failure mode the pinning design exists to avoid.
- **Treating the reproduction checkpoint as *the* reference (phase1 rows)** — rejected: rows 0–8980 would no longer be the author reference, breaking v1 parity and the DToxS benchmark gate.
- **Halting the feature** — rejected after the hybrid was chosen; the hybrid preserves the parity-critical guarantees and produces FDA rows from a frozen, deterministic model.

## Consequences

- `drug_emb_model/hg` is untouched (the model library stays an author-given artifact; the workflow owns the pinned-inference script and the resulting table).
- v2 `input/` is a symlink to v1's `input/`; one source of truth for the frozen UniMol array, SMILES mapping, DToxS, and FDA list.
- The residual mismatch to accept and document: training-DTI edge provenance is unknown (FDA DTI matches by source choice, not identity); the reproduction model's space differs from the author model's (no author model exists).
