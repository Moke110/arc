# AGENTS.md

The **HG-aligned MIRI screening workflow** (`miri_reposition/xpert-fda-v2/`): the same eleven numbered steps as the active `xpert-fda-v1`, but the FDA-extension HG embeddings (rows 8981–11103 of the 11,104-drug table) are produced by **author-pinned frozen inference** instead of low-LR finetune. See `docs/adr/0006` for the decision and the domain glossary in the root `AGENTS.md` for the vocabulary.

## Why v2 exists

v1 produces the FDA-extension HG rows by finetuning the HG model on the extended graph and patching scale with an undocumented norm-mapping. v2 embeds FDA drugs with a frozen, deterministic model over the extended graph while keeping the L1000 side completely frozen, then aligns the FDA rows to the author-embedding distribution. The author HG model was never released (figshare ships only XPert checkpoints + the HG embedding); the local checkpoint `drug_emb_model/hg/model/20260518_102121/` reproduces its own phase1 reference, not the author reference — gate A is repurposed accordingly.

## Structure

- Same 11 step directories (`01_drug_list/`…`11_integrate/`) and `_lib/` helpers as v1; steps 01 and 03–11 run unchanged.
- `input/` is a **symlink to v1's `input/`** — one source of truth for the frozen UniMol array, SMILES mapping, DToxS, and FDA drug list. Never store workflow-owned artifacts there.
- `02_drug_embeddings/` — v2's HG step, workflow-owned:
  - `02_drug_embeddings.py` — idempotent entry (reuses an intact table that passes gates; `--force` regenerates; `--skip-gates` to override a reviewed C/D borderline).
  - `hg_aligned_infer.py` — author-pinned frozen inference + full-CDF quantile alignment (drives the model's layer modules; the library is untouched).
  - `compute_dds_fda_tanimoto.py` — v2-owned FDA DDS edges with the pretraining recipe (`Tanimoto > 0.5`).
  - `author_reproduction_gate.md` — ticket-03 evidence record (documented blocker).
  - `node_feature_recipe_verdict.json` — ticket-02 verdict (fresh UniMol cls_repr does NOT reproduce the L1000 node features, cos 0.16; the frozen array's mol-token is cos 0.97 close).
- `_lib/verify_gates.py` — the step-02 verification seam (gates A–D).
- `_lib/verify_parity_v1.py` — toolchain-parity seam (gate E): steps 01/03/08 outputs + the DToxS benchmark reproduce v1.
- `_lib/compare_v1.py` — regenerates `11_integrate/divergence_vs_v1.md` (gate F, observational).

## Key invariants (gates)

- **A** — the frozen checkpoint is deterministic and re-usable for pinning (reproduces its own phase1 reference within float tolerance).
- **B** — table rows 0–8980 are byte-identical to the author reference (`embeddings/article/HG_drug_embeddings.npy`).
- **C** — FDA/train std ratio ≈ 1.0 (full-CDF quantile alignment; no v1-style norm-mapping of L1000 rows).
- **D** — FDA DTI drug–gene cosine vs random ratio ≥ 2× on the aligned rows (currently 3.53×).
- Any gate failure aborts (exit 1); a researcher may override a reviewed borderline manually via `--skip-gates`.

## Entry points

```bash
# Regenerate the HG table + run gates (idempotent; --force to regenerate)
python 02_drug_embeddings/02_drug_embeddings.py
# Gate harness alone
python _lib/verify_gates.py
# Toolchain parity vs v1 (gate E) and divergence report (gate F)
python _lib/verify_parity_v1.py
python _lib/compare_v1.py
```

## Gotchas

- Steps 04/05 re-run the XPert DToxS finetune/benchmark on the aligned table — GPU-fp-nondeterministic across environments (parity is asserted at the benchmark level, not bit-level).
- Predict (step 06) indexes the precomputed UniMol array by drug index — never regenerate UniMol features from SMILES.
- Large artifacts (`.npy`, `.pth`, `.csv`, `.h5ad`) and `.env` are gitignored — never commit them.
