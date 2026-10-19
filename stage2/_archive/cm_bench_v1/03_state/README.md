# 03_state — State adapter (ticket 11)

**Status: implemented, scored.** Entrant `state` in the cold-drug and cold-dose
leaderboards. Runs in the pinned env `arc-stage2-v1-state` (ADR 0021) over the
vendored repo `pert_model/state/` (pinned commit `9bbfe78a43`, `.git` removed).

- `entry.py` — runs inside the State env: loads the ST-HVG-Tahoe checkpoint in
  `pert_model/state/weights/ST-HVG-Tahoe/`, optionally finetunes, predicts each
  test well, and writes `pred.npz` = Δ aligned to the common gene universe.
- `adapter.py` — `SubprocessAdapter("state", "arc-stage2-v1-state", entry.py)`.
- `hvg_gene_names_recovered.json` — the 1483/2000 recovered HVG identities.
- registered in `_lib/model_registry.py` as `MODEL_ADAPTERS["state"]`.

## How it works

The checkpoint (ST-HVG-Tahoe) is a State transition model trained on Tahoe-100M:
input = 2000 HVG expression per cell + a one-hot over 1138 Tahoe drug-dose
combos; output = predicted perturbed HVG expression. The entry:

1. maps the benchmark baseline (log2-CPM over the 20570-gene universe) into the
   model's 2000-dim HVG space (unmapped/out-of-universe dims = 0);
2. maps each CM drug to the closest Tahoe dose one-hot for the same base name;
3. optionally finetunes (domain adaptation on the fold's control baselines, the
   checkpoint's own test-time-finetune path);
4. predicts, takes Δ = pred − baseline, and scatters back to the universe.

## Documented approximations (read before interpreting the score)

- **HVG identities unpublished.** The public ST-HVG-Tahoe artifacts do not name
  the 2000 HVGs. They were recovered by exact value-matching against the
  `arcinstitute/State-Tahoe-Filtered` training shards: two independent runs
  (8000-row direct hash and full 30394-row block hash) agree on every recovered
  column (**zero conflicts**), covering **1483/2000** dims, **1398 of which are
  in the benchmark universe**. Recovery, not fabrication.
- **Input scale mismatch.** State trained on `log1p(raw counts)`; the benchmark
  only froze `log2(CPM+1)`. The benchmark baseline is fed as-is.
- **One-hot perturbation / cold-drug limitation.** Only **17 of 79** CM drugs
  appear in Tahoe (Afatinib, Bortezomib, Crizotinib, Dabrafenib, Erlotinib,
  Furosemide, Gefitinib, Ponatinib, Regorafenib, Rosiglitazone, Sildenafil,
  Sunitinib, Tofacitinib, Trametinib, Vandetanib, Vemurafenib, Verapamil). The
  other 62 fall back to the DMSO control one-hot, so State predicts ≈no-change
  for them. The model also has no cardiomyocyte cell type (50 cancer lines).
  The resulting cold-drug macro AUPRC (0.0865) is therefore ≈ the no-change
  baseline (0.0844) — reported as-is, never padded with fabricated rows.
