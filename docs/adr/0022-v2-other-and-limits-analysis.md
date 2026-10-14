# v2 non-numbered `_other` stage + `13_limits` limitations analysis

**Status: accepted**

Two linked changes to `workflow/xpert-fda-miri-v2/`, both about *what belongs at
the top level of a workflow*.

**`13_other/` is renamed to `_other/`.** The directory holds ad-hoc analyses that
sit outside the numbered screen (the GO/proliferation study and the top-20
structure gallery); it is not a screen *step*. A leading underscore marks it as
non-numbered workspace material, name-adjacent to `_lib/`, mirroring v1's
archived, non-numbered `report/` stage. Only the directory name and the strings
that referred to it change (script docstrings/usage lines, the generated
`top20_structure.html`, and the workflow `AGENTS.md` entry) — no analysis
behaviour changes.

**`13_limits/` is added as a numbered *analysis* step.** It is a read-only
diagnostic that quantifies three structural limitations of the v2 screen, one
figure each, and writes a deterministic `report.md`:

- `landmark_coverage.png` — for each of the 40 scoring pathways, total genes
  vs the genes actually measured by the 978 L1000-landmark panel.
- `critical_pathway_overlap.png` — pairwise gene overlap of the 8 curated
  critical pathways, row-normalized (`|i∩j| / |i|`, asymmetric), annotated with
  raw counts.
- `drug_embedding_pca.png` — two-panel PCA of every drug passing the step-10
  filter: left = the model's own drug input
  `[HG(256) | mean-pooled projected UniMol tokens(256)]`, right = the HG
  component alone, with the 10 direct α/β-adrenergic ligands highlighted.

Decisions within the step:

- **The step is read-only and consumes existing artifacts** (step 06 finetuned
  checkpoint's `drug_emb.linear`, step 09 `scores.csv`, step 10
  `path_filtered.csv`, plus the frozen inputs and `input/MIRI_pathways.json`).
  It changes no screen artifact and is not a screening step.
- **Figures sit at the step top level, not under `output/`.** `.gitignore`
  ignores every `output/` directory, so deliverables placed there are not
  tracked; this follows the existing convention (`10_filter/pyramid-filters.png`,
  `11_integrate/heatmap_rank.png`, …). CSVs (`landmark_coverage.csv`,
  `critical_pathway_overlap.csv`) do go under `output/`.
- **The adrenergic classification is a pharmacological one**: the direct
  α/β-receptor ligands among the filtered drugs (agonists DOBUTAMINE, MIDODRINE,
  METHOXAMINE, MEPHENTERMINE; antagonists PHENOXYBENZAMINE, PROPRANOLOL,
  SOTALOL, DILEVALOL, LABETALOL, MOXISYLYTE), hard-coded in `limits.py` with a
  comment. Indirect sympathomimetics, reuptake inhibitors, MAOIs and
  α2-antagonist antidepressants are deliberately excluded. Same spirit as the
  adrenergic group in `_other/top20_structures.py`.
- **The drug embedding for PCA** reproduces the tensor the finetuned transformer
  sees — the `drug_emb.linear` 512→256 projection with token 0 overwritten by HG
  — then mean-pools the valid UniMol tokens. Position embeddings and LayerNorm
  are deliberately omitted: they are not drug-identity signal. The vector is
  per-dimension z-scored before PCA so the large-magnitude HG block and the
  small UniMol block are not silently weighted by scale.

**Finding (recorded because it reframes limitation 3).** In the joint input the
model actually consumes, the adrenergic ligands are *not* clustered
(adrenergic-vs-rest Cohen's d on PC1 ≈ 0.34; agonists and antagonists
intermixed). The class signal is real but lives in the HG branch alone (HG-only
d ≈ 1.37) and is diluted by the 512 UniMol-token dimensions in the joint input
(UniMol-only d ≈ 0.33). The two-panel figure is kept precisely so the negative
result and its cause are both visible; `report.md` states it as a negative
result rather than a cluster claim.

**Consequences:** the workflow `AGENTS.md` gains a `13_limits/` entry and the
`_other/` line; `CONTEXT.md` gains the *Limitations analysis* and
*Non-numbered stage* terms. `build_screen_notebook.py` is not touched (step 13
is a diagnostic, not part of the screen walkthrough). The three figures and
`report.md` are tracked; the step's CSVs are under the gitignored `output/`.
