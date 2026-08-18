# DToxS-LINCS testbed workflow (xpert-fda-miri-v1.1)

`workflow/xpert-fda-miri-v1/` is the canonical finalized 11-step screen, and
its frozen DToxS finetune data (`input/dtoxs.h5ad`, 521 samples / 36 drugs) was
built outside the repo from the now-deleted well-level deposit. Rather than
changing v1 (whose benchmark/report numbers are the established baseline), we
created a **testbed sibling** `workflow/xpert-fda-miri-v1.1/` that mirrors v1's
pipeline on the 12-step layout (same as v2) and rebuilds the DToxS finetune
data from the re-homed LINCS LDS-1593 wells (step `03_cm_finetune_data`).

The fresh LINCS build yields **662 samples / 45 drugs** — a legitimate superset
of v1's 521/36 (the raw GEO deposit lacks the non-KI perturbations) with 5
compounds re-indexed by the author canonical-SMILES rule (cyclosporin 73→72,
crizotinib 4786→1176, estradiol 5674→1291, paroxetine 6743→6737, ruxolitinib
6798→6797). Because sample-for-sample identity with v1's frozen data is not
reproducible (and not desirable), the **parity gate** requires **schema parity**
(the rebuilt h5ad is structure-identical to v1's `input/dtoxs.h5ad` in the
fields the XPert loader consumes) and emits an **auditable delta report**
(`03_cm_finetune_data/output/parity_delta.md`). v1 stays frozen; v1.1 tests
the fresh build + renumbering, and its outputs are for comparison only.
