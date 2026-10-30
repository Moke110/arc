# Stage 2 setup

Nine deployed models, one run card each. The goal is to prove the models run in
their own environments, not to compare them.

- **HTML report for non-computational readers:** [`model_architecture_weights_benchmarks.html`](model_architecture_weights_benchmarks.html) —
  architecture (drug encoder vs. main predictor), weights, and published benchmark claims per model.

| model | status | run | artifact |
|---|---|---|---|
| PRnet | ✅ passed (978 demo, setup-only) | `bash stage2/01_setup/prnet/run.sh` | `prnet/output/demo_split_y_pre_array.csv` (661 × 978) |
| CPA | ✅ passed (Kang, IFN-β) | `bash stage2/01_setup/cpa/run.sh` | `cpa/output/cpa_pred.npy` (2,000 × 5,000) |
| chemCPA | ⚠️ load-only forward | `bash stage2/01_setup/chemcpa/run.sh` | `chemcpa/output/pred.npy` (8 × 977) |
| biolord | ✅ passed (3-epoch train) | `bash stage2/01_setup/biolord/run.sh` | `biolord/output/pred_mean.npy` (640 × 4,999) |
| PerturbNet | ✅ checkpoint + latent forward | `bash stage2/01_setup/perturbnet/run.sh` | `perturbnet/output/meta.json`, `pred_latent_*.npy` (2 × 10) |
| scVIDR | ✅ 2,000-gene train + dose predict | `bash stage2/01_setup/scvidr/run.sh` | `scvidr/output/pred_doses.npz` (1/10/100 µM) |
| Squidiff | ✅ 2,000-gene forward | `bash stage2/01_setup/squidiff/run.sh` | `squidiff/output/pred.npy` (2 × 2,000) |
| CMonge | ✅ official demo train + eval | `bash stage2/01_setup/cmonge/run.sh` | `cmonge/output/meta.json`, `demo_logs.yml` |
| DECANT | ✅ 2,000-gene forward, two doses | `bash stage2/01_setup/decant/run.sh` | `decant/output/pred_gene_dose{1,100}.npy` (4 × 2,000) |

Excluded from setup:

- **preprints**: X-Pert, DORA, BaiZe, PerturbLDM — no Stage 2 vendor clone or
  run card; `pert_model/XPert/` is the Stage 1 XPert library, not a Stage 2 entrant;
- **gate-fail references**: WAVE (no dose), ExPO (978-only).

Common conventions:

- environments: one per model (see each run card); `arc` itself stays clean;
- outputs: `<model>/output/` (gitignored);
- logs: `<model>/output/run.log`;
- no upstream model code is edited; adapters live in the run-card directory.

**Key finding:** the public chemCPA weights are L1000 **977/2,000-gene** `ComPert`
state dicts (188 drugs, 3 covariates), not 5,000/20,570-gene models. The chemCPA
run card therefore uses the load-only fallback; real-data inference on the
downloaded 5,000-gene `combo_sciplex` data is blocked without retraining.

PRnet's 978-gene LINCS demo is **setup-only** and is quarantined from any
Stage 2 evaluation.
