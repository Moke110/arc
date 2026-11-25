# Architecture Decision Records

Index of the architecture decision records in `docs/adr/`.

**Convention.** One decision per file, named `NNNN-<slug>.md`. Each record opens with a
`**Status:**` line (`accepted (YYYY-MM-DD)`, `superseded by NNNN`, …) and closes with a
`See` section of repo-relative pointers. Records are historical once accepted: later
decisions amend them, they are not silently rewritten. Read a record together with the
files its `See` section names.

**Numbering note: 0023 and 0050 are gaps, not missing files.** Neither number was ever
assigned (the sequence runs 0022 → 0024 and 0049 → 0051 deliberately); there is nothing to
restore.

The index covers every record in the directory, 0001–0077 (75 files; 0023 and 0050 are the
gaps). Two early records are marked superseded in place, and 0026 is marked
amended by 0033. Stage 3's initial scope in 0034 is refined by 0035–0037.
Its formal workspace/corpus design subsequently replaces the specified wayfinder
details in 0038–0041. ADR 0043 reopens the corpus as source measurements without endpoint gates; ADR 0044 subsequently replaces universal-corpus-first planning with source-first corpora.
ADR 0077 implements the paired 5k-only axis-local expression release (#243–#245); current verification is in the Bases axis-local acceptance receipt.
ADR 0073 remains the disease/four-axis scientific release: full-source collection, tier-first DEG1000/f1000 and four axis sizes and historical full-axis views (#220–#225). ADR 0072 historically amended the disease eligibility
and TOX roster of ADR 0071; its measured release is superseded by 0073. ADR 0071 is the base gene-axis contract and supersedes ADR 0066's *Source eligibility*,
*Priority* and *Candidate table and the two-level rank-sum fusion*, and supersedes ADR
0069/0070 (its `AF` removal is reversed). Each of those four records carries a status-line
pointer to 0071; their bodies are unchanged, per the convention above.

| # | Decision |
|---|---|
| [0001](0001-category-taxonomy.md) | Category taxonomy for the repository root |
| [0002](0002-components-as-libraries.md) | Components are libraries; workflows generate their own artifacts |
| [0003](0003-conda-env-arc.md) | Conda environment renamed from xpert to arc |
| [0004](0004-freeze-and-rebuild-fda-miri-vs.md) | fda_miri_vs frozen as xpert-fda-v1_history and rebuilt — *superseded by 0007* |
| [0005](0005-docs-version-controlled.md) | docs/ and domain docs are version-controlled |
| [0006](0006-hg-aligned-fda-extension.md) | HG-aligned FDA extension via author-pinned frozen inference — *superseded by 0007* |
| [0007](0007-finalize-xpert-fda-miri.md) | Finalize the MIRI screening workflow as xpert-fda-miri-v1 |
| [0008](0008-components-as-pure-static-libraries.md) | Components as pure static libraries (extends ADR 0002) |
| [0009](0009-xpert-fda-miri-v2-union-cm.md) | Union-CM finetune: xpert-fda-miri-v2 |
| [0010](0010-align-union-finetune-to-author-semantics.md) | Align the union-CM finetune data to the author L1000 semantics (xpert-fda-miri-v2) |
| [0011](0011-flatten-pert-data-rehome-lincs-lds1593.md) | Data-lake restructure: flatten pert_data, re-home LINCS LDS-1593, delete dtoxs-xpert |
| [0012](0012-xpert-fda-miri-v1.1-dtoxs-lincs-testbed.md) | DToxS-LINCS testbed workflow (xpert-fda-miri-v1.1) |
| [0013](0013-release-lincs-lds1593-freeze-v1-dtoxs.md) | dtoxs_primary: delete raw GEO deposit, release LINCS LDS-1593 flat, freeze v1 data |
| [0014](0014-unified-benchmark-folds.md) | Unified benchmark fold scheme across v1 / v1.1 / v2 |
| [0015](0015-align-v2-to-v11-primary-ipsc.md) | Align v2 finetune to v1.1 — primary + iPSC-CM, real-dose bins, imported folds |
| [0016](0016-per-plate-normalization.md) | Per-plate matched-control and QN grouping in the rebuilt finetune data |
| [0017](0017-xpert-fda-miri-v3.md) | xpert-fda-miri-v3 — all-CM finetune sibling (per-plate normalized) |
| [0018](0018-dtoxs-primary-bench-and-data-ablation.md) | Step-05 benchmark restructure: dtoxs-primary bench + data-ablation experiment |
| [0019](0019-v2-step05-full-and-data-ablation.md) | v2 step-05 restructure: full/ bench + drug-partitioned data-ablation |
| [0020](0020-stage2-whole-transcriptome-pivot.md) | Stage 2 pivots to whole-transcriptome prediction (no 978-landmark subset) |
| [0021](0021-per-model-environments-stage2.md) | Stage 2 models run in per-model conda environments |
| [0022](0022-v2-other-and-limits-analysis.md) | v2 non-numbered `_other` stage + `13_limits` limitations analysis |
| [0024](0024-stage2-bringup-978-quarantine.md) | Stage 2 bring-up quarantines the 978-landmark PRnet demo |
| [0025](0025-stage1-stage2-directory-split.md) | Split stage1 and stage2, rename setup, and exclude preprint models |
| [0026](0026-stage2-cm-data-contract.md) | Stage 2 CM data — frozen corpus and filtered views — *amended by 0033* |
| [0027](0027-stage2-benchmark-contract.md) | Stage 2 whole-transcriptome benchmark contract |
| [0028](0028-stage2-5k-hvg-benchmark.md) | Stage 2 5,000-gene HVG sensitivity benchmark |
| [0029](0029-tahoe-pseudobulk-pretrain-corpus.md) | Tahoe-100M enters Stage 2 as a pseudobulk pretrain corpus |
| [0030](0030-external-data-provenance-and-step-layout.md) | External data provenance and the Stage 2 step layout |
| [0031](0031-tahoe-cm-schema-view-and-prnet-pretrain.md) | A CM-schema Tahoe pretrain view and the PRnet pretrain bring-up |
| [0032](0032-prnet-formal-transfer.md) | PRnet formal transfer: paired Tahoe-pretrained vs random-init CM study |
| [0033](0033-stage2-dataset-corpus.md) | Stage 2 dataset corpus — native-axis corpora, views downstream, 04 folded into 02 |
| [0034](0034-stage3-phenotype-data-boundary.md) | Stage 3 starts with drug–phenotype data; predicted transcriptomes remain outside the first data-planning effort |
| [0035](0035-stage3-stratified-phenotype-curation.md) | Stage 3 curates separate experimental strata and distinguishes retention, qualified counting, and label readiness |
| [0036](0036-stage3-continuous-evidence-before-direction-labels.md) | Stage 3 retains continuous evidence and unresolved directions instead of manufacturing three-class labels |
| [0037](0037-stage3-phenotype-release-layout.md) | Stage 3 uses existing provenance items and typed CSV/JSON releases with separate qualification and labeling gates |
| [0038](0038-dedicated-phenotype-data-and-stage3-corpus.md) | Dedicated pheno_data category and Stage 3 01_data_corpus workspace boundary |
| [0039](0039-condition-corpus-and-model-view-separation.md) | Fine condition corpus, shared MTT/CCK8 phenotype, endpoint-inclusive sample identity, retained readouts and later model views |
| [0040](0040-numeric-endpoints-and-source-readout-sequences.md) | Source-inclusive sample keys, retained readout sequences and numeric magnitude-plus-significance endpoints |
| [0041](0041-formal-corpus-statistics-export-and-coverage.md) | Statistical comparison/BH families, Parquet corpus files and coverage without per-file 500 filtering |
| [0042](0042-tenv-model-benchmark-redesign.md) | Tenv model comparison with shared Morgan inputs, trainable adapters and three-fold three-seed evaluation |
| [0043](0043-stage3-source-measurement-corpus.md) | Source-measurement samples without significance/endpoint admission gates, revised key and source annotations |
| [0044](0044-stage3-source-first-corpora.md) | Select the best source per condition/phenotype, build its corpus, then convert compatible additions |
| [0045](0045-stage3-morgan-direction-baseline.md) | Morgan baseline and shared shallow-MLP comparison with future transcriptomic Δ phenotype predictors |
| [0046](0046-perturbdiff-standalone-training.md) | One PerturbDiff hyperparameter step with separate Tahoe and tenv studies |
| [0047](0047-perturbdiff-tahoe-transfer-to-tenv.md) | Tahoe transfer comparison across three tenv split regimes |
| [0048](0048-perturbdiff-transfer-diagnosis.md) | Paired controls for the Tahoe transfer gap and Step 06 output placement |
| [0049](0049-perturbdiff-transfer-optimization-comparisons.md) | Five Step 06 Transfer comparisons on cold-cm-drug with one dose-stratified figure |
| [0051](0051-latent-xpert.md) | Latent XPert autoencoder and latent perturbation predictor |
| [0052](0052-stage2-autoencoder-reconstruction-boundary.md) | Step 07 owns transcriptome reconstruction, and Step 08 owns perturbation prediction |
| [0053](0053-step07-compact-matrix-autoencoder.md) | Compare a 64 × 8 matrix latent representation with the Step 07 reconstruction baselines |
| [0054](0054-step07-separate-reconstruction-and-nuisance-tokens.md) | Separate expression reconstruction and nuisance tokens in a Step 07 attention autoencoder |
| [0055](0055-step07-matrix-mlp-reconstruction.md) | Compare the vector MLP and an AAE-shaped matrix MLP in Step 07 |
| [0056](0056-step07-shared-expression-representation.md) | Compare five matrix reconstruction models with a shared frozen expression representation interface |
| [0057](0057-step09-xpert-landmark-reference.md) | Use an internal landmark reference and a shared scoring panel for the Step 09 XPert ablation |
| [0058](0058-step09-direction-matched-pr.md) | Direction-matched significance retrieval and fixed-group pooled PR ending at attainable recall |
| [0059](0059-stage2-significance-retrieval.md) | Three full-CM, significance-and-direction matched retrieval metrics and one shared Stage 2 evaluator |

| [0060](0060-step09-round1-rerun.md) | Rerun Step 09 Round 1 with two models, four runs per model and validation early stopping |
| [0061](0061-stage2-shared-bases-and-indexed-corpora.md) | Shared Stage 2 bases, indexed gene identity and separate Tahoe/tenv training corpora |

| [0062](0062-stage2-accepted-products-and-reproducible-figures.md) | Stage 2 stable accepted products, CSV-reproducible figures and separate audit evidence |
| [0063](0063-tenv-comparison-direct-and-indirect-metrics.md) | Rebuilt tenv comparison: direct and indirect metrics, fold/seed index tables, stage palette, PRnet — *indirect metrics and figures superseded by 0064* |
| [0064](0064-tenv-reliability-threshold-scan.md) | Score the tenv comparison by a predicted-effect threshold scan and a reliability score |
| [0065](0065-stage2-gene-axis-decision-preparation.md) | Prepare the Stage 2 gene-axis decision: public disease sources, full-transcriptome DEGs and an unweighted rank-sum fusion, without deciding the axis — *superseded by 0066* |
| [0066](0066-stage2-corpus-view-and-gene-axis.md) | Stage 2 01_bases: raw counts, a raw corpus, gene-axis views and the released 2,000/5,000-gene axes |
| [0067](0067-errata-to-0066-gse270788-arms-and-view-shares.md) | Errata to ADR 0066: GSE270788 arm sizes and released-axis CP10K shares |
| [0068](0068-fusion-amendment-degenerate-sources-and-rank-norm.md) | Fusion amendment: non-informative sources excluded from fusion/`pN`, unmeasured `rank_norm`, and two removals |
| [0069](0069-remove-af-category-from-01-bases.md) | Remove the `AF` category from Stage 2 `01_bases`: drop the non-informative disease-DEG product and its gene-axis column/row surface |
| [0070](0070-af-stays-out-of-the-gene-axis.md) | AF stays out of the gene axis: the clinical CM-specific audit; informativeness is necessary but not sufficient (*amends 0068, confirms 0069*) |
| [0071](0071-gene-axis-v2-deg-tiered-sources-and-f2000.md) | Gene-axis v2: tiered source collection, DEG-only rank fusion, `DEG2000`/`f2000` and the tenv-HVG fallback (*supersedes 0066 eligibility/fusion, 0069/0070 AF removal*) |
| [0072](0072-disease-eligibility-and-unified-tox.md) | Disease eligibility at 100 axis-resolved DEGs, retained AF evidence and dataset-level unified TOX (*measured successor accepted 2026-11-16*) |
| [0073](0073-full-source-tier-first-deg1000-and-four-axes.md) | Full-source collection, highest DEG source tier, tier-first DEG1000/f1000, four nested axes/views and portrait category-split rule DOCX (*coherently measured release, #220–#225*) |
| [0074](0074-single-dose-tenv-5k-native-comparison.md) | Single-dose G_5000 native model comparison, shared splits and parent-level Δ scoring (*supersedes 0063/0064 current evaluation duties; old release retained until acceptance*) |
| [0075](0075-transigen-all-dose-input-selection.md) | TranSiGen all-dose inputs: seven fixed candidates, nested cold-drug selection and auxiliary held-dose evaluation |

| [0076](0076-dose-input-round1-exploration-boundary.md) | Dose-input round1 exploration, run-level report spread and retirement of final weight publication (*amends 0075*) |

| [0077](0077-five-k-only-axis-local-normalization.md) | 5k-only expression views with axis-local normalization; retain four named axes and old downstream input generations (*implemented paired expression successor, #244–#245*) |

| [0078](0078-dual-source-axis-local-dose-round1.md) | Rebuild dose-input round1 on the paired axis-local 5k views with sample-equal evaluation and cold-drug-dose |

| [0079](0079-stage2-drug-representation-step.md) | Isolate the drug representation comparison and migrate dose and transfer to Steps 04 and 05 |

| [0080](0080-tenv-kpgt-dose-round1.md) | Compare dose inputs on tenv with frozen Step 03 KPGT representations |

| [0081](0081-transigen-width-and-dose-comparison.md) | Compare TranSiGen width and dose encoding on the frozen tenv cold-drug-dose folds |

| [0082](0082-view-variants-tenv-tahoe-replicate-structure.md) | Test whether within-condition replicate structure explains the tenv–Tahoe gap with eight expression view variants (*study step `stage2/tsg_view`*) |
| [0083](0083-stage2-corpus-axis-view-boundaries.md) | Stage 2 corpus, prediction axes and expression view publication |
| [0084](0084-retire-stage2-bases-and-reserve-loss-step.md) | Retire Stage 2 Bases and reserve the loss comparison step |
| [0085](0085-stage2-raw-and-training-data.md) | Publish source-native raw data and plate-matched training data |
| [0086](0086-tenv-axis-matched-plate-control.md) | Condition the tenv axis comparison on matched plate controls in every partition |
| [0087](0087-stage2-active-steps-and-history.md) | Keep Steps 01–04 active and retain Stage 2 history in one local archive |
| [0088](0088-transigen-delta-loss-comparison.md) | Reopen Step 05 to compare TranSiGen Delta losses and weights on tenv |
| [0089](0089-transigen-delta-metric-weight.md) | Vary the TranSiGen Delta loss coefficient independently in Step 05.1 |
| [0090](0090-transigen-mse-l1-sweep.md) | Add a Delta L1 penalty to the native TranSiGen objective in Step 05.2 |
| [0091](0091-transigen-reconstruction-mse-sweep.md) | Vary the shared reconstruction MSE weight in Step 05.3 |
| [0092](0092-transigen-bounded-delta-similarity-sweeps.md) | Sweep bounded Delta similarity loss weights in Steps 05.4 and 05.5 |
| [0093](0093-stage2-shared-plotting-and-in-place-figures.md) | Share Stage 2 plot primitives and replace presentation-only figures at their accepted paths |
| [0094](0094-transigen-shape-round1.md) | Compare TranSiGen expression widths and depths in Stage 2 Step 06.1 |
