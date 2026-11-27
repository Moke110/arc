# Spec 24 documentation authors

## Leaf documents

The 58 document assignments for the CPST and CPSP leaf rewrite are recorded below.
Each author reviewed the previous document at `9f5511b9` and the current owning
sources and evidence. Document links use the canonical paths. The final column
identifies the accepted material and interfaces checked for each assignment.

CPST evidence can be traced through the linked studies' `log/` receipts, the
[shared verification records](../../cpst/util/README.md) and the
[release migration inventory](../../cpst/util/execution/385-release-migration.json).
CPSP evidence is recorded in the
[corpus manifest](../../cpsp/01_corpus/output/manifest.json) and the
[Morgan release record](../../cpsp/02_base_morgan/execution/release_record.json).
The original authors completed their own follow-up document edits. For Round 3,
the coordinator independently confirmed that the README retained its scientific
text and headings, the guide's four Bash blocks passed syntax checks, and the
original verifier resolves at `28f5eea4`. Both Round 3 authors are credited
below for their document and scientific evidence work.

| Document | Author | Evidence checked |
| --- | --- | --- |
| [`cpst/02_gene_axis/02.1_hvg-deg/AGENTS.md`](../../cpst/02_gene_axis/02.1_hvg-deg/AGENTS.md) | `/root/coordinator/doc388_01` | Axis builder, projection API and axis manifest hashes, native corpus coverage, #342 figure receipt. |
| [`cpst/02_gene_axis/02.1_hvg-deg/README.md`](../../cpst/02_gene_axis/02.1_hvg-deg/README.md) | `/root/coordinator/doc388_02` | Four released HVG-DEG axes, nesting and coverage, #342 figures, downstream Step 03 views and Steps 04 and 09 comparisons. |
| [`cpst/02_gene_axis/02.2_denoise-deg/AGENTS.md`](../../cpst/02_gene_axis/02.2_denoise-deg/AGENTS.md) | `/root/coordinator/doc388_03` | Build commands and five public axis reads, disease evidence tests, axis and native corpus manifests, accepted ranking rebuild and four figure PNG replays. |
| [`cpst/02_gene_axis/02.2_denoise-deg/README.md`](../../cpst/02_gene_axis/02.2_denoise-deg/README.md) | `/root/coordinator/doc388_04` | Published ranking and ordered five-axis release, #342 figure hashes, matched native corpus manifest, source for diagnostic and axis criteria. |
| [`cpst/03_data_train/03.1_tenv_hvg-5k/AGENTS.md`](../../cpst/03_data_train/03.1_tenv_hvg-5k/AGENTS.md) | `/root/coordinator/doc388_05` | tenv builder and readers, ordered axis and fold ledgers, #312 manifest and replay, bounded control and Delta reads. |
| [`cpst/03_data_train/03.1_tenv_hvg-5k/README.md`](../../cpst/03_data_train/03.1_tenv_hvg-5k/README.md) | `/root/coordinator/doc388_06` | #312 exact replay and promotion, eight installed payload entries, reader and builder interfaces, Steps 04 and 09 consumers. |
| [`cpst/03_data_train/03.2_tahoe_hvg-5k/AGENTS.md`](../../cpst/03_data_train/03.2_tahoe_hvg-5k/AGENTS.md) | `/root/coordinator/doc388_07` | Tahoe publisher and reader, #312 release and combined-view pins, all three installed H5AD headers. |
| [`cpst/03_data_train/03.2_tahoe_hvg-5k/README.md`](../../cpst/03_data_train/03.2_tahoe_hvg-5k/README.md) | `/root/coordinator/doc388_08` | #312 and #323 release and smoke pins, Tahoe H5AD and CSV headers, publisher and owner verification routes. |
| [`cpst/03_data_train/03.3_tenv_hvg-3k/AGENTS.md`](../../cpst/03_data_train/03.3_tenv_hvg-3k/AGENTS.md) | `/root/coordinator/doc388_09` | #313 axis, fold, native corpus and reference pins, eight installed files, float32 array headers and named readers. |
| [`cpst/03_data_train/03.3_tenv_hvg-3k/README.md`](../../cpst/03_data_train/03.3_tenv_hvg-3k/README.md) | `/root/coordinator/doc388_10` | #313 exact replay and fold parity, all eight installed files, named reader and builder routes, Steps 04 and 09. |
| [`cpst/03_data_train/03.4_tenv_hvg-4k/AGENTS.md`](../../cpst/03_data_train/03.4_tenv_hvg-4k/AGENTS.md) | `/root/coordinator/doc388_11` | #313 exact replay, installed release, frozen 5k reference and fold manifests, reader and replay options. |
| [`cpst/03_data_train/03.4_tenv_hvg-4k/README.md`](../../cpst/03_data_train/03.4_tenv_hvg-4k/README.md) | `/root/coordinator/doc388_12` | #313 population and replay receipt, installed manifest and payload paths, axis and fold provenance, downstream reports. |
| [`cpst/03_data_train/03.5_tenv_denoise-3k/AGENTS.md`](../../cpst/03_data_train/03.5_tenv_denoise-3k/AGENTS.md) | `/root/coordinator/doc388_13` | #313 exact replay, installed manifest and payload paths, Denoise-DEG axis and fold pins, reader and builder code. |
| [`cpst/03_data_train/03.5_tenv_denoise-3k/README.md`](../../cpst/03_data_train/03.5_tenv_denoise-3k/README.md) | `/root/coordinator/doc388_14` | #313 replay and population receipt, Denoise-DEG axis and installed release manifests, reader and downstream reports. |
| [`cpst/03_data_train/03.6_tenv_denoise-4k/AGENTS.md`](../../cpst/03_data_train/03.6_tenv_denoise-4k/AGENTS.md) | `/root/coordinator/doc388_15` | #313 exact replay, installed Denoise-DEG 4k release and input pins, publication and 256-row reader routes. |
| [`cpst/03_data_train/03.6_tenv_denoise-4k/README.md`](../../cpst/03_data_train/03.6_tenv_denoise-4k/README.md) | `/root/coordinator/doc388_16` | #313 replay and installed manifest, eight payload paths, axis and frozen fold records, Steps 04 and 09 use. |
| [`cpst/03_data_train/03.7_tenv_denoise-5k/AGENTS.md`](../../cpst/03_data_train/03.7_tenv_denoise-5k/AGENTS.md) | `/root/coordinator/doc388_17` | #313 formal replay, installed Denoise-DEG 5k manifest, frozen source and reference pins, publisher and reader. |
| [`cpst/03_data_train/03.7_tenv_denoise-5k/README.md`](../../cpst/03_data_train/03.7_tenv_denoise-5k/README.md) | `/root/coordinator/doc388_18` | #313 and #371 manifest pins, axis and split, installed payload paths, downstream Steps 04–10. |
| [`cpst/03_data_train/03.8_tenv_hvg-2k/AGENTS.md`](../../cpst/03_data_train/03.8_tenv_hvg-2k/AGENTS.md) | `/root/coordinator/doc388_19` | #366 replay and independent counts, installed HVG-DEG 2k release, builder and reader routes. |
| [`cpst/03_data_train/03.8_tenv_hvg-2k/README.md`](../../cpst/03_data_train/03.8_tenv_hvg-2k/README.md) | `/root/coordinator/doc388_20` | #366 and #371 matched manifest and counts, eight installed files, separate replay and count audit. |
| [`cpst/03_data_train/03.9_tenv_denoise-2k/AGENTS.md`](../../cpst/03_data_train/03.9_tenv_denoise-2k/AGENTS.md) | `/root/coordinator/doc388_21` | #366 replay and independent counts, installed Denoise-DEG 2k payload sizes and frozen reference pins. |
| [`cpst/03_data_train/03.9_tenv_denoise-2k/README.md`](../../cpst/03_data_train/03.9_tenv_denoise-2k/README.md) | `/root/coordinator/doc388_22` | #366 and #371 matching release digest and counts, eight installed paths, named reader and audit. |
| [`cpst/05_tsg_loss/05.1_loss_metrics/AGENTS.md`](../../cpst/05_tsg_loss/05.1_loss_metrics/AGENTS.md) | `/root/coordinator/doc388_23` | #329 and #330 receipts, native anchors and frozen inputs, current figure hashes, sealed runner pins at `997ca0f1`. |
| [`cpst/05_tsg_loss/05.1_loss_metrics/README.md`](../../cpst/05_tsg_loss/05.1_loss_metrics/README.md) | `/root/coordinator/doc388_24` | #330 formal and curve receipts, 36-fit ledger, 324 fold scores, metric optima and 4,920 epoch rows, report and figure hashes. |
| [`cpst/05_tsg_loss/05.2_mse_l1/AGENTS.md`](../../cpst/05_tsg_loss/05.2_mse_l1/AGENTS.md) | `/root/coordinator/doc388_25` | #332 and #333 receipts, 27 registrations and three anchors, frozen inputs, pinned replay at `4d85563e`. |
| [`cpst/05_tsg_loss/05.2_mse_l1/README.md`](../../cpst/05_tsg_loss/05.2_mse_l1/README.md) | `/root/coordinator/doc388_26` | #333 reports, 243 fold scores and 81 aggregates, metric optima, Step 05.3 reference and figure hashes. |
| [`cpst/05_tsg_loss/05.3_mse-l1_mse-x/AGENTS.md`](../../cpst/05_tsg_loss/05.3_mse-l1_mse-x/AGENTS.md) | `/root/coordinator/doc388_27` | #335 and #336 receipts, 27 registrations and three anchors, fold scores and figures, sealed source at `fa7e9a0e`. |
| [`cpst/05_tsg_loss/05.3_mse-l1_mse-x/README.md`](../../cpst/05_tsg_loss/05.3_mse-l1_mse-x/README.md) | `/root/coordinator/doc388_28` | #336 report, contrasts and optima, 27 registration/manifest pairs, Step 05.2 anchor and figure hashes. |
| [`cpst/05_tsg_loss/05.4_mse-l1_pcc/AGENTS.md`](../../cpst/05_tsg_loss/05.4_mse-l1_pcc/AGENTS.md) | `/root/coordinator/doc388_29` | #338 and #339 receipts, 18 PCC registrations and three Step 05.3 anchors, source and figure hashes at `6e903171`. |
| [`cpst/05_tsg_loss/05.4_mse-l1_pcc/README.md`](../../cpst/05_tsg_loss/05.4_mse-l1_pcc/README.md) | `/root/coordinator/doc388_30` | #339 PCC report, 18 replay receipts and metric optima, figure hashes, all three smoke renderer inputs. |
| [`cpst/05_tsg_loss/05.5_mse-l1_cos/AGENTS.md`](../../cpst/05_tsg_loss/05.5_mse-l1_cos/AGENTS.md) | `/root/coordinator/doc388_31` | #338 and #339 receipts, 18 cosine registrations and three Step 05.3 anchors, CLI and figure hashes, pinned replay. |
| [`cpst/05_tsg_loss/05.5_mse-l1_cos/README.md`](../../cpst/05_tsg_loss/05.5_mse-l1_cos/README.md) | `/root/coordinator/doc388_32` | #339 cosine report and 18 replays, shared bandwidth and anchor provenance, figures and three smoke CSVs. |
| [`cpst/06_tsg_shape/06.1_round1/AGENTS.md`](../../cpst/06_tsg_shape/06.1_round1/AGENTS.md) | `/root/coordinator/doc388_33` | #347 and #350 installed receipts, 30 registrations and replay records, Step 05.2 parity, verification at `fc680914`. |
| [`cpst/06_tsg_shape/06.1_round1/README.md`](../../cpst/06_tsg_shape/06.1_round1/README.md) | `/root/coordinator/doc388_34` | #347 and #350 ledgers, 270 scores and 90 aggregates, parity and figure hashes, owner reproduction links. |
| [`cpst/06_tsg_shape/06.2_round2/AGENTS.md`](../../cpst/06_tsg_shape/06.2_round2/AGENTS.md) | `/root/coordinator/doc388_35` | #352 and #353 receipts, 21 replay records and three reference anchors, figure hashes, formal verification at `87564cda`. |
| [`cpst/06_tsg_shape/06.2_round2/README.md`](../../cpst/06_tsg_shape/06.2_round2/README.md) | `/root/coordinator/doc388_36` | #353 paired means and SDs, 216 fold scores, 72 aggregates, all fit source pins at `f187a4dc`, publication at `87564cda`. |
| [`cpst/06_tsg_shape/06.3_round3/AGENTS.md`](../../cpst/06_tsg_shape/06.3_round3/AGENTS.md) | `/root/coordinator/doc388_37` | #356 matrix and report manifests, all 72 registration/seal pairs and input pins, source and current figure hashes. |
| [`cpst/06_tsg_shape/06.3_round3/README.md`](../../cpst/06_tsg_shape/06.3_round3/README.md) | `/root/coordinator/doc388_38` | #356 failed and passing attempts, 72 ledger rows and 12 recomputed contrasts, pinned source and figure hashes. |
| [`cpst/06_tsg_shape/06.4_round4/AGENTS.md`](../../cpst/06_tsg_shape/06.4_round4/AGENTS.md) | `/root/coordinator/doc388_39` | #358 and #359 receipts, 27 exact replay seals, frozen inputs, 24 pinned code digests at `f387f2c3`, publication at `09da1ead`. |
| [`cpst/06_tsg_shape/06.4_round4/README.md`](../../cpst/06_tsg_shape/06.4_round4/README.md) | `/root/coordinator/doc388_40` | #359 primary means and SDs, eight contrasts, 27 seals and publication artifacts, original source and figure pins. |
| [`cpsp/01_corpus/AGENTS.md`](../../cpsp/01_corpus/AGENTS.md) | `/root/coordinator2/doc388_41` | #386 corpus acceptance, three code digests, five installed source files, eight baseline/H-R payload pins at `2fd12768`. |
| [`cpsp/01_corpus/README.md`](../../cpsp/01_corpus/README.md) | `/root/coordinator2/doc388_42` | #386 rebuilt and installed corpus manifests, both arms' coverage and payload hashes, desalt review and Morgan corpus verification. |
| [`cpsp/02_base_morgan/AGENTS.md`](../../cpsp/02_base_morgan/AGENTS.md) | `/root/coordinator2/doc388_43` | #386 Morgan installation verify, 60 run records and checkpoint symlinks, three aggregate hashes and workflow CLI. |
| [`cpsp/02_base_morgan/README.md`](../../cpsp/02_base_morgan/README.md) | `/root/coordinator2/doc388_44` | Morgan split and release records, original formal replay, #386 verification of 60 fits, ADR 0045 and the CPST delta comparison scope. |
| [`cpst/01_data_raw/AGENTS.md`](../../cpst/01_data_raw/AGENTS.md) | `/root/coordinator2/doc388_45` | Native v2 release and migration receipts #303/#306, 17 installed paths, ten compact output hashes and current reader CLI. |
| [`cpst/01_data_raw/README.md`](../../cpst/01_data_raw/README.md) | `/root/coordinator2/doc388_46` | Native v2 replay and relocation, seven dataset manifests, membership ledger and eight downstream tenv manifests. |
| [`cpst/04_tsg_axis/AGENTS.md`](../../cpst/04_tsg_axis/AGENTS.md) | `/root/coordinator2/doc388_47` | #318, #319, #322 and #323 receipts, 18 registrations, inputs, report and figure digests, seven CLI parsers. |
| [`cpst/04_tsg_axis/README.md`](../../cpst/04_tsg_axis/README.md) | `/root/coordinator2/doc388_48` | #319 and #322 metrics and reports, 78 formal result cells, 18 seal pairs, #323 bounded smoke and figure hashes. |
| [`cpst/07_cpst_alpha/AGENTS.md`](../../cpst/07_cpst_alpha/AGENTS.md) | `/root/coordinator2/doc388_49` | #362 alpha receipt, nine registrations, 30 alpha and 24 Round 4 source hashes, public adapter and bounded tests. |
| [`cpst/07_cpst_alpha/README.md`](../../cpst/07_cpst_alpha/README.md) | `/root/coordinator2/doc388_50` | Alpha release and #362 receipts, nine checkpoint/prediction pins, exact replay receipts, Steps 06.4 and 08–10. |
| [`cpst/08_cpst_chem-emb/AGENTS.md`](../../cpst/08_cpst_chem-emb/AGENTS.md) | `/root/coordinator2/doc388_51` | #363 and #364 receipts, nine historical reuse records and 54 fresh fits, 567 scores, paired threshold and CLI checks. |
| [`cpst/08_cpst_chem-emb/README.md`](../../cpst/08_cpst_chem-emb/README.md) | `/root/coordinator2/doc388_52` | #364 release, 63 coordinates and seven summaries, nine UniMol pins reused by Step 09, figure hashes. |
| [`cpst/09_cpst_axis/AGENTS.md`](../../cpst/09_cpst_axis/AGENTS.md) | `/root/coordinator2/doc388_53` | #367 gate, #368 reuse and #371 formal receipts, 63 new fits and nine reuses, 720 scores and 33,192 condition rows. |
| [`cpst/09_cpst_axis/README.md`](../../cpst/09_cpst_axis/README.md) | `/root/coordinator2/doc388_54` | #371 eight axis summaries and closest paired contrast, fresh/reused ledger, installed publication and figure hashes. |
| [`cpst/10_cpst_dose/AGENTS.md`](../../cpst/10_cpst_dose/AGENTS.md) | `/root/coordinator2/doc388_55` | #377 historical gate, #378 formal and figure receipts, 90 exact replays and 900 scores, 11 parsed CLI routes. |
| [`cpst/10_cpst_dose/README.md`](../../cpst/10_cpst_dose/README.md) | `/root/coordinator2/doc388_56` | #378 90-fit and 900-score ledgers, partition and summary counts, dose arm rankings, publication and figure hashes. |
| [`cpst/util/AGENTS.md`](../../cpst/util/AGENTS.md) | `/root/coordinator2/doc388_57` | Shared API symbols, gene identity and chemical representation pins, both runtime verification records and 36 stage receipt hashes. |
| [`cpst/util/README.md`](../../cpst/util/README.md) | `/root/coordinator2/doc388_58` | Shared API contracts and 28 utility tests, both representation manifests, 36 stage receipts and 14 figure receipts. |
