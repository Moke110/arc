# TranSiGen training and replay evidence

[Step 04 axis study](../../04_tsg_axis/README.md) ·
[Step 05 module study](../../05_tsg_module/README.md) ·
[CarPS-Transcriptome model verification](cpst-verification.md)

Use the study's recorded data digests, source revision and replay commands.
The current Step 04 formal results and Step 05 bounded results use the
frozen five-fold tenv cold-drug protocol. Resolve older registrations
against the Git revision they record. Copy checkpoint folders before
running a replay that writes receipts or temporary products.

## Historical TranSiGen loss and module studies

The [pre-cutover CPST revision](https://github.com/Moke110/carps/tree/325be5dd45e2b670f1e01882e225aa5debbeaab7/cpst)
contains the former Step 05 loss comparison and all four Step 06 module
rounds, including the Round 1 width extension. Each study's report and
guide document its training and replay. The
[historical verification inventory](https://github.com/Moke110/carps/blob/325be5dd45e2b670f1e01882e225aa5debbeaab7/cpst/util/training/verification.md)
lists the registration revisions and formal payload layout. The
[historical figure gate](https://github.com/Moke110/carps/blob/325be5dd45e2b670f1e01882e225aa5debbeaab7/cpst/util/figure/script/verify_report_figures.py)
recreates the 45 PNGs from that revision's sealed plot tables.

The materialized workspace holds the preserved loss products at
`archived/spec27-transigen/loss/<study>/output/` and the module products at
`archived/spec27-transigen/round<n>/output/`. These are local archives.
The loss products and module fits were linked into this archive from the
prior installed outputs on the same ext4 filesystem. Their 74 symlinks
still resolve to their original targets outside the retired active roots.
Treat the preserved products as read-only. The Step 05 bounded report and
receipts are also preserved under `archived/spec27-transigen/spec27-step05-bounded/`.

The archived `round4/` additionally includes the accepted publication
receipt, fit ledger, bandwidths and fold scores. Pass that directory to
the existing `--round-four` option in alpha and chemical-reference replay
when running from a worktree without a local copy of the archive. The
[original fit revision](https://github.com/Moke110/carps/blob/f387f2c370eb87cd27c02502aa98f165851491a7/stage2/06_tsg_shape/06.4_round4/README.md)
provides the original Step 06.4 fit entry. The
[alpha release manifest](https://github.com/Moke110/carps/blob/a7613b74fc23f9f48d8ea0d4479542c008d701db/cpst/07_cpst_alpha/output/release.json)
pins the source and receipt hashes for the historical reference.
