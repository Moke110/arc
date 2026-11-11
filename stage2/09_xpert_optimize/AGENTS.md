# AGENTS.md — stage2/09_xpert_optimize

- Study contract and reproduce commands live in [README.md](README.md). Read
  [ADR 0057](../../docs/adr/0057-step09-xpert-landmark-reference.md) before changing
  gene axes, XPert arms, or the validation and test boundary.
- Prepare landmark inputs with `execution/01_corpus/build.py`. Raw test row joins
  and test expression stay under the evaluator release.
- Prepare the three fold-local HVG panels and the Round 1 G-axis views with
  `execution/01_corpus/gene_axes.py`. Freeze the shared panel before fitting.
- Fit and replay with `execution/script/xpert.py`. Give workers the predictor
  directory and selected checkpoint.
- Original Round 1 test dispatch and label rebuilding validate the committed
  score-free `execution/round1_record/promotion-provenance.json` pointer to the
  original Git seal. Rescoring uses the shared evaluator and creates no new
  promotion or checkpoint selection.
- Before real-data jobs, read the [WSL resource rules](../../docs/agents/execution-environment.md)
  and use `execution/script/run_capped.sh`. Run long fits serially.

- Active retrieval scoring and figures use [`../evaluation/README.md`](../evaluation/README.md), with `AUPRC`, `AUPRC_positive`, and `AUPRC_negative` on the full CM significance axis.
