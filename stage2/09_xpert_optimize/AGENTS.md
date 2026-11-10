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
- Seal validation promotions with `execution/script/evaluate.py seal`. Commit the
  seal on `dev` before building test labels or scoring the test split.
- Before real-data jobs, read the [WSL resource rules](../../docs/agents/execution-environment.md)
  and use `execution/script/run_capped.sh`. Run long fits serially.
