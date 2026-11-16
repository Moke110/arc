# Single-dose tenv 5k native comparison

Accepted 2026-11-16 under [#226](https://github.com/Moke110/arc/issues/226).
The successor `stage2/02_models_tenv-5k` compares randomly initialized native
PRnet and TranSiGen with Mean, Mean_cell, Mean_drug and the project MLP_Morgan
on the released G_5000, 1 µM tenv population (1,402 wells, 459 parents).
This isolates architecture from pretraining and dose adaptation. All learners
share Morgan2048, cold-drug/warm-start five-fold identities and seeds 17/42/2026;
complete fixed-randomness validation treated-expression MSE selects the strictly
lowest checkpoint with patience 50 and at most 1,000 epochs.

This supersedes ADR 0063/0064's current model-comparison and evaluation duties:
primary scoring now aggregates repeated-well Δ within parent before computing
PCC, with parent-paired bootstrap and explicit constant/failure accounting.
Historical evidence retains its original scientific context. The previous step
remains accepted until the successor's complete 90-run matrix, reload, metric,
CSV-only figure and downstream acceptance passes; only then may #231 retire it.
ADR 0062's candidate publication rules still apply. Stage 3 remains frozen.
