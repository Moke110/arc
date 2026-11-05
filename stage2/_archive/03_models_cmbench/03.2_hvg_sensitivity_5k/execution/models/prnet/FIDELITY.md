# PRnet fidelity ledger (ticket 13)

> Derived from `stage2/03_models_cmbench/03.1_cold_drug_18k/models/prnet/FIDELITY.md`; the gene-space rows are re-declared for the 5,000-HVG sensitivity axis (ADR 0028).

Classes follow ticket 09: `config_only`, `adapter_extension`,
`recipe_omission`, `recipe_substitution`. Every row states an expected effect;
the target is zero *undeclared* deviations.

| component | author setting | run setting | class | reason | expected effect | ticket |
|---|---|---|---|---|---|---|
| gene space | 978 LINCS / 5,000 Sci-Plex | 5,000 rank-sum HVG (corpus subset) | recipe_substitution | fixed benchmark gene universe (ADR 0028); width matches Sci-Plex | small: the output width matches the widest author run; gene identity differs | 09 |
| corpus scale | L1000 (1.3 M bulk) / Sci-Plex | 3,600 wells/fold, 276 train drugs | recipe_substitution | the curated CM view is the benchmark | large negative: fewer perturbagens and cell contexts | 09 |
| initialisation | released LINCS/Sci-Plex checkpoints for the demo | random init, no released weights | recipe_omission | ticket 09 forbids starting from public weights | large negative at fixed epochs | 09 |
| expression space | `normalize_total` (median library) + `log1p` | frozen `log1p_cp10k` | recipe_substitution | test wells expose no library size, so an exact inverse transform is impossible; Δ is defined in cp10k log1p | small, monotone per-cell rescaling; no effect on direction ranking | 09 |
| control input | all dose-0 cells appended to every split, paired per well | matched plate DMSO mean only from the fold's training role | adapter_extension | the author pool leaks valid/test controls across folds (ticket 13 item 6) | correctness fix; no expected score effect | 13 |
| training rows | treated + pooled dose-0 cells | training-role treated wells; controls enter as inputs | adapter_extension | fold isolation; controls are never held out | none expected | 04 |
| prediction | sample the Gaussian (NB for LINCS) | Gaussian mean head, zero noise | config_only | deterministic point prediction for one-shot test scoring | slightly lower variance than sampling; no expected bias | 09 |
| scheduler patience | 5 in `train_sciplex.py` | 10 (ticket-09 frozen table) | config_only | ticket 09 froze patience 10 | negligible | 09 |
| drug encoding | rFCFP4 1024 x log10(dose+1) | identical | config_only | author function reused verbatim | none | 09 |
| dose | raw µM | raw µM | config_only | benchmark doses are 0.2/1/10 µM | none | 09 |
| data preprocessing fit | per-cell, no fit | none fitted | config_only | PRnet has no dataset-level preprocessing | none | 04 |
