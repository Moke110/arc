# MLP_Morgan bounded acceptance — 2026-11-16

Issue #230; project configuration of the XPert-defined baseline, not an
assertion of reproducing unavailable complete author baseline training code.

The real frozen `/tmp/arc-tenv5k-inputs` reader drove CPU cold-drug and warm-start
fold0 seed17 candidates, each with smoke_rows8 and max_epochs2. Both terminal
statuses are complete and explicitly smoke=true; weights updated, best epoch2
was restored, and a new Python process reproduced prediction.npy exactly.
Validation MSE was 0.2326799410146369 / 0.22070924950581766 respectively.
Predictions, identity, matched-control delta scores, learning curves, environment,
checkpoint and hashed manifest exist in `/tmp/arc-tenv5k-230-smoke-{cold,warm}`.
These candidates are not formal 1000-epoch-budget matrix entries.

Both ran serially under systemd MemoryMax10G/MemorySwapMax8G and `/usr/bin/time -v`,
threads2, with 41GiB available before launch. Peak RSS including fresh replay was
953528 / 953840 KiB, wall time4.90 /4.78s, no swaps, exit0. The coordinator held
the single heavy-task lock and monitored resources; lock released after replay.

Approved prediction-seam tests first failed for missing model_mlp then passed:
no treated-label/control input is accessible, repeated Morgan inputs yield the
same ordered prediction, inference seed does not alter deterministic output,
negative linear outputs remain, and training learns expression with a singleton
last batch. Shared strict stopping/best-restore tests remain owned by learner.
