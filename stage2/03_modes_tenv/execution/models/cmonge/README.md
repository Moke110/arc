# CMonge runtime admission

Owned by [ticket #50](https://github.com/Moke110/arc/issues/50). The
[FIDELITY.md](FIDELITY.md) ledger explains the native direct-space transport,
shared Morgan adapter, objective and recipe deviations. Source compatibility
and small synthetic tests do not establish formal admission.

The isolated interpreter is `stage2/_envs/cmonge/bin/python` (existing CPU JAX).
The adapter command is that interpreter followed by
`stage2/03_modes_tenv/execution/script/cmonge_adapter.py --config <config>`;
`benchmark.run_one` appends the safe `--request` path. The scheduler owns the
full-space qualification and formal dispatch; direct uncapped jobs are prohibited.
`qualification.json` is a bounded train/save/reload/predict probe, while
`formal.json` retains the 10000-iteration, 512-sample author recipe.

The synthetic artifact test checks joint adapter weight updates, finite native
training/validation objective, baseline and dose participation, exact checkpoint
reload and binary checkpoint checksum validation. Run from repository root:

```bash
systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=8G -- /usr/bin/time -v \
  python -m unittest discover -s stage2/03_modes_tenv/execution/tests -p test_cmonge.py -v
```

Measured synthetic test: PASS, 3.812 seconds unittest time, 642540 KiB peak RSS
(4.26 seconds wall including test runner). An earlier fixture with entirely zero
chemical fingerprints correctly failed adapter-kernel participation; the fixture
now supplies nonzero binary features. No result from this probe is a benchmark
score or evidence of full-axis admission. Real-data admission evidence is recorded
separately after the scheduler-owned full-axis probe.
