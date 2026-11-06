# Morgan integration checks — 2026-11-06

Feature baseline: `c17686dcf1409d8b069c8540685ddf28a3a9a1e3`.
Implementation is on `dev`; unrelated Stage 2 commits are outside this change.

The complete step suite passed 28 tests in 26.000 seconds (27.28 seconds wall),
with maximum RSS 853,960 KiB and zero swaps. The command used the `arc` Python,
two library threads, hidden CUDA devices, `/usr/bin/time -v`, and a systemd scope
with `MemoryMax=10G` and `MemorySwapMax=8G`:

```bash
python -m unittest discover -s stage3/02_base_morgan/execution/tests -v
```

After strengthening detection of external training entry points, all five
resource tests and all three benchmark tests passed again. The benchmark tests
took 6.941 seconds (7.39 seconds wall), maximum RSS 841,784 KiB, with the portable
10 GiB address-space cap and two threads. These include all 60 actual tiny-model
smoke fits with explicitly simulated host telemetry, plus a genuinely capped
lightweight child. Production commands expose no telemetry bypass.

Real frozen-input preparation for all three settings passed in 5.33 seconds,
maximum RSS 242,092 KiB. It preserved 5,488 admitted measurements, 663 parents,
ten splits and 60 expected identities. This preparation check performed no fits.

Both repository routing/tracking gates passed:

```bash
bash script/check_paths.sh
bash script/check_tracked.sh
```

Formal execution, resource evidence and prediction-based replay are recorded
separately in `formal-execution.md` and `../release_record.json`. The final
independent review and broad offline verification are recorded at delivery.
