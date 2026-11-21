# tenv KPGT formal acceptance

On 2026-11-21 (Asia/Shanghai), [#275](https://github.com/Moke110/arc/issues/275)
completed the ADR 0080 comparison with 42 TranSiGen fits, six KPGT_MLP fits
and six Mean_cell evaluations. The [acceptance receipt](acceptance.json)
binds the frozen inputs, fit manifests and replay receipts, validation locks,
scope and timing evidence, CSV tables and all three PNGs.

| Driver | Successful phases | Wall time | Peak RSS | GPU/WSL samples |
|---|---:|---:|---:|---:|
| TranSiGen | 168 | 42:40.35 | 1,929,988 KiB | 241 |
| KPGT_MLP and Mean_cell | 24 | 2:03.83 | 1,608,452 KiB | 13 |

The drivers ran singly in one persistent 10 GiB RAM and 8 GiB swap capped
service. The lowest sampled WSL RAM availability was 40,087,676 KiB and
the highest sampled GPU usage was 3,656 MiB. No fit phase failed in
the accepted retry. Exact fresh-process replay covered both partitions
for each fit and the test partition for each Mean_cell evaluation. The
publisher also reproduced all Mean_cell predictions with the independent
Step 02 reference.

The initial scope stopped after its telemetry sampler failed. Its partial
artifacts remain under `execution/checkpoints/tenv-kpgt-v1/`. The accepted
retry started afresh at 14:47 local time and retains all fit, selection, scope,
timing and monitor bytes in the `retry/` subdirectory. The input freeze is in
the sibling `inputs/` directory. The accepted training plan pins the original
fit-time driver at commit `c0c1af6`. Its verification wrapper now retains the
nested origin identity in results. The wrapper accepts that exact former
driver digest only during read verification and compares every other frozen
field. The receipt separately binds the corrected verifier.

The corrected formal verifier checked all 42 origins in 4.33 seconds.
Scientific preparation took 12.76 seconds. After visual inspection of all
three plots, acceptance rechecked the source fits, controls, resource window
and all CSV digests and reproduced byte-identical PNGs in 9.30 seconds.
The installed successor and both sealed predecessor inventories matched
their receipts.
