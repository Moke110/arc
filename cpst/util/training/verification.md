# Steps 04–06 consolidation verification

The [published receipts](formal-receipts.sha256) retain the exact experiment
inputs, commands, release manifests and result hashes. Verify this receipt
inventory from the repository root with
`sha256sum -c cpst/util/training/formal-receipts.sha256`.
The numbered study READMEs give the full report and figure verification commands.

| Study | Installed formal products under the study's `output/` | Registration code revision | Checked registrations |
| --- | --- | --- | ---: |
| 04 | `formal-matched-control/runs`, `report-matched-control` | `25000482dfa9c30a109318048711ea59f08d740b` | 18 |
| 05.1 | `formal/matrix`, `formal/report` | `997ca0f1ef5fd592b5faeb0979cb03cd5bf16f0c` | 36 |
| 05.2 | `formal/matrix`, `formal/report` | `4d85563ea6d973f1573cfebaac182cb301e857ae` | 27 |
| 05.3 | `formal/matrix`, `formal/report` | `fa7e9a0e580d06a63210e93bff82c655d37b5ae8` | 27 |
| 05.4 and 05.5 | Joint `formal/matrix` and `formal/reports` | `6e9031714000534714304b012959a20992feb214` | 36 shared |
| 06.1 | `formal/prior-347/matrix`, `formal/matrix`, `formal/report` | `1f166d2f0dee53444131eee1dd83481c8dedce4f`, `2a5ac0a7a6c7dba527820564455f638d7e834a81` | 24 prior, 6 extension |
| 06.2 | `formal/matrix`, `formal/report` | `f187a4dc0dfaac43eec0e37f409b3ffe4d2ee0b6` | 21 new |
| 06.3 | `formal/matrix`, `formal/report`, `formal/prior-unnormalized` | `d5c0554a2f9c0b2208c805c6a372bee1953dde58`, `29c7e1680e79b32d40de18950260ff57ad5615bd` | 36 normalized, 27 prior |
| 06.4 | `formal/matrix`, `formal/report` | `f387f2c370eb87cd27c02502aa98f165851491a7` | 27 new |

For all 285 installed registrations in the inventory, the manifest's
registration digest matched and all path-keyed code digests matched
`git cat-file blob <revision>:<registered-path>` from the listed revision.
Reused source coordinates retain their predecessor receipts. The alpha verifier
independently resolves its 30 runtime digests at
`825f48c96a0968c11b0bffa9a8b06686b0d16e03` and checks its Step 06.4 source
revision from the unchanged release.

## Registration parity

Opened the active ordered-axis and tenv training releases with the frozen
KPGT release and vendor source. Bound the smoke and formal rosters through
the original `f0d5bb673435b1a0341301c802ceac5094708a54` entry points and their
refactored equivalents. Removed the path-keyed `code_sha256` and checked
the new explicit `experiment` field against each bound model, ordered axis,
dose transform, objective, fold and seed before comparing the original
configuration. All 343 normalized records were byte-identical as sorted JSON.

| Study | Coordinates exercised |
| --- | ---: |
| 04 | 36 |
| 05.1 | 40 |
| 05.2 | 29 |
| 05.3 | 29 |
| 05.4 and 05.5 | 40 |
| 06.1 | 33 |
| 06.2 | 29 |
| 06.3 | 60 |
| 06.4 | 47 |

## CUDA and scored artifacts

Ran separate `register --smoke --fold 0 --device cuda`, `train`, `predict`
and `replay` processes on the Denoise-DEG 5k release. For each row the original
and current run produced equal 51 checkpoint tensors, bit-exact validation,
prediction and signed Delta arrays, identical learning traces, condition and
gene tables, and exact condition and fold scientific metrics. The signed
architecture also matched its KL diagnostic CSV.

| Numbered runner | Options in addition to the shared flag set | Fold metric CSV SHA-256 |
| --- | --- | --- |
| `04_tsg_axis/execution/run.py` | `--axis denoise-5k` | `ac42fe759beecd426ca0fed09066cded232c310cdbc41c8fec08ea0b4bb1a34a` |
| `05_tsg_loss/similarity/run_similarity.py` | `--metric pcc --alpha 0.5` | `c40d12f24adbcce7bc167ec1518a5985788f645146e10602da273557097fda4a` |
| `06_tsg_shape/06.4_round4/execution/run_measured.py` | `--shape XΔ-XΔ --seed 17` | `c87d38b55c6d1eb36acc6bb1ab76dacc5a546d2066fe8781bd30cc74b27098df` |

Also replayed a scratch copy of the sealed `XΔ-XΔ`, fold 0, seed 17
checkpoint under its original Step 06.4 Git revision using the
[documented command](../../06_tsg_shape/06.4_round4/README.md#replay-a-sealed-step-064-fit).
Validation objective, held-out predictions and signed Delta arrays all
replayed exactly. The current and original evaluators agreed exactly with each
other and with all nine published fold metrics and their eligible counts.
The scoring comparison recalculated the fold's published HVG-DEG 5k training
bandwidth, as required for the MMD calibration.

The environment was Python 3.13.13, PyTorch 2.12.0+cu130, CUDA 13.0 and NumPy
2.4.6 with two Torch threads on the RTX 4070 Ti SUPER. The serial six-fit
scope peaked at 1,986,040 KiB RSS and completed in 72.91 seconds. Its scope
used a 12 GiB memory cap and an 8 GiB swap cap.

All 89 regression cases passed across the relocated test groups, including
the two opt-in Step 06.1 predecessor cases with `TSG_PRIOR_MATRIX` set to
the installed `output/formal/prior-347/matrix`. The Stage 2 retired-path
audit and tracked-file check passed.
