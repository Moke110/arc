# Stage 2 raw data relocation

[Ticket #306](https://github.com/Moke110/arc/issues/306) relocates the accepted
source-native release. Readers and build commands use the new root. Relative
paths in the release inventories and the recorded source pins stay the same.

| Previous location | Current location |
|---|---|
| `stage2/01_corpus/output/` | `stage2/01_data_raw/output/` |
| `stage2/01_corpus/historical/native-v1/` | `stage2/01_data_raw/historical/native-v1/` |
| `stage2/01_corpus/log/` | `stage2/01_data_raw/log/` |

The per-dataset release `output/manifest.json` is sealed by SHA-256
`91a54af534bb34df97619814c38ea676a074e1ec826753095c6f55b068165a90`.
The archived v1 `historical/native-v1/manifest.json` remains the generation
referenced by the full 5k release. Its SHA-256 is
`a4cd6c42165569f77765a9c65de0d9d7e7292c192fdc7d25b455b78fdb50b4e4`.
The [original acceptance](acceptance.json) and [v2 conversion](migration-303.json)
retain their original path strings and evidence.

## Verification on 2026-11-24

- `build_corpus.py verify` at the renamed default output reproduced the v2
  manifest and every payload digest from the sealed v1 release. The bounded
  equivalence checks covered counts, observations, native genes and common-gene
  mappings. The primary-cardiomyocyte ledger retained 842 observations and all
  956 source memberships. The run used 858,768 KiB peak RSS and 4:05 wall time.
- `build_axes.py verify --corpus stage2/01_data_raw/output` reproduced the four
  ordered axes and read a 5,000-gene integer row from tenv, Tahoe bulk and sctox.
- The raw publication downstream gate verified the first and last sealed full
  view rows for real tenv and Tahoe bulk against projected counts and checked the
  control and Delta layers. It resolved the original corpus digest from v1 at
  the new location.
- All 20 focused raw data, prediction-axis and full-view tests passed. The
  Stage 2 path and Git tracking audits passed.
