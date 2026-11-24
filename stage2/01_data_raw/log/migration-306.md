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
