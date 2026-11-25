# Stage 2 raw data

Read ADRs 0083 and 0085 and the [step README](README.md). Owns immutable, source-granular
raw integer counts, observation metadata, matched vehicles and gene coverage.
Retain unmatched observations. The training step owns pooled expression,
Delta and training exclusions. The prediction-axis step owns gene selection.

- `output/api/native_corpus.py` provides the bounded corpus API.
- `output/manifest.json` inventories seven per-dataset releases beside `output/api/`.
- `execution/build_corpus.py` publishes and independently replays candidates.
  Use `migrate --predecessor <sealed-v1> --output <candidate>` to rehearse
  migration away from the accepted output.
- Pass `--predecessor` to replay an archived migration explicitly.
- `execution/checkpoints/` holds disposable build state.
- `log/` records build and replay acceptance evidence.

Follow the repository resource caps for Tahoe source scans. Never load the complete
matrix. Run one heavy job at a time. External source files are read-only.
