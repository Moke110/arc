# Stage 2 source-native corpus

Read ADR 0083 and the [step README](README.md). Owns immutable, source-granular
raw integer counts, observation metadata, matched vehicles and gene coverage.
Retain unmatched observations. Paired expression, dataset exclusions from model
views, axis selection and normalized matrices belong to their consuming steps.

- `output/api/native_corpus.py` provides the bounded corpus API.
- `execution/build_corpus.py` publishes and independently replays candidates.
- `execution/checkpoints/` holds disposable build state.
- `log/` records build and replay acceptance evidence.

Follow the repository resource caps for Tahoe source scans. Never load the complete
matrix. Run one heavy job at a time. External source files are read-only.
