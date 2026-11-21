# Separate the Stage 2 drug representation comparison

**Status:** accepted design (2026-11-21). Implementation pending.

Insert `03_tsg_drug_emb` before dose and transfer, which move to
`04_tsg_dose` and `05_tsg_transfer`. This gives the representation
comparison its own accepted products and report while keeping the dose and
transfer results independently reproducible. Migrate active references and
ignored payload locations with an explicit path mapping for sealed evidence.

Compare representations on complete tenv cold-drug folds, holding the dose
method fixed to FiLM. Vary only the input width of the first drug-encoder
layer to accept each representation at its native width. Its output stays
400-dimensional and the rest of the network shares architecture and
initial weights across arms. This avoids introducing a separate projection
or compression rule for each representation while retaining the existing
FiLM interface.
