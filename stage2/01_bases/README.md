# Stage 2 shared foundations

ADR 0061 assigns shared gene identity, training corpora, HVG, metrics and report
rendering to this boundary. The identity release is the first implemented
interface; existing consumers remain operational while the later tickets migrate.

`output/identity-v1/genes.csv` adds `gene_index` 0–18163 to the exact frozen
four-CM-source axis, without sorting. Its manifest pins the original axis,
source crosswalk, HGNC aliases and independently replayed native mappings.
FAM153B has index 5091: tenv ENSG00000182230 and Tahoe ENSG00000289731/token
61482 retain their own source evidence.

## Identity interface

Add `execution/` to the caller's import path and import `GeneIdentity` and
`native_mapping_rows` from `bases.identity`. `GeneIdentity.load(axis, crosswalk,
annotation, native_mappings=native_mapping_rows(shared))` validates the pinned
axis bytes and size. `indices(ids)` resolves global HGNC/canonical/alias IDs;
`indices(ids, source=...)` resolves native IDs only within that source.
`labels(indices, field='symbol')` converts integers for display or HGNC-native
boundaries; `native_ids(index, source)` returns all evidenced source identities.
`temp_env` and `screen_ctox` evidence is exposed as the corpus names `tenv` and
`sctox`; original evidence records retain their source names. Unknown, off-axis
and ambiguous identities raise `IdentityError`. Duplicate corroborating evidence
is safe; competing identities, including off-axis candidates, fail explicitly.

Build and replay from the repository root:

```bash
python stage2/01_bases/execution/build_identity.py
python stage2/01_bases/execution/build_identity.py --verify
python -m unittest discover -s stage2/01_bases/execution/tests -p test_identity.py
```

The builder accepts explicit `--native`, `--annotation`, `--shared` and `--output`
locations for isolated worktrees. It loads only annotation/crosswalk tables,
never expression matrices. The output stores the indexed axis and hashes rather
than duplicate source data.

## Compatibility locations

- Native model setup and pinned environments: [01_setup](../01_setup/README.md)
  and [_envs](../_envs/README.md).
- Native corpora and bounded shared expression build/replay:
  [02_dataset_corpus](../02_dataset_corpus/README.md),
  `execution/cmdata/shared_cm.py` and `execution/build_corpus.py --shared-cm`.
- Existing significance/retrieval public entry points:
  [evaluation](../evaluation/README.md). Its old labels remain historical;
  this identity ticket does not assert they satisfy the new axis-only family.

These are live compatibility entry points, with their original code and seals
unchanged. New training corpora, new true labels and consumer migration are
separate acceptance boundaries under #145.
