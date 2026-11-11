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

## Indexed training corpora

`output/training-v1/{tahoe_bulk,tenv}.h5ad` are separate treated-well
releases. `X` contains full-axis treated expression; `layers/control` is the
mean of independently normalized matched controls, and `layers/delta` is
`X - control`. `var/gene_index` is the fixed 0-based axis. No statistical
truth or significance labels are present. Tahoe rows are well-by-cell-line
pseudobulks; tenv's `cell_identity` is `hiPSC-CM`.

`obs/sample_id` is stable and source-prefixed. Source row/profile, plate, well,
administered compound name/structure, full stereochemistry-retaining desalted
parent InChIKey, micromolar dose, native dose value/unit and conversion evidence
are retained. Parent identity uses RDKit Cleanup, FragmentParent and Uncharger,
with stereochemistry retained. `drug_parent` is the normalized identity;
`source_parent_key` and the original `drug_key` retain the native identity,
including its original protonation state. A resolvable identity conversion is
recorded rather than excluding the treatment. `dose_original_value` refers to the pinned native corpus input;
its upstream `dose_source` is retained when available. Native corpora and their
raw-source provenance remain unchanged. Unknown Tahoe time/vehicle metadata is
not invented. Tenv parent–dose keys with missing or incompatible time/vehicle
conditions fail before publication. `group_id` hashes parent and dose, and is
not a partition assignment: membership belongs to an explicit partition/cohort.

`matched_control_rows` and `matched_control_ids` are exact JSON lists in native
row space. Tenv uses same-plate **DMSO** (Media controls are excluded); Tahoe
uses same-plate, same-cell-line **DMSO_TF**. The exclusion CSVs account for
technical blanks and excluded treatments with explicit reasons. Source control
profiles stay in the native corpus and are referenced rather than published as
training observations.

Run serially under the execution-environment memory caps:

```bash
python stage2/01_bases/execution/build_corpora.py --record stage2/01_bases/execution/build_record/training_v1_build.json
python stage2/01_bases/execution/build_corpora.py --verify --record stage2/01_bases/execution/build_record/training_v1_verify.json
python -m unittest discover -s stage2/01_bases/execution/tests -p test_corpora.py
```

The CLI accepts explicit `--native`, `--shared`, `--drugs` and `--output` paths
for worktrees. It verifies the compatible `shared-cm-v1` input/payload hashes
and completely replays its raw-count projection before reuse. The independent
release verifier recomputes every included treated profile and matched-control
mean from native raw counts, then checks all three output arrays. Bounded
64-row windows avoid loading the full observation-by-axis matrix. The native
expression layers are never used. The manifest binds source evidence, metadata,
exclusions and output file hashes. `bases.corpora.training_header` reads small
training metadata using h5py; it does not discover or open a labels release.

Historical inputs/results retain their existing normalization semantics; this
release does not retrain or retroactively relabel them.

The real build (not yet accepted) retains 65,318 Tahoe treated pseudobulks and all 4,206 tenv
treated wells. Tahoe excludes 300 wells without resolvable source structures;
tenv excludes its 32 technical blanks. All 1,400 Tahoe DMSO_TF controls and 192
tenv DMSO controls contribute to matched means (80 tenv Media controls do not).
The final tenv parent normalization records 81 native-key conversions and 1,377
parent–dose groups. The conversion table is
`execution/build_record/tenv_parent_identity.csv`; each h5ad also preserves the
per-well source and conversion evidence.

Acceptance is blocked: a subsequent user instruction removed the native/shared
Stage 2 inputs during the final source replay. The build completed, but no
full-release PASS is claimed. See
`execution/build_record/training_v1_acceptance_blocked.json` and the exact
verification failure record. Deleted resources have not been restored.
