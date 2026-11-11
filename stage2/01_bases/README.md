# Stage 2 shared foundations

This is the standalone foundation retained after the user removed the other
Stage 2 steps. Its public code and required identity tables live inside this
step. Immutable raw inputs remain under `pert_data/`. Deleted historical
experiments and their checkpoints are outside this build.

## Identity

`output/identity-v1/genes.csv` fixes `gene_index` 0–18163 in the original
four-CM-source order. `output/native-v1/tables/genes.csv` retains the exact
original axis bytes; the crosswalk retains source evidence. FAM153B is index
5091: tenv ENSG00000182230 and Tahoe ENSG00000289731/token61482.

Add `execution/` to the import path and use `bases.identity.GeneIdentity`.
`indices(ids, source=...)` resolves native namespaces with explicit ambiguity
errors; `labels` and `native_ids` cross the human/native boundary. The shared
mapping files are derived solely from raw gene metadata and pinned HGNC evidence.

```bash
python stage2/01_bases/execution/build_identity.py
python stage2/01_bases/execution/build_identity.py --verify
```

## Native raw counts

`build_native.py` streams the 16 GSE262419 plate files one plate at a time,
sums duplicate probes by HGNC, and writes `output/native-v1/tenv.h5ad` with
4510 raw wells and 19613 native genes. It retains controls and blanks and has
no normalized expression or significance layers. The manifest seals every
immutable plate, hash table, HGNC table, compound registry and output file.
Independent verification sums the source probes directly for every raw well.

```bash
python stage2/01_bases/execution/build_native.py
python stage2/01_bases/execution/build_native.py --verify
```

The CLI accepts `--data-root`, `--output` and `--record` for worktrees.
`bases.raw.tahoe_axis_mapping(axis, crosswalk, gene_metadata)` exposes the
original token-to-index evidence without loading expression. Tahoe's new raw
preparation boundary projects directly onto the frozen 18164 genes; it does
not require recreating the deleted 62710-gene derived matrix.

## Training files and source replay

`output/training-v1/{tahoe_bulk,tenv}.h5ad` keep separate treated wells. `X`
is natural log1p of raw counts projected onto the full fixed axis and scaled
to 10000. `layers/control` averages independently normalized controls;
`layers/delta` is treated minus control. Tenv uses same-plate DMSO, excludes
Media controls and has cell identity `hiPSC-CM`. Tahoe uses same-plate,
same-cell-line DMSO_TF. Training files contain no truth labels.

Samples retain stable source profile IDs, source salt structures/keys,
stereochemistry-aware normalized parent keys, micromolar doses and exact
matched-control IDs. Parent identity uses RDKit Cleanup, FragmentParent and
Uncharger. The tenv release contains 4206 wells and 1377 parent–dose groups;
81 source-key conversions are explicit. Blank and missing-structure exclusions
are recorded. Group IDs are independent of partition membership.

```bash
python stage2/01_bases/execution/verify_training_source.py --source tenv
python stage2/01_bases/execution/verify_training_source.py --source tahoe
```

Replay joins immutable profile/control IDs, so rebuilding a raw matrix in a
different row order cannot change identity. Every treated/control/delta value
is recomputed from raw counts in bounded windows. Equivalent tautomer SMILES
serializations are audited with mandatory exact InChIKey and stereo-aware
parent-key equality; original training strings remain unchanged. New raw
files have new byte seals and are not claimed byte-identical to deleted files.

`build_corpora.py` builds the pair after both bases-native inputs and the
bounded shared expression release exist. `bases.raw.build`/`verify` own the
shared projection/normalization; explicit source subsets support incremental
preparation. All defaults point inside this step. `training_header` opens
small h5ad metadata only, never truth files.

## Measured truth

`build_significance.py --native <bases-native-tenv.h5ad>` builds separate
parent–dose truth from raw counts on exactly 18164 genes. The bases-owned
DESeq2 worker uses a two-sided Wald test, plate effects for multi-plate
contrasts and within-contrast BH. The partition-authorized reader exposes
only requested groups. Untestable and tested-empty are different states.
Truth fitting and full dual-corpus source acceptance have their own records;
a successful training build alone does not establish either result.

Run heavy commands serially under the repository memory caps. Tests are
individual `unittest discover` files under `execution/tests/`.
