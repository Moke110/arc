# 07_lxpert — whole-transcriptome Latent XPert

Step 07 trains a transcriptome autoencoder and uses the refit checkpoint in a
Latent XPert treatment smoke. The development split and full-data refit serve
reconstruction. The treatment branch exercise checks integration.
[ADR 0051](../../docs/adr/0051-latent-xpert.md)

## Expression input release

From the repository root, use the `arc` Python environment. Run each real-data
command serially under the [repository resource limits](../../docs/agents/execution-environment.md).

```bash
python stage2/07_lxpert/execution/script/lxpert.py prepare
python stage2/07_lxpert/execution/script/lxpert.py verify
```

The default source is the Step 02 v2 corpus. `prepare --corpus-dir PATH` selects
another copy of the same corpus and checks it against the frozen source digests.
Both commands take `--release PATH`, which defaults to `output/input/` within
this step. Prepare requires a new destination and publishes the directory after
the build succeeds. Verify checks the release files, expression values, masks,
memberships and declared source hashes. Add `verify --corpus-dir PATH` to rehash
the original files as well.

The release includes 67,018 Tahoe and 4,478 tenv profiles. The 32 tenv technical
blanks are excluded. Raw counts are projected onto the ordered 18,164-gene CM
axis, then normalized by the projected library size with CP10K and `log1p`.
Tahoe maps by Ensembl ID and tenv by HGNC ID. Tahoe's unobserved FAM153B is stored
as zero and masked from the Tahoe reconstruction loss. Tenv keeps the measured
gene. Both sources use the complete ordered decoder axis.

For each source, `prepare` ranks plate IDs by a hash of the seed, source and
plate. The first fifth of the plate list, rounded down with a minimum of one
plate, is held out for development validation. The refit uses every released
profile. The default seed is 20261109.

| Release file | Contract |
| --- | --- |
| `genes.csv` | Frozen ordered CM gene axis |
| `{tahoe,tenv}_expression.npy` | Row-major float32 matrices readable via `np.load(..., mmap_mode="r")` |
| `{tahoe,tenv}_target_mask.npy` | Boolean arrays in gene order. True identifies observed reconstruction targets |
| `{tahoe,tenv}_profiles.csv` | One row per expression row, with the original zero-based HDF5 row, profile ID, plate, control flag and development role |
| `manifest.json` | Source, axis, file and membership SHA-256 digests, plus normalization and population counts |

The [verified input release record](execution/input_record/release.json) captures
the published manifest and resource measurements for the real-data build and
verification.

Run bounded CLI tests with
`python -m unittest discover -s stage2/07_lxpert/execution/tests -p test_inputs.py`.
