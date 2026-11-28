# AGENTS.md — CarPS-Transcriptome shared utilities

Add `cpst/util/output/api` to the Python import path and reuse `cpst_util`
through the [public API](README.md). Follow the parent
[CarPS-Transcriptome workflow](../AGENTS.md).

- **Identity and chemical inputs.** Read the frozen gene order through
  `load_gene_table()` and check its [manifest](output/gene-identity/manifest.json)
  when changing the identity reader. Use `RepresentationTable` with the
  matching `tenv` or `tahoe_bulk` release. Use the installed source-specific
  manifest and its payload, or pass their directory explicitly. Verify parent
  structures with `verified_parent_keys` and select vectors in caller order,
  preserving repeated conditions. Use `open_ecfp4` for the binary fingerprint
  release, with its source manifest when verifying a publication.
- **Training and evaluation.** Start with the numbered study runner and use
  the [shared runtime APIs](README.md#transigen-training-and-evaluation).
  Materialize the study's pinned training, axis, fold, chemical and source inputs.
  Bind new registrations to the scientific configuration and current runtime.
  Keep fold validation, calibration, aggregation and report tables with the study.
  Use the original source revisions for sealed registrations and preserve the
  recorded files. Follow the receipt and replay routes for
  [Steps 04–05](training/verification.md) and [Steps 07–10](training/cpst-verification.md).
- **Derivation and publication.** Reuse `cpst_util.derivation` for DEG, ranking,
  R-worker resolution and manifest verification. Build and verify a candidate
  outside the accepted output before replacing it under the destination lock.
  See the [accepted-product contract](../../docs/adr/0062-stage2-accepted-products-and-reproducible-figures.md).
- **Figures.** Use the [comparison CSV schema](README.md#comparison-figures)
  for plot inputs and the [plotting API](README.md#plotting) for shared marks.
  Run the [current figure verification](../README.md#verify-the-current-figures)
  after changes to shared rendering. The gate checks receipt hashes and both
  independent renders from the published plot data.

## Verify the utility APIs

Run from the repository root in the `arc` environment with the pinned Git
revisions available. Keep test scratch and the Matplotlib cache under `tmp/`.
Run the shared tests before publishing an API change.

```sh
mkdir -p "$PWD/tmp"
UTILITY_TEST_WORK=$(mktemp -d "$PWD/tmp/cpst-util-tests.XXXXXX")
export TMPDIR="$UTILITY_TEST_WORK"
export MPLCONFIGDIR="$UTILITY_TEST_WORK/matplotlib"
PYTHONPATH=cpst:cpst/util/output/api python -m unittest discover -s cpst/util/execution/tests
```

Run the affected training or evaluation test groups one at a time. For example

```sh
PYTHONPATH=cpst:cpst/util/output/api python -m unittest discover \
  -s cpst/05_tsg_module/execution -p test_module_objective.py
python -m unittest discover -s cpst/util/eval/tests/04_tsg_axis
```
