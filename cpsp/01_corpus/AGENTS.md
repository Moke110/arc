# CPSP 01_corpus agent instructions

Read the [step README](README.md) for the Heliste 2021 baseline and H-R source contract. Use [ADR 0044](../../docs/adr/0044-stage3-source-first-corpora.md) for corpus decisions and the [CPSP release boundary](../../docs/adr/0106-cpst-and-cpsp-public-release-boundaries.md) for migration of sealed artifacts.

## Build and verify

- **Inputs.** In the `arc` environment, use the dependency versions recorded in the [corpus manifest](output/manifest.json). Stage the full frozen Heliste `source/` inventory in the checkout where you run the workflow. The [source manifest](../../pheno_data/Heliste-2021/manifest.sha256) pins Dataset S1 and its companion files. The [author log](../../pheno_data/Heliste-2021/author_log.json) declares the tracked local evidence, including the v1 reference and v2 desalted-parent snapshots. Use the [phenotype archive workflow](../../pheno_data/README.md) for acquisition and source verification.
- **Execution.** Run `python cpsp/01_corpus/execution/corpus.py build` from the repository root. The [offline CLI](execution/corpus.py) stages a new release and verifies an existing destination in place. Pass `--output` with a fresh directory for an intentional rebuild.
- **Acceptance.** Run `python cpsp/01_corpus/execution/corpus.py verify` to check the release inventory and hashes and compare both arm exports with an independent offline rebuild. For changes to the build or verification contract, run `python -m unittest discover -s cpsp/01_corpus/execution/tests -v` ([tests](execution/tests/test_corpus.py)). For public path changes, run `bash script/check_cpsp_paths.sh` ([path audit](../../script/check_cpsp_paths.sh)).

## Corpus invariants

Keep baseline and H-R separate. Account for every Dataset S1 cell once in either its arm's corpus or pending-identity report. Preserve the source label, dose, cell reference, original readout, and null replicate even when labels share a desalted parent. Preserve the versioned magnitude-only direction rule and its same-arm vehicle normalization assumption in the annotation. Follow the [Morgan step](../02_base_morgan/README.md) for downstream preparation and evaluation.

Keep the [superseded release record](execution/historical/release_record.json) and its local payloads sealed.
