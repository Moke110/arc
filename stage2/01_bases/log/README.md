# log/ — bases acceptance evidence (post-0066)

ADR 0066 applies the repository rule strictly here: `log/` is execution-run
intermediate output and verification evidence, **not** a history of directory
changes (git is that history). The pre-0066 record tree was deleted at the rebuild
and is not archived; this directory was regenerated at acceptance (#192). The one
restored exception is `raw_counts_provenance.json` (#194): the #186 rehome ruling
requires the recorded source-replay evidence to travel with the regenerated record.

| path | contents |
|---|---|
| [`acceptance.md`](acceptance.md) | the acceptance record: per-product replay results and digests, resource peaks, "no retired path read" evidence, old-to-new mapping, recorded errata and the #173 decision answers |
| [`inputs.json`](inputs.json) | machine-readable input/output digests per product, replay commands and resource measurements, the retired evidence digests preserved across deletion, the figure-reproduction result, the #208 v2 republish gate (`v2_republish`) and every superseded digest (`superseded_digests`) |
| [`raw_counts_provenance.json`](raw_counts_provenance.json) | the two pre-0066 source-replay acceptance records and the native tenv build manifest restored verbatim from git history, digest-pinned by the raw-count release manifest |
| `raw_counts_verify.json` | `build_raw_counts.py verify` replay output (manifest digest + counts) |
| `corpus_verify.json` | `build_corpus.py verify` replay output |
| `disease_deg_verify.json` | `build_disease_deg.py verify --all` replay output (the four registry categories: nHF, ICM, AF, genetic_cardiomyopathy) |
| `cardiotox_verify.json` | `build_disease_deg_cardiotox.py verify` replay output (three products) |
| `disease_fusion_verify.json` | `build_disease_fusion.py verify` replay output (six diseases, v2 DEG2000) |
| `gene_axis_verify.json` | `build_gene_axis.py --verify` replay output |
| `view_publish_verify.json` | `build_view.py publish` result (materialises the gitignored view payloads; manifests byte-identical to the tracked ones) |
| `view_verify.json` | `build_view.py verify` replay output for both views |

Disposable build state lives in `execution/checkpoints/`, never inside `output/`.

## Consequential removals

The #192 rebuild removed two QA helpers that had outlived the evidence they
reported. Neither removal is a new science decision: both produced the retired
HVG population-isolation report, and both lost their inputs with the evaluation
chain (recorded as ADR 0068 (c)).

| removed | commit | why | where the evidence lives now |
|---|---|---|---|
| `execution/verify_hvg_metrics.py` (with `hvg_metrics.md`) | `ee88e3c` (2026-11-13) | it recomputed the retired population-isolation metrics from `bases.metrics`; after the evaluation chain was removed it had no surviving input or consumer | listed in ADR 0066 *Removals* and named as "the already-retired `verify_hvg_metrics` report" in ADR 0068 (c); the retained `bases.hvg` release tests in `execution/tests/test_hvg.py` are the surviving regression coverage; code in git history |
| `output/api/bases/hvg_qa.py` | `6078e9f` (2026-11-13) | it imported `score_groups`/`summarize_runs` from the deleted `bases.metrics`, so it could not import once the evaluation chain was removed; it was the last reporter of that retired QA | ADR 0068 (c); its population-isolation assertion was dropped from `test_hvg.py` (the release tests keep their two cases); code in git history |

No other #192 removal is consequential: the one-off migration/audit scripts and
the superseded build chain had finished their job, and the retired truth/
evaluation products have no surviving consumer (the old→new mapping is in
[`acceptance.md`](acceptance.md)).
