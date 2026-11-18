# log/ — bases acceptance evidence

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
| `tox_candidate_verify.json` | ADR 0072 #216 staged twelve-dataset TOX measurement candidate: fresh source replay, counts, digests and resource peaks; not the accepted successor release |
| `disease_fusion_verify.json` | historical ADR 0071 six-disease fusion replay; superseded by current full-source acceptance |
| `gene_axis_verify.json` | `build_gene_axis.py --verify` replay output |
| `view_publish_verify.json` | `build_view.py publish` result (materialises the gitignored view payloads; manifests byte-identical to the tracked ones) |
| `view_verify.json` | historical two-view replay; superseded by the historical ADR0073 four-view receipt |

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

## ADR 0072 coherent acceptance

`successor_acceptance.json` records #219 final source/integrated replays, all
affected public-contract tests, three offline health checks, exact human artifact
reproduction and main-checkout byte preservation. The successor sections of
`acceptance.md` and `inputs.json` supersede historical six-category records.

## Historical ADR 0073 coherent acceptance

- `full_source_acceptance.json`: #225 complete five-category/fusion/four-axis/eight-payload release, exact manifests and unchanged corpus bytes, affected tests, offline health and source/payload preservation.
- `source_collection_candidate.json` and `tier_first_fusion_candidate.json`: full-source and tier-first fusion candidate evidence from #221/#222.
- `source_derivation_replay.json`: sixteen genuine per-item refits, current R-worker pins, identical scientific CSV/sample bytes and deterministic gzip intermediates.
- `four_axes_candidate.json`: #223 independent upstream and per-view replays, bounded public-reader checks, exact axes/nesting/fallback, payload map and recovery evidence.
- `gene_axis_document_verify.json`: root DOCX and both figure hashes, canonical input/tool pins, deterministic byte replay and every one of thirteen inspected pages.

The older entries/sections above describe their own dated releases. They keep
historical names and counts; current contracts/paths are in `../README.md`.

## Current ADR0077 expression acceptance

- `axis_local_acceptance.json`: paired 5k-local release, real independent replay,
  selected totals, populations, exact preservation and new document QA.
- `axis_local_preservation.json`: streamed unchanged upstream and accepted
  downstream scientific/frozen-input/weight/result/sealed-log inventory.
- `axis_local_candidate.json`: independently built and replayed #244 candidate.

Earlier receipts remain byte-identical history. Disease and all four gene axes
remain the ADR0073 scientific generation; only expression normalization changes.
