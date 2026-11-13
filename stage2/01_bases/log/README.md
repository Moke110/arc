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
| [`inputs.json`](inputs.json) | machine-readable input/output digests per product, replay commands and resource measurements, the retired evidence digests preserved across deletion, and the figure-reproduction result |
| [`raw_counts_provenance.json`](raw_counts_provenance.json) | the two pre-0066 source-replay acceptance records and the native tenv build manifest restored verbatim from git history, digest-pinned by the raw-count release manifest |
| `raw_counts_verify.json` | `build_raw_counts.py verify` replay output (manifest digest + counts) |
| `corpus_verify.json` | `build_corpus.py verify` replay output |
| `disease_deg_verify.json` | `build_disease_deg.py verify --all` replay output (four categories) |
| `cardiotox_verify.json` | `build_disease_deg_cardiotox.py verify` replay output (three products) |
| `gene_axis_verify.json` | `build_gene_axis.py --verify` replay output |
| `view_publish_verify.json` | `build_view.py publish` result (materialises the gitignored view payloads; manifests byte-identical to the tracked ones) |
| `view_verify.json` | `build_view.py verify` replay output for both views |

Disposable build state lives in `execution/checkpoints/`, never inside `output/`.
