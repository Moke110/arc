# Stage 2 shared utilities

Add `stage2/util/output/api` to the Python import path. The public
`stage2_util` package provides source-specific drug identities and representations,
indexed gene identities, seed selection, generic metrics, CSV rendering and derivation utilities.

| Module | Entry point |
| --- | --- |
| `identity` | `load_gene_table()` verifies and reads the accepted 18,164 rows in index order. `GeneIdentity` resolves global and source-specific identifiers. |
| `representations` | `RepresentationTable(source="tenv")` or `RepresentationTable(source="tahoe_bulk")` verifies the source manifest, reads the parent identity table, and selects vectors by ordered parent InChIKeys. `ordered_parents(observations)` and `verified_parent_keys(parents)` validate profile and accepted representation identities. |
| `ecfp4` | `open_ecfp4(path)` verifies the accepted binary fingerprint and parent order. `publish_ecfp4(source, candidate)` preserves frozen input provenance. |
| `seeds` | `read_seeds()` returns the accepted seed tuple. Pass a count to select a prefix. |
| `metrics` | `vector_metrics(truth, prediction)` scores paired finite vectors, including the accepted constant-vector policy. |
| `plotting`, `palette` | `read_rows(path)`, `save_figure(figure, path, **kwargs)` and the Stage 2 color and role APIs support rendering solely from plotted CSV values. |
| `publication` | `publication_lock(destination)` and `replace_release(...)` coordinate verified replacement of accepted releases. |
| `derivation` | `bulk_deg`, `disease_deg`, and `deg_rank` provide common DEG, rank, and worker interfaces. `r_env` locates R resources and `release` provides derivation manifest and candidate verification. |

The default readers open the accepted gene and feature products under util.
For a separate materialized checkout, pass the directory explicitly to `load_gene_table` or `RepresentationTable`. The
feature reader checks the shipped source-specific manifest pin unless an
explicit expected manifest digest is provided. The manifest binds the frozen
input generation, extraction provenance, parent structures and vector files.
Studies retain fold coverage verification, aggregation and report tables.

From the repository root:

```sh
PYTHONPATH=stage2/util/output/api python -m unittest discover -s stage2/util/execution/tests
```
