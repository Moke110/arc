# tsg_view — execution workflow (first draft; ticket 08 finalizes)

Read [README.md](README.md), the [Stage 2 routing](../AGENTS.md), [ADR 0082](../../docs/adr/0082-view-variants-tenv-tahoe-replicate-structure.md),
[ADR 0077](../../docs/adr/0077-five-k-only-axis-local-normalization.md) and
[ADR 0078](../../docs/adr/0078-dual-source-axis-local-dose-round1.md).

This is a **non-numbered study step**: it publishes comparison evidence and no model
weights. It must not edit Steps 01–06 products or their sealed records.

## What to reuse, and what stays independent

- Read the integer corpus, the 5k gene axis and the accepted `view-5k` release through the
  Base step's public API (`stage2/01_bases/output/api/bases/{corpus,gene_axis,view,release}.py`).
  Normalize only with `bases.view.library_totals` / `bases.view.log1p_cp10k`.
- Reuse the sealed Step 04 folds (`tsg_round1.protocol.make_splits`), the sealed freeze's
  `splits.csv`/`parents.csv`/`morgan.npy` and the native `no-dose` learner mechanisms.
- Own the view product manifests (`tsg-view-product-v1`), the view-aware freeze
  (`tsg-view-freeze-v1`), the view-aware verifier (`execution/inputs.py`) and the run
  manifests (`tsg-view-run-v1`). Do not route a derived view through `bases.view.verify`
  or `tsg_round1.inputs.Inputs`; both hard-assert the released full population.

## Build and verify

- `execution/build.py --view <view> --source <tenv|tahoe_bulk> --output <views root>
  --staging <dir> [--repo-root <materialized checkout>]` auto-discovers
  `execution/views/<view>.py::build(context)` and installs the product under
  `<views root>/<view>/`.
- `execution/freeze.py` converts one product into a Step 04-shaped freeze under
  `execution/freezes/<view>/<source>/`; `execution/inputs.py::Inputs` verifies it and
  cross-checks the sealed Step 04 fold mapping.
- Run every build/freeze/fit singly under the repository caps
  (`../01_bases/execution/capped.sh /usr/bin/time -v …`), with 10-second RAM/GPU sampling.
  Write disposable payloads beneath `execution/checkpoints/` or `execution/freezes/`;
  keep manifests, parents, splits and CSV summaries tracked.

## Tests and health

- `cd stage2/tsg_view && /home/chang/miniconda3/envs/arc/bin/python -m unittest discover
  -s execution/tests -p 'test_*.py' -v` — hermetic, no real payloads.
- Health from the repo root: `bash script/check_paths.sh` and `bash script/check_tracked.sh`.
- Leave the sealed Step 01–06 products and their logs intact.
