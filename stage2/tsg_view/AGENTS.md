# tsg_view — execution workflow

Read [README.md](README.md), the [Stage 2 routing](../AGENTS.md), the
[step contract](../AGENTS.md#step-contract), [ADR 0082](../../docs/adr/0082-view-variants-tenv-tahoe-replicate-structure.md),
[ADR 0077](../../docs/adr/0077-five-k-only-axis-local-normalization.md) and
[ADR 0078](../../docs/adr/0078-dual-source-axis-local-dose-round1.md).

This is a **non-numbered study step**: it publishes comparison evidence and **no model
weights**, and it must not edit Steps 01–06 products or their sealed records.

## What to reuse, and what stays independent

- Read the integer corpus, the 5k gene axis and the accepted `view-5k` release through the
  Base step's public API (`stage2/01_bases/output/api/bases/{corpus,gene_axis,view,release}.py`).
  Normalize only with `bases.view.library_totals` / `bases.view.log1p_cp10k`.
- Reuse the sealed Step 04 folds (`tsg_round1.protocol.make_splits`), the sealed freeze's
  `splits.csv`/`parents.csv`/`morgan.npy` and the native `no-dose` learner mechanisms.
- Own the view product manifests (`tsg-view-product-v1`), the view-aware freeze
  (`tsg-view-freeze-v1`), the view-aware verifier (`execution/inputs.py`), the frozen truth
  and diagnostics (`tsg-view-truth-v1`) and the run manifests (`tsg-view-run-v1`). Do not
  route a derived view through `bases.view.verify` or `tsg_round1.inputs.Inputs`; both
  hard-assert the released full population.

## Build, freeze, truth, matrix

One heavy job at a time, under the repository caps (10 GiB memory, 8 GiB swap, 10-second
RAM/GPU telemetry). Use `../01_bases/execution/capped.sh` or `systemd-run --user --scope -p
MemoryMax=10G -p MemorySwapMax=8G`; check `free -g` (≥ 12 GiB available) first and never
run a fit beside `script/verify_all.sh`. Disposable payloads go under
`execution/checkpoints/` or `execution/freezes/`; keep manifests, parents, splits and CSV
summaries tracked.

```bash
PY=/home/chang/miniconda3/envs/arc/bin/python
cd /home/chang/arc/stage2/tsg_view

$PY execution/build.py --view <full|merge|solo|solo-deepest|full-pooled-label> \
    --source <tenv|tahoe_bulk> --output output/views --staging execution/checkpoints \
    --repo-root /home/chang/arc --record log/build-0X.json
$PY execution/freeze.py --view <view> --source <source> --views output/views \
    --destination execution/freezes/<view>/<source> --repo-root /home/chang/arc
$PY execution/truth.py --source all --output output/truth --repo-root /home/chang/arc
$PY execution/formal.py plan   --candidate execution/checkpoints/formal
$PY execution/formal.py run    --candidate execution/checkpoints/formal   # 45 fits, resumable
$PY execution/formal.py verify --candidate execution/checkpoints/formal   # writes formal_verify.json
```

`build.py` auto-discovers `execution/views/<view>.py::build(context)`; `freeze.py` and
`execution/inputs.py` cross-check the sealed `(parent, dose) → partition` mapping so a view
cannot move folds. `formal.py` seals its plan before the first fit and re-execs under the
capped scope; re-running the same command resumes only the missing phases, and
`--only-source` / `--only-fold` keep each fold group whole.

## Report, figures and acceptance

```bash
$PY execution/report.py --runs output/runs --truth output/truth \
    --report output/report --figure-data figure/data
$PY figure/script/render.py --data figure/data --output figure --suite tenv   # pass one
$PY figure/script/render.py --data figure/data --output figure --suite all    # pass two
```

The renderer reads the `figure/data/*.csv` only (no matrices, checkpoints or truth) and is
deterministic; the tenv suite is published first and the Tahoe/cross-source panels are added
afterwards. Keep `execution/render.py`, `execution/report.py` and
`figure/script/render.py` as the only report code paths.

Re-run the sealed receipt from the materialised checkout, **alone** (it hashes the accepted
payloads and the 49 GB Steps 01–06 artifact inventory):

```bash
$PY execution/acceptance.py receipt --root /home/chang/arc --phase final --verify-binary \
    --preregistration log/acceptance/pre_registration.md \
    --health-dir log/acceptance/health --output log/acceptance/receipt.json
```

The receipt fails closed on a missing artifact, a digest mismatch, a cap violation, a
non-identical CSV-only re-render, a missing replay token or health transcript, or a changed
Steps 01–06. Its Steps 01–06 proof reads the pinned per-path baseline
`.scratch/tsg-view/research/baseline/steps01-06-per-path.sha256` (11,468 files), the planning
ref in `refs.txt` and the dirty working tree; all three must be clean.

## Tests and health

```bash
cd stage2/tsg_view
$PY -m unittest discover -s execution/tests -p test_acceptance.py -v
$PY -m unittest discover -s execution/tests -p test_docs.py -v
$PY -m unittest discover -s execution/tests -p 'test_*.py' -v

cd /home/chang/arc
bash script/check_paths.sh
bash script/check_tracked.sh
bash script/verify_all.sh     # alone; heavy (52 GB hashed, 9.7 GiB page cache)
```

The tests are hermetic (synthetic fixtures, no payloads, GPU or network). Health checks run
from the repo root: a git worktree lacks the gitignored external payloads, so
`check_paths.sh` fails there by design and `verify_all.sh` must run in the main checkout.

## Staging and publication

Publish only the minimal accepted products under `output/report/`, `figure/` and the view
manifests under `output/views/`; payloads, freezes and run checkpoints stay gitignored.
Commit acceptance evidence under `log/acceptance/` and the ticket evidence under
`log/build-08.json`. Leave the sealed Step 01–06 products and their logs intact.
