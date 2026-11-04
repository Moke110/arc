# Frozen v1 view-construction code (archival, not an entry point)

These files are a read-only **copy** of the code that produced the five v1
views in [`../datasets-1/`](../datasets-1/). They are kept here because ADR 0033
decision 8 retires the producer (`stage2/02_cm_data/`) but the frozen
benchmark and the concluded transfer study still name its products as input
evidence: if the code disappeared with the step, nobody could tell how
`dtoxs_cm` / `screen_ctox` / `temp_env_*` were filtered or why a row is in a
view.

Origin: `stage2/02_cm_data/execution/cmdata/{views,paths,qa}.py` at commit
`7e7d689` (2026-11-04), before any ticket-#14 edit of `02_cm_data`.

**Not executable as a step.** Nothing imports these modules, no build wires
them in, and no entry point calls them. `views.py` keeps its intra-package
imports (`from .drugs import ...`, `from .genes import ...`) on purpose — the
package context they were written for is retired, so they are evidence, not
runnable code. Do not "fix" them, do not import them into a build, and do not
edit anything in this directory.

| file | origin | what it holds |
|---|---|---|
| `views.py` | `cmdata/views.py`, verbatim behind a banner | the single view writer: tier filter (`_select_rows`), 18,164-gene projection (`map_counts`), the 24-column `obs` contract, the `log1p_cp10k` layer formula and the per-view `uns["control_mean"]` |
| `paths_v1.py` | `cmdata/paths.py`, view/tier constants only | `EXPECTED_VIEWS`, `EXPECTED_TIER_COUNTS`, `TIER_DATASETS`, `VIEWS`, `TEMP_ENV_VIEWS`, `DEFAULT_BUNDLE`, `CORPUS_VERSION` |
| `qa_tier_nesting.py` | `cmdata/qa.py`, cross-view gates only | `CONTRACT_OBS`, `GateResult`, `cross_view_gates()` — shared var axis, `atc-moa ⊆ atc-moa-mesh ⊆ all-smiles` tier nesting and shared `temp_env` controls |

The `content_sha256` of the five views is unchanged by the move; see
[`../README.md`](../README.md) for the verified table.
