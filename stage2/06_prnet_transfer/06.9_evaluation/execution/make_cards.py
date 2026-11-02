#!/usr/bin/env python
"""Generate the 30 model cards with the ticket-#4 writer (never a local copy).

    python execution/make_cards.py --phase test

For each evaluated unit the card is ``ptlib.modelcard.write_model_card(run_dir,
eval_dir=<unit eval dir>, out_path=output/model_cards/<unit_id>.md)``: the card already
carries data identity, structure/recipe, seed, init provenance, stop reason, validation
curve summary and the held-out evaluation block, so ticket #9 only has to point the writer
at the right directories and record the resulting hashes.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import common  # noqa: E402
from common import DEFAULT_OUT, EvaluationError, load_json, now_utc  # noqa: E402
from units import load_unit_manifest  # noqa: E402

SCHEMA = "stage2/06.9_evaluation/model_cards/v1"


def _modelcard_module():
    if str(common.PIPELINE_EXECUTION) not in sys.path:
        sys.path.insert(0, str(common.PIPELINE_EXECUTION))
    from ptlib import modelcard  # noqa: WPS433 - the pipeline's writer

    return modelcard


def _relativize_card(path: Path) -> None:
    """Strip any checkout prefix from a generated card (ADR 0030 §4).

    The shared pipeline writer renders the run manifest's ``command`` and ``run_dir``
    verbatim, and #6/#7's manifests were produced in their own per-ticket worktrees.
    Those checkout prefixes are meaningless from a clone, so they are rewritten
    repo-relative here; no other line of the card moves. The longest root is replaced
    first, so a worktree root never degrades into the main checkout's suffix.
    """
    roots = sorted({str(root) for root in common.repo_roots()}, key=len, reverse=True)
    text = path.read_text(encoding="utf-8")
    normalized = text
    for root in roots:
        normalized = normalized.replace(root + "/", "")
    if normalized != text:
        common.write_text_atomic(path, normalized)


def build_cards(
    *,
    units_path: Path,
    out: Path,
    phase: str,
    cards_dir: Path,
    allow_subset: bool,
) -> dict:
    payload = load_unit_manifest(units_path, require_complete=(phase == "test" and not allow_subset))
    modelcard = _modelcard_module()
    cards_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    missing: list[str] = []
    for unit in payload["units"]:
        unit_id = unit["unit_id"]
        eval_dir = out / "units" / unit_id / "eval"
        run_dir = common.repo_path(unit["run_dir"])
        if not (eval_dir / "aggregate.json").exists():
            missing.append(unit_id)
            continue
        if not (run_dir / "manifest.json").exists():
            raise EvaluationError(f"unit {unit_id}: candidate run dir vanished: {run_dir}")
        target = cards_dir / f"{unit_id}.md"
        modelcard.write_model_card(run_dir, eval_dir=eval_dir, out_path=target)
        _relativize_card(target)
        rows.append(
            {
                "unit_id": unit_id,
                "phase": phase,
                "fold": int(unit["fold"]),
                "seed": int(unit["seed"]),
                "init": unit["init_name"],
                "lr": float(unit["lr"]),
                "run_dir": common.relpath(run_dir),
                "eval_dir": common.relpath(eval_dir),
                "card": target.name,
                "card_sha256": common.sha256_file(target),
                "card_bytes": target.stat().st_size,
            }
        )
    if missing and not allow_subset:
        raise EvaluationError(f"no evaluation output for {len(missing)} unit(s): {missing[:5]}")
    index = {
        "schema": SCHEMA,
        "phase": phase,
        "created_utc": now_utc(),
        "generator": "stage2/06_prnet_transfer/06.4_pipeline/execution/ptlib/modelcard.py",
        "units_manifest": common.relpath(units_path),
        "units_manifest_sha256": common.sha256_file(units_path),
        "n_cards": len(rows),
        "n_missing": len(missing),
        "missing": missing,
        "cards": rows,
    }
    common.write_json(cards_dir / "index.json", index)
    header = "unit_id,phase,fold,seed,init,lr,card,card_sha256,card_bytes\n"
    lines = [header.rstrip("\n")]
    for row in rows:
        lines.append(
            f"{row['unit_id']},{phase},{row['fold']},{row['seed']},{row['init']},"
            f"{row['lr']},{row['card']},{row['card_sha256']},{row['card_bytes']}"
        )
    common.write_text_atomic(cards_dir / "index.csv", "\n".join(lines) + "\n")
    return index


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate the ticket-#9 model cards from the pipeline writer.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--units", type=Path, default=DEFAULT_OUT / "record" / "units.json")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--phase", choices=("valid", "test"), default="test")
    parser.add_argument(
        "--cards-dir", type=Path, default=None,
        help="default output/model_cards for the test phase, output/dry_run/model_cards otherwise",
    )
    parser.add_argument("--allow-subset", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    ns = build_parser().parse_args(argv)
    cards_dir = Path(ns.cards_dir) if ns.cards_dir is not None else (
        DEFAULT_OUT / ("model_cards" if ns.phase == "test" and not ns.allow_subset else "dry_run/model_cards")
    )
    try:
        index = build_cards(
            units_path=ns.units,
            out=ns.out,
            phase=ns.phase,
            cards_dir=cards_dir,
            allow_subset=ns.allow_subset,
        )
    except EvaluationError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(
        {
            "cards_dir": str(cards_dir),
            "n_cards": index["n_cards"],
            "n_missing": index["n_missing"],
            "total_bytes": sum(row["card_bytes"] for row in index["cards"]),
        },
        indent=2,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
