"""Verify the accepted learning figure against the sealed epoch and run CSVs."""

import argparse
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

STEP = Path(__file__).resolve().parents[1]
SCRIPT_DIR = STEP / "figure/script"
sys.path.insert(0, str(SCRIPT_DIR))

from render_learning_curves import EPOCH_SOURCE, FIGURE_NAME, read_curves
from stage2_util.derivation.release import read_json, sha256
from stage2_util.plotting import read_rows


def audit(publication, scratch):
    formal = publication / "log/330-formal.json"
    runs_path = publication / "log/runs.csv"
    source = publication / "figure/data" / EPOCH_SOURCE
    figure = publication / "figure" / FIGURE_NAME
    receipt = read_json(formal)
    source_hash, runs_hash = sha256(source), sha256(runs_path)
    if (receipt["status"] != "PASS" or receipt["report_outputs_sha256"][EPOCH_SOURCE] != source_hash
            or receipt["report_outputs_sha256"]["runs.csv"] != runs_hash):
        raise ValueError("epoch and run CSVs must match the accepted formal report")

    curves = read_curves(source)
    runs = read_rows(runs_path)
    records = {(row["loss"], row["alpha"], int(row["fold"])): row for row in runs}
    if len(runs) != len(curves) or records.keys() != curves.keys():
        raise ValueError("learning curves must match the complete accepted run roster")
    evidence = []
    for key, points in curves.items():
        row = records[key]
        identity, _, _, validation_loss, validation_mse = points[int(row["best_epoch"]) - 1]
        expected = row["run_sha256"], int(row["seed"]), int(row["best_epoch"]), int(row["stopping_epoch"])
        if identity != expected or (validation_loss, validation_mse) != (
                float(row["validation_loss"]), float(row["validation_mse"])):
            raise ValueError(f"learning curve differs from its accepted checkpoint record: {key}")
        evidence.append({name: row[name] for name in ("run_sha256", "checkpoint_sha256", "trace_sha256")})

    renderer = SCRIPT_DIR / "render_learning_curves.py"
    figure_hash = sha256(figure)
    with TemporaryDirectory(prefix="step05-curves-", dir=scratch) as work:
        replay = Path(work)
        subprocess.run([sys.executable, str(renderer), "--data", str(source.parent), "--output", work], check=True)
        if sha256(replay / FIGURE_NAME) != figure_hash:
            raise ValueError("CSV-only learning-curve re-render differs")

    return {
        "status": "PASS",
        "ticket": 330,
        "accepted_formal_report_manifest_sha256": receipt["report_manifest_sha256"],
        "source_csv_sha256": source_hash,
        "accepted_runs_csv_sha256": runs_hash,
        "figure_sha256": figure_hash,
        "renderer_sha256": sha256(renderer),
        "run_count": len(curves),
        "panel_count": len({key[:2] for key in curves}),
        "epoch_rows": sum(map(len, curves.values())),
        "runs": evidence,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--publication", type=Path, default=STEP)
    parser.add_argument("--scratch", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(audit(args.publication, args.scratch), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
