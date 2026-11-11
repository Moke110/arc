"""Publish ``figure/data/diagnostic_tenv_pr.csv`` and prove the CSV-only PNG replay (#156).

The accepted diagnostic evidence under ``log/diagnostic/`` is the scientific
source: ``diagnostic_pr.csv`` holds the plotted curve coordinates and
``diagnostic_pr.json`` holds the displayed equal-group macro AP values and the
fixed-group definition. That evidence is read-only here. This publisher writes a
single complete plotted-values table, renders the PNG from the table alone with
``figure/script/render_diagnostic_tenv_pr.py``, and compares the rendered bytes
with the accepted PNG (sha256 35d14266…).

    python3 stage2/01_bases/execution/publish_diagnostic_figure.py \
        --step-root /home/chang/arc/stage2/01_bases \
        --figure /home/chang/arc/.scratch/worktrees/stage2-bases-156/stage2/01_bases/figure \
        --scratch /home/chang/arc/.scratch/stage2-bases/runtime/v3/figure-render \
        --record /home/chang/arc/.scratch/worktrees/stage2-bases-156/stage2/01_bases/log/156-figure-reproduction.json
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

SCHEMA = "bases-156-diagnostic-figure-v1"
ACCEPTED_PNG = "35d14266283b296c1585af3d58f0bacb15a3d7f854cde4972c433375bbdd911f"
SEQUENCE = ("AUPRC", "AUPRC_positive", "AUPRC_negative")
PANEL_TITLES = {"AUPRC": "Combined",
                "AUPRC_positive": "Upregulated",
                "AUPRC_negative": "Downregulated"}
FIELDS = ("record", "panel", "panel_title", "series", "title", "recall", "precision",
          "contributing_groups", "eligible_groups", "macro_ap", "group_definition",
          "uncertainty", "reference")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def plotted_rows(evidence):
    """Read the accepted curve coordinates without changing their values."""
    with open(evidence / "diagnostic_pr.csv", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or {row["metric"] for row in rows} != set(SEQUENCE):
        raise ValueError("accepted plotted CSV does not cover the three metrics")
    return rows


def build_table(evidence):
    rows = plotted_rows(evidence)
    report = json.loads((evidence / "diagnostic_pr.json").read_text())
    summaries = {(item["label"], metric): item[metric]
                 for item in report["summaries"] for metric in SEQUENCE}
    labels = []
    for row in rows:
        if row["label"] not in labels:
            labels.append(row["label"])
    if len(labels) < 2 or len(summaries) != len(labels) * len(SEQUENCE):
        raise ValueError("accepted panel summaries and curve labels disagree")
    table = [{"record": "figure", "title": report["title"],
              "group_definition": report["curve_semantics"],
              "uncertainty": "none - one deterministic replay; no interval is plotted",
              "reference": "none - no external reference curve is plotted"}]
    for metric in SEQUENCE:
        table.append({"record": "panel", "panel": metric,
                      "panel_title": PANEL_TITLES[metric]})
    for label in labels:
        for metric in SEQUENCE:
            part = [row for row in rows if row["label"] == label and row["metric"] == metric]
            if not part:
                raise ValueError(f"accepted CSV has no curve for {label}/{metric}")
            table.append({"record": "summary", "panel": metric, "series": label,
                          "macro_ap": repr(float(summaries[(label, metric)])),
                          "eligible_groups": part[0]["eligible_groups"],
                          "group_definition": report["curve_semantics"]})
            for row in part:
                table.append({"record": "curve", "panel": metric, "series": label,
                              "recall": row["recall"], "precision": row["precision"],
                              "contributing_groups": row["contributing_groups"],
                              "eligible_groups": row["eligible_groups"]})
    return table


def serialize(table):
    from io import StringIO
    stream = StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
    writer.writeheader()
    for row in table:
        writer.writerow({field: row.get(field, "") for field in FIELDS})
    return stream.getvalue().encode()


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--step-root", type=Path,
                        default=Path(__file__).resolve().parents[1])
    parser.add_argument("--figure", type=Path, help="figure root (defaults to step root/figure)")
    parser.add_argument("--scratch", type=Path, required=True,
                        help="candidate rendering directory outside the published figure tree")
    parser.add_argument("--record", type=Path)
    args = parser.parse_args()
    step = args.step_root.resolve()
    figure = (args.figure or (step / "figure")).resolve()
    evidence = step / "log/diagnostic"
    scratch = args.scratch.resolve()
    scratch.mkdir(parents=True, exist_ok=True)

    table = build_table(evidence)
    payload = serialize(table)
    values = figure / "data/diagnostic_tenv_pr.csv"
    values.parent.mkdir(parents=True, exist_ok=True)
    if values.exists():
        current = values.read_bytes()
        accepted = (evidence / "diagnostic_pr.csv").read_bytes()
        if current not in (payload, accepted):
            raise ValueError("published plotted-values CSV is neither the accepted source "
                             "nor this ticket's complete table; refusing to overwrite evidence")
    values.write_bytes(payload)

    candidate = scratch / "diagnostic_tenv_pr.png"
    if candidate.exists():
        candidate.unlink()
    renderer = figure / "script/render_diagnostic_tenv_pr.py"
    subprocess.run([sys.executable, str(renderer), "--csv", str(values),
                    "--png", str(candidate)], check=True,
                   env={"MPLCONFIGDIR": str(scratch / "matplotlib"), "PATH": "/usr/bin:/bin",
                        "HOME": str(scratch)})
    rendered = digest(candidate)
    accepted = digest(evidence / "diagnostic_pr.png")
    published = figure / "diagnostic_tenv_pr.png"
    if rendered != accepted or accepted != ACCEPTED_PNG:
        raise ValueError(f"rendered PNG {rendered} differs from accepted {accepted}")
    if published.exists() and digest(published) != rendered:
        raise ValueError("published figure PNG is not the accepted panel")
    published.write_bytes(candidate.read_bytes())

    rows = [row for row in table if row["record"] == "curve"]
    series = sorted({row["series"] for row in table if row["record"] == "summary"})
    record = {
        "schema": SCHEMA, "ticket": 156, "parent": 154,
        "result": "PASS",
        "figure": str(published), "values": str(values),
        "renderer": str(renderer),
        "renderer_reads_only": "figure/data/diagnostic_tenv_pr.csv",
        "command": [sys.executable, str(renderer), "--csv", str(values), "--png", str(candidate)],
        "digests": {"png": rendered, "csv": hashlib.sha256(payload).hexdigest(),
                    "accepted_png": accepted,
                    "evidence_csv": digest(evidence / "diagnostic_pr.csv"),
                    "evidence_json": digest(evidence / "diagnostic_pr.json")},
        "counts": {"curve_rows": len(rows), "summary_rows": len(SEQUENCE) * len(series),
                   "panels": len(SEQUENCE), "series": len(series)},
        "values_include": ["recall/precision curve coordinates",
                           "fixed contributing_groups and eligible_groups per point",
                           "displayed equal-group macro AP per panel and series",
                           "fixed-group pooled curve definition",
                           "uncertainty and reference-line statement"],
        "byte_identical_png": rendered == accepted,
        "environment": {"python": sys.version.split()[0],
                        "matplotlib": __import__("matplotlib").__version__,
                        "pandas": __import__("pandas").__version__},
    }
    if args.record:
        args.record.parent.mkdir(parents=True, exist_ok=True)
        args.record.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps(record, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
