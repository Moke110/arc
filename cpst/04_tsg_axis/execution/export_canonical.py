"""Export the accepted matched-control result tables with current axis IDs."""
import argparse
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

STEP = Path(__file__).resolve().parents[1]
STAGE = STEP.parent
sys.path[:0] = [str(STAGE / "util/output/api"), str(STAGE / "03_data_train/output/api")]

from cpst_util.derivation.release import check_outputs, read_json, sha256, write_json
from cpst_util.publication import publication_lock
from training_data import TENV_AXIS_STEPS

RECEIPT = STEP / "log/319-formal.json"
AXIS_TABLES = (
    "condition-disease-f1.csv", "condition-metrics.csv", "fold-disease-f1.csv",
    "fold-metrics.csv", "figure/data/delta-metrics.csv", "figure/data/disease-f1.csv",
)
AXIS_IDS = {axis_id.replace("-", "-deg-", 1): axis_id for axis_id in TENV_AXIS_STEPS}


def export_table(source, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    ids = {old.encode(): new.encode() for old, new in AXIS_IDS.items()}
    seen = set()
    with source.open("rb") as original, target.open("xb") as output:
        header = original.readline()
        if not header.startswith(b"axis,"):
            raise ValueError(f"sealed table lacks its leading axis column: {source}")
        output.write(header)
        for line in original:
            axis, separator, remaining = line.partition(b",")
            if not separator or axis not in ids:
                raise ValueError(f"unexpected sealed axis ID in {source}: {axis!r}")
            output.write(ids[axis] + separator + remaining)
            seen.add(axis)
    if seen != ids.keys():
        raise ValueError(f"sealed table lacks one or more accepted axes: {source}")


def export_report(report, output):
    """Verify the sealed report, then publish an independent axis-name export."""
    report = Path(report).resolve()
    output = Path(output).resolve()
    if output.is_relative_to(report) or report.is_relative_to(output):
        raise ValueError("canonical export must be separate from the sealed report")
    receipt = read_json(RECEIPT)["report"]
    manifest_path = report / "manifest.json"
    if sha256(manifest_path) != receipt["report_manifest_sha256"]:
        raise ValueError("report is not the accepted matched-control release")
    sealed = check_outputs(report, read_json(manifest_path))
    if sealed != receipt["output_sha256"]:
        raise ValueError("report output inventory differs from the accepted receipt")
    with publication_lock(output) as destination:
        if destination.exists():
            raise ValueError("canonical export requires an unused output directory")
        with TemporaryDirectory(prefix="canonical-report-", dir=destination.parent) as scratch:
            candidate = Path(scratch) / "result"
            for name in AXIS_TABLES:
                export_table(report / name, candidate / name)
            record = {
                "axis_ids": AXIS_IDS,
                "sealed_report_manifest_sha256": receipt["report_manifest_sha256"],
                "sealed_table_sha256": {name: sealed[name] for name in AXIS_TABLES},
                "outputs": {name: sha256(candidate / name) for name in AXIS_TABLES},
            }
            write_json(candidate / "manifest.json", record)
            candidate.rename(destination)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, default=STEP / "output/report-matched-control")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(export_report(args.report, args.output), sort_keys=True))


if __name__ == "__main__":
    main()
