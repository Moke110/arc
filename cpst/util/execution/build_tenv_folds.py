"""Freeze and authenticate the CPST tenv evaluation ledger from all eight axis releases."""

import argparse
import csv
from hashlib import file_digest, sha256
from io import StringIO
import json
from pathlib import Path

from cpst_util.seeds import read_seeds
from cpst_util.tenv_folds import (
    COLD_DOSE_DRUG, COLD_DRUG, COLD_SEEDS, CONDITION_COLUMNS, CONDITION_FILE,
    HELD_DOSES, HISTORICAL_COLUMNS, LEDGER_COLUMNS, SCHEMA, SPLIT_SEED, Condition,
    compute_partition_digests, open_tenv_folds,
)

ROOT = Path(__file__).resolve().parents[3]
TRACKED_STEP03 = ROOT / "cpst/03_data_train"
AXES = tuple(f"{family}-{width}k" for family in ("hvg", "denoise") for width in range(2, 6))
SOURCE_COLUMNS = (*CONDITION_COLUMNS, "control_row", "member_count", *HISTORICAL_COLUMNS, *(
    f"cold_dose_drug_{group}{suffix}" for group in range(3) for suffix in ("", "_reason")
))


def digest(payload):
    return sha256(payload).hexdigest()


def axis_root(root, axis):
    path, = root.glob(f"03.*_tenv_{axis}/output")
    return path


def read_axis(training, axis):
    source = axis_root(training, axis)
    tracked = (axis_root(TRACKED_STEP03, axis) / "manifest.json").read_bytes()
    manifest = (source / "manifest.json").read_bytes()
    if tracked != manifest:
        raise ValueError(f"tenv axis release manifest differs from the tracked input: {axis}")
    metadata = json.loads(manifest)
    outputs = metadata["outputs"]
    for name, expected in outputs.items():
        with (source / name).open("rb") as handle:
            if file_digest(handle, "sha256").hexdigest() != expected:
                raise ValueError(f"tenv axis payload differs: {axis}/{name}")
    payload = (source / "treatments.csv").read_bytes()
    rows = tuple(tuple(row[column] for column in SOURCE_COLUMNS)
                 for row in csv.DictReader(StringIO(payload.decode())))
    members = (source / "members.csv").read_bytes()
    return rows, {
        "manifest_sha256": digest(manifest),
        "treatments_sha256": digest(payload),
        "members_sha256": digest(members),
    }


def verify_historical_strips(rows):
    for group, dose in enumerate(HELD_DOSES):
        label = f"cold_dose_drug_{group}"
        for row in rows:
            old = row[f"cold_drug_{group}"]
            held_parent = old == "test"
            held_dose = float(row["dose_um"]) == dose
            expected = ("test", "") if held_parent and held_dose else (old, "")
            if held_parent != held_dose:
                expected = ("", "held_parent" if held_parent else "held_dose")
            if (row[label], row[f"{label}_reason"]) != expected:
                raise ValueError(f"historical tenv dose routing differs: {group}")


def build(training, output):
    provenance = {}
    reference, provenance[AXES[0]] = read_axis(training, AXES[0])
    for axis in AXES[1:]:
        observed, provenance[axis] = read_axis(training, axis)
        if observed != reference or provenance[axis]["members_sha256"] != provenance[AXES[0]]["members_sha256"]:
            raise ValueError(f"tenv condition order, membership or historical roles differ: {axis}")
    rows = tuple(dict(zip(SOURCE_COLUMNS, row, strict=True)) for row in reference)
    verify_historical_strips(rows)
    parents = {row["drug_parent"] for row in rows}
    ranked = sorted(parents, key=lambda parent: (digest(f"{SPLIT_SEED}:{parent}".encode()), parent))
    groups = {parent: index % 5 for index, parent in enumerate(ranked)}
    stream = StringIO()
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(LEDGER_COLUMNS)
    for row in rows:
        condition = Condition.from_row(row)
        writer.writerow((condition.plate, condition.drug_parent, str(condition.dose_um),
                         groups[condition.drug_parent], *(row[name] for name in HISTORICAL_COLUMNS)))
    payload = stream.getvalue().encode()
    manifest = {
        "schema": SCHEMA,
        "split_seed": SPLIT_SEED,
        "parent_group_rule": "Sort by SHA-256 of seed:parent, then by parent; assign rank modulo five.",
        "cold_drug_rotation": "Test group equals fold; validation group equals (fold + 1) modulo five.",
        "dose_crossing": {"held_doses_um": list(HELD_DOSES), "fold_order": "parent-group-major",
                          "development_roles": "Use original cold_drug roles where neither holdout matches."},
        "protocols": {COLD_DRUG: {"folds": 5, "seeds": list(COLD_SEEDS)},
                      COLD_DOSE_DRUG: {"folds": 9, "seeds": list(read_seeds())}},
        "partition_sha256": compute_partition_digests(payload),
        "axis_inputs": provenance,
        "outputs": {CONDITION_FILE: digest(payload)},
    }
    manifest_bytes = (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode()
    output.mkdir(parents=True)
    (output / CONDITION_FILE).write_bytes(payload)
    (output / "manifest.json").write_bytes(manifest_bytes)
    release = open_tenv_folds(output, manifest_sha256=digest(manifest_bytes))
    return release.release_sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training", type=Path, default=TRACKED_STEP03)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(build(args.training, args.output))


if __name__ == "__main__":
    main()
