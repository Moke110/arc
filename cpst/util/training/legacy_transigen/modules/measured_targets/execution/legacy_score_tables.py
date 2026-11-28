"""Read and verify aggregates for downstream CPST comparison tables."""

import csv
import math
import statistics


def read_table(path, keys):
    rows = {}
    with path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            key = tuple(int(row[name]) if name in {"fold", "seed"} else row[name] for name in keys)
            if key in rows:
                raise ValueError(f"duplicate {path.name} row {key}")
            rows[key] = row
    return rows


def value(row, name="value"):
    raw = row[name]
    if not raw:
        return None
    parsed = float(raw)
    if not math.isfinite(parsed):
        raise ValueError(f"nonfinite defined value in {name}")
    return parsed


def same(found, expected):
    if found is None or expected is None:
        return found is expected
    return math.isclose(found, expected, rel_tol=1e-12, abs_tol=1e-12)


def check_statistics(row, values, label, field="value", sd="sample_sd", count="n_defined"):
    defined = [number for number in values if number is not None]
    mean = statistics.mean(defined) if defined else None
    sample_sd = statistics.stdev(defined) if len(defined) > 1 else None
    if (int(row[count]) != len(defined) or int(row["n_coordinates"]) != len(values)
            or not same(value(row, field), mean) or not same(value(row, sd), sample_sd)):
        raise ValueError(f"invalid aggregate {label}")
