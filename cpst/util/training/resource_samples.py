"""Read and verify timed resource samples."""

import json

import numpy as np


def read_samples(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def gaps_summary(samples):
    gaps = [later["monotonic_seconds"] - earlier["monotonic_seconds"]
            for earlier, later in zip(samples, samples[1:])
            if earlier.get("monitor_session") == later.get("monitor_session")]
    if not gaps or min(gaps) <= 0 or not np.isfinite(gaps).all():
        raise ValueError("resource audit must record monotonic samples")
    return {"samples": len(samples), "sample_interval_seconds": 10,
            "maximum_gap_seconds": round(max(gaps), 3), "gaps_over_12_seconds": sum(gap > 12 for gap in gaps)}
