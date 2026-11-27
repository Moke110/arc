"""Read the accepted Stage 2 seeds in their stable order."""

import csv
import hashlib
from io import StringIO
from numbers import Integral
from pathlib import Path

SEEDS_SHA256 = "21e14ff3aff48ee3ddb624c344f728c9a08505e37b714150da9089c4a1760041"


def read_seeds(count=None):
    """Return the entire accepted seed tuple or the requested prefix."""
    path = Path(__file__).resolve().parents[2] / "tables" / "seeds.csv"
    payload = path.read_bytes()
    if hashlib.sha256(payload).hexdigest() != SEEDS_SHA256:
        raise ValueError("Stage 2 seed table digest mismatch")
    seeds = tuple(int(row["seed"]) for row in csv.DictReader(StringIO(payload.decode())))
    if count is None:
        return seeds
    if isinstance(count, bool) or not isinstance(count, Integral) or not 1 <= count <= len(seeds):
        raise ValueError(f"expected a seed count between 1 and {len(seeds)}")
    return seeds[:count]
