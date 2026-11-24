"""Ordered Stage 2 prediction axes and bounded projection of source-native counts."""
from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from hashlib import file_digest
from pathlib import Path

SCHEMA = "stage2-hvg-deg-axes-v1"
SIZES = (2000, 3000, 4000, 5000)
COMMON_SIZE = 18164
COMMON_SHA256 = "19ce44acfd615a601d5bc6617074941e13d51cd6e527c8d519174e5085385f42"
AXIS_RELEASES = {
    "hvg-deg": ("02.1_hvg-deg", SCHEMA, SIZES),
    "denoise-deg": ("02.2_denoise-deg", "stage2-denoise-deg-axes-v2", (1000, 2000, 3000, 4000, 5000)),
}


def sha256(path):
    with Path(path).open("rb") as handle:
        return file_digest(handle, "sha256").hexdigest()


@dataclass(frozen=True)
class PredictionAxis:
    gene_indices: tuple[int, ...]
    symbols: tuple[str, ...]
    sha256: str

    def project(self, corpus, source_name):
        """Bind this axis to a corpus source with complete common-gene coverage."""
        if corpus.manifest["common_axis"] != {"size": COMMON_SIZE, "sha256": COMMON_SHA256}:
            raise ValueError("corpus uses a different common gene axis")
        source = corpus.source(source_name)
        if not source.shared_axis_complete:
            raise ValueError(f"{source_name} does not cover the fixed common gene axis")
        return Projection(source, source.native_columns(self.gene_indices))


@dataclass(frozen=True)
class Projection:
    """A bounded, unnormalized view of one native count matrix in axis order."""

    source: object
    columns: tuple[int, ...]

    def counts(self, start, stop):
        return self.source.counts(start, stop)[:, self.columns]

    def batches(self, batch_size=256):
        for start, block in self.source.batches(batch_size):
            yield start, block[:, self.columns]


class Axes:
    def __init__(self, root, *, family="hvg-deg"):
        if family not in AXIS_RELEASES:
            raise ValueError(f"unknown prediction-axis family: {family}")
        step, schema, self.sizes = AXIS_RELEASES[family]
        self.family = family
        default = Path(__file__).resolve().parents[3] / step / "output/release"
        self.root = default if root is None else Path(root)
        path = self.root / "manifest.json"
        self.manifest = json.loads(path.read_text())
        self.manifest_sha256 = sha256(path)
        if self.manifest.get("schema") != schema:
            raise ValueError(f"unsupported {family} axes release")
        if self.manifest["common_axis"] != {"size": COMMON_SIZE, "sha256": COMMON_SHA256}:
            raise ValueError(f"{family} release must use the fixed Stage 2 common axis")
        if set(self.manifest["axes"]) != {str(size) for size in self.sizes}:
            raise ValueError(f"{family} release does not contain the accepted axis sizes")

    def axis(self, size):
        if size not in self.sizes:
            raise ValueError(f"unknown {self.family} prediction axis size: {size}")
        path = self.root / f"gene_axis_{size}.csv"
        digest = sha256(path)
        if digest != self.manifest["axes"][str(size)]:
            raise ValueError(f"prediction axis digest differs: {path}")
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames != ["gene_index", "symbol"]:
                raise ValueError(f"unexpected prediction axis columns: {path}")
            rows = list(reader)
        indices = tuple(int(row["gene_index"]) for row in rows)
        symbols = tuple(row["symbol"] for row in rows)
        if (len(indices) != size or indices != tuple(sorted(set(indices)))
                or any(not 0 <= gene < COMMON_SIZE for gene in indices)
                or not all(symbols)):
            raise ValueError(f"invalid ordered prediction axis: {path}")
        return PredictionAxis(indices, symbols, digest)


def open_axes(root=None, *, family="hvg-deg"):
    """Open a published axis family, defaulting to the accepted HVG-DEG release."""
    return Axes(root, family=family)
