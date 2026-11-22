"""Rendering primitives for figures whose plotted values are recorded in CSV."""

import csv
from pathlib import Path


def read_rows(path):
    """Read the CSV's rows without loading any study execution assets."""
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def save_figure(figure, path, **kwargs):
    """Write a figure to a destination, creating its parent directory if needed."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, **kwargs)
