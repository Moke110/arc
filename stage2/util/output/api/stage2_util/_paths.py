"""Resolve the repository hosting this Stage 2 utility API."""
from pathlib import Path


def repository():
    return next(parent for parent in Path(__file__).resolve().parents
                if (parent / ".git").exists())
