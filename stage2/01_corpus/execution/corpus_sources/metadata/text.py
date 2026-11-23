"""Small deterministic text helpers shared by the loaders and the corpus."""

from __future__ import annotations

import re

# Sentence-like names in the GEO deposits use underscores and stray spaces.
_WS = re.compile(r"\s+")


def normalize_name(value: object) -> str:
    """Normalize a source-visible compound name.

    Underscores become spaces, leading/trailing whitespace is stripped and
    repeated whitespace collapses.  This is the *only* name canonicalization:
    it never changes chemical identity (e.g. ``Acetylsalicylic_acid`` ->
    ``Acetylsalicylic acid``; ``Acetaminophen `` -> ``Acetaminophen``).
    """
    if value is None:
        return ""
    text = str(value).replace("_", " ").strip()
    return _WS.sub(" ", text)
