"""The Stage 2 palette: six hue families crossed with nine shades.

Every Stage 2 comparison figure selects its colours from this one table so that a
hue means the same entry or split across steps (ADR 0063).  The table is a pure
function of the recipe in :data:`RECIPE`: no file is read to compute a colour, so
a renderer that has only recorded plot values and this public implementation can
reproduce a figure without the step's ``output/`` tree.

Layout.  Shade index 0 is the darkest, most saturated cell of a family and index
8 is the lightest; each row is one shade across the six families, in
:data:`FAMILIES` order.  A figure draws an entry with :data:`LINE_SHADE` for its
line, bar edge and dots, and fills its band and bar with the lighter
:data:`FILL_SHADE` of the same family.  A learning-curve split uses the blue
family for training and the green family for validation at one shade per
``(fold, seed)`` run (:func:`run_shade`).
"""
import colorsys
import hashlib
import json
from pathlib import Path

SCHEMA = "stage2-tenv-stage-palette-v1"
FAMILIES = ("gray", "blue", "green", "orange", "purple", "magenta")
SHADE_COUNT = 9
# Base hue in degrees; ``None`` is the achromatic gray family.
FAMILY_HUES = {"gray": None, "blue": 214.0, "green": 137.0,
               "orange": 28.0, "purple": 276.0, "magenta": 318.0}
LIGHTNESS_START, LIGHTNESS_STEP = 0.22, 0.07
SATURATION_START, SATURATION_STEP = 0.90, -0.06

# A figure entry draws in one shade and fills with a lighter one of its family.
LINE_SHADE = 2
FILL_SHADE = 6
# Entry order used when a caller does not name hues; the deterministic reference
# (gray) comes last so trained entrants take the chromatic families first.
ENTRY_HUES = ("blue", "orange", "green", "purple", "magenta", "gray")
TRAIN_FAMILY = "blue"
VALIDATION_FAMILY = "green"
# One shade per (fold, seed) run; four runs use 2, 4, 6, 8.
RUN_SHADE_ORDER = (2, 4, 6, 8, 1, 3, 5, 7, 0)

# Role -> colour rule.  ``family: 'entry'`` means the entry's own hue and
# ``shade: 'run'`` means that run's shade in the named family.
_ROLE_SPECS = {
    "text": {"family": "gray", "shade": 0},
    "grid": {"family": "gray", "shade": 8},
    "undefined_marker": {"family": "gray", "shade": 4},
    "undefined_label": {"family": "gray", "shade": 4},
    "line": {"family": "entry", "shade": LINE_SHADE},
    "bar_edge": {"family": "entry", "shade": LINE_SHADE},
    "dot": {"family": "entry", "shade": LINE_SHADE},
    "band": {"family": "entry", "shade": FILL_SHADE},
    "bar_fill": {"family": "entry", "shade": FILL_SHADE},
    "train": {"family": TRAIN_FAMILY, "shade": "run"},
    "validation": {"family": VALIDATION_FAMILY, "shade": "run"},
}
ROLES = tuple(_ROLE_SPECS)

RECIPE = {
    "schema": "stage2-tenv-stage-palette-recipe-v1",
    "algorithm": (
        "Per shade s in 0..8: L = 0.22 + 0.07*s, S = 0.90 - 0.06*s (S = 0 for "
        "gray). Convert (H/360, L, S) with colorsys.hls_to_rgb, then each "
        "channel is int(round(channel*255 + 0.5)) and formatted as lowercase "
        "#rrggbb."
    ),
    "families": [{"name": name, "hue_degrees": FAMILY_HUES[name]} for name in FAMILIES],
    "shades": SHADE_COUNT,
    "lightness": {"start": LIGHTNESS_START, "step": LIGHTNESS_STEP},
    "saturation": {"start": SATURATION_START, "step": SATURATION_STEP, "gray": 0.0},
    "shade_roles": {
        "line": LINE_SHADE,
        "fill": FILL_SHADE,
        "train_family": TRAIN_FAMILY,
        "validation_family": VALIDATION_FAMILY,
        "run_shade_order": list(RUN_SHADE_ORDER),
        "entry_hue_order": list(ENTRY_HUES),
    },
}


def _check_family(family):
    if family not in FAMILY_HUES:
        raise ValueError(f"unknown palette family: {family!r}")
    return family


def _check_shade(shade):
    if not isinstance(shade, int) or isinstance(shade, bool) or not 0 <= shade < SHADE_COUNT:
        raise ValueError(f"palette shade must be an integer in 0..{SHADE_COUNT - 1}, got {shade!r}")
    return shade


def _hex(hue, saturation, lightness):
    red, green, blue = colorsys.hls_to_rgb((hue or 0.0) / 360.0, lightness, saturation)
    return "#%02x%02x%02x" % tuple(int(channel * 255 + 0.5) for channel in (red, green, blue))


def color(family, shade):
    """Return the ``#rrggbb`` hex of one table cell."""
    family = _check_family(family)
    shade = _check_shade(shade)
    saturation = 0.0 if family == "gray" else round(SATURATION_START + SATURATION_STEP * shade, 6)
    lightness = round(LIGHTNESS_START + LIGHTNESS_STEP * shade, 6)
    return _hex(FAMILY_HUES[family], saturation, lightness)


def table():
    """Return the palette as ``SHADE_COUNT`` rows of one hex per family."""
    return tuple(tuple(color(family, shade) for family in FAMILIES)
                 for shade in range(SHADE_COUNT))


def entry_hue(index):
    """Return the family assigned to the ``index``-th figure entry."""
    if not isinstance(index, int) or isinstance(index, bool) or not 0 <= index < len(ENTRY_HUES):
        raise ValueError(f"entry index must be an integer in 0..{len(ENTRY_HUES) - 1}, got {index!r}")
    return ENTRY_HUES[index]


def run_shade(index):
    """Return the shade that identifies the ``index``-th ``(fold, seed)`` run."""
    if not isinstance(index, int) or isinstance(index, bool) or not 0 <= index < len(RUN_SHADE_ORDER):
        raise ValueError(f"run index must be an integer in 0..{len(RUN_SHADE_ORDER) - 1}, got {index!r}")
    return RUN_SHADE_ORDER[index]


def role_colors(*, entry_hue=None, run_index=None):
    """Resolve every colour role the two figures need, for one entry and run."""
    if entry_hue is not None:
        _check_family(entry_hue)
    elif any(spec["family"] == "entry" for spec in _ROLE_SPECS.values()):
        raise ValueError("entry hue is required for the entry colour roles")
    if run_index is None:
        raise ValueError("run index is required for the train/validation colour roles")
    resolved = {}
    for role, spec in _ROLE_SPECS.items():
        family = entry_hue if spec["family"] == "entry" else spec["family"]
        shade = run_shade(run_index) if spec["shade"] == "run" else spec["shade"]
        resolved[role] = color(family, shade)
    return resolved


def rows():
    """Return the table as long-form ``shade, family, hex`` records."""
    return [{"shade": shade, "family": family, "hex": color(family, shade)}
            for shade in range(SHADE_COUNT) for family in FAMILIES]


def manifest():
    """Return the product manifest; digests are filled by :func:`publish`."""
    return {
        "schema": SCHEMA,
        "families": list(FAMILIES),
        "shade_count": SHADE_COUNT,
        "recipe": RECIPE,
        "roles": {role: dict(spec) for role, spec in _ROLE_SPECS.items()},
        "module": {"path": "output/api/bases/palette.py",
                   "sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
    }


def publish(directory):
    """Materialise the table and manifest under ``directory``; return the manifest."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    lines = ["shade," + ",".join(FAMILIES)]
    for shade, row in enumerate(table()):
        lines.append(",".join((str(shade),) + row))
    payload = ("\n".join(lines) + "\n").encode()
    (directory / "palette.csv").write_bytes(payload)
    value = manifest()
    value["table"] = {"file": "palette.csv", "rows": SHADE_COUNT,
                      "columns": len(FAMILIES), "bytes": len(payload),
                      "sha256": hashlib.sha256(payload).hexdigest()}
    (directory / "manifest.json").write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    return value
