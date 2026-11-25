"""Stage 2 hue-depth coordinates and transitional family-shade colors.

Use :func:`color_pair` for the active palette (ADR 0093).

The original Stage 2 palette has six hue families crossed with nine shades.

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
import math
from numbers import Integral
from pathlib import Path
from typing import NamedTuple


HUES = tuple(range(20, 341, 40))
DEPTHS = tuple(range(5))
DEPTH_LEVELS = ((.24, .80), (.33, .76), (.42, .72), (.52, .66), (.62, .60))
MAIN_OPACITY = .90
SECONDARY_OPACITY = .40


class ColorPair(NamedTuple):
    """RGBA colors for main marks and secondary fills at one coordinate."""

    main: tuple[float, float, float, float]
    secondary: tuple[float, float, float, float]


def color_pair(hue, depth=2):
    """Resolve an approved hue in degrees and zero-based depth into RGBA colors."""
    if not isinstance(hue, Integral) or hue not in HUES:
        raise ValueError(f"palette hue must be one of {HUES}, got {hue!r}")
    if isinstance(depth, bool) or not isinstance(depth, Integral) or depth not in DEPTHS:
        raise ValueError(f"palette depth must be one of {DEPTHS}, got {depth!r}")
    lightness, saturation = DEPTH_LEVELS[depth]
    rgb = colorsys.hls_to_rgb(hue / 360, lightness, saturation)
    return ColorPair((*rgb, MAIN_OPACITY), (*rgb, SECONDARY_OPACITY))


def sequential_colors(hue):
    """Return the five main colors of one hue from light to dark."""
    return tuple(color_pair(hue, depth).main for depth in reversed(DEPTHS))


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
# The thin reference line and the threshold marker sit between the line and the
# fill, and darker than the line, so they stay distinguishable inside the entry's
# family without leaving the stage palette.
REFERENCE_SHADE = 4
MARKER_SHADE = 0
# The three-class confusion heat map is a sequential gray ramp, ordered from a
# low count (shade 8) to a high count (shade 0), so a categorical entry hue never
# competes with it.
HEATMAP_FAMILY = "gray"
HEATMAP_SHADES = tuple(reversed(range(SHADE_COUNT)))
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
    "dot_edge": {"family": "gray", "shade": 8},
    "band": {"family": "entry", "shade": FILL_SHADE},
    "bar_fill": {"family": "entry", "shade": FILL_SHADE},
    "reference": {"family": "entry", "shade": REFERENCE_SHADE},
    "marker": {"family": "entry", "shade": MARKER_SHADE},
    "train": {"family": TRAIN_FAMILY, "shade": "run"},
    "validation": {"family": VALIDATION_FAMILY, "shade": "run"},
}
ROLES = tuple(_ROLE_SPECS)

RECIPE = {
    "schema": "stage2-tenv-stage-palette-recipe-v1",
    "algorithm": (
        "Per shade s in 0..8: L = 0.22 + 0.07*s, S = 0.90 - 0.06*s (S = 0 for "
        "gray). Convert (H/360, L, S) with colorsys.hls_to_rgb, then each "
        "channel is floor(channel*255 + 0.5) and formatted as lowercase #rrggbb."
    ),
    "rounding": {"operation": "floor(channel*255 + 0.5)", "scale": 255},
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
        "reference": REFERENCE_SHADE,
        "marker": MARKER_SHADE,
    },
    "heatmap": {"family": HEATMAP_FAMILY,
                "shades": list(HEATMAP_SHADES),
                "order": "low count -> high count"},
}


def heatmap_colors():
    """Return the sequential heat-map ramp, low count first, from the table."""
    return tuple(color(HEATMAP_FAMILY, shade) for shade in HEATMAP_SHADES)


def heatmap_color(fraction):
    """Return the ramp cell at a normalized ``fraction`` in ``[0, 1]``."""
    if isinstance(fraction, bool) or not isinstance(fraction, (int, float)):
        raise ValueError("a heat-map fraction must be a number in [0, 1]")
    if not 0 <= float(fraction) <= 1:
        raise ValueError("a heat-map fraction must be a number in [0, 1]")
    ramp = heatmap_colors()
    index = min(int(float(fraction) * len(ramp)), len(ramp) - 1)
    return ramp[index]


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
    return "#%02x%02x%02x" % tuple(math.floor(channel * 255 + 0.5)
                                  for channel in (red, green, blue))


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
        "module": {"path": "output/api/stage2_util/palette.py",
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
