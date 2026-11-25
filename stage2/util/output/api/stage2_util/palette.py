"""Stage 2 hue-depth coordinates and mark colors (ADR 0093)."""
import colorsys
from numbers import Integral
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
