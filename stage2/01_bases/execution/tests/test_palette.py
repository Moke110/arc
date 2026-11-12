"""Stage palette contract: deterministic hex table, colour roles, published product.

The two Stage 2 comparison figures select every colour from this palette, so the
table must be reproducible from a recorded recipe and the published product must
be byte-identical to what the public implementation generates.
"""
import colorsys
import hashlib
import json
import math
import re
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases import palette

STEP_ROOT = Path(__file__).resolve().parents[2]
PRODUCT = STEP_ROOT / "output" / "palette"
HEX = re.compile(r"^#[0-9a-f]{6}$")

# Golden cells computed from the documented HLS recipe, independent of the module.
GOLDEN = {
    ("gray", 0): "#383838",
    ("blue", 2): "#1452a3",
    ("green", 2): "#14a33d",
    ("orange", 2): "#a35714",
    ("purple", 2): "#6a14a3",
    ("magenta", 2): "#a31478",
    ("blue", 6): "#729dd5",
    ("magenta", 8): "#deafd0",
}


def channels(value):
    """Split '#rrggbb' into its three integer channels."""
    return tuple(int(value[index:index + 2], 16) for index in (1, 3, 5))


def relative_luminance(value):
    """WCAG relative luminance of a hex colour."""
    parts = []
    for channel in channels(value):
        scaled = channel / 255.0
        parts.append(scaled / 12.92 if scaled <= 0.03928 else ((scaled + 0.055) / 1.055) ** 2.4)
    return 0.2126 * parts[0] + 0.7152 * parts[1] + 0.0722 * parts[2]


def recipe_table(recipe):
    """Regenerate the whole table from the recorded recipe alone."""
    assert recipe["rounding"]["operation"] == "floor(channel*255 + 0.5)"
    scale = recipe["rounding"]["scale"]
    families = [item["name"] for item in recipe["families"]]
    hues = {item["name"]: item["hue_degrees"] for item in recipe["families"]}
    rows = []
    for shade in range(recipe["shades"]):
        row = []
        for family in families:
            saturation = (0.0 if hues[family] is None
                          else round(recipe["saturation"]["start"]
                                     + recipe["saturation"]["step"] * shade, 6))
            lightness = round(recipe["lightness"]["start"]
                              + recipe["lightness"]["step"] * shade, 6)
            red, green, blue = colorsys.hls_to_rgb(
                (hues[family] or 0.0) / 360.0, lightness, saturation)
            row.append("#%02x%02x%02x" % tuple(math.floor(channel * scale + 0.5)
                                              for channel in (red, green, blue)))
        rows.append(tuple(row))
    return tuple(rows)


def hue_degrees(value):
    """HSL hue in degrees of a hex colour (None for an achromatic colour)."""
    red, green, blue = (channel / 255.0 for channel in channels(value))
    high, low = max(red, green, blue), min(red, green, blue)
    if high == low:
        return None
    span = high - low
    if high == red:
        hue = ((green - blue) / span) % 6.0
    elif high == green:
        hue = (blue - red) / span + 2.0
    else:
        hue = (red - green) / span + 4.0
    return hue * 60.0


class PaletteTableTests(unittest.TestCase):
    def test_table_is_six_families_by_nine_shades_of_distinct_hex(self):
        self.assertEqual(palette.FAMILIES,
                         ("gray", "blue", "green", "orange", "purple", "magenta"))
        self.assertEqual(palette.SHADE_COUNT, 9)
        table = palette.table()
        self.assertEqual(len(table), palette.SHADE_COUNT)
        for row in table:
            self.assertEqual(len(row), len(palette.FAMILIES))
            for cell in row:
                self.assertRegex(cell, HEX)
        flat = [cell for row in table for cell in row]
        self.assertEqual(len(set(flat)), palette.SHADE_COUNT * len(palette.FAMILIES))

    def test_each_family_gets_lighter_from_shade_zero_to_eight(self):
        table = palette.table()
        for column, family in enumerate(palette.FAMILIES):
            values = [row[column] for row in table]
            luminances = [relative_luminance(value) for value in values]
            self.assertEqual(luminances, sorted(luminances), family)
            self.assertEqual(len(set(luminances)), palette.SHADE_COUNT, family)

    def test_gray_is_achromatic_and_every_family_keeps_its_declared_hue(self):
        for shade in range(palette.SHADE_COUNT):
            red, green, blue = channels(palette.color("gray", shade))
            self.assertEqual((red, green, blue), (red, red, red))
        bounds = {"blue": (190, 240), "green": (100, 160), "orange": (10, 50),
                  "purple": (250, 300), "magenta": (300, 340)}
        for family, (low, high) in bounds.items():
            for shade in range(palette.SHADE_COUNT):
                hue = hue_degrees(palette.color(family, shade))
                self.assertIsNotNone(hue, (family, shade))
                self.assertTrue(low <= hue <= high, (family, shade, hue))

    def test_recorded_recipe_reproduces_every_cell(self):
        recipe = palette.RECIPE
        self.assertEqual(recipe["schema"], "stage2-tenv-stage-palette-recipe-v1")
        self.assertEqual(tuple(item["name"] for item in recipe["families"]), palette.FAMILIES)
        self.assertEqual(recipe["shades"], palette.SHADE_COUNT)
        self.assertTrue(recipe["algorithm"])
        for (family, shade), expected in GOLDEN.items():
            self.assertEqual(palette.color(family, shade), expected, (family, shade))

    def test_recorded_recipe_regenerates_all_fifty_four_cells(self):
        self.assertEqual(recipe_table(palette.RECIPE), palette.table())
        committed = tuple(tuple(line.split(",")[1:])
                          for line in (PRODUCT / "palette.csv").read_text().splitlines()[1:])
        self.assertEqual(recipe_table(palette.RECIPE), committed)

    def test_unknown_family_shade_and_run_index_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "family"):
            palette.color("red", 0)
        with self.assertRaisesRegex(ValueError, "shade"):
            palette.color("blue", palette.SHADE_COUNT)
        with self.assertRaisesRegex(ValueError, "shade"):
            palette.color("blue", -1)
        with self.assertRaisesRegex(ValueError, "run"):
            palette.run_shade(palette.SHADE_COUNT)


class PaletteRoleTests(unittest.TestCase):
    def test_entry_line_edge_and_dot_share_the_hue_and_lighter_fill_fills(self):
        for hue in palette.FAMILIES:
            colors = palette.role_colors(entry_hue=hue, run_index=0)
            self.assertEqual(colors["line"], palette.color(hue, palette.LINE_SHADE))
            for role in ("bar_edge", "dot"):
                self.assertEqual(colors[role], colors["line"])
            for role in ("band", "bar_fill"):
                self.assertEqual(colors[role], palette.color(hue, palette.FILL_SHADE))
            self.assertGreater(palette.FILL_SHADE, palette.LINE_SHADE)
            self.assertNotEqual(colors["line"], colors["band"])

    def test_train_and_validation_share_one_shade_per_run_in_blue_and_green(self):
        for run_index in range(4):
            colors = palette.role_colors(entry_hue="blue", run_index=run_index)
            shade = palette.run_shade(run_index)
            self.assertEqual(colors["train"], palette.color("blue", shade))
            self.assertEqual(colors["validation"], palette.color("green", shade))
        self.assertEqual(len({palette.run_shade(index) for index in range(4)}), 4)

    def test_every_role_the_two_figures_need_resolves_from_the_table(self):
        colors = palette.role_colors(entry_hue="blue", run_index=0)
        self.assertEqual(set(colors), set(palette.ROLES))
        self.assertEqual(len(palette.ROLES), 12)
        self.assertEqual(colors["dot_edge"], palette.color("gray", 8))
        published = {cell for row in palette.table() for cell in row}
        for role, value in colors.items():
            self.assertIn(value, published, role)

    def test_role_colors_require_the_identity_their_roles_depend_on(self):
        with self.assertRaisesRegex(ValueError, "entry hue"):
            palette.role_colors(run_index=0)
        with self.assertRaisesRegex(ValueError, "run"):
            palette.role_colors(entry_hue="blue")


class PaletteProductTests(unittest.TestCase):
    def test_published_product_is_the_generated_table_plus_recipe(self):
        manifest = json.loads((PRODUCT / "manifest.json").read_text())
        self.assertEqual(manifest["schema"], "stage2-tenv-stage-palette-v1")
        self.assertEqual(tuple(manifest["families"]), palette.FAMILIES)
        self.assertEqual(manifest["shade_count"], palette.SHADE_COUNT)
        self.assertEqual(manifest["recipe"], palette.RECIPE)
        table = (PRODUCT / "palette.csv").read_text().splitlines()
        self.assertEqual(table[0], "shade," + ",".join(palette.FAMILIES))
        self.assertEqual(len(table), palette.SHADE_COUNT + 1)
        for index, line in enumerate(table[1:]):
            cells = line.split(",")
            self.assertEqual(int(cells[0]), index)
            self.assertEqual(tuple(cells[1:]), palette.table()[index])
        payload = (PRODUCT / "palette.csv").read_bytes()
        self.assertEqual(manifest["table"]["sha256"], hashlib.sha256(payload).hexdigest())
        self.assertEqual(manifest["table"]["bytes"], len(payload))

    def test_publishing_reproduces_the_committed_product_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            palette.publish(Path(tmp))
            for name in ("palette.csv", "manifest.json"):
                self.assertEqual((Path(tmp) / name).read_bytes(),
                                 (PRODUCT / name).read_bytes(), name)

    def test_product_roles_name_the_entry_and_split_colours(self):
        manifest = json.loads((PRODUCT / "manifest.json").read_text())
        roles = manifest["roles"]
        self.assertEqual(roles["line"], {"family": "entry", "shade": palette.LINE_SHADE})
        self.assertEqual(roles["band"], {"family": "entry", "shade": palette.FILL_SHADE})
        self.assertEqual(roles["train"], {"family": "blue", "shade": "run"})
        self.assertEqual(roles["validation"], {"family": "green", "shade": "run"})
        self.assertEqual(set(roles), set(palette.ROLES))


if __name__ == "__main__":
    unittest.main()
