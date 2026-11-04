#!/usr/bin/env python3
"""Shared helpers for editing the 开题 deck with python-pptx.

Import in a build/edit script (run from the repo root):

    import sys; sys.path.insert(0, "presentation/figures_src")
    from deck_edit import (by_name, remove_shape, replace_picture,
                           add_picture_at, fit_box)

All bboxes are (left, top, width, height) EMU tuples. `fit_box` preserves the
image aspect ratio and centres it inside the bbox, matching the convention used
by apply_deck_edits.py.
"""

from pathlib import Path

from PIL import Image
from pptx import Presentation  # noqa: F401  (re-exported for convenience)

_IN = 914400
ROOT = Path("/home/chang/arc")
DECK = ROOT / "presentation/开题3.pptx"


def open_deck(path=DECK):
    return Presentation(str(path))


def save_deck(prs, path=DECK):
    prs.save(str(path))
    return path


def by_name(slide, name):
    for shp in slide.shapes:
        if shp.name == name:
            return shp
    raise KeyError(name)


def remove_shape(shape):
    shape._element.getparent().remove(shape._element)


def fit_box(bbox, iw, ih):
    l, t, w, h = bbox
    s = min(w / iw, h / ih)
    nw, nh = int(iw * s), int(ih * s)
    return l + (w - nw) // 2, t + (h - nh) // 2, nw, nh


def replace_picture(slide, name, image_path, bbox=None):
    """Replace the picture named `name` with image_path (keeps bbox unless given)."""
    shape = by_name(slide, name)
    with Image.open(image_path) as im:
        iw, ih = im.size
    if bbox is None:
        bbox = (shape.left, shape.top, shape.width, shape.height)
    l, t, w, h = fit_box(bbox, iw, ih)
    remove_shape(shape)
    return slide.shapes.add_picture(str(image_path), l, t, w, h)


def add_picture_at(slide, image_path, bbox):
    """Add image_path, aspect-fit and centred inside bbox."""
    with Image.open(image_path) as im:
        iw, ih = im.size
    l, t, w, h = fit_box(bbox, iw, ih)
    return slide.shapes.add_picture(str(image_path), l, t, w, h)


def inch(x):
    return int(x * _IN)
