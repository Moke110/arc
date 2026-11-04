#!/usr/bin/env python3
"""Apply the two 【补】 edits carried in presentation/开题3r.pptx.

Run from the repo root (python-pptx lives in the base conda env):

    /home/chang/miniconda3/bin/python presentation/figures_src/apply_deck_edits_3r.py

Writes presentation/开题4.pptx. Deterministic: always starts from the frozen
开题3r.pptx, never edits in place.

The two notes
-------------
1. slide 15 (【补】 将第 3 部分……标题格式与位置改为与本页相同，无大纲标题)
   Slide 15's title box is the reference: a single 32 pt line, no 24 pt
   grey "3.x　…" outline label, box at (1.05", 0.57", 11.78", 0.58").
   Slides 16-22 (the rest of part 3: 3.2 药效预测 ×5, 3.3 湿实验验证 ×2)
   still carry the 24 pt outline label; drop it and snap the box to the
   slide-15 geometry.

2. slide 24 (【补】 v2 benchmark 中，warm-start 场景下的 PCC, Spearman,
   Pos P@20, Neg P@20 …)
   "v2 benchmark" = our own union-CM finetune benchmark, run by
   stage1/xpert-fda-miri-v2/05_xpert_finetune_bench/full/ on the DToxS-CM
   data. The 【补】 asks for its warm-start PCC / Spearman / Pos P@20 / Neg P@20,
   so a bar chart of *our* model's 10 warm-start test folds is added (generated
   by make_v2_warmstart_figure.py from full/output/benchmark.csv:
   Pearson_deg, Spearman_deg, Precision@20 Positive_deg,
   Precision@20 Negative_deg). The XPert-paper Fig. 2a crop and its legend stay
   on the slide as the literature reference; both are re-laid-out into a
   left column so the two figures sit side by side.
"""

from pathlib import Path

from PIL import Image
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt

ROOT = Path("/home/chang/arc")
SRC = ROOT / "presentation/开题3r.pptx"
OUT = ROOT / "presentation/开题4.pptx"
FIG = ROOT / "presentation/figures"

_IN = 914400

# Slide-15 title geometry (EMU) — the format/position every part-3 title takes.
REF_TITLE_BOX = (959485, 520065, 10768330, 533400)

PART3_SLIDES = range(16, 23)  # 1-based: slides 16..22

# slide 24 — our own v2 (DToxS-CM) finetune benchmark chart
CHART_PNG = FIG / "v2_dtoxS_warmstart.png"

GREY = RGBColor(0x80, 0x80, 0x80)


def by_name(slide, name):
    for shp in slide.shapes:
        if shp.name == name:
            return shp
    raise KeyError(name)


def remove_shape(shape):
    shape._element.getparent().remove(shape._element)


def add_picture_at(slide, image_path, bbox):
    """Add image_path, aspect-fit and centred inside bbox (inches)."""
    l, t, w, h = (int(v * _IN) for v in bbox)
    with Image.open(image_path) as im:
        iw, ih = im.size
    s = min(w / iw, h / ih)
    nw, nh = int(iw * s), int(ih * s)
    return slide.shapes.add_picture(str(image_path),
                                    l + (w - nw) // 2, t + (h - nh) // 2, nw, nh)


def resize_keep_ratio(shape, left, top, width):
    """Move/resize an existing picture, preserving its current displayed aspect."""
    ratio = shape.height / shape.width
    shape.left, shape.top = Inches(left), Inches(top)
    shape.width = Inches(width)
    shape.height = int(Inches(width) * ratio)


def fit_picture_into(slide, image_path, bbox):
    """Add image_path, aspect-fit and centred inside bbox (inches)."""
    l, t, w, h = (int(v * _IN) for v in bbox)
    with Image.open(image_path) as im:
        iw, ih = im.size
    s = min(w / iw, h / ih)
    nw, nh = int(iw * s), int(ih * s)
    return slide.shapes.add_picture(str(image_path),
                                    l + (w - nw) // 2, t + (h - nh) // 2, nw, nh)


def add_label(slide, text, box, size=11, bold=False, color=None,
              align=PP_ALIGN.LEFT):
    tb = slide.shapes.add_textbox(Inches(box[0]), Inches(box[1]),
                                  Inches(box[2]), Inches(box[3]))
    tf = tb.text_frame
    tf.word_wrap = True
    for i, line in enumerate(text.split("\n")):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        r = p.add_run()
        r.text = line
        r.font.size = Pt(size)
        r.font.bold = bold
        if color is not None:
            r.font.color.rgb = color
    return tb


def strip_outline_title(slide):
    """Drop the 24 pt outline-label paragraph and adopt the slide-15 geometry."""
    tb = by_name(slide, "文本框 2")
    tf = tb.text_frame
    paras = list(tf.paragraphs)
    sizes = [r.font.size.pt for r in paras[0].runs if r.font.size is not None]
    if paras and sizes and max(sizes) <= 24:
        paras[0]._p.getparent().remove(paras[0]._p)
    l, t, w, h = REF_TITLE_BOX
    tb.left, tb.top, tb.width, tb.height = l, t, w, h


def set_v2_benchmark_chart(slide):
    """Our own v2 benchmark bar chart, beside the retained XPert-paper figure."""
    remove_shape(by_name(slide, "文本框 8"))  # the 【补】 todo note

    # left column — the XPert-paper Fig. 2a crop + its 10-model legend
    resize_keep_ratio(by_name(slide, "图片 6"), left=0.85, top=2.34, width=7.40)
    resize_keep_ratio(by_name(slide, "图片 9"), left=0.85, top=4.70, width=2.62)

    add_label(slide, "文献 benchmark：XPert 论文 Fig. 2a（warm-start / cold-drug / cold-cell）",
              (0.85, 2.00, 7.40, 0.28))
    add_label(slide, "数据来源：Guo Y, et al. Nat Mach Intell 2026, Fig. 2a",
              (0.85, 5.58, 7.40, 0.28), size=10, color=GREY)

    # right column — our own DToxS-CM warm-start test result
    add_label(slide, "我们自己微调模型的 v2 benchmark\nDToxS-CM · warm-start · 10 折",
              (8.62, 1.98, 3.95, 0.58), size=11, bold=True)
    fit_picture_into(slide, CHART_PNG, (8.62, 2.58, 3.95, 3.95))


def main():
    prs = Presentation(str(SRC))
    S = list(prs.slides)

    for idx in PART3_SLIDES:
        strip_outline_title(S[idx - 1])

    set_v2_benchmark_chart(S[23])  # slide 24

    prs.save(str(OUT))
    print("Saved", OUT.relative_to(ROOT))


if __name__ == "__main__":
    main()
