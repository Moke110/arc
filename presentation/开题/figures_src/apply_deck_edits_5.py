#!/usr/bin/env python3
"""Apply the three 【补】 notes carried in presentation/开题4r.pptx and unify the
document's fonts, writing presentation/开题5.pptx.

    /home/chang/miniconda3/envs/dargus_e1/bin/python \
        presentation/figures_src/apply_deck_edits_5.py

Deterministic: always starts from the frozen 开题4r.pptx, never edits in place.

The three notes
---------------
1. slide 5 (【补】添加 ACC/AHA 急性冠脉综合征指南引用索引)
   Append the guideline citation to the footer: "Rao SV, et al. Circulation
   2025" (2025 ACC/AHA/ACEP/NAEMSP/SCAI Guideline for the Management of
   Patients With Acute Coronary Syndromes, Circulation 151(13), 2025,
   doi:10.1161/CIR.0000000000001309; first author Rao SV).

2. slide 7 (【补】Tahoe 100M 论文总摘要图 / 【补】这张 l1000 图换成分辨率更高的)
   * replace the 386x396 CMap/L1000 infographic with a 4x upscale (1544x1584);
   * add the Tahoe-100M overview figure (bioRxiv 2025.02.20.639398, Figure 1)
     under the text in the left column.

3. slide 9 (【补】Perillen 降尿酸文章引用索引)
   The perillen / hyperuricemia in-vivo validation is Figure 4c-d of the very
   paper already cited (Zhu J, et al. Nat Biotechnol 2021); make that explicit.

Font normalisation (page 3 is the reference)
--------------------------------------------
Page 3 states the standard: the small section label ("1.1　研究背景") is 黑体
24 pt, while the large headline is the theme major font (+mj-lt = 微软雅黑)
32 pt; body text is 等线.  Accordingly:
  * every 方正粗黑宋简体 run -> +mj-lt;
  * the "目录" large titles -> +mj-lt;
  * the cover title -> +mj-lt;
  * footers drop every explicit typeface / size / colour (they inherit, exactly
    like page 3) and use page 3's citation semantics: journal italic, year bold.
"""

import re
from pathlib import Path

from PIL import Image
from pptx import Presentation
from pptx.enum.text import PP_ALIGN
from pptx.oxml import parse_xml
from pptx.oxml.ns import nsdecls, qn
from pptx.util import Inches, Pt

ROOT = Path("/home/chang/arc")
SRC = ROOT / "presentation/开题4r.pptx"
OUT = ROOT / "presentation/开题5.pptx"
FIG = ROOT / "presentation/figures"

L1000_HI = FIG / "l1000_cmap_hi.png"
TAHOE = FIG / "tahoe100m_fig1.png"

THEME_FACE = "+mj-lt"          # theme major font (微软雅黑) — page 3's headline
FANGZHENG = "方正粗黑宋简体"

# every non-page-3 typeface -> its page-3 equivalent
FACE_MAP = {
    "方正粗黑宋简体": THEME_FACE,
    "思源黑体 CN Light": THEME_FACE,
}

# left-column box under the slide-7 text where the Tahoe figure goes (inches)
TAHOE_BOX = (0.95, 3.45, 5.95, 3.02)
TAHOE_CAPTION = (0.95, 6.53, 5.95, 0.28)
TAHOE_CAPTION_TEXT = "Tahoe-100M 数据集总览（bioRxiv 2025, Fig. 1）"


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def by_name(slide, name):
    for shp in slide.shapes:
        if shp.name == name:
            return shp
    raise KeyError(name)


def remove_shape(shape):
    shape._element.getparent().remove(shape._element)


def shape_has(slide, name):
    return any(shp.name == name for shp in slide.shapes)


def fit_picture(slide, image_path, box):
    """aspect-fit image_path, centred inside box=(l,t,w,h) inches."""
    l, t, w, h = (int(v * 914400) for v in box)
    with Image.open(image_path) as im:
        iw, ih = im.size
    s = min(w / iw, h / ih)
    nw, nh = int(iw * s), int(ih * s)
    return slide.shapes.add_picture(str(image_path),
                                    l + (w - nw) // 2, t + (h - nh) // 2, nw, nh)


def run_rpr(run):
    return run._r.find(qn("a:rPr"))


def set_face(run, face):
    """Repoint a run's latin/ea typeface (the children already exist)."""
    rPr = run_rpr(run)
    if rPr is None:
        return False
    hit = False
    for tag in ("a:latin", "a:ea"):
        el = rPr.find(qn(tag))
        if el is not None and el.get("typeface") is not None:
            el.set("typeface", face)
            for a in ("panose", "charset", "pitchFamily"):
                el.attrib.pop(a, None)
            hit = True
    return hit


def face_of(run):
    rPr = run_rpr(run)
    if rPr is None:
        return None
    el = rPr.find(qn("a:latin"))
    return el.get("typeface") if el is not None else None


def shape_runs(shape):
    """Yield every run inside a shape, descending into groups."""
    if shape.shape_type == 6:  # GROUP
        for sub in shape.shapes:
            yield from shape_runs(sub)
        return
    if shape.has_text_frame:
        for para in shape.text_frame.paragraphs:
            for run in para.runs:
                yield run


def all_shapes(shapes):
    for shp in shapes:
        if shp.shape_type == 6:
            yield shp
            yield from all_shapes(shp.shapes)
        else:
            yield shp


def strip_run_format(run, bold, italic):
    """Footer normalisation: keep only lang/altLang + b/i (page-3 footer style)."""
    r = run._r
    rPr = r.find(qn("a:rPr"))
    if rPr is None:
        rPr = parse_xml('<a:rPr %s/>' % nsdecls("a"))
        r.insert(0, rPr)
    for child in list(rPr):
        rPr.remove(child)
    for a in list(rPr.attrib):
        if a not in ("lang", "altLang"):
            del rPr.attrib[a]
    rPr.set("b", "1" if bold else "0")
    rPr.set("i", "1" if italic else "0")


def is_year(text):
    return bool(re.fullmatch(r"\s*\d{4}\s*", text))


def fix_default_faces(root):
    """Catch paragraph-mark / list-style defaults missed by the run-level pass."""
    for rpr in root.iter(qn("a:rPr"), qn("a:endParaRPr"), qn("a:defRPr")):
        for tag in ("a:latin", "a:ea", "a:cs"):
            el = rpr.find(qn(tag))
            if el is None:
                continue
            face = el.get("typeface")
            if face in FACE_MAP:
                el.set("typeface", FACE_MAP[face])
                for a in ("panose", "charset", "pitchFamily"):
                    el.attrib.pop(a, None)


def strip_footer_para_marks(shape):
    """Footer paragraph marks must inherit page-3 style too (no font/size/colour)."""
    txBody = shape.text_frame._txBody
    for rpr in txBody.iter(qn("a:endParaRPr")):
        for child in list(rpr):
            rpr.remove(child)
        for a in list(rpr.attrib):
            if a not in ("lang", "altLang"):
                del rpr.attrib[a]
        rpr.set("b", "1")
        rpr.set("i", "0")


# --------------------------------------------------------------------------- #
# footers
# --------------------------------------------------------------------------- #
def footer_of(slide):
    for shp in slide.shapes:
        try:
            if shp.is_placeholder and "FOOTER" in str(shp.placeholder_format.type):
                return shp
        except Exception:
            pass
    return None


def set_footer_segments(shape, segments):
    """Rewrite a footer as one centred paragraph of (text, italic, bold) runs."""
    tf = shape.text_frame
    tf.clear()
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.CENTER
    for text, ital, bold in segments:
        run = p.add_run()
        run.text = text
        strip_run_format(run, bold=bold, italic=ital)


def normalise_footer(shape):
    """Page-3 footer semantics: journal italic, year bold (no explicit font)."""
    for para in shape.text_frame.paragraphs:
        for run in para.runs:
            t = run.text
            yr = is_year(t)
            strip_run_format(run, bold=yr,
                             italic=(run.font.italic is True) and not yr)


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main():
    prs = Presentation(str(SRC))
    S = list(prs.slides)

    # ---- slide 5: ACC/AHA ACS guideline citation -------------------------- #
    set_footer_segments(footer_of(S[4]), [
        ("Front Physiol.", True, False), (" ", False, False),
        ("2025", False, True), ("；", False, False),
        ("Rao SV, et al. ", False, False), ("Circulation", True, False),
        (" ", False, False), ("2025", False, True),
    ])

    # ---- slide 7: higher-res L1000 + Tahoe-100M overview ------------------ #
    s7 = S[6]
    old = by_name(s7, "Picture 18")
    l, t, w, h = old.left, old.top, old.width, old.height
    remove_shape(old)
    s7.shapes.add_picture(str(L1000_HI), l, t, w, h)

    fit_picture(s7, TAHOE, TAHOE_BOX)
    cl, ct, cw, ch = (int(v * 914400) for v in TAHOE_CAPTION)
    cap = s7.shapes.add_textbox(cl, ct, cw, ch)
    cap.text_frame.word_wrap = True
    cp = cap.text_frame.paragraphs[0]
    cp.alignment = PP_ALIGN.CENTER
    crun = cp.add_run()
    crun.text = TAHOE_CAPTION_TEXT
    crun.font.size = Pt(11)

    for nm in ("文本框 7", "文本框 8"):        # the two 【补】 notes
        if shape_has(s7, nm):
            remove_shape(by_name(s7, nm))

    # ---- slide 9: explicit perillen citation ------------------------------ #
    set_footer_segments(footer_of(S[8]), [
        ("Zhu J, et al. ", False, False), ("Nat Biotechnol", True, False),
        (" ", False, False), ("2021", False, True),
        ("（perillen 降尿酸，Fig. 4c–d）", False, False),
    ])

    # ---- font normalisation (page 3 = reference) -------------------------- #
    for idx, slide in enumerate(S, start=1):
        for shp in all_shapes(slide.shapes):
            for run in shape_runs(shp):
                # (a) 方正粗黑宋简体 -> the theme major font
                if face_of(run) == FANGZHENG:
                    set_face(run, THEME_FACE)
                # (b) the large "目录" titles follow the headline font
                if run.text.strip() == "目录":
                    set_face(run, THEME_FACE)
            # (c) cover title -> headline font
            if idx == 1 and shp.name == "文本框 6":
                for run in shape_runs(shp):
                    set_face(run, THEME_FACE)

    # paragraph marks / shape list styles
    for slide in S:
        fix_default_faces(slide._element)

    # ---- footers: page-3 format everywhere -------------------------------- #
    for slide in S:
        ftr = footer_of(slide)
        if ftr is not None:
            normalise_footer(ftr)
            strip_footer_para_marks(ftr)

    prs.save(str(OUT))
    print("Saved", OUT.relative_to(ROOT))


if __name__ == "__main__":
    main()
