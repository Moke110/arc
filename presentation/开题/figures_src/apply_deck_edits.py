#!/usr/bin/env python3
"""Apply the 开题1.pptx 【补】 edits.

Run from the untouched backup each time so the result is deterministic:

    python presentation/figures_src/apply_deck_edits.py

Covers the 【补】 notes on slides 4, 6, 8, 9, 14, 15, 16, 17 plus a footer
normalisation (first author + journal italic + year bold, footer-default
font/size/colour).
"""

from pathlib import Path

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.util import Emu, Pt

ROOT = Path("/home/chang/arc")
DECK = ROOT / "presentation/开题1.pptx"
BAK = ROOT / "presentation/开题1.pptx.bak"
FIG = ROOT / "presentation/figures"
_IN = 914400

FOOTER_BOX = (int(0.90 * _IN), int(6.85 * _IN), int(11.53 * _IN), int(0.50 * _IN))


# --------------------------------------------------------------------------- #
# generic helpers
# --------------------------------------------------------------------------- #
def remove_shape(shape):
    shape._element.getparent().remove(shape._element)


def by_name(slide, name):
    for shp in slide.shapes:
        if shp.name == name:
            return shp
    raise KeyError(name)


def _fit(bbox, iw, ih):
    l, t, w, h = bbox
    s = min(w / iw, h / ih)
    nw, nh = int(iw * s), int(ih * s)
    return l + (w - nw) // 2, t + (h - nh) // 2, nw, nh


def replace_picture(slide, name, image_path, bbox=None):
    from PIL import Image
    shape = by_name(slide, name)
    with Image.open(image_path) as im:
        iw, ih = im.size
    if bbox is None:
        bbox = (shape.left, shape.top, shape.width, shape.height)
    l, t, w, h = _fit(bbox, iw, ih)
    remove_shape(shape)
    return slide.shapes.add_picture(str(image_path), l, t, w, h)


def add_picture_at(slide, image_path, bbox):
    from PIL import Image
    with Image.open(image_path) as im:
        iw, ih = im.size
    l, t, w, h = _fit(bbox, iw, ih)
    return slide.shapes.add_picture(str(image_path), l, t, w, h)


def set_text(shape, lines, size=14, align=PP_ALIGN.LEFT):
    tf = shape.text_frame
    tf.clear()
    tf.word_wrap = True
    for i, ln in enumerate(lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        run = p.add_run()
        run.text = ln
        run.font.size = Pt(size)


def set_footer(slide, segments, size=10):
    """segments: list of (text, italic, bold). Rebuilds the footer placeholder."""
    ftr = None
    for shp in slide.shapes:
        try:
            if shp.is_placeholder and str(shp.placeholder_format.type) == "FOOTER (15)":
                ftr = shp
                break
        except Exception:
            pass
    if ftr is None:
        return
    ftr.left, ftr.top, ftr.width, ftr.height = FOOTER_BOX
    tf = ftr.text_frame
    tf.clear()
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.CENTER
    for text, ital, bold in segments:
        r = p.add_run()
        r.text = text
        r.font.size = Pt(size)
        r.font.italic = ital
        r.font.bold = bold


def C(author, journal, year):
    """Convenience: citation tuple list for one work."""
    return [(author + " ", False, False), (journal, True, False),
            (" ", False, False), (year, True, True)]


def join_cites(*works):
    out = []
    for k, w in enumerate(works):
        if k:
            out.append(("；", False, False))
        out += w
    return out


# --------------------------------------------------------------------------- #
# slide 6 — inverted pyramid (native freeform shapes)
# --------------------------------------------------------------------------- #
def add_pyramid(slide):
    cx = 6.6 * _IN
    top = 1.62 * _IN
    h = 0.86 * _IN
    gap = 0.07 * _IN
    widths = [9.4, 7.7, 6.0, 4.3, 2.6]
    colors = ["#20456B", "#2E6DA4", "#4F91C6", "#8FB8DC", "#C9DDEF"]
    tcolors = ["#FFFFFF", "#FFFFFF", "#FFFFFF", "#173A57", "#173A57"]
    labels = ["疾病 × 药物反应空间（近乎无限）", "计算模拟", "体外试验",
              "体内试验", "临床试验"]
    for i in range(5):
        wt = widths[i] * _IN
        wb = (widths[i + 1] * _IN) if i + 1 < len(widths) else 0.7 * _IN
        yt = int(top + i * (h + gap))
        yb = int(yt + h)
        pts = [(int(cx - wt / 2), yt), (int(cx + wt / 2), yt),
               (int(cx + wb / 2), yb), (int(cx - wb / 2), yb)]
        fb = slide.shapes.build_freeform(pts[0][0], pts[0][1], scale=1.0)
        fb.add_line_segments(pts[1:], close=True)
        shp = fb.convert_to_shape()
        shp.fill.solid()
        shp.fill.fore_color.rgb = RGBColor.from_string(colors[i].lstrip("#"))
        shp.line.fill.background()
        shp.shadow.inherit = False
        tf = shp.text_frame
        tf.word_wrap = False
        tf.vertical_anchor = MSO_ANCHOR.MIDDLE
        p = tf.paragraphs[0]
        p.alignment = PP_ALIGN.CENTER
        r = p.add_run()
        r.text = labels[i]
        r.font.size = Pt(20)
        r.font.bold = True
        r.font.color.rgb = RGBColor.from_string(tcolors[i].lstrip("#"))


# --------------------------------------------------------------------------- #
# text content
# --------------------------------------------------------------------------- #
SLIDE8_LIST = [
    "代表性虚拟细胞药物扰动模型（按发表年份）",
    "scGen — 2019（VAE + 潜空间向量运算；跨细胞类型/物种外推）",
    "DLEPS — 2021（化学结构→L1000谱；疾病签名反推候选）",
    "DeepCE — 2021（GNN + 多头注意力；结构→差异表达）",
    "CIGER — 2022（化学诱导基因表达排序学习）",
    "chemCPA — 2022（CPA + 化学编码；预测未见药物）",
    "CPA — 2023（组合式扰动自编码器；剂量/细胞/组合外推）",
    "SAMS-VAE — 2023（稀疏加性机制偏移 VAE）",
    "TranSiGen — 2024（自监督表示学习；化学诱导谱重建）",
    "biolord — 2024（解耦表示；跨状态虚拟迁移）",
    "scGPT — 2024（生成式 Transformer 基础模型）",
    "State — 2026（多样本集合表征；跨情境扰动预测）",
    "XPert — 2026（生物学先验双分支 Transformer；解耦细胞态与药物诱导调控）",
]

SLIDE9_TITLE = "AI 增强“单细胞+虚拟细胞”范式预测内皮源 MIRI 治疗候选并体内验证（Zhang 2026）"
SLIDE9_DESC = [
    "基于小鼠心肌缺血再灌注（MIRI）内皮单细胞图谱，经细胞注释、hdWGCNA 与差异分析锁定内皮特异基因 S100A8；",
    "构建蛋白互作网络并用 AlphaFold 预测构象，进行虚拟基因敲除模拟与分子对接，筛得候选药穿心莲内酯（andrographolide, AG）。",
    "在小鼠 MIRI 模型中，AG 缩小梗死面积、改善 LVEF/LVFS 并减轻 IL-17/S100A8 炎症，完成体内阳性验证。",
]

SLIDE15_TEXT = [
    "将药物上下调基因集与疾病上下调基因集比较；若两者变化方向相反，提示药物可能具有心肌保护作用。",
    "1 WCS（Weighted Connectivity Score，加权连接性评分）",
    "在药物扰动谱上做 GSEA 式富集：分别计算疾病上调基因富集分 ES_up 与下调基因富集分 ES_down；",
    "仅当二者方向相反时 WCS =（ES_up − ES_down）/2，否则记为 0。WCS > 0 表示药物逆转疾病签名。",
]

SLIDE16_RGES = [
    "2 RGES（Reverse Gene Expression Score，逆转基因表达评分）",
    "RGES = |ES_up| − |ES_down|。与 WCS 的差异：WCS 在 ES_up、ES_down 同号时置 0，",
    "RGES 取绝对值之差、不置零，对变化幅度更敏感，可区分同向变化的强弱。",
]

SLIDE16_ZHANG = [
    "3 Zhang（Zhang & Gant 相似性评分，ssCMap 类）",
    "对药物差异基因按 |变化| 排秩 R，与疾病签名的有符号秩 s 逐基因相乘求和 C = Σ R·s，",
    "再除以理论上限 C_max 得归一化系数 c ∈ [−1, 1]；c < 0 表示逆转疾病签名。",
]

SLIDE17_TEXT = [
    "4 关键通路 RGES（pathway RGES）",
    "对 MIRI 每条关键通路，将其基因集在药物扰动谱上做 GSEA 富集得 ES，并方向对齐：",
    "UP 通路取 −ES、DOWN 通路取 +ES，得疗效分 efficacy；efficacy > 0 表示药物逆转该通路疾病方向。",
    "☆ 通路逆转率 = 被逆转通路数 / 通路总数（40 条）",
    "8 条关键通路：氧化应激、炎症反应、炎症反应调节、凋亡、坏死性凋亡、细胞死亡调节、钙信号、线粒体基因表达。",
    "直接评价 MIRI 关键机制轴，可解释性最强。",
]


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main():
    prs = Presentation(str(BAK))
    S = list(prs.slides)

    # ---- slide 4: CV drug-development lag figure ----
    s4 = S[3]
    remove_shape(by_name(s4, "文本框 3"))
    add_picture_at(s4, FIG / "cv_pipeline_by_area.png",
                   (int(1.30 * _IN), int(1.70 * _IN), int(10.70 * _IN), int(4.85 * _IN)))

    # ---- slide 6: inverted pyramid ----
    s6 = S[5]
    remove_shape(by_name(s6, "文本框 9"))
    add_pyramid(s6)

    # ---- slide 8: model list ----
    set_text(by_name(S[7], "Content Placeholder 3"), SLIDE8_LIST, size=12)

    # ---- slide 9: cell-specific model with in-vivo validation ----
    s9 = S[8]
    tb = by_name(s9, "文本框 8")
    tb.left, tb.top, tb.width, tb.height = (int(0.87 * _IN), int(1.62 * _IN),
                                            int(11.58 * _IN), int(0.45 * _IN))
    set_text(tb, [SLIDE9_TITLE], size=15)
    cp = by_name(s9, "Content Placeholder 3")
    cp.left, cp.top, cp.width, cp.height = (int(0.87 * _IN), int(2.12 * _IN),
                                            int(11.58 * _IN), int(1.25 * _IN))
    set_text(cp, SLIDE9_DESC, size=13)
    replace_picture(s9, "图片 7", FIG / "zhang2026_workflow.jpg",
                    (int(1.20 * _IN), int(3.45 * _IN), int(4.60 * _IN), int(3.35 * _IN)))
    replace_picture(s9, "Picture 18", FIG / "zhang2026_invivo.jpg",
                    (int(7.00 * _IN), int(3.45 * _IN), int(4.60 * _IN), int(3.35 * _IN)))

    # ---- slide 14: benchmark redraws ----
    s14 = S[13]
    replace_picture(s14, "Picture 18", FIG / "bench_groups_2.5x.png")
    replace_picture(s14, "Picture 19", FIG / "ablation_2.5x.png")
    remove_shape(by_name(s14, "文本框 8"))

    # ---- slide 15: WCS ----
    s15 = S[14]
    set_text(by_name(s15, "Content Placeholder 3"), SLIDE15_TEXT, size=14)
    replace_picture(s15, "Picture 20", FIG / "wcs_distribution.png",
                    (int(1.00 * _IN), int(3.85 * _IN), int(5.00 * _IN), int(2.90 * _IN)))

    # ---- slide 16: RGES + Zhang ----
    s16 = S[15]
    cps = [shp for shp in s16.shapes if shp.has_text_frame
           and shp.name.startswith("Content Placeholder")]
    cps.sort(key=lambda s: s.left)
    set_text(cps[0], SLIDE16_RGES, size=12)
    set_text(cps[1], SLIDE16_ZHANG, size=12)
    pics = [shp for shp in s16.shapes if shp.shape_type == 13]
    pics.sort(key=lambda s: s.left)
    replace_picture(s16, pics[0].name, FIG / "rges_distribution.png")
    replace_picture(s16, pics[1].name, FIG / "zhang_distribution.png")

    # ---- slide 17: pathway reversal + 8 critical pathways ----
    s17 = S[16]
    cp17 = by_name(s17, "Content Placeholder 3")
    cp17.left, cp17.top, cp17.width, cp17.height = (int(0.93 * _IN), int(1.50 * _IN),
                                                    int(5.05 * _IN), int(5.30 * _IN))
    set_text(cp17, SLIDE17_TEXT, size=12)
    remove_shape(by_name(s17, "Picture 20"))
    for nm in ("文本框 7", "文本框 8"):
        remove_shape(by_name(s17, nm))
    grid = ["pathway_reversal_distribution.png",
            "pathway_Inflammation.png", "pathway_Inflamm_Reg.png", "pathway_Apoptosis.png",
            "pathway_Oxid_Stress.png", "pathway_Necroptosis.png", "pathway_Cell_Death_Reg.png",
            "pathway_Ca_Signaling.png", "pathway_Mito_Gene.png"]
    gx0, gy0 = 6.20 * _IN, 1.42 * _IN
    gw, gh = 7.00 * _IN, 5.40 * _IN
    cw, ch = gw / 3, gh / 3
    for k, fn in enumerate(grid):
        r, c = divmod(k, 3)
        add_picture_at(s17, FIG / fn,
                       (int(gx0 + c * cw), int(gy0 + r * ch), int(cw), int(ch)))

    # ---- footer normalisation ----
    set_footer(S[3], join_cites(C("Van Norman GA, et al.", "JACC Basic Transl Sci", "2017"),
                                [("ClinicalTrials.gov", False, False), (" ", False, False),
                                 ("2026", True, True)]))
    set_footer(S[7], join_cites(
        C("Lotfollahi M, et al.", "Nat Methods", "2019"),
        C("Zhu J, et al.", "Nat Biotechnol", "2021"),
        C("Pham T, et al.", "Nat Mach Intell", "2021"),
        C("Hetzel L, et al.", "NeurIPS", "2022"),
        C("Lotfollahi M, et al.", "Mol Syst Biol", "2023"),
        C("Piran Z, et al.", "Nat Biotechnol", "2024"),
        C("Cui H, et al.", "Nat Methods", "2024"),
        C("Guo Y, et al.", "Nat Mach Intell", "2026")))
    set_footer(S[8], join_cites(
        C("Zhang Q, et al.", "Int J Mol Sci", "2026")))
    set_footer(S[13], join_cites(
        C("Guo Y, et al.", "Nat Mach Intell", "2026"),
        C("van Hasselt JGC, et al.", "Nat Commun", "2020"),
        C("Hansen J, et al.", "Nat Commun", "2024")))
    set_footer(S[14], join_cites(
        C("Lamb J, et al.", "Science", "2006"),
        C("Subramanian A, et al.", "Cell", "2017")))
    set_footer(S[15], join_cites(
        [("WCS：", False, False)] + C("Lamb J, et al.", "Science", "2006") +
        [(" / ", False, False)] + C("Subramanian A, et al.", "Cell", "2017"),
        [("RGES：", False, False)] + C("Chen B, et al.", "Nat Commun", "2017"),
        [("Zhang：", False, False)] + C("Zhang SD & Gant TW", "BMC Bioinformatics", "2008")))
    set_footer(S[16], join_cites(
        C("Chen B, et al.", "Nat Commun", "2017"),
        C("Subramanian A, et al.", "Cell", "2017"),
        C("Samart P, et al.", "Brief Bioinform", "2021")))
    set_footer(S[22], join_cites(
        C("Lotfollahi M, et al.", "Nat Methods", "2019"),
        C("Lotfollahi M, et al.", "Mol Syst Biol", "2023"),
        C("Hetzel L, et al.", "NeurIPS", "2022"),
        C("Piran Z, et al.", "Nat Biotechnol", "2024"),
        C("Adduri AK, et al.", "Cell", "2026")))

    prs.save(str(DECK))
    print("Saved", DECK)


if __name__ == "__main__":
    main()
