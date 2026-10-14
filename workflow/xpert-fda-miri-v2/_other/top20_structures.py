#!/usr/bin/env python
"""Render the step-11 top-20 drugs' molecular structures into a self-contained HTML.

Two classifications are produced in one page:
  1. By neurotransmitter system, ordered so pharmacologically related systems sit
     next to each other (monoamines first: adrenergic -> dopaminergic -> serotonergic/
     histaminergic; then the amino-acid glutamatergic system; then the distinct
     cholinergic system; then non-neurotransmitter drugs).
  2. By molecular scaffold (appended at the end): the shared-scaffold families found
     empirically with RDKit Murcko scaffolds + pairwise maximum-common-substructure
     (MCS), followed by the drugs with a unique scaffold.

Every card carries the drug's Chinese name. Reads the authoritative step-11 ranking
(rank + integrated score + SMILES) and renders each SMILES locally with RDKit (no
external structure lookup). Writes `top20_structure.html` (structures embedded as
base64 PNG, portable) next to this script, plus the per-drug PNGs under
`output/top20_structures/`.

Run:  /home/chang/miniconda3/envs/arc/bin/python _other/top20_structures.py
"""
from __future__ import annotations

import base64
import csv
import html
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

HERE = Path(__file__).resolve().parent
WF = HERE.parent
RANKED = WF / "11_integrate" / "output" / "ranked_drugs.csv"
OUT_HTML = HERE / "top20_structure.html"
OUT_PNG_DIR = HERE / "output" / "top20_structures"

TOP_N = 20
IMG_W, IMG_H = 360, 280

# Standard Chinese INN names for the 20 drugs.
CN_NAME = {
    "DOBUTAMINE": "多巴酚丁胺",
    "MIDODRINE": "米多君",
    "MOLINDONE": "吗茚酮",
    "ROPINIROLE": "罗匹尼罗",
    "METHOXAMINE": "甲氧明",
    "TROPICAMIDE": "托吡卡胺",
    "SOTALOL": "索他洛尔",
    "DILEVALOL": "地来洛尔",
    "PHENOXYBENZAMINE": "酚苄明",
    "CYCLOPENTOLATE": "环喷托酯",
    "CARISOPRODOL": "卡立普多",
    "REMOXIPRIDE": "瑞莫必利",
    "MECAMYLAMINE": "美卡拉明",
    "LABETALOL": "拉贝洛尔",
    "DOXYLAMINE": "多西拉敏",
    "PRAMIPEXOLE": "普拉克索",
    "PROPRANOLOL": "普萘洛尔",
    "MIRTAZAPINE": "米氮平",
    "DIETHYLCARBAMAZINE": "乙胺嗪",
    "AMANTADINE": "金刚烷胺",
}

# --- Classification 1: neurotransmitter system -----------------------------
# Ordering keeps related systems adjacent: monoamines (adrenergic, dopaminergic,
# serotonergic/histaminergic) together, then the amino-acid glutamatergic system
# (antiparkinsonian, functionally tied to dopamine), then cholinergic, then other.
SYSTEMS = [
    {
        "cn": "肾上腺素能系统",
        "en": "Adrenergic system (α / β)",
        "color": "#2563eb",
        "note": "去甲肾上腺素 / 肾上腺素受体（α、β）激动或拮抗",
        "drugs": [
            ("DOBUTAMINE", "β₁ 受体激动剂", "正性肌力药（强心）", "β1-adrenergic agonist"),
            ("MIDODRINE", "α₁ 受体激动剂", "升压药（前药 desglymidodrine）", "α1-adrenergic agonist"),
            ("METHOXAMINE", "α₁ 受体激动剂", "血管收缩 / 升压", "α1-adrenergic agonist"),
            ("PHENOXYBENZAMINE", "α 受体拮抗剂（不可逆）", "α₁/α₂ 阻断，嗜铬细胞瘤", "irreversible α-adrenergic antagonist"),
            ("PROPRANOLOL", "β 受体拮抗剂（非选择性）", "β₁/β₂ 阻断，抗心律失常", "non-selective β-adrenergic antagonist"),
            ("SOTALOL", "β 受体拮抗剂 + K⁺ 通道阻滞", "Ⅲ 类抗心律失常药", "β-adrenergic antagonist / class III antiarrhythmic"),
            ("DILEVALOL", "β 受体拮抗剂（β₁ 部分激动）", "拉贝洛尔活性代谢物", "β-adrenergic antagonist (β1 partial agonist)"),
            ("LABETALOL", "α₁ + β 受体拮抗剂", "联合阻断，妊娠期高血压", "combined α1 + β-adrenergic antagonist"),
        ],
    },
    {
        "cn": "多巴胺能系统",
        "en": "Dopaminergic system (D₂ / D₃)",
        "color": "#ea580c",
        "note": "多巴胺 D₂/D₃ 受体激动或拮抗",
        "drugs": [
            ("ROPINIROLE", "D₂/D₃ 受体激动剂", "帕金森病 / 不宁腿综合征", "dopamine D2/D3 agonist"),
            ("PRAMIPEXOLE", "D₂/D₃ 受体激动剂", "帕金森病 / 不宁腿综合征", "dopamine D2/D3 agonist"),
            ("MOLINDONE", "D₂ 受体拮抗剂", "抗精神病药（吲哚酮类）", "dopamine D2 antagonist antipsychotic"),
            ("REMOXIPRIDE", "D₂ 受体拮抗剂", "抗精神病药（苯甲酰胺类）", "dopamine D2 antagonist antipsychotic"),
        ],
    },
    {
        "cn": "谷氨酸能系统",
        "en": "Glutamatergic system (NMDA)",
        "color": "#7c3aed",
        "note": "NMDA 受体拮抗；兼具多巴胺能 / 抗胆碱作用（抗帕金森）",
        "drugs": [
            ("AMANTADINE", "NMDA 受体拮抗剂", "抗帕金森 / 抗病毒（兼促多巴胺释放）", "NMDA antagonist (antiparkinsonian, dopaminergic)"),
        ],
    },
    {
        "cn": "5-羟色胺能 / 组胺能系统",
        "en": "Serotonergic / Histaminergic system",
        "color": "#0d9488",
        "note": "5-HT 受体、H₁ 组胺受体拮抗（抗抑郁 / 抗组胺）",
        "drugs": [
            ("MIRTAZAPINE", "α₂ + 5-HT₂/5-HT₃ 拮抗剂", "去甲肾上腺素能和特异性 5-羟色胺能抗抑郁药（NaSSA）", "α2 / 5-HT2 / 5-HT3 antagonist (NaSSA)"),
            ("DOXYLAMINE", "H₁ 组胺受体拮抗剂", "第一代抗组胺药，兼显著抗胆碱作用", "H1 antihistamine (with antimuscarinic action)"),
        ],
    },
    {
        "cn": "胆碱能系统",
        "en": "Cholinergic system (muscarinic / nicotinic)",
        "color": "#16a34a",
        "note": "毒蕈碱（M）或烟碱（N）型乙酰胆碱受体拮抗",
        "drugs": [
            ("TROPICAMIDE", "毒蕈碱受体拮抗剂", "睫状肌麻痹 / 散瞳（抗胆碱）", "muscarinic antagonist (mydriatic)"),
            ("CYCLOPENTOLATE", "毒蕈碱受体拮抗剂", "睫状肌麻痹 / 散瞳（抗胆碱）", "muscarinic antagonist (mydriatic)"),
            ("MECAMYLAMINE", "烟碱受体拮抗剂", "神经节阻滞剂（降压）", "nicotinic antagonist (ganglionic blocker)"),
        ],
    },
    {
        "cn": "其他（非神经递质靶点）",
        "en": "Other / non-neurotransmitter targets",
        "color": "#6b7280",
        "note": "主要作用不直接经经典神经递质受体",
        "drugs": [
            ("CARISOPRODOL", "中枢性骨骼肌松弛剂", "经代谢物甲丙氨酯作用于 GABA_A（镇静）", "centrally-acting skeletal muscle relaxant (GABA-A related)"),
            ("DIETHYLCARBAMAZINE", "抗寄生虫药", "丝虫病（淋巴丝虫）", "antiparasitic (antifilarial)"),
        ],
    },
]

# --- Classification 2: molecular scaffold ----------------------------------
# Families found empirically (RDKit Murcko scaffold + pairwise MCS). The note on
# each member records the observed scaffold relationship.
SCAFFOLDS = [
    {
        "cn": "芳基烷醇胺 / 芳基烷基胺骨架",
        "en": "Arylethanolamine / arylethylamine (adrenergic ligand scaffold)",
        "color": "#0f766e",
        "note": "Ar–C(OH)–C–N 肾上腺素能配体骨架：6 个含 β-羟基的芳基烷醇胺 + 1 个儿茶酚苯乙胺。"
                "其中 MIDODRINE↔METHOXAMINE（MCS 0.93）、DILEVALOL↔LABETALOL（MCS 1.00）为近同一骨架。",
        "members": [
            ("MIDODRINE", "2,5-二甲氧基苯乙醇胺；与 METHOXAMINE 近同一骨架"),
            ("METHOXAMINE", "2,5-二甲氧基苯乙醇胺；与 MIDODRINE 近同一骨架"),
            ("DOBUTAMINE", "儿茶酚苯乙胺（无 β-OH）；与 LABETALOL 骨架高度重叠"),
            ("DILEVALOL", "4-苯基丁-2-基氨基乙醇；与 LABETALOL 骨架一致（MCS 1.00）"),
            ("LABETALOL", "4-苯基丁-2-基氨基乙醇；与 DILEVALOL 为立体同系"),
            ("SOTALOL", "芳基乙醇胺–甲磺酰苯胺"),
            ("PROPRANOLOL", "芳氧丙醇胺（aryloxypropanolamine）"),
        ],
    },
    {
        "cn": "苯基乙酸酰胺 / 酯骨架",
        "en": "Phenylacetic acid amide / ester (antimuscarinic mydriatic scaffold)",
        "color": "#b45309",
        "note": "Ph–CH–C(=O)–[N/O] 抗胆碱散瞳药骨架（托吡卡胺为酰胺、环喷托酯为酯）。",
        "members": [
            ("TROPICAMIDE", "苯基乙酰胺（–C(=O)N–）"),
            ("CYCLOPENTOLATE", "苯基乙酸酯（–C(=O)O–）"),
        ],
    },
    {
        "cn": "吲哚酮类稠合氮杂环酮骨架",
        "en": "Indolone-type fused N-heterocyclic ketone",
        "color": "#9333ea",
        "note": "稠合吡咯/吲哚环 + 环酮；两者为同类骨架但非同一骨架，且均为多巴胺能药物。",
        "members": [
            ("ROPINIROLE", "吲哚-2-酮（oxindole）骨架"),
            ("MOLINDONE", "4-氧代四氢吲哚（tetrahydroindol-4-one）骨架"),
        ],
    },
    {
        "cn": "独特骨架（本组内无共同骨架）",
        "en": "Unique scaffolds (no shared scaffold within the top-20)",
        "color": "#6b7280",
        "note": "其余 9 个药物各自为独立骨架，未在 top-20 内找到共同骨架伙伴。",
        "members": [
            ("PRAMIPEXOLE", "氨基苯并噻唑 / 噻唑-二胺"),
            ("PHENOXYBENZAMINE", "β-卤代烷胺（芳氧丙胺型氮芥类 α 阻断剂）"),
            ("MECAMYLAMINE", "双环脂肪仲胺"),
            ("DOXYLAMINE", "苯基–吡啶基乙醇胺醚"),
            ("CARISOPRODOL", "氨基甲酸酯（甲丙氨酯类）"),
            ("REMOXIPRIDE", "苯甲酰胺–吡咯烷"),
            ("MIRTAZAPINE", "吡啶并苯并二氮杂䓬 / 哌嗪稠环"),
            ("DIETHYLCARBAMAZINE", "哌嗪甲酰胺"),
            ("AMANTADINE", "金刚烷笼状脂肪胺"),
        ],
    },
]


def read_top20() -> dict[str, dict]:
    """Return {drug_name: {rank, integrated, integrated_zs, smiles, drug_idx}}."""
    with RANKED.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    out: dict[str, dict] = {}
    for i, r in enumerate(rows[:TOP_N], start=1):
        out[r["drug_name"].strip().upper()] = {
            "rank": i,
            "integrated": float(r["integrated"]),
            "integrated_zs": float(r["integrated_zs"]),
            "smiles": r["smiles"],
            "drug_idx": r["drug_idx"],
        }
    return out


def render_png(smiles: str) -> bytes:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"unparseable SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    AllChem.Compute2DCoords(mol)
    mol = Chem.RemoveHs(mol)
    d = rdMolDraw2D.MolDraw2DCairo(IMG_W, IMG_H)
    opts = d.drawOptions()
    opts.bondLineWidth = 2
    opts.padding = 0.06
    opts.minFontSize = 11
    opts.maxFontSize = 16
    opts.clearBackground = True
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    return d.GetDrawingText()


def name_block(name: str, rank: int, color: str) -> str:
    cn = CN_NAME[name]
    return (f'<div class="name"><span class="cn">{html.escape(cn)}</span>'
            f'<span class="en">{html.escape(name)}</span>'
            f'<span class="rank" style="color:{color}">#{rank}</span></div>')


def main() -> None:
    ranked = read_top20()
    OUT_PNG_DIR.mkdir(parents=True, exist_ok=True)

    # Render each structure exactly once, saving the PNG and caching base64.
    png_b64: dict[str, str] = {}
    for name, info in ranked.items():
        png = render_png(info["smiles"])
        (OUT_PNG_DIR / f"{info['rank']:02d}_{name.lower()}.png").write_bytes(png)
        png_b64[name] = base64.b64encode(png).decode("ascii")

    def thumb(name: str) -> str:
        return (f'<div class="thumb"><img alt="{html.escape(name)} structure" '
                f'src="data:image/png;base64,{png_b64[name]}"></div>')

    # ---- System sections ---------------------------------------------------
    seen: set[str] = set()
    system_sections = []
    for s in SYSTEMS:
        cards = []
        for name, cls_cn, use_cn, mech_en in s["drugs"]:
            if name not in ranked:
                raise KeyError(f"{name} not in top-20 ranking")
            seen.add(name)
            info = ranked[name]
            cards.append(
                f'<article class="card" style="--c:{s["color"]}">'
                f'{thumb(name)}'
                f'<div class="body">'
                f'{name_block(name, info["rank"], s["color"])}'
                f'<div class="score">integrated = {info["integrated"]:.3f} '
                f'<span class="zs">(z = {info["integrated_zs"]:.2f})</span></div>'
                f'<div class="row"><span class="tag">{html.escape(cls_cn)}</span></div>'
                f'<div class="use">{html.escape(use_cn)}</div>'
                f'<div class="mech">{html.escape(mech_en)}</div>'
                f'</div></article>'
            )
        system_sections.append(
            f'<section class="system" style="--c:{s["color"]}">'
            f'<header><h2>{html.escape(s["cn"])}'
            f'<span class="en">{html.escape(s["en"])}</span>'
            f'<span class="count">{len(cards)} 个药物</span></h2>'
            f'<p class="sysnote">{html.escape(s["note"])}</p></header>'
            f'<div class="grid">{"".join(cards)}</div></section>'
        )

    missing = set(ranked) - seen
    if missing:
        raise KeyError(f"unclassified drugs: {sorted(missing)}")

    # ---- Scaffold sections -------------------------------------------------
    scaf_seen: set[str] = set()
    scaffold_sections = []
    for s in SCAFFOLDS:
        cards = []
        for name, rel in s["members"]:
            if name not in ranked:
                raise KeyError(f"{name} not in top-20 ranking")
            scaf_seen.add(name)
            info = ranked[name]
            cards.append(
                f'<article class="card" style="--c:{s["color"]}">'
                f'{thumb(name)}'
                f'<div class="body">'
                f'{name_block(name, info["rank"], s["color"])}'
                f'<div class="score">integrated = {info["integrated"]:.3f} '
                f'<span class="zs">(z = {info["integrated_zs"]:.2f})</span></div>'
                f'<div class="rel">{html.escape(rel)}</div>'
                f'</div></article>'
            )
        scaffold_sections.append(
            f'<section class="system" style="--c:{s["color"]}">'
            f'<header><h2>{html.escape(s["cn"])}'
            f'<span class="en">{html.escape(s["en"])}</span>'
            f'<span class="count">{len(cards)} 个药物</span></h2>'
            f'<p class="sysnote">{html.escape(s["note"])}</p></header>'
            f'<div class="grid">{"".join(cards)}</div></section>'
        )

    scaf_missing = set(ranked) - scaf_seen
    if scaf_missing:
        raise KeyError(f"unassigned to a scaffold group: {sorted(scaf_missing)}")

    legend = "".join(
        f'<span class="legend-item"><i style="background:{s["color"]}"></i>'
        f'{html.escape(s["cn"])} <small>{html.escape(s["en"])}</small></span>'
        for s in SYSTEMS
    )
    scaf_legend = "".join(
        f'<span class="legend-item"><i style="background:{s["color"]}"></i>'
        f'{html.escape(s["cn"])} <small>{len(s["members"])} 个药物</small></span>'
        for s in SCAFFOLDS
    )

    doc = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>MIRI Top-20 药物分子结构 — 按神经递质系统 / 分子骨架分类</title>
<style>
:root {{ color-scheme: light; }}
* {{ box-sizing: border-box; }}
body {{ margin: 0; padding: 32px 28px 56px; background: #f5f6f8;
  font-family: -apple-system, "Segoe UI", "Noto Sans CJK SC", "PingFang SC",
    "Microsoft YaHei", Roboto, sans-serif; color: #1f2430; }}
h1 {{ font-size: 26px; margin: 0 0 6px; }}
.sub {{ color: #5b6472; font-size: 14px; margin: 0 0 22px; }}
h2.part {{ font-size: 21px; margin: 40px 0 4px; padding-top: 18px;
  border-top: 2px solid #d8dce3; }}
p.partnote {{ color: #5b6472; font-size: 13px; margin: 0 0 18px; }}
.legend {{ display: flex; flex-wrap: wrap; gap: 10px 18px; margin-bottom: 26px;
  padding: 14px 16px; background: #fff; border: 1px solid #e3e6ec;
  border-radius: 12px; }}
.legend-item {{ font-size: 13px; display: inline-flex; align-items: center;
  gap: 7px; }}
.legend-item i {{ width: 12px; height: 12px; border-radius: 3px;
  display: inline-block; }}
.legend-item small {{ color: #7a8290; }}
.system {{ margin: 0 0 30px; }}
.system > header {{ border-left: 5px solid var(--c); padding: 2px 0 2px 12px;
  margin-bottom: 14px; }}
.system h2 {{ font-size: 18px; margin: 0; display: flex; align-items: baseline;
  flex-wrap: wrap; gap: 10px; }}
.system h2 .en {{ font-size: 13px; font-weight: 500; color: #5b6472; }}
.system h2 .count {{ font-size: 12px; font-weight: 500; color: #fff;
  background: var(--c); border-radius: 999px; padding: 2px 10px; }}
.sysnote {{ margin: 5px 0 0; font-size: 12.5px; color: #7a8290; }}
.grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(330px, 1fr));
  gap: 14px; }}
.card {{ display: flex; gap: 12px; background: #fff; border: 1px solid #e3e6ec;
  border-top: 3px solid var(--c); border-radius: 10px; padding: 12px;
  box-shadow: 0 1px 2px rgba(16,24,40,.04); }}
.thumb {{ flex: 0 0 150px; width: 150px; height: 118px; background: #fff;
  border: 1px solid #eef0f4; border-radius: 8px; display: flex;
  align-items: center; justify-content: center; overflow: hidden; }}
.thumb img {{ width: 100%; height: 100%; object-fit: contain; }}
.body {{ min-width: 0; flex: 1; }}
.name {{ display: flex; align-items: baseline; flex-wrap: wrap; gap: 6px; }}
.name .cn {{ font-weight: 700; font-size: 15.5px; letter-spacing: .01em; }}
.name .en {{ font-size: 11.5px; color: #5b6472; letter-spacing: .02em; }}
.rank {{ font-size: 11px; font-weight: 600; border: 1px solid currentColor;
  border-radius: 5px; padding: 0 5px; }}
.score {{ font-size: 12px; color: #40485a; margin: 3px 0 6px;
  font-variant-numeric: tabular-nums; }}
.score .zs {{ color: #8a92a1; }}
.row {{ margin-bottom: 5px; }}
.tag {{ font-size: 11.5px; font-weight: 600; color: #fff; background: var(--c);
  border-radius: 5px; padding: 2px 7px; }}
.use {{ font-size: 12.5px; color: #2b3242; }}
.mech {{ font-size: 11.5px; color: #8a92a1; margin-top: 2px;
  font-style: italic; }}
.rel {{ font-size: 12px; color: #2b3242; line-height: 1.45; }}
footer {{ margin-top: 34px; font-size: 12px; color: #8a92a1;
  border-top: 1px solid #e3e6ec; padding-top: 14px; }}
</style>
</head>
<body>
<h1>MIRI Top-20 药物分子结构</h1>
<p class="sub">来源：<code>workflow/xpert-fda-miri-v2/11_integrate/output/ranked_drugs.csv</code>
（integrated score 排名前 20）。结构由 RDKit 依据步骤 11 的 SMILES 本地渲染。</p>

<h2 class="part">一、按神经递质系统 / 药理学分类</h2>
<p class="partnote">同一神经递质系统的药物相邻；单胺系统（肾上腺素能 → 多巴胺能 → 5-羟色胺能/组胺能）
相邻，其后为氨基酸类谷氨酸能系统、胆碱能系统，最后为非神经递质靶点药物。</p>
<div class="legend">{legend}</div>
{"".join(system_sections)}

<h2 class="part">二、按分子骨架分类</h2>
<p class="partnote">骨架由 RDKit 的 Murcko 骨架 + 两两最大公共子结构（MCS）实证判定：
仅列出在本组内确证存在共同骨架的药物家族，其余药物各自为独特骨架。每个卡片标注其骨架关系。</p>
<div class="legend">{scaf_legend}</div>
{"".join(scaffold_sections)}

<footer>共 {len(seen)} 个药物 · 神经递质系统 {len(SYSTEMS)} 类 · 分子骨架 {len(SCAFFOLDS)} 组
（含独特骨架组）· 图像由 RDKit 本地生成（<code>_other/top20_structures.py</code>），
未使用外部结构数据库。</footer>
</body>
</html>
"""
    OUT_HTML.write_text(doc, encoding="utf-8")
    print(f"wrote {OUT_HTML}")
    print(f"wrote {len(seen)} PNGs under {OUT_PNG_DIR}")


if __name__ == "__main__":
    main()
