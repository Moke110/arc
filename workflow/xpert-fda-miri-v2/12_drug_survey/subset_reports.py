#!/usr/bin/env python3
"""Step 12 subset — 14-drug heatmap + filtered bilingual HTML report.

Restricts the step-12 evidence survey to a hand-picked subset of 14 drugs from
the top-20 candidates, then emits:

  12_drug_survey/heatmap_scores@14.jpg — score-percentile heatmap (same style
                                         as 11_integrate/heatmap_rank@20.jpg,
                                         only the 14 selected drugs)
  12_drug_survey/reports@14.html       — self-contained bilingual report with
                                         the 14 dossiers + embedded heatmap

Reuses the existing renderers rather than duplicating them:
  - heatmap: render_rank_heatmap from 11_integrate/integrate.py
  - markdown/score rendering: helpers from 12_drug_survey/assemble_reports.py

Non-destructive: the other 6 top-20 dossiers and reports@20.html are untouched.

Usage:
    python 12_drug_survey/subset_reports.py
"""

import base64
import importlib.util
import sys
from datetime import date
from pathlib import Path

import pandas as pd

STEP_DIR = Path(__file__).resolve().parent
WORKFLOW_ROOT = STEP_DIR.parent

sys.path.insert(0, str(WORKFLOW_ROOT / "_lib"))
from scoring import load_critical_names  # noqa: E402

SELECTED = [
    "DOBUTAMINE", "MIDODRINE", "ROPINIROLE", "METHOXAMINE", "TROPICAMIDE",
    "SOTALOL", "PHENOXYBENZAMINE", "CYCLOPENTOLATE", "LABETALOL", "DOXYLAMINE",
    "PRAMIPEXOLE", "PROPRANOLOL", "MIRTAZAPINE", "AMANTADINE",
]

HEATMAP = STEP_DIR / "heatmap_scores@14.jpg"
OUT = STEP_DIR / "reports@14.html"

KEYS = ["drug", "rank", "integrated", "chinese_name", "fda", "nmpa",
        "classification", "score_miri", "score_iri"]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ASSEMBLE = load_module("assemble_step12", STEP_DIR / "assemble_reports.py")


def render_heatmap(ranked, integrated, critical):
    """Reuse step-11's percentile heatmap renderer, restricted to the subset."""
    integrate = load_module("integrate_step11", WORKFLOW_ROOT / "11_integrate" / "integrate.py")
    merge_cols = (["drug_name", "WCS", "RGES", "Zhang_c", "pathway_reversal", "integrated"]
                  + [f"{n}_efficacy" for n in critical])
    subset = (ranked[ranked["drug_name"].isin(SELECTED)][merge_cols]
              .sort_values("integrated", ascending=False).reset_index(drop=True))
    integrate.STEP_DIR = STEP_DIR
    integrate.render_rank_heatmap(
        subset, integrated, critical,
        f"MIRI Drug Ranking — Selected 14 drugs ({len(subset)} drugs)",
        HEATMAP.name)
    return subset


def build_html(rows):
    today = date.today().isoformat()

    legend_rows = "".join(
        f'<tr><td class="sc"><span class="badge" style="background:{ASSEMBLE.SCORE_COLOR[s]}">'
        f'{s:+d}</span></td><td>{cn} · {en}</td></tr>\n'
        for s, (cn, en) in sorted(ASSEMBLE.SCORE_MEANING.items(), reverse=True))

    sum_rows, sections = [], []
    for i, r in enumerate(rows, start=1):
        drug, meta, issues = r["drug"], r["meta"], r["issues"]
        issues_txt = (' <span class="flag">' + " · ".join(issues) + "</span>"
                      if issues else "")
        zh = meta.get("chinese_name", "").strip() or "?"
        fda = meta.get("fda", "").strip() or "?"
        nmpa = meta.get("nmpa", "").strip() or "?"
        cls = meta.get("classification", "").strip() or "?"
        miri = ASSEMBLE.score_badge("score_miri", meta, issues)
        iri = ASSEMBLE.score_badge("score_iri", meta, issues)
        link = f'<a href="drugs/{drug}.md">{drug}</a>'
        sum_rows.append(
            f'<tr><td class="num">{i}</td>'
            f'<td class="num">{r["step11_rank"]}</td>'
            f'<td><b>{link}</b>{issues_txt}</td><td>{zh}</td><td>{fda}</td>'
            f'<td>{nmpa}</td><td>{cls}</td><td class="sc">{miri}</td>'
            f'<td class="sc">{iri}</td></tr>\n')

        body = ASSEMBLE.render_md(r["body"]) if r["body"].strip() else \
            '<p class="flag">(dossier body empty)</p>'
        sections.append(
            f'<section id="drug-{i}">\n'
            f'<h2>{i}. {drug} · {zh}</h2>\n{body}\n</section>\n')

    img_b64 = base64.b64encode(HEATMAP.read_bytes()).decode("ascii")

    return f"""<!DOCTYPE html>
<html lang="zh">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>MIRI Drug Repositioning — 14-Drug Evidence Survey (reports@14)</title>
<style>
:root {{ --ink:#1c2733; --line:#d8dee6; }}
* {{ box-sizing:border-box; }}
body {{ font-family:"Segoe UI",-apple-system,"PingFang SC","Microsoft YaHei",
  sans-serif; color:var(--ink); margin:0; background:#f5f7fa; line-height:1.55; }}
.wrap {{ max-width:1180px; margin:0 auto; padding:28px 20px 80px; }}
h1 {{ font-size:26px; margin:0 0 4px; }}
.sub {{ color:#5b6675; font-size:13.5px; margin:0 0 18px; }}
h2 {{ font-size:20px; margin-top:34px; padding-bottom:6px; border-bottom:2px
  solid var(--line); }}
h3 {{ font-size:16px; margin-top:20px; }}
table {{ border-collapse:collapse; width:100%; background:#fff; font-size:13.5px; }}
th, td {{ border:1px solid var(--line); padding:6px 8px; text-align:left;
  vertical-align:top; }}
th {{ background:#eef1f5; }}
td.sc, th.sc {{ text-align:center; }}
tr td.num {{ text-align:center; }}
.badge {{ display:inline-block; min-width:26px; padding:1px 6px; border-radius:9px;
  color:#fff; font-weight:700; text-align:center; font-size:12.5px; }}
.flag {{ color:#b3401f; font-weight:600; }}
section {{ background:#fff; border:1px solid var(--line); border-radius:6px;
  padding:14px 20px; margin-top:14px; }}
section h2 {{ border:none; margin:0 0 10px; }}
code {{ background:#eef1f5; padding:0 4px; border-radius:4px; font-size:12.5px; }}
a {{ color:#0b5394; }}
table.tbl {{ margin:10px 0 4px; }}
.heatmap {{ background:#fff; border:1px solid var(--line); border-radius:6px;
  padding:12px; text-align:center; }}
.heatmap img {{ max-width:100%; height:auto; }}
</style>
</head>
<body><div class="wrap">

<h1>MIRI 药物重定位 — 选定 14 味候选药调研报告</h1>
<p class="sub">MIRI Drug Repositioning — Selected 14-Drug Evidence Survey
(<b>reports@14</b>) · xpert-fda-miri-v2 step 12 · generated {today} ·
subset of <code>11_integrate/output/ranked_drugs.csv</code> top-20 · per-drug
biomni agent investigation · scores: MIRI/AMI evidence + other
ischemia-reperfusion evidence, each from the legend below. Full dossier per
drug: <code>12_drug_survey/drugs/{{DRUG}}.md</code>.</p>

<h2>综合评分热图 · Integrated score heatmap</h2>
<div class="heatmap">
<img alt="MIRI score-percentile heatmap for the 14 selected drugs"
 src="data:image/jpeg;base64,{img_b64}">
</div>

<h2>评分说明 · Score legend</h2>
<table class="tbl"><thead><tr><th style="width:70px">Score</th><th>Meaning</th>
</tr></thead><tbody>
{legend_rows}
</tbody></table>

<h2>汇总表 · Summary ({len(sum_rows)} drugs)</h2>
<table>
<thead><tr><th>#</th><th>Step-11<br>rank</th><th>Drug</th><th>中文名</th>
<th>FDA</th><th>NMPA</th><th>药理分类 Classification</th>
<th class="sc">心梗/MIRI<br>score</th>
<th class="sc">其他缺血再灌注<br>score_iri</th></tr></thead>
<tbody>
{''.join(sum_rows)}
</tbody></table>

<h2>单药详情 · Drug dossiers</h2>
{''.join(sections)}

</div></body></html>
"""


def main():
    ranked = pd.read_csv(WORKFLOW_ROOT / "11_integrate" / "output" / "ranked_drugs.csv")
    integrated = pd.read_csv(WORKFLOW_ROOT / "11_integrate" / "output" / "integrated_scores.csv")
    critical = load_critical_names(WORKFLOW_ROOT / "input" / "MIRI_pathways.json")

    render_heatmap(ranked, integrated, critical)

    step11_rank = {name: i for i, name in enumerate(ranked["drug_name"], start=1)}
    ordered = [n for n in ranked["drug_name"] if n in SELECTED]

    rows, n_bad = [], 0
    for drug in ordered:
        md_path = STEP_DIR / "drugs" / f"{drug}.md"
        if not md_path.exists():
            rows.append({"drug": drug, "step11_rank": step11_rank[drug],
                         "meta": {}, "body": "", "issues": ["missing dossier"]})
            n_bad += 1
            continue
        text = md_path.read_text(encoding="utf-8")
        meta = ASSEMBLE.parse_header(text)
        issues = [f"missing {k}" for k in KEYS if not meta.get(k, "").strip()]
        if drug not in meta.get("drug", ""):
            issues.append("header drug mismatch")
        rows.append({"drug": drug, "step11_rank": step11_rank[drug], "meta": meta,
                     "body": ASSEMBLE.strip_header(text), "issues": issues})
        if issues:
            n_bad += 1

    OUT.write_text(build_html(rows), encoding="utf-8")
    print(f"  Saved: {OUT.relative_to(WORKFLOW_ROOT)} ({len(rows)} drugs, {n_bad} flagged)")
    for r in rows:
        if r["issues"]:
            print(f"    FLAG {r['drug']}: {r['issues']}")
    return 1 if n_bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
