#!/usr/bin/env python3
"""Step 12 (12_drug_survey) — assemble reports@20.html from the per-drug dossiers.

Reads the top-20 drug ranking from step 11 (`11_integrate/output/
ranked_drugs.csv`, integrated descending) and each bilingual dossier in
`drugs/{DRUG}.md`, then emits a single self-contained bilingual HTML report:

  reports@20.html
    - scope/date header + evidence-score legend
    - summary table (rank, drug, 中文名, FDA, NMPA, 药理分类, score_miri,
      score_iri, link to the per-drug .md)
    - embedded full per-drug dossier sections (markdown rendered to HTML)

Each `drugs/{DRUG}.md` is expected to begin with a parseable key-value header
block inside an HTML comment:

    <!-- dossier
    drug: DOBUTAMINE
    rank: 1
    integrated: +1.7241
    chinese_name: 多巴酚丁胺
    fda: approved (1978; ...)
    nmpa: approved (...)
    classification: ...
    score_miri: 2
    score_iri: 1
    -->

Missing or invalid header values are flagged in the table rather than silently
dropped. Idempotent — safe to re-run after editing any dossier.

Usage:
    python 12_drug_survey/assemble_reports.py
"""

import re
import sys
from datetime import date
from pathlib import Path

import bleach
import markdown
import pandas as pd

STEP_DIR = Path(__file__).resolve().parent
WORKFLOW_ROOT = STEP_DIR.parent
RANKED = WORKFLOW_ROOT / "11_integrate" / "output" / "ranked_drugs.csv"
DRUGS_DIR = STEP_DIR / "drugs"
OUT = STEP_DIR / "reports@20.html"

TOP_N = 20

KEYS = ["drug", "rank", "integrated", "chinese_name", "fda", "nmpa",
        "classification", "score_miri", "score_iri"]

SCORE_MEANING = {
    2: ("临床有效", "clinical efficacy"),
    1: ("临床前有效（无临床报道）", "preclinical efficacy, no clinical report"),
    0: ("无报道 / 临床前报道矛盾", "no report / conflicting preclinical"),
    -1: ("临床前无效（无临床报道）", "preclinical ineffective, no clinical report"),
    -2: ("临床无效 / 有害", "clinically ineffective / harmful"),
}

SCORE_COLOR = {2: "#1a9850", 1: "#a6d96a", 0: "#969696",
               -1: "#fdb863", -2: "#d73027"}

ALLOWED_TAGS = [
    "p", "br", "h1", "h2", "h3", "h4", "h5", "h6", "strong", "em", "ul", "ol",
    "li", "table", "thead", "tbody", "tr", "th", "td", "code", "pre",
    "blockquote", "a", "span", "hr",
]


def parse_header(text):
    """Parse the leading `<!-- dossier ... -->` comment into a dict."""
    m = re.search(r"<!--\s*dossier(.*?)-->", text, re.DOTALL)
    if not m:
        return {}
    out = {}
    for line in m.group(1).splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if ":" in line:
            k, v = line.split(":", 1)
            out[k.strip()] = v.strip()
    return out


def strip_header(text):
    """Return the dossier body with the header comment removed."""
    return re.sub(r"<!--\s*dossier.*?-->\s*", "", text, count=1, flags=re.DOTALL)


def render_md(text):
    html = markdown.markdown(text, extensions=["tables", "fenced_code",
                                               "sane_lists"])
    return bleach.clean(html, tags=ALLOWED_TAGS, attributes={"a": ["href"]},
                        strip=True)


def score_badge(key, meta, issues):
    raw = meta.get(key, "").strip()
    try:
        s = int(raw)
    except (TypeError, ValueError):
        issues.append(f"missing/invalid {key}")
        return '<span class="badge" style="background:#ccc">?</span>'
    color = SCORE_COLOR.get(s, "#ccc")
    meaning = SCORE_MEANING.get(s, ("?", "?"))
    return (f'<span class="badge" title="{meaning[1]}" '
            f'style="background:{color}">{s:+d}</span>')


def build_html(rows):
    today = date.today().isoformat()
    legend_rows = "".join(
        f'<tr><td class="sc"><span class="badge" style="background:{SCORE_COLOR[s]}">'
        f'{s:+d}</span></td><td>{cn} · {en}</td></tr>\n'
        for s, (cn, en) in sorted(SCORE_MEANING.items(), reverse=True))

    sum_rows = []
    sections = []
    for r in rows:
        drug = r["drug"]
        meta = r["meta"]
        issues = r["issues"]
        issues_txt = (' <span class="flag">' + " · ".join(issues) + "</span>"
                      if issues else "")
        zh = meta.get("chinese_name", "").strip() or "?"
        fda = meta.get("fda", "").strip() or "?"
        nmpa = meta.get("nmpa", "").strip() or "?"
        cls = meta.get("classification", "").strip() or "?"
        miri = score_badge("score_miri", meta, issues)
        iri = score_badge("score_iri", meta, issues)
        link = f'<a href="drugs/{drug}.md">{drug}</a>'
        sum_rows.append(
            f'<tr><td class="num">{r["rank"]}</td><td><b>{link}</b>'
            f'{issues_txt}</td><td>{zh}</td><td>{fda}</td><td>{nmpa}</td>'
            f'<td>{cls}</td><td class="sc">{miri}</td><td class="sc">{iri}</td>'
            f'</tr>\n')

        body = render_md(r["body"]) if r["body"].strip() else \
            '<p class="flag">(dossier body empty)</p>'
        sections.append(
            f'<section id="drug-{r["rank"]}">\n'
            f'<h2>{r["rank"]}. {drug} · {zh}</h2>\n'
            f'{body}\n</section>\n')

    return f"""<!DOCTYPE html>
<html lang="zh">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>MIRI Drug Repositioning — Top-20 Evidence Survey (reports@20)</title>
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
</style>
</head>
<body><div class="wrap">

<h1>MIRI 药物重定位 — 前 20 候选药调研报告</h1>
<p class="sub">MIRI Drug Repositioning — Top-20 Evidence Survey
(<b>reports@20</b>) · xpert-fda-miri-v2 step 12 · generated {today} ·
rank source <code>11_integrate/output/ranked_drugs.csv</code> · per-drug
biomni agent investigation · scores: MIRI/AMI evidence + other
ischemia-reperfusion evidence, each from the legend below. Full dossier per
drug: <code>12_drug_survey/drugs/{{DRUG}}.md</code>.</p>

<h2>评分说明 · Score legend</h2>
<table class="tbl"><thead><tr><th style="width:70px">Score</th><th>Meaning</th>
</tr></thead><tbody>
{legend_rows}
</tbody></table>

<h2>汇总表 · Summary ({len(sum_rows)} drugs)</h2>
<table>
<thead><tr><th>Rank</th><th>Drug</th><th>中文名</th><th>FDA</th><th>NMPA</th>
<th>药理分类 Classification</th><th class="sc">心梗/MIRI<br>score</th>
<th class="sc">其他缺血再灌注<br>score_iri</th></tr></thead>
<tbody>
{''.join(sum_rows)}
</tbody></table>

<h2>单药详情 · Drug dossiers</h2>
{''.join(sections)}

</div></body></html>
"""


def main():
    if not RANKED.exists():
        sys.exit(f"missing step-11 ranking: {RANKED}")
    ranked = pd.read_csv(RANKED).sort_values("integrated", ascending=False)
    top = ranked.head(TOP_N).reset_index(drop=True)

    DRUGS_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    n_bad = 0
    for i, (_, drug_row) in enumerate(top.iterrows(), start=1):
        drug = drug_row["drug_name"]
        md_path = DRUGS_DIR / f"{drug}.md"
        if not md_path.exists():
            rows.append({"rank": i, "drug": drug, "meta": {},
                         "body": "", "issues": ["missing dossier"]})
            n_bad += 1
            continue
        text = md_path.read_text(encoding="utf-8")
        meta = parse_header(text)
        issues = []
        for key in KEYS:
            if key not in meta or not meta[key].strip():
                issues.append(f"missing {key}")
        if drug not in meta.get("drug", ""):
            issues.append("header drug mismatch")
        body = strip_header(text)
        rows.append({"rank": i, "drug": drug, "meta": meta,
                     "body": body, "issues": issues})
        if issues:
            n_bad += 1

    OUT.write_text(build_html(rows), encoding="utf-8")
    print(f"Wrote {OUT} ({len(rows)} drugs, {n_bad} flagged)")
    for r in rows:
        if r["issues"]:
            print(f"  FLAG {r['rank']:>2} {r['drug']}: {r['issues']}")
    return 1 if n_bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
