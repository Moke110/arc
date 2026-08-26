#!/usr/bin/env python
"""Generate 12_report/output/report_cn.md — Chinese translation of the drug report.

Changes vs report.md:
  * Full Chinese translation.
  * First (summary) table drops the NMPA column and adds an AMI/MIRI column
    formatted "n/m": n = drug AMI-improvement evidence, m = drug MIRI-improvement
    evidence, with 2=clinical improvement, 1=preclinical improvement,
    0=contradictory reports/ineffective, -1=preclinical worsening, -2=clinical worsening.

Unlike the older hand-authored CN/SCORES tables (which were frozen to one run's
final-drug set), this version is fully solution-driven: every field is built
from the cached A1 fact sheets (12_report/cache/) written by
run_drug_agents.py, exactly like assemble_report.py. The AMI/MIRI evidence
scores are a documented heuristic over the A1 "MIRI/AMI clinical association"
and "preclinical association" fields (see score_evidence()).

Covers the 12_report final-drug set (10_filter survivors with integrated-score
percentile > 0.5), read via run_drug_agents.load_final_drugs().
"""
import json
import re
from datetime import date
from pathlib import Path

import pandas as pd

from run_drug_agents import load_final_drugs

WORKFLOW_ROOT = Path(__file__).resolve().parent.parent
REPORT_DIR = WORKFLOW_ROOT / "12_report"
CACHE_DIR = REPORT_DIR / "cache"
OUT = REPORT_DIR / "output" / "report_cn.md"

FIELDS = ["ENGLISH_NAME", "CHINESE_NAME", "US_FDA_APPROVAL", "CHINA_NMPA_APPROVAL",
          "CLASSIFICATION", "CLINICAL_SCOPE", "MIRI_AMI_CLINICAL", "MIRI_AMI_PRECLINICAL"]

# Standard pharmacopeia transliterations for drugs the A1 agents left "not found"
# (biomni tools carry no Chinese-name source). Marked with ★ in the table.
CHINESE_OVERRIDE = {
    "LEVALLORPHAN": "左啡诺 (烯丙左吗喃)",
    "LOXAPINE": "洛沙平",
    "CYCLIZINE": "赛克利嗪 (环利嗪)",
    "DESIPRAMINE": "地昔帕明 (去甲丙咪嗪)",
    "ALPIDEM": "阿吡坦",
    "IBUDILAST": "异丁司特",
    "CYCLOPENTOLATE": "环喷托酯",
    "DIPHEMANIL": "地芬马尼",
    "BRETYLIUM": "溴苄胺",
    "INFIGRATINIB": "英菲格拉替尼",
}

SCORE_MEANING = {
    2: "临床改善", 1: "临床前改善", 0: "矛盾报道/无效", -1: "临床前恶化", -2: "临床恶化",
}

DETAIL_LABELS = {
    "ENGLISH_NAME": "英文名", "CHINESE_NAME": "中文名",
    "US_FDA_APPROVAL": "美国 FDA 审批", "CHINA_NMPA_APPROVAL": "中国 NMPA 审批",
    "CLASSIFICATION": "药理分类", "CLINICAL_SCOPE": "临床应用范围",
    "MIRI_AMI_CLINICAL": "MIRI/AMI 临床关联", "MIRI_AMI_PRECLINICAL": "MIRI/AMI 临床前关联",
}

# --- AMI/MIRI evidence scoring heuristic (documented) ------------------------
# n (AMI evidence, clinical) and m (MIRI evidence, preclinical) are derived
# from the A1 fields:
#   +2 clinical improvement   — clinical evidence of AMI/MIRI benefit
#   +1 preclinical improvement — animal/mechanistic evidence of benefit
#    0 no/contradictory evidence
#   -1 preclinical worsening  — animal evidence of harm
#   -2 clinical worsening     — clinical contraindication/worsening
_IMPROVE = re.compile(r"(改善|降低|减少|缩小|保护|获益|阳性|有益|缩梗死|提高|减轻|下降)", re.I)
_WORSEN = re.compile(r"(恶化|有害|禁忌|加重|禁用|增加死亡率|升高风险|反指|无获益|未改善|不推荐|不利|风险上升|毒性)", re.I)
_NONE = re.compile(r"(无|未|零|0 结果|未入前|阴性|未检索到|未发现|无证据|not found)", re.I)


def score_evidence(text, scale="clinical"):
    """Heuristic AMI/MIRI evidence score for a field string.

    The regex polarity gives a coarse signal in {-2, 0, 2}; the clinical scale
    keeps it (2=clinical improvement, -2=clinical worsening) and the
    preclinical scale halves it (1 / -1), so all five SCORE_MEANING levels are
    reachable and match the documented meaning.
    """
    if not text:
        return 0
    if _NONE.search(text) and not _IMPROVE.search(text):
        return 0
    has_imp = bool(_IMPROVE.search(text))
    has_harm = bool(_WORSEN.search(text))
    if has_imp and has_harm:
        return 0  # contradictory
    raw = 2 if has_imp else (-2 if has_harm else 0)
    if scale == "preclinical" and raw != 0:
        raw //= 2
    return raw


def parse_solution(text):
    if not text:
        return {k: "" for k in FIELDS}
    fields = {}
    for k in FIELDS:
        m = re.search(rf"(?m)^{re.escape(k)}\s*:\s*(.*)$", text)
        fields[k] = m.group(1).strip() if m else ""
    return fields


def main():
    ranked = load_final_drugs()
    ranked["integrated"] = ranked["integrated"].round(4)
    n_drugs = len(ranked)

    rows = []
    quick = []
    details = []
    n_missing = {}
    for _, r in ranked.iterrows():
        drug = r["drug_name"]
        path = CACHE_DIR / f"{drug}.json"
        if not path.exists():
            continue
        d = json.load(open(path))
        sol = parse_solution(d.get("solution") or "")
        zh = sol["CHINESE_NAME"]
        if (not zh or zh.lower() == "not found") and drug in CHINESE_OVERRIDE:
            zh = CHINESE_OVERRIDE[drug] + " ★"
        n = score_evidence(sol["MIRI_AMI_CLINICAL"], scale="clinical")
        m = score_evidence(sol["MIRI_AMI_PRECLINICAL"], scale="preclinical")
        for k in FIELDS:
            if not sol[k] or sol[k].lower() == "not found":
                n_missing[k] = n_missing.get(k, 0) + 1
        rows.append((drug, sol["ENGLISH_NAME"], zh, sol["US_FDA_APPROVAL"],
                     sol["CLASSIFICATION"], sol["CLINICAL_SCOPE"], f"{n}/{m}"))
        quick.append((drug, zh, sol["MIRI_AMI_CLINICAL"], sol["MIRI_AMI_PRECLINICAL"]))
        details.append((drug, sol, n, m))

    # ---------- summary table ----------
    head = ("| 排名 | 药物 | 英文名 | 中文名 | 美国 FDA | 药理分类 | 临床用途 | "
            "AMI/MIRI |\n| ---: | --- | --- | --- | --- | --- | --- | --- |")
    lines = [head]

    def c(x, k=80):
        s = re.sub(r"\s+", " ", str(x)).strip().replace("|", "/")
        return s[:k] + ("…" if len(s) > k else "") or "not found"

    for i, (drug, en, zh, us, cls, scope, score) in enumerate(rows, 1):
        lines.append(f"| {i} | {c(drug)} | {c(en)} | {c(zh)} | {c(us)} | "
                     f"{c(cls)} | {c(scope)} | {score} |")

    # ---------- MIRI/AMI quick view ----------
    lines.append("")
    lines.append("## MIRI/AMI 关联线索速览")
    lines.append("| 药物 | 临床证据 | 临床前证据 |")
    lines.append("|---|---|---|")
    for drug, zh, clin, pre in quick:
        lines.append(f"| {c(drug)} | {c(clin, 200)} | {c(pre, 200)} |")

    # ---------- per-drug details ----------
    lines.append("")
    lines.append("## 单药详情")
    for i, (drug, sol, n, m) in enumerate(details, 1):
        lines.append(f"\n### {i}. {drug}  —  AMI/MIRI 评分：{n}/{m}"
                     f"（{SCORE_MEANING[n]} / {SCORE_MEANING[m]}）")
        for k in FIELDS:
            v = sol[k].strip()
            if k == "CHINESE_NAME" and (not v or v.lower() == "not found") \
                    and drug in CHINESE_OVERRIDE:
                v = CHINESE_OVERRIDE[drug] + "（标准音译，A1 未在数据库中查到）"
            if not v:
                v = "not found"
            lines.append(f"- **{DETAIL_LABELS[k]}**：{v}")

    out = [
        f"# xpert-fda-miri-v2 MIRI 药物重定位 — {n_drugs} 候选药调研报告（中文）",
        "",
        "> 数据来源：`11_integrate/output/ranked_drugs.csv`（step 11 整合评分，"
        "integrated percentile > 0.5），每种药物由 biomni A1 agent 独立调研"
        "（KEGG/ChEMBL/openFDA/OpenTargets/PubMed/TxGNN），"
        f"{date.today().isoformat()} 生成。",
        "",
        "> AMI/MIRI 列 'n/m'：n = AMI 改善证据（临床），m = MIRI 改善证据"
        "（临床前），由 A1 的 MIRI/AMI 临床与临床前关联字段启发式判定："
        "2=临床改善，1=临床前改善，0=矛盾/无效，-1=临床前恶化，-2=临床恶化。",
        "",
        "## 汇总表",
    ] + lines + [""]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(out))
    print(f"Wrote {OUT} ({len(rows)} drugs, {len(details)} details)")
    print("fields empty/not found:", n_missing)


if __name__ == "__main__":
    main()
