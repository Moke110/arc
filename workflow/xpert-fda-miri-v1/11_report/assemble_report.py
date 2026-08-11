#!/usr/bin/env python
"""Assemble 11_report/output/report.md from the per-drug A1 fact sheets.

Covers the 11_report final-drug set (09_filter survivors with integrated-score
percentile > 0.5; 41 drugs in the current data). Reads the cached per-drug JSONs
(11_report/cache/) written by run_drug_agents.py.
"""
import json
import re
from datetime import date
from pathlib import Path

import pandas as pd

from run_drug_agents import load_final_drugs

WORKFLOW_ROOT = Path(__file__).resolve().parent.parent
REPORT_DIR = WORKFLOW_ROOT / "11_report"
CACHE_DIR = REPORT_DIR / "cache"
OUT = REPORT_DIR / "output" / "report.md"

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
}


def parse_solution(text):
    if not text:
        return {k: "" for k in FIELDS}
    fields = {}
    for k in FIELDS:
        m = re.search(rf"(?m)^{re.escape(k)}\s*:\s*(.*)$", text)
        if m:
            v = m.group(1).strip()
            # strip any trailing content on the same line only
            fields[k] = v
        else:
            fields[k] = ""
    return fields


def main():
    ranked = load_final_drugs()
    ranked["integrated"] = ranked["integrated"].round(4)
    n_drugs = len(ranked)

    rows = []
    detail = []
    n_missing = {}
    for _, r in ranked.iterrows():
        drug = r["drug_name"]
        path = CACHE_DIR / f"{drug}.json"
        if not path.exists():
            rows.append({"Drug": drug, "Integrated": r["integrated"]})
            detail.append((drug, "MISSING RESULT", ""))
            continue
        d = json.load(open(path))
        fields = parse_solution(d.get("solution") or "")
        chinese = fields["CHINESE_NAME"]
        if (not chinese or chinese.lower() == "not found") and drug in CHINESE_OVERRIDE:
            chinese = CHINESE_OVERRIDE[drug] + " ★"
        rows.append({
            "Drug": drug,
            "Integrated": r["integrated"],
            "English name": fields["ENGLISH_NAME"],
            "Chinese name": chinese,
            "US FDA": fields["US_FDA_APPROVAL"],
            "China NMPA": fields["CHINA_NMPA_APPROVAL"],
            "Classification": fields["CLASSIFICATION"],
            "Clinical scope": fields["CLINICAL_SCOPE"],
        })
        for k in FIELDS:
            if not fields[k] or fields[k].lower() == "not found":
                n_missing[k] = n_missing.get(k, 0) + 1
        detail.append((drug, r["integrated"], d.get("solution") or ""))

    table = "| Rank | Drug | English name | 中文名 | US FDA | China NMPA | 药理分类 | 临床用途 |\n"
    table += "|---:|---|---|---|---|---|---|---|\n"
    for i, row in enumerate(rows, 1):
        # keep each cell on a single line, escape pipes, truncate long cells
        def c(x):
            s = re.sub(r"\s+", " ", str(x)).strip().replace("|", "/")
            return s[:160] + ("…" if len(s) > 160 else "")
        table += f"| {i} | {c(row['Drug'])} | {c(row['English name'])} | {c(row['Chinese name'])} | {c(row['US FDA'])} | {c(row['China NMPA'])} | {c(row['Classification'])} | {c(row['Clinical scope'])} |\n"

    parts = []
    parts.append(f"# xpert-fda-miri-v1 MIRI 药物重定位 — {n_drugs} 候选药调研报告\n")
    parts.append("> 数据来源：`10_integrate/output/ranked_drugs.csv`（step 10 整合评分，"
                 "integrated percentile > 0.5），每种药物由 biomni A1 agent 独立调研"
                 "（KEGG/ChEMBL/openFDA/OpenTargets/PubMed/TxGNN），"
                 f"{date.today().isoformat()} 生成。\n")
    parts.append("> 说明：中文名后带 ★ 者为按药典标准音译补充（A1 数据库无中文名来源）；"
                 "审批/关联字段中 \"not found\" 表示所查公开数据库未检索到记录，并非确认不存在。\n")
    parts.append("## 汇总表\n")
    parts.append(table)
    parts.append("")
    parts.append("## MIRI/AMI 关联线索速览")
    parts.append("| Drug | 临床证据 | 临床前证据 |")
    parts.append("|---|---|---|")
    for i, (drug, integ, sol) in enumerate(detail, 1):
        if not sol:
            continue
        f = parse_solution(sol)
        def c(x):
            s = re.sub(r"\s+", " ", str(x)).strip().replace("|", "/")
            return s[:200] + ("…" if len(s) > 200 else "") or "not found"
        parts.append(f"| {drug} | {c(f['MIRI_AMI_CLINICAL'])} | {c(f['MIRI_AMI_PRECLINICAL'])} |")

    parts.append("")
    parts.append("## 单药详情\n")
    for i, (drug, integ, sol) in enumerate(detail, 1):
        f = parse_solution(sol)
        if not f or not any(f.values()):
            parts.append(f"\n### {i}. {drug} — (无结果)\n")
            continue
        parts.append(f"\n### {i}. {drug}\n")
        for k in FIELDS:
            label = {"ENGLISH_NAME": "英文名", "CHINESE_NAME": "中文名",
                     "US_FDA_APPROVAL": "美国 FDA 审批",
                     "CHINA_NMPA_APPROVAL": "中国 NMPA 审批",
                     "CLASSIFICATION": "药理分类", "CLINICAL_SCOPE": "临床应用范围",
                     "MIRI_AMI_CLINICAL": "MIRI/AMI 临床关联",
                     "MIRI_AMI_PRECLINICAL": "MIRI/AMI 临床前关联"}[k]
            v = f[k].strip()
            if k == "CHINESE_NAME" and (not v or v.lower() == "not found") \
                    and drug in CHINESE_OVERRIDE:
                v = CHINESE_OVERRIDE[drug] + "（标准音译，A1 未在数据库中查到）"
            if not v:
                v = "not found"
            parts.append(f"- **{label}**：{v}\n")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as fp:
        fp.write("\n".join(parts))
    print(f"Wrote {OUT} ({len(rows)} drugs, {len(detail)} details)")
    print("fields empty/not found:", n_missing)


if __name__ == "__main__":
    main()
