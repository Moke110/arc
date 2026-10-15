#!/usr/bin/env python3
"""R1 (slide 4): FDA novel-drug approvals by year (2008-2024) x therapeutic area.

The right-hand companion to ``cv_pipeline_by_area.png`` (ClinicalTrials.gov
interventional-study counts by area). Same four areas, same years, same colours.

PROVENANCE / METHOD (real, citable data; no fabricated numbers)
--------------------------------------------------------------
Sources (retrieved 2026-10-13):

1. FDA openFDA "Drugs@FDA" bulk export  (authoritative FDA data)
   https://download.open.fda.gov/drug/drugsfda/drug-drugsfda-0001-of-0001.json.zip
   Each application record carries its ``submissions``. We take applications with
   an ORIGinal submission that was APproved in 2008-2024 whose
   ``submission_class_code`` is "TYPE 1" (New Molecular Entity) or "TYPE 1/4"
   (NME + new combination). That set is FDA's "novel drug" set: the per-year
   counts this produces (25,26,21,29,38,27,41,45,22,47,60,49,54,51,38,55,52 for
   2008..2024) match FDA's published CDER novel-drug-approval counts (24,26,21,
   30,39,27,41,45,22,46,59,48,53,50,37,55,50) to within +/-2 for every year.

2. WHO ATC codes joined via ChEMBL (EBI), max_phase=4, first_approval 2008-2024
   https://www.ebi.ac.uk/chembl/api/data/molecule.json?max_phase=4
   ATC level-1/chemical-subgroup classes are mapped to the four areas:
     L01* / L02* -> Oncology ; C* -> Cardiovascular ; A10* -> Diabetes ;
     J*  -> Infectious disease.
   This is the primary therapeutic-area classifier (ATC is the reference
   therapeutic classification).

3. FDA drug-label indications text (openFDA label export), used only when a drug
   has no ATC code in the four areas. Downloaded from
   https://download.open.fda.gov/drug/label/drug-label-000N-of-0014.json.zip
   and cached (per drug, longest ``indications_and_usage`` string) in
   ``figures_src/fda_approvals_label_cache.json``. A documented keyword lexicon
   maps the indication text to an area.

4. INN stem fallback (e.g. -tinib, -ciclib, -gliflozin, -glutide, -xaban,
   -cycline, -conazole, -vir ...) for the few drugs with neither ATC nor label.

The classification is a heuristic therapeutic-area mapping (option (b) of the
requirement). Drugs whose ATC/label place them outside the four areas (CNS,
immunology, ophthalmology, rare metabolic disease, diagnostics, ...) are
intentionally NOT counted.

Outputs
-------
  presentation/figures/fda_approvals_by_area.csv   (per-year x area counts)
  presentation/figures/fda_approvals_by_area.png   (the figure)

Run from the repo root:  python3 presentation/figures_src/make_fda_approvals_figure.py
"""

import collections
import io
import json
import re
import sys
import time
import urllib.request
import zipfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/home/chang/arc")
OUT = ROOT / "presentation/figures"
SRC = ROOT / "presentation/figures_src"
CACHE = SRC / ".fda_cache"
LABEL_CACHE = SRC / "fda_approvals_label_cache.json"

YEARS = list(range(2008, 2025))
AREAS = [
    ("Oncology", "#C0392B"),
    ("Cardiovascular", "#1F4E79"),
    ("Diabetes", "#E8A13C"),
    ("Infectious disease", "#4A8C5A"),
]

# FDA's published CDER novel-drug-approval counts (2008-2024), used as a
# spot-check of our Drugs@FDA-derived per-year totals. Source: FDA "New Drug
# Therapy Approvals" annual reports (fda.gov/drugs).
FDA_PUBLISHED_NOVEL = {
    2008: 24, 2009: 26, 2010: 21, 2011: 30, 2012: 39, 2013: 27, 2014: 41,
    2015: 45, 2016: 22, 2017: 46, 2018: 59, 2019: 48, 2020: 53, 2021: 50,
    2022: 37, 2023: 55, 2024: 50,
}

SALTS = {
    "hydrochloride", "hydrobromide", "dihydrochloride", "mesylate", "sodium",
    "potassium", "calcium", "magnesium", "acetate", "fumarate", "sulfate",
    "sulphate", "tromethamine", "besylate", "tosylate", "phosphate", "maleate",
    "succinate", "citrate", "bitartrate", "tartrate", "monohydrate", "dihydrate",
    "trihydrate", "ethanolate", "nitrate", "bromide", "chloride", "iodide",
    "ethyl", "meglumine", "arginine", "lysine", "salicylate", "stearate",
    "palmitate", "decanoate", "enantate", "undecanoate", "valerate",
    "dipropionate", "propionate", "pivoxil", "medoxomil", "kamedoxomil",
    "olamine", "fosamil", "anhydrous", "monosodium", "disodium", "hydrate",
    "oxide", "carbonate", "heparin",
}

STRONG = {
    "Oncology": [
        "neoplasm", "cancer", "carcinoma", "melanoma", "lymphoma", "leukemia",
        "leukaemia", "sarcoma", "myeloma", "glioma", "glioblastoma",
        "neuroblastoma", "blastoma", "myelodysplastic", "myelofibrosis",
        "oncolog", "malignan", "metastatic", "metastasis", "adenocarcinoma",
        "urothelial", "hepatocellular", "renal cell", "non-small cell",
        "colorectal", "gastrointestinal stromal", "basal cell", "merkel cell",
        "kaposi", "polycythemia vera", "myeloproliferative", "tumor", "tumour",
    ],
    "Cardiovascular": [
        "hypertension", "hypercholesterolemia", "hyperlipidemia", "dyslipidemia",
        "heart failure", "atrial fibrillation", "coronary", "myocardial", "angina",
        "venous thromboembolism", "deep vein thrombosis", "pulmonary embolism",
        "thromboembolism", "anticoagulant", "platelet", "cholesterol",
        "blood pressure", "atherosclerosis", "hypertriglyceridemia",
        "peripheral artery", "stroke",
    ],
    "Diabetes": ["diabetes", "glycemic control", "blood glucose", "hyperglycemia"],
    "Infectious disease": [
        "antibacterial", "antibiotic", "antiviral", "antifungal", "antimicrobial",
        "antimalarial", "malaria", "human immunodeficiency virus", "hepatitis b",
        "hepatitis c", "influenza", "tuberculosis", "anthrax", "plague",
        "smallpox", "candidiasis", "aspergillosis", "covid",
        "respiratory syncytial virus", "cytomegalovirus", "varicella", "gonorrhea",
        "typhoid", "pneumonia", "bacterial infection", "clostridium difficile",
    ],
}
STEM = {
    "Oncology": ["tinib", "ciclib", "parib", "fenib", "denib", "lisib", "sertib",
                 "zutifan", "zomib", "limab", "lumab", "ubicin", "trexed"],
    "Cardiovascular": ["sartan", "statin", "olol", "pril", "xaban", "grel",
                       "dipine", "apide", "paxar", "renone", "fiban"],
    "Diabetes": ["gliptin", "gliflozin", "glutide", "aglutide", "insulin",
                 "glitazone", "tirzepatide"],
    "Infectious disease": ["vir", "cycline", "oxacin", "floxacin", "conazole",
                           "fungin", "penem", "bactam", "mycin"],
}
PRIORITY = ["Oncology", "Infectious disease", "Diabetes", "Cardiovascular"]

DFDA_URL = ("https://download.open.fda.gov/drug/drugsfda/"
            "drug-drugsfda-0001-of-0001.json.zip")
LABEL_URL = ("https://download.open.fda.gov/drug/label/"
             "drug-label-{part:04d}-of-0014.json.zip")
CHEMBL = ("https://www.ebi.ac.uk/chembl/api/data/molecule.json"
          "?first_approval__gte=2008&first_approval__lte=2024&max_phase=4"
          "&limit=1000&offset={off}")


def norm(s):
    s = s.lower()
    s = re.sub(r"\(.*?\)", " ", s)
    s = re.sub(r"[^a-z0-9\- ]+", " ", s)
    return re.sub(r"\s+", " ", s).strip()


def strip_salt(s):
    toks = s.split()
    while len(toks) > 1 and toks[-1] in SALTS:
        toks = toks[:-1]
    return " ".join(toks)


def get_json(url, cache_name, retries=4):
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / cache_name
    if path.exists():
        return json.load(open(path))
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "research"})
            with urllib.request.urlopen(req, timeout=180) as r:
                data = json.load(r)
            json.dump(data, open(path, "w"))
            return data
        except Exception as e:  # noqa: BLE001
            print(f"  retry {i+1}/{retries} for {url}: {e}", file=sys.stderr)
            time.sleep(3 * (i + 1))
    raise RuntimeError(f"could not fetch {url}")


def download(url, cache_name, retries=4):
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / cache_name
    if path.exists():
        return path
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "research"})
            with urllib.request.urlopen(req, timeout=300) as r:
                path.write_bytes(r.read())
            return path
        except Exception as e:  # noqa: BLE001
            print(f"  retry {i+1}/{retries} for {url}: {e}", file=sys.stderr)
            time.sleep(3 * (i + 1))
    raise RuntimeError(f"could not download {url}")


# --------------------------------------------------------------------------- #
# 1. FDA novel approvals (Drugs@FDA bulk)
# --------------------------------------------------------------------------- #
def fetch_novel():
    zpath = download(DFDA_URL, "drugsfda.json.zip")
    with zipfile.ZipFile(zpath) as z:
        inner = z.namelist()[0]
        data = json.loads(z.read(inner))
    novel = []
    for r in data["results"]:
        for s in r.get("submissions", []) or []:
            if (s.get("submission_type") == "ORIG"
                    and s.get("submission_status") == "AP"
                    and s.get("submission_status_date")):
                y = int(s["submission_status_date"][:4])
                if 2008 <= y <= 2024 and s.get("submission_class_code") in ("TYPE 1", "TYPE 1/4"):
                    of = r.get("openfda", {})
                    gens = list(of.get("generic_name", []) or [])
                    if not gens:
                        for p in r.get("products", []) or []:
                            for a in p.get("active_ingredients", []) or []:
                                gens.append(a.get("name", ""))
                    brands = list(of.get("brand_name", []) or [])
                    if not brands:
                        brands = [p.get("brand_name", "") for p in r.get("products", []) or []]
                    novel.append({
                        "year": y,
                        "app": r["application_number"],
                        "generic": [g for g in gens if g],
                        "brand": [b for b in brands if b],
                    })
    return novel


# --------------------------------------------------------------------------- #
# 2. ATC map (ChEMBL / WHO ATC)
# --------------------------------------------------------------------------- #
def fetch_atc_map():
    path = CACHE / "chembl_name2atc.json"
    if path.exists():
        return json.load(open(path))
    mols, off = [], 0
    while True:
        d = get_json(CHEMBL.format(off=off), f"chembl_{off}.json")
        mols += d.get("molecules", [])
        off += 1000
        if off >= d["page_meta"]["total_count"]:
            break
    name2atc = {}
    for m in mols:
        atc = m.get("atc_classifications") or []
        if not atc:
            continue
        names = set()
        if m.get("pref_name"):
            names.add(norm(m["pref_name"]))
        for syn in m.get("molecule_synonyms", []) or []:
            if syn.get("molecule_synonym"):
                names.add(norm(syn["molecule_synonym"]))
        for n in names:
            for k in (n, strip_salt(n), n.split("-")[0].strip()):
                if k and k not in name2atc:
                    name2atc[k] = atc
    json.dump(name2atc, open(path, "w"))
    return name2atc


# --------------------------------------------------------------------------- #
# 3. label indications (openFDA labels), shipped cache or built from bulk
# --------------------------------------------------------------------------- #
def fetch_label_indications(targets):
    if LABEL_CACHE.exists():
        return json.load(open(LABEL_CACHE))
    print("Building label cache from openFDA label bulk (14 partitions)...")
    cache = {}
    dec = json.JSONDecoder()
    for part in range(1, 15):
        zpath = download(LABEL_URL.format(part=part), f"label{part}.zip")
        with zipfile.ZipFile(zpath) as z:
            text = z.read(z.namelist()[0]).decode("utf-8")
        i = text.find('"results"')
        i = text.find("[", i) + 1
        while True:
            while i < len(text) and text[i] in " \n\r\t,":
                i += 1
            if i >= len(text) or text[i] == "]":
                break
            rec, j = dec.raw_decode(text, i)
            i = j
            of = rec.get("openfda", {})
            cand = [norm(x) for x in (of.get("generic_name", []) or [])]
            cand += [norm(x) for x in (of.get("substance_name", []) or [])]
            cand += [norm(x) for x in (of.get("brand_name", []) or [])]
            hit = [c for c in cand if c in targets]
            if not hit:
                continue
            txt = " ".join(rec.get("indications_and_usage", []) or [])
            for c in hit:
                if len(txt) > len(cache.get(c, "")):
                    cache[c] = txt
        print(f"  label part {part}: cache={len(cache)}")
    json.dump(cache, open(LABEL_CACHE, "w"))
    return cache


# --------------------------------------------------------------------------- #
# classification
# --------------------------------------------------------------------------- #
def keys_of(n):
    ks = set()
    for g in n["generic"]:
        g2 = norm(g)
        for k in (g2, strip_salt(g2), g2.split("-")[0].strip()):
            if k:
                ks.add(k)
    for b in n["brand"]:
        if norm(b):
            ks.add(norm(b))
    return ks


def atc_area(n, atc_map):
    for k in keys_of(n):
        if k in atc_map:
            codes = atc_map[k]
            if any(c.startswith(("L01", "L02")) for c in codes):
                return "Oncology"
            if any(c.startswith("A10") for c in codes):
                return "Diabetes"
            if any(c.startswith("J") for c in codes):
                return "Infectious disease"
            if any(c.startswith("C") for c in codes):
                return "Cardiovascular"
            # ATC B01 = antithrombotic agents (anticoagulants/antiplatelets),
            # clinically cardiovascular (apixaban, ticagrelor, dabigatran, ...).
            if any(c.startswith("B01") for c in codes):
                return "Cardiovascular"
            return "__other__"
    return None


def label_area(n, ind_cache):
    best = ""
    for k in keys_of(n):
        if len(ind_cache.get(k, "")) > len(best):
            best = ind_cache[k]
    if not best:
        return None
    t = best.lower()
    for bad in ("non-cancer", "non cancer", "tumor necrosis factor",
                "tumour necrosis factor", "insulin-like"):
        t = t.replace(bad, " ")
    hits = {a: sum(1 for x in ks if x in t) for a, ks in STRONG.items()}
    b = max(hits, key=lambda a: (hits[a], -PRIORITY.index(a)))
    return b if hits[b] > 0 else None


def stem_area(n):
    h = {}
    for a, stems in STEM.items():
        c = 0
        for g in n["generic"]:
            g2 = strip_salt(norm(g))
            if any(g2.endswith(s) for s in stems):
                c += 1
        for b in n["brand"]:
            if any(norm(b).endswith(s) for s in stems):
                c += 1
        h[a] = c
    b = max(h, key=lambda a: (h[a], -PRIORITY.index(a)))
    return b if h[b] > 0 else None


def classify(n, atc_map, ind_cache):
    a = atc_area(n, atc_map)
    if a == "__other__":
        return None, "atc-other"
    if a:
        return a, "ATC"
    a = label_area(n, ind_cache)
    if a:
        return a, "label"
    a = stem_area(n)
    if a:
        return a, "stem"
    return None, "none"


def main():
    print("== R1 FDA novel-drug approvals by area ==")
    novel = fetch_novel()
    print(f"FDA Drugs@FDA novel approvals (TYPE 1, 2008-2024): {len(novel)}")

    atc_map = fetch_atc_map()
    targets = set()
    for n in novel:
        targets |= keys_of(n)
    ind_cache = fetch_label_indications(targets)

    per = collections.defaultdict(collections.Counter)
    total = collections.Counter()
    src = collections.Counter()
    assigned = collections.defaultdict(list)
    for n in novel:
        total[n["year"]] += 1
        a, s = classify(n, atc_map, ind_cache)
        src[s] += 1
        if a:
            per[n["year"]][a] += 1
            assigned[a].append(f"{n['year']} {(n['generic'] or n['brand'] or ['?'])[0]} [{s}]")

    # ---- print the per-year x area table (deliverable) ----
    print("\nClassification-tier usage:", dict(src))
    print("\nYEAR  Oncology  Cardiovascular  Diabetes  Infectious  assigned/total  (FDA novel)")
    for y in YEARS:
        c = per[y]
        row = [c[a] for a, _ in AREAS]
        print(f"{y}  {row[0]:>8d}  {row[1]:>14d}  {row[2]:>8d}  {row[3]:>10d}  "
              f"{sum(row):>8d}/{total[y]:<6d} ({FDA_PUBLISHED_NOVEL[y]})")
    print("\nArea totals:", {a: len(assigned[a]) for a, _ in AREAS})

    # ---- spot checks ----
    print("\nSpot-check vs FDA published novel-drug totals (2008-2024):")
    ok = sum(1 for y in YEARS if abs(total[y] - FDA_PUBLISHED_NOVEL[y]) <= 2)
    print(f"  our totals within +/-2 of FDA published count in {ok}/{len(YEARS)} years"
          f" (e.g. 2015 {total[2015]} vs {FDA_PUBLISHED_NOVEL[2015]},"
          f" 2018 {total[2018]} vs {FDA_PUBLISHED_NOVEL[2018]},"
          f" 2023 {total[2023]} vs {FDA_PUBLISHED_NOVEL[2023]})")
    for a in ("Oncology", "Infectious disease"):
        ex = [x for x in assigned[a] if x.startswith("2018 ")][:10]
        print(f"  2018 {a} examples: {'; '.join(ex)}")

    # ---- CSV ----
    csv = OUT / "fda_approvals_by_area.csv"
    with open(csv, "w") as f:
        f.write("year," + ",".join(a for a, _ in AREAS) + ",assigned,total_novel\n")
        for y in YEARS:
            c = per[y]
            row = [c[a] for a, _ in AREAS]
            f.write(f"{y}," + ",".join(str(v) for v in row) +
                    f",{sum(row)},{total[y]}\n")
    print(f"\nSaved table: {csv.relative_to(ROOT)}")

    # ---- figure (same styling as cv_pipeline_by_area.png) ----
    plt.rcParams.update({"figure.dpi": 200, "savefig.dpi": 200,
                         "font.size": 16, "axes.labelsize": 18,
                         "xtick.labelsize": 14, "ytick.labelsize": 14,
                         "legend.fontsize": 15, "axes.linewidth": 1.3})
    fig, ax = plt.subplots(figsize=(10, 6))
    for label, color in AREAS:
        ax.plot(YEARS, [per[y][label] for y in YEARS], "-o",
                color=color, lw=2.4, ms=6, label=label)
    ax.set_xlabel("Approval year")
    ax.set_ylabel("FDA novel drug approvals")
    ax.set_xticks(YEARS[::2])
    ax.grid(False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False)
    fig.tight_layout()
    png = OUT / "fda_approvals_by_area.png"
    fig.savefig(png, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved figure: {png.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
