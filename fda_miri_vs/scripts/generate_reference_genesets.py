#!/usr/bin/env python3
"""
Generate MIRI reference gene sets from processed GEO datasets.

Workflow:
  1. GSE160516 (Mouse) — uses pre-mapped mouse→human L1000 CSVs
  2. GSE108940 (Rat)    — uses R homologene-mapped CSVs (geo/GSE108940_human_l1000.csv)
  3. GSE137293 (Pig)    — uses direct symbol matching (pig symbols ≅ human symbols)
                          supplemented by R homologene where available

Outputs JSON files to reference-genesets/:
  - GSE160516_MIRI.json
  - GSE108940_MIRI.json
  - GSE137293_scar_vs_border.json
  - MIRI_consensus.json        (cross-species, ≥2 agree)
  - MIRI_pathways.json         (16 curated MIRI pathway gene sets)
"""

import json
import sys
from pathlib import Path

import pandas as pd

PROJECT_DIR = Path(__file__).resolve().parents[1]
GEO_DIR = PROJECT_DIR / "geo"
OUT_DIR = PROJECT_DIR / "reference-genesets"

OUT_DIR.mkdir(parents=True, exist_ok=True)


# ── L1000 landmark genes ─────────────────────────────────────────────

def load_l1000_genes():
    l1000_path = "/home/chang/arc/xpert/metadata/l1000_gene_info_978.csv"
    df = pd.read_csv(l1000_path)
    genes = df["gene_symbol"].dropna().unique().tolist()
    print(f"L1000 landmark genes: {len(genes)}")
    return set(genes)


# ── Save helper ──────────────────────────────────────────────────────

def save_json(gene_set, fname):
    path = OUT_DIR / fname
    with open(path, 'w') as f:
        json.dump(gene_set, f, indent=2, ensure_ascii=False)
    print(f"  → Saved: {path}")
    return gene_set


# ══════════════════════════════════════════════════════════════════════
# Dataset 1: GSE160516 (Mouse) — already mapped
# ══════════════════════════════════════════════════════════════════════

def generate_gse160516(l1000_set):
    """Load pre-mapped mouse→human L1000 gene lists."""
    up_df = pd.read_csv(GEO_DIR / "GSE160516_up_l1000.csv")
    down_df = pd.read_csv(GEO_DIR / "GSE160516_down_l1000.csv")

    up_genes = up_df["human_gene"].dropna().unique().tolist()
    down_genes = down_df["human_gene"].dropna().unique().tolist()

    # Filter to L1000 + deduplicate
    up_in_l1000 = sorted(set(up_genes) & l1000_set)
    down_in_l1000 = sorted(set(down_genes) & l1000_set)

    # Remove conflicts
    conflicts = set(up_in_l1000) & set(down_in_l1000)
    for g in sorted(conflicts):
        # Keep direction with larger |logFC|
        up_fc = abs(up_df[up_df["human_gene"] == g]["logFC"].values[0])
        down_fc = abs(down_df[down_df["human_gene"] == g]["logFC"].values[0])
        if up_fc >= down_fc:
            down_in_l1000.remove(g)
        else:
            up_in_l1000.remove(g)
    if conflicts:
        print(f"  Resolved {len(conflicts)} conflicts by |logFC|")

    gs = {
        "name": "GSE160516_MIRI",
        "species": "Mus musculus",
        "source_geo": "GSE160516",
        "platform": "Affymetrix Clariom S",
        "comparison": "MIRI (pooled 6h/24h/72h IR) vs Sham",
        "description": (
            "Mouse MIRI model: 45min LAD ischemia + reperfusion (6h/24h/72h). "
            "Pooled IR vs Sham DEGs (limma, P<0.01, |logFC|>0.2), "
            "mouse→human ortholog via homologene, intersected with L1000 978."
        ),
        "up_genes": up_in_l1000,
        "down_genes": down_in_l1000,
        "n_up_l1000": len(up_in_l1000),
        "n_down_l1000": len(down_in_l1000),
    }
    print(f"GSE160516: {len(up_in_l1000)} UP + {len(down_in_l1000)} DOWN in L1000")
    return save_json(gs, "GSE160516_MIRI.json")


# ══════════════════════════════════════════════════════════════════════
# Dataset 2: GSE108940 (Rat) — R homologene mapped
# ══════════════════════════════════════════════════════════════════════

def generate_gse108940(l1000_set):
    """Load R homologene-mapped rat→human L1000 genes."""
    csv_path = GEO_DIR / "GSE108940_human_l1000.csv"
    if not csv_path.exists():
        print("GSE108940_human_l1000.csv not found — run geo/map_orthologs.R first")
        return None

    df = pd.read_csv(csv_path)
    if len(df) == 0:
        print("GSE108940: 0 genes mapped to L1000")
        return None

    up_genes = df[df["direction"] == "UP"]["human_gene"].dropna().unique().tolist()
    down_genes = df[df["direction"] == "DOWN"]["human_gene"].dropna().unique().tolist()

    up_in_l1000 = sorted(set(up_genes) & l1000_set)
    down_in_l1000 = sorted(set(down_genes) & l1000_set)

    # Remove conflicts
    conflicts = set(up_in_l1000) & set(down_in_l1000)
    for g in sorted(conflicts):
        up_fc = abs(df[(df["human_gene"] == g) & (df["direction"] == "UP")]["logFC"].values[0])
        down_fc = abs(df[(df["human_gene"] == g) & (df["direction"] == "DOWN")]["logFC"].values[0])
        if up_fc >= down_fc:
            down_in_l1000.remove(g)
        else:
            up_in_l1000.remove(g)

    gs = {
        "name": "GSE108940_MIRI",
        "species": "Rattus norvegicus",
        "source_geo": "GSE108940",
        "platform": "Agilent SurePrint G3 GE 8x60K",
        "comparison": "WT MIRI (IR) vs WT Sham",
        "description": (
            "Rat MIRI model: LAD ischemia-reperfusion. "
            "WT IR vs WT Sham DEGs (limma, P<0.05, |logFC|>0.3), "
            "rat→human ortholog via homologene R package (inTax=10116, outTax=9606), "
            "filtered to L1000 978 landmark genes."
        ),
        "up_genes": up_in_l1000,
        "down_genes": down_in_l1000,
        "n_up_l1000": len(up_in_l1000),
        "n_down_l1000": len(down_in_l1000),
    }
    print(f"GSE108940: {len(up_in_l1000)} UP + {len(down_in_l1000)} DOWN in L1000")
    if up_in_l1000:
        print(f"  UP: {up_in_l1000}")
    if down_in_l1000:
        print(f"  DOWN: {down_in_l1000}")
    return save_json(gs, "GSE108940_MIRI.json")


# ══════════════════════════════════════════════════════════════════════
# Dataset 3: GSE137293 (Pig) — direct symbol matching + R homologene
# ══════════════════════════════════════════════════════════════════════

def generate_gse137293(l1000_set):
    """
    Load pig DEG results, map to human, intersect with L1000.

    Strategy:
      - Human-readable pig gene symbols → direct match (most match human)
      - ENSSCG genes → attempt homologene mapping (via R output)
      - Comparison: Scar vs Border zone (NOT MIRI vs sham — caveat!)
    """
    deg_path = GEO_DIR / "GSE137293_deg_results.csv"
    r_map_path = GEO_DIR / "GSE137293_human_l1000.csv"

    if not deg_path.exists():
        print("GSE137293_deg_results.csv not found — skipping")
        return None

    deg = pd.read_csv(deg_path)

    # Filter: P<0.05, |log2FC|>0.3
    sig = deg[(deg["P_Value"] < 0.05) & (abs(deg["log2FC"]) > 0.3)].copy()
    if len(sig) < 10:
        sig = deg[deg["P_Value"] < 0.05].copy()
        print(f"  Relaxed to P<0.05 (any |logFC|): {len(sig)} genes")

    up = sig[sig["log2FC"] > 0]
    down = sig[sig["log2FC"] < 0]

    # Direct symbol matching for human-readable genes
    up_readable = [g for g in up["gene"] if not str(g).startswith("ENS")]
    down_readable = [g for g in down["gene"] if not str(g).startswith("ENS")]
    n_ens_up = sum(1 for g in up["gene"] if str(g).startswith("ENS"))
    n_ens_down = sum(1 for g in down["gene"] if str(g).startswith("ENS"))

    # Direct match to L1000 (case-insensitive)
    l1000_upper = {g.upper(): g for g in l1000_set}
    up_in_l1000 = set()
    down_in_l1000 = set()
    for g in up_readable:
        if g in l1000_set:
            up_in_l1000.add(g)
        elif g.upper() in l1000_upper:
            up_in_l1000.add(l1000_upper[g.upper()])
    for g in down_readable:
        if g in l1000_set:
            down_in_l1000.add(g)
        elif g.upper() in l1000_upper:
            down_in_l1000.add(l1000_upper[g.upper()])

    # Supplement with R homologene results (if any)
    if r_map_path.exists():
        r_map = pd.read_csv(r_map_path)
        if len(r_map) > 0:
            r_up = r_map[r_map["direction"] == "UP"]["human_gene"].dropna().tolist()
            r_down = r_map[r_map["direction"] == "DOWN"]["human_gene"].dropna().tolist()
            up_in_l1000.update(r_up)
            down_in_l1000.update(r_down)

    # Remove conflicts
    conflicts = up_in_l1000 & down_in_l1000
    up_in_l1000 = sorted(up_in_l1000 - conflicts)
    down_in_l1000 = sorted(down_in_l1000 - conflicts)

    gs = {
        "name": "GSE137293_scar_vs_border",
        "species": "Sus scrofa",
        "source_geo": "GSE137293",
        "platform": "Illumina NovaSeq 6000 (RNA-seq)",
        "comparison": "Scar (infarct core) vs Border zone (peri-infarct), 4wk post-I/R",
        "description": (
            "Pig I/R model: P30 pigs, 60min LAD occlusion + 4wk reperfusion. "
            "Scar vs Border zone comparison — CAVEAT: no sham control; "
            "comparison captures injury severity gradient, not disease vs healthy. "
            "DEG: scipy Welch t-test on log2(CPM+1), P<0.05, |log2FC|>0.3. "
            "Pig→human: direct gene symbol matching + homologene (inTax=9823). "
            f"Dropped {n_ens_up + n_ens_down} ENSSCG genes without clear ortholog."
        ),
        "up_genes": up_in_l1000,
        "down_genes": down_in_l1000,
        "n_up_l1000": len(up_in_l1000),
        "n_down_l1000": len(down_in_l1000),
        "n_ensscg_dropped": n_ens_up + n_ens_down,
        "_note": "Scar vs border — use with caution; not a MIRI-vs-sham comparison.",
    }
    print(f"GSE137293: {len(up_in_l1000)} UP + {len(down_in_l1000)} DOWN in L1000 "
          f"(dropped {n_ens_up + n_ens_down} ENSSCG)")
    if up_in_l1000:
        print(f"  UP: {up_in_l1000}")
    if down_in_l1000:
        print(f"  DOWN: {down_in_l1000}")
    return save_json(gs, "GSE137293_scar_vs_border.json")


# ══════════════════════════════════════════════════════════════════════
# Cross-species consensus
# ══════════════════════════════════════════════════════════════════════

def generate_consensus(datasets, l1000_set):
    """
    Cross-species consensus: genes dysregulated in same direction in ≥2 datasets.
    Excludes genes with conflicting directions across datasets.
    """
    valid = [d for d in datasets if d is not None]
    if len(valid) < 2:
        print("Need ≥2 valid datasets for consensus")
        return None

    from collections import Counter

    all_up = [set(d["up_genes"]) for d in valid]
    all_down = [set(d["down_genes"]) for d in valid]
    names = [d["name"] for d in valid]

    up_counter = Counter()
    down_counter = Counter()
    for u in all_up:
        up_counter.update(u)
    for d in all_down:
        down_counter.update(d)

    # Consensus: ≥2 datasets, no conflict
    consensus_up = set()
    consensus_down = set()
    for gene, count in up_counter.items():
        if count >= 2 and down_counter.get(gene, 0) == 0:
            consensus_up.add(gene)
    for gene, count in down_counter.items():
        if count >= 2 and up_counter.get(gene, 0) == 0:
            consensus_down.add(gene)

    # Strong consensus (all datasets)
    n = len(valid)
    strong_up = {g for g, c in up_counter.items() if c >= n and down_counter.get(g, 0) == 0}
    strong_down = {g for g, c in down_counter.items() if c >= n and up_counter.get(g, 0) == 0}

    gs = {
        "name": "MIRI_consensus",
        "species": "Multi-species consensus",
        "source_geo": " + ".join(names),
        "comparison": "MIRI/IR vs Sham/Control (cross-species)",
        "description": (
            f"Cross-species consensus from {n} datasets: {', '.join(names)}. "
            f"Genes dysregulated in same direction in ≥2/{n} datasets. "
            f"Strong consensus ({n}/{n} datasets): {len(strong_up)} UP + {len(strong_down)} DOWN."
        ),
        "contributing_datasets": names,
        "n_datasets": n,
        "consensus_threshold": 2,
        "up_genes": sorted(consensus_up),
        "down_genes": sorted(consensus_down),
        "n_up_l1000": len(consensus_up),
        "n_down_l1000": len(consensus_down),
        "strong_up_genes": sorted(strong_up),
        "strong_down_genes": sorted(strong_down),
        "n_strong_up": len(strong_up),
        "n_strong_down": len(strong_down),
    }

    print(f"\nConsensus (≥2/{n} datasets):")
    print(f"  UP: {len(consensus_up)} (strong={len(strong_up)}): {sorted(consensus_up)[:20]}...")
    print(f"  DOWN: {len(consensus_down)} (strong={len(strong_down)}): {sorted(consensus_down)[:20]}...")

    return save_json(gs, "MIRI_consensus.json")


# ══════════════════════════════════════════════════════════════════════
# MIRI Pathway Gene Sets (curated from KEGG + MIRI literature)
# ══════════════════════════════════════════════════════════════════════

def generate_pathway_genesets(l1000_set):
    """Curated MIRI-relevant pathways with direction labels + L1000 intersection."""

    pathways = [
        # ── UP-regulated in MIRI (activated by I/R) ──
        ("Apoptosis", "UP", {
            "BAX", "BCL2", "CASP3", "CASP8", "CASP9", "TP53", "CYCS",
            "APAF1", "BID", "BAD", "BAK1", "XIAP", "BIRC2", "BIRC3", "BIRC5",
            "DIABLO", "HTRA2", "ENDOG", "AIFM1", "PARP1", "DFFA", "DFFB",
            "ROCK1", "NFKB1", "RELA", "JUN", "FOS", "MAPK8", "MAPK9",
            "AKT1", "PIK3CA", "MTOR", "PTEN", "FAS", "FASLG",
            "TNF", "TNFRSF1A", "TRADD", "FADD", "CFLAR",
        }),
        ("TNF_signaling", "UP", {
            "TNF", "TNFRSF1A", "TRADD", "TRAF2", "RIPK1", "BIRC2", "BIRC3",
            "NFKB1", "RELA", "IKBKB", "IKBKG", "CHUK",
            "JUN", "FOS", "CXCL8", "CCL2", "CCL5", "CCL20",
            "CXCL1", "CXCL2", "CXCL3", "CXCL10", "IL6", "IL1B",
            "ICAM1", "VCAM1", "SELE", "MMP3", "MMP9",
            "MAPK1", "MAPK3", "MAPK8", "MAPK14",
            "PIK3CA", "AKT1", "MTOR",
            "PTGS2", "CREB1", "ATF2", "FOSL1", "EDN1", "VEGFA",
            "LIF", "CSF1", "CSF2", "SOCS3",
            "TAB1", "TAB2", "TAB3",
        }),
        ("HIF1_signaling", "UP", {
            "HIF1A", "ARNT", "EP300", "CREBBP", "VHL",
            "EGLN1", "EGLN2", "EGLN3",
            "VEGFA", "FLT1", "KDR", "SERPINE1", "EDN1",
            "NOS2", "NOS3", "EPO", "TF", "TFRC",
            "HK1", "HK2", "PFKL", "PFKP", "ALDOA", "ALDOC",
            "GAPDH", "PGK1", "PGAM1", "ENO1", "ENO2", "PKM",
            "LDHA", "LDHB", "PDK1", "BNIP3",
            "MTOR", "RPS6", "EIF4E", "EIF4EBP1",
            "AKT1", "PIK3CA", "MAPK1", "MAPK3", "STAT3",
            "CDKN1A", "HMOX1", "TIMP1", "MMP9", "IL6", "CXCL8",
            "IGF1R", "EGFR", "ERBB2", "INSR", "NPPA",
        }),
        ("NLRP3_inflammasome", "UP", {
            "NLRP3", "PYCARD", "CASP1", "IL1B", "IL18",
            "NFKB1", "RELA", "TNF", "CXCL8", "CCL2", "IL6",
            "TXNIP", "NEK7", "P2RX7", "PANX1",
            "MAPK8", "MAPK14", "HSP90AA1", "HSP90AB1", "HSPA1A", "HSPD1",
            "RIPK1", "RIPK3", "MLKL", "TAB1", "TAB2", "MAP3K7",
            "IRAK1", "IRAK4", "MYD88", "TLR2", "TLR4", "TICAM1",
            "ATG5", "ATG7", "BECN1", "SQSTM1",
            "STAT1", "JUN", "FOS",
        }),
        ("Ferroptosis", "UP", {
            "GPX4", "ACSL4", "LPCAT3", "ALOX5", "ALOX12", "ALOX15",
            "TFRC", "FTH1", "FTL", "SLC40A1", "SLC11A2",
            "HMOX1", "NCOA4", "SAT1", "GLS2",
            "SLC7A11", "SLC3A2", "NFE2L2", "KEAP1",
            "TXN", "TXNRD1", "GSS", "GSR", "GCLC", "GCLM",
            "GPX1", "CAT", "SOD1", "SOD2", "PRDX1", "PRDX6",
            "VDAC2", "VDAC3", "MAP1LC3B", "ATG5", "ATG7",
            "HSPB1", "HSPA5", "DDIT3", "XBP1", "ERN1", "ATF4",
            "AIFM1", "ENDOG", "PARP1",
            "CASP3", "BAX", "BCL2", "TP53", "CDKN1A", "MDM2",
            "PTEN", "PIK3CA", "AKT1", "MTOR",
        }),
        ("Autophagy", "UP", {
            "MTOR", "ULK1", "ULK2", "ATG13", "RB1CC1", "ATG101",
            "BECN1", "PIK3C3", "PIK3R4", "ATG14", "AMBRA1",
            "ATG3", "ATG4B", "ATG4C", "ATG4D", "ATG5", "ATG7",
            "ATG10", "ATG12", "ATG16L1", "MAP1LC3A", "MAP1LC3B",
            "GABARAP", "GABARAPL1", "GABARAPL2", "SQSTM1",
            "OPTN", "NBR1", "CALCOCO2", "TAX1BP1",
            "BNIP3", "BNIP3L", "FUNDC1", "PINK1", "PRKN",
            "LAMP1", "LAMP2", "CTSD", "CTSB", "CTSL",
            "RAB7A", "PLEKHM1", "VPS11", "VPS16", "VPS18",
            "STX17", "SNAP29", "VAMP8",
            "ITPR1", "TMBIM6", "BCL2", "BCL2L1",
            "HIF1A", "FOXO1", "FOXO3", "TFEB", "TFE3", "ZKSCAN3",
        }),
        ("Calcium_signaling", "UP", {
            "RYR2", "CASQ2", "PLN", "SLC8A1", "ATP2A2",
            "CACNA1C", "CACNA1D", "CACNA1G", "CACNA1H",
            "CACNB1", "CACNB2", "CACNB3", "CACNB4",
            "CALM1", "CALM2", "CALM3",
            "CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G", "CAMK4",
            "PPP3CA", "PPP3CB", "PPP3CC",
            "NFATC1", "NFATC2", "NFATC3", "NFATC4",
            "ITPR1", "ITPR2", "ITPR3",
            "PRKCA", "PRKCB", "PRKCG",
            "PLCB1", "PLCB2", "PLCB3", "PLCB4",
            "GNAQ", "GNA11", "GNAS",
            "ADCY1", "ADCY5", "ADCY6",
            "PDE1A", "PDE1B", "PDE1C",
            "EDN1", "EDNRA", "EDNRB",
            "AGTR1", "ADRB1", "ADRB2",
            "CHRM1", "CHRM2", "CHRM3",
            "HRC", "TRDN", "CASQ1",
        }),
        ("ECM_receptor_interaction", "UP", {
            "COL1A1", "COL1A2", "COL2A1", "COL3A1",
            "COL4A1", "COL4A2", "COL5A1", "COL5A2",
            "COL6A1", "COL6A2", "COL6A3",
            "FN1", "VTN", "TNC", "TNN", "TNR",
            "THBS1", "THBS2",
            "ITGA1", "ITGA2", "ITGA3", "ITGA4", "ITGA5", "ITGA6",
            "ITGAV", "ITGB1", "ITGB3", "ITGB4", "ITGB5", "ITGB6",
            "LAMA1", "LAMA2", "LAMA3", "LAMA4", "LAMA5",
            "LAMB1", "LAMB2", "LAMB3", "LAMC1", "LAMC2",
            "HSPG2", "SDC1", "SDC2", "SDC4", "AGRN", "DAG1",
            "CD44", "CD47", "COMP", "IBSP", "SPP1", "BGLAP",
        }),
        ("TGFbeta_signaling", "UP", {
            "TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2",
            "SMAD1", "SMAD2", "SMAD3", "SMAD4", "SMAD5", "SMAD6", "SMAD7", "SMAD9",
            "BMP2", "BMP4", "BMP6", "BMP7",
            "BMPR1A", "BMPR1B", "BMPR2",
            "ACVR1", "ACVR1B", "ACVR2A", "ACVR2B",
            "INHBA", "INHBB", "FST", "NODAL",
            "ID1", "ID2", "ID3", "ID4",
            "MYC", "JUN", "FOS", "SP1", "E2F4", "E2F5",
            "CDKN1A", "CDKN2B", "GADD45A", "GADD45B",
            "CTNNB1", "LEF1", "TCF7", "TCF7L1", "TCF7L2",
            "FN1", "COL1A1", "COL1A2", "COL3A1", "SERPINE1",
            "MMP2", "MMP9", "TIMP1", "TIMP2",
            "LTBP1", "LTBP2", "LTBP3", "THBS1", "THBS2",
            "NOG", "CHRD", "GREM1", "GREM2",
            "SKI", "SKIL", "TGIF1", "TGIF2",
        }),
        ("mTOR_signaling", "UP", {
            "MTOR", "RPTOR", "MLST8", "AKT1S1", "DEPTOR",
            "AKT1", "AKT2", "AKT3", "PIK3CA", "PIK3CB",
            "PIK3R1", "PIK3R2", "PIK3R3", "PTEN", "PDK1",
            "TSC1", "TSC2", "RHEB",
            "RRAGA", "RRAGB", "RRAGC", "RRAGD",
            "LAMTOR1", "LAMTOR2", "LAMTOR3", "LAMTOR4", "LAMTOR5",
            "RPS6KB1", "RPS6", "EIF4E", "EIF4EBP1", "EIF4B",
            "EIF4G1", "EIF4G2", "EIF4A1", "EIF4A2",
            "ULK1", "ULK2", "ATG13", "RB1CC1", "ATG101",
            "HIF1A", "SREBF1", "SREBF2", "PPARG", "FOXO1", "FOXO3",
            "INSR", "IGF1R", "IRS1", "IRS2",
            "MAPK1", "MAPK3", "TNF", "IKBKB", "NFKB1", "RELA",
            "RICTOR", "PRKAA1", "PRKAA2", "STK11", "STRADA", "STRADB",
            "DDIT4", "DDIT4L", "BNIP3", "GADD45A", "CDKN1A",
        }),
        ("p53_signaling", "UP", {
            "TP53", "MDM2", "MDM4", "CDKN1A", "CDKN2A",
            "BAX", "BCL2", "BCL2L1", "BBC3", "PMAIP1",
            "BIK", "BID", "BAD", "APAF1", "CASP3", "CASP8", "CASP9",
            "CYCS", "AIFM1", "ENDOG", "HTRA2", "DIABLO",
            "GADD45A", "GADD45B", "GADD45G",
            "CCNB1", "CCNB2", "CCND1", "CCNE1", "CDK1", "CDK2", "CDK4",
            "RB1", "E2F1", "E2F2", "E2F3",
            "SESN1", "SESN2", "SESN3", "TSC2", "PTEN",
            "TP53I3", "TIGAR", "SCO2", "FDXR", "RRM2B",
            "ZMAT3", "SHISA5", "SERPINB5", "THBS1", "IGFBP3",
            "PERP", "FAS", "TNFRSF10A", "TNFRSF10B",
            "DDB2", "XPC", "RRM2", "PCNA",
            "ATM", "ATR", "CHEK1", "CHEK2", "SFN", "GTSE1",
        }),
        ("JAK_STAT_signaling", "UP", {
            "IL6", "IL6R", "IL6ST",
            "JAK1", "JAK2", "JAK3", "TYK2",
            "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6",
            "IFNG", "IFNGR1", "IFNGR2",
            "IFNA1", "IFNAR1", "IFNAR2",
            "IL10", "IL10RA", "IL10RB",
            "CSF1", "CSF1R", "CSF2", "CSF2RA", "CSF3", "CSF3R",
            "EPO", "EPOR", "GH1", "GHR", "PRL", "PRLR", "LEP", "LEPR",
            "OSM", "OSMR", "LIF", "LIFR", "CTF1",
            "SOCS1", "SOCS2", "SOCS3", "SOCS4", "SOCS5", "SOCS6", "SOCS7", "CISH",
            "PIAS1", "PIAS2", "PIAS3", "PIAS4",
            "PTPN1", "PTPN2", "PTPN6", "PTPN11",
            "MYC", "JUN", "FOS", "CCND1", "BCL2", "BCL2L1",
            "PIK3CA", "AKT1", "MTOR", "MAPK1", "MAPK3",
            "SPRY1", "SPRY2", "SPRY4", "STAM", "STAM2",
            "CXCL8", "CCL2", "CXCL10", "ICAM1",
        }),
        ("MAPK_signaling", "UP", {
            "MAPK1", "MAPK3", "MAPK8", "MAPK9", "MAPK10",
            "MAPK11", "MAPK12", "MAPK13", "MAPK14",
            "MAP2K1", "MAP2K2", "MAP2K3", "MAP2K4",
            "MAP2K5", "MAP2K6", "MAP2K7",
            "MAP3K1", "MAP3K2", "MAP3K3", "MAP3K4", "MAP3K5",
            "MAP3K7", "MAP3K8", "RAF1", "BRAF",
            "RASGRF1", "RASGRF2", "RASGRP1", "RASGRP2", "RASGRP3",
            "RAPGEF2", "RAP1A", "RAP1B",
            "FGF1", "FGF2", "FGFR1", "FGFR2",
            "PDGFA", "PDGFB", "PDGFRA", "PDGFRB",
            "EGF", "EGFR", "TGFA", "ERBB2",
            "INS", "IGF1", "IGF1R", "NGF", "BDNF", "NTRK1", "NTRK2",
            "VEGFA", "KDR", "FLT1",
            "TNF", "TNFRSF1A", "IL1B", "IL1R1",
            "HSPA1A", "HSPA1B", "HSPA8", "HSPB1",
            "MKNK1", "MKNK2", "RPS6KA1", "RPS6KA2", "RPS6KA3",
            "ATF2", "ATF4", "JUN", "FOS", "ELK1", "ELK4",
            "MYC", "SRF", "TP53", "DDIT3", "NFATC1", "NFATC3",
            "MAX", "MEF2C", "CREB1", "ATF1", "CREB3", "CREB3L1", "CREB3L4",
            "NFKB1", "NFKB2", "RELA", "RELB",
            "TAB1", "TAB2", "TRAF2", "TRAF6",
            "PPP3CA", "PPP3CB", "PPP3CC",
            "DUSP1", "DUSP2", "DUSP3", "DUSP4", "DUSP5", "DUSP6",
            "DUSP7", "DUSP8", "DUSP9", "DUSP10", "DUSP16",
            "PPP5C", "PPM1A", "PPM1B",
            "AKT1", "AKT2", "AKT3", "PIK3CA",
            "CDC42", "RAC1", "RAC2", "RAC3", "RHOA",
            "PAK1", "PAK2", "CRK", "CRKL",
            "ARRB1", "ARRB2", "GRB2", "SOS1", "SOS2",
            "GNA12", "GNA13", "GNG12",
            "RASA1", "RASA2", "NF1",
        }),
        ("PI3K_Akt_signaling", "UP", {
            "PIK3CA", "PIK3CB", "PIK3R1", "PIK3R2", "PIK3R3",
            "AKT1", "AKT2", "AKT3", "PTEN", "PDK1", "PDPK1",
            "MTOR", "RPTOR", "RICTOR", "MLST8",
            "TSC1", "TSC2", "RHEB",
            "GSK3A", "GSK3B", "FOXO1", "FOXO3", "FOXO4",
            "BAD", "BCL2", "BCL2L1", "CASP9",
            "CDKN1A", "CDKN1B", "MDM2", "TP53",
            "CREB1", "NFKB1", "RELA",
            "MYC", "CCND1", "CCNE1", "CDK2", "CDK4", "RB1", "E2F1",
            "INSR", "IGF1R", "IRS1", "IRS2",
            "EGFR", "ERBB2", "ERBB3", "ERBB4",
            "FGFR1", "FGFR2", "FGFR3",
            "PDGFRA", "PDGFRB", "KIT", "FLT3", "CSF1R",
            "VEGFA", "KDR", "FLT1", "ANGPT1", "TEK",
            "NGF", "BDNF", "NTRK1", "NTRK2",
            "ITGA1", "ITGA2", "ITGA3", "ITGA4", "ITGA5", "ITGA6",
            "ITGAV", "ITGB1", "ITGB3", "ITGB4", "ITGB5",
            "LAMA1", "LAMA2", "LAMA3", "LAMA4", "LAMA5",
            "LAMB1", "LAMB2", "LAMB3", "LAMC1", "LAMC2",
            "FN1", "COL1A1", "COL1A2", "COL4A1", "COL4A2", "COL6A1",
            "THBS1", "THBS2", "COMP", "SPP1", "TNC",
            "PPP2R1A", "PPP2R1B", "PPP2R2A", "PPP2R2B", "PPP2R2C",
            "PPP2R3A", "PPP2R5A", "PPP2R5B", "PPP2R5C", "PPP2R5D",
            "PPP2CA", "PPP2CB",
            "HSP90AA1", "HSP90AB1", "HSP90B1",
            "CDC37", "YWHAB", "YWHAZ", "YWHAH", "YWHAQ",
            "PRKAA1", "PRKAA2", "STK11", "STRADA",
            "MAPK1", "MAPK3", "JAK1", "JAK2", "JAK3",
            "RPS6KB1", "RPS6", "EIF4E", "EIF4EBP1", "EIF4B",
            "GYS1", "GYS2", "PHLPP1", "PHLPP2",
            "IL2RG", "IL6", "CSF1", "CSF3", "EPO", "GH1", "PRL",
            "INS", "IGF1", "IGF2", "FGF1", "FGF2", "EGF", "TGFA",
            "HGF", "MET", "PCK1", "PCK2", "G6PC", "G6PC2", "G6PC3",
            "HK1", "HK2", "PFKL", "PFKP", "PKM",
            "SYK", "BTK", "CD19", "BLNK",
            "F2R", "CHRM1", "CHRM2",
            "LPAR1", "LPAR2", "LPAR3",
            "GNB1", "GNG2", "GNG3", "GNG4", "GNG5",
        }),

        # ── DOWN-regulated in MIRI (suppressed by I/R) ──
        ("Fatty_acid_oxidation", "DOWN", {
            "ACADM", "ACADS", "ACADVL", "ACADL", "ACAA2",
            "HADHA", "HADHB", "HADH", "ECHS1", "HSD17B10", "EHHADH",
            "CPT1A", "CPT1B", "CPT2", "SLC25A20",
            "ACOX1", "ACOX2", "ACOX3", "ACAA1",
            "PPARA", "PPARD", "PPARG", "PPARGC1A",
            "MLYCD", "ACACB", "FASN",
            "CD36", "FABP1", "FABP3", "FABP4", "LPL", "VLDLR",
            "SLC27A1", "SLC27A2", "SLC27A4", "SLC27A6",
            "DECR1", "ECI1", "ECI2", "IVD", "GCDH",
            "ETFA", "ETFB", "ETFDH", "ETHE1",
        }),
        ("Cardiac_muscle_contraction", "DOWN", {
            "MYH6", "MYH7", "MYL2", "MYL3", "MYL4",
            "ACTC1", "TPM1", "TNNT2", "TNNI3", "TNNC1",
            "TTN", "ACTN2", "DES", "VIM", "NEBL", "MYBPC3",
            "MYLK2", "MYLK3", "MYL6", "MYL7", "MYL9",
            "RYR2", "CASQ2", "PLN", "ATP2A2", "SLC8A1",
            "CACNA1C", "CACNB2", "SCN5A", "KCNH2", "KCNQ1",
            "KCNE1", "KCNJ2", "KCNJ3", "KCNJ5", "KCNJ11",
            "GJA1", "GJA5", "DSP", "PKP2", "DSG2", "DSC2",
            "CTNNA3", "CTNNB1", "CDH2", "JUP",
        }),
        ("Oxidative_phosphorylation", "DOWN", {
            "NDUFA1", "NDUFA2", "NDUFA3", "NDUFA4",
            "NDUFS1", "NDUFS2", "NDUFS3",
            "COX4I1", "COX5A", "COX5B", "COX6A1", "COX6B1",
            "COX7A2", "COX7C", "COX8A",
            "CYCS", "UQCRB", "UQCRC1", "UQCRC2", "UQCRH", "UQCRQ",
            "ATP5A1", "ATP5B", "ATP5C1", "ATP5D", "ATP5E", "ATP5F1",
            "ATP5G1", "ATP5G2", "ATP5G3", "ATP5H", "ATP5I", "ATP5J",
            "ATP5L", "ATP5O",
            "SDHA", "SDHB", "SDHC", "SDHD",
        }),
    ]

    valid = []
    for name, direction, genes in pathways:
        in_l1000 = sorted(genes & l1000_set)
        if len(in_l1000) >= 3:
            valid.append({
                "name": name,
                "direction": direction,
                "description": f"KEGG/Reactome pathway: {name}. "
                               f"{'Activated' if direction == 'UP' else 'Suppressed'} in MIRI.",
                "genes": in_l1000,
                "n_genes": len(in_l1000),
            })
            print(f"  {name} ({direction}): {len(in_l1000)} genes")
        else:
            print(f"  {name}: SKIPPED ({len(in_l1000)} genes < 3)")

    pw_set = {
        "name": "MIRI_KEGG_pathways",
        "source": "KEGG PATHWAY + Reactome + MIRI literature curation",
        "description": (
            f"Curated MIRI-relevant pathways with direction labels. "
            f"UP = pathway activated in MIRI (drug should DOWN-regulate). "
            f"DOWN = pathway suppressed in MIRI (drug should UP-regulate). "
            f"Genes intersected with L1000 978 landmark genes."
        ),
        "n_pathways_total": len(pathways),
        "n_pathways_valid": len(valid),
        "pathways": valid,
    }

    print(f"\nPathways: {len(valid)}/{len(pathways)} valid (≥3 L1000 genes)")
    return save_json(pw_set, "MIRI_pathways.json")


# ══════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════

def main():
    l1000_set = load_l1000_genes()

    print("\n═══ Generating individual dataset gene sets ═══\n")
    gse160516 = generate_gse160516(l1000_set)
    gse108940 = generate_gse108940(l1000_set)
    gse137293 = generate_gse137293(l1000_set)

    print("\n═══ Generating cross-species consensus ═══\n")
    consensus = generate_consensus([gse160516, gse108940, gse137293], l1000_set)

    print("\n═══ Generating pathway gene sets ═══\n")
    pathways = generate_pathway_genesets(l1000_set)

    print("\n" + "=" * 60)
    print("Summary:")
    print(f"  GSE160516 (mouse):      {gse160516['n_up_l1000']} UP + {gse160516['n_down_l1000']} DOWN")
    if gse108940:
        print(f"  GSE108940 (rat):        {gse108940['n_up_l1000']} UP + {gse108940['n_down_l1000']} DOWN")
    if gse137293:
        print(f"  GSE137293 (pig):        {gse137293['n_up_l1000']} UP + {gse137293['n_down_l1000']} DOWN")
    if consensus:
        print(f"  MIRI consensus (≥2/3):  {consensus['n_up_l1000']} UP + {consensus['n_down_l1000']} DOWN")
    print(f"  Pathways:               {pathways['n_pathways_valid']} valid")
    print("\n✓ All reference gene sets generated!")


if __name__ == "__main__":
    main()
