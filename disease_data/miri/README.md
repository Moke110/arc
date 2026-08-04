# disease_data/miri — MIRI 疾病共享数据

MIRI（心肌缺血再灌注损伤）疾病特异性共享输入。原始 GEO 数据采集与各数据集（GSE）的 DEG 列表在此按 GSE 组织，可被未来的工作流复用。

## 结构

```
miri/
├── scripts/                  # 共享直系同源映射工具 (R)
│   ├── map_to_human.R        # 小鼠/大鼠/猪 → 人类 L1000 映射函数库
│   └── map_orthologs.R       # CLI: 大鼠/猪 DEG 结果 → 人类 L1000 列表
├── GSE160516/                # 小鼠 MIRI (Affymetrix Clariom S)，主要参考
│   ├── data/                 # 原始 GEO 数据 (limma DEG 结果)
│   └── deg/                  # per-GSE DEG 列表 + build_deg.R
├── GSE108940/                # 大鼠 MIRI (Agilent GE2-v5)
│   ├── data/
│   └── deg/
└── GSE137293/                # 猪 scar vs border (RNA-seq，无 sham 对照)
    ├── data/
    └── deg/
```

## 每个 GSE 目录

- `data/` — 原始 GEO 采集（DEG 结果、原始计数）
- `deg/` — per-GSE DEG 列表（映射到人类 L1000）+ 该 GSE 的构建脚本 `build_deg.R`

## 使用方法

```bash
conda activate arc
# 重建某个 GSE 的 DEG 列表
cd disease_data/miri/GSE160516 && Rscript deg/build_deg.R

# 映射大鼠/猪 DEG 到人类 L1000
cd disease_data/miri && Rscript scripts/map_orthologs.R --deg ... --out ... --species rat|pig
```

## 消费方

工作流阶段 07（`miri_reposition/xpert-fda-v1/scripts/07_build_miri_deg_reference.py`）消费各 GSE 的 `deg/` 列表，生成集成 MIRI 参考基因集。

## 说明

- 原单文件 3-GSE R 分析脚本已拆分为各 GSE 的 `build_deg.R`。
- 猪（GSE137293）为 scar vs border 比较（无 sham 对照），使用需谨慎。
- 直系同源映射使用 conda env `arc` 中的 R + homologene。
