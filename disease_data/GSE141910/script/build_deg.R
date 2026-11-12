#!/usr/bin/env Rscript
#
# build_deg.R — full-transcriptome HF-vs-NF differential expression for
# GSE141910 (human left-ventricle bulk RNA-seq, MAGNet consortium).
#
# What is modelled
#   Contrast: HF (DCM+HCM+PPCM, n=200) vs NF (non-failing donor hearts, n=166).
#   Only the case/control group enters the design. Etiology (166 DCM / 28 HCM /
#   6 PPCM among the HF samples), sex, age and race are NOT modelled as
#   covariates: the ticket asks for the pooled HF-vs-NF contrast, and the
#   per-sample etiology is recorded in script/lib_gse141910.R only. The
#   deposited expression values were already SVA-processed by the submitters
#   (GEO `Sample_data_processing`: "Raw gene counts were transformed using SVA
#   and the 'LIMMA-Voom' method"); the surrogate variables are not deposited, so
#   that preprocessing is inherited, not re-modelled.
#
# Method
#   limma on the deposited log2-scale voom-transformed abundance
#   (20,781 Ensembl genes, as shipped in GSE141910_RAW.tar). Raw counts were
#   never deposited, so DESeq2 / limma-voom weighting cannot be re-run from
#   source/; see author_log.json "notes". Genes are collapsed to HGNC symbols
#   (max mean expression per symbol); genes whose Ensembl v67 ID is retired in
#   current org.Hs.eg.db keep the Ensembl ID as `gene_symbol` so the table stays
#   full-transcriptome. eBayes(trend = TRUE, robust = TRUE).
#
# Outputs
#   local/deg/GSE141910_deg.csv   full-transcriptome DEG table (gitignored)
#   local/deg/GSE141910_deg.json  small tracked summary (declared in
#                                 author_log.json local_additions[])
#
# Usage:
#   cd disease_data/GSE141910
#   ulimit -v 6G
#   OMP_NUM_THREADS=2 /home/chang/miniconda3/envs/biomni_e1/bin/Rscript script/build_deg.R

suppressPackageStartupMessages({
  library(limma)
  library(jsonlite)
})

SCRIPT_DIR <- dirname(normalizePath(
  sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1])
))
ITEM_DIR <- dirname(SCRIPT_DIR)
SOURCE_DIR <- file.path(ITEM_DIR, "source")
LOCAL_DEG_DIR <- file.path(ITEM_DIR, "local", "deg")
dir.create(LOCAL_DEG_DIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(SCRIPT_DIR, "lib_gse141910.R"))

PADJ_CUT <- 0.05
LFC_CUT <- 0.5

sha256_file <- function(path) {
  out <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = FALSE)
  sub("[[:space:]].*$", "", out)
}

cat("\n===== GSE141910 - human LV bulk RNA-seq, HF vs NF =====\n\n")

dat <- load_gse141910_expression(SOURCE_DIR)
expr <- dat$expr
group <- dat$group
etiology <- dat$etiology
cat(sprintf("Deposited matrix: %d genes x %d samples (log2 voom abundance)\n",
            nrow(expr), ncol(expr)))
cat(sprintf("  HF: %d, NF: %d\n", sum(group == "HF"), sum(group == "NF")))
cat(sprintf("  HF etiology: %s\n",
            paste(names(table(etiology[group == "HF"])),
                  as.integer(table(etiology[group == "HF"])), sep = "=", collapse = ", ")))

cat("Mapping Ensembl v67 IDs to HGNC symbols (org.Hs.eg.db)...\n")
symbols <- map_ensg_to_symbol(rownames(expr))
n_fallback <- sum(symbols == rownames(expr))
cat(sprintf("  mapped: %d, Ensembl-ID fallback: %d, duplicated symbols collapsed: %d\n",
            length(symbols) - n_fallback, n_fallback,
            sum(duplicated(symbols))))
expr <- collapse_to_symbols(expr, symbols)
cat(sprintf("Gene-level matrix: %d symbols x %d samples\n", nrow(expr), ncol(expr)))

cat("limma fit: HF vs NF ...\n")
design <- model.matrix(~ 0 + factor(group, levels = c("NF", "HF")))
colnames(design) <- c("NF", "HF")
fit <- lmFit(expr, design)
fit2 <- eBayes(contrasts.fit(fit, makeContrasts(HF_vs_NF = HF - NF, levels = design)),
               trend = TRUE, robust = TRUE)
tt <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH", sort.by = "P")

deg <- data.frame(
  gene_symbol = rownames(tt),
  contrast = "HF_vs_NF",
  log2FC = tt$logFC,
  stat = tt$t,
  p_value = tt$P.Value,
  padj = tt$adj.P.Val,
  mean_expr = tt$AveExpr,
  n_case = sum(group == "HF"),
  n_control = sum(group == "NF"),
  stringsAsFactors = FALSE
)
deg <- deg[order(deg$padj, deg$p_value, -abs(deg$log2FC)), ]
rownames(deg) <- NULL

n_up <- sum(deg$padj < PADJ_CUT & deg$log2FC > LFC_CUT, na.rm = TRUE)
n_down <- sum(deg$padj < PADJ_CUT & deg$log2FC < -LFC_CUT, na.rm = TRUE)
cat(sprintf("DEG: %d genes tested; padj<%.2f & |log2FC|>%.1f: %d up, %d down\n",
            nrow(deg), PADJ_CUT, LFC_CUT, n_up, n_down))

out_csv <- file.path(LOCAL_DEG_DIR, "GSE141910_deg.csv")
write.csv(deg, out_csv, row.names = FALSE, quote = FALSE)
cat(sprintf("wrote %s (%.1f MiB)\n", out_csv, file.info(out_csv)$size / 1024^2))

src_digests <- list(
  `GSE141910_RAW.tar` = sha256_file(file.path(SOURCE_DIR, "GSE141910_RAW.tar")),
  `filelist.txt` = sha256_file(file.path(SOURCE_DIR, "filelist.txt")),
  `lib_gse141910.R` = sha256_file(file.path(SCRIPT_DIR, "lib_gse141910.R"))
)

summary <- list(
  accession = "GSE141910",
  species = "Homo sapiens",
  platform = "Illumina HiSeq 2500 (GPL16791), hg19 / Ensembl v67 gene counts",
  tissue_model = "human left-ventricular free-wall bulk tissue (MAGNet consortium)",
  contrast = "HF_vs_NF",
  contrast_description = paste0(
    "heart failure (DCM 166, HCM 28, PPCM 6; n=200) vs non-failing donor hearts ",
    "(n=166); positive log2FC = higher in HF"),
  method = paste0(
    "limma on the deposited log2-scale voom-transformed abundance ",
    "(GEO sample_data_processing: raw counts transformed using SVA and LIMMA-Voom); ",
    "HGNC symbols collapsed by max mean expression, retired Ensembl v67 IDs kept as ",
    "Ensembl-ID fallback; lmFit(~0+group) + eBayes(trend=TRUE, robust=TRUE)"),
  covariates = paste0(
    "group only. Etiology (DCM/HCM/PPCM) recorded (see lib_gse141910.R) but NOT modelled; ",
    "sex, age, race and submitter SVA surrogate variables are not modelled"),
  thresholds = list(padj = PADJ_CUT, abs_log2FC = LFC_CUT, p_adjust_method = "BH"),
  n_genes_tested = nrow(deg),
  n_up = n_up,
  n_down = n_down,
  n_case = sum(group == "HF"),
  n_control = sum(group == "NF"),
  case_etiology = as.list(table(etiology[group == "HF"])),
  expression_matrix = "local/expr/GSE141910_log_expr.tsv.gz",
  samples = "local/expr/GSE141910_samples.csv",
  source_digests = src_digests,
  deg_csv = "local/deg/GSE141910_deg.csv",
  reproducibility = "deterministic; re-running build_deg.R reproduces the CSV sha256 exactly"
)
out_json <- file.path(LOCAL_DEG_DIR, "GSE141910_deg.json")
writeLines(jsonlite::toJSON(summary, pretty = TRUE, auto_unbox = TRUE, digits = 10),
           out_json)
cat(sprintf("wrote %s\n", out_json))
cat(sprintf("deg_csv_sha256: %s\n", sha256_file(out_csv)))
cat("\nGSE141910 DEG done\n")
