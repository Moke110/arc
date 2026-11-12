#!/usr/bin/env Rscript
#
# build_deg.R — full-transcriptome differential expression for GSE116250
# (human left-ventricle bulk RNA-seq; 50 HF = 37 DCM + 13 ICM, 14 NF).
#
# What is modelled
#   Primary contrast HF_vs_NF = (DCM + ICM)/2 - NF  (50 vs 14).
#   The two HF arms are additionally fitted as their own contrasts in the same
#   linear model: DCM_vs_NF (37 vs 14) and ICM_vs_NF (13 vs 14). Donor sex, age
#   and other clinical covariates are not modelled; the ticket asks for the
#   pooled HF contrast with the arms recorded separately, and the arms ARE
#   recorded in the same table through the `contrast` column.
#
# Method
#   limma on log2(RPKM + 1) (RPKM as deposited; 57,974 Ensembl annotations
#   collapsed to one row per HGNC symbol by max mean expression). Genes are
#   pre-filtered to those expressed (RPKM > 0) in >= 25% of samples, documented
#   in the JSON; no L1000 filtering. The 999999999 sentinel cells are NA (see
#   script/lib_gse116250.R) and limma fits each gene on its non-missing samples.
#   lmFit(~0+group) + eBayes(trend = TRUE, robust = TRUE).
#
# Outputs
#   local/deg/GSE116250_deg.csv   full-transcriptome DEG table, 3 contrasts
#                                 (gitignored)
#   local/deg/GSE116250_deg.json  small tracked summary (declared in
#                                 author_log.json local_additions[])
#
# Usage:
#   cd disease_data/GSE116250
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

source(file.path(SCRIPT_DIR, "lib_gse116250.R"))

PADJ_CUT <- 0.05
LFC_CUT <- 0.5

sha256_file <- function(path) {
  out <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = FALSE)
  sub("[[:space:]].*$", "", out)
}

cat("\n===== GSE116250 - human LV bulk RNA-seq, HF vs NF =====\n\n")
dat <- load_gse116250_expression(SOURCE_DIR)
cat(sprintf("Deposited RPKM matrix: %d annotations x %d samples; sentinel cells -> NA: %d\n",
            nrow(dat$expr), ncol(dat$expr), dat$n_sentinel))
cat(sprintf("  groups: %s\n",
            paste(names(table(dat$group)), as.integer(table(dat$group)),
                  sep = "=", collapse = ", ")))

expr <- filter_expressed(dat$expr, min_frac = 0.25)
cat(sprintf("Expression filter (RPKM>0 in >=25%% of samples): %d of %d genes kept\n",
            nrow(expr), nrow(dat$expr)))

group <- factor(dat$group, levels = c("NF", "DCM", "ICM"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

cat("limma fit: HF_vs_NF (primary), DCM_vs_NF, ICM_vs_NF ...\n")
fit <- lmFit(expr, design)
contr <- makeContrasts(
  HF_vs_NF = (DCM + ICM) / 2 - NF,
  DCM_vs_NF = DCM - NF,
  ICM_vs_NF = ICM - NF,
  levels = design
)
fit2 <- eBayes(contrasts.fit(fit, contr), trend = TRUE, robust = TRUE)

contrast_spec <- list(
  HF_vs_NF = c(n_case = sum(group %in% c("DCM", "ICM")), n_control = sum(group == "NF")),
  DCM_vs_NF = c(n_case = sum(group == "DCM"), n_control = sum(group == "NF")),
  ICM_vs_NF = c(n_case = sum(group == "ICM"), n_control = sum(group == "NF"))
)

blocks <- list()
for (cn in names(contrast_spec)) {
  tt <- topTable(fit2, coef = cn, number = Inf, adjust.method = "BH", sort.by = "P")
  blocks[[cn]] <- data.frame(
    gene_symbol = rownames(tt),
    contrast = cn,
    log2FC = tt$logFC,
    stat = tt$t,
    p_value = tt$P.Value,
    padj = tt$adj.P.Val,
    mean_expr = tt$AveExpr,
    n_case = contrast_spec[[cn]][["n_case"]],
    n_control = contrast_spec[[cn]][["n_control"]],
    stringsAsFactors = FALSE
  )
}
deg <- do.call(rbind, blocks)
rownames(deg) <- NULL

count_sig <- function(d) {
  c(n_up = sum(d$padj < PADJ_CUT & d$log2FC > LFC_CUT, na.rm = TRUE),
    n_down = sum(d$padj < PADJ_CUT & d$log2FC < -LFC_CUT, na.rm = TRUE))
}
primary <- blocks$HF_vs_NF
counts <- count_sig(primary)
cat(sprintf("DEG: %d genes tested; HF_vs_NF padj<%.2f & |log2FC|>%.1f: %d up, %d down\n",
            nrow(primary), PADJ_CUT, LFC_CUT, counts[["n_up"]], counts[["n_down"]]))

out_csv <- file.path(LOCAL_DEG_DIR, "GSE116250_deg.csv")
write.csv(deg, out_csv, row.names = FALSE, quote = FALSE)
cat(sprintf("wrote %s (%.1f MiB, %d contrast rows)\n", out_csv,
            file.info(out_csv)$size / 1024^2, nrow(deg)))

arm_counts <- lapply(blocks[-1], count_sig)
src_digests <- list(
  `GSE116250_rpkm.txt.gz` = sha256_file(file.path(SOURCE_DIR, "GSE116250_rpkm.txt.gz")),
  `lib_gse116250.R` = sha256_file(file.path(SCRIPT_DIR, "lib_gse116250.R"))
)

summary <- list(
  accession = "GSE116250",
  species = "Homo sapiens",
  platform = "Illumina HiSeq 2000 (GPL16791), hg19 / Ensembl GRCh37.82 Cufflinks RPKM",
  tissue_model = "human left-ventricular bulk tissue",
  contrast = "HF_vs_NF",
  contrast_description = paste0(
    "primary HF_vs_NF = (DCM + ICM)/2 - NF, 50 HF (37 DCM + 13 ICM) vs 14 NF; ",
    "DCM_vs_NF and ICM_vs_NF fitted in the same model and recorded as extra ",
    "`contrast` rows; positive log2FC = higher in disease"),
  arm_contrasts = list(
    DCM_vs_NF = list(n_case = 37L, n_control = 14L,
                     n_up = arm_counts$DCM_vs_NF[["n_up"]],
                     n_down = arm_counts$DCM_vs_NF[["n_down"]]),
    ICM_vs_NF = list(n_case = 13L, n_control = 14L,
                     n_up = arm_counts$ICM_vs_NF[["n_up"]],
                     n_down = arm_counts$ICM_vs_NF[["n_down"]])
  ),
  method = paste0(
    "limma on log2(RPKM+1); Ensembl annotations collapsed to HGNC symbols by max ",
    "mean expression; 999999999 sentinel cells -> NA; expression filter RPKM>0 in ",
    ">=25% of samples; lmFit(~0+group) with group in {NF,DCM,ICM} + ",
    "eBayes(trend=TRUE, robust=TRUE)"),
  covariates = paste0(
    "group only (three levels); sex, age and other clinical covariates not modelled",
    "; the DCM and ICM arms are reported as separate contrasts"),
  thresholds = list(padj = PADJ_CUT, abs_log2FC = LFC_CUT, p_adjust_method = "BH"),
  n_genes_tested = nrow(primary),
  n_up = counts[["n_up"]],
  n_down = counts[["n_down"]],
  n_case = 50L,
  n_control = 14L,
  group_counts = list(NF = 14L, DCM = 37L, ICM = 13L),
  normalization = "log2(RPKM + 1), RPKM as deposited",
  expression_filter = "expressed (RPKM > 0) in >= 25% of 64 samples",
  expression_matrix = "local/expr/GSE116250_log_expr.tsv.gz",
  samples = "local/expr/GSE116250_samples.csv",
  source_digests = src_digests,
  deg_csv = "local/deg/GSE116250_deg.csv",
  reproducibility = "deterministic; re-running build_deg.R reproduces the CSV sha256 exactly"
)
out_json <- file.path(LOCAL_DEG_DIR, "GSE116250_deg.json")
writeLines(jsonlite::toJSON(summary, pretty = TRUE, auto_unbox = TRUE, digits = 10),
           out_json)
cat(sprintf("wrote %s\n", out_json))
cat(sprintf("deg_csv_sha256: %s\n", sha256_file(out_csv)))
cat("\nGSE116250 DEG done\n")
