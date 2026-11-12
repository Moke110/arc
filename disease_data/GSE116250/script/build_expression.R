#!/usr/bin/env Rscript
#
# build_expression.R — normalized, gene-level, log-scale expression matrix for
# GSE116250 (human left-ventricle bulk RNA-seq; 50 HF / 14 NF).
#
# This is the only interface the Stage 2 gene-axis step reads from this item.
# The matrix is log2(RPKM + 1) of the deposited Cufflinks RPKM values
# (GEO sample_data_processing), with the 999999999 sentinel cells set to NA.
# RPKM is already library-size normalized, so no further normalization is
# applied; mitochondrial genes were removed by the submitters (57,974
# annotations). Annotations are collapsed to one row per HGNC symbol by max mean
# expression. No expression filtering here (full measured transcriptome); the
# DEG step applies its documented filter.
#
# Outputs
#   local/expr/GSE116250_log_expr.tsv.gz  genes x samples (+ header line)
#   local/expr/GSE116250_samples.csv      sample,group (NF/DCM/ICM)
#
# Usage:
#   cd disease_data/GSE116250
#   ulimit -v 6G
#   OMP_NUM_THREADS=2 /home/chang/miniconda3/envs/biomni_e1/bin/Rscript script/build_expression.R

suppressPackageStartupMessages({
  library(data.table)
})

SCRIPT_DIR <- dirname(normalizePath(
  sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1])
))
ITEM_DIR <- dirname(SCRIPT_DIR)
SOURCE_DIR <- file.path(ITEM_DIR, "source")
LOCAL_EXPR_DIR <- file.path(ITEM_DIR, "local", "expr")
dir.create(LOCAL_EXPR_DIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(SCRIPT_DIR, "lib_gse116250.R"))

cat("\n===== GSE116250 - expression matrix =====\n\n")
dat <- load_gse116250_expression(SOURCE_DIR)
expr <- round(dat$expr, 6)
cat(sprintf("Gene-level log2 matrix: %d symbols x %d samples; NA cells: %d\n",
            nrow(expr), ncol(expr), sum(is.na(expr))))

out_tsv <- file.path(LOCAL_EXPR_DIR, "GSE116250_log_expr.tsv.gz")
con <- gzfile(out_tsv, "w")
writeLines(sprintf(
  "# genes=%d samples=%d normalization=log2(RPKM+1)_as_deposited,sentinel_999999999_as_NA contrast=HF_vs_NF",
  nrow(expr), ncol(expr)), con)
write.table(data.frame(gene_symbol = rownames(expr), expr, check.names = FALSE),
            con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)
cat(sprintf("wrote %s (%.1f MiB)\n", out_tsv, file.info(out_tsv)$size / 1024^2))

samples <- data.frame(sample = dat$sample, group = dat$group, stringsAsFactors = FALSE)
out_csv <- file.path(LOCAL_EXPR_DIR, "GSE116250_samples.csv")
write.csv(samples, out_csv, row.names = FALSE, quote = FALSE)
cat(sprintf("wrote %s (%d rows: %s)\n", out_csv, nrow(samples),
            paste(names(table(samples$group)), as.integer(table(samples$group)),
                  sep = "=", collapse = ", ")))
cat("\nGSE116250 expression done\n")
