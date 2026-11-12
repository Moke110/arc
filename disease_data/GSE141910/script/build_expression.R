#!/usr/bin/env Rscript
#
# build_expression.R — normalized, gene-level, log-scale expression matrix for
# GSE141910 (human left-ventricle bulk RNA-seq, MAGNet consortium).
#
# This is the only interface the Stage 2 gene-axis step reads from this item.
# The matrix is the deposited GSE141910 per-sample voom-transformed abundance
# (already log2-scale and library-size normalized by the submitters; GEO's
# sample_data_processing says "Raw gene counts were transformed using SVA and
# the 'LIMMA-Voom' method"), collapsed to one row per HGNC symbol by max mean
# expression. It is NOT re-normalized here: the raw counts are not deposited, so
# any additional normalization would be a transformation of an already
# log2-scale matrix, not a count-based normalization. Values are rounded to 6
# decimal places for file size; no filtering is applied (full measured
# transcriptome, 20,781 Ensembl genes -> symbol-level rows).
#
# Outputs
#   local/expr/GSE141910_log_expr.tsv.gz  genes x samples (+ header line)
#   local/expr/GSE141910_samples.csv      sample,group
#
# Usage:
#   cd disease_data/GSE141910
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

source(file.path(SCRIPT_DIR, "lib_gse141910.R"))

cat("\n===== GSE141910 - expression matrix =====\n\n")
dat <- load_gse141910_expression(SOURCE_DIR)
expr <- dat$expr
cat(sprintf("Deposited matrix: %d genes x %d samples\n", nrow(expr), ncol(expr)))

symbols <- map_ensg_to_symbol(rownames(expr))
expr <- collapse_to_symbols(expr, symbols)
expr <- round(expr, 6)
cat(sprintf("Gene-level log2 matrix: %d symbols x %d samples\n", nrow(expr), ncol(expr)))

out_tsv <- file.path(LOCAL_EXPR_DIR, "GSE141910_log_expr.tsv.gz")
con <- gzfile(out_tsv, "w")
writeLines(sprintf(
  "# genes=%d samples=%d normalization=voom_log2_abundance_as_deposited(SVA_preprocessed_by_submitters) contrast=HF_vs_NF",
  nrow(expr), ncol(expr)), con)
write.table(data.frame(gene_symbol = rownames(expr), expr, check.names = FALSE),
            con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)
cat(sprintf("wrote %s (%.1f MiB)\n", out_tsv, file.info(out_tsv)$size / 1024^2))

samples <- data.frame(sample = dat$sample, group = dat$group, stringsAsFactors = FALSE)
out_csv <- file.path(LOCAL_EXPR_DIR, "GSE141910_samples.csv")
write.csv(samples, out_csv, row.names = FALSE, quote = FALSE)
cat(sprintf("wrote %s (%d rows: %s)\n", out_csv, nrow(samples),
            paste(names(table(samples$group)), as.integer(table(samples$group)),
                  sep = "=", collapse = ", ")))
cat("\nGSE141910 expression done\n")
