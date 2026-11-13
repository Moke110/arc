#!/usr/bin/env Rscript
# Pooled class-vs-vehicle edgeR contrasts from the DToxS GSE217421 read counts.
#
# The deposit ships its own differential-comparison pipeline (Program-Source-
# Codes.tar.gz). This script runs that *pinned* code on a pooling design: for one
# mechanism class and one (cell line, plate) stratum, the class members' wells
# are pooled into one Disease group ("CLASS") and compared with the same
# stratum's DMSO control wells ("CTRL"), using the deposit's own edgeR settings
# (dispersion "trended", calcNormFactors TMM, remove_small filter, exactTest).
#
# The shipped Calc tables only exist per drug, so the pooling refit is this
# step's declared analysis; the method and the inputs are recorded by the
# Python driver.
#
# Usage:
#   Rscript cardiotox_edgeR.R COUNTS_TSV SAMPLES_TSV PIPELINE_LIBRARY OUT_TSV
#
# COUNTS_TSV: genes x wells (first column gene symbol, tab separated)
# SAMPLES_TSV: well, stratum, group (group in CTRL/CLASS)
# OUT_TSV: stratum, gene_symbol, logFC, logCPM, PValue, FDR, n_case, n_control

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: cardiotox_edgeR.R COUNTS_TSV SAMPLES_TSV PIPELINE_LIBRARY OUT_TSV")
}
counts_path <- args[1]
samples_path <- args[2]
func_dir <- args[3]
out_path <- args[4]

suppressPackageStartupMessages(require(edgeR))
# edgeR >= 4 removed decideTestsDGE, which the pinned DToxS pipeline calls; the
# two-group exactTest decision is decideTests(method = "separate") with the same
# BH adjustment.
if (!exists("decideTestsDGE")) {
  decideTestsDGE <- function(object, adjust.method = "BH", p.value = 0.05, ...)
    decideTests(object, method = "separate", adjust.method = adjust.method,
                p.value = p.value, ...)
}
source(file.path(func_dir, "compareGeneExpDiffSubset.R"), local = TRUE)

counts <- read.delim(counts_path, check.names = FALSE, stringsAsFactors = FALSE)
genes <- counts[[1]]
matrix <- as.matrix(counts[, -1, drop = FALSE])
storage.mode(matrix) <- "double"
rownames(matrix) <- genes
samples <- read.delim(samples_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!identical(colnames(matrix), samples$well)) {
  stop("the counts columns and the samples table are not in the same order")
}
if (!setequal(unique(samples$group), c("CTRL", "CLASS"))) {
  stop("the samples table needs both CTRL and CLASS groups")
}

results <- list()
for (stratum in sort(unique(samples$stratum))) {
  keep <- samples$stratum == stratum
  groups <- samples$group[keep]
  if (length(unique(groups)) != 2L) {
    stop("stratum ", stratum, " lacks one arm")
  }
  if (min(table(groups)) < 2L) {
    stop("stratum ", stratum, " needs at least two wells per arm")
  }
  subset <- matrix[, keep, drop = FALSE]
  response <- compareGeneExpDiffSubset(
    read_counts = subset,
    group_names = groups,
    group_pair = c("CTRL", "CLASS"),
    subset_name = stratum,
    dispersion = "trended",
    min_samples = 2,
    min_samples_req = "all",
    fdr = 0.1,
    deg_only = FALSE,
    normalize = TRUE,
    norm_base = 1e+6,
    remove_small = TRUE,
    verbose = FALSE,
    func_dir = func_dir)
  if (is.null(response)) {
    stop("the pinned edgeR pipeline returned no result for stratum ", stratum)
  }
  results[[length(results) + 1L]] <- data.frame(
    stratum = stratum,
    gene_symbol = rownames(response),
    logFC = response[, "logFC"],
    logCPM = response[, "logCPM"],
    PValue = response[, "PValue"],
    FDR = response[, "FDR"],
    n_case = sum(samples$group[keep] == "CLASS"),
    n_control = sum(samples$group[keep] == "CTRL"),
    stringsAsFactors = FALSE)
}
output <- do.call(rbind, results)
format_number <- function(values) {
  out <- rep("", length(values))
  finite <- is.finite(values)
  out[finite] <- sprintf("%.10g", values[finite])
  out
}
output$logFC <- format_number(output$logFC)
output$logCPM <- format_number(output$logCPM)
output$PValue <- format_number(output$PValue)
output$FDR <- format_number(output$FDR)
write.table(output, out_path, sep = "\t", quote = FALSE, row.names = FALSE,
            na = "")
cat("edgeR pooled class contrasts:", length(results), "strata,", nrow(output),
    "rows\n")
