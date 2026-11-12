#!/usr/bin/env Rscript
#
# build_expression.R — normalized gene-level log expression for GSE282618
# (human hypertensive heart disease [HHD] vs control; post-mortem LV+RV bulk
# RNA-seq). ** EXPLORATORY ITEM ** — n = 5 subjects.
#
# Input (source/, ADR 0030 — never modified):
#   source/GSE282618_counts.txt.gz — 60,646 Ensembl gene IDs x 10 integer-count
#   libraries (5 subjects x LV/RV; 3 hypertensive: NDRI1, NDRI3, NDRI5;
#   2 control: NDRI4, NDRI7).
#
# Normalization (documented because the gene-axis step reads this matrix):
#   log2(CPM + 1) per GEO library, where CPM = count / (library size / 1e6)
#   and the library size is the column sum over every Ensembl gene in the
#   upstream counts table (before any symbol mapping). This matrix keeps all 10
#   GEO libraries; build_deg.R pools the paired LV+RV libraries per subject for
#   the 3-vs-2 contrast.
#
# Gene symbols: identical Ensembl -> symbol mapping and duplicate-symbol
# collapse as build_deg.R (KEEP IN SYNC) so both files share one gene set:
# org.Hs.eg.db ENSEMBL -> SYMBOL, unmapped Ensembl IDs dropped, duplicate
# symbols collapsed to the Ensembl ID with the highest mean CPM across the 10
# libraries.
#
# Sample IDs: the GEO series-matrix titles are "NDRI4 LV" (with a space); the
# matrix header cannot carry that safely, so the shipped sample ID is
# `NDRI4_LV` and the verbatim GEO title is kept in the sample sheet's
# `geo_title` column.
#
# Outputs (local/, never source/):
#   local/expr/GSE282618_log_expr.tsv.gz  genes x 10 libraries, log2(CPM+1);
#                                         first line is a metadata comment
#                                         '# genes=... samples=...
#                                         normalization=... contrast=...'
#                                         (same convention as the sibling
#                                         gene-axis items)
#   local/expr/GSE282618_samples.csv      sample, group, subject, site, geo_title
#
# Usage (after build_deg.R, which is independent of this script):
#   cd disease_data/GSE282618
#   /home/chang/miniconda3/envs/biomni_e1/bin/Rscript script/build_expression.R

suppressPackageStartupMessages({
  library(digest)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

SCRIPT_DIR <- dirname(normalizePath(
  sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1])
))
ITEM_DIR <- dirname(SCRIPT_DIR)
SOURCE_DIR <- file.path(ITEM_DIR, "source")
LOCAL_EXPR_DIR <- file.path(ITEM_DIR, "local", "expr")
dir.create(LOCAL_EXPR_DIR, recursive = TRUE, showWarnings = FALSE)

COUNTS_FILE <- file.path(SOURCE_DIR, "GSE282618_counts.txt.gz")
EXPR_TSV <- file.path(LOCAL_EXPR_DIR, "GSE282618_log_expr.tsv.gz")
SAMPLES_CSV <- file.path(LOCAL_EXPR_DIR, "GSE282618_samples.csv")

# GEO series matrix: NDRI1/NDRI3/NDRI5 hypertensive, NDRI4/NDRI7 control.
HHD_SUBJECTS <- c("NDRI1", "NDRI3", "NDRI5")
CONTROL_SUBJECTS <- c("NDRI4", "NDRI7")

sha256 <- function(path) digest::digest(file = path, algo = "sha256")

cat("\n===== GSE282618 - log2(CPM+1) expression matrix (EXPLORATORY) =====\n\n")

if (!file.exists(COUNTS_FILE)) {
  stop("missing source/", basename(COUNTS_FILE),
       " — run 'bash script/fetch.sh' first")
}

raw <- read.delim(gzfile(COUNTS_FILE), check.names = FALSE,
                  stringsAsFactors = FALSE)
counts <- as.matrix(raw[, -1, drop = FALSE])
storage.mode(counts) <- "double"
ensg <- raw[[1]]
stopifnot(!anyDuplicated(ensg), all(grepl("^ENSG", ensg)))
rownames(counts) <- ensg

sample <- colnames(counts)
subject <- sub(" .*$", "", sample)
site <- sub("^.* ", "", sample)
if (!setequal(unique(subject), c(HHD_SUBJECTS, CONTROL_SUBJECTS))) {
  stop("subject IDs in the counts table do not match the GEO series-matrix subjects")
}
if (!setequal(site, c("LV", "RV"))) {
  stop("unexpected ventricular-site labels in counts columns")
}
group <- ifelse(subject %in% HHD_SUBJECTS, "HHD", "control")
if (sum(group == "HHD") != 6L || sum(group == "control") != 4L) {
  stop("expected 6 HHD and 4 control libraries")
}

# --- Ensembl ID -> human symbol, collapse duplicate symbols ---------------
# KEEP IN SYNC with script/build_deg.R (same annotation and collapse rule).
symbol <- suppressWarnings(mapIds(
  org.Hs.eg.db, keys = ensg, column = "SYMBOL", keytype = "ENSEMBL",
  multiVals = "first"
))
mapped <- !is.na(symbol) & nzchar(symbol)
lib_size <- colSums(counts)
cpm <- sweep(counts, 2, lib_size / 1e6, "/")
mean_cpm <- rowMeans(cpm)
sel <- data.frame(
  ensg = ensg[mapped],
  symbol = as.character(symbol[mapped]),
  mean_cpm = mean_cpm[mapped],
  stringsAsFactors = FALSE
)
sel <- sel[order(sel$symbol, -sel$mean_cpm, sel$ensg, method = "radix"), ]
n_collapsed <- sum(duplicated(sel$symbol))
sel <- sel[!duplicated(sel$symbol), ]

counts <- counts[sel$ensg, , drop = FALSE]
rownames(counts) <- sel$symbol

# --- log2(CPM + 1), library size = sum over all upstream genes --------------
# lib_size was computed above on the full upstream matrix (before mapping), so
# subsetting the gene rows does not change the library size.
cpm <- sweep(counts, 2, lib_size / 1e6, "/")
log_expr <- log2(cpm + 1)
# The GEO titles are "NDRI4 LV" (space); the matrix header and the sample-sheet
# `sample` column share the same space-free ID, while `geo_title` keeps the
# verbatim label.
colnames(log_expr) <- gsub(" ", "_", sample)

genes <- rownames(counts)
ord <- order(genes, method = "radix")
log_expr <- log_expr[ord, , drop = FALSE]
out <- data.frame(gene_symbol = genes[ord], log_expr, check.names = FALSE)
# Leading metadata comment line, the same convention as the sibling gene-axis
# items: '# genes=<n> samples=<n> normalization=<desc> contrast=<contrast>'.
# It is written through the gzip connection so a re-run reproduces the file
# byte-for-byte.
con <- gzfile(EXPR_TSV, "w")
writeLines(sprintf(
  paste0("# genes=%d samples=%d normalization=",
         "log2(CPM+1)_library_sizes_from_all_source_genes contrast=HHD_vs_control"),
  nrow(out), ncol(out) - 1L), con)
write.table(out, con, sep = "\t", quote = FALSE, row.names = FALSE,
            col.names = TRUE)
close(con)

# --- Sample sheet (clean IDs + verbatim GEO titles) -----------------------
sample_sheet <- data.frame(
  sample = gsub(" ", "_", sample),
  group = group,
  subject = subject,
  site = site,
  geo_title = sample,
  stringsAsFactors = FALSE
)
write.csv(sample_sheet, SAMPLES_CSV, row.names = FALSE)

cat(sprintf("log_expr: %d genes x %d samples; samples.csv: %d rows\n",
            nrow(out), ncol(out) - 1L, nrow(sample_sheet)))
cat(sprintf("groups: %s\n",
            paste(names(table(sample_sheet$group)),
                  table(sample_sheet$group), sep = "=", collapse = ", ")))
cat(sprintf("duplicate symbols collapsed: %d\n", n_collapsed))
cat(sprintf("log_expr.tsv.gz sha256: %s\n", sha256(EXPR_TSV)))
cat(sprintf("samples.csv sha256: %s\n", sha256(SAMPLES_CSV)))
