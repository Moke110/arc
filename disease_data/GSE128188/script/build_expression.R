#!/usr/bin/env Rscript
#
# build_expression.R — normalized gene-level log expression for GSE128188
# (human atrial fibrillation vs sinus rhythm; paired left/right atrial
# appendage bulk RNA-seq).
#
# Input (source/, ADR 0030 — never modified):
#   source/GSE128188_CountsV2_EdgeR.txt.gz — Gene x 20 integer-count columns
#   (5 SR donors + 5 AF donors, each with paired LA and RA libraries).
#
# Normalization (documented because the gene-axis step reads this matrix):
#   log2(CPM + 1) per GEO library, where CPM = count / (library size / 1e6)
#   and the library size is the column sum over every gene in the upstream
#   counts table (no additional filtering). This is the same library-level
#   matrix the series ships; the DESeq2 contrast in build_deg.R pools the
#   paired LA+RA libraries per donor, so the DEG file's n_case/n_control are
#   donors (5/5) while this matrix keeps all 20 GEO libraries.
#
# Outputs (local/, never source/):
#   local/expr/GSE128188_log_expr.tsv.gz  genes x 20 libraries, log2(CPM+1)
#   local/expr/GSE128188_samples.csv      sample, group, donor, tissue
#
# Usage (after build_deg.R, which is independent of this script):
#   cd disease_data/GSE128188
#   /home/chang/miniconda3/envs/biomni_e1/bin/Rscript script/build_expression.R

suppressPackageStartupMessages({
  library(digest)
})

SCRIPT_DIR <- dirname(normalizePath(
  sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1])
))
ITEM_DIR <- dirname(SCRIPT_DIR)
SOURCE_DIR <- file.path(ITEM_DIR, "source")
LOCAL_EXPR_DIR <- file.path(ITEM_DIR, "local", "expr")
dir.create(LOCAL_EXPR_DIR, recursive = TRUE, showWarnings = FALSE)

COUNTS_FILE <- file.path(SOURCE_DIR, "GSE128188_CountsV2_EdgeR.txt.gz")
EXPR_TSV <- file.path(LOCAL_EXPR_DIR, "GSE128188_log_expr.tsv.gz")
SAMPLES_CSV <- file.path(LOCAL_EXPR_DIR, "GSE128188_samples.csv")

# GEO series-matrix sample titles give the group of each donor (see
# build_deg.R; KEEP IN SYNC). Counts columns are `<donor>-<tissue>`.
SR_DONORS <- c("7", "12b", "20", "27", "30")
AF_DONORS <- c("5", "14", "25", "28", "31")

sha256 <- function(path) digest::digest(file = path, algo = "sha256")

cat("\n===== GSE128188 - log2(CPM+1) expression matrix =====\n\n")

if (!file.exists(COUNTS_FILE)) {
  stop("missing source/", basename(COUNTS_FILE),
       " — run 'bash script/fetch.sh' first")
}

raw <- read.delim(gzfile(COUNTS_FILE), check.names = FALSE,
                  stringsAsFactors = FALSE)
counts <- as.matrix(raw[, -1, drop = FALSE])
storage.mode(counts) <- "double"
rownames(counts) <- raw[[1]]

samples <- colnames(counts)
donor <- sub("-[^-]+$", "", samples)
tissue <- sub("^.*-", "", samples)
if (!setequal(unique(donor), c(SR_DONORS, AF_DONORS))) {
  stop("donor IDs in the counts table do not match the GEO series-matrix donors")
}
if (!setequal(tissue, c("LA", "RA"))) {
  stop("unexpected tissue labels in counts columns")
}
group <- ifelse(donor %in% SR_DONORS, "SR", "AF")
if (sum(group == "SR") != 10L || sum(group == "AF") != 10L) {
  stop("expected 10 SR and 10 AF libraries")
}

# --- log2(CPM + 1), library size = sum over all source genes ----------------
lib_size <- colSums(counts)
if (any(lib_size <= 0)) stop("empty library among the counts columns")
cpm <- sweep(counts, 2, lib_size / 1e6, "/")
log_expr <- log2(cpm + 1)

# --- Write the genes x samples matrix --------------------------------------
genes <- rownames(counts)
ord <- order(genes, method = "radix")
log_expr <- log_expr[ord, , drop = FALSE]
out <- data.frame(gene_symbol = genes[ord], log_expr, check.names = FALSE)
con <- gzfile(EXPR_TSV, "w")
write.table(out, con, sep = "\t", quote = FALSE, row.names = FALSE,
            col.names = TRUE)
close(con)

# --- Write the sample sheet ------------------------------------------------
sample_sheet <- data.frame(
  sample = samples,
  group = group,
  donor = donor,
  tissue = tissue,
  stringsAsFactors = FALSE
)
write.csv(sample_sheet, SAMPLES_CSV, row.names = FALSE)

cat(sprintf("log_expr: %d genes x %d samples; samples.csv: %d rows\n",
            nrow(out), ncol(out) - 1L, nrow(sample_sheet)))
cat(sprintf("groups: %s\n",
            paste(names(table(sample_sheet$group)),
                  table(sample_sheet$group), sep = "=", collapse = ", ")))
cat(sprintf("log_expr.tsv.gz sha256: %s\n", sha256(EXPR_TSV)))
cat(sprintf("samples.csv sha256: %s\n", sha256(SAMPLES_CSV)))
