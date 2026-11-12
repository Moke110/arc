#!/usr/bin/env Rscript
#
# build_expression.R — normalized log expression matrix for GSE243674
# (human hiPSC-CM chemotherapy cardiotoxicity, 6 healthy female donors).
#
# Input (source/, untouched upstream download, ADR 0030):
#   source/GSE243674_Counts_RNA_ERMatthews.txt.gz
#   ENTREZID x 72 wells, featureCounts gene-level counts on hg38.
#
# Outputs (declared in author_log.json; both regenerable and gitignored):
#   local/expr/GSE243674_log_expr.tsv.gz
#     genes x 72 wells, `gene_symbol` + one column per well (the source's own
#     `<DRUG>.<donor>.<time>` labels). Normalization: log2(CP10K + 1), i.e.
#     count / library_size * 1e4 then log2(x + 1), rounded to 4 decimals; the
#     library size is the column sum over the mapped genes retained in this
#     matrix. A leading comment line records genes=/samples=/normalization=.
#   local/expr/GSE243674_samples.csv
#     sample,drug,time_h,donor,group — one row per well, in matrix column order.
#     `drug` is the lower-case source token (dnr/dox/epi/mtx/trz/veh); the JSON
#     summary in local/deg/ reports the expanded drug names and classes.
#     `group` is treated/vehicle.
#
# Gene level: ENTREZID -> HGNC symbol with the pinned
# `pert_data/Annotation/source/hgnc_complete_set.txt` table, identical to
# script/build_deg.R (keep the two mapping blocks in sync; both use the same
# deterministic one-to-many rule and the same duplicate-symbol collapse). The
# resulting gene set is the same 26,867 symbols the DEG table is keyed by.
#
# Usage:
#   cd disease_data/GSE243674
#   Rscript script/build_expression.R
#
# Pinned R: /home/chang/miniconda3/envs/biomni_e1/bin/Rscript.
# Cap the run: ulimit -v 6G; OMP_NUM_THREADS=2. Peak RSS is recorded in
# author_log.json / the ticket record.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

options(stringsAsFactors = FALSE, datatable.verbose = FALSE)

# --- Locate the item from the script path (paths do not encode depth) --------
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("run this script as: Rscript script/build_expression.R")
SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", script_arg)))
find_repo_root <- function(start_dir) {
  dir <- normalizePath(start_dir, mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(dir, ".git"))) return(dir)
    parent <- dirname(dir)
    if (identical(parent, dir)) {
      stop("cannot find the repository root (no .git directory above ", start_dir, ")")
    }
    dir <- parent
  }
}
REPO_ROOT <- find_repo_root(SCRIPT_DIR)
ITEM_DIR <- dirname(SCRIPT_DIR)                       # disease_data/GSE243674
SOURCE_DIR <- file.path(ITEM_DIR, "source")
LOCAL_EXPR_DIR <- file.path(ITEM_DIR, "local", "expr")
dir.create(LOCAL_EXPR_DIR, recursive = TRUE, showWarnings = FALSE)

SOURCE_FILE <- file.path(SOURCE_DIR, "GSE243674_Counts_RNA_ERMatthews.txt.gz")
HGNC_TABLE <- Sys.getenv(
  "ARC_HGNC_TABLE",
  unset = file.path(REPO_ROOT, "pert_data", "Annotation", "source", "hgnc_complete_set.txt")
)
OUT_TSV <- file.path(LOCAL_EXPR_DIR, "GSE243674_log_expr.tsv.gz")
OUT_SAMPLES <- file.path(LOCAL_EXPR_DIR, "GSE243674_samples.csv")

if (!file.exists(SOURCE_FILE)) {
  stop("missing ", SOURCE_FILE, " — run `bash script/fetch.sh` first")
}
if (!file.exists(HGNC_TABLE)) {
  stop("pinned HGNC table not found: ", HGNC_TABLE,
       "\nIt lives in the gitignored pert_data/Annotation/source/ tree; in a",
       " checkout without it, point ARC_HGNC_TABLE at a materialised copy.")
}
sha256_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)

# --- Read the count matrix (same validation as script/build_deg.R) ----------
cat("GSE243674: reading", basename(SOURCE_FILE), "\n")
counts_dt <- fread(cmd = paste("gzip -cd", shQuote(SOURCE_FILE)), sep = " ", header = TRUE)
if (ncol(counts_dt) != 73L || names(counts_dt)[1] != "ENTREZID") {
  stop("unexpected matrix shape: ", nrow(counts_dt), " x ", ncol(counts_dt),
       " (expected ENTREZID + 72 wells)")
}
wells <- names(counts_dt)[-1]
counts <- as.matrix(counts_dt[, -1, with = FALSE])
rownames(counts) <- as.character(counts_dt[[1]])
gene_ids <- rownames(counts)
if (anyDuplicated(gene_ids)) stop("duplicated ENTREZID rows in the source matrix")

parsed <- regmatches(wells, regexec("^([A-Z]+)\\.([0-9]+)\\.([0-9]+)h$", wells))
if (any(lengths(parsed) != 4L)) {
  stop("cannot parse well labels: ", paste(wells[lengths(parsed) != 4L], collapse = ", "))
}
well_dt <- data.table(
  sample = wells,
  token = vapply(parsed, `[`, character(1), 2L),
  donor = vapply(parsed, `[`, character(1), 3L),
  time_h = as.integer(vapply(parsed, `[`, character(1), 4L))
)
if (!setequal(unique(well_dt$token), c("DNR", "DOX", "EPI", "MTX", "TRZ", "VEH"))) {
  stop("unexpected drug tokens in the source header: ",
       paste(sort(unique(well_dt$token)), collapse = ", "))
}
well_dt[, group := ifelse(token == "VEH", "vehicle", "treated")]
well_dt[, drug := tolower(token)]
if (nrow(well_dt) != 72L || length(unique(well_dt$donor)) != 6L ||
    !setequal(unique(well_dt$time_h), c(3L, 24L)) ||
    sum(well_dt$group == "vehicle") != 12L) {
  stop("unbalanced design in the source header: ", nrow(well_dt), " wells")
}

# --- ENTREZID -> HGNC symbol (identical rule to script/build_deg.R) ---------
hg <- fread(HGNC_TABLE, sep = "\t", header = TRUE, quote = "",
            select = c("hgnc_id", "symbol", "status", "locus_group", "entrez_id"))
hg <- hg[status == "Approved" & !is.na(entrez_id) & entrez_id != ""]
hg[, entrez_id := as.character(entrez_id)]
hg[, locus_rank := match(locus_group,
                         c("protein-coding gene", "non-coding RNA", "pseudogene"),
                         nomatch = 99L)]
hg[, hgnc_num := as.integer(sub("^HGNC:", "", hgnc_id))]
setorder(hg, entrez_id, locus_rank, hgnc_num)
map <- hg[!duplicated(entrez_id), .(entrez_id, symbol)]

hit <- map$symbol[match(gene_ids, map$entrez_id)]
mapped <- !is.na(hit)
counts <- counts[mapped, , drop = FALSE]
symbols <- hit[mapped]
if (any(duplicated(symbols))) {
  counts <- rowsum(counts, group = symbols)
  symbols <- rownames(counts)
}
rownames(counts) <- symbols
n_genes <- length(symbols)

# --- log2(CP10K + 1) --------------------------------------------------------
lib_size <- colSums(counts)
if (any(lib_size <= 0)) stop("empty library in ", paste(wells[lib_size <= 0], collapse = ", "))
cpm <- sweep(counts, 2, lib_size, "/") * 1e4
log_expr <- log2(cpm + 1)

header <- sprintf(
  "# genes=%d samples=%d normalization=log2(CP10K+1) contrast=none source=%s",
  n_genes, ncol(log_expr), basename(SOURCE_FILE)
)
con <- gzfile(OUT_TSV, "w")
on.exit(close(con), add = TRUE)
writeLines(header, con)
writeLines(paste(c("gene_symbol", colnames(log_expr)), collapse = "\t"), con)
chunk <- 2000L
for (start in seq(1L, n_genes, by = chunk)) {
  idx <- start:min(start + chunk - 1L, n_genes)
  vals <- matrix(formatC(log_expr[idx, , drop = FALSE], format = "f", digits = 4),
                 nrow = length(idx))
  body <- do.call(paste, c(list(rownames(log_expr)[idx]),
                           as.list(as.data.frame(vals, stringsAsFactors = FALSE)),
                           list(sep = "\t")))
  writeLines(body, con)
}
close(con)
on.exit()

# --- samples.csv ------------------------------------------------------------
samples <- well_dt[, .(sample, drug, time_h, donor, group)]
fwrite(samples, OUT_SAMPLES, sep = ",", quote = FALSE)

cat(sprintf("  mapped %d/%d ENTREZIDs; %d genes x %d wells\n",
            sum(mapped), length(gene_ids), n_genes, ncol(log_expr)))
cat(sprintf("  library sizes (mapped genes): %d..%d\n", min(lib_size), max(lib_size)))
cat(sprintf("GSE243674: wrote %s (%.1f MB, sha256 %s)\n", OUT_TSV,
            file.info(OUT_TSV)$size / 1e6, sha256_file(OUT_TSV)))
cat(sprintf("GSE243674: wrote %s (%d rows, sha256 %s)\n", OUT_SAMPLES,
            nrow(samples), sha256_file(OUT_SAMPLES)))
cat("GSE243674 build_expression done\n")
