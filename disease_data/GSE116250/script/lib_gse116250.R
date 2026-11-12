#!/usr/bin/env Rscript
#
# Per-item sample annotation + expression loader for GSE116250 (human
# left-ventricle bulk RNA-seq, 50 HF / 14 NF) — shared by script/build_deg.R and
# script/build_expression.R so the labels and the normalization live in one
# place.
#
# SAMPLE LABELS. GSE116250_rpkm.txt.gz has two annotation columns (`Gene` =
# Ensembl ID, `Common_name` = HGNC symbol) followed by 64 sample columns whose
# names already carry the group prefix: NF1..NF15 (14 non-failing donors),
# DCM2..DCM78 (37 dilated cardiomyopathy) and ICM47..ICM64 (13 ischaemic
# cardiomyopathy). 50 HF = 37 DCM + 13 ICM. This matches the GEO series matrix
# (Sample_title / characteristics `disease:`), which is not part of source/ under
# ticket #179's download scope; the group is parsed from the column names, i.e.
# from source/ itself.
#
# NORMALIZATION. The deposited values are linear-scale RPKM (Cufflinks against
# Ensembl GRCh37.82, mitochondrial genes removed -> 57,974 annotations; GEO
# sample_data_processing). This item works in log2(RPKM + 1), the standard
# log-scale transform for RPKM data and the scale limma assumes; no
# library-size renormalization is applied because RPKM already is one.
#
# SENTINEL. 656 cells across 65 genes carry the value 999999999, an obvious
# invalid-RPKM fill/overflow code (max valid RPKM in the file is 687,869; the
# sentinel genes are otherwise ordinary, e.g. MYH7). They are set to NA before
# the log transform and stay NA in the exported matrix; limma fits each gene on
# its non-missing samples (see build_deg.R), so no gene is dropped for a sentinel.

GSE116250_SENTINEL <- 999999999

# Load and collapse the RPKM matrix. Returns
#   list(expr = log2(RPKM+1) genes x samples matrix with NA for sentinels,
#        sample = sample column names,
#        group = "NF"/"DCM"/"ICM",
#        n_sentinel = number of sentinel cells replaced)
load_gse116250_expression <- function(source_dir) {
  f <- file.path(source_dir, "GSE116250_rpkm.txt.gz")
  if (!file.exists(f)) stop("GSE116250_rpkm.txt.gz not found in ", source_dir,
                            " — run script/fetch.sh first")
  dt <- data.table::fread(f, header = TRUE, data.table = FALSE, showProgress = FALSE)
  if (!identical(names(dt)[1:2], c("Gene", "Common_name"))) {
    stop("unexpected first two columns: ", paste(names(dt)[1:2], collapse = ", "))
  }
  ensg <- as.character(dt[[1]])
  sym <- as.character(dt[[2]])
  samples <- names(dt)[-(1:2)]
  x <- as.matrix(dt[, -(1:2), drop = FALSE])
  storage.mode(x) <- "double"
  n_sentinel <- sum(x == GSE116250_SENTINEL, na.rm = TRUE)
  x[x == GSE116250_SENTINEL] <- NA
  logx <- log2(x + 1)

  # Collapse to one row per Ensembl gene and then one row per HGNC symbol,
  # keeping the row with the highest mean log-expression (deterministic).
  ord <- order(ensg, -rowMeans(logx, na.rm = TRUE), seq_along(ensg))
  logx <- logx[ord, , drop = FALSE]
  ensg <- ensg[ord]
  sym <- sym[ord]
  keep <- !duplicated(ensg)
  logx <- logx[keep, , drop = FALSE]
  sym <- sym[keep]
  ord <- order(sym, -rowMeans(logx, na.rm = TRUE), seq_along(sym))
  logx <- logx[ord, , drop = FALSE]
  sym <- sym[ord]
  keep <- !duplicated(sym)
  logx <- logx[keep, , drop = FALSE]
  sym <- sym[keep]
  rownames(logx) <- sym
  colnames(logx) <- samples

  group <- sub("[0-9]+$", "", samples)
  if (!all(group %in% c("NF", "DCM", "ICM"))) {
    stop("unexpected sample names: ", paste(samples[!group %in% c("NF", "DCM", "ICM")],
                                            collapse = ", "))
  }
  list(expr = logx, sample = samples, group = group, n_sentinel = n_sentinel)
}

# Keep genes that are expressed (RPKM > 0, i.e. log2(RPKM+1) > 0) in at least
# `min_frac` of the samples; NA counts as not expressed (documented in the JSON).
filter_expressed <- function(expr, min_frac = 0.25) {
  keep <- rowMeans(!is.na(expr) & expr > 0) >= min_frac
  expr[keep, , drop = FALSE]
}
