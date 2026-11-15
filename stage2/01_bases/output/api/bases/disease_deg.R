#!/usr/bin/env Rscript
# disease_deg.R — the one live DEG implementation for the disease-side sources.
#
# ADR 0066 / spec D9: every published disease DEG is produced by one
# implementation, appropriate to the deposition. This is NOT the retired
# evaluation chain (`deseq2_contrasts.R` fitted tenv parent-dose truth); it is a
# new limma-based fitter written for the CM-pseudobulk and normalised-expression
# sources:
#
#   * integer counts      -> limma-voom (TMM size factors, voom precision
#                            weights, moderated t). A count model.
#   * normalised expression (e.g. log1p(CP10K) pseudobulk) -> limma on the
#                            supplied log2-scale matrix (moderated t). A
#                            documented linear model.
#
# One joint design is always fitted once; every contrast of that fit is written
# as a row set of the long output table. A shared-control source therefore does
# not get separate case-control fits.
#
# Usage
#   R_HOME=/path/to/arc/lib/R /path/to/arc/lib/R/bin/exec/R --no-save --slave \
#     -f disease_deg.R --args --matrix M.tsv.gz --samples S.csv \
#     --design D.json --out O.csv
#
# Inputs
#   --matrix  TSV (optionally gzipped), header `source_gene_id`, `gene_symbol`,
#             then one column per sample in the sample sheet's row order.
#   --samples CSV with `sample` and the design factor column named by the design
#             spec (`factor`, default `group`). Row order must equal the matrix
#             column order; a mismatch fails.
#   --design  JSON: {"kind":"counts"|"expression", "formula":"~0 + group",
#             "factor":"group", "contrasts":[{"label":"NICM-Donor",
#             "numerator":"NICM","denominator":"Donor"}, ...]}
#             Optional, both off by default:
#               "block":"donor"  a sample-sheet column naming a blocking factor;
#                   one consensus intra-block correlation is estimated with
#                   duplicateCorrelation() and given to voom() and lmFit(), for
#                   designs where one subject contributes repeated profiles
#                   (e.g. the Kuppe-MI donor x zone pseudobulk). For a repeated-
#                   measures design the two calls must share the correlation.
#               "filter":{"min_count":10,"min_total_count":15,"large_n":10,
#                   "min_prop":0.7}  a gene filter on the counts, applying the
#                   edgeR::filterByExpr rule (minimum group size, min.count CPM
#                   cutoff at the median library size, min.total.count) without
#                   edgeR: the pinned `arc` R env ships limma but not edgeR, so
#                   the rule is written out here rather than imported. Genes
#                   with an all-zero row are always dropped; the filter only
#                   removes more, and the dropped count is reported.
#   --out     long CSV: source_gene_id, gene_symbol, contrast, log2fc, stat,
#             p_value, padj, base_mean (one row per tested gene per contrast).
#
# The last stdout line is `DEG_SUMMARY <json>` with the fit description
# (limma version, kind, design, n_genes_input/tested, coefficient names); the
# caller records it in the item's tracked DEG summary.
suppressPackageStartupMessages(library(limma))

fail <- function(...) {
  message("disease_deg.R: ", paste0(...))
  quit(status = 1)
}

parse_args <- function(argv) {
  values <- list()
  i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--") || i == length(argv)) {
      fail("bad argument: ", key)
    }
    values[[substring(key, 3)]] <- argv[[i + 1]]
    i <- i + 2
  }
  values
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
for (required in c("matrix", "samples", "design", "out")) {
  if (is.null(args[[required]])) fail("missing --", required)
}

read_design <- function(path) {
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    return(jsonlite::fromJSON(path, simplifyVector = FALSE))
  }
  fail("jsonlite is required to read the design spec")
}
design_spec <- read_design(args$design)
kind <- design_spec$kind
if (!kind %in% c("counts", "expression")) fail("kind must be counts or expression")
factor_name <- if (is.null(design_spec$factor)) "group" else design_spec$factor
formula <- design_spec$formula
if (is.null(formula)) fail("design spec needs a formula")
block_name <- design_spec$block
filter_spec <- design_spec$filter
if (!is.null(filter_spec) && identical(filter_spec, FALSE)) filter_spec <- NULL

# The edgeR::filterByExpr rule, written out because the pinned `arc` R env has no
# edgeR (see the module header). lib.size is the raw column sum, which is what
# voom's own normalisation uses, so the filter and the model agree on the library
# size; the rule's parameters are the edgeR defaults unless the design overrides
# them.
count_filter <- function(values, group, spec) {
  min_count <- if (is.null(spec$min_count)) 10 else spec$min_count
  min_total_count <- if (is.null(spec$min_total_count)) 15 else spec$min_total_count
  large_n <- if (is.null(spec$large_n)) 10 else spec$large_n
  min_prop <- if (is.null(spec$min_prop)) 0.7 else spec$min_prop
  sizes <- tabulate(as.integer(factor(group)))
  sizes <- sizes[sizes > 0L]
  min_sample_size <- min(sizes)
  if (min_sample_size > large_n) {
    min_sample_size <- large_n + (min_sample_size - large_n) * min_prop
  }
  lib_size <- colSums(values)
  cutoff <- min_count / median(lib_size) * 1e6
  tol <- 1e-14
  cpm <- t(t(values) / lib_size) * 1e6
  (rowSums(cpm >= cutoff) >= (min_sample_size - tol)) &
    (rowSums(values) >= min_total_count - tol)
}

matrix_sep <- if (grepl("\\.csv(\\.gz)?$", args$matrix)) "," else "\t"
matrix <- read.delim(args$matrix, check.names = FALSE, stringsAsFactors = FALSE,
                     sep = matrix_sep)
if (ncol(matrix) < 3) fail("matrix needs two id columns and at least one sample")
gene_id <- as.character(matrix[[1]])
gene_symbol <- as.character(matrix[[2]])
values <- as.matrix(matrix[, -(1:2), drop = FALSE])
storage.mode(values) <- "double"
rownames(values) <- gene_id
n_genes_matrix <- length(gene_id)
consensus_correlation <- NULL

samples <- read.csv(args$samples, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("sample", factor_name) %in% colnames(samples))) {
  fail("sample sheet lacks sample/", factor_name)
}
if (nrow(samples) != ncol(values)) fail("sample sheet and matrix disagree on sample count")
if (!identical(as.character(samples$sample), colnames(values))) {
  fail("sample sheet row order differs from matrix column order")
}
samples[[factor_name]] <- factor(samples[[factor_name]])
block <- NULL
if (!is.null(block_name)) {
  if (!is.character(block_name) || length(block_name) != 1L ||
      !block_name %in% colnames(samples)) {
    fail("design block must name one sample-sheet column")
  }
  block <- factor(samples[[block_name]])
  if (nlevels(block) < 2L) fail("the block factor has fewer than two levels")
}
design <- model.matrix(as.formula(formula), data = samples)
if (nrow(design) != nrow(samples)) fail("design matrix row count mismatch")
coefficient_names <- colnames(design)

if (identical(kind, "counts")) {
  if (any(!is.finite(values)) || any(values < 0) || any(values != round(values))) {
    fail("counts mode requires finite non-negative integer counts")
  }
  # The pinned `arc` R env ships limma but not edgeR, so TMM is not available and
  # is deliberately not re-implemented here: voom's own library-size
  # normalisation (total counts + 0.5 offset, CPM scaling) is the documented
  # count model for these pseudobulk profiles.
  keep <- rowSums(values) > 0
  if (!is.null(filter_spec)) {
    keep <- keep & count_filter(values, samples[[factor_name]], filter_spec)
  }
  values <- values[keep, , drop = FALSE]
  gene_id <- gene_id[keep]
  gene_symbol <- gene_symbol[keep]
  n_genes_dropped <- n_genes_matrix - length(gene_id)
  if (nrow(values) == 0) fail("no gene with a positive count")
  voom_fit <- voom(values, design, plot = FALSE)
  if (is.null(block)) {
    fit <- lmFit(voom_fit, design)
  } else {
    consensus_correlation <- duplicateCorrelation(voom_fit, design, block = block)
    consensus_correlation <- consensus_correlation$consensus.correlation
    voom_fit <- voom(values, design, block = block,
                     correlation = consensus_correlation, plot = FALSE)
    fit <- lmFit(voom_fit, design, block = block, correlation = consensus_correlation)
  }
  fit <- eBayes(fit)
  tested <- rownames(voom_fit$E)
  base_mean <- voom_fit$E
} else {
  if (any(!is.finite(values))) {
    keep <- rowSums(is.finite(values)) == ncol(values)
    values <- values[keep, , drop = FALSE]
    gene_id <- gene_id[keep]
    gene_symbol <- gene_symbol[keep]
    if (nrow(values) == 0) fail("no gene is finite in every sample")
  }
  n_genes_dropped <- n_genes_matrix - length(gene_id)
  if (is.null(block)) {
    fit <- lmFit(values, design)
  } else {
    consensus_correlation <- duplicateCorrelation(values, design, block = block)
    consensus_correlation <- consensus_correlation$consensus.correlation
    fit <- lmFit(values, design, block = block, correlation = consensus_correlation)
  }
  fit <- eBayes(fit)
  tested <- rownames(values)
  base_mean <- values
}

blocks <- list()
for (contrast in design_spec$contrasts) {
  numerator <- paste0(factor_name, contrast$numerator)
  denominator <- paste0(factor_name, contrast$denominator)
  if (!numerator %in% coefficient_names) fail("unknown coefficient: ", numerator)
  if (!denominator %in% coefficient_names) fail("unknown coefficient: ", denominator)
  expression <- paste0(numerator, " - ", denominator)
  contrast_matrix <- makeContrasts(contrasts = expression, levels = design)
  fit_contrast <- eBayes(contrasts.fit(fit, contrast_matrix))
  table <- topTable(fit_contrast, coef = 1, number = Inf, sort.by = "none",
                    adjust.method = "BH")
  symbol <- gene_symbol[match(rownames(table), gene_id)]
  blocks[[length(blocks) + 1]] <- data.frame(
    source_gene_id = rownames(table),
    gene_symbol = ifelse(is.na(symbol), "", symbol),
    contrast = contrast$label,
    log2fc = table$logFC,
    stat = table$t,
    p_value = table$P.Value,
    padj = table$adj.P.Val,
    base_mean = table$AveExpr,
    stringsAsFactors = FALSE
  )
}
result <- do.call(rbind, blocks)
if (is.null(result) || nrow(result) == 0) fail("fit produced no DEG rows")
write.csv(result, args$out, row.names = FALSE, quote = FALSE, na = "NA")

summary <- list(
  implementation = "output/api/bases/disease_deg.R",
  kind = kind,
  formula = formula,
  factor = factor_name,
  coefficient_names = I(coefficient_names),
  contrasts = I(vapply(design_spec$contrasts, function(x) x$label, character(1))),
  n_genes_input = length(gene_id),
  n_genes_matrix = n_genes_matrix,
  n_genes_dropped = n_genes_dropped,
  filter = if (is.null(filter_spec)) NULL else filter_spec,
  block = if (is.null(block_name)) NULL else block_name,
  n_blocks = if (is.null(block)) NULL else nlevels(block),
  consensus_correlation = consensus_correlation,
  n_genes_tested = length(tested),
  n_rows = nrow(result),
  limma_version = as.character(packageVersion("limma")),
  r_version = R.version.string
)
cat("DEG_SUMMARY ", jsonlite::toJSON(summary, auto_unbox = TRUE), "\n", sep = "")
