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

matrix_sep <- if (grepl("\\.csv(\\.gz)?$", args$matrix)) "," else "\t"
matrix <- read.delim(args$matrix, check.names = FALSE, stringsAsFactors = FALSE,
                     sep = matrix_sep)
if (ncol(matrix) < 3) fail("matrix needs two id columns and at least one sample")
gene_id <- as.character(matrix[[1]])
gene_symbol <- as.character(matrix[[2]])
values <- as.matrix(matrix[, -(1:2), drop = FALSE])
storage.mode(values) <- "double"
rownames(values) <- gene_id

samples <- read.csv(args$samples, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("sample", factor_name) %in% colnames(samples))) {
  fail("sample sheet lacks sample/", factor_name)
}
if (nrow(samples) != ncol(values)) fail("sample sheet and matrix disagree on sample count")
if (!identical(as.character(samples$sample), colnames(values))) {
  fail("sample sheet row order differs from matrix column order")
}
samples[[factor_name]] <- factor(samples[[factor_name]])
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
  values <- values[keep, , drop = FALSE]
  gene_id <- gene_id[keep]
  gene_symbol <- gene_symbol[keep]
  if (nrow(values) == 0) fail("no gene with a positive count")
  voom_fit <- voom(values, design, plot = FALSE)
  fit <- lmFit(voom_fit, design)
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
  fit <- lmFit(values, design)
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
  n_genes_tested = length(tested),
  n_rows = nrow(result),
  limma_version = as.character(packageVersion("limma")),
  r_version = R.version.string
)
cat("DEG_SUMMARY ", jsonlite::toJSON(summary, auto_unbox = TRUE), "\n", sep = "")
