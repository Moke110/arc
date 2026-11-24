# Exploratory parallel contrasts against the DMSO baseline of one plate.
suppressPackageStartupMessages(library(limma))
args <- commandArgs(trailingOnly = TRUE)
metadata <- read.csv(args[[2]], stringsAsFactors = FALSE)
n_genes <- as.integer(args[[3]])
con <- file(args[[1]], "rb")
counts <- matrix(readBin(con, integer(), n_genes * nrow(metadata), size = 4,
                         endian = "little"), nrow = n_genes)
close(con)

group <- relevel(factor(metadata$group), ref = "DMSO")
design <- model.matrix(~group)
observed <- rowSums(counts) > 0L
expression <- voom(counts[observed, , drop = FALSE], design,
                   lib.size = metadata$depth, normalize.method = "none")
fit <- eBayes(lmFit(expression, design))
ids <- which(observed) - 1L
tables <- lapply(seq.int(2L, ncol(design)), function(coef_index) {
  p <- p.adjust(fit$p.value[, coef_index], method = "BH")
  lfc <- fit$coefficients[, coef_index]
  selected <- which(is.finite(p) & p < 0.05 & abs(lfc) > 1.5)
  data.frame(group = rep(levels(group)[coef_index], length(selected)),
             gene_index = ids[selected], log2fc = lfc[selected], padj = p[selected])
})
write.csv(do.call(rbind, tables), args[[4]], row.names = FALSE)
