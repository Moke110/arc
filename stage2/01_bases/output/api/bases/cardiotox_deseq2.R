#!/usr/bin/env Rscript
# Pooled class-vs-vehicle DESeq2 contrasts from the GSE243674 read counts.
#
# The item's own build_deg.R fits DESeq2 per drug x time against the same-time
# vehicle wells with the donor-paired design `~ donor + group`. Ticket #189 pools
# a mechanism class's member wells into one Disease group, so this script runs
# the same model and mapping rule on the pooled counts: for each time point the
# class members' wells are the CLASS arm and the vehicle wells are the control
# arm. The ENTREZID -> HGNC symbol mapping mirrors the item's documented rule
# (approved rows only; protein-coding < non-coding RNA < pseudogene < other, then
# lowest hgnc_id; several Entrez IDs mapping to one symbol are summed).
#
# Usage:
#   Rscript cardiotox_deseq2.R COUNTS_TSV SAMPLES_TSV HGNC_TSV OUT_TSV
#
# COUNTS_TSV: ENTREZID x wells (first column ENTREZID, tab separated)
# SAMPLES_TSV: well, donor, group, time_h
# OUT_TSV: time_h, gene_symbol, log2FC, stat, p_value, padj, baseMean,
#          n_case, n_control

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: cardiotox_deseq2.R COUNTS_TSV SAMPLES_TSV HGNC_TSV OUT_TSV")
}
counts_path <- args[1]
samples_path <- args[2]
hgnc_path <- args[3]
out_path <- args[4]

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
})

counts_dt <- fread(counts_path, sep = "\t", header = TRUE)
gene_ids <- as.character(counts_dt[[1]])
counts <- as.matrix(counts_dt[, -1, with = FALSE])
rownames(counts) <- gene_ids
samples <- fread(samples_path, sep = "\t", header = TRUE)
if (!identical(colnames(counts), samples$well)) {
  stop("the counts columns and the samples table are not in the same order")
}

# --- ENTREZID -> HGNC symbol with the pinned table (the item's rule) --------
hg <- fread(hgnc_path, sep = "\t", header = TRUE, quote = "",
            select = c("hgnc_id", "symbol", "status", "locus_group", "entrez_id"))
hg <- hg[status == "Approved" & !is.na(entrez_id) & entrez_id != ""]
hg[, entrez_id := as.character(entrez_id)]
hg[, locus_rank := match(locus_group,
                         c("protein-coding gene", "non-coding RNA", "pseudogene"),
                         nomatch = 99L)]
hg[, hgnc_num := as.integer(sub("^HGNC:", "", hgnc_id))]
setorder(hg, entrez_id, locus_rank, hgnc_num)
mapping <- hg[!duplicated(entrez_id), .(entrez_id, symbol)]
symbols <- mapping$symbol[match(gene_ids, mapping$entrez_id)]
keep <- !is.na(symbols)
counts <- counts[keep, , drop = FALSE]
symbols <- symbols[keep]
duplicated_symbols <- unique(symbols[duplicated(symbols)])
if (length(duplicated_symbols) > 0L) {
  counts <- rowsum(counts, group = symbols)
  symbols <- rownames(counts)
}
rownames(counts) <- symbols

results <- list()
for (time_h in sort(unique(samples$time_h))) {
  if (length(unique(samples$group[samples$time_h == time_h])) != 2L) {
    next
  }
  part <- samples[time_h == time_h]
  matrix <- counts[, part$well, drop = FALSE]
  coldata <- data.frame(
    donor = factor(part$donor),
    group = factor(part$group, levels = c("vehicle", "CLASS")),
    row.names = part$well)
  dds <- DESeqDataSetFromMatrix(matrix, coldata, design = ~ donor + group)
  dds <- DESeq(dds, quiet = TRUE)
  result <- results(dds, contrast = c("group", "CLASS", "vehicle"))
  results[[length(results) + 1L]] <- data.frame(
    time_h = time_h,
    gene_symbol = rownames(result),
    log2FC = result$log2FoldChange,
    stat = result$stat,
    p_value = result$pvalue,
    padj = result$padj,
    baseMean = result$baseMean,
    n_case = sum(part$group == "CLASS"),
    n_control = sum(part$group == "vehicle"),
    stringsAsFactors = FALSE)
}
output <- do.call(rbind, results)
format_number <- function(values) {
  out <- rep("", length(values))
  finite <- is.finite(values)
  out[finite] <- sprintf("%.10g", values[finite])
  out
}
output$log2FC <- format_number(output$log2FC)
output$stat <- format_number(output$stat)
output$p_value <- format_number(output$p_value)
output$padj <- format_number(output$padj)
output$baseMean <- format_number(output$baseMean)
write.table(output, out_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat("DESeq2 pooled class contrasts:", length(results), "time points,", nrow(output),
    "rows\n")
