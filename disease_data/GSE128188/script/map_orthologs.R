#!/usr/bin/env Rscript
#
# Ortholog mapping for rat (GSE108940) and pig (GSE137293) DEG results.
#
# Reads a per-GSE DEG results CSV from disease_data/{GSE}/source/ and writes
# the human L1000-mapped per-GSE DEG list to disease_data/{GSE}/local/deg/.
#
# Usage:
#   Rscript script/map_orthologs.R \
#     --gse GSE108940 --deg source/GSE108940_WT_deg_results.csv \
#     --out local/deg/GSE108940_human_l1000.csv --species rat
#   Rscript script/map_orthologs.R \
#     --gse GSE137293 --deg source/GSE137293_deg_results.csv \
#     --out local/deg/GSE137293_human_l1000.csv --species pig
#
# Column conventions differ per source:
#   rat:  Gene, logFC, P.Value      (limma output, run_deg_analysis.R)
#   pig:  gene, log2FC, P_Value     (scipy Welch t-test output)

suppressPackageStartupMessages({
  library(optparse)
})

SCRIPT_DIR <- dirname(normalizePath(
  sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1])
))
source(file.path(SCRIPT_DIR, "map_to_human.R"))

parse_args <- function() {
  opt <- make_option(c("--deg"), type = "character", help = "Path to per-GSE DEG results CSV")
  opt <- c(opt, make_option(c("--out"), type = "character", help = "Output CSV path (human L1000 list)"))
  opt <- c(opt, make_option(c("--species"), type = "character",
                            help = "Species: rat or pig"))
  parser <- OptionParser(option_list = opt, description = "Map rat/pig DEG results to human L1000")
  parse_args2(parser, args = commandArgs(trailingOnly = TRUE))
}

main <- function() {
  args <- parse_args()
  for (need in c("deg", "out", "species")) {
    if (is.null(args[[need]])) stop(sprintf("Missing required argument: --%s", need))
  }
  if (!args$species %in% c("rat", "pig")) stop("--species must be rat or pig")

  deg <- read.csv(args$deg, stringsAsFactors = FALSE)
  mapped <- if (args$species == "rat") {
    map_deg_to_l1000(deg, gene_col = "Gene", logfc_col = "logFC",
                     p_col = "P.Value", species = "rat")
  } else {
    map_deg_to_l1000(deg, gene_col = "gene", logfc_col = "log2FC",
                     p_col = "P_Value", species = "pig")
  }

  dir.create(dirname(args$out), recursive = TRUE, showWarnings = FALSE)
  write.csv(mapped, args$out, row.names = FALSE)
  cat(sprintf("Saved %s: %d UP + %d DOWN (human L1000)\n",
              args$out, sum(mapped$direction == "UP"), sum(mapped$direction == "DOWN")))
}

main()
