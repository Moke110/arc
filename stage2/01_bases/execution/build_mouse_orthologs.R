#!/usr/bin/env Rscript
#
# Generate the pinned mouse -> human ortholog table for the Stage 2 gene-axis
# evidence product (ticket #183, spec section C).
#
# The repository has no axis-level mouse->human table. This script derives one
# deterministically from the two mouse MIRI DEG items with the pinned `arc`
# environment's homologene package (the same `mouse2human` definition used by the
# shared per-item ortholog helper in each disease item's `script/map_to_human.R`),
# writes the small table to the gene-axis evidence product and prints the pinned
# package versions. The
# table is then a tracked, digest-sealed input of the gene-axis build; re-running
# this script reproduces it byte for byte.
#
# Run with the pinned `arc` R entry point (the env's `R` wrapper is mis-patched):
#   R_HOME=/home/chang/miniconda3/envs/arc/lib/R \
#     /home/chang/miniconda3/envs/arc/lib/R/bin/exec/R --no-save --slave \
#     --file stage2/01_bases/execution/build_mouse_orthologs.R

suppressPackageStartupMessages({
  library(homologene)
  library(org.Mm.eg.db)
})

find_repo_root <- function(start_dir) {
  dir <- normalizePath(start_dir, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(dir, ".git"))) return(dir)
    parent <- dirname(dir)
    if (identical(parent, dir)) {
      stop("cannot find the repository root (no .git directory above ", start_dir, ")")
    }
    dir <- parent
  }
}

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (!length(script_arg)) stop("run this script with R --file=...")
SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
REPO_ROOT <- find_repo_root(SCRIPT_DIR)

mouse_sources <- c(
  file.path(REPO_ROOT, "disease_data/GSE160516/source/GSE160516_deg_results.csv"),
  file.path(REPO_ROOT, "disease_data/GSE108940/source/GSE108940_WT_deg_results.csv")
)
missing <- mouse_sources[!file.exists(mouse_sources)]
if (length(missing)) stop("missing mouse DEG source: ", paste(missing, collapse = ", "))

mouse_symbols <- unique(unlist(lapply(mouse_sources, function(path) {
  as.character(read.csv(path, stringsAsFactors = FALSE)$Gene)
})))
mouse_symbols <- sort(mouse_symbols)

mapped <- homologene::mouse2human(mouse_symbols)
table <- unique(data.frame(
  mouse_symbol = as.character(mapped$mouseGene),
  human_symbol = as.character(mapped$humanGene),
  stringsAsFactors = FALSE
))
table <- table[order(table$mouse_symbol, table$human_symbol), , drop = FALSE]
table <- table[table$mouse_symbol != "" & table$human_symbol != "", , drop = FALSE]

output <- file.path(REPO_ROOT, "stage2/01_bases/output/gene-axis/mouse_human_orthologs.csv")
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.csv(table, output, row.names = FALSE, quote = FALSE)

cat("mouse symbols:", length(mouse_symbols), "\n")
cat("mapping rows:", nrow(table), "\n")
cat("distinct mapped mouse symbols:", length(unique(table$mouse_symbol)), "\n")
cat("R:", R.version.string, "\n")
cat("homologene:", as.character(packageVersion("homologene")), "\n")
cat("org.Mm.eg.db:", as.character(packageVersion("org.Mm.eg.db")), "\n")
cat("wrote:", output, "\n")
