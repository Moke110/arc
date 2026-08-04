#!/usr/bin/env Rscript
#
# Shared ortholog-mapping helper library for disease_data/miri.
#
# Provides:
#   - load_landmark_genes()   : read the 978 L1000 landmark gene symbols
#   - mouse2human()           : mouse symbol -> human symbol (homologene)
#   - rat2human()             : rat symbol -> human symbol (Entrez via homologene)
#   - species2human()         : generic taxid -> human symbol (homologene)
#   - map_deg_to_l1000()      : full pipeline: DEG results -> human L1000 up/down lists
#
# Sourced by the per-GSE build_deg.R scripts. All paths are resolved relative
# to this file so the helper works from any checkout location.

suppressPackageStartupMessages({
  library(org.Mm.eg.db)
  library(org.Rn.eg.db)
  library(org.Hs.eg.db)
  library(homologene)
  library(annotate)
})

SCRIPT_DIR <- dirname(normalizePath(
  sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1])
))
MIRI_DIR <- dirname(SCRIPT_DIR)   # disease_data/miri
L1000_CSV <- Sys.getenv(
  "ARC_L1000_CSV",
  unset = normalizePath(file.path(MIRI_DIR, "..", "..", "pert_model", "xpert", "metadata", "l1000_gene_info_978.csv"))
)

# `%||%` — null-coalescing operator (base R equivalent)
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# L1000 landmark genes
# ---------------------------------------------------------------------------

load_landmark_genes <- function() {
  if (!file.exists(L1000_CSV)) {
    stop(sprintf("L1000 gene list not found: %s\nSet ARC_L1000_CSV to the path of l1000_gene_info_978.csv", L1000_CSV))
  }
  l1000 <- read.csv(L1000_CSV, stringsAsFactors = FALSE)
  unique(l1000$gene_symbol)
}

# ---------------------------------------------------------------------------
# Species -> human symbol mapping
# ---------------------------------------------------------------------------

mouse2human <- function(mouse_genes) {
  if (length(mouse_genes) == 0 || all(is.na(mouse_genes))) {
    return(setNames(rep(NA_character_, length(mouse_genes)), mouse_genes))
  }
  m <- homologene::mouse2human(mouse_genes)
  if (is.data.frame(m)) {
    setNames(m[[2]], m[[1]])
  } else {
    m
  }
}

rat2human <- function(rat_genes) {
  result <- setNames(rep(NA_character_, length(rat_genes)), rat_genes)
  if (length(rat_genes) == 0 || all(is.na(rat_genes))) return(result)
  rat_entrez <- tryCatch(
    mapIds(org.Rn.eg.db, keys = rat_genes, column = "ENTREZID",
           keytype = "SYMBOL", multiVals = "first"),
    error = function(e) setNames(rep(NA_character_, length(rat_genes)), rat_genes)
  )
  valid <- !is.na(rat_entrez)
  if (sum(valid) == 0) return(result)

  he <- tryCatch(
    homologene(as.integer(rat_entrez[valid]), inTax = 10116, outTax = 9606),
    error = function(e) data.frame()
  )
  if (nrow(he) == 0) return(result)

  # homologene columns: 10116 (rat symbol), 9606 (human symbol)
  rat_sym <- as.character(he[[1]])
  human_sym <- as.character(he[[2]])
  for (i in seq_len(nrow(he))) {
    idx <- which(rat_genes == rat_sym[i])
    if (length(idx) > 0) result[idx[1]] <- human_sym[i]
  }
  result
}

species2human <- function(genes, in_tax) {
  # Exclude Ensembl-style IDs (e.g. ENSSCG...) which homologene cannot resolve
  valid <- genes[!grepl("^ENS", genes)]
  mapping <- tryCatch(
    homologene(valid, inTax = in_tax, outTax = 9606),
    error = function(e) data.frame()
  )
  result <- setNames(rep(NA_character_, length(genes)), genes)
  if (nrow(mapping) == 0) return(result)
  in_col <- as.character(in_tax)
  out_col <- "9606"
  if (!in_col %in% names(mapping) || !out_col %in% names(mapping)) return(result)
  hit <- match(genes, mapping[[in_col]])
  result[!is.na(hit)] <- as.character(mapping[[out_col]][hit[!is.na(hit)]])
  result
}

# ---------------------------------------------------------------------------
# Full DEG -> human L1000 mapping
# ---------------------------------------------------------------------------

map_deg_to_l1000 <- function(deg, gene_col = "Gene", logfc_col = "logFC",
                             p_col = "P.Value", species = c("mouse", "rat", "pig"),
                             p_threshold = 0.05, lfc_threshold = 0.3) {
  species <- match.arg(species)
  landmark <- load_landmark_genes()

  sig <- deg[!is.na(deg[[p_col]]) & deg[[p_col]] < p_threshold & abs(deg[[logfc_col]]) > lfc_threshold, ]
  if (nrow(sig) == 0) {
    sig <- deg[!is.na(deg[[p_col]]) & deg[[p_col]] < p_threshold, ]
  }

  mapper <- switch(species,
    mouse = function(g) mouse2human(g),
    rat = function(g) rat2human(g),
    pig = function(g) species2human(g, 9823)
  )

  up <- sig[sig[[logfc_col]] > 0, ]
  down <- sig[sig[[logfc_col]] < 0, ]

  up_map <- mapper(up[[gene_col]])
  down_map <- mapper(down[[gene_col]])
  if (is.null(names(up_map))) names(up_map) <- up[[gene_col]]
  if (is.null(names(down_map))) names(down_map) <- down[[gene_col]]

  up_df <- data.frame(
    species_gene = up[[gene_col]],
    human_gene = unname(up_map[up[[gene_col]]]),
    logFC = up[[logfc_col]],
    P.Value = up[[p_col]],
    stringsAsFactors = FALSE
  )
  down_df <- data.frame(
    species_gene = down[[gene_col]],
    human_gene = unname(down_map[down[[gene_col]]]),
    logFC = down[[logfc_col]],
    P.Value = down[[p_col]],
    stringsAsFactors = FALSE
  )

  to_l1000 <- function(df, dir_) {
    df <- df[!is.na(df$human_gene) & df$human_gene %in% landmark, ]
    if (nrow(df) == 0) return(df)
    df <- df[order(abs(df$logFC), decreasing = TRUE), ]
    df <- df[!duplicated(df$human_gene), ]
    df$direction <- dir_
    df[, c("human_gene", "direction", "logFC", "P.Value")]
  }

  rbind(to_l1000(up_df, "UP"), to_l1000(down_df, "DOWN"))
}
