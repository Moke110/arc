#!/usr/bin/env Rscript
#
# Ortholog mapping: Rat (GSE108940) + Pig (GSE137293) → Human gene symbols
# Uses homologene R package for proper ortholog inference.
#
# Outputs:
#   geo/GSE108940_human_l1000.csv   — rat DEGs mapped to human L1000
#   geo/GSE137293_human_l1000.csv   — pig DEGs mapped to human L1000

suppressPackageStartupMessages({
  library(homologene)
})

# Load L1000 landmark genes
l1000_csv <- "/home/chang/arc/xpert/metadata/l1000_gene_info_978.csv"
l1000 <- read.csv(l1000_csv, stringsAsFactors = FALSE)
landmark_human <- unique(l1000$gene_symbol)
cat(sprintf("L1000 landmark genes: %d\n", length(landmark_human)))

# ============================================================
# Helper: map species genes to human gene symbols
# ============================================================
map_to_human <- function(genes, in_tax) {
  # homologene() takes gene symbols, returns data.frame with homologues
  result <- homologene(genes, inTax = in_tax, outTax = 9606)

  if (nrow(result) == 0) {
    return(data.frame(
      input_gene = genes,
      human_gene = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  # result has columns: input gene info + homolog info
  # The homolog symbol is in a column — need to check exact structure
  # Typically: `1`=input gene ID, `2`=input gene symbol, `3`=input tax,
  #            `4`=homolog gene ID, `5`=homolog symbol, `6`=homolog tax

  # result has columns named by tax IDs: e.g. "10116" (input symbol), "9606" (homolog symbol)
  in_col <- as.character(in_tax)
  out_col <- "9606"
  mapping <- result[, c(in_col, out_col)]
  names(mapping) <- c("input_gene", "human_gene")

  # Merge back to original gene list
  df <- data.frame(input_gene = genes, stringsAsFactors = FALSE)
  df <- merge(df, mapping, by = "input_gene", all.x = TRUE)

  # Deduplicate: keep first mapping per input gene
  df <- df[!duplicated(df$input_gene), ]

  return(df)
}

# ============================================================
# Process GSE108940 (Rat, taxid=10116)
# ============================================================
cat("\n--- GSE108940 (Rat → Human) ---\n")

d2 <- read.csv("/home/chang/arc/MIRI_fda_xes/geo/GSE108940_WT_deg_results.csv")

# Filter: P < 0.05, |logFC| > 0.3
g2 <- d2[d2$P.Value < 0.05 & abs(d2$logFC) > 0.3, ]
cat(sprintf("Rat DEGs filtered: %d (P<0.05, |logFC|>0.3)\n", nrow(g2)))

g2_up <- g2[g2$logFC > 0, ]
g2_down <- g2[g2$logFC < 0, ]
cat(sprintf("  Up: %d, Down: %d\n", nrow(g2_up), nrow(g2_down)))

# Map to human
cat("Mapping rat up genes to human...\n")
rat_up_genes <- g2_up$Gene
up_map <- map_to_human(rat_up_genes, in_tax = 10116)
cat(sprintf("  Mapped: %d/%d (%.1f%%)\n",
            sum(!is.na(up_map$human_gene)), nrow(up_map),
            100 * sum(!is.na(up_map$human_gene)) / nrow(up_map)))

cat("Mapping rat down genes to human...\n")
rat_down_genes <- g2_down$Gene
down_map <- map_to_human(rat_down_genes, in_tax = 10116)
cat(sprintf("  Mapped: %d/%d (%.1f%%)\n",
            sum(!is.na(down_map$human_gene)), nrow(down_map),
            100 * sum(!is.na(down_map$human_gene)) / nrow(down_map)))

# Build result data frames with logFC
up_result <- data.frame(
  rat_gene = rat_up_genes,
  human_gene = up_map$human_gene[match(rat_up_genes, up_map$input_gene)],
  logFC = g2_up$logFC,
  P.Value = g2_up$P.Value,
  stringsAsFactors = FALSE
)
down_result <- data.frame(
  rat_gene = rat_down_genes,
  human_gene = down_map$human_gene[match(rat_down_genes, down_map$input_gene)],
  logFC = g2_down$logFC,
  P.Value = g2_down$P.Value,
  stringsAsFactors = FALSE
)

# Filter to L1000
up_l1000 <- up_result[!is.na(up_result$human_gene) & up_result$human_gene %in% landmark_human, ]
down_l1000 <- down_result[!is.na(down_result$human_gene) & down_result$human_gene %in% landmark_human, ]

cat(sprintf("Rat up in L1000: %d\n", nrow(up_l1000)))
cat(sprintf("Rat down in L1000: %d\n", nrow(down_l1000)))

# Remove duplicates (keep highest |logFC|)
up_l1000 <- up_l1000[order(abs(up_l1000$logFC), decreasing = TRUE), ]
up_l1000 <- up_l1000[!duplicated(up_l1000$human_gene), ]
down_l1000 <- down_l1000[order(abs(down_l1000$logFC), decreasing = TRUE), ]
down_l1000 <- down_l1000[!duplicated(down_l1000$human_gene), ]

cat(sprintf("After dedup: up=%d, down=%d\n", nrow(up_l1000), nrow(down_l1000)))

# Remove conflicts (gene assigned to both up and down)
conflict_genes <- intersect(up_l1000$human_gene, down_l1000$human_gene)
if (length(conflict_genes) > 0) {
  cat(sprintf("Conflict genes (in both up & down): %d\n", length(conflict_genes)))
  up_l1000 <- up_l1000[!up_l1000$human_gene %in% conflict_genes, ]
  down_l1000 <- down_l1000[!down_l1000$human_gene %in% conflict_genes, ]
}

# Write combined CSV
combined <- rbind(
  data.frame(human_gene = up_l1000$human_gene, direction = "UP",
             logFC = up_l1000$logFC, P.Value = up_l1000$P.Value,
             stringsAsFactors = FALSE),
  data.frame(human_gene = down_l1000$human_gene, direction = "DOWN",
             logFC = down_l1000$logFC, P.Value = down_l1000$P.Value,
             stringsAsFactors = FALSE)
)
write.csv(combined, "/home/chang/arc/MIRI_fda_xes/geo/GSE108940_human_l1000.csv",
          row.names = FALSE)
cat(sprintf("Saved: GSE108940_human_l1000.csv (%d genes)\n", nrow(combined)))


# ============================================================
# Process GSE137293 (Pig, taxid=9823)
# ============================================================
cat("\n--- GSE137293 (Pig → Human) ---\n")

# Read the pig DEG results generated by Python
pig_deg_path <- "/home/chang/arc/MIRI_fda_xes/geo/GSE137293_deg_results.csv"
if (file.exists(pig_deg_path)) {
  d3 <- read.csv(pig_deg_path)
  cat(sprintf("Pig DEGs loaded: %d\n", nrow(d3)))

  # Filter
  g3 <- d3[d3$P_Value < 0.05 & abs(d3$log2FC) > 0.3, ]
  if (nrow(g3) < 10) {
    g3 <- d3[d3$P_Value < 0.05, ]
    cat(sprintf("  Relaxed: %d (P<0.05, any |logFC|)\n", nrow(g3)))
  }

  g3_up <- g3[g3$log2FC > 0, ]
  g3_down <- g3[g3$log2FC < 0, ]
  cat(sprintf("  Up: %d, Down: %d\n", nrow(g3_up), nrow(g3_down)))

  # Map to human — skip ENSSCG genes (homologene needs gene symbols)
  pig_up_genes <- g3_up$gene
  pig_down_genes <- g3_down$gene

  # Filter out ENS* IDs for homologene (it needs valid gene symbols)
  pig_up_valid <- pig_up_genes[!grepl("^ENS", pig_up_genes)]
  pig_down_valid <- pig_down_genes[!grepl("^ENS", pig_down_genes)]

  n_ens_up <- sum(grepl("^ENS", pig_up_genes))
  n_ens_down <- sum(grepl("^ENS", pig_down_genes))
  cat(sprintf("  ENSSCG genes dropped: %d up + %d down\n", n_ens_up, n_ens_down))

  if (length(pig_up_valid) > 0) {
    cat(sprintf("Mapping %d pig up genes to human...\n", length(pig_up_valid)))
    pig_up_map <- map_to_human(pig_up_valid, in_tax = 9823)
    cat(sprintf("  Mapped: %d/%d\n", sum(!is.na(pig_up_map$human_gene)), nrow(pig_up_map)))
  } else {
    pig_up_map <- data.frame(input_gene = character(0), human_gene = character(0))
  }

  if (length(pig_down_valid) > 0) {
    cat(sprintf("Mapping %d pig down genes to human...\n", length(pig_down_valid)))
    pig_down_map <- map_to_human(pig_down_valid, in_tax = 9823)
    cat(sprintf("  Mapped: %d/%d\n", sum(!is.na(pig_down_map$human_gene)), nrow(pig_down_map)))
  } else {
    pig_down_map <- data.frame(input_gene = character(0), human_gene = character(0))
  }

  # Build result data frames
  up_result_p <- data.frame(
    pig_gene = pig_up_valid,
    human_gene = pig_up_map$human_gene[match(pig_up_valid, pig_up_map$input_gene)],
    logFC = g3_up$log2FC[match(pig_up_valid, g3_up$gene)],
    P.Value = g3_up$P_Value[match(pig_up_valid, g3_up$gene)],
    stringsAsFactors = FALSE
  )
  down_result_p <- data.frame(
    pig_gene = pig_down_valid,
    human_gene = pig_down_map$human_gene[match(pig_down_valid, pig_down_map$input_gene)],
    logFC = g3_down$log2FC[match(pig_down_valid, g3_down$gene)],
    P.Value = g3_down$P_Value[match(pig_down_valid, g3_down$gene)],
    stringsAsFactors = FALSE
  )

  # Filter to L1000
  up_p_l1000 <- up_result_p[!is.na(up_result_p$human_gene) & up_result_p$human_gene %in% landmark_human, ]
  down_p_l1000 <- down_result_p[!is.na(down_result_p$human_gene) & down_result_p$human_gene %in% landmark_human, ]

  cat(sprintf("Pig up in L1000: %d\n", nrow(up_p_l1000)))
  cat(sprintf("Pig down in L1000: %d\n", nrow(down_p_l1000)))

  # Dedup
  up_p_l1000 <- up_p_l1000[order(abs(up_p_l1000$logFC), decreasing = TRUE), ]
  up_p_l1000 <- up_p_l1000[!duplicated(up_p_l1000$human_gene), ]
  down_p_l1000 <- down_p_l1000[order(abs(down_p_l1000$logFC), decreasing = TRUE), ]
  down_p_l1000 <- down_p_l1000[!duplicated(down_p_l1000$human_gene), ]

  # Remove conflicts
  conflict_p <- intersect(up_p_l1000$human_gene, down_p_l1000$human_gene)
  if (length(conflict_p) > 0) {
    cat(sprintf("  Conflicts: %d, removing\n", length(conflict_p)))
    up_p_l1000 <- up_p_l1000[!up_p_l1000$human_gene %in% conflict_p, ]
    down_p_l1000 <- down_p_l1000[!down_p_l1000$human_gene %in% conflict_p, ]
  }

  # Build combined — handle empty case
  rows_list <- list()
  if (nrow(up_p_l1000) > 0) {
    rows_list[[1]] <- data.frame(human_gene = up_p_l1000$human_gene, direction = "UP",
               logFC = up_p_l1000$logFC, P.Value = up_p_l1000$P.Value,
               stringsAsFactors = FALSE)
  }
  if (nrow(down_p_l1000) > 0) {
    rows_list[[length(rows_list) + 1]] <- data.frame(human_gene = down_p_l1000$human_gene, direction = "DOWN",
               logFC = down_p_l1000$logFC, P.Value = down_p_l1000$P.Value,
               stringsAsFactors = FALSE)
  }
  if (length(rows_list) > 0) {
    combined_p <- do.call(rbind, rows_list)
  } else {
    combined_p <- data.frame(human_gene = character(0), direction = character(0),
                             logFC = numeric(0), P.Value = numeric(0))
  }
  write.csv(combined_p, "/home/chang/arc/MIRI_fda_xes/geo/GSE137293_human_l1000.csv",
            row.names = FALSE)
  cat(sprintf("Saved: GSE137293_human_l1000.csv (%d genes)\n", nrow(combined_p)))
} else {
  cat(sprintf("Pig DEG results not found at %s — skipping pig mapping\n", pig_deg_path))
  # Create empty file
  empty_df <- data.frame(human_gene = character(0), direction = character(0),
                         logFC = numeric(0), P.Value = numeric(0))
  write.csv(empty_df, "/home/chang/arc/MIRI_fda_xes/geo/GSE137293_human_l1000.csv",
            row.names = FALSE)
}

cat("\nOrtholog mapping complete!\n")
