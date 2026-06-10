#!/usr/bin/env Rscript
#
# Map mouse/rat DEGs to human L1000 landmark genes via orthologs
# Uses homologene for mouse/rat-to-human conversion

suppressPackageStartupMessages({
  library(org.Mm.eg.db)
  library(org.Rn.eg.db)
  library(org.Hs.eg.db)
  library(homologene)
  library(annotate)
})

# Load L1000 landmark genes
l1000_csv <- "/home/chang/arc/xpert/metadata/l1000_gene_info_978.csv"
l1000 <- read.csv(l1000_csv, stringsAsFactors = FALSE)
landmark_human <- unique(l1000$gene_symbol)  # column 5
cat(sprintf("L1000 landmark genes: %d\n", length(landmark_human)))

# ============================================================
# Helper: map mouse genes to human gene symbols
# ============================================================
mouse2human_symbol <- function(mouse_genes) {
  # Use homologene: mouse symbol -> human symbol
  human_symbols <- mouse2human(mouse_genes)
  return(human_symbols)
}

# ============================================================
# Helper: map rat genes to human gene symbols
# homologene package doesn't have rat2human - use Entrez-based mapping
# ============================================================
rat2human_symbol <- function(rat_genes) {
  # Step 1: rat gene symbol -> rat Entrez ID
  rat_entrez_all <- mapIds(org.Rn.eg.db, keys=rat_genes, column="ENTREZID",
                           keytype="SYMBOL", multiVals="first")

  # Remove NAs
  valid <- !is.na(rat_entrez_all)
  if (sum(valid) == 0) {
    return(setNames(rep(NA_character_, length(rat_genes)), rat_genes))
  }

  rat_e_valid <- rat_entrez_all[valid]
  rat_genes_valid <- names(rat_e_valid)

  # Step 2: Use homologene to get human homologs from rat Entrez
  tryCatch({
    he <- homologene(as.integer(rat_e_valid), inTax=10116, outTax=9606)
  }, error = function(e) {
    cat("homologene error:", conditionMessage(e), "\n")
    # Fallback: use Ensembl or case-insensitive matching
    return(setNames(rep(NA_character_, length(rat_genes)), rat_genes))
  })

  if (nrow(he) == 0) {
    return(setNames(rep(NA_character_, length(rat_genes)), rat_genes))
  }

  # Step 3: human Entrez -> human gene symbol
  human_entrez_ids <- as.character(he[, 2])
  human_symbols <- mapIds(org.Hs.eg.db,
                          keys=human_entrez_ids,
                          column="SYMBOL", keytype="ENTREZID",
                          multiVals="first")

  # Step 4: Build mapping from original rat genes
  result <- setNames(rep(NA_character_, length(rat_genes)), rat_genes)

  # Match rat Entrez back to original genes
  rat_e_in_map <- as.integer(he[, 1])
  for (i in seq_len(nrow(he))) {
    rat_e <- rat_e_in_map[i]
    idx <- which(rat_e_valid == rat_e)
    if (length(idx) > 0) {
      result[rat_genes_valid[idx[1]]] <- human_symbols[i]
    }
  }

  return(result)
}

# ============================================================
# Process GSE160516 (Mouse)
# ============================================================
d1 <- read.csv("/home/chang/arc/drug_ranking/xpert/miri_efficacy/geo/GSE160516_deg_results.csv")

# Filter: P < 0.01, |logFC| > 0.2 (broader threshold for more genes)
g1 <- d1[d1$P.Value < 0.01 & abs(d1$logFC) > 0.2, ]
cat(sprintf("\nGSE160516 DEGs: %d total filtered\n", nrow(g1)))

g1_up <- g1[g1$logFC > 0, ]
g1_down <- g1[g1$logFC < 0, ]
cat(sprintf("  Up: %d, Down: %d\n", nrow(g1_up), nrow(g1_down)))

# Map to human (homologene may return data.frame or named vector)
cat("Mapping mouse up genes to human...\n")
mouse_up_genes <- g1_up$Gene
human_up_raw <- mouse2human_symbol(mouse_up_genes)
# Normalize to named character vector
if (is.data.frame(human_up_raw)) {
  human_up_vec <- setNames(human_up_raw[, 2], human_up_raw[, 1])
} else {
  human_up_vec <- as.character(human_up_raw)
}
cat("Mapping mouse down genes to human...\n")
mouse_down_genes <- g1_down$Gene
human_down_raw <- mouse2human_symbol(mouse_down_genes)
if (is.data.frame(human_down_raw)) {
  human_down_vec <- setNames(human_down_raw[, 2], human_down_raw[, 1])
} else {
  human_down_vec <- as.character(human_down_raw)
}

# Build data frame with proper alignment using named vector lookup
up_map <- data.frame(
  mouse_gene = mouse_up_genes,
  human_gene = human_up_vec[mouse_up_genes],
  logFC = g1_up$logFC,
  P.Value = g1_up$P.Value,
  stringsAsFactors = FALSE
)
down_map <- data.frame(
  mouse_gene = mouse_down_genes,
  human_gene = human_down_vec[mouse_down_genes],
  logFC = g1_down$logFC,
  P.Value = g1_down$P.Value,
  stringsAsFactors = FALSE
)

cat(sprintf("Mouse up mapped: %d/%d (%.1f%%)\n",
            sum(!is.na(up_map$human_gene)), nrow(up_map),
            100*sum(!is.na(up_map$human_gene))/nrow(up_map)))
cat(sprintf("Mouse down mapped: %d/%d (%.1f%%)\n",
            sum(!is.na(down_map$human_gene)), nrow(down_map),
            100*sum(!is.na(down_map$human_gene))/nrow(down_map)))

# Filter to L1000
up_map_lm <- up_map[!is.na(up_map$human_gene) & up_map$human_gene %in% landmark_human, ]
down_map_lm <- down_map[!is.na(down_map$human_gene) & down_map$human_gene %in% landmark_human, ]

cat(sprintf("Mouse up in L1000: %d\n", nrow(up_map_lm)))
cat(sprintf("Mouse down in L1000: %d\n", nrow(down_map_lm)))

# Remove duplicates (keep highest |logFC|)
up_map_lm <- up_map_lm[order(abs(up_map_lm$logFC), decreasing = TRUE), ]
up_map_lm <- up_map_lm[!duplicated(up_map_lm$human_gene), ]
down_map_lm <- down_map_lm[order(abs(down_map_lm$logFC), decreasing = TRUE), ]
down_map_lm <- down_map_lm[!duplicated(down_map_lm$human_gene), ]

set1_up <- unique(up_map_lm$human_gene)
set1_down <- unique(down_map_lm$human_gene)
cat(sprintf("GSE160516 L1000 unique: up=%d, down=%d\n", length(set1_up), length(set1_down)))

# Save
write.csv(up_map_lm, "/home/chang/arc/drug_ranking/xpert/miri_efficacy/geo/GSE160516_up_l1000.csv", row.names=FALSE)
write.csv(down_map_lm, "/home/chang/arc/drug_ranking/xpert/miri_efficacy/geo/GSE160516_down_l1000.csv", row.names=FALSE)

# ============================================================
# Process GSE108940 (Rat)
# ============================================================
d2 <- read.csv("/home/chang/arc/drug_ranking/xpert/miri_efficacy/geo/GSE108940_WT_deg_results.csv")

# Filter: P < 0.05, |logFC| > 0 (broadest)
g2 <- d2[d2$P.Value < 0.05 & abs(d2$logFC) > 0, ]
cat(sprintf("\nGSE108940 DEGs: %d total filtered\n", nrow(g2)))

g2_up <- g2[g2$logFC > 0, ]
g2_down <- g2[g2$logFC < 0, ]
cat(sprintf("  Up: %d, Down: %d\n", nrow(g2_up), nrow(g2_down)))

# Map to human
cat("Mapping rat up genes to human...\n")
rat_up_genes <- g2_up$Gene
human_up2_map <- rat2human_symbol(rat_up_genes)
cat("Mapping rat down genes to human...\n")
rat_down_genes <- g2_down$Gene
human_down2_map <- rat2human_symbol(rat_down_genes)

# Normalize to named character vector
if (is.data.frame(human_up2_map)) {
  human_up2_vec <- setNames(human_up2_map[, 2], human_up2_map[, 1])
} else {
  human_up2_vec <- as.character(human_up2_map)
}
if (is.data.frame(human_down2_map)) {
  human_down2_vec <- setNames(human_down2_map[, 2], human_down2_map[, 1])
} else {
  human_down2_vec <- as.character(human_down2_map)
}

up_map2 <- data.frame(
  rat_gene = rat_up_genes,
  human_gene = human_up2_vec[rat_up_genes],
  logFC = g2_up$logFC,
  P.Value = g2_up$P.Value,
  stringsAsFactors = FALSE
)
down_map2 <- data.frame(
  rat_gene = rat_down_genes,
  human_gene = human_down2_vec[rat_down_genes],
  logFC = g2_down$logFC,
  P.Value = g2_down$P.Value,
  stringsAsFactors = FALSE
)

cat(sprintf("Rat up mapped: %d/%d (%.1f%%)\n",
            sum(!is.na(up_map2$human_gene)), nrow(up_map2),
            100*sum(!is.na(up_map2$human_gene))/nrow(up_map2)))
cat(sprintf("Rat down mapped: %d/%d (%.1f%%)\n",
            sum(!is.na(down_map2$human_gene)), nrow(down_map2),
            100*sum(!is.na(down_map2$human_gene))/nrow(down_map2)))

# Filter to L1000 (remove any LMNTX prefix from homologene)
up_map2$human_gene <- gsub("^LMNTX", "", up_map2$human_gene)
down_map2$human_gene <- gsub("^LMNTX", "", down_map2$human_gene)

up_map2_lm <- up_map2[!is.na(up_map2$human_gene) & up_map2$human_gene %in% landmark_human, ]
down_map2_lm <- down_map2[!is.na(down_map2$human_gene) & down_map2$human_gene %in% landmark_human, ]

cat(sprintf("Rat up in L1000: %d\n", nrow(up_map2_lm)))
cat(sprintf("Rat down in L1000: %d\n", nrow(down_map2_lm)))

# Remove duplicates
up_map2_lm <- up_map2_lm[order(abs(up_map2_lm$logFC), decreasing = TRUE), ]
up_map2_lm <- up_map2_lm[!duplicated(up_map2_lm$human_gene), ]
down_map2_lm <- down_map2_lm[order(abs(down_map2_lm$logFC), decreasing = TRUE), ]
down_map2_lm <- down_map2_lm[!duplicated(down_map2_lm$human_gene), ]

set2_up <- unique(up_map2_lm$human_gene)
set2_down <- unique(down_map2_lm$human_gene)
cat(sprintf("GSE108940 L1000 unique: up=%d, down=%d\n", length(set2_up), length(set2_down)))

write.csv(up_map2_lm, "/home/chang/arc/drug_ranking/xpert/miri_efficacy/geo/GSE108940_up_l1000.csv", row.names=FALSE)
write.csv(down_map2_lm, "/home/chang/arc/drug_ranking/xpert/miri_efficacy/geo/GSE108940_down_l1000.csv", row.names=FALSE)

# ============================================================
# Cross-dataset consensus
# ============================================================
cons_up <- intersect(set1_up, set2_up)
cons_down <- intersect(set1_down, set2_down)

# Check conflicts
conf_up_dn <- intersect(set1_up, set2_down)
conf_dn_up <- intersect(set1_down, set2_up)

cat(sprintf("\n=== 2-DATASET CONSENSUS ===\n"))
cat(sprintf("Consensus UP: %d\n", length(cons_up)))
cat(sprintf("Consensus DOWN: %d\n", length(cons_down)))
cat(sprintf("TOTAL: %d\n", length(cons_up) + length(cons_down)))
cat(sprintf("\nConflicts (DS1 up vs DS2 down): %d\n", length(conf_up_dn)))
cat(sprintf("Conflicts (DS1 down vs DS2 up): %d\n", length(conf_dn_up)))

cat("\nConsensus UP genes:\n")
cat(paste(sort(cons_up), collapse = ", "), "\n")
cat("\nConsensus DOWN genes:\n")
cat(paste(sort(cons_down), collapse = ", "), "\n")

if (length(conf_up_dn) > 0) {
  cat("\nConflict UP(vsDOWN):\n")
  cat(paste(sort(conf_up_dn), collapse = ", "), "\n")
}
if (length(conf_dn_up) > 0) {
  cat("\nConflict DOWN(vsUP):\n")
  cat(paste(sort(conf_dn_up), collapse = ", "), "\n")
}

# Save consensus results
writeLines(sort(cons_up), "/home/chang/arc/drug_ranking/xpert/miri_efficacy/geo/consensus_up.txt")
writeLines(sort(cons_down), "/home/chang/arc/drug_ranking/xpert/miri_efficacy/geo/consensus_down.txt")
writeLines(c(sort(set1_up), "", sort(set1_down)),
           "/home/chang/arc/drug_ranking/xpert/miri_efficacy/geo/GSE160516_l1000_genes.txt")
writeLines(c(sort(set2_up), "", sort(set2_down)),
           "/home/chang/arc/drug_ranking/xpert/miri_efficacy/geo/GSE108940_l1000_genes.txt")

cat("\nMapping complete!\n")
