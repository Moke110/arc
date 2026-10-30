#!/usr/bin/env Rscript
#
# Build per-GSE DEG list for GSE160516 (mouse MIRI, Affymetrix Clariom S).
#
# This is the mouse portion of the original fda_miri_vs/geo/run_deg_analysis.R,
# split into a per-GSE script. Steps:
#   1. Read CEL files from source/ (GEO acquisition)
#   2. Quantile-normalize PM intensities, collapse probes to genes (max probe)
#   3. limma DEG: pooled IR vs Sham
#   4. Map mouse -> human L1000 (shared helper map_to_human.R)
#   5. Write per-GSE DEG list (human L1000 CSV) + MIRI JSON to local/deg/
#
# Usage:
#   cd disease_data/GSE160516
#   Rscript script/build_deg.R
#
# Prerequisite R packages: oligo, limma, pd.clariom.s.mouse,
#   clariomsmousetranscriptcluster.db, homologene, org.Mm.eg.db, org.Hs.eg.db

suppressPackageStartupMessages({
  library(oligo)
  library(limma)
  library(pd.clariom.s.mouse)
  library(clariomsmousetranscriptcluster.db)
  library(jsonlite)
})

options(affyio.threads = 1)

SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1])))
source(file.path(SCRIPT_DIR, "map_to_human.R"))
# Item layout (ADR 0030): the shared helper defines REPO_ROOT / ITEM_DIR /
# SOURCE_DIR / LOCAL_DEG_DIR; the DEG result CSV is read from source/ and the
# derived human L1000 lists are written to local/deg/.
DEG_DIR <- LOCAL_DEG_DIR
DATA_DIR <- SOURCE_DIR
dir.create(DEG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("\n===== GSE160516 - Mouse MIRI (Affymetrix Clariom S) =====\n\n")

# --- Read CEL files ---
cel_files <- list.files(DATA_DIR, pattern = "\\.CEL\\.gz$", full.names = TRUE)
if (length(cel_files) == 0) {
  stop("No CEL.gz files found in ", DATA_DIR,
       " — download the GSE160516 CEL files first (GEO acquisition).")
}
cat(sprintf("Reading %d CEL files...\n", length(cel_files)))
raw <- read.celfiles(cel_files)

pm_intens <- pm(raw)
cat(sprintf("PM matrix: %d probes x %d samples\n", nrow(pm_intens), ncol(pm_intens)))

cat("Log2 transforming and quantile normalizing...\n")
expr_raw <- log2(pm_intens + 1)
expr_norm <- normalizeBetweenArrays(expr_raw, method = "quantile")

# --- Probe-to-gene mapping via pd.clariom.s.mouse ---
cat("Mapping probes to genes...\n")
con <- db(pd.clariom.s.mouse)
map_query <- "
  SELECT pm.fid, fs.man_fsetid
  FROM pmfeature pm
  JOIN featureSet fs ON pm.fsetid = fs.fsetid
  WHERE fs.type = 1
"
probe_map <- dbGetQuery(con, map_query)
tc_ids <- unique(probe_map$man_fsetid)
tc_to_gene <- mapIds(clariomsmousetranscriptcluster.db, keys = tc_ids,
                     column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
tc_to_gene_vec <- setNames(as.character(tc_to_gene), names(tc_to_gene))
probe_map$gene <- tc_to_gene_vec[probe_map$man_fsetid]
probe_map <- probe_map[!is.na(probe_map$gene), ]

pm_fids <- as.integer(rownames(expr_norm))
expr_idx <- match(probe_map$fid, pm_fids)
valid <- !is.na(expr_idx)
probe_map <- probe_map[valid, ]
expr_idx <- expr_idx[valid]

# Collapse to gene level (max probe per gene)
genes <- unique(probe_map$gene)
gene_expr <- matrix(NA_real_, nrow = length(genes), ncol = ncol(expr_norm))
rownames(gene_expr) <- genes
colnames(gene_expr) <- colnames(expr_norm)
for (i in seq_along(genes)) {
  g <- genes[i]
  idx <- expr_idx[probe_map$gene == g]
  if (length(idx) == 1) {
    gene_expr[i, ] <- expr_norm[idx, ]
  } else {
    means <- rowMeans(expr_norm[idx, , drop = FALSE])
    gene_expr[i, ] <- expr_norm[idx[which.max(means)], ]
  }
}
cat(sprintf("Gene-level matrix: %d genes x %d samples\n", nrow(gene_expr), ncol(gene_expr)))

# --- limma DEG: Sham vs IR ---
sample_names <- colnames(gene_expr)
group <- ifelse(grepl("Con", sample_names), "Sham", "IR")
cat(sprintf("  Sham: %d, IR: %d\n", sum(group == "Sham"), sum(group == "IR")))

design <- model.matrix(~ 0 + factor(group, levels = c("Sham", "IR")))
colnames(design) <- c("Sham", "IR")
fit <- lmFit(gene_expr, design)
contrast <- makeContrasts(IR_vs_Sham = IR - Sham, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
deg <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")
deg$Gene <- rownames(deg)
deg <- deg[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
cat(sprintf("DEG: %d genes, adj.P<0.05 & |logFC|>0.5: %d\n",
            nrow(deg), sum(deg$adj.P.Val < 0.05 & abs(deg$logFC) > 0.5)))

write.csv(deg, file.path(DATA_DIR, "GSE160516_deg_results.csv"), row.names = FALSE, quote = FALSE)

# --- Map mouse -> human L1000 ---
cat("\nMapping mouse DEGs to human L1000...\n")
mapped <- map_deg_to_l1000(deg, species = "mouse", p_threshold = 0.01, lfc_threshold = 0.2)

# Split up/down CSVs (historical format: mouse_gene, human_gene, logFC, P.Value)
up <- mapped[mapped$direction == "UP", c("species_gene", "human_gene", "logFC", "P.Value")]
down <- mapped[mapped$direction == "DOWN", c("species_gene", "human_gene", "logFC", "P.Value")]
names(up)[1] <- "mouse_gene"
names(down)[1] <- "mouse_gene"
write.csv(up, file.path(DEG_DIR, "GSE160516_up_l1000.csv"), row.names = FALSE)
write.csv(down, file.path(DEG_DIR, "GSE160516_down_l1000.csv"), row.names = FALSE)

# --- Write per-GSE MIRI JSON ---
up_genes <- sort(unique(mapped$human_gene[mapped$direction == "UP"]))
down_genes <- sort(unique(mapped$human_gene[mapped$direction == "DOWN"]))
gs <- list(
  name = "GSE160516_MIRI",
  species = "Mus musculus",
  source_geo = "GSE160516",
  platform = "Affymetrix Clariom S",
  comparison = "MIRI (pooled 6h/24h/72h IR) vs Sham",
  description = paste0("Mouse MIRI model: 45min LAD ischemia + reperfusion (6h/24h/72h). ",
                       "Pooled IR vs Sham DEGs (limma, P<0.01, |logFC|>0.2), ",
                       "mouse-to-human ortholog via homologene, intersected with L1000 978."),
  up_genes = up_genes,
  down_genes = down_genes,
  n_up_l1000 = length(up_genes),
  n_down_l1000 = length(down_genes)
)
writeLines(jsonlite::toJSON(gs, pretty = TRUE, auto_unbox = TRUE),
           file.path(DEG_DIR, "GSE160516_MIRI.json"))

cat(sprintf("\nGSE160516 done: %d UP + %d DOWN (human L1000)\n",
            length(up_genes), length(down_genes)))
