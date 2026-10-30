#!/usr/bin/env Rscript
#
# Build per-GSE DEG list for GSE108940 (rat MIRI, Agilent GE2-v5).
#
# This is the rat portion of the original fda_miri_vs/geo/run_deg_analysis.R,
# split into a per-GSE script. Steps:
#   1. Parse Agilent FEATURES files from source/ (GEO acquisition)
#   2. Log2 + quantile normalize, collapse probes to genes (max probe)
#   3. limma DEG: WT IR vs WT Sham
#   4. Map rat -> human L1000 (shared helper map_to_human.R)
#   5. Write per-GSE DEG list (human L1000 CSV) + MIRI JSON to local/deg/
#
# Usage:
#   cd disease_data/GSE108940
#   Rscript script/build_deg.R
#
# Prerequisite R packages: limma, homologene, org.Rn.eg.db, org.Hs.eg.db

suppressPackageStartupMessages({
  library(limma)
  library(jsonlite)
})

SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1])))
source(file.path(SCRIPT_DIR, "map_to_human.R"))
# Item layout (ADR 0030): the shared helper defines REPO_ROOT / ITEM_DIR /
# SOURCE_DIR / LOCAL_DEG_DIR; the DEG result CSV is read from source/ and the
# derived human L1000 lists are written to local/deg/.
DEG_DIR <- LOCAL_DEG_DIR
DATA_DIR <- SOURCE_DIR
dir.create(DEG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("\n===== GSE108940 - Rat MIRI (Agilent GE2-v5) =====\n\n")

# --- Parse Agilent files ---
agilent_files <- list.files(DATA_DIR, pattern = "\\.txt\\.gz$", full.names = TRUE)
if (length(agilent_files) == 0) {
  stop("No Agilent .txt.gz files found in ", DATA_DIR,
       " — download the GSE108940 Agilent files first (GEO acquisition).")
}
cat(sprintf("Found %d Agilent files\n", length(agilent_files)))

parse_agilent <- function(fp) {
  lines <- readLines(fp)
  feat_line <- which(grepl("^FEATURES", lines))[1]
  if (is.na(feat_line)) return(NULL)
  header <- unlist(strsplit(lines[feat_line], "\t"))[-1]
  data_lines <- lines[(feat_line + 1):length(lines)]
  blank <- which(data_lines == "" | (grepl("^[A-Z]", data_lines) & !grepl("^DATA", data_lines)))
  if (length(blank) > 0) data_lines <- data_lines[1:(blank[1] - 1)]
  data_lines <- sub("^DATA\t", "", data_lines)
  if (length(data_lines) == 0) return(NULL)
  con <- textConnection(paste(data_lines, collapse = "\n"))
  dat <- read.delim(con, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  close(con)
  colnames(dat) <- header
  pn_col <- which(header == "ProbeName")
  gn_col <- which(header == "GeneName")
  gs_col <- which(header == "gProcessedSignal")
  data.frame(
    ProbeName = as.character(dat[, pn_col]),
    GeneName = as.character(dat[, gn_col]),
    Signal = suppressWarnings(as.numeric(as.character(dat[, gs_col]))),
    stringsAsFactors = FALSE
  )
}

all_data <- list()
for (f in agilent_files) {
  bname <- sub("\\.txt\\.gz$", "", basename(f))
  parsed <- parse_agilent(f)
  if (!is.null(parsed)) all_data[[bname]] <- parsed
}
sample_names <- names(all_data)
cat(sprintf("Parsed %d files\n", length(all_data)))

all_probes <- unique(unlist(lapply(all_data, function(x) x$ProbeName)))
ag_expr <- matrix(NA_real_, nrow = length(all_probes), ncol = length(sample_names))
rownames(ag_expr) <- all_probes
colnames(ag_expr) <- sample_names
for (i in seq_along(sample_names)) {
  df <- all_data[[sample_names[i]]]
  df_u <- df[!duplicated(df$ProbeName), ]
  idx <- match(df_u$ProbeName, all_probes)
  ag_expr[idx, i] <- df_u$Signal
}

is_control <- grepl("^(GE_|DarkCorner|BrightCorner|dmr|NC|gps|Corner|ETG|E1A|Pro25G|eQC|AFFX|ERCC|RNA|Spike)",
                    rownames(ag_expr), ignore.case = TRUE)
ag_expr_filt <- ag_expr[!is_control, ]
ag_expr_filt[ag_expr_filt < 1] <- 1
ag_expr_log <- log2(ag_expr_filt)
ag_expr_norm <- normalizeBetweenArrays(ag_expr_log, method = "quantile")

gene_names_raw <- all_data[[1]][match(rownames(ag_expr_norm), all_data[[1]]$ProbeName), "GeneName"]
gene_names_raw[gene_names_raw == ""] <- NA
ag_genes <- gene_names_raw
valid_ag <- !is.na(ag_genes)
ag_expr_valid <- ag_expr_norm[valid_ag, ]
ag_genes_valid <- ag_genes[valid_ag]

unique_ag_genes <- unique(ag_genes_valid)
ag_gene_expr <- matrix(NA_real_, nrow = length(unique_ag_genes), ncol = ncol(ag_expr_valid))
rownames(ag_gene_expr) <- unique_ag_genes
colnames(ag_gene_expr) <- colnames(ag_expr_valid)
for (i in seq_along(unique_ag_genes)) {
  g <- unique_ag_genes[i]
  probes <- which(ag_genes_valid == g)
  if (length(probes) == 1) {
    ag_gene_expr[i, ] <- ag_expr_valid[probes, ]
  } else {
    means <- rowMeans(ag_expr_valid[probes, , drop = FALSE])
    ag_gene_expr[i, ] <- ag_expr_valid[probes[which.max(means)], ]
  }
}
cat(sprintf("Gene-level: %d genes x %d samples\n", nrow(ag_gene_expr), ncol(ag_gene_expr)))

# --- WT Sham vs WT IR ---
ag_group <- ifelse(grepl("WT_SHAM", colnames(ag_gene_expr)), "WT_Sham",
            ifelse(grepl("WT_IR", colnames(ag_gene_expr)), "WT_IR",
            ifelse(grepl("TG_SHAM", colnames(ag_gene_expr)), "TG_Sham", "TG_IR")))
wt_idx <- ag_group %in% c("WT_Sham", "WT_IR")
wt_expr <- ag_gene_expr[, wt_idx]
wt_grp <- factor(ag_group[wt_idx], levels = c("WT_Sham", "WT_IR"))

d_wt <- model.matrix(~ 0 + wt_grp)
colnames(d_wt) <- c("WT_Sham", "WT_IR")
f_wt <- lmFit(wt_expr, d_wt)
c_wt <- makeContrasts(IR_vs_Sham = WT_IR - WT_Sham, levels = d_wt)
f2_wt <- contrasts.fit(f_wt, c_wt)
f2_wt <- eBayes(f2_wt, trend = TRUE, robust = TRUE)
deg <- topTable(f2_wt, coef = 1, number = Inf, adjust.method = "BH")
deg$Gene <- rownames(deg)
deg <- deg[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
cat(sprintf("DEG: %d genes, adj.P<0.05 & |logFC|>0.5: %d\n",
            nrow(deg), sum(deg$adj.P.Val < 0.05 & abs(deg$logFC) > 0.5)))

write.csv(deg, file.path(DATA_DIR, "GSE108940_WT_deg_results.csv"), row.names = FALSE, quote = FALSE)

# --- Map rat -> human L1000 ---
cat("\nMapping rat DEGs to human L1000...\n")
mapped <- map_deg_to_l1000(deg, species = "rat", p_threshold = 0.05, lfc_threshold = 0.3)
write.csv(mapped, file.path(DEG_DIR, "GSE108940_human_l1000.csv"), row.names = FALSE)

up_genes <- sort(unique(mapped$human_gene[mapped$direction == "UP"]))
down_genes <- sort(unique(mapped$human_gene[mapped$direction == "DOWN"]))
gs <- list(
  name = "GSE108940_MIRI",
  species = "Rattus norvegicus",
  source_geo = "GSE108940",
  platform = "Agilent SurePrint G3 GE 8x60K",
  comparison = "WT MIRI (IR) vs WT Sham",
  description = paste0("Rat MIRI model: LAD ischemia-reperfusion. ",
                       "WT IR vs WT Sham DEGs (limma, P<0.05, |logFC|>0.3), ",
                       "rat-to-human ortholog via homologene, filtered to L1000 978."),
  up_genes = up_genes,
  down_genes = down_genes,
  n_up_l1000 = length(up_genes),
  n_down_l1000 = length(down_genes)
)
writeLines(jsonlite::toJSON(gs, pretty = TRUE, auto_unbox = TRUE),
           file.path(DEG_DIR, "GSE108940_MIRI.json"))

cat(sprintf("\nGSE108940 done: %d UP + %d DOWN (human L1000)\n",
            length(up_genes), length(down_genes)))
