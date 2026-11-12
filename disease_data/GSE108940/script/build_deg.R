#!/usr/bin/env Rscript
#
# Build per-GSE DEG list for GSE108940 (mouse MIRI, GPL7202 Agilent Whole Mouse
# Genome 4x44K; NCBI GEO taxon Mus musculus, 24 samples).
#
# This is the GSE108940 portion of the original fda_miri_vs/geo/run_deg_analysis.R,
# split into a per-GSE script. Steps:
#   1. Parse Agilent FEATURES files from source/ (GEO acquisition: the GSE108940
#      RAW tar, i.e. 24 GSM*_..._GE2-v5_*.txt.gz files)
#   2. Log2 + quantile normalize, collapse probes to genes (max probe)
#   3. limma DEG: WT IR vs WT Sham (12 wild-type samples: 6 IR + 6 sham)
#   4. Map mouse -> human L1000 (shared helper map_to_human.R, homologene
#      inTax=10090, outTax=9606)
#   5. Write per-GSE DEG list (human L1000 CSV) + MIRI JSON to local/deg/
#
# The raw FEATURES files are NOT redistributed: source/ keeps only the limma DEG
# table (GSE108940_WT_deg_results.csv). When source/ has no .txt.gz files this
# script therefore takes the DEG table as given and runs steps 4-5 only; it
# cannot rebuild the DEG table itself without re-downloading GSE108940_RAW.tar
# (380,979,200 bytes) and re-running steps 1-3. The map-only path is the
# reproducible way to regenerate local/deg/ from the tracked-source DEG table.
#
# Usage:
#   cd disease_data/GSE108940
#   Rscript script/build_deg.R
#
# Prerequisite R packages: limma, jsonlite, homologene, org.Mm.eg.db, org.Hs.eg.db

suppressPackageStartupMessages({
  library(limma)
  library(jsonlite)
  library(org.Mm.eg.db)
})

SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1])))
source(file.path(SCRIPT_DIR, "map_to_human.R"))
# Item layout (ADR 0030): the shared helper defines REPO_ROOT / ITEM_DIR /
# SOURCE_DIR / LOCAL_DEG_DIR; the DEG result CSV is read from source/ and the
# derived human L1000 lists are written to local/deg/.
DEG_DIR <- LOCAL_DEG_DIR
DATA_DIR <- SOURCE_DIR
DEG_CSV <- file.path(DATA_DIR, "GSE108940_WT_deg_results.csv")
dir.create(DEG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("\n===== GSE108940 - Mouse MIRI (GPL7202 Agilent Whole Mouse Genome 4x44K) =====\n\n")

# --- DEG table: rebuild from raw FEATURES if present, else reuse the local one ---
agilent_files <- list.files(DATA_DIR, pattern = "\\.txt\\.gz$", full.names = TRUE)

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

if (length(agilent_files) > 0) {
  cat(sprintf("Found %d Agilent FEATURES files — rebuilding the DEG table from raw data\n",
              length(agilent_files)))

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

  write.csv(deg, DEG_CSV, row.names = FALSE, quote = FALSE)
} else {
  if (!file.exists(DEG_CSV)) {
    stop("No Agilent .txt.gz files and no ", basename(DEG_CSV), " in ", DATA_DIR,
         " — download GSE108940_RAW.tar (380,979,200 bytes) and extract it first.")
  }
  cat("No raw Agilent FEATURES files in source/ — mapping the existing DEG table only\n")
  cat(sprintf("(cannot rebuild the DEG table without GSE108940_RAW.tar; using %s)\n",
              basename(DEG_CSV)))
  deg <- read.csv(DEG_CSV, stringsAsFactors = FALSE)
}

# --- Map mouse -> human L1000 ---
cat("\nMapping mouse DEGs to human L1000...\n")
mapped <- map_deg_to_l1000(deg, species = "mouse", p_threshold = 0.05, lfc_threshold = 0.3)
write.csv(mapped, file.path(DEG_DIR, "GSE108940_human_l1000.csv"), row.names = FALSE)

up_genes <- sort(unique(mapped$human_gene[mapped$direction == "UP"]))
down_genes <- sort(unique(mapped$human_gene[mapped$direction == "DOWN"]))
gs <- list(
  name = "GSE108940_MIRI",
  species = "Mus musculus",
  source_geo = "GSE108940",
  platform = "Agilent-014868 Whole Mouse Genome Microarray 4x44K G4122F (GPL7202)",
  comparison = "WT MIRI (IR) vs WT Sham",
  description = paste0("Mouse MIRI model: LAD ischemia-reperfusion (MuRF1 WT + I/R vs sham). ",
                       "WT IR vs WT Sham DEGs (limma, P<0.05, |logFC|>0.3), ",
                       "mouse-to-human ortholog via homologene R package (inTax=10090, outTax=9606), ",
                       "filtered to L1000 978 landmark genes."),
  up_genes = up_genes,
  down_genes = down_genes,
  n_up_l1000 = length(up_genes),
  n_down_l1000 = length(down_genes)
)
writeLines(jsonlite::toJSON(gs, pretty = TRUE, auto_unbox = TRUE),
           file.path(DEG_DIR, "GSE108940_MIRI.json"))

cat(sprintf("\nGSE108940 done: %d UP + %d DOWN (human L1000)\n",
            length(up_genes), length(down_genes)))
