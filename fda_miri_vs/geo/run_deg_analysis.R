#!/usr/bin/env Rscript
#
# Complete DEG analysis for MIRI reference gene set
# GSE160516: Mouse Clariom S microarray (4 sham vs 12 IR)
# GSE108940: Rat Agilent GE2-v5 (WT sham vs WT IR)
#
# Output: *_deg_results.csv files for downstream Python consensus

suppressPackageStartupMessages({
  library(oligo)
  library(limma)
  library(pd.clariom.s.mouse)
  library(clariomsmousetranscriptcluster.db)
})

# Force single-threaded operation
options(affyio.threads = 1)

cat("\n")
cat("########################################################################\n")
cat("#  DATASET 1: GSE160516 - Mouse MIRI (Affymetrix Clariom S)\n")
cat("########################################################################\n\n")

# --- Read CEL files ---
cel_files <- list.files("GSE160516", pattern = "\\.CEL\\.gz$", full.names = TRUE)
cat(sprintf("Reading %d CEL files...\n", length(cel_files)))
raw <- read.celfiles(cel_files)

# --- Get PM intensities ---
pm_intens <- pm(raw)
cat(sprintf("PM matrix: %d probes x %d samples\n", nrow(pm_intens), ncol(pm_intens)))

# --- Log2 transform + quantile normalize ---
cat("Log2 transforming and quantile normalizing...\n")
expr_raw <- log2(pm_intens + 1)
expr_norm <- normalizeBetweenArrays(expr_raw, method = "quantile")

# --- Probe-to-gene mapping via pd.clariom.s.mouse ---
cat("Mapping probes to genes...\n")
con <- db(pd.clariom.s.mouse)

# Get mapping: pm feature index -> fsetid -> man_fsetid (transcript cluster ID)
map_query <- "
  SELECT pm.fid, fs.man_fsetid
  FROM pmfeature pm
  JOIN featureSet fs ON pm.fsetid = fs.fsetid
  WHERE fs.type = 1
"
probe_map <- dbGetQuery(con, map_query)
cat(sprintf("Probe-to-TC mapping: %d probes -> %d unique TCs\n",
            nrow(probe_map), length(unique(probe_map$man_fsetid))))

# Map TC IDs to gene symbols
tc_ids <- unique(probe_map$man_fsetid)
cat(sprintf("Mapping %d TC IDs to gene symbols...\n", length(tc_ids)))
tc_to_gene <- mapIds(
  clariomsmousetranscriptcluster.db,
  keys = tc_ids,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)
cat(sprintf("Mapped: %d TCs have gene symbols\n", sum(!is.na(tc_to_gene))))

# Build probe_id -> gene_symbol lookup
tc_to_gene_vec <- setNames(as.character(tc_to_gene), names(tc_to_gene))
probe_map$gene <- tc_to_gene_vec[probe_map$man_fsetid]
probe_map <- probe_map[!is.na(probe_map$gene), ]
cat(sprintf("Usable probe-gene pairs: %d\n", nrow(probe_map)))

# --- Collapse probes to gene level ---
# Align probe index (fid) with PM matrix rows
pm_fids <- as.integer(rownames(expr_norm))
expr_idx <- match(probe_map$fid, pm_fids)
valid <- !is.na(expr_idx)
probe_map <- probe_map[valid, ]
expr_idx <- expr_idx[valid]
cat(sprintf("Probes in expression matrix: %d\n", length(expr_idx)))

# Build gene expression matrix (max probe per gene)
genes <- unique(probe_map$gene)
n_genes <- length(genes)
cat(sprintf("Unique genes: %d\n", n_genes))

gene_expr <- matrix(NA_real_, nrow = n_genes, ncol = ncol(expr_norm))
rownames(gene_expr) <- genes
colnames(gene_expr) <- colnames(expr_norm)

pb <- txtProgressBar(min = 0, max = n_genes, style = 3)
for (i in seq_len(n_genes)) {
  g <- genes[i]
  idx <- expr_idx[probe_map$gene == g]
  if (length(idx) == 1) {
    gene_expr[i, ] <- expr_norm[idx, ]
  } else {
    means <- rowMeans(expr_norm[idx, , drop = FALSE])
    best <- idx[which.max(means)]
    gene_expr[i, ] <- expr_norm[best, ]
  }
  if (i %% 5000 == 0) setTxtProgressBar(pb, i)
}
close(pb)

cat(sprintf("Gene-level matrix: %d genes x %d samples\n", nrow(gene_expr), ncol(gene_expr)))

# --- limma DEG: Sham vs IR ---
sample_names <- colnames(gene_expr)
group <- ifelse(grepl("Con", sample_names), "Sham", "IR")

cat("\nSample groups:\n")
cat(sprintf("  Sham: %d, IR: %d\n", sum(group == "Sham"), sum(group == "IR")))

design <- model.matrix(~ 0 + factor(group, levels = c("Sham", "IR")))
colnames(design) <- c("Sham", "IR")

fit <- lmFit(gene_expr, design)
contrast <- makeContrasts(IR_vs_Sham = IR - Sham, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

deg160516 <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")
deg160516$Gene <- rownames(deg160516)
deg160516 <- deg160516[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]

cat(sprintf("\nGSE160516 DEG summary:\n"))
cat(sprintf("  Total genes: %d\n", nrow(deg160516)))
cat(sprintf("  adj.P < 0.05: %d\n", sum(deg160516$adj.P.Val < 0.05)))
cat(sprintf("  adj.P < 0.05 & |logFC| > 0.5: %d\n",
            sum(deg160516$adj.P.Val < 0.05 & abs(deg160516$logFC) > 0.5)))
cat(sprintf("  adj.P < 0.05 & |logFC| > 0.3: %d\n",
            sum(deg160516$adj.P.Val < 0.05 & abs(deg160516$logFC) > 0.3)))
cat(sprintf("  Up-regulated (adj.P < 0.05 & logFC > 0.5): %d\n",
            sum(deg160516$adj.P.Val < 0.05 & deg160516$logFC > 0.5)))
cat(sprintf("  Down-regulated (adj.P < 0.05 & logFC < -0.5): %d\n",
            sum(deg160516$adj.P.Val < 0.05 & deg160516$logFC < -0.5)))

write.csv(deg160516, "GSE160516_deg_results.csv", row.names = FALSE, quote = FALSE)
cat("\nSaved: GSE160516_deg_results.csv\n")

# --- Time-point specific analyses ---
for (tp in c("6h", "24h", "72h")) {
  tp_samp <- colnames(gene_expr)[grepl(paste0("IR", tp), colnames(gene_expr))]
  sham_samp <- colnames(gene_expr)[grepl("Con", colnames(gene_expr))]
  sub_expr <- gene_expr[, c(sham_samp, tp_samp)]
  sub_group <- factor(c(rep("Sham", length(sham_samp)), rep("IR", length(tp_samp))),
                      levels = c("Sham", "IR"))

  d_tp <- model.matrix(~ 0 + sub_group)
  colnames(d_tp) <- c("Sham", "IR")
  f_tp <- lmFit(sub_expr, d_tp)
  c_tp <- makeContrasts(IR_vs_Sham = IR - Sham, levels = d_tp)
  f2_tp <- contrasts.fit(f_tp, c_tp)
  f2_tp <- eBayes(f2_tp, trend = TRUE, robust = TRUE)

  res_tp <- topTable(f2_tp, coef = 1, number = Inf, adjust.method = "BH")
  write.csv(data.frame(Gene = rownames(res_tp), res_tp[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]),
            sprintf("GSE160516_deg_%s.csv", tp),
            row.names = FALSE, quote = FALSE)
  cat(sprintf("  %s: adj.P<0.05 & |logFC|>0.5 = %d (up=%d, down=%d)\n", tp,
              sum(res_tp$adj.P.Val < 0.05 & abs(res_tp$logFC) > 0.5),
              sum(res_tp$adj.P.Val < 0.05 & res_tp$logFC > 0.5),
              sum(res_tp$adj.P.Val < 0.05 & res_tp$logFC < -0.5)))
}

########################################################################
# DATASET 2: GSE108940
########################################################################
cat("\n")
cat("########################################################################\n")
cat("#  DATASET 2: GSE108940 - Rat MIRI (Agilent GE2-v5)\n")
cat("########################################################################\n\n")

agilent_files <- list.files("GSE108940", pattern = "\\.txt\\.gz$", full.names = TRUE)
cat(sprintf("Found %d Agilent files\n", length(agilent_files)))

# Parse Agilent files
parse_agilent <- function(fp) {
  lines <- readLines(fp)
  # Find FEATURES section (header is ON the FEATURES line itself)
  feat_line <- which(grepl("^FEATURES", lines))[1]
  if (is.na(feat_line)) return(NULL)

  # Parse header from the FEATURES line
  header <- unlist(strsplit(lines[feat_line], "\t"))
  # Remove "FEATURES" prefix (first element), keep column names
  header <- header[-1]

  # Read data lines (starting from line after FEATURES header)
  data_lines <- lines[(feat_line + 1):length(lines)]
  # Find end of data: blank line, or a section header that isn't "DATA"
  blank <- which(data_lines == "" | (grepl("^[A-Z]", data_lines) & !grepl("^DATA", data_lines)))
  if (length(blank) > 0) data_lines <- data_lines[1:(blank[1] - 1)]
  # Remove the "DATA" prefix from each line (Agilent format has "DATA\t" at start)
  data_lines <- sub("^DATA\t", "", data_lines)

  if (length(data_lines) == 0) return(NULL)

  # Parse tab-separated data
  con <- textConnection(paste(data_lines, collapse = "\n"))
  dat <- read.delim(con, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  close(con)
  colnames(dat) <- header

  # Find relevant columns
  pn_col <- which(header == "ProbeName")
  gn_col <- which(header == "GeneName")
  gs_col <- which(header == "gProcessedSignal")

  result <- data.frame(
    ProbeName = as.character(dat[, pn_col]),
    GeneName = as.character(dat[, gn_col]),
    Signal = suppressWarnings(as.numeric(as.character(dat[, gs_col]))),
    stringsAsFactors = FALSE
  )
  return(result)
}

# Parse all files
all_data <- list()
for (f in agilent_files) {
  bname <- sub("\\.txt\\.gz$", "", basename(f))
  parsed <- parse_agilent(f)
  if (!is.null(parsed)) {
    all_data[[bname]] <- parsed
  }
}
sample_names_ag <- names(all_data)
cat(sprintf("Parsed %d files\n", length(all_data)))

# Build expression matrix
all_probes <- unique(unlist(lapply(all_data, function(x) x$ProbeName)))
cat(sprintf("Total unique probes: %d\n", length(all_probes)))

ag_expr <- matrix(NA_real_, nrow = length(all_probes), ncol = length(sample_names_ag))
rownames(ag_expr) <- all_probes
colnames(ag_expr) <- sample_names_ag

for (i in seq_along(sample_names_ag)) {
  df <- all_data[[sample_names_ag[i]]]
  df_u <- df[!duplicated(df$ProbeName), ]
  idx <- match(df_u$ProbeName, all_probes)
  ag_expr[idx, i] <- df_u$Signal
}

# Filter control probes and log2 transform
is_control <- grepl("^(GE_|DarkCorner|BrightCorner|dmr|NC|gps|Corner|ETG|E1A|Pro25G|eQC|AFFX|ERCC|RNA|Spike)",
                    rownames(ag_expr), ignore.case = TRUE)
cat(sprintf("Control probes: %d\n", sum(is_control)))
ag_expr_filt <- ag_expr[!is_control, ]

# Log2 transform (Agilent signals can be negative; set floor)
ag_expr_filt[ag_expr_filt < 1] <- 1
ag_expr_log <- log2(ag_expr_filt)

# Quantile normalize
ag_expr_norm <- normalizeBetweenArrays(ag_expr_log, method = "quantile")

# Get gene names
gene_names_raw <- all_data[[1]][match(rownames(ag_expr_norm), all_data[[1]]$ProbeName), "GeneName"]
gene_names_raw[gene_names_raw == ""] <- NA
ag_genes <- gene_names_raw

# Remove probes without gene symbol
valid_ag <- !is.na(ag_genes)
ag_expr_valid <- ag_expr_norm[valid_ag, ]
ag_genes_valid <- ag_genes[valid_ag]
cat(sprintf("Probes with gene symbols: %d\n", sum(valid_ag)))

# Collapse to gene level
unique_ag_genes <- unique(ag_genes_valid)
n_ag_genes <- length(unique_ag_genes)
ag_gene_expr <- matrix(NA_real_, nrow = n_ag_genes, ncol = ncol(ag_expr_valid))
rownames(ag_gene_expr) <- unique_ag_genes
colnames(ag_gene_expr) <- colnames(ag_expr_valid)

for (i in seq_len(n_ag_genes)) {
  g <- unique_ag_genes[i]
  probes <- which(ag_genes_valid == g)
  if (length(probes) == 1) {
    ag_gene_expr[i, ] <- ag_expr_valid[probes, ]
  } else {
    means <- rowMeans(ag_expr_valid[probes, , drop = FALSE])
    best <- probes[which.max(means)]
    ag_gene_expr[i, ] <- ag_expr_valid[best, ]
  }
}
cat(sprintf("Gene-level: %d genes x %d samples\n", nrow(ag_gene_expr), ncol(ag_gene_expr)))

# Determine groups from sample names
ag_group <- ifelse(grepl("WT_SHAM", colnames(ag_gene_expr)), "WT_Sham",
            ifelse(grepl("WT_IR", colnames(ag_gene_expr)), "WT_IR",
            ifelse(grepl("TG_SHAM", colnames(ag_gene_expr)), "TG_Sham", "TG_IR")))

cat("\nAgilent sample groups:\n")
print(table(ag_group))

# WT Sham vs WT IR DEG
wt_idx <- ag_group %in% c("WT_Sham", "WT_IR")
wt_expr <- ag_gene_expr[, wt_idx]
wt_grp <- factor(ag_group[wt_idx], levels = c("WT_Sham", "WT_IR"))

cat(sprintf("\nWT samples for DEG: %d\n", sum(wt_idx)))

d_wt <- model.matrix(~ 0 + wt_grp)
colnames(d_wt) <- c("WT_Sham", "WT_IR")
f_wt <- lmFit(wt_expr, d_wt)
c_wt <- makeContrasts(IR_vs_Sham = WT_IR - WT_Sham, levels = d_wt)
f2_wt <- contrasts.fit(f_wt, c_wt)
f2_wt <- eBayes(f2_wt, trend = TRUE, robust = TRUE)

deg108940 <- topTable(f2_wt, coef = 1, number = Inf, adjust.method = "BH")
deg108940$Gene <- rownames(deg108940)
deg108940 <- deg108940[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]

cat(sprintf("\nGSE108940 WT DEG summary:\n"))
cat(sprintf("  Total genes: %d\n", nrow(deg108940)))
cat(sprintf("  adj.P < 0.05: %d\n", sum(deg108940$adj.P.Val < 0.05)))
cat(sprintf("  adj.P < 0.05 & |logFC| > 0.5: %d\n",
            sum(deg108940$adj.P.Val < 0.05 & abs(deg108940$logFC) > 0.5)))
cat(sprintf("  adj.P < 0.05 & |logFC| > 0.3: %d\n",
            sum(deg108940$adj.P.Val < 0.05 & abs(deg108940$logFC) > 0.3)))
cat(sprintf("  Up-regulated (adj.P < 0.05 & logFC > 0.5): %d\n",
            sum(deg108940$adj.P.Val < 0.05 & deg108940$logFC > 0.5)))
cat(sprintf("  Down-regulated (adj.P < 0.05 & logFC < -0.5): %d\n",
            sum(deg108940$adj.P.Val < 0.05 & deg108940$logFC < -0.5)))

write.csv(deg108940, "GSE108940_WT_deg_results.csv", row.names = FALSE, quote = FALSE)
cat("\nSaved: GSE108940_WT_deg_results.csv\n")

cat("\nDEG analysis complete!\n")
