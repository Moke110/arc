#!/usr/bin/env Rscript
#
# Build per-GSE DEG list for GSE137293 (pig scar vs border, Illumina NovaSeq 6000).
#
# The original pig DEG analysis was performed with a scipy Welch t-test on
# log2(CPM+1) counts (there is no limma step — the GEO series lacks a sham
# control). This per-GSE script reproduces that analysis and maps the results
# to human L1000:
#   1. Read raw counts from data/GSE137293_COUNTS.txt.gz (GEO acquisition)
#   2. log2(CPM+1), Welch t-test scar vs border (8 vs 8)
#   3. BH-adjusted p-values
#   4. Map pig -> human L1000 (shared helper map_to_human.R; drops ENSSCG IDs)
#   5. Write per-GSE DEG list (human L1000 CSV) + MIRI JSON to deg/
#
# Usage:
#   cd disease_data/miri/GSE137293
#   Rscript deg/build_deg.R
#
# Prerequisite R packages: scran (for .gz read) or base gzfile, homologene.

suppressPackageStartupMessages({
  library(jsonlite)
})

SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1])))
DEG_DIR <- SCRIPT_DIR
DATA_DIR <- normalizePath(file.path(SCRIPT_DIR, "..", "data"))
source(file.path(dirname(dirname(SCRIPT_DIR)), "scripts", "map_to_human.R"))

cat("\n===== GSE137293 - Pig scar vs border =====\n\n")

counts_path <- file.path(DATA_DIR, "GSE137293_COUNTS.txt.gz")
if (!file.exists(counts_path)) stop("Counts file not found: ", counts_path)

counts <- read.delim(gzfile(counts_path), stringsAsFactors = FALSE)
cat(sprintf("Counts: %d genes x %d samples\n", nrow(counts), ncol(counts) - 1))

samples <- colnames(counts)[-1]
scar <- samples[grepl("scar", samples, ignore.case = TRUE)]
border <- samples[grepl("border", samples, ignore.case = TRUE)]
cat(sprintf("Scar: %d, Border: %d\n", length(scar), length(border)))

# log2(CPM+1) normalization (matches original scipy pipeline)
count_mat <- as.matrix(counts[, -1])
totals <- colSums(count_mat)
totals[totals == 0] <- 1
cpm <- sweep(count_mat, 2, totals, "/") * 1e6
log2_cpm <- log2(cpm + 1)

# Welch t-test per gene (equivalent to scipy.stats.ttest_ind equal_var=False)
p_vals <- apply(log2_cpm, 1, function(x) {
  sx <- x[scar]
  bx <- x[border]
  t.test(sx, bx, alternative = "two.sided", var.equal = FALSE)$p.value
})
p_adj <- p.adjust(p_vals, method = "BH")

mean_scar <- rowMeans(log2_cpm[, scar, drop = FALSE])
mean_border <- rowMeans(log2_cpm[, border, drop = FALSE])
log2fc <- mean_scar - mean_border

deg <- data.frame(
  gene = counts[[1]],
  log2FC = log2fc,
  mean_scar = mean_scar,
  mean_border = mean_border,
  P_Value = p_vals,
  adj_P_Val = p_adj,
  stringsAsFactors = FALSE
)
deg <- deg[order(abs(deg$log2FC), decreasing = TRUE), ]
cat(sprintf("DEG: %d genes, adj.P<0.05 & |log2FC|>0.3: %d\n",
            nrow(deg), sum(deg$adj_P_Val < 0.05 & abs(deg$log2FC) > 0.3)))

write.csv(deg, file.path(DATA_DIR, "GSE137293_deg_results.csv"), row.names = FALSE)

# --- Map pig -> human L1000 ---
# Faithful to the original Python pipeline (generate_reference_genesets.py):
#   1. Direct gene-symbol matching for human-readable symbols (most match human)
#   2. Supplement with homologene results for the remaining genes
cat("\nMapping pig DEGs to human L1000...\n")
landmark <- load_landmark_genes()
l1000_upper <- toupper(landmark)

sig <- deg[!is.na(deg$P_Value) & deg$P_Value < 0.05 & abs(deg$log2FC) > 0.3, ]
if (nrow(sig) < 10) sig <- deg[!is.na(deg$P_Value) & deg$P_Value < 0.05, ]

mapped <- do.call(rbind, lapply(split(sig, sig$log2FC > 0), function(grp) {
  dir_ <- if (all(grp$log2FC > 0)) "UP" else "DOWN"
  readable <- grp$gene[!grepl("^ENS", grp$gene)]
  direct <- readable[toupper(readable) %in% l1000_upper]
  # homologene supplement for the rest
  homo <- species2human(setdiff(readable, direct), 9823)
  human <- c(direct, as.character(homo[!is.na(homo)]))
  human <- unique(human[human %in% landmark])
  if (length(human) == 0) return(NULL)
  data.frame(human_gene = human, direction = dir_, stringsAsFactors = FALSE)
}))
if (is.null(mapped)) mapped <- data.frame(human_gene = character(0), direction = character(0))
write.csv(mapped, file.path(DEG_DIR, "GSE137293_human_l1000.csv"), row.names = FALSE)

up_genes <- sort(unique(mapped$human_gene[mapped$direction == "UP"]))
down_genes <- sort(unique(mapped$human_gene[mapped$direction == "DOWN"]))
n_ens <- sum(grepl("^ENS", deg$gene))
gs <- list(
  name = "GSE137293_scar_vs_border",
  species = "Sus scrofa",
  source_geo = "GSE137293",
  platform = "Illumina NovaSeq 6000 (RNA-seq)",
  comparison = "Scar (infarct core) vs Border zone (peri-infarct), 4wk post-I/R",
  description = paste0("Pig I/R model: P30 pigs, 60min LAD occlusion + 4wk reperfusion. ",
                       "Scar vs Border zone comparison — CAVEAT: no sham control; ",
                       "comparison captures injury severity gradient, not disease vs healthy. ",
                       "DEG: Welch t-test on log2(CPM+1), P<0.05, |log2FC|>0.3. ",
                       "Pig-to-human: direct gene symbol matching + homologene (inTax=9823). ",
                       "Dropped ", n_ens, " ENSSCG genes without clear ortholog."),
  up_genes = up_genes,
  down_genes = down_genes,
  n_up_l1000 = length(up_genes),
  n_down_l1000 = length(down_genes),
  n_ensscg_dropped = n_ens,
  note = "Scar vs border — use with caution; not a MIRI-vs-sham comparison."
)
writeLines(jsonlite::toJSON(gs, pretty = TRUE, auto_unbox = TRUE),
           file.path(DEG_DIR, "GSE137293_scar_vs_border.json"))

cat(sprintf("\nGSE137293 done: %d UP + %d DOWN (human L1000)\n",
            length(up_genes), length(down_genes)))
