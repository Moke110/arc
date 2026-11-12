#!/usr/bin/env Rscript
#
# build_deg.R — full-transcriptome DESeq2 DEG for GSE128188
# (human atrial fibrillation vs sinus rhythm; paired left/right atrial
# appendage bulk RNA-seq).
#
# Input (source/, ADR 0030 — never modified):
#   source/GSE128188_CountsV2_EdgeR.txt.gz — Gene x 20 integer-count columns,
#   5 SR donors (7, 12b, 20, 27, 30) + 5 AF donors (5, 14, 25, 28, 31), each
#   donor contributing one LA and one RA appendage library.
#
# Design decision — the paired structure is POOLED, not blocked (stated here
# because the ticket asks which and why):
#   Disease status is a donor-level property: each donor is either AF or SR and
#   contributes both appendages, so every donor falls in exactly one group. A
#   donor fixed effect is therefore perfectly confounded with the donor-level
#   contrast (the donor indicators span the same subject-level space as the AF
#   indicator), which makes `~ donor + group` non-identifiable; DESeq2 has no
#   random-effect term either. The two remaining DESeq2 options are
#   (a) treat the 20 libraries as independent with a `~ tissue + group` design,
#       which pseudo-replicates each donor (10 donors counted twice) and
#       inflates type-I error, or
#   (b) pool the paired LA+RA libraries per donor into 10 independent
#       biological units and run the unpaired 5 AF vs 5 SR contrast.
#   We take (b). Counts are summed across each donor's two appendages; every
#   donor contributes one LA and one RA library, so tissue composition is
#   identical in both groups and the LA/RA baseline cancels out of the group
#   contrast, while no donor enters the model twice. `n_case`/`n_control` in
#   the CSV are therefore donors (5/5); the shipped expression matrix keeps all
#   20 GEO libraries (script/build_expression.R).
#
#   Group labels are not in the supplementary table; they are taken from the
#   GEO series matrix (!Sample_title: SR_LA_7, SR_RA_7, ... AF_LA_5, AF_RA_5,
#   ...). The donor lists are pinned below and validated against the counts
#   columns, so a silent label drift fails loudly.
#
# Outputs (both in local/, never in source/):
#   local/deg/GSE128188_deg.csv   full transcriptome, signed AF-vs-SR effect
#   local/deg/GSE128188_deg.json  small TRACKED summary + input/output digests
#
# Usage:
#   cd disease_data/GSE128188
#   /home/chang/miniconda3/envs/biomni_e1/bin/Rscript script/build_deg.R
#
# Prerequisite R packages: DESeq2, jsonlite, digest.

suppressPackageStartupMessages({
  library(DESeq2)
  library(jsonlite)
  library(digest)
})

SCRIPT_DIR <- dirname(normalizePath(
  sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1])
))
ITEM_DIR <- dirname(SCRIPT_DIR)
SOURCE_DIR <- file.path(ITEM_DIR, "source")
LOCAL_DEG_DIR <- file.path(ITEM_DIR, "local", "deg")
dir.create(LOCAL_DEG_DIR, recursive = TRUE, showWarnings = FALSE)

COUNTS_FILE <- file.path(SOURCE_DIR, "GSE128188_CountsV2_EdgeR.txt.gz")
DEG_CSV <- file.path(LOCAL_DEG_DIR, "GSE128188_deg.csv")
DEG_JSON <- file.path(LOCAL_DEG_DIR, "GSE128188_deg.json")

ACCESSION <- "GSE128188"
CONTRAST <- "AF_vs_SR"
N_CASE_EXPECTED <- 5L
N_CONTROL_EXPECTED <- 5L
THRESH_PADJ <- 0.05
THRESH_LFC <- 0.5

# GEO series-matrix sample titles give the group of each donor. The counts
# table column names are `<donor>-<tissue>`; both appendages of every donor are
# present.
SR_DONORS <- c("7", "12b", "20", "27", "30")
AF_DONORS <- c("5", "14", "25", "28", "31")

sha256 <- function(path) digest::digest(file = path, algo = "sha256")

cat("\n===== GSE128188 - human AF vs SR (paired LA+RA atrial appendage) =====\n\n")

if (!file.exists(COUNTS_FILE)) {
  stop("missing source/", basename(COUNTS_FILE),
       " — run 'bash script/fetch.sh' first")
}

# --- Read counts -----------------------------------------------------------
raw <- read.delim(gzfile(COUNTS_FILE), check.names = FALSE,
                  stringsAsFactors = FALSE)
counts <- as.matrix(raw[, -1, drop = FALSE])
storage.mode(counts) <- "integer"
rownames(counts) <- raw[[1]]
stopifnot(!anyDuplicated(rownames(counts)))
cat(sprintf("Counts: %d genes x %d libraries\n", nrow(counts), ncol(counts)))

# --- Sample metadata -------------------------------------------------------
samples <- colnames(counts)
donor <- sub("-[^-]+$", "", samples)
tissue <- sub("^.*-", "", samples)
if (!setequal(tissue, c("LA", "RA"))) {
  stop("unexpected tissue labels in counts columns: ",
       paste(sort(unique(tissue)), collapse = ", "))
}
if (!setequal(unique(donor), c(SR_DONORS, AF_DONORS))) {
  stop("donor IDs in the counts table do not match the GEO series-matrix donors")
}
group <- ifelse(donor %in% SR_DONORS, "SR", "AF")
meta <- data.frame(sample = samples, donor = donor, tissue = tissue,
                   group = group, stringsAsFactors = FALSE)
if (any(table(meta$donor) != 2L)) {
  stop("every donor must contribute exactly one LA and one RA library")
}
cat(sprintf("Libraries: SR %d, AF %d; donors: %d\n",
            sum(meta$group == "SR"), sum(meta$group == "AF"),
            length(unique(meta$donor))))

# --- Pool the paired LA+RA libraries per donor -----------------------------
donors <- sort(unique(meta$donor))
pooled <- vapply(donors, function(d) {
  rowSums(counts[, meta$donor == d, drop = FALSE])
}, numeric(nrow(counts)))
colnames(pooled) <- donors
donor_group <- ifelse(donors %in% SR_DONORS, "SR", "AF")
coldata <- data.frame(
  group = factor(donor_group, levels = c("SR", "AF")),
  row.names = donors
)
stopifnot(sum(coldata$group == "SR") == N_CONTROL_EXPECTED,
          sum(coldata$group == "AF") == N_CASE_EXPECTED)
cat(sprintf("Donor-pooled matrix: %d genes x %d donors (SR %d, AF %d)\n",
            nrow(pooled), ncol(pooled),
            sum(coldata$group == "SR"), sum(coldata$group == "AF")))

# --- DESeq2: AF vs SR ------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = round(pooled), colData = coldata,
                              design = ~ group)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("group", "AF", "SR"))
cat(sprintf("DESeq2: %d genes tested; %d with finite padj\n",
            nrow(res), sum(!is.na(res$padj))))

# --- Full-transcriptome DEG table -----------------------------------------
deg <- data.frame(
  gene_symbol = rownames(res),
  contrast = CONTRAST,
  log2FC = res$log2FoldChange,
  stat = res$stat,
  p_value = res$pvalue,
  padj = res$padj,
  baseMean = res$baseMean,
  n_case = as.integer(sum(coldata$group == "AF")),
  n_control = as.integer(sum(coldata$group == "SR")),
  stringsAsFactors = FALSE
)
deg <- deg[order(deg$gene_symbol, method = "radix"), , drop = FALSE]
rownames(deg) <- NULL
write.csv(deg, DEG_CSV, row.names = FALSE, quote = FALSE, na = "NA")

sig <- !is.na(deg$padj) & deg$padj < THRESH_PADJ &
  abs(deg$log2FC) >= THRESH_LFC
n_up <- sum(sig & deg$log2FC > 0, na.rm = TRUE)
n_down <- sum(sig & deg$log2FC < 0, na.rm = TRUE)
cat(sprintf("DEG calls (padj<%.2f & |log2FC|>=%.1f): %d up / %d down\n",
            THRESH_PADJ, THRESH_LFC, n_up, n_down))

# --- Tracked JSON summary --------------------------------------------------
desc <- paste0(
  "Full-transcriptome DESeq2 AF-vs-SR differential expression on donor-pooled ",
  "counts (paired LA+RA appendages summed per donor, 20 libraries -> 10 donors; ",
  "5 AF / 5 SR). Positive log2FC = higher in AF."
)
summary_list <- list(
  accession = ACCESSION,
  title = "Differentially Expressed Genes for Atrial Fibrillation Identified using RNA Sequencing from Paired Human Left and Right Atrial Appendages.",
  species = "Homo sapiens",
  platform = "GPL18573 (Illumina NextSeq 500, Homo sapiens)",
  tissue_model = "human paired left and right atrial appendage, bulk RNA-seq",
  contrast = CONTRAST,
  contrast_description = "positive log2FC = higher expression in atrial fibrillation than in sinus rhythm",
  method = paste0(
    "DESeq2 ", as.character(packageVersion("DESeq2")), " Wald test; ",
    "pairing pooled: counts summed across each donor's two appendages ",
    "(20 libraries -> 10 donors), design ~ group, BH adjustment; ",
    "DESeq2 independent filtering leaves padj NA for filtered genes"
  ),
  pairing = "pooled",
  pairing_rationale = paste0(
    "disease status is donor-level (each donor is AF or SR and contributes both ",
    "LA and RA), so a donor fixed block is confounded with the contrast; ",
    "DESeq2 has no random effects; pooling the paired libraries avoids ",
    "pseudo-replication"
  ),
  thresholds = list(padj = THRESH_PADJ, abs_log2FC = THRESH_LFC),
  n_genes_tested = as.integer(nrow(deg)),
  n_genes_padj = as.integer(sum(!is.na(deg$padj))),
  n_up = as.integer(n_up),
  n_down = as.integer(n_down),
  n_case = as.integer(sum(coldata$group == "AF")),
  n_control = as.integer(sum(coldata$group == "SR")),
  n_libraries = as.integer(ncol(counts)),
  exploratory = FALSE,
  description = desc,
  digests = list(),
  generated_by = "script/build_deg.R"
)
summary_list$digests <- list()
summary_list$digests[["source/GSE128188_CountsV2_EdgeR.txt.gz"]] <- sha256(COUNTS_FILE)
summary_list$digests[["local/deg/GSE128188_deg.csv"]] <- sha256(DEG_CSV)
writeLines(jsonlite::toJSON(summary_list, pretty = TRUE, auto_unbox = TRUE),
           DEG_JSON)

cat(sprintf("\nWrote %s (%d rows)\n", sub(paste0(ITEM_DIR, "/"), "", DEG_CSV),
            nrow(deg)))
cat(sprintf("Wrote %s\n", sub(paste0(ITEM_DIR, "/"), "", DEG_JSON)))
cat(sprintf("deg.csv sha256: %s\n", sha256(DEG_CSV)))
