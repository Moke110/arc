#!/usr/bin/env Rscript
#
# build_deg.R — full-transcriptome DESeq2 DEG for GSE282618
# (human hypertensive heart disease [HHD] vs control; post-mortem LV+RV bulk
# RNA-seq). ** EXPLORATORY ITEM ** — see the note below.
#
# Input (source/, ADR 0030 — never modified):
#   source/GSE282618_counts.txt.gz — 60,646 Ensembl gene IDs x 10 integer-count
#   libraries: 5 subjects (3 hypertensive: NDRI1, NDRI3, NDRI5; 2 control:
#   NDRI4, NDRI7), each contributing one LV and one RV library.
#
# EXPLORATORY, by ticket and by data: 3 HHD vs 2 control subjects (n = 5
# subjects) cannot support a stable transcriptome-wide DEG list. The gene list
# is a small-n reference for the gene-axis evidence, never a claim about the
# HHD transcriptome. `exploratory` is true in the tracked JSON and the
# author_log, and the tiny-n caveat is repeated in the JSON description.
#
# Design decision (same principle as GSE128188):
#   Disease status is a subject-level property and each subject contributes both
#   ventricles, so a subject fixed block is perfectly confounded with the HHD
#   contrast and DESeq2 has no random-effect term. We therefore pool the paired
#   LV+RV libraries per subject into 5 independent biological units and run the
#   unpaired 3 HHD vs 2 control DESeq2. `n_case`/`n_control` in the CSV are
#   subjects (3/2); the shipped expression matrix keeps all 10 GEO libraries
#   (script/build_expression.R).
#
# Gene symbols: the upstream matrix is Ensembl-ID-keyed (60,646 IDs). Symbols
# come from org.Hs.eg.db (ENSEMBL -> SYMBOL, multiVals = "first"); Ensembl IDs
# with no symbol in the pinned annotation database are dropped from both the
# DEG and the expression matrix (they cannot join a human-symbol axis; the drop
# is recorded in the JSON, not silent). Duplicate symbols are collapsed to the
# Ensembl ID with the highest mean CPM across the 10 libraries (ties broken by
# Ensembl ID) so the output is one row per gene symbol.
#
# Group labels are not in the supplementary table; they come from the GEO
# series matrix (!Sample_characteristics_ch1 "disease: hypertension|control").
# The subject lists are pinned below and validated against the counts columns.
#
# Outputs (both in local/, never in source/):
#   local/deg/GSE282618_deg.csv   full transcriptome (symbol-mapped), signed
#                                 HHD-vs-control effect
#   local/deg/GSE282618_deg.json  small TRACKED summary + input/output digests
#
# Usage:
#   cd disease_data/GSE282618
#   /home/chang/miniconda3/envs/biomni_e1/bin/Rscript script/build_deg.R
#
# Prerequisite R packages: DESeq2, jsonlite, digest, org.Hs.eg.db,
#   AnnotationDbi.

suppressPackageStartupMessages({
  library(DESeq2)
  library(jsonlite)
  library(digest)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

SCRIPT_DIR <- dirname(normalizePath(
  sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1])
))
ITEM_DIR <- dirname(SCRIPT_DIR)
SOURCE_DIR <- file.path(ITEM_DIR, "source")
LOCAL_DEG_DIR <- file.path(ITEM_DIR, "local", "deg")
dir.create(LOCAL_DEG_DIR, recursive = TRUE, showWarnings = FALSE)

COUNTS_FILE <- file.path(SOURCE_DIR, "GSE282618_counts.txt.gz")
DEG_CSV <- file.path(LOCAL_DEG_DIR, "GSE282618_deg.csv")
DEG_JSON <- file.path(LOCAL_DEG_DIR, "GSE282618_deg.json")

ACCESSION <- "GSE282618"
CONTRAST <- "HHD_vs_control"
N_CASE_EXPECTED <- 3L
N_CONTROL_EXPECTED <- 2L
THRESH_PADJ <- 0.05
THRESH_LFC <- 0.5

# GEO series matrix: NDRI1/NDRI3/NDRI5 are hypertensive, NDRI4/NDRI7 control.
HHD_SUBJECTS <- c("NDRI1", "NDRI3", "NDRI5")
CONTROL_SUBJECTS <- c("NDRI4", "NDRI7")

sha256 <- function(path) digest::digest(file = path, algo = "sha256")

cat("\n===== GSE282618 - human HHD vs control (LV+RV, EXPLORATORY n=5) =====\n\n")

if (!file.exists(COUNTS_FILE)) {
  stop("missing source/", basename(COUNTS_FILE),
       " — run 'bash script/fetch.sh' first")
}

# --- Read counts (Ensembl-ID-keyed) ---------------------------------------
raw <- read.delim(gzfile(COUNTS_FILE), check.names = FALSE,
                  stringsAsFactors = FALSE)
counts <- as.matrix(raw[, -1, drop = FALSE])
storage.mode(counts) <- "double"
ensg <- raw[[1]]
stopifnot(!anyDuplicated(ensg), all(grepl("^ENSG", ensg)))
rownames(counts) <- ensg
cat(sprintf("Counts: %d Ensembl genes x %d libraries\n",
            nrow(counts), ncol(counts)))

# --- Sample metadata -------------------------------------------------------
sample <- colnames(counts)
subject <- sub(" .*$", "", sample)
site <- sub("^.* ", "", sample)
if (!setequal(site, c("LV", "RV"))) {
  stop("unexpected ventricular-site labels in counts columns: ",
       paste(sort(unique(site)), collapse = ", "))
}
if (!setequal(unique(subject), c(HHD_SUBJECTS, CONTROL_SUBJECTS))) {
  stop("subject IDs in the counts table do not match the GEO series-matrix subjects")
}
group <- ifelse(subject %in% HHD_SUBJECTS, "HHD", "control")
meta <- data.frame(sample = sample, subject = subject, site = site,
                   group = group, stringsAsFactors = FALSE)
if (any(table(meta$subject) != 2L)) {
  stop("every subject must contribute exactly one LV and one RV library")
}
cat(sprintf("Libraries: HHD %d, control %d; subjects: %d (HHD %d, control %d)\n",
            sum(meta$group == "HHD"), sum(meta$group == "control"),
            length(unique(meta$subject)),
            length(unique(meta$subject[meta$group == "HHD"])),
            length(unique(meta$subject[meta$group == "control"]))))

# --- Ensembl ID -> human symbol, collapse duplicate symbols ---------------
# KEEP IN SYNC with script/build_expression.R (same annotation and collapse
# rule so the DEG table and the expression matrix share one gene set).
cat("Mapping Ensembl IDs to symbols via org.Hs.eg.db...\n")
symbol <- suppressWarnings(mapIds(
  org.Hs.eg.db, keys = ensg, column = "SYMBOL", keytype = "ENSEMBL",
  multiVals = "first"
))
mapped <- !is.na(symbol) & nzchar(symbol)
n_unmapped <- sum(!mapped)
cat(sprintf("  mapped: %d / %d (dropped %d without a symbol)\n",
            sum(mapped), length(ensg), n_unmapped))

# Mean CPM across libraries, used only to choose one Ensembl ID per symbol.
lib_size <- colSums(counts)
cpm <- sweep(counts, 2, lib_size / 1e6, "/")
mean_cpm <- rowMeans(cpm)
sel <- data.frame(
  ensg = ensg[mapped],
  symbol = as.character(symbol[mapped]),
  mean_cpm = mean_cpm[mapped],
  stringsAsFactors = FALSE
)
sel <- sel[order(sel$symbol, -sel$mean_cpm, sel$ensg, method = "radix"), ]
dup_symbol <- duplicated(sel$symbol)
n_collapsed <- sum(dup_symbol)
sel <- sel[!dup_symbol, ]
cat(sprintf("  symbol genes: %d (collapsed %d duplicate symbols)\n",
            nrow(sel), n_collapsed))

counts <- counts[sel$ensg, , drop = FALSE]
rownames(counts) <- sel$symbol

# --- Pool the paired LV+RV libraries per subject ---------------------------
subjects <- sort(unique(meta$subject))
pooled <- vapply(subjects, function(s) {
  rowSums(counts[, meta$subject == s, drop = FALSE])
}, numeric(nrow(counts)))
colnames(pooled) <- subjects
subject_group <- ifelse(subjects %in% HHD_SUBJECTS, "HHD", "control")
coldata <- data.frame(
  group = factor(subject_group, levels = c("control", "HHD")),
  row.names = subjects
)
stopifnot(sum(coldata$group == "control") == N_CONTROL_EXPECTED,
          sum(coldata$group == "HHD") == N_CASE_EXPECTED)
cat(sprintf("Subject-pooled matrix: %d genes x %d subjects (HHD %d, control %d)\n",
            nrow(pooled), ncol(pooled),
            sum(coldata$group == "HHD"), sum(coldata$group == "control")))

# --- DESeq2: HHD vs control ------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = round(pooled), colData = coldata,
                              design = ~ group)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("group", "HHD", "control"))
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
  n_case = as.integer(sum(coldata$group == "HHD")),
  n_control = as.integer(sum(coldata$group == "control")),
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
  "EXPLORATORY — full-transcriptome DESeq2 HHD-vs-control differential ",
  "expression on subject-pooled counts (paired LV+RV libraries summed per ",
  "subject, 10 libraries -> 5 subjects; 3 HHD / 2 control). n=5 subjects ",
  "cannot support a stable DEG list; this is a small-n reference only. ",
  "Positive log2FC = higher in HHD. Hypertensive heart disease is a ",
  "borderline >5% prevalence category (docs/research/",
  "stage2-stage3-gene-axis-datasets.md, 'Disease-prevalence taxonomy')."
)
summary_list <- list(
  accession = ACCESSION,
  title = "Multiomics Investigation of the Hypertensive Human Heart",
  species = "Homo sapiens",
  platform = "GPL30173 (NextSeq 2000, Homo sapiens)",
  tissue_model = "human post-mortem left and right ventricle, bulk RNA-seq",
  contrast = CONTRAST,
  contrast_description = "positive log2FC = higher expression in hypertensive heart disease than in control",
  method = paste0(
    "DESeq2 ", as.character(packageVersion("DESeq2")), " Wald test; ",
    "pairing pooled: counts summed across each subject's LV+RV libraries ",
    "(10 libraries -> 5 subjects), design ~ group, BH adjustment; ",
    "DESeq2 independent filtering leaves padj NA for filtered genes"
  ),
  pairing = "pooled",
  pairing_rationale = paste0(
    "disease status is subject-level (each subject is HHD or control and ",
    "contributes both LV and RV), so a subject fixed block is confounded with ",
    "the contrast; DESeq2 has no random effects; pooling the paired libraries ",
    "avoids pseudo-replication"
  ),
  symbol_mapping = paste0(
    "org.Hs.eg.db ", as.character(packageVersion("org.Hs.eg.db")),
    " ENSEMBL -> SYMBOL (multiVals = first); unmapped Ensembl IDs dropped; ",
    "duplicate symbols collapsed to the Ensembl ID with the highest mean CPM"
  ),
  thresholds = list(padj = THRESH_PADJ, abs_log2FC = THRESH_LFC),
  n_genes_tested = as.integer(nrow(deg)),
  n_genes_padj = as.integer(sum(!is.na(deg$padj))),
  n_up = as.integer(n_up),
  n_down = as.integer(n_down),
  n_case = as.integer(sum(coldata$group == "HHD")),
  n_control = as.integer(sum(coldata$group == "control")),
  n_subjects = as.integer(length(unique(meta$subject))),
  n_libraries = as.integer(ncol(counts)),
  source_genes = as.integer(length(ensg)),
  dropped_no_symbol = as.integer(n_unmapped),
  duplicate_symbols_collapsed = as.integer(n_collapsed),
  exploratory = TRUE,
  exploratory_reason = paste0(
    "3 hypertensive vs 2 control subjects (n=5); the dataset survey records ",
    "the series as exploratory because n=5 subjects cannot support a stable ",
    "transcriptome-wide DEG list"
  ),
  description = desc,
  digests = list(),
  generated_by = "script/build_deg.R"
)
summary_list$digests <- list()
summary_list$digests[["source/GSE282618_counts.txt.gz"]] <- sha256(COUNTS_FILE)
summary_list$digests[["local/deg/GSE282618_deg.csv"]] <- sha256(DEG_CSV)
writeLines(jsonlite::toJSON(summary_list, pretty = TRUE, auto_unbox = TRUE),
           DEG_JSON)

cat(sprintf("\nWrote %s (%d rows)\n", sub(paste0(ITEM_DIR, "/"), "", DEG_CSV),
            nrow(deg)))
cat(sprintf("Wrote %s\n", sub(paste0(ITEM_DIR, "/"), "", DEG_JSON)))
cat(sprintf("deg.csv sha256: %s\n", sha256(DEG_CSV)))
