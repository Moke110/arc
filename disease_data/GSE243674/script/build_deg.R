#!/usr/bin/env Rscript
#
# build_deg.R — full-transcriptome DESeq2 for GSE243674 (human hiPSC-CM
# chemotherapy cardiotoxicity, 6 healthy female donors, 0.5 uM, 3 h / 24 h).
#
# Input (source/, untouched upstream download, ADR 0030):
#   source/GSE243674_Counts_RNA_ERMatthews.txt.gz
#   ENTREZID x 72 wells, space-delimited, counts from featureCounts on hg38.
#   Column labels are `<DRUG>.<donor>.<time>` with DRUG in
#   {DNR, DOX, EPI, MTX, TRZ, VEH}. The drug names actually present are the
#   abbreviations in that header; GEO series-matrix metadata expands them to
#   daunorubicin, doxorubicin, epirubicin, mitoxantrone, trastuzumab and vehicle
#   (verified against GSE243674's `!Sample_characteristics_ch1` treatment rows;
#   see author_log.json). Nothing is assumed: an unknown header token aborts
#   this script instead of being silently relabelled.
#
# Analysis:
#   One DESeq2 fit per contrast (drug x time vs vehicle at the same time point),
#   on the matched 12 wells of that time (6 drug + 6 vehicle, same 6 donors).
#   Design `~ donor + group` — each donor contributes a matched vehicle and drug
#   well, so the donor term removes hiPSC-line baseline differences and the
#   coefficient is the donor-adjusted drug-vs-vehicle effect. No LFC shrinkage:
#   `log2FC` is the DESeq2 MLE and `stat` the Wald statistic. BH padj with
#   DESeq2's default independent filtering (alpha = 0.05). No L1000 filtering.
#
# Mapping (ENTREZID -> HGNC symbol):
#   The pinned `pert_data/Annotation/source/hgnc_complete_set.txt` table is used
#   (no download, no org.Hs.eg.db round trip). Approved rows only; an Entrez ID
#   with several approved HGNC rows is resolved deterministically
#   (protein-coding < non-coding RNA < pseudogene < other, then lowest hgnc_id)
#   and every such case, plus every lost Entrez ID and its read share, is
#   recorded in the JSON. Set ARC_HGNC_TABLE to override the path in checkouts
#   (e.g. a worktree) where the gitignored pert_data source is not materialised.
#
# Outputs (declared in author_log.json):
#   local/deg/GSE243674_deg.csv   long table, one row per contrast x gene:
#                                 contrast, gene_symbol, log2FC (drug vs vehicle),
#                                 stat, p_value, padj, baseMean, n_case, n_control
#   local/deg/GSE243674_deg.json  tracked compact summary: per-contrast counts,
#                                 species/platform/cell model/method/thresholds,
#                                 mapping evidence (file + sha256 + losses),
#                                 drug-class labels and input/output digests
#
# Usage:
#   cd disease_data/GSE243674
#   Rscript script/build_deg.R
#
# Pinned R: /home/chang/miniconda3/envs/biomni_e1/bin/Rscript (DESeq2 1.x).
# Cap the run: ulimit -v 6G; OMP_NUM_THREADS=2. Peak RSS is recorded in
# author_log.json / the ticket record.

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(digest)
  library(jsonlite)
})

options(stringsAsFactors = FALSE, datatable.verbose = FALSE)

# --- Locate the item from the script path (paths do not encode depth) --------
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("run this script as: Rscript script/build_deg.R")
SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", script_arg)))
find_repo_root <- function(start_dir) {
  dir <- normalizePath(start_dir, mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(dir, ".git"))) return(dir)
    parent <- dirname(dir)
    if (identical(parent, dir)) {
      stop("cannot find the repository root (no .git directory above ", start_dir, ")")
    }
    dir <- parent
  }
}
REPO_ROOT <- find_repo_root(SCRIPT_DIR)
ITEM_DIR <- dirname(SCRIPT_DIR)                       # disease_data/GSE243674
SOURCE_DIR <- file.path(ITEM_DIR, "source")
LOCAL_DEG_DIR <- file.path(ITEM_DIR, "local", "deg")
dir.create(LOCAL_DEG_DIR, recursive = TRUE, showWarnings = FALSE)

SOURCE_FILE <- file.path(SOURCE_DIR, "GSE243674_Counts_RNA_ERMatthews.txt.gz")
HGNC_TABLE <- Sys.getenv(
  "ARC_HGNC_TABLE",
  unset = file.path(REPO_ROOT, "pert_data", "Annotation", "source", "hgnc_complete_set.txt")
)
OUT_CSV <- file.path(LOCAL_DEG_DIR, "GSE243674_deg.csv")
OUT_JSON <- file.path(LOCAL_DEG_DIR, "GSE243674_deg.json")

if (!file.exists(SOURCE_FILE)) {
  stop("missing ", SOURCE_FILE, " — run `bash script/fetch.sh` first")
}
if (!file.exists(HGNC_TABLE)) {
  stop("pinned HGNC table not found: ", HGNC_TABLE,
       "\nIt lives in the gitignored pert_data/Annotation/source/ tree; in a",
       " checkout without it, point ARC_HGNC_TABLE at a materialised copy.")
}

sha256_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
source_sha256 <- sha256_file(SOURCE_FILE)
source_bytes <- file.info(SOURCE_FILE)$size
hgnc_sha256 <- sha256_file(HGNC_TABLE)

# --- The drug arms actually present in the source header ---------------------
# Full names come from the GEO series-matrix metadata, not from guesswork; the
# abbreviation is what this item's own source file carries.
DRUG_META <- data.table(
  token = c("DNR", "DOX", "EPI", "MTX", "TRZ", "VEH"),
  name = c("daunorubicin", "doxorubicin", "epirubicin", "mitoxantrone",
           "trastuzumab", "vehicle"),
  drug_class = c("anthracycline TOP2 inhibitor", "anthracycline TOP2 inhibitor",
                 "anthracycline TOP2 inhibitor", "anthraquinone TOP2 inhibitor",
                 "HER2-targeted antibody", "vehicle")
)

# --- Read the count matrix ---------------------------------------------------
cat("GSE243674: reading", basename(SOURCE_FILE), "\n")
counts_dt <- fread(cmd = paste("gzip -cd", shQuote(SOURCE_FILE)), sep = " ", header = TRUE)
if (ncol(counts_dt) != 73L || names(counts_dt)[1] != "ENTREZID") {
  stop("unexpected matrix shape: ", nrow(counts_dt), " x ", ncol(counts_dt),
       " (expected ENTREZID + 72 wells)")
}
wells <- names(counts_dt)[-1]
counts <- as.matrix(counts_dt[, -1, with = FALSE])
rownames(counts) <- as.character(counts_dt[[1]])
gene_ids <- rownames(counts)
if (anyDuplicated(gene_ids)) stop("duplicated ENTREZID rows in the source matrix")

# --- Decode wells: <DRUG>.<donor>.<time>h -----------------------------------
parsed <- regmatches(wells, regexec("^([A-Z]+)\\.([0-9]+)\\.([0-9]+)h$", wells))
if (any(lengths(parsed) != 4L)) {
  stop("cannot parse well labels: ", paste(wells[lengths(parsed) != 4L], collapse = ", "))
}
well_dt <- data.table(
  sample = wells,
  token = vapply(parsed, `[`, character(1), 2L),
  donor = vapply(parsed, `[`, character(1), 3L),
  time_h = as.integer(vapply(parsed, `[`, character(1), 4L))
)
unknown <- setdiff(unique(well_dt$token), DRUG_META$token)
if (length(unknown) > 0L) {
  stop("unknown drug token(s) in the source header: ", paste(unknown, collapse = ", "),
       "\nVerify the arms against GEO GSE243674 metadata and extend DRUG_META.")
}
well_dt <- merge(well_dt, DRUG_META, by = "token", sort = FALSE)
setorder(well_dt, donor, time_h, token)
well_dt[, group := ifelse(token == "VEH", "vehicle", "treated")]
if (nrow(well_dt) != 72L || length(unique(well_dt$donor)) != 6L ||
    !setequal(unique(well_dt$time_h), c(3L, 24L)) ||
    sum(well_dt$group == "vehicle") != 12L) {
  stop("unbalanced design in the source header: ", nrow(well_dt), " wells")
}
cat(sprintf("  %d wells x %d ENTREZIDs; %d donors; times %s; drugs %s\n",
            length(wells), nrow(counts), length(unique(well_dt$donor)),
            paste(sort(unique(well_dt$time_h)), collapse = "/"),
            paste(sort(unique(well_dt$token[well_dt$token != "VEH"])), collapse = "/")))

# --- ENTREZID -> HGNC symbol with the pinned table --------------------------
hg <- fread(HGNC_TABLE, sep = "\t", header = TRUE, quote = "",
            select = c("hgnc_id", "symbol", "status", "locus_group", "entrez_id"))
if (!all(c("hgnc_id", "symbol", "status", "entrez_id") %in% names(hg))) {
  stop("pinned HGNC table lacks the expected columns: ", paste(names(hg), collapse = ", "))
}
hg <- hg[status == "Approved" & !is.na(entrez_id) & entrez_id != ""]
hg[, entrez_id := as.character(entrez_id)]
hg[, locus_rank := match(locus_group,
                         c("protein-coding gene", "non-coding RNA", "pseudogene"),
                         nomatch = 99L)]
hg[, hgnc_num := as.integer(sub("^HGNC:", "", hgnc_id))]
setorder(hg, entrez_id, locus_rank, hgnc_num)
multi <- hg[, .N, by = entrez_id][N > 1L & entrez_id %in% gene_ids]
multi_detail <- hg[entrez_id %in% multi$entrez_id,
                   .(entrez_id, symbol, locus_group, hgnc_id)]
# deterministic first row per Entrez ID after the ordering above
map <- hg[!duplicated(entrez_id), .(entrez_id, symbol)]

hit <- map$symbol[match(gene_ids, map$entrez_id)]
mapped <- !is.na(hit)
n_input <- length(gene_ids)
n_mapped_ids <- sum(mapped)
n_unmapped_ids <- n_input - n_mapped_ids
total_reads <- sum(counts)
unmapped_read_share <- if (n_unmapped_ids > 0L) {
  sum(counts[!mapped, , drop = FALSE]) / total_reads
} else 0

counts <- counts[mapped, , drop = FALSE]
symbols <- hit[mapped]
n_genes_all_zero <- sum(rowSums(counts) == 0L)
# Collapse several Entrez IDs that map to one symbol (sum counts, the usual
# gene-level rule). Recorded even when zero so the mapping evidence is complete.
dup_symbols <- unique(symbols[duplicated(symbols)])
if (length(dup_symbols) > 0L) {
  counts <- rowsum(counts, group = symbols)
  symbols <- rownames(counts)
}
rownames(counts) <- symbols  # gene_symbol is the gene-level index from here on
if (anyDuplicated(symbols)) stop("gene-symbol collapse left duplicates")
n_genes <- length(symbols)
cat(sprintf("  mapped %d/%d ENTREZIDs to HGNC symbols (%d lost, %.3f%% of reads); %d genes; %d all-zero\n",
            n_mapped_ids, n_input, n_unmapped_ids,
            100 * unmapped_read_share, n_genes, n_genes_all_zero))

# --- Per-contrast DESeq2 ----------------------------------------------------
contrasts <- unique(well_dt[token != "VEH",
                            .(drug_token = token, drug = name, drug_class, time_h)])
contrasts[, contrast := sprintf("%s_%dh_vs_vehicle", tolower(drug_token), time_h)]
setorder(contrasts, drug_token, time_h)

rows <- vector("list", nrow(contrasts))
summary_rows <- vector("list", nrow(contrasts))

for (i in seq_len(nrow(contrasts))) {
  cst <- contrasts[i]
  sel <- well_dt[(token == cst$drug_token & time_h == cst$time_h) |
                   (token == "VEH" & time_h == cst$time_h)]
  # `group` must have vehicle as the reference level so log2FC = drug - vehicle.
  sel[, group := factor(group, levels = c("vehicle", "treated"))]
  sel[, donor := factor(donor)]
  mat <- counts[, sel$sample, drop = FALSE]
  meta <- data.frame(row.names = sel$sample, group = sel$group, donor = sel$donor)

  dds <- DESeqDataSetFromMatrix(countData = mat, colData = meta, design = ~ donor + group)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("group", "treated", "vehicle"),
                 independentFiltering = TRUE, alpha = 0.05)
  res_dt <- data.table(
    contrast = cst$contrast,
    gene_symbol = rownames(res),
    log2FC = res$log2FoldChange,
    stat = res$stat,
    p_value = res$pvalue,
    padj = res$padj,
    baseMean = res$baseMean,
    n_case = sum(meta$group == "treated"),
    n_control = sum(meta$group == "vehicle")
  )
  # deterministic display order: most significant first, NA padj last
  res_dt[, padj_na := is.na(padj)]
  setorder(res_dt, padj_na, padj, gene_symbol)
  res_dt[, padj_na := NULL]
  rows[[i]] <- res_dt

  sig <- !is.na(res_dt$padj) & res_dt$padj < 0.05
  summary_rows[[i]] <- list(
    contrast = cst$contrast,
    drug = cst$drug,
    drug_token = cst$drug_token,
    drug_class = cst$drug_class,
    time_h = as.integer(cst$time_h),
    n_case = sum(meta$group == "treated"),
    n_control = sum(meta$group == "vehicle"),
    n_genes_tested = nrow(res_dt),
    n_padj_not_na = sum(!is.na(res_dt$padj)),
    n_up = sum(sig & res_dt$log2FC > 0.5),
    n_down = sum(sig & res_dt$log2FC < -0.5),
    n_significant = sum(sig)
  )
  cat(sprintf("  %-22s up %5d  down %5d  tested %6d  (padj<0.05 & |log2FC|>0.5)\n",
              cst$contrast, summary_rows[[i]]$n_up, summary_rows[[i]]$n_down,
              summary_rows[[i]]$n_genes_tested))
}

deg <- rbindlist(rows)
stopifnot(identical(names(deg),
                    c("contrast", "gene_symbol", "log2FC", "stat", "p_value",
                      "padj", "baseMean", "n_case", "n_control")))
fwrite(deg, OUT_CSV, sep = ",", quote = FALSE, na = "NA")
cat("GSE243674: wrote", OUT_CSV, sprintf("(%.1f MB, %d rows)\n",
                                         file.info(OUT_CSV)$size / 1e6, nrow(deg)))

# --- Tracked compact summary ------------------------------------------------
# Counts are taken from the table that was just written, so JSON and CSV agree
# by construction. No timestamps: rebuilding on the pinned R reproduces the JSON
# byte-for-byte (the CSV digest is recorded and is stable on that R too).
csv_sha256 <- sha256_file(OUT_CSV)
deg_json <- list(
  item = "GSE243674",
  accession = "GSE243674",
  accession_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243674",
  title = "Anthracyclines induce cardiotoxicity through a shared gene expression response signature",
  species = "Homo sapiens",
  species_taxid = 9606L,
  platform = "GPL21697 NextSeq 550 (Homo sapiens); featureCounts gene-level counts on hg38",
  cell_model = "hiPSC-CM monolayer, 6 healthy female donors (UCSD iPSC lines), 0.5 uM, 3 h / 24 h",
  contrast_direction = "drug vs vehicle (positive log2FC = higher in drug-treated wells)",
  design = list(
    n_wells = nrow(well_dt),
    n_donors = length(unique(well_dt$donor)),
    timepoints_h = sort(unique(well_dt$time_h)),
    vehicle_wells = sum(well_dt$group == "vehicle"),
    arms = lapply(seq_len(nrow(DRUG_META)), function(k) list(
      token = DRUG_META$token[k],
      drug = DRUG_META$name[k],
      drug_class = DRUG_META$drug_class[k],
      wells = sum(well_dt$token == DRUG_META$token[k])
    ))
  ),
  method = list(
    deg = "DESeq2 (DESeq2::DESeq + results) per drug x time vs vehicle at the same time point",
    model = "~ donor + group; donor is a fixed effect because each donor contributes matched vehicle and drug wells",
    test = "Wald test, unshrunken MLE log2FoldChange, BH padj, DESeq2 default independent filtering (alpha = 0.05)",
    genes = "all mapped genes (no L1000 or expression-level filter); all-zero genes carry NA statistics",
    normalization = "DESeq2 median-of-ratios size factors within each 12-well contrast",
    drug_names_source = "column tokens DNR/DOX/EPI/MTX/TRZ/VEH are expanded using GSE243674 series-matrix metadata (!Sample_characteristics_ch1 treatment)"
  ),
  thresholds = list(padj = 0.05, abs_log2FC = 0.5),
  contrasts = summary_rows,
  mapping = list(
    gene_id_type = "ENTREZID",
    table = "pert_data/Annotation/source/hgnc_complete_set.txt",
    table_sha256 = hgnc_sha256,
    table_bytes = file.info(HGNC_TABLE)$size,
    table_rows_approved = nrow(hg),
    input_entrez = n_input,
    mapped_entrez = n_mapped_ids,
    unmapped_entrez = n_unmapped_ids,
    unmapped_read_share = round(unmapped_read_share, 6),
    one_to_many_entrez = nrow(multi),
    one_to_many_resolution = lapply(seq_len(nrow(multi)), function(k) list(
      entrez_id = multi$entrez_id[k],
      candidates = multi_detail[entrez_id == multi$entrez_id[k], symbol],
      chosen = map[entrez_id == multi$entrez_id[k], symbol],
      rule = "protein-coding < non-coding RNA < pseudogene < other, then lowest hgnc_id"
    )),
    duplicate_symbols_collapsed = length(dup_symbols),
    genes = n_genes,
    genes_all_zero = n_genes_all_zero
  ),
  inputs = list(
    path = "source/GSE243674_Counts_RNA_ERMatthews.txt.gz",
    bytes = source_bytes,
    sha256 = source_sha256,
    manifest = "manifest.sha256"
  ),
  outputs = list(
    deg_csv = "local/deg/GSE243674_deg.csv",
    deg_csv_sha256 = csv_sha256,
    deg_csv_rows = nrow(deg)
  ),
  generated_by = "script/build_deg.R"
)
writeLines(toJSON(deg_json, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null"),
           OUT_JSON)
cat("GSE243674: wrote", OUT_JSON, "\n")
cat("GSE243674 build_deg done\n")
