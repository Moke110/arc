#!/usr/bin/env Rscript
#
# Per-item sample annotation + expression loader for GSE141910 (human LV bulk
# RNA-seq, MAGNet consortium) — shared by script/build_deg.R and
# script/build_expression.R so the sample labels and the normalization are
# defined in exactly one place.
#
# SAMPLE LABELS. GSE141910's 366 per-sample supplementary CSVs are named
# `GSM<id>_<title>.csv.gz`; the case/control label is the `Sample_title`
# (= the per-sample file basename) mapped to `Sample_characteristics_ch1:
# etiology`, extracted from the GEO series matrix
# (https://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141910/matrix/GSE141910_series_matrix.txt.gz)
# on 2026-11-13. The series matrix is NOT part of source/ (ticket #179 scope:
# exactly GSE141910_RAW.tar + filelist.txt), so the extracted map is pinned here
# as the tracked input; it reproduces from this file every run. HF = DCM + HCM +
# PPCM (200 samples), NF = non-failing donor hearts (166 samples).
#
# EXPRESSION VALUES. Each per-sample CSV holds one column of log2-scale
# normalized abundance for 20,781 Ensembl gene IDs. GEO's own
# `Sample_data_processing` says exactly what they are: "Raw gene counts were
# transformed using SVA and the 'LIMMA-Voom' method" and
# "Supplementary_files_format_and_content: Voom transformed abundance
# measurements". So the deposited matrix is already the log2-scale,
# library-size-normalized, submitter-SVA-processed voom transform; this item
# does NOT re-normalize it. The HF-vs-NF limma fit in build_deg.R therefore runs
# on these deposited values (limma on log2-scale expression) — the raw counts
# were never deposited, so DESeq2/voom weighting cannot be re-run from source/.

GSE141910_SAMPLES <- data.frame(
  sample = c("C00039", "C00055", "C00074", "C00085", "C00105", "C00120", "C00132", "C00158", "C00168", "C00186", "C00188", "C00190", "C00192", "C00193", "C00195", "C00196", "C00197", "C00198", "C00199", "C00202", "C01413", "C01418", "C01631", "C01776", "C01783", "C01794", "C01809", "C01814", "C01816", "C01862", "C01875", "C01881", "C01895", "C01902", "C01904", "C01905", "C01910", "C01912", "C01914", "C01926", "C01929", "C01931", "C01945", "C01952", "C01957", "C01960", "C01997", "C02014", "C02018", "C02022", "C02029", "C02045", "C02050", "C02054", "C02056", "C02058", "C02065", "C02067", "C02069", "C02080", "C02081", "C02085", "C02095", "C02096", "C02101", "C02108", "C02109", "C02116", "C02164", "C02172", "C02190", "C02205", "C02210", "C02233", "C02242", "C02272", "C02276", "C02280", "C02297", "C02298", "C02301", "C02303", "C02305", "C02312", "C02318", "C02327", "C02328", "C02329", "C02340", "C02344", "C02349", "C02350", "C02351", "C02352", "C02354", "C02355", "C02356", "C02360", "C02363", "C02364", "C02374", "C02375", "C02380", "C02387", "C02397", "C02420", "C02452", "C02458", "C02464", "C02468", "C02475", "C02491", "C02496", "C02510", "C02516", "C02529", "C02530", "C02546", "C02592", "C02645", "C02660", "C02665", "C02698", "C02720", "C02724", "C02728", "P01007", "P01027", "P01032", "P01044", "P01045", "P01052", "P01054", "P01067", "P01070", "P01071", "P01076", "P01077", "P01095", "P01097", "P01100", "P01101", "P01104", "P01106", "P01109", "P01110", "P01116", "P01117", "P01118", "P01119", "P01121", "P01122", "P01124", "P01127", "P01134", "P01135", "P01137", "P01140", "P01146", "P01151", "P01155", "P01163", "P01164", "P01165", "P01167", "P01168", "P01171", "P01173", "P01179", "P01183", "P01184", "P01185", "P01197", "P01199", "P01202", "P01204", "P01205", "P01206", "P01207", "P01208", "P01209", "P01211", "P01215", "P01217", "P01221", "P01222", "P01224", "P01226", "P01228", "P01233", "P01238", "P01241", "P01245", "P01246", "P01247", "P01248", "P01252", "P01254", "P01259", "P01262", "P01264", "P01266", "P01275", "P01279", "P01282", "P01287", "P01288", "P01289", "P01290", "P01292", "P01293", "P01294", "P01296", "P01298", "P01300", "P01304", "P01305", "P01307", "P01309", "P01310", "P01318", "P01320", "P01322", "P01323", "P01325", "P01329", "P01330", "P01334", "P01337", "P01339", "P01340", "P01347", "P01348", "P01350", "P01352", "P01354", "P01356", "P01357", "P01358", "P01359", "P01360", "P01362", "P01365", "P01367", "P01369", "P01371", "P01374", "P01375", "P01377", "P01379", "P01380", "P01381", "P01383", "P01385", "P01387", "P01388", "P01396", "P01397", "P01400", "P01402", "P01404", "P01406", "P01412", "P01413", "P01414", "P01420", "P01421", "P01422", "P01425", "P01430", "P01433", "P01436", "P01437", "P01440", "P01444", "P01447", "P01452", "P01454", "P01457", "P01462", "P01463", "P01465", "P01466", "P01469", "P01471", "P01472", "P01473", "P01479", "P01481", "P01485", "P01488", "P01490", "P01492", "P01497", "P01498", "P01500", "P01503", "P01504", "P01507", "P01508", "P01509", "P01510", "P01513", "P01515", "P01516", "P01521", "P01528", "P01532", "P01535", "P01537", "P01538", "P01539", "P01540", "P01543", "P01545", "P01547", "P01549", "P01556", "P01558", "P01559", "P01561", "P01562", "P01563", "P01567", "P01568", "P01573", "P01578", "P01580", "P01581", "P01582", "P01584", "P01585", "P01588", "P01590", "P01591", "P01592", "P01593", "P01595", "P01596", "P01597", "P01600", "P01602", "P01603", "P01605", "P01606", "P01608", "P01609", "P01610", "P01611", "P01615", "P01616", "P01617", "P01618", "P01621", "P01622", "P01623", "P01626", "P01627", "P01628", "P01629", "P01630", "P01631", "P01634", "P01635", "P01639", "P01640"),
  group = c("NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "HF", "NF", "HF", "HF", "NF", "HF", "NF", "NF", "HF", "NF", "HF", "HF", "HF", "NF", "NF", "HF", "NF", "HF", "HF", "NF", "NF", "HF", "NF", "NF", "NF", "NF", "HF", "NF", "HF", "NF", "NF", "NF", "NF", "NF", "HF", "HF", "NF", "NF", "HF", "NF", "NF", "HF", "NF", "NF", "NF", "NF", "NF", "NF", "HF", "NF", "NF", "NF", "NF", "HF", "NF", "NF", "HF", "NF", "NF", "NF", "NF", "HF", "NF", "NF", "NF", "NF", "HF", "HF", "NF", "NF", "HF", "NF", "HF", "HF", "HF", "HF", "NF", "HF", "NF", "HF", "HF", "HF", "HF", "HF", "NF", "NF", "HF", "HF", "HF", "HF", "NF", "HF", "NF", "HF", "HF", "HF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "HF", "NF", "NF", "HF", "HF", "NF", "NF", "NF", "HF", "NF", "NF", "NF", "NF", "NF", "NF", "HF", "NF", "NF", "NF", "HF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "HF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "HF", "HF", "HF", "HF", "NF", "HF", "NF", "NF", "HF", "NF", "NF", "NF", "HF", "NF", "NF", "HF", "NF", "NF", "HF", "NF", "HF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "HF", "NF", "HF", "HF", "HF", "NF", "NF", "NF", "NF", "NF", "HF", "HF", "HF", "HF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "HF", "HF", "HF", "HF", "NF", "NF", "NF", "NF", "NF", "HF", "NF", "HF", "NF", "HF", "NF", "NF", "HF", "HF", "NF", "HF", "NF", "HF", "HF", "HF", "HF", "NF", "HF", "HF", "NF", "HF", "HF", "NF", "NF", "NF", "NF", "HF", "NF", "HF", "HF", "HF", "HF", "NF", "HF", "HF"),
  etiology = c("NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "PPCM", "PPCM", "HCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "HCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "PPCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "HCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "HCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "HCM", "HCM", "DCM", "DCM", "HCM", "DCM", "DCM", "DCM", "HCM", "DCM", "HCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "HCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "DCM", "HCM", "NF", "DCM", "DCM", "NF", "DCM", "NF", "NF", "DCM", "NF", "DCM", "DCM", "DCM", "NF", "NF", "DCM", "NF", "DCM", "DCM", "NF", "NF", "DCM", "NF", "NF", "NF", "NF", "DCM", "NF", "DCM", "NF", "NF", "NF", "NF", "NF", "DCM", "DCM", "NF", "NF", "PPCM", "NF", "NF", "DCM", "NF", "NF", "NF", "NF", "NF", "NF", "DCM", "NF", "NF", "NF", "NF", "HCM", "NF", "NF", "DCM", "NF", "NF", "NF", "NF", "HCM", "NF", "NF", "NF", "NF", "PPCM", "DCM", "NF", "NF", "DCM", "NF", "DCM", "DCM", "DCM", "DCM", "NF", "DCM", "NF", "HCM", "DCM", "DCM", "DCM", "PPCM", "NF", "NF", "DCM", "DCM", "DCM", "DCM", "NF", "DCM", "NF", "HCM", "HCM", "DCM", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "HCM", "NF", "NF", "DCM", "DCM", "NF", "NF", "NF", "DCM", "NF", "NF", "NF", "NF", "NF", "NF", "DCM", "NF", "NF", "NF", "HCM", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "DCM", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "HCM", "HCM", "DCM", "DCM", "NF", "DCM", "NF", "NF", "HCM", "NF", "NF", "NF", "HCM", "NF", "NF", "DCM", "NF", "NF", "DCM", "NF", "HCM", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "DCM", "NF", "HCM", "DCM", "HCM", "NF", "NF", "NF", "NF", "NF", "DCM", "DCM", "DCM", "DCM", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF", "DCM", "DCM", "DCM", "DCM", "NF", "NF", "NF", "NF", "NF", "DCM", "NF", "DCM", "NF", "DCM", "NF", "NF", "DCM", "DCM", "NF", "HCM", "NF", "DCM", "DCM", "DCM", "DCM", "NF", "DCM", "DCM", "NF", "DCM", "DCM", "NF", "NF", "NF", "NF", "DCM", "NF", "DCM", "HCM", "HCM", "DCM", "NF", "DCM", "DCM"),
  stringsAsFactors = FALSE
)

# ENSG -> HGNC symbol. org.Hs.eg.db (pinned for this repository's external
# items) is a current-release database; GSE141910 was quantified against
# Ensembl v67 (hg19), so retired/merged v67 IDs do not map. Those genes are kept
# and labelled with the Ensembl ID as the `gene_symbol`, so the DEG table stays
# full-transcriptome instead of silently dropping ~18% of the measured genes.
map_ensg_to_symbol <- function(ensg) {
  suppressPackageStartupMessages({
    library(org.Hs.eg.db)
    library(AnnotationDbi)
  })
  sym <- mapIds(org.Hs.eg.db, keys = ensg, column = "SYMBOL",
                keytype = "ENSEMBL", multiVals = "first")
  sym <- as.character(sym)
  sym[is.na(sym) | sym == ""] <- ensg[is.na(sym) | sym == ""]
  sym
}

# Read the 366 per-sample CSVs out of source/GSE141910_RAW.tar and return
#   list(expr = genes x samples matrix (log2 abundance, as deposited),
#        sample = sample titles in column order,
#        group = "HF"/"NF", etiology = DCM/HCM/PPCM/NF,
#        ensg = rownames(expr))
load_gse141910_expression <- function(source_dir) {
  tarfile <- file.path(source_dir, "GSE141910_RAW.tar")
  if (!file.exists(tarfile)) {
    stop("GSE141910_RAW.tar not found in ", source_dir,
         " — run script/fetch.sh first")
  }
  td <- tempfile("GSE141910_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  members <- untar(tarfile, list = TRUE)
  members <- members[grepl("\\.csv\\.gz$", members)]
  untar(tarfile, files = members, exdir = td, tar = "internal")
  files <- sort(list.files(td, pattern = "\\.csv\\.gz$", full.names = TRUE))
  if (length(files) != nrow(GSE141910_SAMPLES)) {
    stop("expected ", nrow(GSE141910_SAMPLES), " per-sample CSVs, found ", length(files))
  }
  titles <- sub("^GSM[0-9]+_", "", sub("\\.csv\\.gz$", "", basename(files)))
  idx <- match(titles, GSE141910_SAMPLES$sample)
  if (anyNA(idx)) stop("unlabelled sample(s): ", paste(titles[is.na(idx)], collapse = ", "))

  mat <- NULL
  ensg <- NULL
  for (i in seq_along(files)) {
    dt <- data.table::fread(files[i], header = TRUE, data.table = FALSE,
                            showProgress = FALSE)
    if (ncol(dt) != 2L) stop("unexpected column count in ", files[i], ": ", ncol(dt))
    if (is.null(ensg)) {
      ensg <- as.character(dt[[1]])
      mat <- matrix(NA_real_, nrow = length(ensg), ncol = length(files),
                    dimnames = list(ensg, titles))
    } else if (!identical(as.character(dt[[1]]), ensg)) {
      stop("gene order differs in ", files[i])
    }
    mat[, i] <- as.numeric(dt[[2]])
  }
  list(expr = mat, sample = titles, group = GSE141910_SAMPLES$group[idx],
       etiology = GSE141910_SAMPLES$etiology[idx], ensg = ensg)
}

# Collapse an Ensembl-keyed matrix to one row per symbol: when several Ensembl
# genes map to the same HGNC symbol, keep the row with the highest mean
# expression (standard, deterministic collapse; documented in the *_deg.json).
collapse_to_symbols <- function(expr, symbols) {
  ord <- order(symbols, -rowMeans(expr), seq_along(symbols))
  keep <- !duplicated(symbols[ord])
  ord <- sort(ord[keep], method = "radix")
  out <- expr[ord, , drop = FALSE]
  rownames(out) <- symbols[ord]
  out
}
