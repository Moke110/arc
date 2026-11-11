#!/usr/bin/env Rscript
# Fit one raw-count contrast at a time, streaming its inputs and results.
suppressPackageStartupMessages(library(DESeq2))

fit_contrast <- function(counts, samples, design) {
    warnings <- character()
    messages <- character()
    if (any(colSums(counts) == 0)) {
        return(list(reason = "zero_library", warnings = "", messages = ""))
    }
    if (nrow(samples) <= ncol(model.matrix(design, samples))) {
        return(list(reason = "no_residual_degrees_of_freedom", warnings = "", messages = ""))
    }
    size_factors <- if (any(rowSums(counts == 0) == 0)) "ratio" else "poscounts"
    result <- tryCatch(
        withCallingHandlers({
            dataset <- DESeqDataSetFromMatrix(counts, samples, design)
            dataset <- DESeq(dataset, test = "Wald", fitType = "parametric",
                             sfType = size_factors, quiet = TRUE, parallel = FALSE)
            contrast <- results(dataset, contrast = c("condition", "treated", "control"),
                                lfcThreshold = 0, altHypothesis = "greaterAbs",
                                independentFiltering = FALSE)
            list(reason = "", values = list(contrast$log2FoldChange, contrast$pvalue,
                                            p.adjust(contrast$pvalue, method = "BH", n = nrow(counts))))
        }, warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        }, message = function(m) {
            messages <<- c(messages, conditionMessage(m))
            invokeRestart("muffleMessage")
        }),
        error = function(e) list(reason = paste0("DESeq2: ", conditionMessage(e)))
    )
    if (!is.null(result$values) && !any(is.finite(result$values[[2]]))) {
        result$reason <- "no_finite_pvalues"
    }
    result$size_factors <- size_factors
    result$warnings <- paste(unique(warnings), collapse = " | ")
    result$messages <- paste(unique(messages), collapse = " | ")
    result
}

main <- function(args) {
    directory <- args[[1]]
    n_genes <- as.integer(args[[2]])
    groups <- read.csv(file.path(directory, "groups.csv"), stringsAsFactors = FALSE)
    samples <- read.csv(file.path(directory, "samples.csv"), stringsAsFactors = FALSE)
    counts_file <- file(file.path(directory, "counts.bin"), "rb")
    names <- c("log2_fold_change", "pvalue", "padj")
    outputs <- lapply(file.path(directory, paste0(names, ".bin")), file, open = "wb")
    on.exit({
        close(counts_file)
        for (output in outputs) close(output)
    })
    reports <- vector("list", nrow(groups))
    offset <- 0L
    for (i in seq_len(nrow(groups))) {
        group <- groups[i, ]
        n_samples <- group$n_controls + group$n_treated
        metadata <- samples[offset + seq_len(n_samples), ]
        offset <- offset + n_samples
        metadata$plate <- factor(metadata$plate)
        metadata$condition <- factor(metadata$condition, levels = c("control", "treated"))
        design <- if (group$n_plates > 1L) ~ plate + condition else ~ condition
        expected <- n_genes * n_samples
        values <- readBin(counts_file, integer(), n = expected, size = 4L, endian = "little")
        stopifnot(length(values) == expected, !anyNA(values))
        fit <- if (group$n_treated < 2L || group$n_controls < 2L) {
            list(reason = "insufficient_replication", warnings = "", messages = "")
        } else {
            fit_contrast(matrix(values, nrow = n_genes), metadata, design)
        }
        if (is.null(fit$values)) fit$values <- rep(list(rep(NA_real_, n_genes)), 3L)
        for (j in seq_along(outputs)) {
            writeBin(as.double(fit$values[[j]]), outputs[[j]], size = 8L, endian = "little")
        }
        n_testable <- sum(is.finite(fit$values[[2]]))
        reports[[i]] <- data.frame(
            status = if (n_testable > 0L) "eligible" else "untestable",
            n_genes_testable = n_testable,
            n_significant = sum(fit$values[[3]] < 0.05, na.rm = TRUE),
            size_factors = if (is.null(fit$size_factors)) "" else fit$size_factors,
            reason = fit$reason,
            warnings = fit$warnings,
            messages = fit$messages
        )
        if (i %% 25L == 0L || i == nrow(groups)) {
            cat(sprintf("DESeq2 contrasts %d/%d\n", i, nrow(groups)))
            flush.console()
        }
    }
    stopifnot(offset == nrow(samples), length(readBin(counts_file, raw(), n = 1L)) == 0L)
    write.csv(do.call(rbind, reports), file.path(directory, "results.csv"), row.names = FALSE)
    writeLines(c(as.character(getRversion()), as.character(packageVersion("DESeq2"))),
               file.path(directory, "versions.txt"))
}

main(commandArgs(trailingOnly = TRUE))
