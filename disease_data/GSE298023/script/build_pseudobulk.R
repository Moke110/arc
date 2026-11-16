#!/usr/bin/env Rscript
# Author CM gate and heart identities; sparse sum before the shared DEG worker.
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(Matrix))
a <- commandArgs(TRUE)
if (length(a) != 2L) stop('usage: build_pseudobulk.R author.rds output-directory')
x <- readRDS(a[1])
meta <- slot(x, 'meta.data')
counts <- slot(slot(x, 'assays')[['RNA']], 'counts')
rm(x); invisible(gc()) # Release normalised/scaled assays before aggregation.
required <- c('ident','heart_number','condition','region','run_number')
if (!all(required %in% names(meta))) stop('missing author identity/gate columns')
if (!inherits(counts,'sparseMatrix')) stop('RNA counts must be sparse')
if (!identical(colnames(counts),rownames(meta))) stop('metadata/count barcode order differs')
if (any(!is.finite(counts@x)) || any(counts@x < 0) || any(counts@x != floor(counts@x))) stop('RNA counts are not nonnegative integers')
if (anyNA(meta[,required]) || any(!meta$condition %in% c('NF','DCM','ICM'))) stop('invalid author donor/condition')
if (any(tapply(as.character(meta$condition),meta$heart_number,function(z) length(unique(z))) != 1L)) stop('one heart has multiple diseases')
keep <- which(as.character(meta$ident) == 'Cardiomyocytes')
if (!length(keep)) stop('author CM gate selects no cells')
hearts <- sort(unique(as.character(meta$heart_number[keep])))
indicator <- sparseMatrix(i=keep,j=match(meta$heart_number[keep],hearts),x=1,dims=c(nrow(meta),length(hearts)))
pb <- as.matrix(counts %*% indicator)
colnames(pb) <- hearts
if (anyDuplicated(rownames(pb)) || anyNA(rownames(pb))) stop('invalid source gene axis')
donors <- do.call(rbind,lapply(hearts,function(h) {
 z <- meta[keep[as.character(meta$heart_number[keep])==h],,drop=FALSE]
 data.frame(sample=h,donor=h,group=as.character(z$condition[1]),cells_cm=nrow(z),
            batch=paste0('run_',paste(sort(unique(as.character(z$run_number))),collapse='_')),
            regions=paste(sort(unique(as.character(z$region))),collapse=';'))
}))
dir.create(a[2],recursive=TRUE,showWarnings=FALSE)
write.table(data.frame(source_gene_id=rownames(pb),gene_symbol=rownames(pb),pb,check.names=FALSE),file.path(a[2],'gse298023_cm_pseudobulk.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.csv(donors,file.path(a[2],'gse298023_donors.csv'),row.names=FALSE,quote=TRUE)
cat('author CM cells:',length(keep),'donors:',length(hearts),'genes:',nrow(pb),'\n');print(donors)
