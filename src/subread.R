#!/usr/bin/env Rscript
library(Rsubread)

args <- commandArgs(trailingOnly=TRUE)


folder_name <- args[1]
filestem <- args[2]

files <- list.files(folder_name,
                    pattern=paste0(filestem, '.r[12].fastq(.gz)+')
)

print(files)

quit(status=0)


fq1 <- paste0(folder_name, '/', filestem, '.r1.fastq.gz')
fq2 <- paste0(folder_name, '/', filestem, '.r2.fastq.gz')


# barcodes_fn <- 'barcodes-cut.tsv'
subjunc(
    index='/gpfs/gsfs9/users/wellerca/scrpio/db/hg38',
    readfile1='01_120405.r1.fastq.gz',
    readfile2='01_120405.r2.fastq.gz',
    input_format='gzFASTQ',
    output_format='BAM',
    useAnnotation = TRUE,
    annot.inbuilt = 'hg38',
    GTF.featureType = "exon",
    GTF.attrType = "gene_id",
    # number of threads
    nthreads = 2
    # other parameters passed to align, subjunc and featureCounts functions 
)

saveRDS(o, file='UMARY-1544.RDS')



quit(status=0)

# Generate sparse matrix for seurat
library(data.table)
library(Matrix)
library(Rsubread)

o <- readRDS('UMARY-1544.RDS')

features <- fread('features.tsv.gz', col.names=c('ENSEMBL','SYMBOL','FEATURE','CHR','START','STOP'))
to_merge <- copy(features[, .SD, .SDcols='ENSEMBL'])


counts <- as.data.table(o$counts[[1]])
cell_umis <- colnames(counts)


counts[, 'ENSEMBL' := o$annotation$GeneID]
counts <- merge(to_merge, counts, by='ENSEMBL', all=TRUE)



fillNA <- function(DT, sdcols, x) {
    # replaces all NA values with x in data.table DT
    for (j in sdcols)
        set(DT,which(is.na(DT[[j]])),j,x)
}

fillNA(counts, cell_umis, 0)

# Reorder to `features` table order
setkey(counts, 'ENSEMBL')
counts <- counts[.(features$ENSEMBL)]

counts[,'ENSEMBL' := NULL]

barcodes <- colnames(counts)
barcodes <- paste0(barcodes, '-1')
writeLines(barcodes, con='barcodes.tsv')
system('gzip barcodes.tsv')




A <- as(as.matrix(counts), "sparseMatrix") 
writeMM(A, file='matrix.mtx')
system('gzip matrix.mtx')


# library("org.Hs.eg.db")

# get_ids <- function(entrez_list) {
#     out <- select(org.Hs.eg.db, 
#         keys = entrez_list,
#         columns = c("ENSEMBL", "ENTREZID", "SYMBOL"),
#         keytype = "ENSEMBL")
#     return(out)
# }
