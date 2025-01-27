library(tximport)
library(data.table)
library(dplyr)
library(readr)

# thank our noodly lord for Love
# http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html

# working dir, biowulf2
working_dir <- '/data/bryanjm/rawData/unified_gene_expression/salmon_counts'
# eyeMac
#working_dir <- '/Volumes/ThunderBay/PROJECTS/mcgaughey/unified_gene_expression/salmon_counts'

setwd(working_dir)

# pull in salmon files
files <- list.files(path=working_dir,recursive=TRUE,pattern='quant.sf')

# Gene TX to name conversion
load('/data/bryanjm/rawData/gencode_v25_annotation.Rdata')
anno <- gencode_v25_annotation %>% select(Transcript.ID, Gene.Name)

# merge tx specific counts to gene level and scale TPM
txi <- tximport(files, type = "salmon", tx2gene = anno, reader = read_tsv, countsFromAbundance = c("no"))

counts <- data.frame(txi$counts)
names <- sapply(files, function(x) strsplit(x,"\\/")[[1]][1])
colnames(counts) <- names

save(counts, file = "/data/bryanjm/rawData/unified_gene_expression/counts_processed.Rdata")