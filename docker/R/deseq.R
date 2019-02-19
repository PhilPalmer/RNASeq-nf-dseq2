#!/opt/conda/envs/RNASeq-nf-salmon-dseq2/bin/Rscript

library("tximport")
library("readr")
library("tximportData")
library("DESeq2")

dir <- system.file("extdata", package="tximportData")
samples <- read.table(args[1], header=TRUE, sep=",")

tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))

files <- file.path(".",samples$samples, "quant.sf")
names(files) <- samples$samples

txi <- tximport(files, type="salmon", tx2gene=tx2gene, ignoreAfterBar=TRUE)

condition <- samples$condition

ddsTxi <- DESeqDataSetFromTximport(txi,
                                colData = samples,
                                design = ~ condition)

save(ddsTxi, file="deseq.RData")