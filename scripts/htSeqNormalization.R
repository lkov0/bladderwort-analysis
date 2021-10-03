#htSeqNormalization.R
# get normalized_counts.txt file from DESeq2 normalized table.

library(DESeq2)

directory <- "~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out"

sampleFiles <- grep("ALL.includingPacBio.txt",list.files(directory),value=TRUE)
sampleCondition <- sub("_htSeq_counts.ALL.includingPacBio.txt","",sampleFiles)
sampleCondition <- sub("[0-9]","",sampleCondition)

sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)

library("DESeq2")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)

#get normalized counts
dds_normalized <- estimateSizeFactors(ddsHTSeq)
normalized_counts <- counts(dds_normalized, normalized=TRUE)
colnames(normalized_counts) <- gsub("_htSeq_counts.ALL.includingPacBio.txt","_geneCounts",colnames(normalized_counts))

names(normalized_counts) <- sub("_htSeq_counts.ALL.includingPacBio.txt","_geneCounts",names(normalized_counts))
write.table(normalized_counts, "~/Xfer/Bladderwort_3pri/3_Quantification/normalized_counts.txt", col.names=T, row.names=T, quote=F, sep="\t")
