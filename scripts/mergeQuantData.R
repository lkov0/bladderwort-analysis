#mergeQuantData.R

library(vegan)
library(ggplot2)
library(DESeq2)

sampleData <- read.table("sampledata.txt", header=T, sep="\t", row.names = 1)
bladder1 <- read.table("1B_htSeq_counts.txt", header=F, sep="\t", row.names=1)
bladder2 <- read.table("2B_htSeq_counts.txt", header=F, sep="\t", row.names=1)
bladder3 <- read.table("3B_htSeq_counts.txt", header=F, sep="\t", row.names=1)
leaf1 <- read.table("1L_htSeq_counts.txt", header = F, sep = "\t", row.names=1)
leaf2 <- read.table("2L_htSeq_counts.txt", header = F, sep = "\t", row.names=1)
leaf3 <- read.table("3L_htSeq_counts.txt", header = F, sep = "\t", row.names=1)
stem1 <- read.table("1S_htSeq_counts.txt", header = F, sep = "\t", row.names=1)
stem2 <- read.table("2S_htSeq_counts.txt", header = F, sep = "\t", row.names=1)
stem3 <- read.table("3S_htSeq_counts.txt", header = F, sep = "\t", row.names=1)
rhizoid1 <- read.table("1R_htSeq_counts.txt", header = F, sep = "\t", row.names=1)
rhizoid2 <- read.table("2R_htSeq_counts.txt", header = F, sep = "\t", row.names=1)
rhizoid3 <- read.table("3R_htSeq_counts.txt", header = F, sep = "\t", row.names=1)

bladder1$V3 <- NULL

names(bladder1) <- ("bladder1")
names(bladder2) <- ("bladder2")
names(bladder3) <- ("bladder3")
names(leaf1) <- ("leaf1")
names(leaf2) <- ("leaf2")
names(leaf3) <- ("leaf3")
names(stem1) <- ("stem1")
names(stem2) <- ("stem2")
names(stem3) <- ("stem3")
names(rhizoid1) <- ("rhizoid1")
names(rhizoid2) <- ("rhizoid2")
names(rhizoid3) <- ("rhizoid3")

expr <- cbind(bladder1, bladder2, bladder3, leaf1, leaf2, leaf3, stem1, stem2, stem3, rhizoid1, rhizoid2, rhizoid3)
expr <- expr[1:29666,]

#using DESeq2 for normalization (vst) 

expr.dds <- DESeqDataSetFromMatrix(countData = expr, colData = sampleData, design = ~ condition)
expr.normalized <- counts(expr.dds.nb, normalized=TRUE)
expr.vsd <- vst(expr.dds, blind = FALSE)
pcaData <- plotPCA(expr.vsd, intgroup=c("condition", "type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  ggtitle("Bladderwort 3' Data - PCA") +
  theme_classic() + scale_color_manual(values=c("black", "red", "blue", "green4")) +
  coord_fixed()

ggsave("~/Dropbox/Figures/bladderwort_3primePCA.png", height=5, width=6)

pcaData.nb <- plotPCA(counts(expr.dds.nb, normalized=TRUE), intgroup=c("condition", "type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData.nb, "percentVar"))
ggplot(pcaData.nb, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  ggtitle("Bladderwort 3' Data - PCA") +
  theme_classic() + scale_color_manual(values=c("black", "red", "blue", "green4")) +
  coord_fixed()


