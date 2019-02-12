#appendFPKM.R

library(argparse)
library(ggplot2)
parser=ArgumentParser()
parser$add_argument("-i", "--inputFile", help="input genePairs file")
parser$add_argument("-f", "--fpkmFile", help="fpkm file from cufflinks - genes.fpkm_tracking")
parser$add_argument("-p", "--plotName", help="name for boxplot")
parser$add_argument("-o", "--outputFile", help="output file name")
args=parser$parse_args()
cat('appending FPKM values to genePairs file.',"\n")

genePairs <- read.table(args$inputFile, header=T, sep="\t") 
fpkms <- read.table(args$fpkmFile, header=T, sep="\t") 

genePairs$tracking_id1 <-  gsub("ID=", "", genePairs$geneId1)
genePairs$tracking_id2 <- gsub("ID=", "", genePairs$geneId2)
fpkms$tracking_id1 <- fpkms$tracking_id
fpkms$tracking_id2 <- fpkms$tracking_id

genePairs.fpkms <- merge(genePairs, fpkms[,c(10,14)], by="tracking_id1")
names(genePairs.fpkms)[15] <- "FPKM1"
genePairs.fpkms <- merge(genePairs.fpkms, fpkms[,c(10,15)], by="tracking_id2")
names(genePairs.fpkms)[16] <- "FPKM2"

genePairs.fpkms$foldChange <- 0

for(i in 1:nrow(genePairs.fpkms)) {
    if(genePairs.fpkms[i,"FPKM1"] == 0 | genePairs.fpkms[i,"FPKM2"] == 0) {
        genePairs.fpkms[i,"foldChange"] = NA
    }
    else {
        if(genePairs.fpkms[i,"FPKM1"] > genePairs.fpkms[i,"FPKM2"]) {
            genePairs.fpkms[i,"foldChange"] = (genePairs.fpkms[i,"FPKM1"] - genePairs.fpkms[i,"FPKM2"]) / genePairs.fpkms[i,"FPKM2"]
        }
        if(genePairs.fpkms[i,"FPKM1"] < genePairs.fpkms[i,"FPKM2"]) {
            genePairs.fpkms[i, "foldChange"] = (genePairs.fpkms[i,"FPKM2"] - genePairs.fpkms[i,"FPKM1"]) / genePairs.fpkms[i,"FPKM1"]
        }
    }
}

genePairs.fpkms <- subset(genePairs.fpkms, !(is.na(foldChange)))

#write table containing only expressed gene pairs
write.table(genePairs.fpkms, file=args$outputFile, row.names=F, col.names=T, quote=F, sep="\t")

png(args$plotName)
boxplot(genePairs.fpkms$foldChange ~ genePairs.fpkms$orientation, ylim=c(0,1000), main="expressed gene pairs - fold change", ylab="Fold change", xlab="orientation")
dev.off()
