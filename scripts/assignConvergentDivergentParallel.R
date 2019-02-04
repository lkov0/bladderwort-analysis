#############################
# Author: Lynsey Kovar
# assignConvergentDivergentParallel.R
# takes input GFF file containing only gene annotations and pulls out adjacent genes, assigns them as convergent, divergent, or parrallel to each other, and gets length of intergenic region between them. 
#############################

library(argparse)
library(ggplot2)
parser=ArgumentParser()
parser$add_argument("-i", "--inputFile", help="input GFF file")
parser$add_argument("-o", "--outputFile", help="output file name")
args=parser$parse_args()
cat('Assigning divergent convergent parallel gene pairs.',"\n")

#load GFF file
gffFile <- read.table(args$inputFile, header=F, sep="\t")

#order GFF file by chromosome and start coordinates of gene
gffFile <- gffFile[order(gffFile$V1, gffFile$V4),]


#find gene pairs. first checking if they are on the same chromosome.
for(i in 1:(nrow(gffFile)-1)) {

    if(gffFile[i,1] == gffFile[i+1,1]) {
        if(!(exists("genePairs"))) {
            genePairs <- data.frame(chr1 = gffFile[i,1], start1=gffFile[i,4], stop1=gffFile[i,5], strand1=gffFile[i,7], geneId1=gffFile[i,9], chr2 = gffFile[i+1,1], start2=gffFile[i+1,4], stop2=gffFile[i+1,5], strand2=gffFile[i+1,7], geneId2=gffFile[i+1,9])
        }
        else {
            genePairs.toMerge <- data.frame(chr1 = gffFile[i,1], start1=gffFile[i,4], stop1=gffFile[i,5], strand1=gffFile[i,7], geneId1=gffFile[i,9], chr2 = gffFile[i+1,1], start2=gffFile[i+1,4], stop2=gffFile[i+1,5], strand2=gffFile[i+1,7], geneId2=gffFile[i+1,9])
            genePairs <- rbind(genePairs, genePairs.toMerge)
        }
    }
}

#change strand1 and strand2 to character vectors
genePairs$strand1 <- as.character(genePairs$strand1)
genePairs$strand2 <- as.character(genePairs$strand2)

#make orientation column
genePairs$orientation <- "blank"

#assign orientation classification based on strand of each gene in a gene pair
for(i in 1:nrow(genePairs)) {
    if((genePairs[i,4] == "+" & genePairs[i,9] == "+") | (genePairs[i,4] == "-" & genePairs[i,9] == "-")) {
        genePairs[i,11] = "parallel"
    }
    if(genePairs[i,4] == "+" & genePairs[i,9] == "-") {
        genePairs[i,11] = "convergent"
    }
    if(genePairs[i,4] == "-" & genePairs[i,9] == "+") {
        genePairs[i,11] = "divergent"
    }
}

#TODO: add distance between gene pairs



cat(table(genePairs$orientation))

write.table(genePairs, args$outputFile, row.names=F, col.names=T, quote=F, sep="\t")
        