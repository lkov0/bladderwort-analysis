#############################
# Author: Lynsey Kovar
# makeIntergenicBed.R
# Makes a bedfile of coordinates for intergenic regions of parallel, convergent, and divergent genes. Input is generated from assignConvergentDivergentParallel.R.
#############################

# input: files generated from assignConvergentDivergentParallel.R

library(argparse)
parser=ArgumentParser()
parser$add_argument("-i", "--inputFile", help="input text file generated from assignConvergentDivergentParallel.R")
parser$add_argument("-o", "--outputStem", help="output file stem")
args=parser$parse_args()
cat('Making bedfiles for intergenic regions of convergent, divergent, and parallel gene pairs.',"\n")

genePairs <- read.table(args$inputFile, header=T, sep="\t")

genePairs.bed <- data.frame(chr = genePairs$chr1, start = genePairs$stop1, stop = genePairs$start2, name = genePairs$geneId1)

genePairs.bed$length <- genePairs.bed$stop - genePairs.bed$start

genePairs.bed <- subset(genePairs.bed, length > 10)
genePairs.bed$length <- NULL

write.table(genePairs.bed, paste(args$outputStem, ".bed", sep=""), row.names=F, col.names=F, sep="\t", quote=F)

### make it so everything below here happens if an additional option (splitOrientation == TRUE)

# genePairs <- read.table(args$inputFile, header=T, sep="\t")
# 
# genePairs.parallel <- subset(genePairs, orientation == "parallel")
# genePairs.convergent <- subset(genePairs, orientation == "convergent")
# genePairs.divergent <- subset(genePairs, orientation == "divergent")
#     
# genePairs.parallelIntergenic <- data.frame(chr = genePairs.parallel$chr1, start = genePairs.parallel$stop1, stop = genePairs.parallel$start2, name=paste(genePairs.parallel$geneId1, genePairs.parallel$geneId2, sep="."))
# genePairs.convergentIntergenic <- data.frame(chr = genePairs.convergent$chr1, start = genePairs.convergent$stop1, stop = genePairs.convergent$start2, name=paste(genePairs.convergent$geneId1, genePairs.convergent$geneId2, sep="."))
# genePairs.divergentIntergenic <- data.frame(chr = genePairs.divergent$chr1, start = genePairs.divergent$stop1, stop = genePairs.divergent$start2, name=paste(genePairs.divergent$geneId1, genePairs.divergent$geneId2, sep="."))
# write.table(genePairs.parallelIntergenic, paste(args$outputStem, ".intergenic.parallel.bed", sep=""), row.names=F, col.names=F, sep="\t", quote=F)
# write.table(genePairs.convergentIntergenic, paste(args$outputStem, ".intergenic.convergent.bed", sep=""), row.names=F, col.names=F, sep="\t", quote=F)
# write.table(genePairs.divergentIntergenic, paste(args$outputStem, ".intergenic.divergent.bed", sep=""), row.names=F, col.names=F, sep="\t", quote=F)


