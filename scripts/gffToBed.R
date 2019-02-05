#############################
# Author: Lynsey Kovar
# gffToBed.R
# takes input GFF file and converts it to a bed file with the feature ID as the name column 
#############################


library(argparse)
parser=ArgumentParser()
parser$add_argument("-i", "--inputFile", help="input GFF file")
parser$add_argument("-o", "--outputFile", help="output bed file name")
args=parser$parse_args()
cat('Converting from GFF to bed format',"\n")

#read in GFF file
inputGFF <- read.table(args$inputFile, header=F, sep="\t")

#generate bed file, leaving out "ID=" for ID names
outputBed <- data.frame(chrom = inputGFF$V1, chromStart = inputGFF$V4, chromEnd = inputGFF$V5, name=gsub("ID=", "", inputGFF$V9), score=".", strand=inputGFF$V7)

#write bed file output
write.table(outputBed, args$outputFile, row.names=F, col.names=F, quote=F, sep="\t")