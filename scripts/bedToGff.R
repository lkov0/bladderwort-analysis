#bedToGff.R

library(argparse)
parser=ArgumentParser()
parser$add_argument("-i", "--inputFile", help="input bed file")
parser$add_argument("-o", "--outputFile", help="output bed file name")
parser$add_argument("-f", "--feature", help="type of feature (gene, CDS, etc.)")
parser$add_argument("-p", "--prefix", help="prefix for gene ids")
parser$add_argument("-s", "--source", help="source of regions generated")
args=parser$parse_args()
cat('Converting bed file to GFF',"\n")

bedFile <- read.table(args$inputFile)
gff <- data.frame(chr = bedFile$V1, source = args$source, feature = args$feature, start = bedFile$V2, end = bedFile$V3, score = ".", strand = bedFile$V6, frame = ".", attribute= paste(args$feature, "_id=", args$prefix, row.names(bedFile), sep=""))

write.table(gff, args$outputFile, row.names=F, col.names=F, sep="\t", quote=F)
