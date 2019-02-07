#############################
# Author: Lynsey Kovar
# blastToBed.R
# Takes input blast file and puts DB regions that were hit in bed format
#############################


library(argparse)
parser=ArgumentParser()
parser$add_argument("-i", "--inputFile", help="input blast output file (generated using option outfmt=6)")
parser$add_argument("-o", "--outputFile", help="output bed file name")
args=parser$parse_args()
cat('Getting DB hits from blast file -> bed file',"\n")

blastOut <- read.table(args$inputFile, header=F, sep="\t")
bedOut <- data.frame(chrom=blastOut$V2, chromStart=blastOut$V9, chromEnd=blastOut$V10, name=".", score=".")
bedOut$strand <- ifelse(bedOut$chromStart > bedOut$chromEnd, "-", "+")
    
    for(i in 1:nrow(bedOut)) {
        x = 0
        if(bedOut[i,2] > bedOut[i,3]) {
            x = bedOut[i,2]
            bedOut[i,2] = bedOut[i,3]
            bedOut[i,3] = x
        }
        
    }
    
write.table(bedOut, args$outputFile, row.names=F, col.names=F, sep="\t", quote=F)