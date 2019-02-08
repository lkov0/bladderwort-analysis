#############################
# Author: Lynsey Kovar
# makeGeneMap.R
# Takes bedtools intersect generated file mapping db Blast hits to genes and Blast output of hits of query genic regions to db genome to generate a map file containing putative identical genes between the two genomes.
#############################


library(argparse)
parser=ArgumentParser()
parser$add_argument("-i", "--bedtoolsFile", help="input bedtools file generated from bedtools intersect with db genome hits and db genome genic regions")
parser$add_argument("-b", "--blastFile", help="input blast_out file (should be format 6)")
parser$add_argument("-o", "--outputFile", help="output map file name")
args=parser$parse_args()
cat('Generating genic map file between genomes',"\n")

inputBedtools <- read.table(args$bedtoolsFile, header=F, sep="\t")
inputBlast <- read.table(args$blastFile, header=F, sep="\t")

#remove hits that did not fall within genic regions
#bedtools coded these as -1 for V8, V9 and V11

inputBedtools <- subset(inputBedtools, V8 != -1)

#found 22445 db genes with hits (new utricularia genome)

#make sure column V9 < column V10 in inputBlast so intervals will match between inputBlast and inputBedtools files

    for(i in 1:nrow(inputBlast)) {
        x = 0
        if(inputBlast[i,9] > inputBlast[i,10]) {
            x = inputBlast[i,9]
            inputBlast[i,9] = inputBlast[i,10]
            inputBlast[i,10] = x
        }
        
    }

inputBedtools$dbGeneHit <- paste(inputBedtools$V1, inputBedtools$V2, inputBedtools$V3, sep="")
inputBlast$dbGeneHit <- paste(inputBlast$V2, inputBlast$V9, inputBlast$V10, sep="")

names(inputBedtools)[10] <- "dbGeneName"
names(inputBlast)[1] <- "queryGeneName"


genomeMap <- merge(inputBedtools[,c(10,13)], inputBlast,  by="dbGeneHit")
inputBlast.DbGeneHits <- subset(inputBlast, dbGeneName != ".")
inputBlast.noDbGeneHits <- subset(inputBlast, dbGeneName == ".")

genomeMap <- genomeMap[,c("dbGeneName","queryGeneName")]

write.table(genomeMap, file=args$outputFile, row.names=F, col.names=T, quote=F, sep="\t")

