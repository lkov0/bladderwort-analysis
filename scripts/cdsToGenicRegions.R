#############################
# Author: Lynsey Kovar
# cdsToGenicRegions.R
# Converts CDS regions to genic regions associated with the coge_fid field in the input gff file. This should be the last attribute in the attribute column.
#############################


library(argparse)
library(dplyr)
library(tidyr)
parser=ArgumentParser()
parser$add_argument("-i", "--inputFile", help="input gff file")
parser$add_argument("-o", "--outputFile", help="Output gff file name")
args=parser$parse_args()
cat('Converting CDS regions to genic regions',"\n")

gffFile <- read.table(args$inputFile, header=F, sep="\t")

#fix names in gff file since there is no header
names(gffFile) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
gffFile <- subset(gffFile, feature == "CDS")

cat("there are ", nrow(gffFile), " CDS regions.\n")

#parse Info column into separate columns using ';' as the delimiter
gffFile.parsed <- gffFile %>% separate(attribute, c("attribute1","coge_fid"), sep="coge_fid=")

#clean up memory
rm(gffFile)
gffFile.parsed$attribute1=NULL
gc()


#find the minimum and maximum positions for genes using their CDS coordinates, structuring new dataframe in gff format
for(gene in unique(gffFile.parsed$coge_fid)) {

    if(gene == unique(gffFile.parsed$coge_fid)[1]) {
        geneInfo <- subset(gffFile.parsed, coge_fid==gene)
        gffFile.new <- data.frame(seqname=geneInfo[1,1], source=geneInfo[1,2], feature="gene", start=min(geneInfo$start), end=max(geneInfo$end), score=geneInfo[1,6], strand=geneInfo[1,7], frame=geneInfo[1,8], attribute=gene)
    }
    
    else {
        geneInfo <- subset(gffFile.parsed, coge_fid==gene)
        gffFile.toAppend <- data.frame(seqname=geneInfo[1,1], source=geneInfo[1,2], feature="gene", start=min(geneInfo$start), end=max(geneInfo$end), score=geneInfo[1,6], strand=geneInfo[1,7], frame=geneInfo[1,8], attribute=gene)
        gffFile.new <- rbind(gffFile.new, gffFile.toAppend)
    }
    
}

#fix names to match bladderwort contig names in genome file, also add ID= to attribute column so that it can be parsed by cufflinks
gffFile.new$seqname <- paste("lcl|", gffFile.new$seqname, sep="")
gffFile.new$attribute <- paste("ID=", gffFile.new$attribute, sep="")

cat("there are ", nrow(gffFile.new), " genes.")

#write new gff file
write.table(gffFile.new, args$outputFile, row.names=F, col.names=F, sep="\t", quote=F)
