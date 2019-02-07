#############################
# Author: Lynsey Kovar
# makeGeneMap.R
# Takes bedtools intersect generated file mapping db Blast hits to genes and Blast output of hits of query genic regions to db genome to generate a map file containing putative identical genes between the two genomes.
#############################


library(argparse)
parser=ArgumentParser()
parser$add_argument("-i", "--bedtoolsFile", help="input blast output file (generated using option outfmt=6)")
parser$add_argument("-b", "--blastFile", help="input blast_out file (should be format 6)")
parser$add_argument("-o", "--outputFile", help="output map file name")
args=parser$parse_args()
cat('Generating genic map file between genomes',"\n")

inputBedtools <- read.table("~/Dropbox/Bladderwort/3_Analysis/u.gibba_NEW.blastHits2genes.txt", header=F, sep="\t")
inputBlast <- read.table("~/Dropbox/Bladderwort/3_Analysis/u.gibba_OLD.u.gibba_NEW.blastOut.txt", header=F, sep="\t")

#remove hits that did not fall within genic regions
#bedtools coded these as -1 for V8, V9 and V11

inputBedtools <- subset(inputBedtools, V8 != -1)

#found 22445 genes with hits