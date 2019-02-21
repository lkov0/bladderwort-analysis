#script to get average coverage of genomic ranges given a file generated from bedtools genome cov with -d option and a file containing coordinites to append to. 

library(argparse)
parser=ArgumentParser()
parser$add_argument("-c", "--coverageFile", help="Output from bedtools genomecov -d. Lists coverage at each genomic position")
parser$add_argument("-i", "--inputFile", help="Input gene pairs file generated from appendFPKM.R")
parser$add_argument("-o", "--outputFile", help="output file name")
args=parser$parse_args()
cat("finding gene pairs shared between ", args$dbFile, " and ", args$queryFile, "\n")

genomeCov <- read.table("~/Dropbox/Bladderwort/3_Analysis/Coverage/SRR768657_U.gibba_NEW.genomecov.txt", header=T, sep="\t")
genePairs <- read.table("~/Dropbox/Bladderwort/3_Analysis/SRR768657_u.gibba_NEW.genePairs.foldChange.ALL.txt", header=T, sep="\t")