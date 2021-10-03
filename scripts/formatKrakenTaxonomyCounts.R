# formatKrakenTaxonomyCounts.R
# Purpose: to format kraken output files for downstream taxonomy-based analyses, with counts of sequences assigned to each taxonomic assignment

library(stringr)
library(argparse)
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-i", "--inputFile", help="Input kraken report - generated with flag --use-mpa-style")
parser$add_argument("-o", "--outputFile", help="Output file name")
args=parser$parse_args()
cat('Formatting kraken output for downstream taxonomy-based analyses',"\n")

taxonomy <- read.delim(args$inputFile, header=FALSE, sep="\t")
taxonomy$V1 <- as.character(taxonomy$V1)

d <- str_extract(taxonomy$V1, "d__[^|]+")
k <- str_extract(taxonomy$V1, "k__[^|]+")
p <- str_extract(taxonomy$V1, "p__[^|]+")
c <- str_extract(taxonomy$V1, "c__[^|]+")
o <- str_extract(taxonomy$V1, "o__[^|]+")
f <- str_extract(taxonomy$V1, "f__[^|]+")
g <- str_extract(taxonomy$V1, "g__[^|]+")
s <- str_extract(taxonomy$V1, "s__[^|]+")

taxonomy.reformatted <- data.frame( readCount = taxonomy$V2, domain = d, kingdom = k, phylum = p, class = c, order = o, family = f, genus = g, species = s )

write.table(taxonomy.reformatted, args$outputFile, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
