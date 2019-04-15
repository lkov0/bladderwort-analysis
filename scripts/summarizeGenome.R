#summarizeGenome.R

# for file in $(ls); do
# grep -v "#" $file | awk '{print $1, $2, $3, $4, $5}' | sed "s/ /\t/g" >  $file.simple
# done

library(argparse)
parser=ArgumentParser()
parser$add_argument("-g", "--gffFile", help="gff file with only the first five columns and no comment lines")
args=parser$parse_args()

genome <- read.table(args$gffFile, sep="\t")

genome_chr <- subset(genome, V3 == "chromosome" & (V1 != "Mt") & (V1 != "Pt"))
genome_genes <- subset(genome, (V1 != "Mt") & (V1 != "Pt") & (V3 == "gene"))

genome_genes$length <- genome_genes$V5 - genome_genes$V4

genome_length <- sum(genome_chr$V5)
genes_length <- sum(genome_genes$length)

intergenic_length <- genome_length - genes_length
intergenic_length <- intergenic_length / 1000000

cat("\n", args$gffFile, ":\n")
cat("intergenic sequence (Mb): ", intergenic_length, "\n\n")
