# intergenic region analysis

# first, check orientation of candidates
library(ggplot2)


geneHits <- read.table("~/Xfer/jwlab/Bladderwort_3pri/1_Assembly/scaffolds_lk_anno_coordinates_in_PacBio.bed", stringsAsFactors = F)
head(geneHits)

uni_cand <- read.table("~/Xfer/jwlab/Bladderwort_3pri/4_CandidateMining/uni_cand.3pri.summary.tsv", header = T, stringsAsFactors = F)
bid_cand <- read.table("~/Xfer/jwlab/Bladderwort_3pri/4_CandidateMining/bid_cand.3pri.summary.tsv", header = T, stringsAsFactors = F)
ter_cand <- read.table("~/Xfer/jwlab/Bladderwort_3pri/4_CandidateMining/term_cand.3pri.summary.tsv", header = T, stringsAsFactors = F)
ins_cand <- read.table("~/Xfer/jwlab/Bladderwort_3pri/4_CandidateMining/ins_cand.3pri.summary.tsv", header = T, stringsAsFactors = F)

isOnSameChrom <- function(bedfile, candfile) {
  
  output <- data.frame(gene1 = character(nrow(candfile)), gene2 = character(nrow(candfile)), onSameChrom = logical(nrow(candfile)), distInPacBio = numeric(nrow(candfile)), distInMyGenome = numeric(nrow(candfile)), pairId = character(nrow(candfile)), stringsAsFactors = F)
  
  for(i in 1:nrow(candfile)) {
    
    output[i,1] = as.character(candfile[i, "geneId1"])
    output[i,2] = as.character(candfile[i, "geneId2"])
    output[i,5] = candfile[i, "distance"]
    output[i,6] = candfile[i, "candidate_name"]
    
    if(!(candfile[i, "geneId1"] %in% bedfile$V4) | !(candfile[i, "geneId2"] %in% bedfile$V4)) {
      output[i,3] = FALSE
      output[i,4] = NA
    }
    else {
      if((bedfile[bedfile$V4 == candfile[i, "geneId1"], 1] == bedfile[bedfile$V4 == candfile[i, "geneId2"], 1])) {
        output[i,3] = TRUE
        if(bedfile[bedfile$V4 == candfile[i, "geneId1"], 2] < bedfile[bedfile$V4 == candfile[i, "geneId2"], 2]) {
          output[i,4] = bedfile[bedfile$V4 == candfile[i, "geneId2"], 2] - bedfile[bedfile$V4 == candfile[i, "geneId1"], 3]
        }
        else {
          output[i,4] = bedfile[bedfile$V4 == candfile[i, "geneId1"], 2] - bedfile[bedfile$V4 == candfile[i, "geneId2"], 3]
        }
      }
      else {
        output[i,3] = FALSE
        output[i,4] = NA
      } 
    }
  }
  
  return(output)
  
}

uni_info <- isOnSameChrom(geneHits, uni_cand)
bid_info <- isOnSameChrom(geneHits, bid_cand)
ter_info <- isOnSameChrom(geneHits, ter_cand)
ins_info <- isOnSameChrom(geneHits, ins_cand)

ter_info$type <- "terminator"
ins_info$type <- "insulator"
bid_info$type <- "bidirectional promoter"
uni_info$type <- "unidirectional promoter"
all_info <- rbind(uni_info, bid_info, ter_info, ins_info)

ggplot(all_info, aes(x = distInPacBio, y = distInMyGenome, col = type)) + 
  geom_point(size = 2) + 
  theme_bw() + 
  xlab("PacBio genome distance (bp)") + 
  ylab("My genome distance (bp)") + 
  ggtitle("Intergenic region length comparison\nMy genome vs. published PacBio") + 
  theme(axis.text = element_text(size = 14, color  = "black"), axis.title = element_text(size = 14), plot.title = element_text(size = 16), legend.title = element_blank(), legend.text = element_text(size = 12)) + 
  scale_color_manual(values = c("red", "blue", "green", "black"))

# plot difference in lengh - pacBio genome vs. my genome for each type of regulatory sequence
ggplot(all_info, aes(x = type, y = differenceInLength, col = type)) + 
  geom_boxplot() +
  theme_bw() + 
  xlab("Element type") + 
  ylab("PacBio genome difference in length (bp)") + 
  ggtitle("Intergenic sequence length difference vs. PacBio genome") + 
  theme(axis.text = element_text(size = 14, color  = "black"), axis.title = element_text(size = 14), plot.title = element_text(size = 16), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  scale_color_manual(values = c("red", "blue", "green", "black"))

ggsave("~/Dropbox/Bladderwort/0_Plots/intergenic_length_comparison_boxplots.png")
head(genePairs)

# compare amount and length of each type of sequence
genePairs_myGenome <- read.table("~/Xfer-work/genomes/utricularia/scaffolds.ugibba_lk.ALL.includingPacBio.genePairs.txt", header=T, sep="\t")
genePairs_pbio <- read.table("~/Xfer-work/genomes/utricularia/u.gibba_NEW.genePairs.txt", header=T, sep="\t")

genePairs_pbio$genome = "PacBio"
genePairs_myGenome$genome = "myGenome"

genePairs_all <- rbind(genePairs_myGenome, genePairs_pbio)
genePairs_all_lt10000 <- subset(genePairs_all, distance < 10000)
goodDists_lt10000 <- subset(genePairs_all_lt10000, distance > 0)

ggplot(goodDists_lt10000, aes(x = orientation, y = distance)) + geom_boxplot() + facet_grid(. ~ genome) +
  theme_bw() + 
  xlab("Element type") + 
  ylab("Intergenic length (bp)") + 
  ggtitle("Intergenic sequence length by orientation") + 
  theme(axis.text = element_text(size = 14, color  = "black"), axis.title = element_text(size = 14), plot.title = element_text(size = 16))
ggsave("~/Dropbox/Bladderwort/0_Plots/intergenic_length_comparison_byOrientation.png")

#plotting expression summaries

allExpr <- read.table("~/Xfer/jwlab/Bladderwort_3pri/3_Quantification/genePairs.3primeFoldChange.tsv", header=T, sep="\t")
all_expr_distGt0 <- subset(allExpr, distance > 0)

ggplot(all_expr_distGt0, aes(x = orientation, y = log(all_tissues_foldChange))) + 
  geom_boxplot() +
  xlab("Element type") + 
  ylab("Avg fold change (log transformed)") + 
  ggtitle("Orientation by log(foldChange) - my genome") + 
  theme(axis.text = element_text(size = 14, color  = "black"), axis.title = element_text(size = 14), plot.title = element_text(size = 16))
ggsave("~/Dropbox/Bladderwort/0_Plots/orientation_by_fold_change.png")

# plotting corr summaries
ggplot(all_expr_distGt0, aes(x = orientation, y = exprCorr)) + 
    geom_violin() + 
    geom_jitter(position=position_jitter(width=.1, height = 0)) +
  xlab("Element type") + 
  ylab("Expression correlation") + 
  ggtitle("Orientation by expression correlation - my genome") + 
  theme(axis.text = element_text(size = 14, color  = "black"), axis.title = element_text(size = 14), plot.title = element_text(size = 16))

ggsave("~/Dropbox/Bladderwort/0_Plots/orientation_by_correlation.png")
