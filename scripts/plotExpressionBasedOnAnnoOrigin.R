# plot expression based on old or new genome. 

library(ggplot2)

expr_1S <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/1S_htSeq_counts.ALL.includingPacBio.txt", sep="\t")
pbioAnnos <- read.table("~/Xfer-work/genomes/utricularia/annos_missing_in_my_genome.bed")
allAnnos_gff <- read.table("~/Xfer-work/genomes/utricularia/scaffolds.ugibba_lk.ALL.includingPacBio.gff")

allAnnos_gff$annoOrigin <- "na"

for(i in 1:nrow(allAnnos_gff)) {
  if(allAnnos_gff[i,"bedKey"] %in% allAnnos$bedKey) {
    allAnnos_gff[i, "annoOrigin"] <- "PacBio Genome"
  }
  else {
    allAnnos_gff[i, "annoOrigin"] <- "Our Genome"
  }
}

pbioAnnos_gff <- merge(pbioAnnos_gff, expr_1S, by.x="V9", by.y="V1")

ggplot(pbioAnnos_gff, aes(x=annoOrigin, y=logExpression)) + geom_violin(col = "black") + 
  geom_boxplot(col = "black", outlier.size = 3) + 
  ylab("log(read counts per gene)") + 
  xlab("Annotation Origin") + 
  ggtitle("Combined annotation gene counts summary") + 
  theme_bw() + theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black"), title = element_text(size = 14))

ggsave("~/Dropbox/Figures/annotationCountsPerGene.png", dpi = 400)
