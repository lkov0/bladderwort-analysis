# runCemiTools.R
# run CemiTools cluster analysis on the 3' expression dataset, get promoter regions associated with each module

library(CEMiTool)

# import raw htSeq count data (no L3 since that's an outlier)

B1 <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/1B_htSeq_counts.ALL.txt", row.names = 1)
L1 <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/1L_htSeq_counts.ALL.txt", row.names = 1)
S1 <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/1S_htSeq_counts.ALL.txt", row.names = 1)
R1 <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/1R_htSeq_counts.ALL.txt", row.names = 1)
B2 <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/2B_htSeq_counts.ALL.txt", row.names = 1)
L2 <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/2L_htSeq_counts.ALL.txt", row.names = 1)
S2 <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/2S_htSeq_counts.ALL.txt", row.names = 1)
R2 <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/2R_htSeq_counts.ALL.txt", row.names = 1)
B3 <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/3B_htSeq_counts.ALL.txt", row.names = 1)
S3 <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/3S_htSeq_counts.ALL.txt", row.names = 1)
R3 <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/3R_htSeq_counts.ALL.txt", row.names = 1)

names(B1) <- c("B1")
names(L1) <- c("L1")
names(S1) <- c("S1")
names(R1) <- c("R1")
names(B2) <- c("B2")
names(L2) <- c("L2")
names(S2) <- c("S2")
names(R2) <- c("R2")
names(B3) <- c("B3")
names(S3) <- c("S3")
names(R3) <- c("R3")

all_counts_pbio <- cbind(B1, L1, S1, R1, B2, L2, S2, R2, B3, S3, R3)
all_counts_pbio <- all_counts_pbio[1:19203,]

cem <- cemitool(all_counts_pbio,  apply_vst = TRUE)
generate_report(cem = cem, force = T)
write_files(cem = cem, force = T)
save_plots(cem = cem, force = T)

# make bedfiles for upstream regions of module genes
module_assignments <- module_genes(cem)

# there are 10 modules
M1 <- subset(module_assignments, modules == "M1") #611
M2 <- subset(module_assignments, modules == "M2") #276
M3 <- subset(module_assignments, modules == "M3") #124
M4 <- subset(module_assignments, modules == "M4") #108
M5 <- subset(module_assignments, modules == "M5") #105
M6 <- subset(module_assignments, modules == "M6") #83
M7 <- subset(module_assignments, modules == "M7") #53
M8 <- subset(module_assignments, modules == "M8") #47
M9 <- subset(module_assignments, modules == "M9") #43
M10 <- subset(module_assignments, modules == "M10") #34

#obtain annotation info
annos <- read.table("~/Work-xfer/genomes/utricularia/scaffolds.ugibba_lk.genic.gff")
annos$genes <- gsub("ID=", "", annos$V9)
annos$genes <- gsub(";.*", "", annos$genes)

M1 <- merge(annos, M1, by = "genes")
M2 <- merge(annos, M2, by = "genes")
M3 <- merge(annos, M3, by = "genes")
M4 <- merge(annos, M4, by = "genes")
M5 <- merge(annos, M5, by = "genes")
M6 <- merge(annos, M6, by = "genes")
M7 <- merge(annos, M7, by = "genes")
M8 <- merge(annos, M8, by = "genes")
M9 <- merge(annos, M9, by = "genes")
M10 <- merge(annos, M10, by = "genes")

# get promoter regions (-200 from TSS, convert from 1 based indexing to 0 based indexing)
getPromoterRegions <- function(gff) {
  gff$promStart = 0
  gff$promStop = 0
  
  for(i in 1:nrow(gff)) {
    if(gff[i,"V7"] == "+")
      if(gff[i,"V4"] < 201) {
        gff[i, "promStart"] = 0
        gff[i, "promStop"] = gff[i, "promStart"]
      }
      else {
        gff[i, "promStart"] = gff[i, "V4"] - 201
        gff[i, "promStop"] = gff[i, "V4"]
      }
    else {
      gff[i, "promStart"] = gff[i, "V5"] - 1
      gff[i, "promStop"] = gff[i, "V5"] + 200
    }
  }
  
  bed <- data.frame(chrom = gff$V1, start = gff$promStart, stop = gff$promStop, name = gff$genes, score = ".", strand = gff$V7)
  
  return(bed)
}

write.table(getPromoterRegions(M1), file = "~/Xfer/Bladderwort_3pri/5_Meme/CEMiTools/module1.cemitool.bed", row.names = F, col.names = F, sep = "\t", quote = F)
write.table(getPromoterRegions(M2), file = "~/Xfer/Bladderwort_3pri/5_Meme/CEMiTools/module2.cemitool.bed", row.names = F, col.names = F, sep = "\t", quote = F)
write.table(getPromoterRegions(M3), file = "~/Xfer/Bladderwort_3pri/5_Meme/CEMiTools/module3.cemitool.bed", row.names = F, col.names = F, sep = "\t", quote = F)
write.table(getPromoterRegions(M4), file = "~/Xfer/Bladderwort_3pri/5_Meme/CEMiTools/module4.cemitool.bed", row.names = F, col.names = F, sep = "\t", quote = F)
write.table(getPromoterRegions(M5), file = "~/Xfer/Bladderwort_3pri/5_Meme/CEMiTools/module5.cemitool.bed", row.names = F, col.names = F, sep = "\t", quote = F)
write.table(getPromoterRegions(M6), file = "~/Xfer/Bladderwort_3pri/5_Meme/CEMiTools/module6.cemitool.bed", row.names = F, col.names = F, sep = "\t", quote = F)
write.table(getPromoterRegions(M7), file = "~/Xfer/Bladderwort_3pri/5_Meme/CEMiTools/module7.cemitool.bed", row.names = F, col.names = F, sep = "\t", quote = F)
write.table(getPromoterRegions(M8), file = "~/Xfer/Bladderwort_3pri/5_Meme/CEMiTools/module8.cemitool.bed", row.names = F, col.names = F, sep = "\t", quote = F)
write.table(getPromoterRegions(M9), file = "~/Xfer/Bladderwort_3pri/5_Meme/CEMiTools/module9.cemitool.bed", row.names = F, col.names = F, sep = "\t", quote = F)
write.table(getPromoterRegions(M10), file = "~/Xfer/Bladderwort_3pri/5_Meme/CEMiTools/module10.cemitool.bed", row.names = F, col.names = F, sep = "\t", quote = F)
