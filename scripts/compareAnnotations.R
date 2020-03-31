#compareAnnotations.r

#get count matrix for new genome + new annotation
S1_new <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/1S_htSeq_counts.txt")
L1_new <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/1L_htSeq_counts.txt")
B1_new <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/1B_htSeq_counts.txt")
R1_new <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/1R_htSeq_counts.txt")
S2_new <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/2S_htSeq_counts.txt")
L2_new <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/2L_htSeq_counts.txt")
B2_new <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/2B_htSeq_counts.txt")
R2_new <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/2R_htSeq_counts.txt")
S3_new <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/3S_htSeq_counts.txt")
L3_new <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/3L_htSeq_counts.txt")
B3_new <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/3B_htSeq_counts.txt")
R3_new <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/3R_htSeq_counts.txt")

names(S1_new) <- c("geneId", "stem1")
names(L1_new) <- c("geneId", "leaf1")
names(B1_new) <- c("geneId", "bladder1")
names(R1_new) <- c("geneId", "rhizoid1")
names(S2_new) <- c("geneId", "stem2")
names(L2_new) <- c("geneId", "leaf2")
names(B2_new) <- c("geneId", "bladder2")
names(R2_new) <- c("geneId", "rhizoid2")
names(S3_new) <- c("geneId", "stem3")
names(L3_new) <- c("geneId", "leaf3")
names(B3_new) <- c("geneId", "bladder3")
names(R3_new) <- c("geneId", "rhizoid3")
htseq_table_new <- data.frame(geneId = S1_new[,1], stem1 = S1_new[,2], leaf1 = L1_new[,2], bladder1 = B1_new[,2], rhizoid1 = R1_new[,2], stem2 = S2_new[,2], leaf2 = L2_new[,2], bladder2 = B2_new[,2], rhizoid2 = R2_new[,2], stem3 = S3_new[,2], leaf3 = L3_new[,2], bladder3 = B3_new[,2], rhizoid3 = R3_new[,2])

#get count matrix for old genome + old annotation
S1_old <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/1S_htSeq_counts_pbioAnnos.gt5k.txt")
L1_old <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/1L_htSeq_counts_pbioAnnos.gt5k.txt")
B1_old <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/1B_htSeq_counts_pbioAnnos.gt5k.txt")
R1_old <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/1R_htSeq_counts_pbioAnnos.gt5k.txt")
S2_old <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/2S_htSeq_counts_pbioAnnos.gt5k.txt")
L2_old <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/2L_htSeq_counts_pbioAnnos.gt5k.txt")
B2_old <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/2B_htSeq_counts_pbioAnnos.gt5k.txt")
R2_old <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/2R_htSeq_counts_pbioAnnos.gt5k.txt")
S3_old <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/3S_htSeq_counts_pbioAnnos.gt5k.txt")
L3_old <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/3L_htSeq_counts_pbioAnnos.gt5k.txt")
B3_old <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/3B_htSeq_counts_pbioAnnos.gt5k.txt")
R3_old <- read.table("~/Xfer/Bladderwort_3pri/2_Alignment/htSeq_out/3R_htSeq_counts_pbioAnnos.gt5k.txt")

names(S1_old) <- c("geneId", "stem1")
names(L1_old) <- c("geneId", "leaf1")
names(B1_old) <- c("geneId", "bladder1")
names(R1_old) <- c("geneId", "rhizoid1")
names(S2_old) <- c("geneId", "stem2")
names(L2_old) <- c("geneId", "leaf2")
names(B2_old) <- c("geneId", "bladder2")
names(R2_old) <- c("geneId", "rhizoid2")
names(S3_old) <- c("geneId", "stem3")
names(L3_old) <- c("geneId", "leaf3")
names(B3_old) <- c("geneId", "bladder3")
names(R3_old) <- c("geneId", "rhizoid3")
htseq_table_old <- data.frame(geneId = S1_old[,1], stem1 = S1_old[,2], leaf1 = L1_old[,2], bladder1 = B1_old[,2], rhizoid1 = R1_old[,2], stem2 = S2_old[,2], leaf2 = L2_old[,2], bladder2 = B2_old[,2], rhizoid2 = R2_old[,2], stem3 = S3_old[,2], leaf3 = L3_old[,2], bladder3 = B3_old[,2], rhizoid3 = R3_old[,2])

#how many genes have hits?

#what proportion of reads were classified?

