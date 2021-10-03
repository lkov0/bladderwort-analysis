# getDiffExpGenePairs_3prime.R
# use gene pair table to extract top candidates for each regulatory element class. 

library(matrixStats)

genePairs <- read.table("~/Xfer/Bladderwort_3pri/3_Quantification/genePairs.3primeFoldChange.tsv", header=T, sep="\t")
#add gene1 and gene2 average expression values

#distance between gene pairs - negative indicates overlapping genes
genePairs$distance <- genePairs$start2 - genePairs$stop1

genePairs$gene1_avg <- rowMeans(genePairs[,c(13:23)])
genePairs$gene2_avg <- rowMeans(genePairs[,c(24:34)])

#get median for gene1 and gene2
gene1_median <- median(subset(genePairs, gene1_avg > 0)$gene1_avg) # 9.83
gene2_median <- median(subset(genePairs, gene2_avg > 0)$gene2_avg) # 9.61

uni_cand <- read.table("~/Xfer/Bladderwort_3pri/4_CandidateMining/uni_cand.3pri.tsv", header=T, sep="\t", stringsAsFactors = F)
bid_cand <- read.table("~/Xfer/Bladderwort_3pri/4_CandidateMining/bid_cand.3pri.tsv", header=T, sep="\t", stringsAsFactors = F)
ter_cand <- read.table("~/Xfer/Bladderwort_3pri/4_CandidateMining/term_cand.3pri.tsv", header=T, sep="\t", stringsAsFactors = F)
ins_cand <- read.table("~/Xfer/Bladderwort_3pri/4_CandidateMining/ins_cand.3pri.tsv", header=T, sep="\t", stringsAsFactors = F)

#####
#unidirectional promoter candidates
#####
# parallel orientation, fold change > 90% quantile, right gene > median expression value, shared between datasets, expression independence
uni_cand <- subset(genePairs, orientation == "parallel") #10574 genes are parallel orientation
quantile(genePairs$all_tissues_foldChange, seq(0,1,0.1), na.rm=TRUE) #90% quantile = 15
uni_cand <- subset(uni_cand, (strand1 == '+' & all_tissues_foldChange >= 15) | (strand1 == '-' & all_tissues_foldChange <= -15)) #157 of these are expressed > 90% quantile
uni_cand.fwd <- subset(uni_cand, strand1 == "+")
uni_cand.rev <- subset(uni_cand, strand2 == "-")
uni_cand.fwd <- subset(uni_cand.fwd, gene2_avg > gene2_median)
uni_cand.rev <- subset(uni_cand.rev, gene1_avg > gene1_median)
uni_cand <- rbind(uni_cand.fwd, uni_cand.rev)
uni_cand <- subset(uni_cand, exprCorr < 0.5 & exprCorr > -0.5) #43
uni_cand$med_all_tissues_foldChange <- abs(rowMedians(as.matrix(uni_cand[,c(45:48)])))
uni_cand <- uni_cand[order(-uni_cand$med_all_tissues_foldChange, uni_cand$distance), ] # order
uni_cand <- subset(uni_cand, distance <=1000 & distance > 0) # 43 candidates
rownames(uni_cand) <- seq(1:nrow(uni_cand))
uni_cand$candidate_name <- paste("uni", rownames(uni_cand), sep="_")
write.table(uni_cand, "~/Xfer/Bladderwort_3pri/4_CandidateMining/uni_cand.3pri.tsv", row.names=F, col.names=T, sep="\t", quote=F)

uni_cand_easyRead <- uni_cand[,c(55,2,1,3,4,5,6,7,8,9,10,11,12,46,47,48,49,50,51,52,53,54)]
write.table(uni_cand_easyRead, "~/Xfer/Bladderwort_3pri/4_CandidateMining/uni_cand.3pri.summary.tsv", row.names=F, col.names=T, sep="\t", quote=F)
write.table(uni_cand_easyRead[1:20,], "~/Xfer/Bladderwort_3pri/4_CandidateMining/uni_cand.3pri.top20.tsv", row.names=F, col.names=T, sep="\t", quote=F)

#####
#bidirectional promoter candidates
#####
# divergent orientation
# fc < median fc in expression
# > median expression value
# shared between datasets
# highly correlated in expression
bid_cand <- subset(genePairs, orientation == "divergent") # 1070 divergent gene pairs
median_fc <- median(subset(genePairs, all_tissues_foldChange > 0)$all_tissues_foldChange)
bid_cand <- subset(bid_cand, all_tissues_foldChange < median_fc) # 705 less than median fold change value
bid_cand <- subset(bid_cand, gene1_avg > gene1_median) # 676 showing expression greater than median expression value
bid_cand <- subset(bid_cand, exprCorr > 0.5 | exprCorr < -0.5) # 274 showing weakly correlated expression
bid_cand$med_all_tissues_foldChange <- abs(rowMedians(as.matrix(bid_cand[,c(45:48)])))
bid_cand <- bid_cand[order(bid_cand$med_all_tissues_foldChange, -bid_cand$exprCorr, bid_cand$distance), ] # order
bid_cand <- subset(bid_cand, distance <=1000 & distance > 0) # 161 candidates
rownames(bid_cand) <- seq(1:nrow(bid_cand))
bid_cand$candidate_name <- paste("bid", rownames(bid_cand), sep="_")
write.table(bid_cand, "~/Xfer/Bladderwort_3pri/4_CandidateMining/bid_cand.3pri.tsv", row.names=F, col.names=T, sep="\t", quote=F)

bid_cand_easyRead <- bid_cand[,c(55,2,1,3,4,5,6,7,8,9,10,11,12,46,47,48,49,50,51,52,53,54)]
write.table(bid_cand_easyRead, "~/Xfer/Bladderwort_3pri/4_CandidateMining/bid_cand.3pri.summary.tsv", row.names=F, col.names=T, sep="\t", quote=F)
write.table(bid_cand_easyRead[1:20,], "~/Xfer/Bladderwort_3pri/4_CandidateMining/bid_cand.3pri.top20.tsv", row.names=F, col.names=T, sep="\t", quote=F)

#terminator candidates
# parallel and convergent
# fc > 90% quantile
# left gene > median expression value
# shared between datasets
term_cand <- subset(genePairs, orientation == "parallel" | orientation == "convergent") #15411 parallel or convergent
term_cand <- subset(term_cand, (strand1 == '-' & all_tissues_foldChange >= 15) | (strand1 == '+' & all_tissues_foldChange <= -15)) #237
term_cand <- subset(term_cand, (strand1 == '-' & gene2_avg > gene2_median) | (strand1 == "+" & gene1_avg > gene1_median)) #237
term_cand <- subset(term_cand, exprCorr < 0.5 & exprCorr > -0.5) #93
term_cand$med_all_tissues_foldChange <- abs(rowMedians(as.matrix(term_cand[,c(45:48)])))
term_cand <- term_cand[order(-term_cand$med_all_tissues_foldChange, term_cand$distance), ] # order
term_cand <- subset(term_cand, distance <=1000 & distance > 0) # 75 candidates
rownames(term_cand) <- seq(1:nrow(term_cand))
term_cand$candidate_name <- paste("ter", rownames(term_cand), sep="_")
write.table(term_cand, "~/Xfer/Bladderwort_3pri/4_CandidateMining/term_cand.3pri.tsv", row.names=F, col.names=T, sep="\t", quote=F)

term_cand_easyRead <- term_cand[,c(55,2,1,3,4,5,6,7,8,9,10,11,12,46,47,48,49,50,51,52,53,54)]
write.table(term_cand_easyRead, "~/Xfer/Bladderwort_3pri/4_CandidateMining/term_cand.3pri.summary.tsv", row.names=F, col.names=T, sep="\t", quote=F)
write.table(term_cand_easyRead[1:20,], "~/Xfer/Bladderwort_3pri/4_CandidateMining/term_cand.3pri.top20.tsv", row.names=F, col.names=T, sep="\t", quote=F)

#insulator candidates
# divergent and convergent orientation
# fold change > 90% quantile
# > median expression value
# shared between datasets
ins_cand <- subset(genePairs, orientation == "divergent" | orientation == "convergent") #9382 parallel or convergent
ins_cand <- subset(ins_cand, (strand1 == '-' & all_tissues_foldChange >= 15) | (strand1 == '+' & all_tissues_foldChange <= -15)) #164
ins_cand <- subset(ins_cand, (strand1 == "-" & gene2_avg > gene2_median) | (strand1 == "+" & gene1_avg > gene1_median)) #164
ins_cand <- subset(ins_cand, exprCorr < 0.5 & exprCorr > -0.5) #63
ins_cand <- subset(ins_cand, distance <=1000 & distance > 0) #43
ins_cand$med_all_tissues_foldChange <- abs(rowMedians(as.matrix(ins_cand[,c(45:48)])))
ins_cand <- ins_cand[order(-ins_cand$med_all_tissues_foldChange, ins_cand$distance), ] # order
rownames(ins_cand) <- seq(1:nrow(ins_cand))
ins_cand$candidate_name <- paste("ins", rownames(ins_cand), sep="_")
write.table(ins_cand, "~/Xfer/Bladderwort_3pri/4_CandidateMining/ins_cand.3pri.tsv", row.names=F, col.names=T, sep="\t", quote=F)

ins_cand_easyRead <- ins_cand[,c(55,2,1,3,4,5,6,7,8,9,10,11,12,46,47,48,49,50,51,52,53,54)]
write.table(ins_cand_easyRead, "~/Xfer/Bladderwort_3pri/4_CandidateMining/ins_cand.3pri.summary.tsv", row.names=F, col.names=T, sep="\t", quote=F)
write.table(ins_cand_easyRead[1:20,], "~/Xfer/Bladderwort_3pri/4_CandidateMining/ins_cand.3pri.top20.tsv", row.names=F, col.names=T, sep="\t", quote=F)

