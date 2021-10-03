# gcContentAnalyis.R
# investigate gc content of candidate regions and ouput plots of results.

library(ggplot2)

ins_gc <- read.table("~/Work/Bladderwort/5_Meme/ins_cand.gccontent.txt", header = T)
term_gc <- read.table("~/Work/Bladderwort/5_Meme/term_cand.gccontent.txt", header = T)
uni_gc <- read.table("~/Work/Bladderwort/5_Meme/uni_cand.gccontent.txt", header = T)
bid_gc <- read.table("~/Work/Bladderwort/5_Meme/bid_cand.gccontent.txt", header = T)
all_gc <- read.table("~/Work/Bladderwort/5_Meme/all_pairs.gccontent.txt", header = T)

ins_gc$type <- "insulator"
term_gc$type <- "terminator"
uni_gc$type <- "unidirectional promoter"
bid_gc$type <- "bidirectional promoter"
all_gc$type <- "all"

gcsummary <- rbind(ins_gc, term_gc, uni_gc, bid_gc)
gcsummary$type <- as.factor(gcsummary$type)

ggplot(gcsummary, aes(y = X.GC, x = type)) + geom_dotplot(binaxis = 'y', stackdir = 'center') + theme_bw() + ylab("GC Content") + xlab("Type") + ggtitle("GC content by putative CRE class")
ggsave("~/Dropbox/0_Plots/GCcontent_CreCandidates.png")
fit <- aov(X.GC ~ type, data = gcsummary)
summary(fit)

