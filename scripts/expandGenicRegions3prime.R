#expandGenicRegions3prime.R

gff <- read.table("Xfer-work/bladderwort/genomes/utricularia/u.gibba_NEW.forStar.gff", header=F, sep="\t")

gff.pos <- subset(gff, V7 == "+")
gff.neg <- subset(gff, V7 == "-")

rownames(gff.pos) = seq(1, nrow(gff.pos))
rownames(gff.neg) = seq(1, nrow(gff.neg))

gff.pos <- gff.pos[order(gff.pos$V1, gff.pos$V4),]
gff.neg <- gff.neg[order(gff.neg$V1, gff.neg$V4),]

gff.pos.new <- gff.pos
gff.neg.old <- gff.neg

for(i in 1:nrow(gff.pos)) {
    if(i == nrow(gff.pos)) {
        gff.pos.new[i, "V5"] =  gff.pos[i, "V5"] + 500
    }
    else {
        if(gff.pos[i, "V1"] == gff.pos[i + 1, "V1"]) {
            if(gff.pos[i + 1, "V4"] - gff.pos[i, "V5"] <= 500) {
                gff.pos.new[i, "V5"] =  (gff.pos[i + 1, "V4"] - 1)
            }
            if(gff.pos[i + 1, "V4"] - gff.pos[i, "V5"] > 500) {
                gff.pos.new[i, "V5"] =  gff.pos[i, "V5"] + 500
            }
        }
        else {
            gff.pos.new[i, "V5"] =  gff.pos[i, "V5"] + 500
        }
    }    
}

for(i in 1:nrow(gff.neg)) {
    if(i == 1) {
        if(gff.neg.new[i, "V4"] - 500 >= 0) {
            gff.neg.new[i, "V4"] =  gff.neg[i, "V4"] - 500
        }
        else {
            gff.neg.new[i, "V4"] = 1
        }
    }
    else {
        if(gff.neg[i, "V1"] == gff.neg[i - 1, "V1"]) {
            if((gff.neg[i, "V4"] - gff.neg[i - 1, "V5"] <= 500) & (gff.neg[i, "V4"] - gff.neg[i - 1, "V5"] > 0)) {
                gff.neg.new[i, "V4"] =  (gff.neg[i - 1, "V5"] + 1)
            }
            if(gff.neg[i, "V4"] - gff.neg[i - 1, "V5"] > 500) {
                gff.neg.new[i, "V4"] =  gff.neg[i, "V4"] - 500
            }
        }
        else {
            if(gff.neg.new[i, "V4"] - 500 >= 0) {
                gff.neg.new[i, "V4"] =  gff.neg[i, "V4"] - 500
        }
            else {
                gff.neg.new[i, "V4"] = 1
            }
        }
    }
}

gff.new <- rbind(gff.pos.new, gff.neg.new)
gff.new <- format(gff.new, scientific=F)
write.table(gff.new, file="/work/jawlab/data/bladderwort/genomes/utricularia/u.gibba_NEW.genes_for_3prime.txt", row.names=F, col.names=F, quote=F, sep="\t")
