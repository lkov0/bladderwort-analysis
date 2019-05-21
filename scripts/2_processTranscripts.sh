#!/bin/bash

#2_processTranscripts.sh



# fix gff file to include genic regions

# use cufflinks to find genic regions showing expression
for bam in $(ls $alignmentDir/*.bam | sed "s/.gsnap.bam//g" | sed "s|$alignmentDir/||g"); do
    for genome in $(ls $dataDir/New_Genome/*.genic.gff); do
        stem=${genome/$dataDir\/*_Genome\/}
        stem=${stem/.genic.gff/}
        echo cufflinks --GTF $genome -o $cuffDir/$bam $alignmentDir/$bam.gsnap.bam
    done
done


#get expression values in FPKM for each genic locus
cufflinks --GTF /home/lynseykovar/Work/Bladderwort/0_Data/New_Genome/u.gibba_NEW.genic.gff -o $cuffDir/SRR094438_U.gibba_NEW $alignmentDir/SRR094438_U.gibba_NEW.gsnap.bam
cufflinks --GTF /home/lynseykovar/Work/Bladderwort/0_Data/New_Genome/u.gibba_NEW.genic.gff -o $cuffDir/SRR768657_U.gibba_NEW $alignmentDir/SRR768657_U.gibba_NEW.gsnap.bam
cufflinks --GTF /home/lynseykovar/Work/Bladderwort/0_Data/Old_Genome/u.gibba_OLD.genic.gff -o $cuffDir/SRR094438_U.gibba_OLD $alignmentDir/SRR094438_U.gibba_OLD.gsnap.bam
cufflinks --GTF /home/lynseykovar/Work/Bladderwort/0_Data/Old_Genome/u.gibba_OLD.genic.gff -o $cuffDir/SRR768657_U.gibba_OLD $alignmentDir/SRR768657_U.gibba_OLD.gsnap.bam
