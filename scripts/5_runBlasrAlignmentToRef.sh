#!/bin/bash

genome=$1
outBam=$2
pbioDir=$3
dataDir=$4
scriptDir=$5
alignDir=$6
procs=$7
# 
# module load BLASR
# module load SAMtools
# module load BEDTools
# module load BLAST
# module load R

prefix=merged_reads.u.gibba_NEW.alignment
reads=$pbioDir/merged_reads.bam

# samtools merge $reads $pbioDir/*.bam

# blasr $reads $genome -nproc $procs -unaligned $prefix.unaligned.fasta --bam --out $alignDir/$prefix.bam
# samtools sort $alignDir/$prefix.bam -o $alignDir/$prefix.sorted.bam
# samtools index $alignDir/$prefix.sorted.bam

# generate genome file (format chrName\tchrSize)
# cat $genome | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > $genome.contigSizes

# run bedtools to get coverage info
# bedtools genomecov -bga -ibam $alignDir/$prefix.sorted.bam -g $genome.contigSizes > $alignDir/$prefix.genomecov.txt
# grep -w 0$ $alignDir/$prefix.genomecov.txt > $alignDir/$prefix.nocov.regions.txt

# get read lengths and fasta files for each genome

# for stem in $(ls $pbioDir/*.bam | sed "s/.subreads.bam//g"); do
#     samtools view -c $stem.subreads.bam
#     bam2fastq -o $stem $stem.subreads.bam
#     bam2fasta -o $stem $stem.subreads.bam
#     gunzip $stem.fastq.gz
#     gunzip $stem.fasta.gz
#     awk '{if(NR%4==2) print length($1)}' $stem.fastq > $stem.readLengths.txt
# done

#Blasting unaligned reads to old and new genome
# makeblastdb -dbtype nucl -in $dataDir/Old_Genome/Utricularia_gibba.4.1.fa -out u.gibba_OLD
# blastn -db $dataDir/New_Genome/u.gibba_NEW -query $prefix.unaligned.fasta -num_threads 8 -perc_identity 95 -num_alignments 5 -out $stem.unaligned.u.gibba_NEW.blastOut.txt -outfmt "6 std qlen"
# blastn -db $dataDir/New_Genome/u.gibba_OLD -query $prefix.unaligned.fasta -num_threads 8 -perc_identity 95 -num_alignments 5 -out $stem.unaligned.u.gibba_OLD.blastOut.txt -outfmt "6 std qlen"

# #blasting raw reads to old and new genome, obtaining DB hits as bed files, then checking overlaps with genes.
# for stem in $(ls $pbioDir/m54193*.bam | sed "s/.subreads.bam//g"); do
# blastn -db $dataDir/New_Genome/u.gibba_NEW -query $stem.fasta -num_threads 8 -perc_identity 95 -num_alignments 5 -out $stem.u.gibba_NEW.blastOut.txt -outfmt "6 std qlen"
# blastn -db $dataDir/Old_Genome/u.gibba_OLD -query $stem.fasta -num_threads 8 -perc_identity 95 -num_alignments 5 -out $stem.u.gibba_OLD.blastOut.txt -outfmt "6 std qlen"
# Rscript $scriptDir/blastToBed.R -i $stem.u.gibba_NEW.blastOut.txt -o $stem.u.gibba_NEW.blastOut.dbhits.bed
# awk -F"\t" '{print $1, $2, $3}' $stem.u.gibba_NEW.blastOut.dbhits.bed | sed "s/ /\t/g" | sed "s/lcl|//g" > $stem.u.gibba_NEW.blastOut.dbhits.bed
# bedtools intersect -u -a $dataDir/New_Genome/u.gibba_NEW.genes.simple.bed -b $stem.u.gibba_NEW.blastOut.dbhits.bed > $stem.genesWithHits_blast.txt
# done

# #checking bedgraph file for overlap with genes
# grep -v 0$ $alignDir/$prefix.genomecov.txt | awk -F"\t" '{print $1, $2, $3}' | sed "s/ /\t/g" | sed "s/lcl|//g" > $alignDir/$prefix.cov.regions.bed
# #writes entries in -a bam that have overlap in -b bam
# awk '{print $1, $2, $3}' $dataDir/New_Genome/u.gibba_NEW.genes.bed | sed "s/ /\t/g" > $dataDir/New_Genome/u.gibba_NEW.genes.simple.bed
# bedtools intersect -u -a $dataDir/New_Genome/u.gibba_NEW.genes.simple.bed -b $alignDir/$prefix.cov.regions.bed > $alignDir/genesWithHits_blasr.txt

