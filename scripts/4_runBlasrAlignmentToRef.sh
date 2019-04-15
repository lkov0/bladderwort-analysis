#!/bin/bash

reads=$1
genome=$2
outBam=$3
procs=$4

# module load BLASR
# module load SamTools
module load BEDTools

prefix=${outBam/.bam/}
# echo $prefix

# blasr $reads $genome --nproc $procs --bam --out $outBam
# samtools sort $outBam -o $prefix.sorted.bam
# samtools index $prefix.sorted.bam

# # generate genome file (format chrName\tchrSize)
# cat $genome | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > $genome.contigSizes
# 
# # run bedtools to get coverage info
# bedtools genomecov -bga -ibam $prefix.sorted.bam -g $genome.contigSizes > $prefix.genomecov.txt
grep -w 0$ $prefix.genomecov.txt > $prefix.nocov.regions.txt
