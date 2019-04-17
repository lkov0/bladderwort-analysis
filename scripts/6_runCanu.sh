#!/bin/bash
#6_runErrorCorrectWithCanu.sh

canuDir=$1
readDir=$2
nProcs=$3

module load canu

#combine reads

cat $readDir/*.fastq > merged.PacBio.fastq
canu -correct -p bladderwort_test -d $canuDir/bladderwort_test genomeSize=100m -pacbio-raw merged.PacBio.fastq