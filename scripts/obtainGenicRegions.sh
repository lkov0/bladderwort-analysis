#!/bin/bash

# Author: Lynsey Kovar
# Date: Jan 24 2019
# Purpose: script for obtaining genic regions from gff file containing only CDS regions

dataDir=$1
scriptDir=$2
gffFile=$3


stem=${gffFile/.gff/}

grep "##" $gffFile > $stem.commentedLines.txt

# Rscript to convert from gff to bed - this worked because the number of genes produced by the script (29666) matched up with the number of genes listed in the coge annotation
Rscript $scriptDir/cdsToGenicRegions.R -i $stem.gff -o $stem.genic.gff.txt

# add commented lines back to produce authentic gff file
cat $stem.commentedLines.txt $stem.genic.gff.txt > $stem.genic.gff