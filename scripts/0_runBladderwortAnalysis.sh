###########################
# Author: Lynsey Kovar
# Date: Jan 16 2019
# Purpose: Parent shell script for bladderwort analysis.
# Software needed: gmap (v2018-07-04), gsnap (v2018-07-04), samtools (v1.6 using htslib 1.6), cufflinks (v2.2.1), ncbi blast (v2.7.1+), meme suite (5.0.4), mummer (v4.0.0.beta2)
# for interactive session: qsub -I -q s_interq -l walltime=10:00:00 -l nodes=1:ppn=16:inter -l mem=16gb
###########################

# parentDir=~/Xfer/Bladderwort
# scriptDir=~/Xfer/Repositories/bladderwort-analysis/scripts
parentDir=/scratch/lk82153/jwlab/Bladderwort_3pri
scriptDir=~/Repositories/bladderwort-analysis/scripts
dataDir=/work/jawlab/data/bladderwort
threeprimeDir=$dataDir/nextseq_rna_3prime
assemblyDir=$parentDir/1_Assembly
alignmentDir=$parentDir/2_Alignment
quantDir=$parentDir/3_Quantification
analysisDir=$parentDir/4_CandidateMining
memeDir=$parentDir/5_Meme
consDir=$parentDir/6_ConservationAnalysis

if [ ! -e $parentDir ]; then mkdir $parentDir; fi
if [ ! -e $assemblyDir ]; then mkdir $assemblyDir; fi
if [ ! -e $alignmentDir ]; then mkdir $alignmentDir; fi
if [ ! -e $quantDir ]; then mkdir $quantDir; fi
if [ ! -e $analysisDir ]; then mkdir $analysisDir; fi
if [ ! -e $memeDir ]; then mkdir $memeDir; fi
if [ ! -e $consDir ]; then mkdir $consDir; fi

refGenome=$dataDir/genomes/utricularia/scaffolds.ugibba_lk.fasta

PROCS=24

# download RNA-Seq reads from NCBI.
# SRR094438: low coverage, average of 31,500 reads for each condition (~820k total)
# SRR768657: higher coverage, but conditions were pooled leading to dimming of tissue-specific expression signal. 

# PATH="/home/lynseykovar/Programs/sratoolkit.2.9.2-centos_linux64/bin/:$PATH"
 
# download SRA data
# for file in SRR094438 SRR768657; do 
# fastq-dump -O $dataDir --gzip $fastqDir/$file
# done
# 
# $scriptDir/2_runAlignment.sh $threeprimeDir $alignmentDir $refGenome $dataDir $PROCS

$scriptDir/3_evaluateGenePairs.sh $parentDir $scriptDir $dataDir $alignmentDir $quantDir

# $scriptDir/4_findMotifs.sh $parentDir $scriptDir $dataDir $analysisDir $memeDir $NPROCS

################
# Aligning PacBio Reads to reference for our genotype
################
#running this on sapelo so the number 64 pertains to threads
# pbioAlignmentDir=$parentDir/5_PacBioAlignment
# 
# if [ ! -e $pbioAlignmentDir ] ; then mkdir $pbioAlignmentDir; fi
# # 
# # #aligning 3.1gb read file first 
# $scriptDir/5_runBlasrAlignmentToRef.sh $dataDir/New_Genome/Utricularia_gibba_v2.faa $pbioAlignmentDir/bladderwort_alignTest.bam $dataDir/PacBio $dataDir $scriptDir $pbioAlignmentDir 48
# # 
# canuDir=$parentDir/6_Canu
# 
# if [ ! -e $canuDir ] ; then mkdir $canuDir; fi

#assemble with CANU on sapelo
# $scriptDir/6_runCanu.sh $canuDir $dataDir/PacBio $dataDir $dataDir/New_Genome/Utricularia_gibba_v2.faa $scriptDir 8
