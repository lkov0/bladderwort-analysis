###########################
# Author: Lynsey Kovar
# Date: Jan 16 2019
# Purpose: Parent shell script for bladderwort analysis.
# Software needed: gmap (v2018-07-04), gsnap (v2018-07-04), samtools (v1.6 using htslib 1.6), cufflinks (v2.2.1), ncbi blast (v2.7.1+) 
###########################

parentDir=$HOME/Dropbox/Bladderwort
scriptDir=$HOME/Work/bladderwort-analysis/scripts
dataDir=$parentDir/0_Data
fastqDir=$dataDir/Fastq
alignmentDir=$parentDir/1_Alignment
cuffDir=$parentDir/2_Cufflinks

if [ ! -e $alignmentDir ]; then mkdir $alignmentDir; fi
if [ ! -e $cuffDir ]; then mkdir $cuffDir; fi

threads=8

# download RNA-Seq reads from NCBI.
# SRR094438: low coverage, average of 31,500 reads for each condition (~820k total)
# SRR768657: higher coverage, but conditions were pooled leading to dimming of tissue-specific expression signal. 

PATH="/home/lynseykovar/Programs/sratoolkit.2.9.2-centos_linux64/bin/:$PATH"
 
# download SRA data
# for file in SRR094438 SRR768657; do 
# fastq-dump -O $dataDir --gzip $fastqDir/$file
# done
# 
# $parentDir/1_runAlignment.sh $fastqDir $alignmentDir
# 
# # for genome in $dataDir/New_Genome/Utricularia_gibba_v2.faa $dataDir/Old_Genome/Utricularia_gibba.4.1.fa; do
# gmap_build -d U.gibba_NEW $dataDir/New_Genome/Utricularia_gibba_v2.faa -D $dataDir/New_Genome/
# gmap_build -d U.gibba_OLD $dataDir/Old_Genome/Utricularia_gibba.4.1.fa -D $dataDir/Old_Genome/

# fix gff file to include genic regions
# $parentDir/obtainGenicRegions.sh $dataDir $scriptDir $dataDir/New_Genome/u.gibba_NEW.gff

# for transcripts in $fastqDir/SRR094438.fastq.gz $fastqDir/SRR768657.fastq.gz; do
#     for genome in $dataDir/New_Genome/U.gibba_NEW $dataDir/Old_Genome/U.gibba_OLD; do
#         stem=${transcripts/$fastqDir\/}
#         stem=${stem/.fastq.gz/}
#         genomeDir=${genome/U.gibba_*/}
#         genomeName=${genome/\/home\/lynseykovar\/Work\/Bladderwort\/0_Data\/*_Genome\//}
#         gsnap -D $genomeDir -d $genomeName --batch=5 --novelsplicing 1 --nthreads $(($threads-2)) --ordered --format sam --output-file $alignmentDir/${stem}_${genomeName}.gsnap.sam --gunzip $transcripts;
#         samtools sort -o $alignmentDir/${stem}_${genomeName}.gsnap.bam $alignmentDir/${stem}_${genomeName}.gsnap.sam
#         samtools index $alignmentDir/${stem}_${genomeName}.gsnap.bam
#         rm $alignmentDir/${stem}_${genomeName}.gsnap.sam
#     done
# done

# use cufflinks to find genic regions showing expression
# for bam in $(ls $alignmentDir/*.bam | sed "s/.gsnap.bam//g" | sed "s|$alignmentDir/||g"); do
#     for genome in $(ls $dataDir/New_Genome/*.genic.gff); do
#         stem=${genome/$dataDir\/*_Genome\/}
#         stem=${stem/.genic.gff/}
#         echo cufflinks --GTF $genome -o $cuffDir/$bam $alignmentDir/$bam.gsnap.bam
#     done
# done


#get expression values in FPKM for each genic locus
# cufflinks --GTF /home/lynseykovar/Work/Bladderwort/0_Data/New_Genome/u.gibba_NEW.genic.gff -o $cuffDir/SRR094438_U.gibba_NEW $alignmentDir/SRR094438_U.gibba_NEW.gsnap.bam
# cufflinks --GTF /home/lynseykovar/Work/Bladderwort/0_Data/New_Genome/u.gibba_NEW.genic.gff -o $cuffDir/SRR768657_U.gibba_NEW $alignmentDir/SRR768657_U.gibba_NEW.gsnap.bam
# cufflinks --GTF /home/lynseykovar/Work/Bladderwort/0_Data/Old_Genome/u.gibba_OLD.genic.gff -o $cuffDir/SRR094438_U.gibba_OLD $alignmentDir/SRR094438_U.gibba_OLD.gsnap.bam
# cufflinks --GTF /home/lynseykovar/Work/Bladderwort/0_Data/Old_Genome/u.gibba_OLD.genic.gff -o $cuffDir/SRR768657_U.gibba_OLD $alignmentDir/SRR768657_U.gibba_OLD.gsnap.bam

#Pull out gene pairs that show high fold change in expression

analysisDir=$parentDir/3_Analysis

if [ ! -e $analysisDir ]; then mkdir $analysisDir; fi

#need to split up divergent convergent parallel before I call diffExp because since the RNA libraries are unstranded cufflinks strand calling is basically meaningless

# for cuff in $(ls $cuffDir); do
# # Rscript $scriptDir/getGenePairs.R -i $cuffDir/$cuff/genes.fpkm_tracking -o $cuffDir/$cuff/$cuff.bed
# # bedtools window -a $cuffDir/$cuff/$cuff.bed -b $cuffDir/$cuff/$cuff.bed > $analysisDir/$cuff.overlaps.txt
# # Rscript $scriptDir/getDiffExpressedNeighbors.R -i $analysisDir/$cuff.overlaps.txt -o $analysisDir/$cuff.diffExp.txt
# done

#pull out convergent, divergent, parallel gene pairs
# Rscript $scriptDir/assignConvergentDivergentParallel.R -i $dataDir/New_Genome/u.gibba_NEW.genic.gff -o $dataDir/New_Genome/u.gibba_NEW.genePairs.txt
# Rscript $scriptDir/assignConvergentDivergentParallel.R -i $dataDir/Old_Genome/u.gibba_OLD.genic.gff -o $dataDir/Old_Genome/u.gibba_OLD.genePairs.txt

#append FPKM to *.genePairs.txt
# Rscript $scriptDir/appendFPKM.R -i $dataDir/New_Genome/u.gibba_NEW.genePairs.txt -f $cuffDir/SRR094438_U.gibba_NEW/genes.fpkm_tracking -o $analysisDir/SRR094438_U.gibba_NEW.genePairs.foldChange.txt -p $analysisDir/SRR094438_u.gibba_NEW.foldChange.png
# Rscript $scriptDir/appendFPKM.R -i $dataDir/New_Genome/u.gibba_NEW.genePairs.txt -f $cuffDir/SRR768657_U.gibba_NEW/genes.fpkm_tracking -o $analysisDir/SRR768657_u.gibba_NEW.genePairs.foldChange.txt -p $analysisDir/SRR768657_u.gibba_NEW.foldChange.png
# Rscript $scriptDir/appendFPKM.R -i $dataDir/Old_Genome/u.gibba_OLD.genePairs.txt -f $cuffDir/SRR094438_U.gibba_OLD/genes.fpkm_tracking -o $analysisDir/SRR094438_u.gibba_OLD.genePairs.foldChange.txt -p $analysisDir/SRR094438_u.gibba_OLD.foldChange.png
# Rscript $scriptDir/appendFPKM.R -i $dataDir/Old_Genome/u.gibba_OLD.genePairs.txt -f $cuffDir/SRR768657_U.gibba_OLD/genes.fpkm_tracking -o $analysisDir/SRR768657_u.gibba_OLD.genePairs.foldChange.txt -p $analysisDir/SRR768657_u.gibba_OLD.foldChange.png

#find shared high fold change in expression gene pairs between RNA-seq runs aligned to the same genome
# Rscript $scriptDir/findSharedExpressedGenePairs.R -r1 $analysisDir/SRR768657_u.gibba_NEW.genePairs.foldChange.txt -r2 $analysisDir/SRR094438_U.gibba_NEW.genePairs.foldChange.txt -o $analysisDir/shared_u.gibba_NEW.genePairs.foldChange.txt
# Rscript $scriptDir/findSharedExpressedGenePairs.R -r1 $analysisDir/SRR768657_u.gibba_OLD.genePairs.foldChange.txt -r2 $analysisDir/SRR094438_u.gibba_OLD.genePairs.foldChange.txt -o $analysisDir/shared_u.gibba_OLD.genePairs.foldChange.txt

#confirm shared regions in old vs. new genome
#first need to make a key of which genes in old genome match genes in new genome - using blast. only taking best hit.

#make blast database using new genome
makeblastdb -dbtype nucl -in $dataDir/New_Genome/Utricularia_gibba_v2.faa -out $dataDir/New_Genome/u.gibba_NEW

#query old genome against new genome



#pull out intergenic regions and search for motifs





