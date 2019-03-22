###########################
# Author: Lynsey Kovar
# Date: Jan 16 2019
# Purpose: Parent shell script for bladderwort analysis.
# Software needed: gmap (v2018-07-04), gsnap (v2018-07-04), samtools (v1.6 using htslib 1.6), cufflinks (v2.2.1), ncbi blast (v2.7.1+), meme suite (5.0.4), mummer (v4.0.0.beta2)
###########################

parentDir=/scratch/lk82153/jwlab/Bladderwort
scriptDir=/scratch/lk82153/jwlab/Repositories/bladderwort-analysis/scripts
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

# PATH="/home/lynseykovar/Programs/sratoolkit.2.9.2-centos_linux64/bin/:$PATH"
 
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
# Rscript $scriptDir/appendFPKM.R -i $dataDir/New_Genome/u.gibba_NEW.genePairs.txt -f $cuffDir/SRR094438_U.gibba_NEW/genes.fpkm_tracking -o $analysisDir/SRR094438_U.gibba_NEW.genePairs.foldChange.ALL.txt -p $analysisDir/SRR094438_u.gibba_NEW.foldChange.ALL.png
# Rscript $scriptDir/appendFPKM.R -i $dataDir/New_Genome/u.gibba_NEW.genePairs.txt -f $cuffDir/SRR768657_U.gibba_NEW/genes.fpkm_tracking -o $analysisDir/SRR768657_u.gibba_NEW.genePairs.foldChange.ALL.txt -p $analysisDir/SRR768657_u.gibba_NEW.foldChange.ALL.png
# Rscript $scriptDir/appendFPKM.R -i $dataDir/Old_Genome/u.gibba_OLD.genePairs.txt -f $cuffDir/SRR094438_U.gibba_OLD/genes.fpkm_tracking -o $analysisDir/SRR094438_u.gibba_OLD.genePairs.foldChange.ALL.txt -p $analysisDir/SRR094438_u.gibba_OLD.foldChange.ALL.png
# Rscript $scriptDir/appendFPKM.R -i $dataDir/Old_Genome/u.gibba_OLD.genePairs.txt -f $cuffDir/SRR768657_U.gibba_OLD/genes.fpkm_tracking -o $analysisDir/SRR768657_u.gibba_OLD.genePairs.foldChange.ALL.txt -p $analysisDir/SRR768657_u.gibba_OLD.foldChange.ALL.png

#find shared high fold change in expression gene pairs between RNA-seq runs aligned to the same genome, minimum 50 fold change
# Rscript $scriptDir/findSharedExpressedGenePairs.R -r1 $analysisDir/SRR768657_u.gibba_NEW.genePairs.foldChange.ALL.txt -r2 $analysisDir/SRR094438_U.gibba_NEW.genePairs.foldChange.ALL.txt -o $analysisDir/shared_u.gibba_NEW.genePairs.foldChange.txt -f 50.0
# Rscript $scriptDir/findSharedExpressedGenePairs.R -r1 $analysisDir/SRR768657_u.gibba_OLD.genePairs.foldChange.ALL.txt -r2 $analysisDir/SRR094438_u.gibba_OLD.genePairs.foldChange.ALL.txt -o $analysisDir/shared_u.gibba_OLD.genePairs.foldChange.txt -f 50.0

#make fasta files of genes for each genome using bedtools - include feature name as name
# Rscript $scriptDir/gffToBed.R -i $dataDir/New_Genome/u.gibba_NEW.genic.gff -o $dataDir/New_Genome/u.gibba_NEW.genic.bed
# Rscript $scriptDir/gffToBed.R -i $dataDir/Old_Genome/u.gibba_OLD.genic.gff -o $dataDir/Old_Genome/u.gibba_OLD.genic.bed
# bedtools getfasta -name -fi $dataDir/New_Genome/Utricularia_gibba_v2.faa -bed $dataDir/New_Genome/u.gibba_NEW.genic.bed -fo $dataDir/New_Genome/u.gibba_NEW.genic.fa
# bedtools getfasta -name -fi $dataDir/Old_Genome/Utricularia_gibba.4.1.fa -bed $dataDir/Old_Genome/u.gibba_OLD.genic.bed -fo $dataDir/Old_Genome/u.gibba_OLD.genic.fa 

# Rscript $scriptDir/

# Make map file for identical genes between genomes #
#####################################################

#make blast database using new genome
# makeblastdb -dbtype nucl -in $dataDir/New_Genome/Utricularia_gibba_v2.faa -out $dataDir/New_Genome/u.gibba_NEW
# 
# #query old genome against new genome
# blastn -db $dataDir/New_Genome/u.gibba_NEW -query $dataDir/Old_Genome/u.gibba_OLD.genic.fa -num_threads 8 -perc_identity 95 -num_alignments 1 -out $analysisDir/u.gibba_OLD.u.gibba_NEW.blastOut.txt -outfmt 6 

#generate bedfile of genomic regions hit by genes
# Rscript $scriptDir/blastToBed.R -i $analysisDir/u.gibba_OLD.u.gibba_NEW.blastOut.txt -o $analysisDir/u.gibba_NEW.blastHits.bed

#find which of the regions in the new genome containing a hit from the old genome correspond to genes
# bedtools intersect -loj -a $analysisDir/u.gibba_NEW.blastHits.bed -b $dataDir/New_Genome/u.gibba_NEW.genic.bed > $analysisDir/u.gibba_NEW.blastHits2genes.txt

#make gene map from oldGenome:newGenome - note: not all genes in old genome were found in new genome and some genes were found 
# Rscript $scriptDir/makeGeneMap.R -i $analysisDir/u.gibba_NEW.blastHits2genes.txt -b $analysisDir/u.gibba_OLD.u.gibba_NEW.blastOut.txt -o $dataDir/geneMap_u.gibba_OLD_u.gibba_NEW.txt

# find gene pairs showing high fold change expression in each genome #
######################################################################

#using output from findSharedExpressedGenePairs.R, find gene pairs in this list shared between genomes

# Rscript $scriptDir/getGenePairsSharedInGenomes.R -d $analysisDir/shared_u.gibba_NEW.genePairs.foldChange.txt -q $analysisDir/shared_u.gibba_OLD.genePairs.foldChange.txt -m $dataDir/geneMap_u.gibba_OLD_u.gibba_NEW.txt -o matched 


###############
# getting different FPKM values for jason's and my alignment (SRR768657 to u.gibba OLD) going to realign as a test and find correlation with my original results
###############
# for transcripts in $fastqDir/SRR768657.fastq.gz; do
#     for genome in $dataDir/Old_Genome/U.gibba_OLD; do
#         stem=${transcripts/$fastqDir\/}
#         stem=${stem/.fastq.gz/}
#         genomeDir=${genome/U.gibba_*/}
#         genomeName=${genome/\/home\/lynseykovar\/Dropbox\/Bladderwort\/0_Data\/*_Genome\//}
#         gsnap -D $genomeDir -d $genomeName --batch=5 --novelsplicing 1 --nthreads $(($threads-2)) --ordered --format sam --output-file $alignmentDir/${stem}_${genomeName}.gsnap.test.sam --gunzip $transcripts;
#         samtools sort -o $alignmentDir/${stem}_${genomeName}.gsnap.test.bam $alignmentDir/${stem}_${genomeName}.gsnap.test.sam
#         samtools index $alignmentDir/${stem}_${genomeName}.gsnap.test.bam
#         rm $alignmentDir/${stem}_${genomeName}.gsnap.test.sam
#     done
# done


# cufflinks --GTF $dataDir/Old_Genome/u.gibba_OLD.genic.gff -o $cuffDir/SRR768657_U.gibba_OLD_test $alignmentDir/SRR768657_U.gibba_OLD.gsnap.test.bam

# Rscript $scriptDir/appendFPKM.R -i $dataDir/Old_Genome/u.gibba_OLD.genePairs.txt -f $cuffDir/SRR768657_U.gibba_OLD_test/genes.fpkm_tracking -o $analysisDir/SRR768657_u.gibba_OLD.genePairs.foldChange.ALL.test.txt -p $analysisDir/SRR768657_u.gibba_OLD.foldChange.ALL.test.png

#FPKMs are perfectly correlated.... hmmm. 

#only difference i can see is that jason is using gtf as input for cufflinks instead of gtf. will use exact GTF file jason used (includes, genes, mRNA, exons, etc)
# gffread $dataDir/Old_Genome/Other_GFF/Utricularia_gibba.4.1.gff3 -T -o $dataDir/Old_Genome/u.gibba_OLD.gtf

# cufflinks --GTF $dataDir/Old_Genome/u.gibba_OLD.gtf -o $cuffDir/SRR768657_U.gibba_OLD_test_gtf $alignmentDir/SRR768657_U.gibba_OLD.gsnap.test.bam
# Rscript $scriptDir/appendFPKM.R -i $dataDir/Old_Genome/u.gibba_OLD.genePairs.txt -f $cuffDir/SRR768657_U.gibba_OLD_test_gtf/genes.fpkm_tracking -o $analysisDir/SRR768657_u.gibba_OLD.genePairs.foldChange.ALL.test.gtf.txt -p $analysisDir/SRR768657_u.gibba_OLD.foldChange.ALL.test.gtf.png

#ok, there is a better correlation between jason's and my results now. See plots directory. There are only ~20 candidate genes he found that I didn't find, and about 40 unique to my dataset. (candidate gene pairs are those with >100fpkm difference in expression AND distance <= 1kb) I think this difference has to do with cufflinks version discrepencies since our alignments essentially looked the same at the loci that were different between our cufflinks output files.

#make bedfile with coverage for all genic positions
covDir=$analysisDir/Coverage

if [ ! -e $covDir ]; then mkdir $covDir; fi

# samtools sort -o $alignmentDir/SRR094438_U.gibba_NEW.gsnap.sorted.bam $alignmentDir/SRR094438_U.gibba_NEW.gsnap.bam 
# samtools sort -o $alignmentDir/SR7R68657_U.gibba_NEW.gsnap.sorted.bam $alignmentDir/SRR768657_U.gibba_NEW.gsnap.bam

# bedtools genomecov -d -ibam $alignmentDir/SRR094438_U.gibba_NEW.gsnap.sorted.bam > $analysisDir/Coverage/SRR094438_U.gibba_NEW.genomecov.txt
# bedtools genomecov -d -ibam $alignmentDir/SRR768657_U.gibba_NEW.gsnap.sorted.bam > $analysisDir/Coverage/SRR768657_U.gibba_NEW.genomecov.txt

#get average coverage in intergenic regions. 

################
# motif finding
################

#####meme

#make bedfile of candidate intergenic sequences
# Rscript $scriptDir/makeIntergenicBed.R -i $analysisDir/shared_u.gibba_NEW.genePairs.foldChange.txt -o $analysisDir/u.gibba_NEW_candidateRegions

#obtain fasta file of intergenic sequences
# bedtools getfasta -name -fi $dataDir/New_Genome/Utricularia_gibba_v2.fa -bed $analysisDir/u.gibba_NEW_candidateRegions.intergenic.divergent.bed -fo $dataDir/New_Genome/u.gibba_NEW_candidateRegions.intergenic.divergent.fasta

# bedtools getfasta -name -fi $dataDir/New_Genome/Utricularia_gibba_v2.fa -bed $analysisDir/u.gibba_NEW_candidateRegions.intergenic.convergent.bed -fo $dataDir/New_Genome/u.gibba_NEW_candidateRegions.intergenic.convergent.fasta
# 
# bedtools getfasta -name -fi $dataDir/New_Genome/Utricularia_gibba_v2.fa -bed $analysisDir/u.gibba_NEW_candidateRegions.intergenic.parallel.bed -fo $dataDir/New_Genome/u.gibba_NEW_candidateRegions.intergenic.parallel.fasta


#using meme to find overrepresented sequences in intergenic sequences - just divergent for now since these should contain enrichment of insulator elements.
# meme $dataDir/New_Genome/u.gibba_NEW_candidateRegions.intergenic.divergent.fasta -oc $analysisDir/meme_divergent -dna -p 8 -nmotifs 20

####mummer

###finding conserved sequences in intergenic regions. using grape, mimulus, papaya, tomato, and arabidopsis
module load BLAST+
blastDir=$analysisDir/Blast
if [ ! -e $blastDir ]; then mkdir $blastDir; fi


makeblastdb -dbtype nucl -in $dataDir/New_Genome/u.gibba_NEW_candidateRegions.intergenic.parallel.fasta -out $dataDir/New_Genome/intergenic.parallel
makeblastdb -dbtype nucl -in $dataDir/New_Genome/u.gibba_NEW_candidateRegions.intergenic.convergent.fasta -out $dataDir/New_Genome/intergenic.convergent
makeblastdb -dbtype nucl -in $dataDir/New_Genome/u.gibba_NEW_candidateRegions.intergenic.divergent.fasta -out $dataDir/New_Genome/intergenic.divergent

for genome in $(ls $dataDir/Genomes | sed "s/.fna//g"); do
blastn -db $dataDir/New_Genome/intergenic.parallel -query $dataDir/Genomes/$genome.fna -num_threads 8 -perc_identity 95 -num_alignments 1 -out $blastDir/$genome.parallel.blastOut.txt -outfmt "6 std qlen"
blastn -db $dataDir/New_Genome/intergenic.convergent -query $dataDir/Genomes/$genome.fna -num_threads 8 -perc_identity 95 -num_alignments 1 -out $blastDir/$genome.convergent.blastOut.txt -outfmt "6 std qlen"
blastn -db $dataDir/New_Genome/intergenic.divergent -query $dataDir/Genomes/$genome.fna -num_threads 8 -perc_identity 95 -num_alignments 1 -out $blastDir/$genome.divergent.blastOut.txt -outfmt "6 std qlen"
done


###############
# comparing blast output
###############

# blastn -db $dataDir/New_Genome/u.gibba_NEW -query $dataDir/Old_Genome/u.gibba_OLD.genic.fa -num_threads 8 -perc_identity 95 -num_alignments 1 -out $analysisDir/u.gibba_OLD.u.gibba_NEW.blastOut.forAnalysis.txt -outfmt "6 std qlen"
