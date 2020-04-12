######################
# Author: Lynsey Kovar
# Date: Jan 16 2019
# Purpose: preliminary bladderwort analysis using old published datasets
# Software needed: R, Blast
######################

parentDir=$1
scriptDir=$2
dataDir=$3
alignDir=$4
analysisDir=$5

myAnnos=$dataDir/genomes/utricularia/scaffolds.ugibba_lk.ALL.includingPacBio.gff
myGenome=$dataDir/genomes/utricularia/scaffolds.ugibba_lk.fasta

module load R

#Pull out gene pairs that show high fold change in expression

#pull out convergent, divergent, parallel gene pairs
Rscript $scriptDir/assignConvergentDivergentParallel.R -i $myAnnos -o $dataDir/genomes/utricularia/scaffolds.ugibba_lk.ALL.includingPacBio.genePairs.txt

# append expression values to genePairs file

# find shared high fold change in expression gene pairs between RNA-seq runs aligned to the same genome, minimum 50 fold change
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
# covDir=$analysisDir/Coverage
# 
# if [ ! -e $covDir ]; then mkdir $covDir; fi

# samtools sort -o $alignmentDir/SRR094438_U.gibba_NEW.gsnap.sorted.bam $alignmentDir/SRR094438_U.gibba_NEW.gsnap.bam 
# samtools sort -o $alignmentDir/SR7R68657_U.gibba_NEW.gsnap.sorted.bam $alignmentDir/SRR768657_U.gibba_NEW.gsnap.bam

# # bedtools genomecov -d -ibam $alignmentDir/SRR094438_U.gibba_NEW.gsnap.sorted.bam > $analysisDir/Coverage/SRR094438_U.gibba_NEW.genomecov.txt
# bedtools genomecov -d -ibam $alignmentDir/SRR768657_U.gibba_NEW.gsnap.sorted.bam > $analysisDir/Coverage/SRR768657_U.gibba_NEW.genomecov.txt

#get average coverage in intergenic regions. 



