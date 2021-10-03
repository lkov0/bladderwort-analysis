######################
# Author: Lynsey Kovar
# Date: Jan 16 2019
# Purpose: to align fastq files to U. gibba genome from 3' mRNA sequencing
# Software needed: FastQC, MultiQC, STAR, Trimmomatic, SAMtools, HTSeq
######################

fastqDir=$1
alignmentDir=$2
refGenome=$3
dataDir=$4
PROCS=$5

genePairs=$dataDir/genomes/utricularia/u.gibba_NEW.genePairs.txt

module load FastQC
module load MultiQC
module load STAR
module load Trimmomatic/0.36-Java-1.8.0_144
module load SAMtools
module load HTSeq

# quality check 3' rna-seq reads
mkdir $fastqDir/fastqc
fastqc -t $PROCS -o $fastqDir/fastqc $fastqDir/*fastq.gz

cd $fastqDir/fastqc/
multiqc $fastqDir/fastqc/* 

# trim reads with trimmomatic
for seq in $(ls $fastqDir | grep .fastq.gz | sed "s/.fastq.gz//g"); do
    java -jar /usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/trimmomatic-0.36.jar SE $fastqDir/$seq.fastq.gz $fastqDir/$seq.trimmed.fq.gz ILLUMINACLIP:$HOME/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -threads $PROCS 
done

# fastqc/multiqc for trimmed reads
if [ ! -e $fastqDir/fastqc_trimmed ] ; then mkdir $fastqDir/fastqc_trimmed; fi
fastqc -t $PROCS -o $fastqDir/fastqc_trimmed $fastqDir/*trimmed.fq.gz
multiqc -f -o $fastqDir/multiqc_trimmed $fastqDir/fastqc_trimmed

# comine lanes of data for each sample
for sample in 1L 2L 3L 1R 2R 3R 1B 2B 3B 1S 2S 3S; do 
    cat ${sample}*.trimmed.fq.gz > ${sample}.combined.trimmed.fq.gz
done

# align 3' transcripts with star
STAR_genome=$dataDir/genomes/STAR_genome_pbio
mkdir $STAR_genome
STAR --runMode genomeGenerate --runThreadN $PROCS --genomeDir $STAR_genome --genomeFastaFiles $refGenome --genomeSAindexNbases 12 --sjdbGTFfile $dataDir/genomes/utricularia/u.gibba_NEW.genes_for_3prime.gff

fcDir=$alignmentDir/featureCounts_out
htDir=$alignmentDir/htSeq_out
macDir=$alignmentDir/macs2_out

if [ ! -e $fcDir ] ; then mkdir $fcDir; fi
if [ ! -e $htDir ] ; then mkdir $htDir; fi
if [ ! -e $macDir ] ; then mkdir $macDir; fi

for seq in $(ls $fastqDir | grep combined.trimmed.fq.gz | sed "s/.combined.trimmed.fq.gz//g"); do
    
    STAR --runThreadN $PROCS --genomeDir $STAR_genome --sjdbOverhang 100 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 --outFilterMismatchNmax 2 --readFilesIn $fastqDir/$seq.combined.trimmed.fq.gz --readFilesCommand zcat --outFileNamePrefix $alignmentDir/${seq}_multimapTroubleshooting --outReadsUnmapped $alignmentDir/unmapped_${seq} --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts
        
    samtools index $alignmentDir/${seq}_Aligned.sortedByCoord.out.bam
    
    samtools index $alignmentDir/${seq}_multimapTroubleshootingAligned.sortedByCoord.out.bam
    
    htseq-count -f bam -t gene -i ID $alignmentDir/PacBio_Genome/${seq}_multimapTroubleshootingAligned.sortedByCoord.out.bam $dataDir/genomes/utricularia/u.gibba_NEW.genes_for_3prime.gff > $alignmentDir/PacBio_Genome/htSeq_out/${seq}_htSeq_counts.txt
    
    htseq-count -f bam -t gene -i ID $alignmentDir/PacBio_Genome/${seq}_multimapTroubleshootingAligned.sortedByCoord.out.bam $dataDir/genomes/utricularia/u.gibba_NEW.genic.gff > $alignmentDir/PacBio_Genome/htSeq_out/${seq}_htSeq_counts_noAdjustment.txt
done

# get regions corresponding to genic region + "3'" for each gene (u.gibba_NEW.genes_for_3prime.txt)
Rscript $scriptDir/expandGenicRegions3prime.R

for file in $(ls $htDir); do
    awk '{print $2}' $file.txt > $file.counts;
done

# make gff
grep "##" $dataDir/New_Genome/u.gibba_NEW.gff > u.gibba_NEW.gffheader.txt 
cat u.gibba_NEW.gffheader.txt u.gibba_NEW.genes_for_3prime.txt > u.gibba_NEW.genes_for_3prime.gff

$parentDir/obtainGenicRegions.sh $dataDir $scriptDir u.gibba_NEW.genes_for_3prime.gff

#######################
# Tests - alignment and quantification (to see which annotations capture the most data from the three full transcript mRNA-Seq samples)
#######################

# map mRNA-seq reads to genome
# STAR_genome=$dataDir/genomes/STAR_genome_ugibba_lk
# myGenome=$dataDir/genomes/utricularia/scaffolds.ugibba_lk.fasta
# myAnnos=$dataDir/genomes/utricularia/scaffolds.ugibba_lk.gff
# htAnnos=$dataDir/genomes/utricularia/scaffolds.ugibba_lk.genic.htSeq.gff

# mkdir $STAR_genome
# STAR --runMode genomeGenerate --runThreadN $PROCS --genomeDir $STAR_genome --genomeFastaFiles $myGenome --genomeSAindexNbases 12 --sjdbGTFfile $myAnnos --sjdbOverhang 149

# STAR_genome_test=$dataDir/genomes/STAR_genome_ugibba_lk_test
# mkdir $STAR_genome_test
# STAR --runMode genomeGenerate --runThreadN $PROCS --genomeDir $STAR_genome_test --genomeFastaFiles $myGenome --genomeSAindexNbases 12
# STAR --runThreadN $PROCS --genomeDir $STAR_genome --sjdbOverhang 149 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 --outFilterMismatchNmax 2 --readFilesIn $dataDir/mRNA-seq/tank1_forward_paired.fq.gz $dataDir/mRNA-seq/tank1_reverse_paired.fq.gz --readFilesCommand zcat --outFileNamePrefix $alignmentDir/tank1_ugibba_lk --outReadsUnmapped $alignmentDir/unmapped_tank1 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts
# STAR --runThreadN $PROCS --genomeDir $STAR_genome_test --sjdbOverhang 149 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 --outFilterMismatchNmax 2 --readFilesIn $dataDir/mRNA-seq/tank1_forward_paired.fq.gz $dataDir/mRNA-seq/tank1_reverse_paired.fq.gz --readFilesCommand zcat --outFileNamePrefix $alignmentDir/tank1_ugibba_lk_test --outReadsUnmapped $alignmentDir/unmapped_tank1_test --outSAMtype BAM SortedByCoordinate
# STAR --runThreadN $PROCS --genomeDir $STAR_genome --sjdbOverhang 149 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 --outFilterMismatchNmax 2 --readFilesIn $dataDir/mRNA-seq/tank2_forward_paired.fq.gz $dataDir/mRNA-seq/tank2_reverse_paired.fq.gz --readFilesCommand zcat --outFileNamePrefix $alignmentDir/tank2_ugibba_lk --outReadsUnmapped $alignmentDir/unmapped_tank2 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts
# STAR --runThreadN $PROCS --genomeDir $STAR_genome --sjdbOverhang 149 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 --outFilterMismatchNmax 2 --readFilesIn $dataDir/mRNA-seq/tank3_forward_paired.fq.gz $dataDir/mRNA-seq/tank3_reverse_paired.fq.gz --readFilesCommand zcat --outFileNamePrefix $alignmentDir/tank3_ugibba_lk --outReadsUnmapped $alignmentDir/unmapped_tank3 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts

# for seq in $(ls $fastqDir | grep combined.trimmed.fq.gz | sed "s/.combined.trimmed.fq.gz//g"); do
#     STAR --runThreadN $PROCS --genomeDir $STAR_genome --sjdbOverhang 149 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 --outFilterMismatchNmax 2 --readFilesIn $fastqDir/$seq.combined.trimmed.fq.gz --readFilesCommand zcat --outFileNamePrefix $alignmentDir/${seq}_ugibba_lk --outReadsUnmapped $alignmentDir/unmapped_${seq} --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts
#     samtools index $alignmentDir/${seq}_ugibba_lkAligned.sortedByCoord.out.bam
#     htseq-count -f bam -t gene -i gene_id $alignmentDir/${seq}_ugibba_lkAligned.sortedByCoord.out.bam $htAnnos > $htDir/${seq}_htSeq_counts.txt
#     htseq-count -f bam -t gene -i gene_id $alignmentDir/${seq}_ugibba_lkAligned.sortedByCoord.out.bam /scratch/lk82153/jwlab/Bladderwort_3pri/1_Assembly/pacbio_anno_coordinates.gff > $htDir/${seq}_htSeq_counts_pbioAnnos.txt
#     htseq-count -f bam -t gene -i gene_id $alignmentDir/${seq}_ugibba_lkAligned.sortedByCoord.out.bam /scratch/lk82153/jwlab/Bladderwort_3pri/1_Assembly/pacbio_anno_coordinates.gt5k.gff > $htDir/${seq}_htSeq_counts_pbioAnnos.gt5k.txt
# done


# multiqc -o multiqc_pbioAnnos *pbioAnnos.txt
# multiqc -o multiqc_ALL *ALL.txt
# multiqc -o multiqc_ALL.includingPacBio *ALL.includingPacBio
