#8_runGenomeAssembly.sh

# parentDir=~/Xfer/Bladderwort
parentDir=/scratch/lk82153/jwlab/Bladderwort
assemblyDir=$parentDir/8_GenomeAssembly
dataDir=$parentDir/0_Data
dnaDir=$dataDir/Bladderwort_Illumina
rnaDir=$dataDir/RNA-seq/mRNA
reference=$dataDir/New_Genome/Utricularia_gibba_v2.fa
index=$dataDir/New_Genome/Utricularia_gibba_v2
PROCS=16

# module load Bowtie2
# module load TopHat
# module load FastQC
# module load MultiQC
# module load seqtk
# module load Kraken2
# module load Trimmomatic/0.36-Java-1.8.0_144
# module load STAR
# module load spades
# module load QUAST
# module load seqtk
# module load Trinity
# module load Java/1.7.0_80
# Module load Maker
# module load BLAST+
# module load BWA
# module load SAMtools
module load OpenMPI
module load Maker


#first, run alignment of DNA and RNA-seq reads to genome

#fastqc
# mkdir $rnaDir/fastqc #because for some reason fastqc can't just make this directory itself. >:(
# fastqc -t $PROCS -o $rnaDir/fastqc $rnaDir/*.fastq.gz
# multiqc $rnaDir/fastqc/* 

#trim reads since fastqc indicated adapter content issues and low quality bases
# this command removes illumina adapters as provided in TruSeq3-PE.fa file. Remove leading and trailing low quality or N bases < quality 3. Cut when average quality per base drops below 15. Drop reads less than 36bp long after these steps.

# java -jar /usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/trimmomatic-0.36.jar PE $rnaDir/1ABladmRNA_S1_L001_R1_001.fastq.gz $rnaDir/1ABladmRNA_S1_L001_R2_001.fastq.gz $rnaDir/tank1_forward_paired.fq.gz $rnaDir/tank1_forward_unpaired.fq.gz $rnaDir/tank1_reverse_paired.fq.gz $rnaDir/tank1_reverse_unpaired.fq.gz ILLUMINACLIP:$HOME/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -threads $PROCS 
# 
# java -jar /usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/trimmomatic-0.36.jar PE $rnaDir/2ABladmRNA_S2_L001_R1_001.fastq.gz $rnaDir/2ABladmRNA_S2_L001_R2_001.fastq.gz $rnaDir/tank2_forward_paired.fq.gz $rnaDir/tank2_forward_unpaired.fq.gz $rnaDir/tank2_reverse_paired.fq.gz $rnaDir/tank2_reverse_unpaired.fq.gz ILLUMINACLIP:$HOME/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -threads $PROCS 
# 
# java -jar /usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/trimmomatic-0.36.jar PE $rnaDir/3ABladmRNA_S3_L001_R1_001.fastq.gz $rnaDir/3ABladmRNA_S3_L001_R2_001.fastq.gz $rnaDir/tank3_forward_paired.fq.gz $rnaDir/tank3_forward_unpaired.fq.gz $rnaDir/tank3_reverse_paired.fq.gz $rnaDir/tank3_reverse_unpaired.fq.gz ILLUMINACLIP:$HOME/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -threads $PROCS 

# # # # Note: bowtie alignment only produced ~2.7% mapping rate, but reads were untrimmed when i tried this. Quality trimming caused most of the reverse reads to be dropped completely, but will try mapping again with starr
# # # # #make bowtie index
# # # # bowtie2-build --threads $PROCS $reference $index
# # # # 
# # # # bowtie2 -p $PROCS --local -x $index -1 $dnaDir/uGibbaTest1_S1_L001_R1_001.fastq.gz -2 $dnaDir/uGibbaTest1_S1_L001_R2_001.fastq.gz -S $dnaDir/illuminaAlignment1.sam
# # # # 
# # # # tophat2 -o $rnaDir/tank1_mRNA -p $PROCS $index $rnaDir/1ABladmRNA_S1_L001_R1_001.fastq.gz $rnaDir/1ABladmRNA_S1_L001_R2_001.fastq.gz
# # # # tophat2 -o $rnaDir/tank2_mRNA -p $PROCS $index $rnaDir/1ABladmRNA_S1_L001_R1_001.fastq.gz $rnaDir/1ABladmRNA_S1_L001_R2_001.fastq.gz
# # # # tophat2 -o $rnaDir/tank3_mRNA -p $PROCS $index $rnaDir/3ABladmRNA_S3_L001_R1_001.fastq.gz $rnaDir/3ABladmRNA_S3_L001_R2_001.fastq.gz

#star alignment
STAR_genome=$dataDir/New_Genome/STAR_genome
# 
if [ ! -e $STAR_genome ]; then mkdir $STAR_genome; fi

# STAR --runMode genomeGenerate --runThreadN $PROCS --genomeDir $STAR_genome --genomeFastaFiles $reference 
# 
# STAR --genomeDir $STAR_genome --runThreadN $PROCS --readFilesIn $rnaDir/tank1_forward_paired.fq.gz $rnaDir/tank1_reverse_paired.fq.gz --readFilesCommand zcat --outFileNamePrefix $rnaDir/tank1_ --outReadsUnmapped $rnaDir/$rnaDir/unmapped_tank1_ --outSAMtype BAM SortedByCoordinate
# 
# STAR --genomeDir $STAR_genome --runThreadN $PROCS --readFilesIn $rnaDir/tank2_forward_paired.fq.gz $rnaDir/tank2_reverse_paired.fq.gz --readFilesCommand zcat --outFileNamePrefix $rnaDir/tank2_ --outReadsUnmapped $rnaDir/$rnaDir/unmapped_tank2_ --outSAMtype BAM SortedByCoordinate
# 
# STAR --genomeDir $STAR_genome --runThreadN $PROCS --readFilesIn $rnaDir/tank3_forward_paired.fq.gz $rnaDir/tank3_reverse_paired.fq.gz --readFilesCommand zcat --outFileNamePrefix $rnaDir/tank3_ --outReadsUnmapped $rnaDir/$rnaDir/unmapped_tank3_ --outSAMtype BAM SortedByCoordinate

# samtools index -@ $PROCS $rnaDir/tank1_Aligned.sortedByCoord.out.bam
# samtools index -@ $PROCS $rnaDir/tank2_Aligned.sortedByCoord.out.bam
# samtools index -@ $PROCS $rnaDir/tank3_Aligned.sortedByCoord.out.bam

#merge bam files
# samtools merge -@ $PROCS $rnaDir/bladderwortMRNA_Aligned.sortedByCoord.out.bam $rnaDir/tank1_Aligned.sortedByCoord.out.bam $rnaDir/tank2_Aligned.sortedByCoord.out.bam $rnaDir/tank3_Aligned.sortedByCoord.out.bam

#classify tophat sequences

KrakenDB=/scratch/lk82153/jwlab/kraken_db2

# seqtk seq -A $rnaDir/1ABladmRNA_S1_L001_R1_001.fastq > $rnaDir/1ABladmRNA_S1_L001_R1_001.fasta
# seqtk seq -A $rnaDir/2ABladmRNA_S2_L001_R1_001.fastq > $rnaDir/2ABladmRNA_S2_L001_R1_001.fasta
# seqtk seq -A $rnaDir/3ABladmRNA_S3_L001_R1_001.fastq > $rnaDir/3ABladmRNA_S3_L001_R1_001.fasta

# kraken2 -db $KrakenDB --output $rnaDir/tank1.kraken.out --report $rnaDir/tank1.kraken.report --threads $PROCS $rnaDir/1ABladmRNA_S1_L001_R1_001.fasta
# kraken2 -db $KrakenDB --output $rnaDir/tank2.kraken.out --report $rnaDir/tank2.kraken.report --threads $PROCS $rnaDir/2ABladmRNA_S2_L001_R1_001.fasta
# kraken2 -db $KrakenDB --output $rnaDir/tank3.kraken.out --report $rnaDir/tank3.kraken.report --threads $PROCS $rnaDir/3ABladmRNA_S3_L001_R1_001.fasta

###############
# Genome assembly
###############

#first merged reads from each flowcell into combined R1 and R2 files 
# cat $dnaDir/uGibbaTest1_S1_L001_R1_001.fastq $dnaDir/uGibbaTest1_2_S1_L001_R1_001.fastq > $dnaDir/uGibbaMerged_S1_L001_R1_001.fastq
# cat $dnaDir/uGibbaTest1_S1_L001_R2_001.fastq $dnaDir/uGibbaTest1_2_S1_L001_R2_001.fastq > $dnaDir/uGibbaMerged_S1_L001_R2_001.fastq

# running trimmomatic - no quality trimming because spades will do that for me. just removing adapter sequences. 
# java -jar /usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/trimmomatic-0.36.jar PE $dnaDir/uGibbaMerged_S1_L001_R1_001.fastq.gz $dnaDir/uGibbaMerged_S1_L001_R2_001.fastq.gz $dnaDir/uGibbaMerged_forward_trimmed_paired.fq.gz $dnaDir/uGibbaMerged_forward_trimmed_unpaired.fq.gz $dnaDir/uGibbaMerged_reverse_trimmed_paired.fq.gz $dnaDir/uGibbaMerged_reverse_trimmed_unpaired.fq.gz ILLUMINACLIP:$HOME/TruSeq3-PE.fa:2:30:10 MINLEN:36 -threads $PROCS

# fastqc -t $PROCS -o $dnaDir/fastqc $dnaDir/uGibbaMerged*.gz

# cd $dnaDir/fastqc/
# multiqc $dnaDir/fastqc/*

# spades.py --continue -1 $dnaDir/uGibbaMerged_forward_trimmed_paired.fq.gz -2 $dnaDir/uGibbaMerged_reverse_trimmed_paired.fq.gz -t $PROCS -o $assemblyDir/spades_assembly1

# spades.py --continue -o $assemblyDir/spades_assembly1

# # run quast to assess first genome assembly

# quast.py -o $assemblyDir/spades_assembly1/quast_all -r $dataDir/New_Genome/Utricularia_gibba_v2.fa -g $dataDir/New_Genome/u.gibba_NEW.genic.gff -t $PROCS --eukaryote -b --min-identity 90.0 --threads $PROCS $assemblyDir/spades_assembly1/scaffolds.fasta $assemblyDir/spades_assembly1/contigs.fasta

# # compare these results to the published illumina assembly mapped to the published pacbio assembly

# quast.py -o $assemblyDir/spades_assembly1/quast_publishedGenomes -r $dataDir/New_Genome/Utricularia_gibba_v2.fa -g $dataDir/New_Genome/u.gibba_NEW.gencd ic.gff -t $PROCS --eukaryote -b --min-identity 90.0 --threads $PROCS $dataDir/Old_Genome/Utricularia_gibba.4.1.fa

# # create blob plot of assembly - blobtools installed on local machine
# # first blast assembly against nt database on sapelo2 (need output format to be seqID, taxID and score using "grep -wFf ~/neededAccessions.list nucl.accession2taxid.txt > contigs.ugibba.accesion2taxid" to make accession to taxid map that was small enough to use for mergining in R.

# blastn -db /db/ncbiblast/nrte/latest/nt -query $assemblyDir/spades_assembly1/contigs.fasta -out $assemblyDir/spades_assembly1/contigs.blobtools.blastout -perc_identity 90 -num_threads $PROCS -max_target_seqs 10 -outfmt 6

# bwa index $assemblyDir/spades_assembly1/contigs.fasta
# bwa mem $assemblyDir/spades_assembly1/contigs.fasta $dnaDir/uGibbaMerged_forward_trimmed_paired.fq.gz $dnaDir/uGibbaMerged_reverse_trimmed_paired.fq.gz > $assemblyDir/spades_assembly1/contigs.blobtools.sorted.sam

# samtools view -S -b $assemblyDir/spades_assembly1/contigs.blobtools.sorted.sam | samtools sort - > $assemblyDir/spades_assembly1/contigs.blobtools.sorted.bam

# had to modify blast output to contain contig, taxid of hit, and bitscore as the first three columns. used accession to taxid map and grep -
# blobtools create -i $assemblyDir/spades_assembly1/contigs.fasta -b $assemblyDir/spades_assembly1/contigs.blobtools.sorted.bam -t $assemblyDir/spades_assembly1/contigs.blobtools.blastout -o $assemblyDir/spades_assembly1/contigs.blobtools && \
# blobtools view -i $assemblyDir/spades_assembly1/contigs.blobtools.blobDB.json && \
# blobtools plot -i $assemblyDir/spades_assembly1/contigs.blobtools.blobDB.json --format svg --notitle

# # run kraken on contigs and scaffolds
# kraken2 -db $KrakenDB --output $assemblyDir/spades_assembly1/contigs.kraken.out --report $assemblyDir/spades_assembly1/contigs.kraken.report --threads $PROCS $assemblyDir/spades_assembly1/contigs.fasta
# kraken2 -db $KrakenDB --output $assemblyDir/spades_assembly1/scaffolds.kraken.out --report $assemblyDir/spades_assembly1/scaffolds.kraken.report --threads $PROCS $assemblyDir/spades_assembly1/scaffolds.fasta
# 
# # # get sequences classified as bladderwort taxid 13748
# awk '{print $1, $2, $3, $4}' $assemblyDir/spades_assembly1/contigs.kraken.out | sed "s/ /\t/g" | awk '$3=="13748"' > $assemblyDir/spades_assembly1/contigs.ugibba.txt
# awk '{print $1, $2, $3, $4}' $assemblyDir/spades_assembly1/scaffolds.kraken.out | sed "s/ /\t/g" | awk '$3=="13748"' > $assemblyDir/spades_assembly1/scaffolds.ugibba.txt

# # get unclassified contigs and sequences
# awk '$1=="U"' $assemblyDir/spades_assembly1/contigs.kraken.out | awk '{print $2}'  > $assemblyDir/spades_assembly1/contigs.unclassified.txt
# seqtk subseq $assemblyDir/spades_assembly1/contigs.fasta $assemblyDir/spades_assembly1/contigs.unclassified.txt > $assemblyDir/spades_assembly1/contigs.unclassified.fasta

# # #get names of sequences classified as bladderwort, extract using seqtk
# awk '{print $2}' $assemblyDir/spades_assembly1/contigs.ugibba.txt > $assemblyDir/spades_assembly1/contig_names.ugibba.txt
# awk '{print $2}' $assemblyDir/spades_assembly1/scaffolds.ugibba.txt > $assemblyDir/spades_assembly1/scaffold_names.ugibba.txt
# 
# seqtk subseq $assemblyDir/spades_assembly1/contigs.fasta $assemblyDir/spades_assembly1/contig_names.ugibba.txt > $assemblyDir/spades_assembly1/contigs.ugibba.fasta
# seqtk subseq $assemblyDir/spades_assembly1/scaffolds.fasta $assemblyDir/spades_assembly1/scaffold_names.ugibba.txt > $assemblyDir/spades_assembly1/scaffolds.ugibba.fasta
# 
# # # run quast on just the u gibba assembly
# quast.py -o $assemblyDir/spades_assembly1/quast_ugibba_output_inclOld -r $dataDir/New_Genome/Utricularia_gibba_v2.fa -g $dataDir/New_Genome/u.gibba_NEW.genic.gff -t $PROCS --eukaryote -b --min-identity 90.0 $assemblyDir/spades_assembly1/scaffolds.ugibba.fasta $assemblyDir/spades_assembly1/contigs.ugibba.fasta

# # #try to scaffold with ragoo instead of spades on the spades ugibba contigs
# # #since all files need to be in the same directory...
# # 
# # ragooDir=$assemblyDir/Ragoo
# # # 
# # # if [ ! -e $ragooDir ]; then mkdir $ragooDir; fi
# # # cp $dataDir/New_Genome/u.gibba_NEW.genic.gff $assemblyDir/contigs.ugibba.fasta $dataDir/New_Genome/Utricularia_gibba_v2.fasta $dnaDir/uGibbaMerged_forward_trimmed_paired.fq.gz $dnaDir/uGibbaMerged_reverse_trimmed_paired.fq.gz $ragooDir
# # # 
# # # cd $ragooDir
# # # ragoo.py -gff u.gibba_NEW.genic.gff -R uGibbaMerged_combined_trimmed_paired.fq -T sr -t 7 -s contigs.ugibba.fasta Utricularia_gibba_v2.fasta
# # 
# # #evaluate ragoo assembly
# # # quast.py -o $ragooDir/quast_results -r $dataDir/New_Genome/Utricularia_gibba_v2.fa -g $dataDir/New_Genome/u.gibba_NEW.genic.gff -t $PROCS --eukaryote -b --min-identity 90.0 $ragooDir/ragoo_output/ragoo.fasta

######################
# Maker annotation on ragoo assembly
######################

#first, trinity assemble transcripts

#  Trinity --seqType fq --JM 50G --left $rnaDir/tank1_forward_paired.fq.gz,$rnaDir/tank2_forward_paired.fq.gz,$rnaDir/tank3_forward_paired.fq.gz --right $rnaDir/tank1_reverse_paired.fq.gz,$rnaDir/tank2_reverse_paired.fq.gz,$rnaDir/tank3_reverse_paired.fq.gz --CPU $PROCS --output $assemblyDir/Trinity

# kraken2 -db $KrakenDB --output $assemblyDir/Trinity/trinity.kraken.out --report $assemblyDir/Trinity/trinity.kraken.report --threads $PROCS $assemblyDir/Trinity/Trinity.fasta
 
# awk '{print $1, $2, $3, $4}' $assemblyDir/Trinity/trinity.kraken.out | sed "s/ /\t/g" | awk '$3=="13748"' > $assemblyDir/Trinity/trinity.ugibba.txt
# awk '{print $2}' $assemblyDir/Trinity/trinity.ugibba.txt > $assemblyDir/Trinity/trinity_names.ugibba.txt
# 
# seqtk subseq $assemblyDir/Trinity/Trinity.fasta $assemblyDir/Trinity/trinity_names.ugibba.txt > $assemblyDir/Trinity/Trinity.ugibba.fasta

# running maker
# edited maker_opts.ctl to include trinity assembled transcripts and asterid proteins downloaded from ensembl PLAZA
makerDir=$assemblyDir/Maker

if [ ! -e $makerDir ]; then mkdir $makerDir; fi
cd $makerDir

maker -base ugibba_rnd1 maker_opts.ctl maker_bopts.ctl maker_exe.ctl

###
# Wait, what are the unclassified transcripts - blastx
###

# cat $assemblyDir/Trinity/trinity.kraken.out | grep "^U" | awk '{print $2}' > $assemblyDir/Trinity/unclassfied_trinity.txt
# seqtk subseq $assemblyDir/spades_assembly1/contigs.fasta $assemblyDir/spades_assembly1/contigs.unclassified.txt > $assemblyDir/spades_assembly1/contigs.unclassified.fasta
# blastx -db /db/ncbiblast/nrte/latest/nr -query $assemblyDir/Trinity/unclassified.fasta -out $assemblyDir/Trinity/unclassified.blastout -num_threads $PROCS -max_target_seqs 10 -outfmt 6

# generate maker contig files
# maker -CTL 
# only changed these options:
# - genome=/scratch/lk82153/jwlab/Bladderwort/0_Data/New_Genome/Utricularia_gibba_v2.fa
# - est= 
# 
# blastn -db $dataDir/New_Genome/u.gibba_NEW.genes -query $assemblyDir/scaffolds.ugibba.fasta -out $assemblyDir/scaffolds.ugibba.geneblast.txt -perc_identity 90 -num_threads $PROCS -max_target_seqs 10 -outfmt 6
