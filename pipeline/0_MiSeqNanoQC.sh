#0_MiSeqNanoQC.sh

parentDir=/scratch/lk82153/jwlab/Bladderwort
qcDir=$parentDir/0_MiSeqNanoQC
fastqDir=$parentDir/0_Data/MiSeqNano_run
rnaFastqDir=$parentDir/0_Data/Fastq

krakenDB=$parentDir/../kraken_db2
refDir=$parentDir/0_Data/New_Genome
refGenome=$refDir/Utricularia_gibba_v2.faa

if [ ! -e $qcDir ]; then mkdir $qcDir; fi

# first transferred files from basespace mount to MiSeqNano_run
module load Kraken2


kraken2 --db $krakenDB --threads 48 --unclassified-out $fastqDir/unclassified_sequences#_kivancDB.fa --output $qcDir/MiSeqNano.kraken_kivancDB.out --report $qcDir/MiSeqNano.kraken_kivancDB.report --paired --gzip-compressed $fastqDir/uGibbaTest1_S1_L001_R1_001.fastq.gz $fastqDir/uGibbaTest1_S1_L001_R2_001.fastq.gz

# #align reads to bladderwort ref: how many align???
bowtie2-build $refGenome $refDir/u.gibba
bowtie2 -p 8 --phred33 -x $refDir/u.gibba -1 $fastqDir/uGibbaTest1_S1_L001_R1_001.fastq.gz -2 $fastqDir/uGibbaTest1_S1_L001_R2_001.fastq.gz -S $qcDir/uGibbaTest1_uGibbaAlignment.sam
samtools view -b $qcDir/uGibbaTest1_uGibbaAlignment.sam -o $qcDir/uGibbaTest1_uGibbaAlignment.bam

#align reads to genome - local
bowtie2 -p 8 --local --phred33 -x $refDir/u.gibba -1 $fastqDir/uGibbaTest1_S1_L001_R1_001.fastq.gz -2 $fastqDir/uGibbaTest1_S1_L001_R2_001.fastq.gz -S $qcDir/uGibbaTest1_uGibbaAlignment_local.sam
samtools view -b $qcDir/uGibbaTest1_uGibbaAlignment_local.sam -o $qcDir/uGibbaTest1_uGibbaAlignment_local.bam
samtools sort $qcDir/uGibbaTest1_uGibbaAlignment_local.bam > $qcDir/uGibbaTest1_uGibbaAlignment_local.sorted.bam
samtools index -b $qcDir/uGibbaTest1_uGibbaAlignment_local.sorted.bam

# # next i want to compare read mapping statistics of our genome with another genotype (ibarra-laclette 2013 illumina reads).
# # download SRA data
for file in SRR768651 SRR768652; do 
  fastq-dump -v --gzip $fastqDir/$file
done

#align ibarra-laclette reads to genome - end to end
fastq_pair $fastqDir/SRR768651.fastq.gz $fastqDir/SRR768652.fastq.gz
bowtie2 -p 8 --phred33 -x $refDir/u.gibba -U $fastqDir/SRR768651.fastq -S $qcDir/ibarra-laclette2013_uGibbaAlignment.sam
samtools view -b $qcDir/ibarra-laclette2013_uGibbaAlignment.sam -o $qcDir/ibarra-laclette2013_uGibbaAlignment.bam
samtools sort $qcDir/ibarra-laclette2013_uGibbaAlignment.bam > $qcDir/ibarra-laclette2013_uGibbaAlignment.sorted.bam
samtools index -b $qcDir/ibarra-laclette2013_uGibbaAlignment.sorted.bam

bowtie2 -p 8 --local --phred33 -x $refDir/u.gibba -U $fastqDir/SRR768651.fastq -S $qcDir/ibarra-laclette2013_uGibbaAlignment.local.sam
samtools view -b $qcDir/ibarra-laclette2013_uGibbaAlignment.local.sam -o $qcDir/ibarra-laclette2013_uGibbaAlignment.local.bam
samtools sort $qcDir/ibarra-laclette2013_uGibbaAlignment.local.bam > $qcDir/ibarra-laclette2013_uGibbaAlignment.local.sorted.bam
samtools index -b $qcDir/ibarra-laclette2013_uGibbaAlignment.local.sorted.bam

# # perform kraken analysis on set of reads from ibarra laclette paper 
kraken2 --db $krakenDB --threads 48 --unclassified-out $fastqDir/unclassified_sequences_SRR768651.fq --output $qcDir/MiSeqNano.kraken_SRR768651.out --report $qcDir/MiSeqNano.kraken_SRR768651.report $fastqDir/SRR768651.fastq

# # perform kraken analysis on set of RNA reads
kraken2 --db $krakenDB --threads 4 --output $qcDir/MiSeqNano.kraken_SRR768657.out --report $qcDir/MiSeqNano.kraken_SRR768657.report $rnaFastqDir/SRR768657.fastq
