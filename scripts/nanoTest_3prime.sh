#nanoTest_3prime.sh

# parentDir=/scratch/lk82153/jwlab/Bladderwort
parentDir=~/Xfer/Bladderwort
dataDir=$parentDir/0_Data
rnaDir=$dataDir/RNA-seq/3primeNano
reference=$dataDir/New_Genome/Utricularia_gibba_v2.fa
index=$dataDir/New_Genome/Utricularia_gibba_v2
nanoTestDir=$parentDir/0_3primeNanoQC
PROCS=48

if [ ! -e $nanoTestDir ] ; then mkdir $nanoTestDir; fi

# module load FastQC
# module load MultiQC
# module load Trimmomatic/0.36-Java-1.8.0_144
# module load STAR
# module load seqtk
# module load Kraken2
module load SAMtools

# module load Biopython

#fastqc
# mkdir $rnaDir/fastqc #because for some reason fastqc can't just make this directory itself. >:(
# if [ ! -e $rnaDir/fastqc ] ; then mkdir $rnaDir/fastqc; fi
# fastqc -t $PROCS -o $rnaDir/fastqc $rnaDir/*.fastq.gz
# multiqc -f -o $rnaDir/multiqc $rnaDir/fastqc/*

#trim reads since fastqc indicated adapter content issues and low quality bases
# this command removes illumina adapters as provided in TruSeq3-PE.fa file. Remove leading and trailing low quality or N bases < quality 3. Cut when average quality per base drops below 15. Drop reads less than 36bp long after these steps.

# for seq in $(ls $rnaDir | grep .fastq.gz | sed "s/.fastq.gz//g"); do
#     java -jar /usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/trimmomatic-0.36.jar SE $rnaDir/$seq.fastq.gz $rnaDir/$seq.trimmed.fq.gz ILLUMINACLIP:$HOME/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -threads $PROCS 
# done

# if [ ! -e $rnaDir/fastqc_trimmed ] ; then mkdir $rnaDir/fastqc_trimmed; fi
# fastqc -t $PROCS -o $rnaDir/fastqc_trimmed $rnaDir/*trimmed.fq.gz
# multiqc -f -o $rnaDir/multiqc_trimmed $rnaDir/fastqc_trimmed

# STAR_genome=$dataDir/New_Genome/STAR_genome

# if [ ! -e $STAR_genome ]; then mkdir $STAR_genome; fi
# 
# STAR --runMode genomeGenerate --runThreadN $PROCS --genomeDir $STAR_genome --genomeFastaFiles $reference 

for seq in $(ls $rnaDir | grep .fastq.gz | sed "s/.fastq.gz//g"); do
#     STAR --genomeDir $STAR_genome --runThreadN $PROCS --readFilesIn $rnaDir/$seq.trimmed.fq.gz --readFilesCommand zcat --outFileNamePrefix $nanoTestDir/${seq}_ --outReadsUnmapped $nanoTestDir/unmapped_${seq} --outSAMtype BAM SortedByCoordinate
samtools index $nanoTestDir/${seq}_Aligned.sortedByCoord.out.bam
done

# KrakenDB=/scratch/lk82153/jwlab/kraken_db2
# 
# # Kraken for read classification
# for seq in $(ls $rnaDir | grep .fastq.gz | sed "s/.fastq.gz//g"); do
#     zcat $rnaDir/$seq.fastq.gz | seqtk seq -A - > $rnaDir/$seq.fasta
#     kraken2 -db $KrakenDB --output $nanoTestDir/$seq.kraken.out --report $nanoTestDir/$seq.kraken.report --threads $PROCS $rnaDir/$seq.fasta
# done

for seq in $(ls $rnaDir | grep .fastq.gz | sed "s/.fastq.gz//g"); do
echo $seq
python3.7 ../../Repositories/my-seq-utils/my-seq-utils/gene_hits.py ${seq}_Aligned.sortedByCoord.out.bam ../0_Data/New_Genome/u.gibba_NEW.genic.bed
done
