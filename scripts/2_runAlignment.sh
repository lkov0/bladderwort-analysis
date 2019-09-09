######################
# Author: Lynsey Kovar
# Date: Jan 16 2018
# Purpose: to align fastq files to U. gibba genome
# Software needed: gsnap
######################

fastqDir=$1
alignmentDir=$2

for genome in $dataDir/New_Genome/Utricularia_gibba_v2.faa $dataDir/Old_Genome/Utricularia_gibba.4.1.fa; do
gmap_build -d U.gibba_NEW $dataDir/New_Genome/Utricularia_gibba_v2.faa -D $dataDir/New_Genome/
gmap_build -d U.gibba_OLD $dataDir/Old_Genome/Utricularia_gibba.4.1.fa -D $dataDir/Old_Genome/

$parentDir/obtainGenicRegions.sh $dataDir $scriptDir $dataDir/New_Genome/u.gibba_NEW.gff

for transcripts in $fastqDir/SRR094438.fastq.gz $fastqDir/SRR768657.fastq.gz; do
    for genome in $dataDir/New_Genome/U.gibba_NEW $dataDir/Old_Genome/U.gibba_OLD; do
        stem=${transcripts/$fastqDir\/}
        stem=${stem/.fastq.gz/}
        genomeDir=${genome/U.gibba_*/}
        genomeName=${genome/\/home\/lynseykovar\/Work\/Bladderwort\/0_Data\/*_Genome\//}
        gsnap -D $genomeDir -d $genomeName --batch=5 --novelsplicing 1 --nthreads $(($threads-2)) --ordered --format sam --output-file $alignmentDir/${stem}_${genomeName}.gsnap.sam --gunzip $transcripts;
        samtools sort -o $alignmentDir/${stem}_${genomeName}.gsnap.bam $alignmentDir/${stem}_${genomeName}.gsnap.sam
        samtools index $alignmentDir/${stem}_${genomeName}.gsnap.bam
        rm $alignmentDir/${stem}_${genomeName}.gsnap.sam
    done
done