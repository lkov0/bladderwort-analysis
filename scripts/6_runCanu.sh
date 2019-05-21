#!/bin/bash
#6_runErrorCorrectWithCanu.sh

canuDir=$1
readDir=$2
dataDir=$3
genome=$4
scriptDir=$5
nProcs=$6

# module load Java/1.8.0_181
# module load canu
# module load BLAST+
# module load QUAST
module load Kraken2

#combine reads

# cat $readDir/*.fastq > $readDir/ merged.PacBio.fastq
# canu -correct stopOnLowCoverage=8 minOverlapLength=200 -minReadLength=500 -p bladderwort_test -d $canuDir/bladderwort_test genomeSize=100m java=/usr/local/apps/eb/Java/1.8.0_144/bin/java -pacbio-raw merged.PacBio.fastq

# canu -trim stopOnLowCoverage=8 minOverlapLength=200 -minReadLength=500 -p bladderwort_test -d $canuDir/bladderwort_test genomeSize=100m java=/usr/local/apps/eb/Java/1.8.0_144/bin/java correctedErrorRate=0.105 -pacbio-corrected $canuDir/bladderwort_test/bladderwort_test.correctedReads.fasta.gz

#for assembly, using slightly higher error rate to adjust for the fact that we have lower coverage (~25x), canu quickstart guide lists an error rate of 0.105 for this
# canu -assemble stopOnLowCoverage=8 minOverlapLength=200 -minReadLength=500 -p bladderwort_test gridOptions="-l nodes=1:ppn=8,walltime=24:00:00" -d $canuDir/bladderwort_test genomeSize=100m java=/usr/local/apps/eb/Java/1.8.0_144/bin/java correctedErrorRate=0.105 -pacbio-corrected $canuDir/bladderwort_test/bladderwort_test.trimmedReads.fasta.gz


# gunzip $canuDir/bladderwort_test/bladderwort_test.correctedReads.fasta.gz

# blast error corrected reads against old and new genome
# blastn -db $dataDir/New_Genome/u.gibba_NEW -query $canuDir/bladderwort_test/bladderwort_test.correctedReads.fasta -num_threads $nProcs -perc_identity 90 -num_alignments 5 -out $canuDir/bladderwort_test.corrected.u.gibba_NEW.blastOut.txt -outfmt "6 std qlen"
# blastn -db $dataDir/New_Genome/u.gibba_OLD -query $canuDir/bladderwort_test/bladderwort_test.correctedReads.fasta -num_threads $nProcs -perc_identity 90 -num_alignments 5 -out $canuDir/bladderwort_test.corrected.u.gibba_OLD.blastOut.txt -outfmt "6 std qlen"

# using QUAST to assess assembly
# quast.py -o $canuDir/quast_results -r $dataDir/New_Genome/Utricularia_gibba_v2.fa -g $dataDir/New_Genome/u.gibba_NEW.genic.gff -t $nProcs --eukaryote -b --min-identity 90.0 $canuDir/bladderwort_test/bladderwort_test.contigs.fasta

# using Kraken to characterize unassembled contigs - first building kraken database - actually using kraken db on my work computer since sapelo seems to have trouble getting files from the NCBI server
# kraken2-build --standard --threads $nProcs --db $canuDir/KrakenDB_standardAndPlant
# kraken2-build --download-taxonomy --threads $nProcs --db $canuDir/KrakenDB_standardAndPlant
# kraken2-build --download-library "plant" --threads $nProcs --db $canuDir/KrakenDB_standardAndPlant
# kraken2-build --build --threads $nProcs --db $canuDir/KrakenDB_standardAndPlant
# kraken2-build --db $canuDir/KrakenDB_standardAndPlant --clean --threads $nProcs

#Using krakenDB from Kivanc, should have standard, plant and fungi
KrakenDB=~/Xfer/kraken_db2
# add bladderwort genome
# sed "s/>/>|kraken:taxid|13748 /g" $genome > $genome.kraken.fa
# kraken2-build --add-to-library $genome.kraken.fa --db $KrakenDB
# kraken2-build --build --threads 8 --db $KrakenDB
#ok, can't add to prebuilt database, and i don't have the library files(presumably they were deleted to save space) will try installing fresh on the cluster

# kraken2-build --download-taxonomy --threads $nProcs --db $KrakenDB
# kraken2-build --download-library human --threads $nProcs --db $KrakenDB
# kraken2-build --download-library plant --threads $nProcs --db $KrakenDB
# kraken2-build --download-library bacteria --threads $nProcs --db $KrakenDB
# kraken2-build --download-library fungi --threads $nProcs --db $KrakenDB
# kraken2-build --add-to-library $genome.kraken.fa --threads $nProcs --db $KrakenDB
# kraken2-build --build --threads $nProcs --db $KrakenDB

# run kraken on error corrected reads - 
# kraken2 -db $KrakenDB --output $canuDir/correctedReads.kraken.out --use-mpa-style --report $canuDir/correctedReads.kraken.report --threads $nProcs $canuDir/bladderwort_test/bladderwort_test.correctedReads.fasta.gz

Rscript $scriptDir/formatKrakenTaxonomyCounts.R -i $canuDir/correctedReads.kraken.report -o $canuDir/correctedReads.kraken.taxonomyCounts.txt;

