#8_runGenomeAssembly.sh

parentDir=/scratch/lk82153/jwlab/Bladderwort
scriptDir=/home/lk82153/Repositories/bladderwort-analysis/scripts
assemblyDir=$parentDir/8_GenomeAssembly
dataDir=/work/jawlab/data/bladderwort
dnaDir=$dataDir/Bladderwort_Illumina
rnaDir=$dataDir/mRNA-seq
reference=$dataDir/New_Genome/Utricularia_gibba_v2.fa
index=$dataDir/New_Genome/Utricularia_gibba_v2
PROCS=12

module load Bowtie2
module load TopHat
module load FastQC
module load MultiQC
module load seqtk
module load Kraken2
module load Trimmomatic/0.36-Java-1.8.0_144
module load STAR
module load spades
module load QUAST
module load seqtk
module load Trinity
module load Java/1.7.0_80
module load BLAST+
module load BWA
module load SAMtools
module load OpenMPI
module load Maker/2.31.10-foss-2016b
module load BEDTools
module load HMMER/2.3-foss-2016b
module load AUGUSTUS/3.2.3-foss-2016b-Python-2.7.14/module load busco/3.0.2
module load R

###############
# Genome assembly
###############

# first merged reads from each flowcell into combined R1 and R2 files 
cat $dnaDir/uGibbaTest1_S1_L001_R1_001.fastq $dnaDir/uGibbaTest1_2_S1_L001_R1_001.fastq > $dnaDir/uGibbaMerged_S1_L001_R1_001.fastq
cat $dnaDir/uGibbaTest1_S1_L001_R2_001.fastq $dnaDir/uGibbaTest1_2_S1_L001_R2_001.fastq > $dnaDir/uGibbaMerged_S1_L001_R2_001.fastq

# running trimmomatic - no quality trimming because spades will do that for me. just removing adapter sequences. 
java -jar /usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/trimmomatic-0.36.jar PE $dnaDir/uGibbaMerged_S1_L001_R1_001.fastq.gz $dnaDir/uGibbaMerged_S1_L001_R2_001.fastq.gz $dnaDir/uGibbaMerged_forward_trimmed_paired.fq.gz $dnaDir/uGibbaMerged_forward_trimmed_unpaired.fq.gz $dnaDir/uGibbaMerged_reverse_trimmed_paired.fq.gz $dnaDir/uGibbaMerged_reverse_trimmed_unpaired.fq.gz ILLUMINACLIP:$HOME/TruSeq3-PE.fa:2:30:10 MINLEN:36 -threads $PROCS

fastqc -t $PROCS -o $dnaDir/fastqc $dnaDir/uGibbaMerged*.gz

cd $dnaDir/fastqc/
multiqc $dnaDir/fastqc/*

spades.py --continue -1 $dnaDir/uGibbaMerged_forward_trimmed_paired.fq.gz -2 $dnaDir/uGibbaMerged_reverse_trimmed_paired.fq.gz -t $PROCS -o $assemblyDir/spades_assembly1
spades.py --continue -o $assemblyDir/spades_assembly1

# # run quast to assess first genome assembly
quast.py -o $assemblyDir/spades_assembly1/quast_all -r $dataDir/New_Genome/Utricularia_gibba_v2.fa -g $dataDir/New_Genome/u.gibba_NEW.genic.gff -t $PROCS --eukaryote -b --min-identity 90.0 --threads $PROCS $assemblyDir/spades_assembly1/scaffolds.fasta $assemblyDir/spades_assembly1/contigs.fasta

# # compare these results to the published illumina assembly mapped to the published pacbio assembly
quast.py -o $assemblyDir/spades_assembly1/quast_publishedGenomes -r $dataDir/New_Genome/Utricularia_gibba_v2.fa -g $dataDir/New_Genome/u.gibba_NEW.gencd ic.gff -t $PROCS --eukaryote -b --min-identity 90.0 --threads $PROCS $dataDir/Old_Genome/Utricularia_gibba.4.1.fa

# # create blob plot of assembly - blobtools installed on local machine
# # first blast assembly against nt database on sapelo2 (need output format to be seqID, taxID and score using "grep -wFf ~/neededAccessions.list nucl.accession2taxid.txt > contigs.ugibba.accesion2taxid" to make accession to taxid map that was small enough to use for mergining in R.
blastn -db /db/ncbiblast/nrte/latest/nt -query $assemblyDir/spades_assembly1/contigs.fasta -out $assemblyDir/spades_assembly1/contigs.blobtools.blastout -perc_identity 90 -num_threads $PROCS -max_target_seqs 10 -outfmt 6

bwa index $assemblyDir/spades_assembly1/contigs.fasta
bwa mem $assemblyDir/spades_assembly1/contigs.fasta $dnaDir/uGibbaMerged_forward_trimmed_paired.fq.gz $dnaDir/uGibbaMerged_reverse_trimmed_paired.fq.gz > $assemblyDir/spades_assembly1/contigs.blobtools.sorted.sam
samtools view -S -b $assemblyDir/spades_assembly1/contigs.blobtools.sorted.sam | samtools sort - > $assemblyDir/spades_assembly1/contigs.blobtools.sorted.bam

# had to modify blast output to contain contig, taxid of hit, and bitscore as the first three columns. used accession to taxid map and grep -
blobtools create -i $assemblyDir/spades_assembly1/contigs.fasta -b $assemblyDir/spades_assembly1/contigs.blobtools.sorted.bam -t $assemblyDir/spades_assembly1/contigs.blobtools.blastout -o $assemblyDir/spades_assembly1/contigs.blobtools && \
blobtools view -i $assemblyDir/spades_assembly1/contigs.blobtools.blobDB.json && \
blobtools plot -i $assemblyDir/spades_assembly1/contigs.blobtools.blobDB.json --format svg --notitle

# # run kraken on contigs and scaffolds
kraken2 -db $KrakenDB --output $assemblyDir/spades_assembly1/contigs.kraken.out --report $assemblyDir/spades_assembly1/contigs.kraken.report --threads $PROCS $assemblyDir/spades_assembly1/contigs.fasta
kraken2 -db $KrakenDB --output $assemblyDir/spades_assembly1/scaffolds.kraken.out --report $assemblyDir/spades_assembly1/scaffolds.kraken.report --threads $PROCS $assemblyDir/spades_assembly1/scaffolds.fasta
 
# # get sequences classified as bladderwort taxid 13748
awk '{print $1, $2, $3, $4}' $assemblyDir/spades_assembly1/contigs.kraken.out | sed "s/ /\t/g" | awk '$3=="13748"' > $assemblyDir/spades_assembly1/contigs.ugibba.txt
awk '{print $1, $2, $3, $4}' $assemblyDir/spades_assembly1/scaffolds.kraken.out | sed "s/ /\t/g" | awk '$3=="13748"' > $assemblyDir/spades_assembly1/scaffolds.ugibba.txt

# # get unclassified contigs and sequences
awk '$1=="U"' $assemblyDir/spades_assembly1/contigs.kraken.out | awk '{print $2}'  > $assemblyDir/spades_assembly1/contigs.unclassified.txt
seqtk subseq $assemblyDir/spades_assembly1/contigs.fasta $assemblyDir/spades_assembly1/contigs.unclassified.txt > $assemblyDir/spades_assembly1/contigs.unclassified.fasta

# # #get names of sequences classified as bladderwort, extract using seqtk
awk '{print $2}' $assemblyDir/spades_assembly1/contigs.ugibba.txt > $assemblyDir/spades_assembly1/contig_names.ugibba.txt
awk '{print $2}' $assemblyDir/spades_assembly1/scaffolds.ugibba.txt > $assemblyDir/spades_assembly1/scaffold_names.ugibba.txt

seqtk subseq $assemblyDir/spades_assembly1/contigs.fasta $assemblyDir/spades_assembly1/contig_names.ugibba.txt > $assemblyDir/spades_assembly1/contigs.ugibba.fasta
seqtk subseq $assemblyDir/spades_assembly1/scaffolds.fasta $assemblyDir/spades_assembly1/scaffold_names.ugibba.txt > $assemblyDir/spades_assembly1/scaffolds.ugibba.fasta

awk '$2 > 1000' scaffolds.ugibba.fasta.fai | cut -f1 > scaffolds.ugibba.fasta.greaterThan1000.txt

# get only scaffolds > 1kb
seqtk subseq $assemblyDir/spades_assembly1/scaffolds.fasta $assemblyDir/spades_assembly1/scaffolds.ugibba.fasta.greaterThan1000.txt > $assemblyDir/spades_assembly1/scaffolds.ugibba.greaterthan1000.fasta

# # # run quast on just the u gibba assembly
quast.py -o $assemblyDir/spades_assembly1/quast_ugibba_output_inclOld -r $dataDir/New_Genome/Utricularia_gibba_v2.fa -g $dataDir/New_Genome/u.gibba_NEW.genic.gff -t $PROCS --eukaryote -b --min-identity 90.0 $assemblyDir/spades_assembly1/scaffolds.ugibba.fasta $assemblyDir/spades_assembly1/contigs.ugibba.fasta

cp $assemblyDir/spades_assembly1/scaffolds.ugibba.greaterthan1000.fasta /work/jawlab/data/bladderwort/genomes/utricularia/scaffolds.ugibba_lk.fasta

# remove chloroplast and mitochondrial contigs (blasting contigs against U. gibba chloroplast and U. reniformis mitochondrion

makeblastdb -dbtype nucl -in /work/jawlab/data/bladderwort/genomes/utricularia/scaffolds.ugibba_lk.fasta -out /work/jawlab/data/bladderwort/genomes/utricularia/scaffolds.ugibba_lk.fasta

blastn -db $dataDir/genomes/utricularia/scaffolds.ugibba_lk.fasta -query $dataDir/genomes/utricularia/Utricularia_gibba_v2.chloroplast_NC_021449.1.fasta -out $assemblyDir/spades_assembly1/scaffolds.ugibba_lk.chloroHits.txt -perc_identity 90 -num_threads $PROCS -max_target_seqs 10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen"
blastn -db $dataDir/genomes/utricularia/scaffolds.ugibba_lk.fasta -query $dataDir/genomes/utricularia/Utricularia_reniformis.mitochondrion_NC_034982.1.fasta -out $assemblyDir/spades_assembly1/scaffolds.ugibba_lk.mitoHits.txt -perc_identity 90 -num_threads $PROCS -max_target_seqs 10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen"

cut -f2 $assemblyDir/spades_assembly1/scaffolds.ugibba_lk.chloroHits.txt | uniq > $assemblyDir/spades_assembly1/organellar_contigs.txt 
cut -f2 $assemblyDir/spades_assembly1/scaffolds.ugibba_lk.mitoHits.txt | uniq >> $assemblyDir/spades_assembly1/organellar_contigs.txt

# get only non-chloro/mito scaffolds, run quast compared to assembly with chloro/mito contigs
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' $dataDir/genomes/utricularia/scaffolds.ugibba_lk.fasta | grep -Ffv $assemblyDir/spades_assembly1/organellar_contigs.txt - | tr "\t" "\n" > $assemblyDir/spades_assembly1/scaffolds.ugibba.noOrganellar.fasta
quast.py -o $assemblyDir/spades_assembly1/quast_ugibba_output_rmChloroMito -r $dataDir/genomes/utricularia/Utricularia_gibba_v2.fa -g $dataDir/genomes/utricularia/u.gibba_NEW.genic.gff -t $PROCS --eukaryote -b --min-identity 90.0 $assemblyDir/spades_assembly1/scaffolds.ugibba.noOrganellar.fasta $dataDir/genomes/utricularia/scaffolds.ugibba_lk.fasta

myGenome=/work/jawlab/data/bladderwort/genomes/utricularia/scaffolds.ugibba_lk.fasta

######################
# Maker annotation on assembly > 1kb
######################

#first, trinity assemble transcripts

Trinity --seqType fq --JM 50G --left $rnaDir/tank1_forward_paired.fq.gz,$rnaDir/tank2_forward_paired.fq.gz,$rnaDir/tank3_forward_paired.fq.gz --right $rnaDir/tank1_reverse_paired.fq.gz,$rnaDir/tank2_reverse_paired.fq.gz,$rnaDir/tank3_reverse_paired.fq.gz --CPU $PROCS --output $assemblyDir/Trinity

# take only U.gibba classified contigs
kraken2 -db $KrakenDB --output $assemblyDir/Trinity/trinity.kraken.out --report $assemblyDir/Trinity/trinity.kraken.report --threads $PROCS $assemblyDir/Trinity/Trinity.fasta
 
awk '{print $1, $2, $3, $4}' $assemblyDir/Trinity/trinity.kraken.out | sed "s/ /\t/g" | awk '$3=="13748"' > $assemblyDir/Trinity/trinity.ugibba.txt
awk '{print $2}' $assemblyDir/Trinity/trinity.ugibba.txt > $assemblyDir/Trinity/trinity_names.ugibba.txt

seqtk subseq $assemblyDir/Trinity/Trinity.fasta $assemblyDir/Trinity/trinity_names.ugibba.txt > $assemblyDir/Trinity/Trinity.ugibba.fasta

# running maker NOTE: maker steps taken from https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2
# edited maker_opts.ctl to include trinity assembled transcripts and asterid proteins downloaded from ensembl PLAZA
makerDir=$assemblyDir/Maker_2ndTry

if [ ! -e $makerDir ]; then mkdir $makerDir; fi
cd $makerDir

maker -base ugibba_rnd1 maker_opts.ctl maker_bopts.ctl maker_exe.ctl

# ########
# # SNAP
# ########
snapDir=$makerDir/Snap
if [ ! -e $snapDir ]; then mkdir $snapDir; fi

mkdir $snapDir/roundOne
cd $snapDir/roundOne

# export 'confident' gene models from MAKER and rename to something meaningful
maker2zff -x 0.25 -l 50 -d $makerDir/ugibba_rnd1.maker.output/ugibba_rnd1_master_datastore_index.log
rename 's/genome/ugibba_rnd1.zff.length50_aed0.25/g' *
# # gather some stats and validate
fathom ugibba_rnd1.zff.length50_aed0.25.ann ugibba_rnd1.zff.length50_aed0.25.dna -gene-stats > gene-stats.log 2>&1
fathom ugibba_rnd1.zff.length50_aed0.25.ann ugibba_rnd1.zff.length50_aed0.25.dna -validate > validate.log 2>&1
# # collect the training sequences and annotations, plus 1000 surrounding bp for training
fathom ugibba_rnd1.zff.length50_aed0.25.ann ugibba_rnd1.zff.length50_aed0.25.dna -categorize 1000 > categorize.log 2>&1
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1
# # create the training parameters
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1
cd ..
# # assemble the HMM
hmm-assembler.pl ugibba_rnd1.zff.length50_aed0.25 params > ugibba_rnd1.zff.length50_aed0.25.hmm

###########
# AUGUSTUS - Round1
###########

gff3_merge -s -d $makerDir/ugibba_rnd1.maker.output/ugibba_rnd1_master_datastore_index.log > $makerDir/ugibba_rnd1.all.maker.gff
fasta_merge -d $makerDir/ugibba_rnd1.maker.output/Bcon_rnd1_master_datastore_index.log
# # GFF w/o the sequences
gff3_merge -n -s -d $makerDir/ugibba_rnd1.maker.output/ugibba_rnd1_master_datastore_index.log > $makerDir/ugibba_rnd1.all.maker.noseq.gff
# 
awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' $makerDir/ugibba_rnd1.all.maker.noseq.gff | \
  awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+500; else print $1, $2-500, $3+500 }' | \
  bedtools getfasta -fi $myGenome -bed - -fo $makerDir/ugibba_rnd1.all.maker.transcripts500.fasta

cd $makerDir
export AUGUSTUS_CONFIG_PATH=~/Augustus/config
python /usr/local/bin/busco/scripts/run_BUSCO.py -f -i ugibba_rnd1.all.maker.transcripts500.fasta -o ugibba_rnd1.maker.output -l ../eudicotyledons_odb10/ \
    -m geno -sp tomato -c 8 --long -z --augustus_parameters='--progress=true'

# performed busco test on included busco data with installation. interestingly, sample data results using busco on cluster did not match expected output. must be something wrong with installation? NOTE: the installation problem did not get resolved so I needed to install locally

rename "s/BUSCO_ugibba_rnd1.maker.output_3427324986/Utricularia_gibba/g" *
sed -i 's/BUSCO_ugibba_rnd1.maker.output_3427324986/Utricularia_gibba/g' Utricularia_gibba_parameters.cfg
sed -i 's/BUSCO_ugibba_rnd1.maker.output_3427324986/Utricularia_gibba/g' Utricularia_gibba_parameters.cfg.orig1

mkdir ~/Xfer-home/Augustus/config/species/Utricularia_gibba
cp Utricularia_gibba_* ~/Xfer-home/Augustus/config/species/Utricularia_gibba/

awk '{ if ($2 == "est2genome") print $0 }' ugibba_rnd1.all.maker.noseq.gff > ugibba_rnd1.maker.output/ugibba_rnd1.all.maker.est2genome.gff
# # protein alignments
awk '{ if ($2 == "protein2genome") print $0 }' ugibba_rnd1.all.maker.noseq.gff > ugibba_rnd1.maker.output/ugibba_rnd1.all.maker.protein2genome.gff
# # repeat alignments
awk '{ if ($2 ~ "repeat") print $0 }' ugibba_rnd1.all.maker.noseq.gff > ugibba_rnd1.maker.output/ugibba_rnd1.all.maker.repeats.gff

# source activate /scratch/lk82153/jwlab/maker-env

makerDirLscratch=/lscratch/"$PBS_JOBID"/maker_run
mkdir $makerDirLscratch
cd $makerDirLscratch

# #need to set RepeatMasker environment variables
export AUGUSTUS_CONFIG_PATH=~/Augustus/config
export REPEATMASKER_LIB_DIR=/scratch/lk82153/jwlab/maker-env/share/RepeatMasker/Libraries
export REPEATMASKER_MATRICES_DIR=/scratch/lk82153/jwlab/maker-env/share/RepeatMasker/Matrices

mpiexec -n 40 maker -base ugibba_rnd2 $makerDir/round2_maker_opts.ctl $makerDir/maker_bopts.ctl $makerDir/maker_exe.ctl

conda deactivate

mv $makerDirLscratch/ugibba_rnd2.maker.output $makerDir

###############
# Snap/Augustus - round2 + maker round 3
###############

gff3_merge -s -d $makerDir/ugibba_rnd2.maker.output/ugibba_rnd2_master_datastore_index.log > $makerDir/ugibba_rnd2.all.maker.gff
fasta_merge -d $makerDir/ugibba_rnd2.maker.output/ugibba_rnd2_master_datastore_index.log 
# # GFF w/o the sequences
gff3_merge -n -s -d $makerDir/ugibba_rnd2.maker.output/ugibba_rnd2_master_datastore_index.log > $makerDir/ugibba_rnd2.all.maker.noseq.gff

mkdir $makerDir/Snap/roundTwo
cd $makerDir/Snap/roundTwo
# # export 'confident' gene models from MAKER and rename to something meaningful
maker2zff -x 0.25 -l 50 -d ../../ugibba_rnd2.maker.output/ugibba_rnd2_master_datastore_index.log
rename 's/genome/ugibba_rnd2.zff.length50_aed0.25/g' *
# # gather some stats and validate
fathom ugibba_rnd2.zff.length50_aed0.25.ann ugibba_rnd2.zff.length50_aed0.25.dna -gene-stats > gene-stats.log 2>&1
fathom ugibba_rnd2.zff.length50_aed0.25.ann ugibba_rnd2.zff.length50_aed0.25.dna -validate > validate.log 2>&1
# # collect the training sequences and annotations, plus 1000 surrounding bp for training
fathom ugibba_rnd2.zff.length50_aed0.25.ann ugibba_rnd2.zff.length50_aed0.25.dna -categorize 1000 > categorize.log 2>&1
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1
# # create the training parameters
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1
cd ..
# # assembly the HMM
hmm-assembler.pl ugibba_rnd2.zff.length50_aed0.25 params > ugibba_rnd2.zff.length50_aed0.25.hmm

awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' ugibba_rnd2.all.maker.noseq.gff | \
  awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+500; else print $1, $2-500, $3+500 }' | \
  bedtools getfasta -fi $myGenome -bed - -fo ugibba_rnd2.all.maker.transcripts500.fasta

cd $makerDir
python2.7 /usr/local/bin/busco/scripts/run_BUSCO.py -f -i ugibba_rnd2.all.maker.transcripts500.fasta -o ugibba_rnd2.maker.output -l ../eudicotyledons_odb10/ -m geno -sp tomato -c $PROCS --long -z --augustus_parameters='--progress=true'

# performed busco test on included busco data with installation. interestingly, sample data results using busco on cluster did not match expected output. must be something wrong with installation? NOTE: the installation problem did not get resolved so I needed to install locally

rename 's/BUSCO_ugibba_rnd2.maker.output_3001126393/Utricularia_gibba_2nd_training/g' *
sed -i 's/BUSCO_ugibba_rnd2.maker.output_3001126393/Utricularia_gibba_2nd_training/g' Utricularia_gibba_2nd_training_parameters.cfg
sed -i 's/BUSCO_ugibba_rnd2.maker.output_3001126393/Utricularia_gibba_2nd_training/g' Utricularia_gibba_2nd_training_parameters.cfg.orig1

mkdir $AUGUSTUS_CONFIG_PATH/species/Utricularia_gibba_2nd_training
cp Utricularia_gibba_* ~/Home-xfer/Augustus/config/species/Utricularia_gibba_2nd_training/

cd $makerDir

# using transcript/protein/repeat gffs from initial round of maker

source activate /scratch/lk82153/jwlab/maker-env

makerDirLscratch=/lscratch/"$PBS_JOBID"/maker_run
mkdir $makerDirLscratch
cd $makerDirLscratch

# #need to set RepeatMasker environment variables
export AUGUSTUS_CONFIG_PATH=~/Augustus/config
export REPEATMASKER_LIB_DIR=/scratch/lk82153/jwlab/maker-env/share/RepeatMasker/Libraries
export REPEATMASKER_MATRICES_DIR=/scratch/lk82153/jwlab/maker-env/share/RepeatMasker/Matrices

mpiexec -n 40 maker -base ugibba_rnd3 $makerDir/round3_maker_opts.ctl $makerDir/maker_bopts.ctl $makerDir/maker_exe.ctl

conda deactivate

mv $makerDirLscratch/ugibba_rnd3.maker.output $makerDir

gff3_merge -s -d $makerDir/ugibba_rnd3.maker.output/ugibba_rnd3_master_datastore_index.log > $makerDir/ugibba_rnd3.all.maker.gff
fasta_merge -d $makerDir/ugibba_rnd3.maker.output/ugibba_rnd3_master_datastore_index.log
# # GFF w/o the sequences
gff3_merge -n -s -d $makerDir/ugibba_rnd3.maker.output/ugibba_rnd3_master_datastore_index.log > $makerDir/ugibba_rnd3.all.maker.noseq.gff

###############
# Snap/Augustus - round3 + maker round 4
###############

mkdir Snap/roundThree
cd Snap/roundThree
# # export 'confident' gene models from MAKER and rename to something meaningful
maker2zff -x 0.25 -l 50 -d ../../ugibba_rnd3.maker.output/ugibba_rnd3_master_datastore_index.log
mv genome.ann ugibba_rnd3.zff.length50_aed0.25.ann
mv genome.dna ugibba_rnd3.zff.length50_aed0.25.dna
# # gather some stats and validate
fathom ugibba_rnd3.zff.length50_aed0.25.ann ugibba_rnd3.zff.length50_aed0.25.dna -gene-stats > gene-stats.log 2>&1
fathom ugibba_rnd3.zff.length50_aed0.25.ann ugibba_rnd3.zff.length50_aed0.25.dna -validate > validate.log 2>&1
# # collect the training sequences and annotations, plus 1000 surrounding bp for training
fathom ugibba_rnd3.zff.length50_aed0.25.ann ugibba_rnd3.zff.length50_aed0.25.dna -categorize 1000 > categorize.log 2>&1
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1
# # create the training parameters
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1
cd ..
# # assemble the HMM
hmm-assembler.pl ugibba_rnd3.zff.length50_aed0.25 params > ugibba_rnd3.zff.length50_aed0.25.hmm

cd $makerDir

awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' ugibba_rnd3.all.maker.noseq.gff | \
  awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+500; else print $1, $2-500, $3+500 }' | \
  bedtools getfasta -fi ../spades_assembly1/scaffolds.ugibba.fasta -bed - -fo ugibba_rnd3.all.maker.transcripts500.fasta

cd $makerDir
python2.7 run_BUSCO.py -f -i ugibba_rnd3.all.maker.transcripts500.fasta -o ugibba_rnd3.maker.output -l eudicotyledons_odb10/ -m geno -sp tomato -c $PROCS --long -z --augustus_parameters='--progress=true'

performed busco test on included busco data with installation. interestingly, sample data results using busco on cluster did not match expected output. must be something wrong with installation? NOTE: the installation problem did not get resolved so I needed to install locally

rename -i 's/BUSCO_ugibba_rnd1.maker.output_2318503959/Utricularia_gibba/g' Utricularia_gibba_parameters.cfg
sed -i 's/BUSCO_ugibba_rnd1.maker.output_2318503959/Utricularia_gibba/g' Utricularia_gibba_parameters.cfg
sed -i 's/BUSCO_ugibba_rnd1.maker.output_2318503959/Utricularia_gibba/g' Utricularia_gibba_parameters.cfg.orig1

mkdir $AUGUSTUS_CONFIG_PATH/species/Utricularia_gibba_2nd_training
cp Utricularia_gibba_* ~/Home-xfer/Augustus/config/species/Utricularia_gibba_2nd_training/

cd $makerDir

export AUGUSTUS_CONFIG_PATH=~/Augustus/config
maker -base ugibba_rnd3 round3_maker_opts.ctl maker_bopts.ctl maker_exe.ctl

gff3_merge -s -d $makerDir/ugibba_rnd3.maker.output/ugibba_rnd3_master_datastore_index.log > $makerDir/ugibba_rnd3.all.maker.gff
fasta_merge -d $makerDir/ugibba_rnd3.maker.output/ugibba_rnd3_master_datastore_index.log
# # GFF w/o the sequences
gff3_merge -n -s -d $makerDir/ugibba_rnd3.maker.output/ugibba_rnd3_master_datastore_index.log > $makerDir/ugibba_rnd3.all.maker.noseq.gff
export AUGUSTUS_CONFIG_PATH=~/Augustus/config
maker -base ugibba_rnd4 round4_maker_opts.ctl maker_bopts.ctl maker_exe.ctl

gff3_merge -s -d $makerDir/ugibba_rnd4.maker.output/ugibba_rnd4_master_datastore_index.log > $makerDir/ugibba_rnd4.all.maker.gff
fasta_merge -d $makerDir/ugibba_rnd4.maker.output/ugibba_rnd4_master_datastore_index.log
# # GFF w/o the sequences
gff3_merge -n -s -d $makerDir/ugibba_rnd4.maker.output/ugibba_rnd4_master_datastore_index.log > $makerDir/ugibba_rnd4.all.maker.noseq.gff


#########################
# Maker results summary
#########################

#get number of genes and average length for each assembly - > 5000

# # >5000
cat ugibba_rnd1.all.maker.noseq.gff | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 -$4) } END { print NR, sum / NR }'
# #22356 - 1580.57
cat ugibba_rnd2.all.maker.noseq.gff | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 -$4) } END { print NR, sum / NR }'
# #15528 - 2283.86
cat ugibba_rnd3.all.maker.noseq.gff | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 -$4) } END { print NR, sum / NR }'
# #16958 2181.74
cat ugibba_rnd4.all.maker.noseq.gff | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 -$4) } END { print NR, sum / NR }'
# # 16909 2188.82

# # > 1000
cat ugibba_rnd1.all.maker.noseq.gff | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 -$4) } END { print NR, sum / NR }'
# # 22209 1595.33
cat ugibba_rnd2.all.maker.noseq.gff | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 -$4) } END { print NR, sum / NR }'
# # 19361 2110.47
cat ugibba_rnd3.all.maker.noseq.gff | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 -$4) } END { print NR, sum / NR }'
# # 19203 2139.83

# # maker AED distribution for >1000
perl ../Maker/AED_cdf_generator.pl -b 0.025 ugibba_rnd1.all.maker.gff > rnd1.AED.txt
perl ../Maker/AED_cdf_generator.pl -b 0.025 ugibba_rnd2.all.maker.gff > rnd2.AED.txt
perl ../Maker/AED_cdf_generator.pl -b 0.025 ugibba_rnd3.all.maker.gff > rnd3.AED.txt

# #fix names of annotations
maker_map_ids --prefix ugibba --justify 5  $makerDir/ugibba_rnd3.all.maker.gff > $makerDir/ugibba_rnd3.all.maker.name.map

map_gff_ids $makerDir/ugibba_rnd3.all.maker.name.map $makerDir/ugibba_rnd3.all.maker.gff
map_gff_ids $makerDir/ugibba_rnd3.all.maker.name.map $makerDir/ugibba_rnd3.all.maker.noseq.gff

map_fasta_ids $makerDir/ugibba_rnd3.all.maker.name.map $makerDir/ugibba_rnd3.all.maker.transcripts.fasta
map_fasta_ids $makerDir/ugibba_rnd3.all.maker.name.map $makerDir/ugibba_rnd3.all.maker.proteins.fasta

# #move annotations to genome folder
cp $makerDir/ugibba_rnd3.all.maker.noseq.gff $dataDir/genomes/utricularia/scaffolds.ugibba_lk.gff


# get genePairs file for downstream analysis
Rscript $scriptDir/assignConvergentDivergentParallel.R -i $dataDir/genomes/utricularia/scaffolds.ugibba_lk.gff -o $dataDir/genomes/utricularia/scaffolds.ugibba_lk.genePairs.txt

########################################
# Other Stuff
########################################

# # compare assemblies - blast
module load BLAST+
makeblastdb -in /work/jawlab/data/bladderwort/genomes/utricularia/Utricularia_gibba_v2.fa -dbtype nucl -out /work/jawlab/data/bladderwort/genomes/utricularia/Utricularia_gibba_v2
blastn -db /work/jawlab/data/bladderwort/genomes/utricularia/Utricularia_gibba_v2 -query /scratch/lk82153/jwlab/Bladderwort/8_GenomeAssembly/spades_assembly1/scaffolds.ugibba.fasta -out /scratch/lk82153/jwlab/Bladderwort/8_GenomeAssembly/spades_assembly1/scaffolds.ugibba_to_pacbio.blastout -perc_identity 90 -num_threads 8 -max_target_seqs 10 -outfmt 6

seqtk subseq /scratch/lk82153/jwlab/Bladderwort/8_GenomeAssembly/spades_assembly1/contigs.fasta /scratch/lk82153/jwlab/Bladderwort/8_GenomeAssembly/spades_assembly1/scaffolds.ugibba.unitig_0.txt > /scratch/lk82153/jwlab/Bladderwort/8_GenomeAssembly/spades_assembly1/scaffolds.ugibba.unitig_0.fasta

# extract genic regions
bedtools getfasta -s -fi $dataDir/genomes/utricularia/Utricularia_gibba_v2.fa -bed $dataDir/genomes/utricularia/u.gibba_NEW.genes.bed -fo $dataDir/genomes/utricularia/u.gibba_NEW.genes.fa

# # compare annotations with blast - pacbio annos mapped to my genome
blastn -db $myGenome -query $dataDir/genomes/utricularia/u.gibba_NEW.genes.fa -out /scratch/lk82153/jwlab/Bladderwort_3pri/1_Assembly/genes_PacBio_to_scaffolds_lk.blastout -perc_identity 90 -num_threads 8 -outfmt 6

# make bed file of gene 
# # compare annotations with blast - my annotations mapped to pacbio genome
Rscript $scriptDir/gffToBed.R -i $dataDir/genomes/utricularia/scaffolds.ugibba_lk.ALL.includingPacBio.gff -o $dataDir/genomes/utricularia/scaffolds.ugibba_lk.ALL.matches
Rscript $scriptDir/blastToBed.R -i /scratch/lk82153/jwlab/Bladderwort_3pri/1_Assembly/genes_PacBio_to_scaffolds_lk.blastout -o /scratch/lk82153/jwlab/Bladderwort_3pri/1_Assembly/pacbio_anno_coordinates.bed

Rscript $scriptDir/bedtogff.R -i /scratch/lk82153/jwlab/Bladderwort_3pri/1_Assembly/pacbio_anno_coordinates.bed -o /scratch/lk82153/jwlab/Bladderwort_3pri/1_Assembly/pacbio_anno_coordinates.gff

bedtools getfasta -s -name -fi $dataDir/genomes/utricularia/scaffolds.ugibba_lk.fasta -bed $dataDir/genomes/utricularia/scaffolds.ugibba_lk.ALL.includingPacBio.bed -fo $dataDir/genomes/utricularia/scaffolds.ugibba_lk.genes.fa

blastn -db $dataDir/genomes/utricularia/Utricularia_gibba_v2 -query $dataDir/genomes/utricularia/scaffolds.ugibba_lk.genes.fa -out /scratch/lk82153/jwlab/Bladderwort_3pri/1_Assembly/genes_scaffolds_lk_to_PacBio.blastout -perc_identity 90 -num_threads 8 -outfmt 6

# make bed file of gene matches
Rscript $scriptDir/blastToBed.R -i /scratch/lk82153/jwlab/Bladderwort_3pri/1_Assembly/genes_scaffolds_lk_to_PacBio.blastout -o /scratch/lk82153/jwlab/Bladderwort_3pri/1_Assembly/scaffolds_lk_anno_coordinates_in_PacBio.bed

Rscript $scriptDir/bedToGff.R -i /scratch/lk82153/jwlab/Bladderwort_3pri/1_Assembly/scaffolds_lk_anno_coordinates_in_PacBio.bed -o /scratch/lk82153/jwlab/Bladderwort_3pri/1_Assembly/scaffolds_lk_anno_coordinates_in_PacBio.gff -s BLAST -f gene -p scaffolds_lk
