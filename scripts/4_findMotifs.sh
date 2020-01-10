#4_findMotifs.sh

parentDir=$1
scriptDir=$2
dataDir=$3
analysisDir=$4
memeDir=$5

NPROCS=$6

# module load R

################
# motif finding
################

#get top 40 candidates in each category, get fasta file of sequences
for file in $(ls $analysisDir/*summary.tsv | sed "s/.tsv//g"); do
    head $file.tsv -n 41 > $file.top40.tsv
    module load R
    Rscript $scriptDir/makeIntergenicBed.R -i $file.top40.tsv -o $memeDir/$file
    module unload R
    module load BEDTools
    bedtools getfasta -name -fi $dataDir/genomes/utricularia/Utricularia_gibba_v2.fa -bed $memeDir/$file.bed -fo $memeDir/$file.fasta
done


# module unload R

# module load BEDTools

# module unload BEDTools

#####meme
# module load OpenMPI
# module load MEME

#use meme on my computer at work.....

# for type in bidirectionalPromoter insulator terminator unidirectionalPromoter; do
#     meme $memeDir/${type}_candidates.fasta -oc $memeDir/meme_${type} -dna -p $NPROCS -nmotifs 10;
# done

motifFile=$memeDir/JASPAR2018_CORE_plants_non-redundant_pfms_meme.txt
for type in bidirectionalPromoter insulator terminator unidirectionalPromoter; do
    fimo --oc $memeDir/fimo_jaspar_${type} $motifFile $memeDir/${type}_candidates.fasta;
done


promoterFile=$memeDir/Transcription_factor_weight_matrix_plantPAN.edited.memeFormat.txt
for type in bidirectionalPromoter insulator terminator unidirectionalPromoter; do
    fimo --oc $memeDir/fimo_promoter_${type} $promoterFile $memeDir/${type}_candidates.fasta;
done

plottingDir=$parentDir/7_Plotting

if [ ! -e $plottingDir ]; then mkdir $plottingDir; done

#get genome file from index

for file in $(ls $alignmentDir/*multimapTroubleshooting*bam | sed "s/.bam//g"); do
#     bedtools bamtobed -i $file.bam > $file.bed
    bedtools genomecov -bga -i $file.bed -g /work/jawlab/data/bladderwort/genomes/utricularia/Utricularia_gibba_v2.genome > $file.bedgraph
done

#perform union on samples of same tissue bedgraph
bedtools unionbedg -i 1L_multimapTroubleshootingAligned.sortedByCoord.out.bedgraph 2L_multimapTroubleshootingAligned.sortedByCoord.out.bedgraph > $plottingDir/all_leaf.bedgraph
bedtools unionbedg -i 1B_multimapTroubleshootingAligned.sortedByCoord.out.bedgraph 2B_multimapTroubleshootingAligned.sortedByCoord.out.bedgraph 3B_multimapTroubleshootingAligned.sortedByCoord.out.bedgraph > $plottingDir/all_bladder.bedgraph
bedtools unionbedg -i 1S_multimapTroubleshootingAligned.sortedByCoord.out.bedgraph 2S_multimapTroubleshootingAligned.sortedByCoord.out.bedgraph 3S_multimapTroubleshootingAligned.sortedByCoord.out.bedgraph > $plottingDir/all_stem.bedgraph
bedtools unionbedg -i 1R_multimapTroubleshootingAligned.sortedByCoord.out.bedgraph 2R_multimapTroubleshootingAligned.sortedByCoord.out.bedgraph 3R_multimapTroubleshootingAligned.sortedByCoord.out.bedgraph > $plottingDir/all_rhizoid.bedgraph

####mummer

###finding conserved sequences in intergenic regions. using grape, mimulus, papaya, tomato, and arabidopsis
# module load BLAST+
# blastDir=$analysisDir/Blast
# if [ ! -e $blastDir ]; then mkdir $blastDir; fi
# 
# 
# makeblastdb -dbtype nucl -in $dataDir/New_Genome/u.gibba_NEW_candidateRegions.intergenic.parallel.fasta -out $dataDir/New_Genome/intergenic.parallel
# makeblastdb -dbtype nucl -in $dataDir/New_Genome/u.gibba_NEW_candidateRegions.intergenic.convergent.fasta -out $dataDir/New_Genome/intergenic.convergent
# makeblastdb -dbtype nucl -in $dataDir/New_Genome/u.gibba_NEW_candidateRegions.intergenic.divergent.fasta -out $dataDir/New_Genome/intergenic.divergent
# 
# for genome in $(ls $dataDir/Genomes | sed "s/.fna//g"); do
# blastn -db $dataDir/New_Genome/intergenic.parallel -query $dataDir/Genomes/$genome.fna -num_threads 8 -perc_identity 70 -num_alignments 1 -out $blastDir/$genome.parallel.blastOut.70.txt -outfmt "6 std qlen"
# blastn -db $dataDir/New_Genome/intergenic.convergent -query $dataDir/Genomes/$genome.fna -num_threads 8 -perc_identity 70 -num_alignments 1 -out $blastDir/$genome.convergent.blastOut.70.txt -outfmt "6 std qlen"
# blastn -db $dataDir/New_Genome/intergenic.divergent -query $dataDir/Genomes/$genome.fna -num_threads 8 -perc_identity 70 -num_alignments 1 -out $blastDir/$genome.divergent.blastOut.70.txt -outfmt "6 std qlen"
# done


###############
# comparing blast output
###############

# blastn -db $dataDir/New_Genome/u.gibba_NEW -query $dataDir/Old_Genome/u.gibba_OLD.genic.fa -num_threads 8 -perc_identity 95 -num_alignments 1 -out $analysisDir/u.gibba_OLD.u.gibba_NEW.blastOut.forAnalysis.txt -outfmt "6 std qlen"

# previously found that there was limited concordance of expression values between genomes... need to check this
# first generating list of all gene pair matches between genomes for dataset SRR768657 (More reads), also using GTF annotation generated FPKM values for old genome since these matched up with jason's FPKM values. 

# Rscript $scriptDir/getGenePairsSharedInGenomes.R -d $analysisDir/SRR768657_u.gibba_NEW.genePairs.foldChange.ALL.txt -q $analysisDir/SRR768657_u.gibba_OLD.genePairs.foldChange.ALL.test.gtf.txt -m $dataDir/geneMap_u.gibba_OLD_u.gibba_NEW.txt -o genomeMatched