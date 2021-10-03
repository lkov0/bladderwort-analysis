#4_findMotifs.sh

parentDir=$1
scriptDir=$2
dataDir=$3
analysisDir=$4
memeDir=$5

NPROCS=$6

module load R

intersectDir=$memeDir/Intersect

if [ ! -e $intersectDir ]; then mkdir $intersectDir; fi

# find intersection between conserved regions and CRE regions
Rscript $scriptDir/makeIntergenicBed.R -i $analysisDir/uni_cand.3pri.pbioAssembly.tsv -o $memeDir/uni_cand
Rscript $scriptDir/makeIntergenicBed.R -i $analysisDir/bid_cand.3pri.pbioAssembly.tsv -o $memeDir/bid_cand
Rscript $scriptDir/makeIntergenicBed.R -i $analysisDir/term_cand.parallel.3pri.pbioAssembly.tsv -o $memeDir/term_cand.parallel
Rscript $scriptDir/makeIntergenicBed.R -i $analysisDir/term_cand.convergent.3pri.pbioAssembly.tsv -o $memeDir/term_cand.convergent
Rscript $scriptDir/makeIntergenicBed.R -i $analysisDir/ins_cand.convergent.3pri.pbioAssembly.tsv -o $memeDir/ins_cand.convergent
Rscript $scriptDir/makeIntergenicBed.R -i $analysisDir/ins_cand.divergent.3pri.pbioAssembly.tsv -o $memeDir/ins_cand.divergent

module unload R
module load BEDTools

for cand in uni_cand bid_cand term_cand.convergent term_cand.parallel ins_cand.convergent ins_cand.divergent; do
    bedtools sort -i $memeDir/$cand.bed > $memeDir/$cand.sorted.bed
    bedtools intersect -a $memeDir/$cand.sorted.bed -b $memeDir/asterids.most-conserved.final.bed > $intersectDir/$cand.asterids.intersect.txt
    bedtools intersect -a $memeDir/$cand.sorted.bed -b $memeDir/rosids.most-conserved.final.bed > $intersectDir/$cand.rosids.intersect.txt
    bedtools intersect -a $memeDir/$cand.sorted.bed -b $memeDir/monocots.most-conserved.final.bed > $intersectDir/$cand.monocots.intersect.txt
done

# also get all intergenic regions - my genome
module load R
Rscript $scriptDir/makeIntergenicBed.R -i $dataDir/genomes/utricularia/scaffolds.ugibba_lk.ALL.includingPacBio.genePairs.txt -o $memeDir/all_pairs
module unload R
module load BEDTools
bedtools getfasta -name -fi $dataDir/genomes/utricularia/scaffolds.ugibba_lk.fasta -bed $memeDir/all_pairs.bed -fo $memeDir/all_pairs.fasta

# also get all intergenic regions -pbio genome
module load R
Rscript $scriptDir/makeIntergenicBed.R -i $dataDir/genomes/utricularia/u.gibba_NEW.genePairs.txt -o $memeDir/PacBio_Genome/all_pairs_pacbio
module unload R
module load BEDTools
bedtools getfasta -name -fi $dataDir/genomes/utricularia/Utricularia_gibba_v2.fa -bed $memeDir/PacBio_Genome/all_pairs_pacbio.bed -fo $memeDir/PacBio_Genome/all_pairs_pacbio.fasta

# intersect hits from pacbio conserved regions
bedtools sort -i $memeDir/PacBio_Genome/all_pairs_pacbio.bed > $memeDir/PacBio_Genome/all_pairs_pacbio.sorted.bed
bedtools intersect -a $memeDir/PacBio_Genome/all_pairs_pacbio.sorted.bed -b $memeDir/PacBio_Genome/asterids.most-conserved.final.bed > $intersectDir/all_pairs.asterids.intersect.txt
bedtools intersect -a $memeDir/PacBio_Genome/all_pairs_pacbio.sorted.bed -b $memeDir/PacBio_Genome/rosids.most-conserved.final.bed > $intersectDir/all_pairs.rosids.intersect.txt
bedtools intersect -a $memeDir/PacBio_Genome/all_pairs_pacbio.sorted.bed -b $memeDir/PacBio_Genome/monocots.most-conserved.final.bed > $intersectDir/all_pairs.monocots.intersect.txt

################
# motif finding 
################

#get top 40 candidates in each category, get fasta file of sequences
for file in $(ls $analysisDir/*summary.tsv | sed "s/.tsv//g"); do
    head $file.tsv -n 41 > $file.top40.tsv
    module load R
    Rscript $scriptDir/makeIntergenicBed.R -i $file.top20.tsv -o $memeDir/$file
    module unload R
    module load BEDTools
    bedtools getfasta -name -fi $dataDir/genomes/utricularia/scaffolds.ugibba_lk.fasta -bed $memeDir/$file.bed -fo $memeDir/$file.fasta
done

#use meme on my computer at work.....

# for type in bidirectionalPromoter insulator terminator unidirectionalPromoter; do
#     meme $memeDir/${type}_candidates.fasta -oc $memeDir/meme_${type} -dna -p $NPROCS -nmotifs 10;
# done
# 
# motifFile=$memeDir/JASPAR2018_CORE_plants_non-redundant_pfms_meme.txt


# for type in bidirectionalPromoter insulator terminator unidirectionalPromoter; do
#     fimo --oc $memeDir/fimo_jaspar_${type} $motifFile $memeDir/${type}_candidates.fasta;
# done
# 
# 
# promoterFile=$memeDir/Transcription_factor_weight_matrix_plantPAN.edited.memeFormat.txt
# for type in bidirectionalPromoter insulator terminator unidirectionalPromoter; do
#     fimo --oc $memeDir/fimo_promoter_${type} $promoterFile $memeDir/${type}_candidates.fasta;
# done
# 
# plottingDir=$parentDir/7_Plotting
# 
# if [ ! -e $plottingDir ]; then mkdir $plottingDir; done
# 
# #get genome file from index
# 
# for file in $(ls $alignmentDir/*multimapTroubleshooting*bam | sed "s/.bam//g"); do
# #     bedtools bamtobed -i $file.bam > $file.bed
#     bedtools genomecov -bga -i $file.bed -g /work/jawlab/data/bladderwort/genomes/utricularia/Utricularia_gibba_v2.genome > $file.bedgraph
# done
# 
# #perform union on samples of same tissue bedgraph
# bedtools unionbedg -i 1L_multimapTroubleshootingAligned.sortedByCoord.out.bedgraph 2L_multimapTroubleshootingAligned.sortedByCoord.out.bedgraph > $plottingDir/all_leaf.bedgraph
# bedtools unionbedg -i 1B_multimapTroubleshootingAligned.sortedByCoord.out.bedgraph 2B_multimapTroubleshootingAligned.sortedByCoord.out.bedgraph 3B_multimapTroubleshootingAligned.sortedByCoord.out.bedgraph > $plottingDir/all_bladder.bedgraph
# bedtools unionbedg -i 1S_multimapTroubleshootingAligned.sortedByCoord.out.bedgraph 2S_multimapTroubleshootingAligned.sortedByCoord.out.bedgraph 3S_multimapTroubleshootingAligned.sortedByCoord.out.bedgraph > $plottingDir/all_stem.bedgraph
# bedtools unionbedg -i 1R_multimapTroubleshootingAligned.sortedByCoord.out.bedgraph 2R_multimapTroubleshootingAligned.sortedByCoord.out.bedgraph 3R_multimapTroubleshootingAligned.sortedByCoord.out.bedgraph > $plottingDir/all_rhizoid.bedgraph

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

#############
# CEMiTool co-expression based promoter finding (taking 200bp before the annotated start position of the gene)
# two modules were found
#############

# for analyses done to get the cemitool.bed files, see runCemiTools.R

cemiDir=$memeDir/CEMiTools
if [ ! -e $cemiDir ]; then mkdir $cemiDir; fi
mv $memeDir/*cemitool.bed $cemiDir

motifFile=$memeDir/PacBio_Genome/JASPAR2018_CORE_plants_non-redundant_pfms_meme.txt
promoterFile=$memeDir/PacBio_Genome/Transcription_factor_weight_matrix_plantPAN.edited.memeFormat.txt
intergenicFasta=$dataDir/


module load BEDTools
module load MEME

for i in {1..10}; do 
    
#     bedtools getfasta -s -name -fi $dataDir/genomes/utricularia/scaffolds.ugibba_lk.fasta -bed $cemiDir/module$i.cemitool.bed -fo $cemiDir/module$i.promoterRegions.cemitool.fasta # NOTE: length zero promoter regions (TSS started on edge of contig) and promoter regions extending beyond the length of the contig were skipped
    
#     meme $cemiDir/module$i.promoterRegions.cemitool.fasta -oc $cemiDir/MEMEout.module$i -dna -p 1 -nmotifs 10;

## use TOMTOM to find matches to JASPAR plant database in our cemitool module motifs
    tomtom -thresh 0.005 -o $cemiDir/MEMEout.module$i/TOMTOMout.jaspar.module$i $cemiDir/MEMEout.module$i/meme.txt $motifFile

    tomtom -thresh 0.005 -o $cemiDir/MEMEout.module$i/TOMTOMout.promoters.module$i $cemiDir/MEMEout.module$i/meme.txt $promoterFile
    
done

# take only significant motifs (no significant motifs in module 9 or 10)
meme-get-motif -ia -id MEME-1 -id MEME-2 -id MEME-3 -id MEME-4 -id MEME-5 -id MEME-6 -id MEME-7 -id MEME-8 -id MEME-9 $cemiDir/MEMEout.module1/meme.txt > $cemiDir/MEMEout.module1/meme.significant.txt

meme-get-motif -ia -id MEME-1 -id MEME-2 -id MEME-3 -id MEME-4 -id MEME-5 -id MEME-8 $cemiDir/MEMEout.module2/meme.txt > $cemiDir/MEMEout.module2/meme.significant.txt

meme-get-motif -ia -id MEME-1 -id MEME-2 -id MEME-3 $cemiDir/MEMEout.module3/meme.txt > $cemiDir/MEMEout.module3/meme.significant.txt

meme-get-motif -ia -id MEME-1 -id MEME-2 $cemiDir/MEMEout.module4/meme.txt > $cemiDir/MEMEout.module4/meme.significant.txt

meme-get-motif -ia -id MEME-1 -id MEME-2 $cemiDir/MEMEout.module5/meme.txt > $cemiDir/MEMEout.module5/meme.significant.txt

meme-get-motif -ia -id MEME-1 -id MEME-2 $cemiDir/MEMEout.module6/meme.txt > $cemiDir/MEMEout.module6/meme.significant.txt

meme-get-motif -ia -id MEME-1 $cemiDir/MEMEout.module7/meme.txt > $cemiDir/MEMEout.module7/meme.significant.txt

meme-get-motif -ia -id MEME-1$cemiDir/MEMEout.module8/meme.txt > $cemiDir/MEMEout.module8/meme.significant.txt

for i in {1..8}; do 

    fimo --o $cemiDir/MEMEout.module$i/FIMO.promoters.module$i $cemiDir/MEMEout.module$i/meme.significant.txt $memeDir/PacBio_Genome/all_pairs_pacbio.fasta

done

#aggregate promoters based on significance in meme and homology to Plant PAN in tomtom

meme-get-motif -ia -id MEME-1 -id MEME-2 -id MEME-3 -id MEME-5 -id MEME-9 $cemiDir/MEMEout.module1/meme.txt > $cemiDir/MEMEout.module1/meme.promoters.txt

meme-get-motif -ia -id MEME-1 -id MEME-2 -id MEME-3 -id MEME-4 $cemiDir/MEMEout.module2/meme.txt > $cemiDir/MEMEout.module2/meme.promoters.txt

meme-get-motif -ia -id MEME-1 $cemiDir/MEMEout.module3/meme.txt > $cemiDir/MEMEout.module3/meme.promoters.txt

meme-get-motif -ia -id MEME-1 -id MEME-2 $cemiDir/MEMEout.module4/meme.txt > $cemiDir/MEMEout.module4/meme.promoters.txt

meme-get-motif -ia -id MEME-1 $cemiDir/MEMEout.module6/meme.txt > $cemiDir/MEMEout.module6/meme.promoters.txt

meme-get-motif -ia -id MEME-1 $cemiDir/MEMEout.module7/meme.txt > $cemiDir/MEMEout.module7/meme.promoters.txt

# combine promoters

meme2meme $cemiDir/MEMEout.module1/meme.promoters.txt \
    $cemiDir/MEMEout.module2/meme.promoters.txt \
    $cemiDir/MEMEout.module3/meme.promoters.txt \
    $cemiDir/MEMEout.module4/meme.promoters.txt \
    $cemiDir/MEMEout.module6/meme.promoters.txt \
    $cemiDir/MEMEout.module7/meme.promoters.txt > $cemiDir/cemiPromoters.txt

# find similarity to putative cre motifs (from our genome)
    
tomtom -thresh 0.005 -o $memeDir/meme.bid_cand/TOMTOM.cemiPromoters $memeDir/meme.bid_cand/meme.txt $cemiDir/cemiPromoters.txt

tomtom -thresh 0.005 -o $memeDir/meme.ins_cand/TOMTOM.cemiPromoters $memeDir/meme.ins_cand/meme.txt $cemiDir/cemiPromoters.txt

tomtom -thresh 0.005 -o $memeDir/meme.term_cand/TOMTOM.cemiPromoters $memeDir/meme.term_cand/meme.txt $cemiDir/cemiPromoters.txt

tomtom -thresh 0.005 -o $memeDir/meme.uni_cand/TOMTOM.cemiPromoters $memeDir/meme.uni_cand/meme.txt $cemiDir/cemiPromoters.txt

#########################
# intergenic regions
# find top 10 motifs in each intergenic region class
#########################

# extract significant motifs with meme get motif

for cand in uni_cand bid_cand term_cand.convergent term_cand.parallel ins_cand.convergent ins_cand.divergent; do
    # get fasta seqs
#     bedtools getfasta -name -fi $dataDir/genomes/utricularia/Utricularia_gibba_v2.fa -bed $memeDir/PacBio_Genome/$cand.bed -fo $memeDir/PacBio_Genome/$cand.fasta
    
    meme $memeDir/PacBio_Genome/$cand.fasta -oc $memeDir/PacBio_Genome/MEMEout.$cand -dna -p 1 -nmotifs 10
done

#extract significant motifs
meme-get-motif -ia -id MEME-1 -id MEME-2 -id MEME-3 $memeDir/PacBio_Genome/MEMEout.bid_cand/meme.txt > $memeDir/PacBio_Genome/MEMEout.bid_cand/meme.significant.txt

# no significant ins_cand.convergent candidates

meme-get-motif -ia -id MEME-1 -id MEME-2 -id MEME-3 $memeDir/PacBio_Genome/MEMEout.ins_cand.divergent/meme.txt > $memeDir/PacBio_Genome/MEMEout.ins_cand.divergent/meme.significant.txt

meme-get-motif -ia -id MEME-1 $memeDir/PacBio_Genome/MEMEout.term_cand.convergent/meme.txt > $memeDir/PacBio_Genome/MEMEout.term_cand.convergent/meme.significant.txt

meme-get-motif -ia -id MEME-1 -id MEME-2 -id MEME-3 -id MEME-4 -id MEME-5 $memeDir/PacBio_Genome/MEMEout.term_cand.parallel/meme.txt > $memeDir/PacBio_Genome/MEMEout.term_cand.parallel/meme.significant.txt

meme-get-motif -ia -id MEME-1 -id MEME-2 -id MEME-3 -id MEME-4 -id MEME-5 $memeDir/PacBio_Genome/MEMEout.uni_cand/meme.txt > $memeDir/PacBio_Genome/MEMEout.uni_cand/meme.significant.txt

# use TOMTOM to find matches to plant PAN and JASPAR databases

# fimo in pacbio intergenic regions
for cand in uni_cand bid_cand term_cand.convergent term_cand.parallel  ins_cand.divergent; do

    tomtom -thresh 0.005 -o $memeDir/PacBio_Genome/MEMEout.$cand/TOMTOMout.jaspar.module$i $memeDir/PacBio_Genome/MEMEout.$cand/meme.significant.txt $motifFile

    tomtom -thresh 0.005 -o $memeDir/PacBio_Genome/MEMEout.$cand/TOMTOMout.promoters.module$i $memeDir/PacBio_Genome/MEMEout.$cand/meme.significant.txt $promoterFile

    fimo --o $memeDir/PacBio_Genome/MEMEout.$cand/FIMO.intergenic.$cand $memeDir/PacBio_Genome/MEMEout.$cand/meme.significant.txt $memeDir/PacBio_Genome/all_pairs_pacbio.fasta

done

#find putative promoters in conserved intergenic sequences

for cand in uni_cand bid_cand term_cand.convergent term_cand.parallel ins_cand.convergent ins_cand.divergent; do

#     bedtools getfasta -name -fi $dataDir/genomes/utricularia/Utricularia_gibba_v2.fa -bed $memeDir/PacBio_Genome/$cand.conservedRegs.bed -fo $memeDir/PacBio_Genome/$cand.conservedRegs.fasta
    
    # find promoters from "our" set of promoters
#     fimo --oc $memeDir/PacBio_Genome/FIMOout.$cand.conservedRegs.promoters $cemiDir/cemiPromoters.txt $memeDir/PacBio_Genome/$cand.conservedRegs.fasta # 
    #look for proportion of sequences with matches to promoters
    #color which conserved sequences had hits to promoters in PlantPAN in R
    
    # find promoters from PlantPAN
    fimo --oc $memeDir/PacBio_Genome/FIMOout.$cand.conservedRegs.promoters $promoterFile $memeDir/PacBio_Genome/$cand.conservedRegs.fasta
done

