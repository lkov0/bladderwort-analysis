#! /bin/bash
# original author - Jason Wallace. Modified for Lynsey's bladderwort analysis

# Make multiple genome alignment and do GERP analysis
# # # ROAST=$HOME/Software/Aligners/multiz-tba.012109/roast
# GERP=$HOME/Software/Aligners/gerp++/
# MULTIZ=$HOME/Software/Aligners/multiz-tba.012109/
# IGVTOOLS="bash /home/jgwall/Software/Browsers/IGVTools/igvtools"

#specify number of processors
# PROCS=48

# module load Biopython/1.68-foss-2016b-Python-3.5.2

# Set up directories
parentdir=~/Xfer/Bladderwort
scriptdir=~/Xfer/Repositories/bladderwort-analysis/scripts
# parentdir=/scratch/lk82153/jwlab/Bladderwort
# scriptdir=/scratch/lk82153/jwlab/Repositories/bladderwort-analysis/scripts
datadir=$parentdir/0_Data
# parentdir=$HOME/Work/Bladderwort/7_multiSpeciesAlignment
# datadir=$HOME/Work/Bladderwort/0_Data
genomedir=$datadir/Genomes/Asterids
sourcedir=$datadir/Genomes/Utricularia
aligndir=$parentdir/7b_PairwiseAlignment
combinedir=$aligndir/bMultiple_Alignment

if [ ! -e $aligndir ]; then mkdir $aligndir; fi
if [ ! -e $combinedir ]; then mkdir $combinedir; fi

# Alignment tree from phyloT (http://phylot.biobyte.de/)
# 13748
# 3625
# 4072
# 49390
# 72917
# 267555
# 192259
# 4232
# 4236
# 4155
# 4097
# 33119
# 212142
# 4081
# 69266
# 13750
TREE="(((VacciniumCorymbosum) ActinidiaChinensis) ((LactucaSativa (HelianthusAnnuus)) (CuscutaAustralis((CapsicumAnnuum SolanumLycopersicum) (PetuniaAxillaris)) (CoffeaCanephora (UtriculariaGibba MimulusGuttatus GenliseaAurea)))))"
# TREE2="((SolanumLycopersicum) (CoffeaCanephora (UtriculariaGibba MimulusGuttatus GenliseaAurea)))"

###########
# Pairwise genome alignments - Split apart for simplicity
###########

tname=UtriculariaGibba
target_genome=$sourcedir/UtriculariaGibba.GCA_002189035.1_U_gibba_v2_genomic.fna
joined_genome=$sourcedir/UtriculariaGibba.GCA_002189035.1_U_gibba_v2_genomic.all_in_one.fna
target_size=1500000000 # 1GB = bigger than any U.gibba chromosome
query_size=1000000  # 1 MB segments for each query

# #Join the U.gibba genome into a single scaffold for convenience
# python3 $scriptdir/1_JoinFastaWithPad.py -i $sourcedir/UtriculariaGibba.GCA_002189035.1_U_gibba_v2_genomic.fna -o $joined_genome \
#   -k $sourcedir/UtriculariaGibba.GCA_002189035.1_U_gibba_v2_genomic.all_in_one.key.txt -n 1000 --newname ${tname}_joined

# for query_genome in $(ls $genomedir/); do
#     qname=$(sed "s/\..*//g" <<<$query_genome)
#     $scriptdir/1a_PairwiseAlignmentComponent.sh $tname $joined_genome $target_size $qname $query_genome $query_size $aligndir $genomedir $scriptdir
# done

# Make multiple alignment file from individual alignments
# NOTE: Needed to modify names in the tree file to match those in the maf file names exactly
# for pairwise in $(ls $aligndir/1e_*/1h_*.multialign.sorted.maf); do
#     echo $pairwise
# #   # Get file name in shape for roast (I hate programs with fixed file name requirements)
#   outfile=$(basename $pairwise)
#   outfile=${outfile/multialign.sorted.maf/toast2.maf}
#   outfile=${outfile/1h_${tname}_/$tname.}
#   echo $outfile
#   cp $pairwise $combinedir/$outfile
# done
# cd $combinedir
# roast + X=2 E=$tname "${TREE}" $combinedir/*.toast2.maf 2_combined.roast.maf
# roast + X=2 E=$tname "${TREE2}" $combinedir/*.toast2.maf 2_combined_asteridSmall.roast.maf



$scriptdir/runPhastCons.sh $aligndir $combinedir/2_combined.roast.maf $joined_genome
