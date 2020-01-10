#! /bin/bash
# original author - Jason Wallace. Modified for Lynsey's bladderwort analysis

# Make multiple genome alignment and do GERP analysis
# # # ROAST=$HOME/Software/Aligners/multiz-tba.012109/roast
# GERP=$HOME/Software/Aligners/gerp++/
# MULTIZ=$HOME/Software/Aligners/multiz-tba.012109/
# IGVTOOLS="bash /home/jgwall/Software/Browsers/IGVTools/igvtools"

#specify number of processors
PROCS=48

module load Biopython/1.68-foss-2016b-Python-3.5.2
module load ucsc/359

# Set up directories
parentdir=/scratch/lk82153/jwlab/Bladderwort_3pri
scriptdir=/home/lk82153/Repositories/bladderwort-analysis/scripts
# parentdir=/scratch/lk82153/jwlab/Bladderwort
# scriptdir=/scratch/lk82153/jwlab/Repositories/bladderwort-analysis/scripts
datadir=/work/jawlab/data/bladderwort
# parentdir=$HOME/Work/Bladderwort/7_multiSpeciesAlignment
# datadir=$HOME/Work/Bladderwort/0_Data
genomedir=$datadir/genomes/monocots
sourcedir=$datadir/genomes/utricularia
aligndir=$parentdir/6_ConservationAnalysis_monocots
combinedir=$aligndir/Multiple_Alignment_monocots

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
TREEMONOCOT="(SpirodelaPolyrhiza ((SetariaItalica SorghumBicolor)(OryzaSativa BrachypodiumDistachyon)) PhalaenopsisEquestris)"
TREEROSID= "(KalanchoeFedtschenkoi (((CaricaPapaya ArabidopsisThaliana) DurioZibethinus) ((CannabisSativa PrunusPersica) (GlycineMax MedicagoTruncatula) BetulaPendula CitrullusLanatus ManihotEsculenta) VitisVinifera) BetaVulgaris)"
# TREE2="((SolanumLycopersicum) (CoffeaCanephora (UtriculariaGibba MimulusGuttatus GenliseaAurea)))"

###########
# Pairwise genome alignments - Split apart for simplicity
###########

tname=UtriculariaGibba
target_genome=$sourcedir/scaffolds.ugibba_lk.fasta
joined_genome=$sourcedir/scaffolds.ugibba_lk.all_in_one.fna
target_size=1500000000 # 1GB = bigger than any U.gibba chromosome
query_size=1000000  # 1 MB segments for each query

#Join the U.gibba genome into a single scaffold for convenience
# python3 $scriptdir/1_JoinFastaWithPad.py -i $sourcedir/scaffolds.ugibba_lk.fasta -o $joined_genome \
#   -k $sourcedir/scaffolds.ugibba_lk.all_in_one.all_in_one.key.txt -n 1000 --newname ${tname}_joined

# for query_genome in $(ls $genomedir/); do
#     qname=$(sed "s/\..*//g" <<<$query_genome)
#     $scriptdir/1a_PairwiseAlignmentComponent.sh $tname $joined_genome $target_size $qname $query_genome $query_size $aligndir $genomedir $scriptdir
# done

# Make multiple alignment file from individual alignments
# NOTE: Needed to modify names in the tree file to match those in the maf file names exactly
for pairwise in $(ls $aligndir/1e_*/1h_*.multialign.sorted.maf); do
    echo $pairwise
#   # Get file name in shape for roast (I hate programs with fixed file name requirements)
  outfile=$(basename $pairwise)
  outfile=${outfile/multialign.sorted.maf/toast2.maf}
  outfile=${outfile/1h_${tname}_/$tname.}
  echo $outfile
  cp $pairwise $combinedir/$outfile
done
cd $combinedir
roast + X=2 E=$tname "${TREEMONOCOT}" $combinedir/*.toast2.maf 2_combined.roast.maf

# # $scriptdir/runPhastCons.sh $aligndir $combinedir/2_combined.roast.maf $joined_genome

# #fix keyfile
# # sed "s/, whole genome shotgun sequence//g" $sourcedir/UtriculariaGibba.GCA_002189035.1_U_gibba_v2_genomic.all_in_one.key.txt | sed "s/\s[NC].*.Umecuaro /\t/g" | sed "s/ mitochondrial//g" | sed "s/ chloroplast//g" | sed "s/chromosome 1/unitig_0/g" | sed "s/chromosome 2/unitig_22/g" | sed "s/chromosome 3/unitig_26/g" | sed "s/chromosome 4/unitig_32/g" > $sourcedir/UtriculariaGibba.GCA_002189035.1_U_gibba_v2_genomic.all_in_one.key.fixed.txt
# 
