#! /bin/bash
# original author - Jason Wallace. Modified for Lynsey's bladderwort analysis

# Make multiple genome alignment and do GERP analysis
# # # ROAST=$HOME/Software/Aligners/multiz-tba.012109/roast
# GERP=$HOME/Software/Aligners/gerp++/
# MULTIZ=$HOME/Software/Aligners/multiz-tba.012109/
# IGVTOOLS="bash /home/jgwall/Software/Browsers/IGVTools/igvtools"

#specify number of processors
PROCS=48

# Set up directories
parentdir=/scratch/lk82153/jwlab/Bladderwort
datadir=$parentdir/0_Data
scriptdir=/scratch/lk82153/jwlab/Repositories/bladderwort-analysis/scripts/
# parentdir=$HOME/Work/Bladderwort/7_multiSpeciesAlignment
# datadir=$HOME/Work/Bladderwort/0_Data
genomedir=$datadir/Genomes/Asterids
sourcedir=$datadir/Genomes/Tomato
aligndir=7_PairwiseAlignment
combinedir=$aligndir/Multiple_Alignment
# gerpdir=3_GERP
# wigdir=3_GERP/3c_WiggleFiles

if [ ! -e $aligndir ]; then mkdir $aligndir; fi
if [ ! -e $combinedir ]; then mkdir $combinedir; fi
# if [ ! -e $gerpdir ]; then mkdir $gerpdir; fi
# if [ ! -e $wigdir ]; then mkdir $wigdir; fi

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
TREE="((((Vaccinium_corymbosum,Vaccinium_macrocarpon)Vaccinium,Actinidia_chinensis)Ericales,((Lactuca_sativa,(Helianthus_annuus,Erigeron_canadensis)Asteroideae)Asteraceae,((Cuscuta_australis,((Capsicum_annuum,Solanum_lycopersicum)Solanoideae,(Nicotiana_tabacum,(Petunia_axillaris,Petunia_integrifolia_subsp._inflata)Petunia))Solanaceae)Solanales,(Coffea_canephora,((Utricularia_gibba,Genlisea_aurea)Lentibulariaceae,Erythranthe_guttata)Lamiales))lamiids))asterids);"


###########
# Pairwise genome alignments - Split apart for simplicity
###########

tname=SolanumLycopersicum
target_genome=$sourcedir/SolanumLycopersicum.GCF_000188115.4_SL3.0_genomic.fna
target_genome_stem=SolanumLycopersicum.GCF_000188115.4_SL3.0_genomic
target_size=1000000000 # 1GB = bigger than any utricularia chromosome
query_size=1000000  # 1 MB segments for each query

# # Join the Utricularia genome into a  single scaffold for convenience
python3 1_JoinFastaWithPad.py -i $target_genome -o $sourcedir/$target_genome_stem.all_in_one.fna -k $sourcedir/$target_genome_stem.all_in_one.key.txt -n 1000 --newname utricularia_joined #--debug

# for query_genome in $(ls $genomedir/ | sed "s/\.fna//g"); do
#     qname=${qname/.*/}
#     
#     ./1a_PairwiseAlignmentComponent.sh $tname $target_genome $target_size $qname $query_genome $query_size $aligndir
#     
# done

# Make multiple alignment file from individual alignments
# PATH=$PATH:$HOME/Software/Aligners/multiz-tba.012109/   # Add utilities to PATH
# for pairwise in $aligndir/1e_*/1h_*.multialign.maf; do
#   # Get file name in shape for roast (I hate programs with fixed file name requirements)
#   outfile=$(basename $pairwise)
#   outfile=${outfile/multialign.maf/toast2.maf}
#   outfile=${outfile/1h_${tname}_/$tname.}
#   echo $outfile
#   cp $pairwise $combinedir/$outfile
# done
# cd $combinedir
# roast + X=2 E=$tname "${TREE}" $combinedir/*.toast2.maf 2_combined.roast.maf


# # # # DEPRECATED - Make phylogenetic tree with 4-fold degenerate sites for neutral evolution model
# # # #  GFF file does not contain start/stop codon locations, just transcripts. Too easy to mess that up, so abandoned this and will just take random sites
# # # # sed -r -e "s/>lcl\|/>/" $orig_genome > $gerpdir/3_genome_names_stripped.fa    # strip "lcl" from beginning of contig names
# # # python3 3a_GetFourfoldDegenerateSites.py -f $gerpdir/3_genome_names_stripped.fa -g $sourcedir/Utricularia_gibba.4.1.windowmasked.short.gff -o $gerpdir/3a_fourfold_sites.txt 

# Make phylogenetic tree with N random sites, using make_phylogeny from QIIME
nsites=10000 # Number of random sites to select
# python3 3a_SelectRandomSitesFromMaf.py -i $combinedir/2_combined.roast.maf -o $gerpdir/3a_combined.random_sites.$nsites.txt \ 
#   --fasta $gerpdir/3a_combined.random_sites.$nsites.fasta --n-taxa 15 --n-sites $nsites --seed 1  #--debug
# make_phylogeny.py -i $gerpdir/3a_combined.random_sites.$nsites.fasta -o $gerpdir/3a_combined.random_sites.$nsites.tre

# NOTE: Had to manually modify the tree file to remove quotes and support numbers so would be parsed by GERP correctly

# # Run GERP
# $GERP/gerpcol -f $combinedir/2_combined.roast.maf -t $gerpdir/3a_combined.random_sites.$nsites.mod.tre -e utricularia   # -j    # -j projects reference sequence; not sure what that means
# $GERP/gerpcol -f $combinedir/2_combined.roast.rates
# mv $combinedir/2_combined.roast.maf.rates $gerpdir/3c_combined.roast.maf.rates
# mv $combinedir/2_combined.roast.maf.rates.elems $gerpdir/3c_combined.roast.maf.rates.elems
# mv $combinedir/2_combined.roast.maf.rates.exclude $gerpdir/3c_combined.roast.maf.rates.exclude

# Convert GERP output for loading into IGV
# offset=$(grep "^s" 2_MultipleAlignment/2_combined.roast.maf | head -n1 | tr -s ' ' | cut -f3 -d' ') # Complex code to get the start position of the first alignment
# offset=$(($offset-1))   # remove 1 because the above gets the _start_ position, when I need the number of nt _before_ that
# python3 3c_ConvertGerpElemsToBed.py -i $gerpdir/3c_combined.roast.maf.rates.elems --merge-key $sourcedir/Utricularia_gibba.4.1.windowmasked.all_in_one.key.txt \
#   -o $gerpdir/3c_combined.roast.maf.rates.elems.bed --offset $offset
# TODO: ~30 elements seem to be missing due to corner cases where they're just barely at the end (1 nt)
# python3 3c_ConvertGerpRatesToWig.py -i $gerpdir/3c_combined.roast.maf.rates --merge-key $sourcedir/Utricularia_gibba.4.1.windowmasked.all_in_one.key.txt \
#   -o $wigdir/3c_ --offset $offset #--debug
# cat $wigdir/3c_*.wig > $gerpdir/3d_gerp_scores.wig

# # # # # # TODO: Run GERP analysis
# # # # # # python3 3b_ConcatenateMafIntoMultifasta.py -i $combinedir/2_combined.roast.maf -o $combinedir/2_combined.roast.aln -t $tname -p '-' --keyfile $combinedir/2b_multifasta_key.txt
# # # # # # $GERP/gerpcol -f $combinedir/2_combined.roast.aln -t  $gerpdir/3a_combined.random_sites.$nsites.mod.tre -e utricularia -a
# # # # # # mv $combinedir/2_combined.roast.aln.rates $gerpdir/3a_combined.roast.aln.rates
# # # # # # $GERP/gerpelem -f 3a_combined.roast.aln.rates
# # # # # 
