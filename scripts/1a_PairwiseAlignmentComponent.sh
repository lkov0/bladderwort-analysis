#! /bin/bash

# Script to do pairwise multiple alignment of two genomes


# Software downloaded from UCSC, from http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/
UCSC=/home/lk82153/jwLab/Programs/UCSC_utils

# arguments
tname=$1
target_genome=$2
target_size=$3
qname=$4
query_genome=$5
query_size=$6
aligndir=$7
qgenomedir=$8
scriptdir=$9

targetdir=$aligndir/1a_${tname}_split
querydir=$aligndir/1a_${qname}_split
pairdir=$aligndir/1b_${tname}_${qname}_pairwise
chaindir=$aligndir/1e_${tname}_${qname}_chains

if [ ! -e $targetdir ]; then mkdir $targetdir; fi
if [ ! -e $querydir ]; then mkdir $querydir; fi
if [ ! -e $pairdir ]; then mkdir $pairdir; fi
if [ ! -e $chaindir ]; then mkdir $chaindir; fi

target_sizes=$chaindir/1e_$tname.combined.sizes
query_sizes=$chaindir/1e_$qname.combined.sizes

# Split target genome by segments
# outprefix=$targetdir/${tname}
# python3 $scriptdir/1a_CutGenomesBySize.py -i $target_genome -o $outprefix --rename $tname --min 1000 --winsize $target_size #--debug
# 
# # Split query genome by segments
# outprefix=$querydir/${qname}
# python3 $scriptdir/1a_CutGenomesBySize.py -i $qgenomedir/$query_genome -o $outprefix --rename $qname --min 1000 --winsize $query_size #--debug

# 
# # Do all pairwise target-query alignments

# lastz_commands=$aligndir/1b_${tname}_${qname}.lastz_commands.txt
# chain_commands=$aligndir/1c_${tname}_${qname}.chain_commands.txt
# echo "echo -e #######\nRunning Pairwise alignments\n#######" > $lastz_commands
# echo "echo -e #######\nConverting lastz output to chains\n#######" > $chain_commands
# i=0
# for target_win in $targetdir/${tname}*.fa; do
#     for query_win in $querydir/${qname}*.fa; do
#     
#       axtfile=$pairdir/1c_${tname}_${qname}.$i.axt 
#       chainfile=$pairdir/1d_${tname}_${qname}.$i.chain
#       
# #       Align with LASTZ and convert to a chain file
#       echo "module load LASTZ; lastz $target_win $query_win --ambiguous=iupac --format=axt --inner=2000 --xdrop=9400 --gappedthresh=3000 --hspthresh=2200 > $axtfile" >> $lastz_commands
#       echo "module load LASTZ; $UCSC/axtChain -linearGap=loose $axtfile -faT $target_win -faQ $query_win $chainfile" >> $chain_commands
#       
#       i=$(($i+1))
#       
#     done
# done
# 
# module load parallel
# 
# cat $lastz_commands | parallel --progress
# cat $chain_commands | parallel --progress

# # ## Get sizes of chromosomes/chromosome segments ##
cat $targetdir/*.fa > $chaindir/1e_$tname.combined.fa
$UCSC/faSize $chaindir/1e_$tname.combined.fa -detailed > $target_sizes
rm $chaindir/1e_$tname.combined.fa

echo ------------ Done with target chr sizes --------------

cat $querydir/*.fa > $chaindir/1e_$qname.combined.fa
$UCSC/faSize $chaindir/1e_$qname.combined.fa -detailed > $query_sizes
rm $chaindir/1e_$qname.combined.fa

echo ------------ Done with query chr sizes --------------

# Sort and filter chains ##    
find $pairdir -maxdepth 1 -name "*.chain" > $chaindir/1e_${tname}_${qname}_all.chains.txt
$UCSC/chainMergeSort -inputList=$chaindir/1e_${tname}_${qname}_all.chains.txt > $chaindir/1e_${tname}_${qname}_all.chain.merge
$UCSC/chainPreNet $chaindir/1e_${tname}_${qname}_all.chain.merge $target_sizes $query_sizes $chaindir/1e_${tname}_${qname}_all.chain.merge.pre

echo ---------- Done with sort and filter chains ----------------

# ## Net the chains ##
$UCSC/chainNet $chaindir/1e_${tname}_${qname}_all.chain.merge.pre -minSpace=1 $target_sizes $query_sizes $chaindir/1f_${tname}_${qname}.target_$tname.net $chaindir/1f_${tname}_${qname}.query_$qname.net

echo ---------- Done with netting of chains ----------------


# ## Add synteny info to nets ##
$UCSC/netSyntenic $chaindir/1f_${tname}_${qname}.target_$tname.net $chaindir/1f_${tname}_${qname}.target_$tname.noClass.net


# Make 2bit files
$UCSC/faToTwoBit $targetdir/${tname}*.fa $chaindir/1g_$tname.seqs.2bit
$UCSC/faToTwoBit $querydir/${qname}*.fa $chaindir/1g_$qname.seqs.2bit

echo ---------- Done with 2bit file making ----------------

# Convert to MAF files
$UCSC/netToAxt $chaindir/1f_${tname}_${qname}.target_$tname.noClass.net $chaindir/1e_${tname}_${qname}_all.chain.merge.pre $chaindir/1g_$tname.seqs.2bit $chaindir/1g_$qname.seqs.2bit \
  $chaindir/1h_${tname}_${qname}.multialign.axt
$UCSC/axtSort $chaindir/1h_${tname}_${qname}.multialign.axt $chaindir/1h_${tname}_${qname}.multialign.sorted.axt
$UCSC/axtToMaf $chaindir/1h_${tname}_${qname}.multialign.sorted.axt $target_sizes $query_sizes $chaindir/1h_${tname}_${qname}.multialign.sorted.maf

echo ---------- Done with MAF making ----------------

echo ---------- Done with species $qname ---------------
