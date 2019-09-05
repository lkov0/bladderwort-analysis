# runPhastCons.sh

alignDir=$1
mafFile=$2
refGenome=$3
phastDir=$alignDir/PhastCons
chunkDir=$phastDir/Chunks
treeDir=$phastDir/Trees
elementDir=$phastDir/Elements
scoreDir=$phastDir/Scores

if [ ! -e $phastDir ] ; then mkdir $phastDir; fi
if [ ! -e $chunkDir ] ; then mkdir $chunkDir; fi
if [ ! -e $treeDir ] ; then mkdir $treeDir; fi
if [ ! -e $scoreDir ] ; then mkdir $scoreDir; fi
if [ ! -e $elementDir ] ; then mkdir $elementDir; fi


# msa_split $mafFile --in-format MAF --refseq $refGenome \
#     --windows 1000000,0 --out-root $chunkDir/asteridAlignment1 --out-format SS \
#     --min-informative 1000 --between-blocks 5000

# phyloFit --msa-format MAF --tree "((HelianthusAnnuus,LactucaSativa),(((((PetuniaAxillaris,(CapsicumAnnuum,SolanumLycopersicum)),CuscutaAustralis),((GenliseaAurea,UtriculariaGibba),MimulusGuttatus)),CoffeaCanephora),(ActinidiaChinensis,VacciniumCorymbosum)));" --out-root $treeDir/phyloFit $mafFile

# treeCommands=$treeDir/treeCommands.txt
# rm $treeCommands
# for file in $(ls $chunkDir/*.ss | sed "s/\.ss//g"); do 
#     echo "phastCons --target-coverage 0.125 --expected-length 20 --gc 0.4 --estimate-trees $file $file.ss $treeDir/phyloFit.mod --no-post-probs" >> $treeCommands
# done
# 
# cat $treeCommands | parallel --progress
# 
# mv $chunkDir/*.mod $treeDir

# cd $phastDir
# ls Trees/*.cons.mod > cons.txt
# phyloBoot --read-mods '*cons.txt' --output-average ave.cons.mod 
# ls $treeDir/*.noncons.mod > noncons.txt
# phyloBoot --read-mods '*noncons.txt' --output-average ave.noncons.mod

rm -f $elementDir/* $scoreDir/*
phastCommands=$phastDir/phastCommands.txt
rm $phastCommands
for file in $(ls $chunkDir/*.ss | sed "s/\.ss//g"); do
    echo "phastCons --target-coverage 0.125 --expected-length 20 --most-conserved $file.bed --score $file.ss $phastDir/ave.cons.mod,$phastDir/ave.noncons.mod > $file.wig" >> $phastCommands
done

cat $phastCommands | parallel --progress

mv $chunkDir/*.bed $elementDir
mv $chunkDir/*.wig $scoreDir
