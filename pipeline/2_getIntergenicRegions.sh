# 2_getIntergenicRegions.sh

parentDir=$1
scriptDir=$2
dataDir=$3
analysisDir=$4
memeDir=$5

NPROCS=$6

intersectDir=$memeDir/Intersect

if [ ! -e $intersectDir ]; then mkdir $intersectDir; fi

module load R

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
