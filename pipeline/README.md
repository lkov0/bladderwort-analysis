## This directory contains wrapper scripts for running analyses on the Bladderwort genome

### File descriptions:

- `0_MiSeqNanoQC.sh` - contains code to run an initial check to make sure the short read data would produce adequate genomic coverage for a short read genome assembly
- `0_runBladderwortAnalysis.sh` - A wrapper script for the analysis performed. input: raw sequencing data, final output: regions of interest and metadata from the subsequent steps + some exploratory analyses on the sequences themselves.
- `1_runAlignment.sh` - Initial QC and alignment of 3' sequencing reads to the scaffolded *Utricularia gibba* assembly.
- `4_findMotifs.sh` - All downstream analyses on the identified candidate intergenic sequences (i.e. MEME suite, intersections with putative conserved regions from runMultigenomeAlignment_*.sh)
- `runGenomeAssembly.sh` - Code used to generate and analyze the *Utricularia gibba* genome and corresponding Maker annotation.
- `runMultigenomeAlignment_*.sh`- Code used to generate multigenome alignments for asterids, eudicots, monocots and rosids.
