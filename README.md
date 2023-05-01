# FISH694

This repository holds script associated with Fish694 - Genomics and Bioinformatics

Included in the /scripts folder are the scripts associated with the course project: comparing eDNA species assignment using two bioinformatic pipelines
This workflow uses data collected from summer 2021 
(water samples of nearshore fish communities in eelgrass and understory kelp beds in southern Southeast Alaska). 
These data were filtered, extracted, and sequenced before and the scripts in this repository start from the sequenced Illumina data. 
The data *are not included* in this repository (to big to be stored here) - please contact author (Lia Domke) for access to data. 

Steps include: 

edna-bioninformatic.sh script 
  - contains the HPC related steps (used UAF Chinook 04)
  - input: demultiplexed fastq files from Illumina MiSeq (files are in .fastq.gz format)
      - process includes using cutadapt to remove primer adaptors
      - reading out .fastq files

01_DADA2_bioinformatics.Rmd script
  - follows DADA2 tutorial (uses R)
  - input .fastq files that are demultiplexed and have primers removed
  - trimms, processes for errors, then identifies amplicon sequence variants. Removes chimeras etc (see script for more details)
  - output is ASV table and ASV table with no chimeras

Then we processed the same ASV table in two different ways: 
BLASTn and Insect

**BLASTn**
BLASTn needs to be done through the HPC or VM - go back to 01_edna_bioinformatic.sh script
Then outputs from this script need to be cleaned manually (custom R scripts) go to -->

02_blastn_taxonomic_filtering.Rmd script
  - input pooled taxonomy that comes from BLASTn and taxonkit
  - script cleans and identifies the top identified identification (species) for each ASV 
  - if there is conflict (i.e. multiple species with high confidence) then we go up a taxonomic level (species --> genus)
  - if the top hit and second hit are more than 2% diff than we go with top hit
 
 **Insect**
 This can be done in R using pre-classifed model specific by primer (MiFish) - see https://shaunpwilkinson.github.io/post/insect-vignette/
 
 03_insect_bioinformatics.Rmd script
  - input is ASV tables (no-chim) from the DADA2 script output
  - this uses the insect classifier to identify taxons from ASVs
  
 **Final comparison step**
 Lets see how everything did! 
 This uses the outputs from steps 2 and 3 and compares ASV by ASV the taxon identification
 
 04_Bioinformatic_comparison_analysis.Rmd 
  - inputs are the outputs from steps 2 and 3
  - outputs are graphs and tables comparing these two pipeliens. 
