# edna-bioinformatics
# this is a script to outline the steps I've taken to set up my comparison of bioinformatic pipelines 
# start by going to the hpc
lkdmac:~ liadomke$ ssh lkdomke@chinook04.alaska.edu

# then create new directories to hold the files
chinook04.alaska.edu % cd /center1/FISH694/lkdomke
chinook04.alaska.edu % mkdir edna

# move raw zipped files from desktop to chinook 
# THIS has to be done from the DESKTOP (local host) to the REMOTE. 
lkdmac:Data liadomke$ pwd
/Users/liadomke/Desktop/Grad School/UAF/Coursework/FISH694_Genom-Bioinform/FISH694/Data # pulling from here
lkdmac:Data liadomke$ scp 20220126_LiaMiFisheDNA.zip lkdomke@chinook04.alaska.edu:/center1/FISH694/lkdomke/edna

# unzip the files, had to use unzip instead of tar -xvzf because it was .zip file
chinook04.alaska.edu % unzip edna/20220126_LiaMiFisheDNA.zip
Archive:  edna/20220126_LiaMiFisheDNA.zip

# move to the edna directory
chinook04.alaska.edu % mv  20220126_LiaMiFisheDNA/ edna

# okay so now we have unzip,  untrimmed files - we need to remove the primer section
# can do this using a conda cutadapt env
# set this up:
chinook04.alaska.edu % conda create -n cutadapt
 chinook04.alaska.edu % conda install -n cutadapt cutadapt
 (base) chinook04.alaska.edu % conda activate cutadapt
# version is 4.3 
 # we need to know the primer sequence for MiFish to remove it: 
 # FWD 5'-GTCGGTAAAACTCGTGCCAGC
 # REV 5'-CATAGTGGGGTATCTAATCCCAGTTTG
 # one thing should be noted - these files are already demultiplexed. There are triplicate presented here


# we have the adaptor sequence from the 5 prime end so we want to make sure to use the -g -G flags rather than -a/-A
# see cutadapt -h if needed, -p indicates the paired end files since we're providing the froward R1 and reverse reads R2
cutadapt -g GTCGGTAAAACTCGTGCCAGC -G CATAGTGGGGTATCTAATCCCAGTTTGG -o e00005-A_S87_L001_R1_trimmed.fastq.gz -p e00005-A_S87_L001_R2_trimmed.fastq.gz e00005-A_S87_L001_R1_001.fastq.gz e00005-A_S87_L001_R2_001.fastq.gz 
#This is cutadapt 4.3 with Python 3.10.10, right now this just put it in the same directory
# okay so we ran that on one and the test looks pretty good

# Processing paired-end reads on 1 core ...
# Done           00:00:00         4,708 reads @  62.1 µs/read;   0.97 M reads/minute
# Finished in 0.306 s (65.080 µs/read; 0.92 M reads/minute).

# === Summary ===

# Total read pairs processed:              4,708
#   Read 1 with adapter:                   4,699 (99.8%)
#   Read 2 with adapter:                   4,693 (99.7%)
# Pairs written (passing filters):         4,708 (100.0%)

# Total basepairs processed:     1,419,452 bp
#   Read 1:       709,153 bp
#   Read 2:       710,299 bp
# Total written (filtered):      1,189,405 bp (83.8%)
#   Read 1:       610,481 bp
#   Read 2:       578,924 bp

# === First read: Adapter 1 ===

# Sequence: GTCGGTAAAACTCGTGCCAGC; Type: regular 5'; Length: 21; Trimmed: 4699 times

# Minimum overlap: 3
# No. of allowed errors:
# 1-9 bp: 0; 10-19 bp: 1; 20-21 bp: 2

# Overview of removed sequences
# length  count   expect  max.err error counts
# 20      20      0.0     2       0 19 1
# 21      4666    0.0     2       4179 443 44
# 22      13      0.0     2       2 7 4

# === Second read: Adapter 2 ===

# Sequence: CATAGTGGGGTATCTAATCCCAGTTTGG; Type: regular 5'; Length: 28; Trimmed: 4693 times

# Minimum overlap: 3
# No. of allowed errors:
# 1-9 bp: 0; 10-19 bp: 1; 20-28 bp: 2

# Overview of removed sequences
# length  count   expect  max.err error counts
# 25      1       0.0     2       0 1
# 27      34      0.0     2       0 3 31
# 28      4650    0.0     2       4 4511 135
# 29      8       0.0     2       0 0 8

# Now we've done it for one - lets create a folder to store the trimmed sequences
# we also need to create a loop to repeat this operation across all sequences
mkdir trimmed
# first we dont have a file that has all the sample names need to do that
mv *_trimmed.fastq.gz trimmed # move any trimmed files (from experiementing) to the new dir
# just kidding - if we have the trimmed directory in the tree with the 20220126_LiaMiFisheDNA dir then it copies the file names
# within trimmed too 
# lets delete it for now
rm -r trimmed
# again jsut kidding this was in the pset7 too
ls *_L001_R1_001.fastq.gz | cut -f '1 2'  -d "_" > sample-names.txt
# this creates a new file in the same step of pulling all the .fastq.gz file names into this sample-names text file
# can check this with cat or nano

# now create loop, note this is run within the 20220126_LiaMiFisheDNA dir 
# MAKE SURE CUTADAPT is activated
for sample in $(cat sample-names.txt) 
do 

    echo "On sample: $sample"

    cutadapt -g GTCGGTAAAACTCGTGCCAGC \
             -G CATAGTGGGGTATCTAATCCCAGTTTGG \
             -o ${sample}_L001_R1_trimmed.fastq.gz -p ${sample}_L001_R2_trimmed.fastq.gz \
             ${sample}_L001_R1_001.fastq.gz ${sample}_L001_R2_001.fastq.gz \
             >> cutadapt_primer_trimming.txt 2>&1

done

# this takes every filename instance in sample-names and inserts it where $sample is 
# then runs the cutadapt with the MiFish primer region
# then puts all the output in a "cutadapt_primer_trimming" txt file, not sure exactly what 2>&1 does
# this is based on the loop in the pset7 metabarcoding assignment

paste sample-names.txt <(grep "passing" cutadapt_primer_trimming.txt | cut -f3 -d "(" | tr -d ")") <(grep "filtered" cutadapt_primer_trimming.txt | cut -f3 -d "(" | tr -d ")")

# this tells you the fraction of reads retained in samples (column 2)
# also the fraction of bps retained in samples (column 3)
# e00005-A_S87    100.0%  83.8%
# e00005-B_S183   100.0%  83.8%
# e00005-C_S279   100.0%  83.8%
# e00007-A_S88    100.0%  83.8%
# e00007-B_S184   100.0%  83.8%
# e00007-C_S280   100.0%  83.8%
# e00010-A_S89    100.0%  83.8%
# e00010-B_S185   100.0%  83.8%
# e00010-C_S281   100.0%  83.8%
# e00011-A_S90    100.0%  83.8%
# e00011-B_S186   100.0%  83.8%
# e00011-C_S282   100.0%  83.8%
# e00012-A_S91    100.0%  83.8%
# e00012-B_S187   100.0%  83.8%
# e00012-C_S283   100.0%  83.8%
# e00013-A_S92    100.0%  83.8%
# e00013-B_S188   100.0%  83.8%
# e00013-C_S284   100.0%  83.8%
# e00014-A_S93    100.0%  83.8%
# e00014-B_S189   100.0%  83.8%
# e00014-C_S285   100.0%  83.8%
# e00015-A_S94    100.0%  83.8%
# e00015-B_S190   100.0%  83.8%
# e00015-C_S286   100.0%  83.8%
# e00016-A_S95    100.0%  83.8%
# e00016-B_S191   100.0%  83.8%
# e00016-C_S287   100.0%  83.8%
# e00017-A_S96    100.0%  95.9%
# e00017-B_S192   100.0%  83.8%
# e00017-C_S288   100.0%  87.3%
# e00289-A_S1     100.0%  83.7%
# e00289-B_S97    100.0%  83.8%
# e00289-C_S193   100.0%  83.7%
# e00290-A_S2     100.0%  83.8%
# e00290-B_S98    100.0%  91.9%
# e00290-C_S194   100.0%  86.5%
# e00291-A_S3     100.0%  83.7%
# e00291-B_S99    100.0%  83.8%
# e00291-C_S195   100.0%  83.7%
# e00292-A_S4     100.0%  83.7%
# e00292-B_S100   100.0%  83.7%
# e00292-C_S196   100.0%  83.8%
# e00293-A_S5     100.0%  83.7%
# e00293-B_S101   100.0%  83.7%
# e00293-C_S197   100.0%  83.7%
# e00294-A_S6     100.0%  83.8%
# e00294-B_S102   100.0%  83.7%
# e00294-C_S198   100.0%  83.7%
# e00295-A_S7     100.0%  83.7%
# e00295-B_S103   100.0%  83.7%
# e00295-C_S199   100.0%  83.7%
# e00296-A_S8     100.0%  83.7%
# e00296-B_S104   100.0%  83.7%
# e00296-C_S200   100.0%  83.8%
# e00297-A_S9     100.0%  83.7%
# e00297-B_S105   100.0%  83.7%
# e00297-C_S201   100.0%  83.7%
# e00298-A_S10    100.0%  83.8%
# e00298-B_S106   100.0%  83.8%
# e00298-C_S202   100.0%  83.7%
# e00299-A_S11    100.0%  83.7%
# e00299-B_S107   100.0%  83.7%
# e00299-C_S203   100.0%  83.7%
# e00300-A_S12    100.0%  83.8%
# e00300-B_S108   100.0%  83.7%
# e00300-C_S204   100.0%  83.7%
# e00301-A_S13    100.0%  83.7%
# e00301-B_S109   100.0%  83.7%
# e00301-C_S205   100.0%  83.7%
# e00302-A_S14    100.0%  83.7%
# e00302-B_S110   100.0%  83.7%
# e00302-C_S206   100.0%  83.7%
# e00303-A_S15    100.0%  83.7%
# e00303-B_S111   100.0%  83.7%
# e00303-C_S207   100.0%  83.8%
# e00304-A_S16    100.0%  83.8%
# e00304-B_S112   100.0%  83.8%
# e00304-C_S208   100.0%  83.7%
# e00305-A_S17    100.0%  83.8%
# e00305-B_S113   100.0%  83.8%
# e00305-C_S209   100.0%  83.8%
# e00306-A_S18    100.0%  83.8%
# e00306-B_S114   100.0%  83.8%
# e00306-C_S210   100.0%  83.8%
# e00307-A_S19    100.0%  83.7%
# e00307-B_S115   100.0%  83.7%
# e00307-C_S211   100.0%  83.7%
# e00308-A_S20    100.0%  83.7%
# e00308-B_S116   100.0%  83.8%
# e00308-C_S212   100.0%  83.8%
# e00309-A_S21    100.0%  83.7%
# e00309-B_S117   100.0%  83.7%
# e00309-C_S213   100.0%  83.7%
# e00310-A_S22    100.0%  83.7%
# e00310-B_S118   100.0%  83.7%
# e00310-C_S214   100.0%  83.7%
# e00311-A_S23    100.0%  83.8%
# e00311-B_S119   100.0%  83.7%
# e00311-C_S215   100.0%  83.7%
# e00312-A_S24    100.0%  83.8%
# e00312-B_S120   100.0%  83.8%
# e00312-C_S216   100.0%  84.4%
# e00313-A_S25    100.0%  83.8%
# e00313-B_S121   100.0%  83.8%
# e00313-C_S217   100.0%  83.8%
# e00314-A_S26    100.0%  83.7%
# e00314-B_S122   100.0%  83.7%
# e00314-C_S218   100.0%  83.8%
# e00315-A_S27    100.0%  83.7%
# e00315-B_S123   100.0%  83.7%
# e00315-C_S219   100.0%  83.8%
# e00316-A_S28    100.0%  83.7%
# e00316-B_S124   100.0%  83.7%
# e00316-C_S220   100.0%  83.8%
# e00317-A_S30    100.0%  96.5%
# e00317-B_S126   100.0%  83.8%
# e00317-C_S222   100.0%  83.8%
# e00318-A_S31    100.0%  83.8%
# e00318-B_S127   100.0%  83.8%
# e00318-C_S223   100.0%  83.8%
# e00319-A_S32    100.0%  83.8%
# e00319-B_S128   100.0%  83.7%
# e00319-C_S224   100.0%  83.7%
# e00320-A_S33    100.0%  83.7%
# e00320-B_S129   100.0%  83.7%
# e00320-C_S225   100.0%  83.7%
# e00321-A_S34    100.0%  83.8%
# e00321-B_S130   100.0%  83.8%
# e00321-C_S226   100.0%  83.9%
# e00322-A_S35    100.0%  83.7%
# e00322-B_S131   100.0%  83.7%
# e00322-C_S227   100.0%  83.7%
# e00323-A_S36    100.0%  83.7%
# e00323-B_S132   100.0%  83.7%
# e00323-C_S228   100.0%  83.8%
# e00324-A_S37    100.0%  83.7%
# e00324-B_S133   100.0%  83.7%
# e00324-C_S229   100.0%  83.7%
# e00325-A_S38    100.0%  83.8%
# e00325-B_S134   100.0%  83.7%
# e00325-C_S230   100.0%  83.8%
# e00326-A_S39    100.0%  83.8%
# e00326-B_S135   100.0%  83.9%
# e00326-C_S231   100.0%  90.1%
# e00327-A_S40    100.0%  83.7%
# e00327-B_S136   100.0%  85.6%
# e00327-C_S232   100.0%  83.8%
# e00328-A_S41    100.0%  83.8%
# e00328-B_S137   100.0%  83.8%
# e00328-C_S233   100.0%  83.8%
# e00329-A_S42    100.0%  83.7%
# e00329-B_S138   100.0%  83.8%
# e00329-C_S234   100.0%  97.1%
# e00330-A_S43    100.0%  83.7%
# e00330-B_S139   100.0%  83.8%
# e00330-C_S235   100.0%  83.7%
# e00331-A_S44    100.0%  83.8%
# e00331-B_S140   100.0%  83.8%
# e00331-C_S236   100.0%  83.8%
# e00332-A_S45    100.0%  100.0%
# e00332-B_S141   100.0%  87.8%
# e00332-C_S237   100.0%  83.8%
# e00333-A_S46    100.0%  86.0%
# e00333-B_S142   100.0%  91.9%
# e00333-C_S238   100.0%  83.8%
# e00334-A_S47    100.0%  83.7%
# e00334-B_S143   100.0%  83.7%
# e00334-C_S239   100.0%  83.7%
# e00335-A_S48    100.0%  83.7%
# e00335-B_S144   100.0%  83.8%
# e00335-C_S240   100.0%  83.8%
# e00336-A_S49    100.0%  83.8%
# e00336-B_S145   100.0%  93.5%
# e00336-C_S241   100.0%  83.7%
# e00337-A_S50    100.0%  83.8%
# e00337-B_S146   100.0%  88.3%
# e00337-C_S242   100.0%  83.7%
# e00338-A_S51    100.0%  83.8%
# e00338-B_S147   100.0%  84.0%
# e00338-C_S243   100.0%  83.8%
# e00339-A_S52    100.0%  83.7%
# e00339-B_S148   100.0%  83.7%
# e00339-C_S244   100.0%  83.7%
# e00340-A_S53    100.0%  83.7%
# e00340-B_S149   100.0%  83.7%
# e00340-C_S245   100.0%  83.8%
# e00341-A_S54    100.0%  83.7%
# e00341-B_S150   100.0%  83.8%
# e00341-C_S246   100.0%  83.8%
# e00342-A_S55    100.0%  85.0%
# e00342-B_S151   100.0%  83.7%
# e00342-C_S247   100.0%  83.7%
# e00343-A_S56    100.0%  83.8%
# e00343-B_S152   100.0%  83.8%
# e00343-C_S248   100.0%  83.8%
# e00344-A_S57    100.0%  83.7%
# e00344-B_S153   100.0%  83.8%
# e00344-C_S249   100.0%  83.7%
# e00345-A_S58    100.0%  83.7%
# e00345-B_S154   100.0%  83.8%
# e00345-C_S250   100.0%  83.7%
# e00346-A_S59    100.0%  83.7%
# e00346-B_S155   100.0%  83.7%
# e00346-C_S251   100.0%  83.7%
# e00347-A_S61    100.0%  83.7%
# e00347-B_S157   100.0%  83.7%
# e00347-C_S253   100.0%  83.7%
# e00348-A_S62    100.0%  83.8%
# e00348-B_S158   100.0%  83.7%
# e00348-C_S254   100.0%  83.7%
# e00349-A_S63    100.0%  83.8%
# e00349-B_S159   100.0%  84.9%
# e00349-C_S255   100.0%  83.7%
# e00350-A_S64    100.0%  83.7%
# e00350-B_S160   100.0%  83.7%
# e00350-C_S256   100.0%  83.7%
# e00351-A_S65    100.0%  83.7%
# e00351-B_S161   100.0%  83.7%
# e00351-C_S257   100.0%  83.7%
# e00352-A_S66    100.0%  83.7%
# e00352-B_S162   100.0%  83.7%
# e00352-C_S258   100.0%  83.7%
# e00353-A_S67    100.0%  83.8%
# e00353-B_S163   100.0%  83.8%
# e00353-C_S259   100.0%  85.8%
# e00354-A_S68    100.0%  84.0%
# e00354-B_S164   100.0%  84.0%
# e00354-C_S260   100.0%  84.0%
# e00355-A_S69    100.0%  83.8%
# e00355-B_S165   100.0%  83.8%
# e00355-C_S261   100.0%  83.8%
# e00356-A_S70    100.0%  83.7%
# e00356-B_S166   100.0%  83.7%
# e00356-C_S262   100.0%  83.8%
# e00357-A_S71    100.0%  83.7%
# e00357-B_S167   100.0%  83.7%
# e00357-C_S263   100.0%  83.8%
# e00358-A_S72    100.0%  83.7%
# e00358-B_S168   100.0%  83.7%
# e00358-C_S264   100.0%  83.7%
# e00359-A_S73    100.0%  83.8%
# e00359-B_S169   100.0%  83.8%
# e00359-C_S265   100.0%  83.6%
# e00360-A_S74    100.0%  83.7%
# e00360-B_S170   100.0%  83.7%
# e00360-C_S266   100.0%  83.7%
# e00361-A_S75    100.0%  83.7%
# e00361-B_S171   100.0%  83.7%
# e00361-C_S267   100.0%  83.7%
# e00362-A_S76    100.0%  83.7%
# e00362-B_S172   100.0%  83.7%
# e00362-C_S268   100.0%  83.7%
# e00363-A_S77    100.0%  83.7%
# e00363-B_S173   100.0%  83.7%
# e00363-C_S269   100.0%  83.7%
# e00364-A_S78    100.0%  83.7%
# e00364-B_S174   100.0%  83.7%
# e00364-C_S270   100.0%  83.7%
# e00365-A_S79    100.0%  83.8%
# e00365-B_S175   100.0%  83.7%
# e00365-C_S271   100.0%  83.7%
# e00366-A_S80    100.0%  100.0%
# e00366-B_S176   100.0%  84.7%
# e00366-C_S272   100.0%  93.0%
# e00367-A_S81    100.0%  83.7%
# e00367-B_S177   100.0%  83.7%
# e00367-C_S273   100.0%  83.7%
# e00368-A_S82    100.0%  83.7%
# e00368-B_S178   100.0%  83.7%
# e00368-C_S274   100.0%  83.7%
# e00369-A_S83    100.0%  83.8%
# e00369-B_S179   100.0%  83.7%
# e00369-C_S275   100.0%  83.7%
# e00370-A_S84    100.0%  84.1%
# e00370-B_S180   100.0%  87.8%
# e00370-C_S276   100.0%  83.8%
# e00371-A_S85    100.0%  100.0%
# e00371-B_S181   100.0%  100.0%
# e00371-C_S277   100.0%  83.7%
# e00372-A_S86    100.0%  83.8%
# e00372-B_S182   100.0%  83.8%
# e00372-C_S278   100.0%  83.8%
# NEGATIVE-Plate1-A_S29   100.0%  83.7%
# NEGATIVE-Plate1-B_S125  100.0%  83.7%
# NEGATIVE-Plate1-C_S221  100.0%  83.7%
# POSITIVE-Plate1-A_S60   100.0%  95.1%
# POSITIVE-Plate1-B_S156
# POSITIVE-Plate1-C_S252
# Undetermined_S0

# unzip files using pigz while they're in the trimmed file
(cutadapt) chinook04.alaska.edu % pigz -d trimmed/*.gz
# now we need to move all the trimmed samples into one file location
# we likely need to create the ASV table within HPC using R 
# note that we'll be pulling this from https://benjjneb.github.io/dada2/tutorial.html
# creating r_env within conda
conda create -n r_env r-essentials r-base

# it said again that there wasn't enough space in the home directory for another envs
# so I'm removing the miniconda3 directory from within the $HOME directory and moving it to the CENTER1
# this worked so we now have r-base and r-essentials loaded into the r_env 
R -Version 
# lets us know we're running R version 4.2.3
# want to install dada2 for use with our trimmed data
install.packages("dada2")
# prompts for selecting our CRAN for this session, 74 is OR in USA
# can't install dada2 with the most up to date version of R
# need to delete the r_env and re install it with an older version of R
conda create -n r_env
conda install -c conda-forge r-base=4.0.0 r-essentials

# okay great now we're in R and it has the version of 4.0.0
# install DADA2
# had to go through devtools
install.packages("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.16")
# which then required updating a LOT of packages (~52)
path <- "~/MiSeq_SOP" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_trimmed.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_trimmed.fastq", full.names = TRUE))
_R1.fastq
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:2]) 
plotQualityProfile(fnRs[1:2]) 
# I'm having issues with this showing up with a popout image. 

## OKAY - so actually turns out my computer DOES have enough computing power to run DADA2 with my samples in R on my desktop rather than the HPC
# so since I have a base script for this from the homework and the tutorial information online lets do that instead (https://benjjneb.github.io/dada2/tutorial.html)
# have to move the trimmed files to the project directory
# need to do this FROM that directory (/Users/liadomke/Desktop/Grad School/UAF/Coursework/FISH694_Genom-Bioinform/FISH694/Data/eDNA_files)
scp -r lkdomke@chinook04.alaska.edu:/center1/FISH694/lkdomke/edna . 

# go to 01_DADA2_bioinformatics.Rmd - complete the steps there BEFORE coming back to this script. 


## MAKE SURE YOU'VE DONE THE OTHER STEPS IN DADA2 BEFORE RUNNING THE SCRIPT BELOW (you won't be able to anyway)

# output of the r scripts were two csv files
# MiFish_ASV_seqtab_nochim.csv and MiFish_ASVtable.csv within a csv_outputs folder
# move folder to hpc
scp -r csv_outputs lkdomke@chinook04.alaska.edu:/center1/FISH694/lkdomke/edna

# which worked
# so now we need to do PIPELINE VERSION ONE
# which requires we use these ASV tables to "blast" or query against known files

# first we need the blastn 
# I chose to put it in my cutadapt env
conda activate cutadapt
conda install blast

# also need to download the nucleotide data base to compare it to nt, this script is within the blast program
update_blastdb.pl --decompress nt

# then you have to appropriately set up the nucleotide database using: 
# since you installed blast in the cutadapt env this should work
# see makeblastdb -help if you're unsure about format
makeblastdb -in [input database] -out [output database] -dbtype nucl 
# fill in the []
makeblastdb -in nt -out nt_db -dbtype nucl

# this is the blastn base script - remember not good to run things in the login node, trouble shoot using the partition = debug
nohup blastn -db nt -query MiFish_ASV_seqtab_nochim.csv -perc_identity 96 -qcov_hsp_perc 100 -out MiFish_ASV_blast.out -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames staxids"

#maybe that should've been run in a script?
# what woudl that look like: 

# first module load slurm
module load slurm
# make file 
nano blastn
#!/bin/bash

#SBATCH --partition=t1standard
#SBATCH --ntasks=24
#SBATCH --tasks-per-node=24
#SBATCH --exclusive
#SBATCH --time=14:00:00
#SBATCH --mem=120G
#SBATCH --mail-user=lkdomke@alaska.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o ouput_blastn.txt


conda activate cutadapt

nohup blastn -remote -db nt -query MiFish_ASV_seqtab_nochim.csv \
-perc_identity 96 -qcov_hsp_perc 100 -out MiFish_ASV_blast.out  \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames staxids"

nohup blastn -db nt -query MiFish_ASV_seqtab_nochim.csv \
 -perc_identity 96 -qcov_hsp_perc 100 -num_threads 10 -out MiFish_ASV_blast.out \
 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames staxids"


# save blastn.sh
# make executable 
chmod +x blastn.sh

# sbatch
sbatch blastn.sh

# the output of the blastn will not have taxonomic lineage so we have to use another tool kit (taxonkit) to get those
conda activate cutadapt 
conda install taxonkit 

cat MiFish_ASV_blast.out | taxonkit lineage -c -i 14 > MiFish_96perc_taxlineage.txt

taxonkit reformat MiFish_96perc_taxlineage.txt -i 16 > MiFish_96perc_shorttaxlineage.txt


# move these files back to the r project
# go to your desired location /Users/liadomke/Desktop/Grad School/UAF/Coursework/FISH694_Genom-Bioinform/FISH694/Data/blast_DADA2_outputs
scp -r lkdomke@chinook04.alaska.edu:/center1/FISH694/lkdomke/edna .