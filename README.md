# Specificity_Calculator

Specificity is a metric that measures how many base pairs of cleaned reads are aligned to target sequences, expressed as a percentage.
These scripts calculate specificity for an entire merged contig (Specificity_Calc_v1.py) or for 
a portion of the contig, such as coding or flanking (Specificity_Calc_Beds_v1.py).  The function for parsing 
the SAM files to count base pairs is also included here.  It worked well with the tests written.

-----------------------------------------------------------------------------------------
Specificity_Calc_v1.py :
Usage: python Specificity_Calc_v1.py [directory with bam files] [directory with associated cleaned reads files]

The 'Sample_indexXX__sorted.bam' files are located in the 'intarget_assemblies' folder resulting 
from script 6 (6-TransExonCapPhylo 'contig' option) of the 
https://github.com/MVZSEQ/denovoTargetCapturePhylogenomics workflow.

The cleaned reads files (see naming scheme below) are created during script 2 (2-ScrubReads).

You may consider copying and pasting each of these sets of files to new directories, as
the output will be written to the directory with bam files.


This script creates a sam file from each bam file, iterates over the lines in the sam file 
to count basepairs (using the defined sam_base_counter to parse out sam codes), then shifts
to cleaned reads and counts bases for files matching bam sample name, then calculates the 
specificity (bam bp / cleaned read bp) and writes to output file in bam directory.

example of bam file names in dir 1:
JMPD002_index1_sorted.bam

example of cleaned file names in dir 2:
JMPD002_index1_1_final
JMPD002_index1_2_final
JMPD002_index1_U_final

removes sam files after each sample (which are large in size, sometimes >1GB!!)

***Make sure for every bam file, there are matching cleaned reads files in the other
directory, or you will crash this script.

------------------------------------------------------------------------------------------
Specificity_Calc_Beds_v1.py:
Usage: python Specificity_Calc_v1.py [directory with bam files] [directory with associated cleaned reads files] [dir with bed files] [output file prefix (ex. Target, Flanking)]

The 'Sample_indexXX__sorted.bam' files are located in the 'intarget_assemblies' folder resulting 
from script 6 (6-TransExonCapPhylo 'contig' option) of the 
https://github.com/MVZSEQ/denovoTargetCapturePhylogenomics workflow.

The cleaned reads files (see naming scheme below) are created during script 2 (2-ScrubReads).

You may consider copying and pasting each of these sets of files to new directories, as
the output will be written to the directory with bam files.

This script creates a sam file from each bam file using appropriate bed file, then
iterates over the lines in the sam file to count basepairs 
(using the defined sam_base_counter to parse out sam codes), then shifts
to cleaned reads and counts bases for files matching bam sample name, then calculates the 
specificity (bam bp / cleaned read bp) and writes to output file in bam directory.

****Only use this version if you want to track specificity for target or flanking,
otherwise the other script will do the target+flanking automatically.

bam file names in dir 1:
JMPD002_index1_sorted.bam

cleaned file names in dir 2:
JMPD002_index1_1_final
JMPD002_index1_2_final
JMPD002_index1_U_final

bed file names in dir 3:
JMPD004_index76_coding.bed
or
JMPD002_index7_flanking_ONLY.bed


removes sam files after each sample (which are large in size, sometimes >1GB)

***Make sure for every bam file, there are matching cleaned reads files in the cleaned reads
directory, or you will crash this script.

***Make sure you only have one bed file per sample in directory (coding or flanking, not both)





##############
DEPENDENCIES:
numpy - Numerical Python
samtools - needs to be in path to call from command line
##############
------------------------
written for Python 2.7.3
Dan Portik
daniel.portik@berkeley.edu
October 2015
------------------------
'''

