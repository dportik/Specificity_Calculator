import sys
import os
import subprocess as sp
import shutil
import numpy as np
import time
'''
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
bam_dir = sys.argv[1]
clean_dir = sys.argv[2]
bed_dir = sys.argv[3]
file_prefix = sys.argv[4]

#define very tricky function for parsing out base pairs in sam file column 6 (2S2M8I56M1H, 99M1H, 74M, etc)
#where numbers preceeding an 'M' should be added together for multiple instances of 'M' on each line
def sam_base_counter(x):
    i = int(0)
    line_bp = int(0)
    for char in x:
        if i <= int(1):
            if char == 'M':
                h = i - int(1)
                if base_str[h].isdigit():
                    bp = int(base_str[h])
                    line_bp += bp
                    
        elif i <= int(2) and i >= int(1):
            if char == 'M':
                h = i - int(1)
                g = i - int(2)
                if base_str[h].isdigit() and base_str[g].isdigit():
                    bp = int( str(base_str[g]) + str(base_str[h]))
                    line_bp += bp
                elif base_str[h].isdigit() and base_str[g].isdigit() == False:
                    bp = int(base_str[h])
                    line_bp += bp
                    
        elif i > int(2):
            if char == 'M':
                h = i - int(1)
                g = i - int(2)
                f = i - int(3)
                if base_str[h].isdigit() and base_str[g].isdigit() and base_str[f].isdigit():
                    bp = int( str(base_str[f]) + str(base_str[g]) + str(base_str[h]))
                    line_bp += bp
                elif base_str[h].isdigit() and base_str[g].isdigit() and base_str[f].isdigit() == False:
                    bp = int(str(base_str[g]) + str(base_str[h]))
                    line_bp += bp
                elif base_str[h].isdigit():
                    bp = int(base_str[h])
                    line_bp += bp
        i += 1
    return int(line_bp)

def perc_calc(x,y):
    percent = int( ( float(x) / float(y) )*100)
    return percent

#create list to store bed file names in
bed_list = []

#go to directory with bed files
os.chdir(bed_dir)
for fl in os.listdir('.'):
    #JMSR008_index56_coding.bed
    if fl.endswith('.bed'):
        #add file names to bed list
        bed_list.append(fl)

#go to directory with bam files
os.chdir(bam_dir)

out_name = "{}_Specificity_Calcs_Output.txt".format(file_prefix)
fh_out = open(out_name, 'a')
fh_out.write("Sample"+'\t'+"Sam_bp_count"+'\t'+"cleaned_bp_count"+'\t'+"perc_Specificity"+'\n')
fh_out.close()

#iterate over files in bam directory
for fl in os.listdir('.'):
    #JMSR008_index56_sorted.bam
    if fl.endswith('_sorted.bam'):
        print "-----------------------------------------------------------------------"
        print time.asctime( time.localtime(time.time()) ), '\n'
        
        
        #split the file name
        flnames = fl.split('_')
        #reconstruct the sample prefix
        #JMSR008_index56
        pref = flnames[0]+'_'+flnames[1]

        #find the matching bed file
        for bed in bed_list:
            bed_split = bed.split('_')
            bed_name = bed_split[0]+'_'+bed_split[1]
            if bed_name == pref:
                bed_file = bed
        bed_path = bed_dir+'/'+bed_file

        print "Found {0} and now converting to sam format using {1}...".format(fl, bed_file)

        #convert bam file into MASSIVE sam file (sometimes >1GB)
        print "system call: samtools view -L {0} {1} > {2}_bed.sam".format(bed_path, fl, pref), '\n'
        sam_string = "samtools view -L {0} {1} > {2}_bed.sam".format(bed_path, fl, pref)
        proc_sam = sp.call(sam_string, shell=True)

        print "Finished writing {}_bed.sam".format(pref), '\n'
        temp_sam = "{}_bed.sam".format(pref)

        #start bp count for this file
        #**** sam count
        file_bp = int(0)
        print "1. Beginning base pair counts in {}_bed.sam...".format(pref)

        #iterate over lines and return bp per line, add to file_bp
        with open(temp_sam, 'r') as fh_sam:
            for line in fh_sam:
                line = line.strip()
                line = line.split()
                base_str = line[5]
                
                line_bp = sam_base_counter(base_str)
                file_bp += line_bp
                
            print '\t', "- There are {0} base pairs in {1}_bed.sam.".format(file_bp, pref), '\n'

        print "2. Changing directory to find matching cleaned reads files."
        os.chdir(clean_dir)
        total_clean_bp = int(0)
        
        for fl in os.listdir('.'):
            if fl.startswith(pref):
                names = fl.split('_')
                truename = names[0]+'_'+names[1]
                if truename == pref:
                    print '\t', "Beginning base pair counts of {}...".format(fl)
                    fh_temp = open(fl, 'r')
                    lines = fh_temp.readlines()
                    
                    i = int(1)
                    cleaned_file_bp = int(0)
                    total_lines = len(lines)

                    while i < total_lines:
                        line_bp = len(lines[i])
                        cleaned_file_bp += line_bp
                        i+=4
                    total_clean_bp += cleaned_file_bp
                    print '\t', "- {0} contains {1} base pairs.".format(fl, cleaned_file_bp), '\n'

        print '\t', "These files for sample {0} total {1} base pairs!".format(truename, total_clean_bp), '\n'

        temp_perc = perc_calc(file_bp,total_clean_bp)
        print "{0} has {1}% specificity.".format(pref, temp_perc), '\n'
        
        os.chdir(bam_dir)
        fh_out = open(out_name, 'a')
        fh_out.write(pref+'\t'+"{}".format(file_bp)+'\t'+"{}".format(total_clean_bp)+'\t'+"{}".format(temp_perc)+'\n')
        fh_out.close()
        for fl in os.listdir('.'):
            if fl.endswith('_bed.sam'):
                os.remove(fl)
                print "Removing {}...".format(fl), '\n'
        print time.asctime( time.localtime(time.time()) )
        print "-----------------------------------------------------------------------", '\n', '\n'

print "All finished, check '{0}' in '{1}' for awesome specificity results.".format(out_name, bam_dir)

fh_out.close()
