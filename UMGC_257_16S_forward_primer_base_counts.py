from Bio import SeqIO
import numpy as np
#import matplotlib.pyplot as plt
from Bio.Restriction import *
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
import itertools
from Bio.SeqRecord import SeqRecord
import re
import os
import sys

path = "~/Seq/UMGC_257"
os.chdir(os.path.expanduser(path))
folder_name = os.getcwd()

R1_primer = "GTGCCAGCAGCCGCGGTAA" #Read 1 primer sequence
R2_primer = "GGACTACCAGGGTATCTAAT" #Read 2 primer sequence
R1_adapter = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
R2_adapter =  "CTGTCTCTTATACACATCTGACGCTGCCGACGA"

os.chdir(folder_name) 
data_file_names = !ls

#Trim and stitch reads
for j in data_file_names:
    if j[-9:] == "001.fastq":
        if j[-11] == '1': #If this file corresponds to read 1
            c_file = j[:-6]
            R1 = c_file + "_trimmed.fastq"
            execute = "cutadapt -a " + R1_adapter + " " + j + " > " + R1
            os.system(execute) #Trim Illumina adapter from end of read
        elif j[-11] == '2': #If this file corresponds to read 2
            c_file = j[:-6]
            R2 = c_file + "_trimmed.fastq"
            execute = "cutadapt -a " + R2_adapter + " " + j + " > " + R2
            os.system(execute) #Trim Illumina adapter from end of read
            execute = "pandaseq -B -f " + R1 + " -r " + R2 + " -O 601 -l 150 -L 350 -F -w " + c_file + "_MERGED.fastq -U " + c_file + "_SINGLES.fastq -g " + c_file + "_log_1.txt"
            os.system(execute) #Stitch read 1 and read 2 with pandaseq

data_file_names = !ls
files = []
for i in data_file_names:
    if i[-12:] == "MERGED.fastq":
        fx = folder_name + "/" + i
        files.append(fx)

#Write out text file
save_name = folder_name + "/Primer_edits_V4F_base_counts.txt"
save_file = open(save_name, "w")
header = ("Sample",'\t',"total reads",'\t',"no mismatch",'\t',"total corrected",'\t',"N_1",'\t',"N_2",'\t',"N_3",'\t',"N_4",'\t',"N_5",'\t',"N_6",'\t',"N_7",'\t',"N_8",'\t',"N_9",'\t',"N_10",'\t',"N_11",'\t',"N_12",'\t',"N_13",'\t',"N_14",'\t',"N_15",'\t',"N_16",'\t',"N_17",'\t',"N_18",'\t',"N_19",'\t',"1_A",'\t',"2_A",'\t',"3_A",'\t',"4_A",'\t',"5_A",'\t',"6_A",'\t',"7_A",'\t',"8_A",'\t',"9_A",'\t',"10_A",'\t',"11_A",'\t',"12_A",'\t',"13_A",'\t',"14_A",'\t',"15_A",'\t',"16_A",'\t',"17_A",'\t',"18_A",'\t',"19_A",'\t',"1_T",'\t',"2_T",'\t',"3_T",'\t',"4_T",'\t',"5_T",'\t',"6_T",'\t',"7_T",'\t',"8_T",'\t',"9_T",'\t',"10_T",'\t',"11_T",'\t',"12_T",'\t',"13_T",'\t',"14_T",'\t',"15_T",'\t',"16_T",'\t',"17_T",'\t',"18_T",'\t',"19_T",'\t',"1_G",'\t',"2_G",'\t',"3_G",'\t',"4_G",'\t',"5_G",'\t',"6_G",'\t',"7_G",'\t',"8_G",'\t',"9_G",'\t',"10_G",'\t',"11_G",'\t',"12_G",'\t',"13_G",'\t',"14_G",'\t',"15_G",'\t',"16_G",'\t',"17_G",'\t',"18_G",'\t',"19_G",'\t',"1_C",'\t',"2_C",'\t',"3_C",'\t',"4_C",'\t',"5_C",'\t',"6_C",'\t',"7_C",'\t',"8_C",'\t',"9_C",'\t',"10_C",'\t',"11_C",'\t',"12_C",'\t',"13_C",'\t',"14_C",'\t',"15_C",'\t',"16_C",'\t',"17_C",'\t',"18_C",'\t',"19_C",'\n')
save_file.write(''.join(map(str, header)))
new_tab = '\t'
new_line = '\n'

##############################################################################

for q in files:
    GG_filename = q
    print(GG_filename)
    #Count number of records in the file
    #Finding V4 forward primer sites
    FP = []
    F_match_count = 0
    V4_F_pos = []
    For_primer = []
    mismatch_rname = []
    mismatch_seq = []
    count = 0
    for record in SeqIO.parse(GG_filename, "fastq"):
        rname = record.description
        V4_F = record.seq.find("GTGCCAGCAGCCGCGGTAA")
        For_primer.append(str(record.seq[0:19]))
        count += 1
        #if V4_F == -1:
        #    V4_F = record.seq.find("GTGCCAGCCGCCGCGGTAA")
        if V4_F == -1:
            V4_F_pos.append(V4_F)
            mismatch_rname.append(rname)
            mismatch_seq.append(record.seq)
        else:
            V4_F_pos.append(V4_F)
            F_match_count +=1
    GG_rec = count
    print("Found " + str(F_match_count) + " total records containing a perfect match to the V4_F primer out of " + str(GG_rec) + " total records")
    print("There were " + str(GG_rec) + " records in the sample")

    bases = ["A","G","T","C","N"]
    A_count = []
    G_count = []
    T_count = []
    C_count = []
    N_count = []
    for i in range(19):
        base_dict = dict.fromkeys(bases,0)
        #base_dict[item] += 1
        for j in For_primer:
            #print(i)
            base_dict[j[i]] += 1
        A_count.append(base_dict["A"])
        G_count.append(base_dict["G"])
        T_count.append(base_dict["T"])
        C_count.append(base_dict["C"])
        N_count.append(base_dict["N"])
            
    save_file.write(GG_filename)
    save_file.write(new_tab)
    save_file.write(str(GG_rec))
    save_file.write(new_tab)
    save_file.write(str(F_match_count))
    save_file.write(new_tab)
    save_file.write("NA")
    save_file.write(new_tab)
    for i, item in enumerate(N_count):
        save_file.write(str(item))
        save_file.write(new_tab)
    for i, item in enumerate(A_count):
        save_file.write(str(item))
        save_file.write(new_tab)
    for i, item in enumerate(T_count):
        save_file.write(str(item))
        save_file.write(new_tab)
    for i, item in enumerate(G_count):
        save_file.write(str(item))
        save_file.write(new_tab)
    for i, item in enumerate(C_count):
        save_file.write(str(item))
        save_file.write(new_tab)    
    save_file.write(new_line)   

save_file.close()