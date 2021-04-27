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

data_file_names = !ls
files = []
for i in data_file_names:
    if i[-12:] == "MERGED.fastq":
        fx = folder_name + "/" + i
        files.append(fx)

#Write out text file
save_name = folder_name + "/Primer_corrections_V4R_thresh3.txt"
save_file = open(save_name, "w")
header = ("Sample",'\t',"total reads",'\t',"no mismatch",'\t',"total corrected",'\t',"1",'\t',"2",'\t',"3",'\t',"4",'\t',"5",'\t',"6",'\t',"7",'\t',"8",'\t',"9",'\t',"10",'\t',"11",'\t',"12",'\t',"13",'\t',"14",'\t',"15",'\t',"16",'\t',"17",'\t',"18",'\t',"19",'\t',"20",'\t',"1_A",'\t',"2_A",'\t',"3_A",'\t',"4_A",'\t',"5_A",'\t',"6_A",'\t',"7_A",'\t',"8_A",'\t',"9_A",'\t',"10_A",'\t',"11_A",'\t',"12_A",'\t',"13_A",'\t',"14_A",'\t',"15_A",'\t',"16_A",'\t',"17_A",'\t',"18_A",'\t',"19_A",'\t',"20_A",'\t',"1_T",'\t',"2_T",'\t',"3_T",'\t',"4_T",'\t',"5_T",'\t',"6_T",'\t',"7_T",'\t',"8_T",'\t',"9_T",'\t',"10_T",'\t',"11_T",'\t',"12_T",'\t',"13_T",'\t',"14_T",'\t',"15_T",'\t',"16_T",'\t',"17_T",'\t',"18_T",'\t',"19_T",'\t',"20_T",'\t',"1_G",'\t',"2_G",'\t',"3_G",'\t',"4_G",'\t',"5_G",'\t',"6_G",'\t',"7_G",'\t',"8_G",'\t',"9_G",'\t',"10_G",'\t',"11_G",'\t',"12_G",'\t',"13_G",'\t',"14_G",'\t',"15_G",'\t',"16_G",'\t',"17_G",'\t',"18_G",'\t',"19_G",'\t',"20_G",'\t',"1_C",'\t',"2_C",'\t',"3_C",'\t',"4_C",'\t',"5_C",'\t',"6_C",'\t',"7_C",'\t',"8_C",'\t',"9_C",'\t',"10_C",'\t',"11_C",'\t',"12_C",'\t',"13_C",'\t',"14_C",'\t',"15_C",'\t',"16_C",'\t',"17_C",'\t',"18_C",'\t',"19_C",'\t',"20_C",'\t','\n')
save_file.write(''.join(map(str, header)))
new_tab = '\t'
new_line = '\n'


for q in files:
    GG_filename = q

    #Count number of records in the file
    count = 0
    for record in SeqIO.parse(GG_filename, "fastq"):
        count += 1
    GG_rec = count
    print("There were " + str(GG_rec) + " records in the sample")

    #Finding V4 forward primer sites
    RP = []
    R_match_count = 0
    V4_R_pos = []
    For_primer = []
    mismatch_rname = []
    mismatch_seq = []
    for record in SeqIO.parse(GG_filename, "fastq"):
        rname = record.description
        rc_record = record.reverse_complement()
        #R2_primer = "GGACTAC[ACT][ACG]GGGT[AT]TCTAAT" 
        V4_R = rc_record.seq.find("GGACTACCAGGGTATCTAAT")
        if V4_R == -1:
            V4_R_pos.append(V4_R)
            mismatch_rname.append(rname)
            mismatch_seq.append(rc_record.seq)
        else:
            V4_R_pos.append(V4_R)
            R_match_count +=1
    print("Found " + str(R_match_count) + " total records containing a perfect match to the V4_R primer out of " + str(GG_rec) + " total records")


    #find positions of corrections
    corr_pos = []
    corr_base = []
    corr_seq = []
    corr_name = []
    large_num_corr = []
    counter = 0
    for i in mismatch_seq:
        c_pos_item = []
        c_base = []
        count = 0
        for j, item in enumerate(i[0:20]):
            if item != "N":
                if item == "GGACTACCAGGGTATCTAAT"[j]:
                    pass
                else:
                    c_pos_item.append(j)
                    c_base.append(item)
                    count += 1
        if count > 0 and count < 4:
            corr_pos.append(c_pos_item)
            corr_base.append(c_base)
            corr_seq.append(i)
            corr_name.append(mismatch_rname[counter])
        else: 
            #print "> 3 corrections"
            large_num_corr.append(i)
        counter += 1
        
    import itertools
    chain_pos = itertools.chain(*corr_pos)
    chain_base = itertools.chain(*corr_base)

    corr_positions = list(chain_pos)
    corr_nucleotide = list(chain_base)
    
    try:
        print("Found " + str(len(corr_pos)) + "/" + str(GG_rec) +" total records containing a likely true mismatch (" + str((len(corr_pos))/float(GG_rec)*100) + "%)")
    except ZeroDivisionError:
        print "none found"
    save_file.write(q.split("/")[-1])
    save_file.write(new_tab)
    save_file.write(str(GG_rec))
    save_file.write(new_tab)
    save_file.write(str(R_match_count))
    save_file.write(new_tab)
    save_file.write(str(len(corr_pos)))
    save_file.write(new_tab)
    
    # plt.clf()
#     if len(corr_pos) > 0:
#         #fig_dir = folder_name + "/Figs"
#         #os.chdir(fig_dir)
#         #Make position histogram - where are these sites?
#         fig = plt.figure()
#         plt.title("Distribution of primer corrections")
#         ax = fig.add_subplot(111)
#         x = corr_positions
#         numBins = 19
#         ax.hist(x,numBins,color='green')
#         ax.set_ylabel("Number of corrections")
#         ax.set_xlabel('Position')
#         #ticks = range(1, (numBins+1))
#         #plt.xticks(ticks)
#         plt_name = GG_filename[:-6] + "_corrections.jpg"
#         plt.savefig(plt_name)
#         
    #Write out file with correction data
    corr_by_pos = []
    A_list = []
    T_list = []
    G_list = []
    C_list = []
    for i in range(20):
        count = 0
        count_A = 0
        count_T = 0
        count_G = 0
        count_C = 0
        for j, item in enumerate(corr_positions):
            if item == i:
                count += 1
                if corr_nucleotide[j] == 'A':
                    count_A +=1
                elif corr_nucleotide[j] == 'T':
                    count_T +=1
                elif corr_nucleotide[j] == 'G':
                    count_G +=1
                elif corr_nucleotide[j] == 'C':
                    count_C +=1
        corr_by_pos.append(count)
        A_list.append(count_A)
        T_list.append(count_T)
        G_list.append(count_G)
        C_list.append(count_C)
                
    for i, item in enumerate(corr_by_pos):
        save_file.write(str(item))
        save_file.write(new_tab)
    for i, item in enumerate(A_list):
        save_file.write(str(item))
        save_file.write(new_tab)
    for i, item in enumerate(T_list):
        save_file.write(str(item))
        save_file.write(new_tab)
    for i, item in enumerate(G_list):
        save_file.write(str(item))
        save_file.write(new_tab)
    for i, item in enumerate(C_list):
        save_file.write(str(item))
        save_file.write(new_tab)    
    save_file.write(new_line)   

save_file.close()


