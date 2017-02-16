import csv
import os
import re
import numpy as np
import pandas as pd
import math

from Bio import SeqIO
from Bio.Seq import Seq

def segment_translate(raw_data_path, clean_data_path, species, outname, min_len):
    #find all files in specified directories
    dict_files = []
    for file in os.listdir(os.path.join(raw_data_path,species,"dict")):
        if file.endswith(".txt"):
            dict_files.append(file)

    seq_files = []
    for file in os.listdir(os.path.join(raw_data_path,species,"seq")):
        if file.endswith(".fa"):
            seq_files.append(file)

    clean_dir = os.path.join(clean_data_path, species)

    if not os.path.exists(clean_dir):
        os.makedirs(clean_dir)

    with open(os.path.join(clean_dir, outname + ".csv"), 'wb') as outfile:
        # Assuming that we process from 3' to 5' and UTR regions are continuous
        # once introns are excluded
        csvwriter = csv.writer(outfile, delimiter=',')
        n_bp = 0

        #Iterate through seq-annotation pairs to generate labeled training set
        for key_file in dict_files:
            search_string = re.sub("\..+$", '', file)
            ct = 0
            for test_file in seq_files:
                if re.sub("\..+$", '', test_file) == search_string:
                    seq_file = seq_files[ct]
                    break
                ct += 1

            #Convert csv file contain coordinates of gene elements to pd dataframe
            seq_index = pd.read_csv(os.path.join(raw_data_path,species,"dict", key_file))
            seq_index.columns = seq_index.columns.str.replace("\\s+","")
            #Filter for protein-coding Transcripts
            seq_index = seq_index[~seq_index['ProteinID'].isnull()]
            print(len(seq_index.index))
            #Remove overlapping proteins. Keep protein that starts first
            seq_agg = seq_index.groupby(['ProteinID']).agg({'ExonChrStart(bp)': min,
                                                            'ExonChrEnd(bp)': max,
                                                            'Strand' : min})
            ct = 0
            len_index = len(seq_agg.index)
            while ct < len_index:
                ind = seq_agg.index[ct]
                strand = seq_agg['Strand'][ind]
                start = seq_agg['ExonChrStart(bp)'][ind]
                end = seq_agg['ExonChrEnd(bp)'][ind]
                seq_agg = seq_agg[(seq_agg['Strand'][ind] != strand) | (seq_agg.index.isin([ind])) | (seq_agg['ExonChrStart(bp)'] >= end) | (seq_agg['ExonChrEnd(bp)'] < start)]
                len_index = len(seq_agg.index)
                ct = 1 + [item for item in range(len(seq_agg.index)) if seq_agg.index[item] == ind][0]

            print(pd.DataFrame.head(seq_agg))
            print(len(seq_agg.index))
            #print(pd.DataFrame.head(seq_index_filt))
            #Keep only names that appear in filtered set

            seq_index = seq_index[seq_index['ProteinID'].isin(seq_agg.index)]
            # print(len(seq_index.index))
            seq_index['Dummy'] = seq_index['Strand']*seq_index['ExonChrStart(bp)']
            seq_index.sort_values(['ProteinID', 'Dummy'])
            print(len(seq_index.index))

            fasta_sequences = SeqIO.parse(open(os.path.join(raw_data_path,species,"seq", seq_file)),'fasta')

            for fasta in fasta_sequences:
                #Initialize string and indexing parameters
                pt_prev = ''
                seq = ''
                len_3 = 0
                len_5 = 0
                last_ind = seq_index.index[-1]
                first_ind = seq_index.index[0]
                for i in seq_index.index:
                    if i < last_ind:
                        strand_id = seq_index['Strand'][i]
                        pt = seq_index['ProteinID'][i]
                        ex_s = seq_index['ExonChrStart(bp)'][i]-1
                        ex_e = seq_index['ExonChrEnd(bp)'][i]

                    if (i != first_ind and pt != pt_prev) or (i == last_ind):
                        len_3 = int(len_3)
                        len_5 = int(len_5)

                        if len_5 > min_len:
                            ind_5 = int(math.floor(len_5 / 3)*3)
                            csvwriter.writerow([str(pt_prev), seq[0:ind_5].translate(), 0])
                        if len_3 > min_len:
                            ind_3 = int(math.floor(len_3 / 3) * 3)
                            csvwriter.writerow([str(pt_prev), seq[-ind_3:].translate(), 0])
                        seq_raw = seq[len_5:len(seq)-len_3]
                        bound = int(math.floor(len(seq_raw)/3)*3)
                        if bound > min_len:
                            csvwriter.writerow([str(pt_prev), seq_raw[0:bound].translate(), 1])

                        n_bp += len(seq)
                        if i == last_ind:
                            break

                        seq = ''
                        len_3 = 0
                        len_5 = 0

                    if i < last_ind:
                        seq_new = Seq(str(fasta.seq[ex_s:ex_e]))
                        if strand_id == -1:
                            seq_new = seq_new.reverse_complement()

                        seq += seq_new
                        if not np.isnan(seq_index["3'UTRStart"][i]) or not np.isnan(seq_index["3'UTREnd"][i]):
                            len_3 += seq_index["3'UTREnd"][i] - seq_index["3'UTRStart"][i] + 1

                        if not np.isnan(seq_index["5'UTRStart"][i]) or not np.isnan(seq_index["5'UTREnd"][i]):
                            len_5 += seq_index["5'UTREnd"][i] - seq_index["5'UTRStart"][i] + 1

                        pt_prev = pt

        training_set = pd.read_csv(os.path.join(clean_dir, outname + ".csv"))
        return(training_set)

if __name__ == "__main__":
    # set path to raw data
    raw_data_path = os.path.join(os.getcwd(), "..", "..", "..", "sequence_data", "raw")
    clean_data_path = os.path.join(os.getcwd(), "..", "..", "..", "sequence_data", "clean")
    species = "human"
    outname = "testfile"
    min_len = 9

    call = segment_translate(raw_data_path, clean_data_path, species, outname, min_len)
    print(pd.DataFrame.head(call))