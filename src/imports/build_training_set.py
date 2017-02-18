import csv
import os

import re
import numpy as np
import pandas as pd
import math
import itertools

from Bio.Data import CodonTable
from Bio import SeqIO
from Bio.Seq import Seq

##Parse coding and noncoding nc sequence from fasta files. Translate.
# Sequences are truncated to nearest multiple of 3 as necessary
# Saves resulting dataframe to file and returns for immediate use in subsequent functions
def segment_translate(raw_data_path, clean_data_path, species, outname, min_len):
    #"raw_data_path": path to directory for raw fasta files
    #"clean_data_path": path to directory for cleaned data (training sets, keys, etc.)
    #"species": specifies subfolder for raw and clean data (sorted by species)
    #"outname": name of set to be saved
    #"min_len": minimum length for parsed training sequnence. Should be kept >= 3
    #find all files in specified directories
    dict_files = []
    seq_files = []
    for file in os.listdir(os.path.join(raw_data_path,species)):
        #Thes directories should have ONLY relevant fasta and dict files
        if file.endswith(".txt"):
            dict_files.append(file)
        if file.endswith(".fa"):
            seq_files.append(file)

    clean_dir = os.path.join(clean_data_path, 'training', 'char', species)

    if not os.path.exists(clean_dir):
        os.makedirs(clean_dir)

    with open(os.path.join(clean_dir, outname + ".csv"), 'wb') as outfile:
        # Assuming that we process from 3' to 5' and UTR regions are continuous once introns are excluded
        csvwriter = csv.writer(outfile, delimiter=',')
        #tracks total bp's saved to final set
        n_bp = 0
        print("Parsing sequences from raw fasta files")
        #Iterate through seq-annotation pairs to generate labeled training set
        for key_file in dict_files:
            search_string = re.sub("\..+$", '', file)
            ct = 0
            seq_file = ''
            for test_file in seq_files:
                if re.sub("\..+$", '', test_file) == search_string:
                    seq_file = seq_files[ct]
                    break
                ct += 1
            print("Processing key file: " + key_file + ", seq file: " + seq_file)
            #Convert csv file contain coordinates of gene elements to pd dataframe
            seq_index = pd.read_csv(os.path.join(raw_data_path,species, key_file))
            seq_index.columns = seq_index.columns.str.replace("\\s+","")
            #Filter for protein-coding Transcripts
            seq_index = seq_index[~seq_index['ProteinID'].isnull()]

            #Remove overlapping proteins. Keep protein that starts first.
            seq_agg = seq_index.groupby(['ProteinID']).agg({'ExonChrStart(bp)': min,
                                                            'ExonChrEnd(bp)': max,
                                                            'Strand' : min})
            seq_agg['DummyStart'] = seq_agg['Strand'] * seq_agg['ExonChrStart(bp)']
            seq_agg['DummyStop']  = seq_agg['Strand'] * seq_agg['ExonChrEnd(bp)']

            seq_agg = seq_agg.sort_values(["DummyStart"])

            ct = 0
            len_index = len(seq_agg.index)
            while ct < len_index - 1:
                ind = seq_agg.index[ct]
                stop = seq_agg['DummyStop'][ind]
                increment = 1
                next_start = seq_agg['DummyStart'][seq_agg.index[ct+increment]]
                drop_list = []
                #Drop all entries that overlap
                while next_start < stop and (ct + increment) < len_index - 1:
                    drop_list.append(seq_agg.index[ct+increment])
                    increment += 1
                    next_start = seq_agg['DummyStart'][seq_agg.index[ct + increment]]

                seq_agg = seq_agg[~seq_agg.index.isin(drop_list)]
                #Update length
                len_index = len(seq_agg.index)
                #find new position of current index in list and increment by 1
                ct = 1 + [item for item in range(len(seq_agg.index)) if seq_agg.index[item] == ind][0]

            #Keep only names that appear in filtered set
            seq_index = seq_index[seq_index['ProteinID'].isin(seq_agg.index)]

            seq_index['Dummy'] = seq_index['Strand']*seq_index['ExonChrStart(bp)']
            seq_index = seq_index.sort_values(['ProteinID', 'Dummy'])

            #Open Seq File. Framework can handle file containing multiple sequences
            fasta_sequences = SeqIO.parse(open(os.path.join(raw_data_path,species, seq_file)),'fasta')

            for fasta in fasta_sequences:
                #Initialize string and indexing parameters
                pt_prev = ''
                seq = ''
                len_3 = 0
                len_5 = 0
                last_ind = seq_index.index[-1]
                first_ind = seq_index.index[0]
                for i in seq_index.index:

                    if i != last_ind:
                        strand_id = seq_index['Strand'][i]
                        pt = seq_index['ProteinID'][i]
                        #Genome indexing starts at 1 in dict, so we must adjust
                        ex_s = seq_index['ExonChrStart(bp)'][i]-1
                        #Dict coordinates are inclusive, so no need to adjust end
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

                    if i != last_ind:
                        seq_new = Seq(str(fasta.seq[ex_s:ex_e]))
                        if strand_id == -1:
                            seq_new = seq_new.reverse_complement()

                        seq += seq_new
                        if not np.isnan(seq_index["3'UTRStart"][i]) and not np.isnan(seq_index["3'UTREnd"][i]):
                            len_3 += seq_index["3'UTREnd"][i] - seq_index["3'UTRStart"][i] + 1

                        if not np.isnan(seq_index["5'UTRStart"][i]) and not np.isnan(seq_index["5'UTREnd"][i]):
                            len_5 += seq_index["5'UTREnd"][i] - seq_index["5'UTRStart"][i] + 1

                        pt_prev = pt

        training_set = pd.read_csv(os.path.join(clean_dir, outname + ".csv"), header=None)

        return(training_set)

#generate efficient AA encoding scheme based on "Amino Acid Encoding Schemes for Machine Learning Methods" results (2011)
def make_aa_vector_key(datapath):

    #Enumerate all possible 3nt nucleic acid combos
    nc_acids = ["".join(seq) for seq in itertools.product("ATCG", repeat=3)]

    aa_nc_key = pd.DataFrame(np.nan, index=range(len(nc_acids)), columns=['AA', 'NC'], dtype='S3')

    #Make NC-AA key
    for i in xrange(len(nc_acids)):
        aa_nc_key.set_value(i, 'AA', str(Seq(nc_acids[i]).translate()))
        aa_nc_key.set_value(i, 'NC', str(nc_acids[i]))

    #Map AA and NC characters to digits for indexing
    nc_num_dict = {'A': 0, 'T': 1, 'C': 2, 'G': 3}

    aa_vec = list(set(aa_nc_key['AA']))
    aa_num_dict = {key: value for (key,value) in zip(aa_vec, range(len(aa_vec)))}

    #Create array to store transition info
    aa_transfers = np.zeros((4*4,len(aa_vec)))

    for i in aa_nc_key.index:
        aa_id = aa_num_dict[aa_nc_key['AA'][i]]
        aa_col = aa_transfers[:,aa_id]
        nc_str = aa_nc_key['NC'][i]

        fr = nc_num_dict[nc_str[0]]
        for j in xrange(1,len(nc_str)):
            to = nc_num_dict[nc_str[j]]
            aa_col[4*to+fr] = 1
            fr = to
        aa_transfers[:, aa_id] = aa_col
    aa_mat = aa_transfers
    aa_transfers = pd.DataFrame(aa_transfers)
    aa_transfers.columns = aa_vec

    aa_transfers.to_csv(os.path.join(datapath, 'clean', 'keys', "aa_vectors.csv"))
    return(aa_mat, aa_num_dict)

#Code to pad AA sequences with "random" seq to bring them up to standardized length for CNN input
#Path to clean data dir

def clean_pad_aa(char_aa_set, clean_dir, dataname, seq_len, min_len, fragment_size, species):
    #char_aa_set = pd.read_csv(os.path.join(clean_dir, "training", "char", species, dataname + ".csv"), header=None)

    char_aa_set.columns = ["ProteinID", "AA","ID"]
    n_rows = len(char_aa_set.index)

    seq_concat = ''.join(list(char_aa_set["AA"]))
    seq_concat = re.sub("X","",seq_concat)
    seg = 0
    segments = [0]
    while seg < len(seq_concat) - 1:
        inc = min(seg + np.random.randint(1,fragment_size+1), len(seq_concat) - 1)
        seg += inc
        segments.append(seg)
    #shuffle sequence
    seg_index = np.arange(len(segments)-1)
    np.random.shuffle(seg_index)
    seq_shuffle = ''

    for fragment in seg_index:
        start = segments[fragment]
        end = segments[fragment+1]
        seq_shuffle += seq_concat[start:end]

    #Make copy of char_aa_set. Will be updated as we go
    char_aa_pad = char_aa_set
    #Keep track of rows that need to be dropped. Will do this at the end
    drop_list= []
    for ind in char_aa_set.index:
        aa_string = char_aa_pad["AA"][ind]
        aa_id = char_aa_pad["ID"][ind]

        if aa_id == 0:
            split_aa = re.split('X+', aa_string)
            outstring = ''
            max_length = 0
            for st in split_aa:
                if len(st) > max_length:
                    outstring = st
                    max_length = len(st)
            aa_string = outstring

        elif aa_id == 1:
            if 'X' in aa_string:
                drop_list.append(ind)
                continue
        else:
            print('Warning: Unidentified Label')
            print(aa_string)
            drop_list.append(ind)
            continue
        n = len(aa_string)

        #Filter out sequences of unsuitable length. Could clip overly long sequences
        if n > seq_len or n < min_len:
            drop_list.append(ind)
            continue
        #Randomly assign start-end coordinates for actual seq
        start = np.random.randint(0, seq_len-n + 1) #inclusive
        end = start + n #exclusive
        #Determine size of padding sequences
        n_left = start
        left_start = np.random.randint(0, len(seq_shuffle)-n_left)
        left_pad = seq_shuffle[left_start:left_start+n_left]
        n_right = seq_len - end
        right_start = np.random.randint(0, len(seq_shuffle) - n_right)
        right_pad = seq_shuffle[right_start:right_start + n_right]

        char_aa_pad.set_value(ind, 'AA', left_pad + aa_string + right_pad)

    char_aa_pad = char_aa_pad[~char_aa_pad.index.isin(drop_list)]
    char_aa_pad.to_csv(os.path.join(clean_dir, "training", "char", species, dataname + "_padded.csv"))

    print(str(len(char_aa_pad.index)) + " kept out of " + str(len(char_aa_set.index)) + ' (' + str(
        sum(char_aa_pad["ID"])) + " coding)")
    return(char_aa_pad)


def vectorize_aa(aa_vector_key, aa_num_dict, char_aa_padded, clean_dir,seq_len):
    #3D Array to store vectorized AA sequences
    aa_vector_set = np.zeros((16,seq_len,len(char_aa_padded.index)))
    label_vec = np.zeros((len(char_aa_padded.index),1))
    ct = 0
    for ind in char_aa_padded.index:
        label_vec[ct] = char_aa_padded["ID"][ind]
        aa_string = char_aa_padded["AA"][ind]
        aa_num = []
        for i in xrange(len(aa_string)):
            aa_num.append(aa_num_dict[aa_string[i]])
        aa_vector_set[:,:,ct] = aa_vector_key[:,aa_num]
        ct += 1
    return(aa_vector_set,label_vec)

if __name__ == "__main__":
    # set path to raw data
    raw_data_path = os.path.join(os.getcwd(), "..", "..", "..", "sequence_data", "raw")
    clean_data_path = os.path.join(os.getcwd(), "..", "..", "..", "sequence_data", "clean")
    species = "human"
    outname = "testfile"
    min_len = 9

    ##Extract set of coding and noncoding sequences from genomic DNA
    char_aa_set = segment_translate(raw_data_path, clean_data_path, species, outname, min_len)
    #Make Amino Acid Key
    datapath = os.path.join(os.getcwd(), "..", "..", "..", "sequence_data")
    (aa_vector_key, aa_num_dict) = make_aa_vector_key(datapath)

    # Set standardized seq length
    seq_len = 300
    # minimum allowed length for sequence
    min_len = 60
    # Set max "fragment" size for random seq compilation
    fragment_size = 10
    #Standardize Sequence Length. Remove sequences with ambiguous AA's ("X")
    char_aa_padded = clean_pad_aa(char_aa_set, clean_data_path, outname, seq_len, min_len, fragment_size, species)

    #Generate Final Vectorized Training Set and LAbel Vector
    (aa_training_set, label_vec) = vectorize_aa(aa_vector_key, aa_num_dict, char_aa_padded, clean_data_path, seq_len)
