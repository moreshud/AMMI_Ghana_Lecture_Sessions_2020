#!/usr/bin/env python

# class fastq(object):
#     def __init__(self,filename):
#         self.filename = filename
#         self.__sequences = {}
        
#     def parse_file(self):
#         with open(self.filename, 'r') as f:
#             content = f.readlines()

#             # Recreate content without lines that start with @ and +
#             content = [line for line in content if not line[0] in '@+']

#             # Now the lines you want are alternating, so you can make a dict
#             # from key/value pairs of lists content[0::2] and content[1::2]
#             data = dict(zip(content[0::2], content[1::2]))

#         return data


# file_name = "sb008.fastz"
# fastq(file_name).parse_file()



# with open('sb008.fastz') as f:
#     lines=f.readlines()
# head=[item[:-1] for item in lines[::4]] #get rid of '\n'
# read=[item[:-1] for item in lines[1::4]]
# qual=[item[:-1] for item in lines[3::4]]
# dict(zip(read, qual))


#imports packages
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

# reading the file
path = './sb008.fastz'
data = open(path,'r')
my_data = {}
my_data_len = {}
for line in data.readlines():
    line = line.rstrip('\n')
    if line[0]=='>':
        key = line[1:]
    else:
        good_line = ''
        for char in line.replace("N", ""):
            if char != '-':
                good_line += char
        my_data_len[key] = len(good_line)
        my_data[key] = good_line
#print(my_data)
data.close

# no of comments, sequence max leng, 4-3-1 kmers, co-occurence matrix, bag of words, similarity check,

print('THE NUMBWER OF COMMENTS: ',len(my_data.keys()))
print('The size of the seq with max length:',max(my_data_len.values()))

# getting the kmer function
def getKmers(sequence, size=6):
    return [sequence[x:x+size] for x in range(len(sequence) - size + 1)]

for seq in my_data.values():
    # print(getKmers(seq))
    break

# Bag of words function - row vector with key - sequence apperance pattern 
def BoG(seq,k):
    kmers_list = getKmers(seq,k)
    bog = {}
    for w in kmers_list:
        if w in bog.keys():
            bog[w] += 1
        else:
            bog[w] = 1
    return bog

for value in my_data.values():
    # print(BoG(value, 5))
    break
# Co-occurence matrix with key rows and unique kmers pattern as columns
def co_matrix(k):
    cols = []
    for key in my_data.keys():
        cols += getKmers(my_data[key],k)
    cols = set(cols)
    co_mtx = []
    for value in my_data.values():
        seq_bog = BoG(value,k)
        row = []
        for col in cols:
            if col in seq_bog.keys():
                row.append(seq_bog[col])
            else:
                row.append(0)
        co_mtx.append(row)
    return np.array(co_mtx)

# print(co_matrix(k=6))

co_maxt = co_matrix(k=6)

# Similarity functions using cosine
from itertools import combinations
ind = list(combinations(np.arange(8),2))
print(ind)

def similarity(comtx):
    sim = {}
    names = list(my_data.keys())
    for indx in ind:
        i, i_ = indx
        A,B = comtx[i],comtx[i_]
        sim_lavue = np.dot(A,B)/(np.linalg.norm(A)*np.linalg.norm(B))
        key = names[i] + ' vs ' + names[i_]
        sim[key] = sim_lavue
    return sim
print(similarity(co_maxt))

# Kmeans clustering of 2 species
# - getting their values  - you can work with index 0, 1----8
sim_val = similarity(co_maxt)
species = ["mouse", "rat"]  # key ordering matters because h,c,m,r,d,c,a,e


def generate_table(species_list):
    data_table = np.zeros((6,len(species)))

    for i in range(len(species)):
        idx = 0
        for key in sim_val: 
            if (species[i] in key) and (species[0] + " vs " + species[1] not in key): # ordering based keys matters say mouse vs chimps
                data_table[idx, i] = sim_val[key]
                idx += 1    
    return data_table

X = generate_table(species)
print(X)


# Kmeans algorithm
from sklearn.cluster import KMeans

model = KMeans(n_clusters=2)
model.fit(X)
predictions = model.predict(X)
print(predictions)

# clustering output and visualization
from collections import defaultdict

groups = list(my_data.keys())
for s in species:
    groups.remove(s)

clusters = defaultdict(list)
for index, value in enumerate(groups):
    clusters[species[predictions[index]]].append(value)

print(list(clusters.items()))

plt.scatter(X[:, 0], X[:, 1], c=predictions, s=100, cmap='plasma')
for i, group in enumerate(groups):
    plt.annotate(group, (X[:, 0][i], X[:, 1][i]), textcoords="offset points", xytext=(0,10),ha='left', va='center_baseline')

centers = model.cluster_centers_
plt.scatter(centers[:, 0], centers[:, 1], c='red', s=200, alpha=1)
for i, specie in enumerate(species):
    plt.annotate(specie, (centers[:, 0][i], centers[:, 1][i]), textcoords="offset points", xytext=(0,10),ha='right', va='center')
# plt.show()
# quick one with experimenting with different kmers sequence

# Blasting using Biopython - becareful for kmers iteration not to overwite or choiceful kmers approach to proceed
def parse_file_write(data):
    file = open("cleaned_sb008_data.fasta", "w+")

    for key, value in data.items():
        file.write(">" + key + "\n" + value + "\n")
    file.close()

parse_file_write(my_data)

# blasting starts
from Bio.Blast import NCBIWWW
from Bio import SeqIO, SearchIO
from Bio.SeqUtils import GC
# taking more than needed
# fasta_string = open("cleaned_sb008_data.fasta")
# result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string.read())
# blast_result = open("blasted_sb008_data_1.xml", "w")
# blast_result.write(result_handle.read())
# fasta_string.close()
# blast_result.close()
# result_handle.close()

# this approach works
# record_dict = SeqIO.to_dict(SeqIO.parse("cleaned_sb008_data.fasta", "fasta"))
# result_handle = NCBIWWW.qblast("blastn", "nt", record_dict)
# blast_result = open("blasted_sb008_data_2.xml", "w")
# blast_result.write(result_handle.read())
# record_dict
# blast_result.close()
# result_handle.close()
# print("Sequence length (bp)", len(record_dict))
# print("GC content (%)", GC(record_dict))
# blast_qresult = SearchIO.read(result_handle, "blast-xml")
# print(blast_qresult)

# print([hit.description for hit in blast_qresult[:5]])
# first_hit = blast_qresult[0]
# print(first_hit.description)
# first_hsp = first_hit[0]
# print(first_hsp.evalue, first_hsp.bitscore)
# print(first_hsp.aln)


# motif assignments -- yet to be test runned
from Bio import motifs
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, generic_dna
from Bio import AlignIO
import os

# alignment - archived
# input_file = './cleaned_sb008_data.fasta'
# records = SeqIO.parse(input_file, 'fasta')
# records = list(records) # make a copy, otherwise our generator
#                         # is exhausted after calculating maxlen
# maxlen = max(len(record.seq) for record in records)

# # pad sequences so that they all have the same length
# for record in records:
#     if len(record.seq) != maxlen:
#         sequence = str(record.seq).ljust(maxlen, 'N')
#         record.seq = Seq(sequence)
# assert all(len(record.seq) == maxlen for record in records)

# # write to temporary file and do alignment
# output_file = '{}_padded.fasta'.format(os.path.splitext(input_file)[0])
# with open(output_file, 'w') as f:
#     SeqIO.write(records, f, 'fasta')
# alignments = AlignIO.read(output_file, "fasta")
# instances = [Seq(str(align.seq).upper(),  IUPAC.IUPACAmbiguousDNA()) for align in alignments]
# print(len(alignments[0].seq.upper()) == len(alignments[1].seq.upper()))
# print(str(alignments[0].seq)[0])


# creating motifs needed instances to be aligned
# sequences = [value for value in my_data.values()]
# aligned_sequences = [seq.ljust(max(my_data_len.values()), "N") for seq in sequences] 
# instances = [aligned.upper() for aligned in aligned_sequences]
# # print(instances)
# m = motifs.create(instances, IUPAC.ambiguous_dna)
# print(m)
# print(len(m))
# print(m.counts)
print(m.alphabet)
print(len(m.consensus))
# print(len(m.anticonsensus))
# print(len(m.degenerate_consensus))
# # print(m.pwm)
# weblogo(m, "mymotifs.png")
# m.weblogo("mymotifs.png") #Sequence logo

# # # jaspar and trasfac

# # # exercise 4 ncbi site
