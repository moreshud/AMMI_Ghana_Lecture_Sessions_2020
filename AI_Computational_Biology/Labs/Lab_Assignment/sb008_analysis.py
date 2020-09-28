#!/usr/bin/env python
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from collections import defaultdict
from itertools import combinations

class SequenceAnalysis(object):
    def __init__(self,filepath):
        self.filepath = filepath
        self.data = {}
        self.sequence_length = []
        self.nucleotide_frequency = {}
        self.counter = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        self.size = 3
        
    def parse_file(self):
        with open(self.filepath, 'r') as content:
            for line in content.readlines():
                # remove the next line effect on each line
                line = line.rstrip("\n")
                
                # check for the comment and make key/value dict altenatively
                if line[0] == ">":
                    key = line[1:]
                else:
                    line = line.translate({ord(i): None for i in 'N-'}).upper() # remove all N and - in the main sequence
                    self.data[key] = line
                    self.sequence_length.append(len(line))
        content.close()
        # return self.data, self.sequence_length

    def nucleotide_statistics(self):
        for key, value in self.data.items():
            countDict ={'A': 0, 'T': 0, 'G': 0, 'C': 0}
            for nucleotide in value:
                countDict[nucleotide] += 1
                self.counter[nucleotide] += 1
            self.nucleotide_frequency[key] = countDict

    def getKmers(self, sequence, size):
        return [sequence[i:i+size] for i in range(len(sequence) - size + 1)]
    
    def bag_of_words(self, sequence, size):
        kmers_words = self.getKmers(sequence,size)
        container = {}
        for word in kmers_words:
            if word in container.keys():
                container[word] += 1
            else:
                container[word] = 1
        return container

    def co_occurence_matrix(self,size):
        matrix, decoded_sequences = [], []

        for val in self.data.values():
            decoded_sequences += self.getKmers(val,size)

        unique_sequences = set(decoded_sequences) # makes the column of the matrix

        for val in self.data.values():
            sequence_bog = self.bag_of_words(val, size)
            row = []

            for col in unique_sequences:
                if col in sequence_bog.keys():
                    row.append(sequence_bog[col])
                else:
                    row.append(0)
            matrix.append(row)
        return np.array(matrix)

    def similarity(self, size = None):
        size = size if size is not None else self.size
        sim = {}
        indices = list(combinations(np.arange(len(self.data.keys())), 2))
        names = list(self.data.keys())
        for idx in indices:
            i, j = idx
            A,B = self.co_occurence_matrix(size)[i], self.co_occurence_matrix(size)[j]
            value = np.dot(A,B)/(np.linalg.norm(A)*np.linalg.norm(B))
            key = names[i] + " vs " + names[j]
            sim[key] = value
        return sim

    def generate_table(self, cluster): # key ordering matters because h,c,m,r,d,c,a,e
        cluster_length = len(cluster)
        sim_value = self.similarity()
        data_table = np.zeros((6,cluster_length))

        for i in range(cluster_length):
            idx = 0
            for key in sim_value: 
                if (cluster[i] in key) and (cluster[0] + " vs " + cluster[1] not in key): # ordering based keys matters say mouse vs chimps
                    data_table[idx, i] = sim_value[key]
                    idx += 1    
        return data_table

    def generate_table(self, species):
        data_table = np.zeros((6,len(species)))
        sim_val = self.similarity()

        for i in range(len(species)):
            idx = 0
            for key in sim_val: 
                if (species[i] in key) and (species[0] + " vs " + species[1] not in key): # ordering based keys matters say mouse vs chimps
                    data_table[idx, i] = sim_val[key]
                    idx += 1    
        return data_table

    def parse_file_write(self):
        file = open("cleaned_sb008_data.fasta", "w+")

        for key, value in self.data.items():
            file.write(">" + key + "\n" + value + "\n")
        file.close()

  
path = './sb008.fastz'
seq_analysis = SequenceAnalysis(path)
seq_analysis.parse_file()
seq_analysis.nucleotide_statistics()

################################################################################################################################
# Rank the species that are more similar to human based on k-mer, similarity and clustering. Use the Sb008.fastz to 
################################################################################################################################

# Ques. 1a: compute statics including the size and frequency of a, c, g,t and n in each sequence
data = seq_analysis.data
length = seq_analysis.sequence_length
species_nucleotide_frequencies = seq_analysis.nucleotide_frequency
nucleotides_total = seq_analysis.counter 
num_comments, max_sequence, min_sequence = len(length), max(length), min(length)
specie_with_max_length = list(data.keys())[length.index(max_sequence)]
specie_with_min_length = list(data.keys())[length.index(min_sequence)]
# =====>
print("i. The data consist of genomic sequces of {} species".format(num_comments))
print("ii. The nucleotides frequecy in each species:")
for key, value in species_nucleotide_frequencies.items():
    print(key, " - ", value)
print("iii. The occurence of each nucleotide in the data is {} respectively".format(nucleotides_total))
print("iv. {} - {} has the maximum length of sequence while {} - {} has the minimum length of sequece".format(specie_with_max_length, max_sequence, specie_with_min_length, min_sequence))

# Ques 1b & c. compute k-mer similarity matrix based on k from 3 to 7 with ranking among species
kmers_data_list = []
# for i in range(3,8):
#     simDict = {"kmer": i}
#     simDict.update(seq_analysis.similarity(i))
#     kmers_data_list.append(simDict)
# kmers_data = kmers_data_list
# print(pd.DataFrame(kmers_data))

# Ques 1d. is their any agreement between the different k?
print("The ranking order is preserved with different k. However, the ranking value decreases as k increases")

# Ques 1e. cluster the species - using the default 3-kmers
species = ["human", "elephant"]
X = seq_analysis.generate_table(species) # cluster table values

# Kmeans algorithm
from sklearn.cluster import KMeans

model = KMeans(n_clusters=2)
model.fit(X)
predictions = model.predict(X)
print(predictions)

# clustering output and visualization
from collections import defaultdict

groups = list(data.keys())
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
plt.show()


######################
# 2. Sequnece decoding
######################
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO, SearchIO
from Bio.SeqUtils import GC
# creating fasta file from the processed data
seq_analysis.parse_file_write()

# Ques 2a. blast the sequence in Sb008 to identify the function of the sequence, the location in human and identify the genes, the Exons, the introns, the transcription factor binding sites
# Are we using the website for this

# Ques 2b. use bio python to automate the process and analyse the 8 sequences. it takes a while to load
# record_dict = SeqIO.to_dict(SeqIO.parse("cleaned_sb008_data.fasta", "fasta"))
# result_handle = NCBIWWW.qblast("blastn", "nt", record_dict)
# blast_result = open("blasted_sb008_data.xml", "w")
# blast_result.write(result_handle.read())
# print(record_dict)
# print("Sequence length (bp)", len(record_dict))
# blast_result.close()
# result_handle.close()

####################################################################################
# 3. finding the motifs of size 15 that are common to human and chimp Sb008.fastz file
# use meme to extract all the motifs of size 15 that is supporting at least 10 genomic 
# region of the human in Sb008.fastz.
#################################################################################### 
from Bio import motifs, SeqUtils
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
import os

# Ques 3a. Provide the list of these motifs using consensus (human), PWM and logo representation. 
# human sequence only for now
with open("meme_human.xml", "r") as file:
    human_meme = motifs.parse(file, "MEME")
file.close()

print("Motifs for human:")
instances = [meme.consensus for meme in human_meme]
m = motifs.create(instances)
print(m)
print(m.consensus)

# Position Weight Matrix (PWM) and logo representation of Motifs
for i in range(len(human_meme)):
    m = human_meme[i]
    print("The PWM of the motif {}".format(i+1))
    print(m.pwm)
    m.weblogo("WebLogo of the motif "+str(i+1)+".png", "PNG")

# Ques 3b. check the property of these regions containing the motifs with blast. What type of elements are you finding?

# Ques 3c. find the motifs that are size 15 common across the 8 mammals in Sb008.fastz
with open("meme_species.xml", "r") as file:
    species_meme = motifs.parse(file, "MEME")
file.close()

print(species_meme.version, species_meme.sequences)
instances = [meme.consensus for meme in species_meme]
m = motifs.create(instances)
print("Motifs for all species:")
print(m)
print(m.consensus)
for meme in species_meme:
    print(SeqUtils.nt_search(value, m.consensus))


#########################
# NCBI
#########################
# 1. Choose a disease of interest (Breast cancer, Cystic fibrosis, alzheimer)
print("Cystic fibrosis")
# 2. Identify the name of the genes associated to the disease
print("CFTR CF transmembrane conductance regulator [ Homo sapiens (human) ]")
# 3. Find the sequence in human and the location of the gene in Human genome
print("Chromosome 7 - NC_000007.14", "27 exon counts")
# 4. Download the DNA sequences of the gene
print("Download attached in the folder")
# 5. Find the known variant associated to gene in relation of the disease
print("GRCh38/hg38 7p22.3-q36.3(chr7:53985-159282531)x1")