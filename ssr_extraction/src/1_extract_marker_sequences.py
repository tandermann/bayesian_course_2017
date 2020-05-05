#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 11:46:16 2020

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import numpy as np
import itertools
import pandas as pd
import re, os, glob
from fuzzywuzzy import fuzz
from sklearn.cluster import DBSCAN
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
import gzip

def check_nucleotide_div(sequence,n_nucleotide_types_required = 4):
    # this can be set to any number , e.g. <= 0 requires that all 4 nucleotides are present, <= 1 requires 3 out of 4 etc...
    return len(np.where(np.array([sequence.count(i) for i in ['A','C','T','G']]) == 0)[0]) <= (4-n_nucleotide_types_required)

def logent(x):  
    if x<=0:     
        return 0  
    else:  
        return -x*np.log(x)  

def entropy(sequence):
    lisfreq=[sequence.count(base)*1.0/len(sequence) for base in ["A","C","G","T"]]
    return sum([logent(elem) for elem in lisfreq])


def compute_similarity(s1, s2):
    return 1.0 - (0.01 * max(
        fuzz.ratio(s1, s2),
        fuzz.token_sort_ratio(s1, s2),
        fuzz.token_set_ratio(s1, s2)))
        
def create_simmatrix(a,b):
    X = np.zeros((len(a), len(b)))
    for i in range(len(a)):
        #print every 100 rows info
        if i > 0 and i % 100 == 0:
            print("Processed %d/%d sequences" % (i, X.shape[0]))
        for j in range(len(b)):
            if X[i, j] == 0.0:        
                X[i, j] = compute_similarity(a[i].lower(), b[j].lower())
                #X[j, i] = X[i, j]
    return X

def smart_cluster_based_on_matrix(scipy_matrix, original_list,sim_value):
    # takes is scipy matrix and the list  and returns a pandas df with the identified clusters
    clust = DBSCAN(eps=float(sim_value), min_samples=1, metric="precomputed")
    clust.fit(scipy_matrix)
    # print cluster report
    preds = clust.labels_
    clabels = np.unique(preds)
    output = []
    for i in range(clabels.shape[0]):
        if clabels[i] < 0:
            continue
        cmem_ids = np.where(preds == clabels[i])[0]
        cmembers = []
        for cmem_id in cmem_ids:
            cmembers.append(original_list[cmem_id])
        #if clabels.shape[0] > 11:
        #    if i < 11:
        #        print("Cluster#%d: %s" % (i, ", ".join(cmembers)))
        #    elif i == 11:
        #        print("...")
        #else:
        #    print("Cluster#%d: %s" % (i, ", ".join(cmembers)))
        output.append([", ".join(cmembers)])
    final_output = pd.DataFrame(output)
    return final_output


motif_repeat = 6
marker_length = 25
buffer_between_marker_and_motif = 10
buffer_after_motif = 10
min_entropy_score = 2.7#3.8
outdir = '/Users/tobias/GitHub/ssr_extraction/data/marker_sequences'
if not os.path.exists(outdir):
    os.makedirs(outdir)

# 1: Define the motifs
# define the rep. patterns to search for
permutations = np.unique(list(map(''.join, itertools.chain(itertools.product(['A','C','G','T'], ['A','C','G','T']), itertools.product(['A','C','G','T'], ['A','C','G','T'])))))
# remove repeats of same nucleotide, e.g. AA
permutations = [i for i in permutations if not i[0] == i[1]]
motifs = [''.join([i]*motif_repeat) for i in permutations]
#motifs = ['ACGTACGTACGTACGTACGT']


# get the fastq sequences and join all into one list
#headers = []
sequences = []
fastq_files = glob.glob('/Users/tobias/GitHub/ssr_extraction/data/fastq_modified/*.gz')
for fastqfile in fastq_files:
    if fastqfile.endswith('fastq.gz') or fastqfile.endswith('fq.gz'):
        print('Processing file', os.path.basename(fastqfile))
        fastq_sequences = [str(record.seq) for record in SeqIO.parse(gzip.open(fastqfile, "rt"), "fastq")]
        #fastq_content = [[record.id,str(record.seq)] for record in SeqIO.parse(gzip.open(fastqfile, "rt"), "fastq")]
        #fastq_headers,fastq_sequences = zip(*fastq_content)
        #headers += fastq_headers
        sequences += fastq_sequences
sequences = np.array(sequences)

# find all reads that contain the motif
indices = [np.flatnonzero(np.core.defchararray.find(sequences,i)!=-1) for i in motifs]

# get the string of nucleotides before the motif as marker
potential_primers = np.array([[index_list[j],re.findall('[ACTG]{%i}%s[ACTG]{%i}'%(marker_length+buffer_between_marker_and_motif,motifs[i],buffer_after_motif),sequence)[0][:marker_length]] for i,index_list in enumerate(indices) for j, sequence in enumerate(sequences[index_list].astype(str)) if len(re.findall('[ACTG]{%i}%s[ACTG]{%i}'%(marker_length+buffer_between_marker_and_motif,motifs[i],buffer_after_motif),sequence)) == 1 ])
sequence_indices = potential_primers[:,0].astype(int)
potential_primer_sequences = potential_primers[:,1].astype(str)

# run the potential primers through nucleotide diversity filter
#target_indices = [i for i,primer in enumerate(potential_primer_sequences) if check_nucleotide_div(primer,4)]
#sequence_indices = sequence_indices[target_indices]
#potential_primer_sequences = potential_primer_sequences[target_indices]

# run the potential primers through entropy filter
target_indices = [i for i,primer in enumerate(potential_primer_sequences) if np.exp(entropy(primer))>=min_entropy_score]
sequence_indices = sequence_indices[target_indices]
potential_primer_sequences = potential_primer_sequences[target_indices]
# remove all duplicates
potential_primer_sequences = np.unique(potential_primer_sequences)

# cluster the primer sequences and select only 1 per cluster
matrix = create_simmatrix(potential_primer_sequences,potential_primer_sequences)
clusters = smart_cluster_based_on_matrix(matrix,potential_primer_sequences,0.1)
# write the identified motif clusters to file
clusters.to_csv(os.path.join(outdir,'motif_clusters.txt'), sep = '\t', index = False, header=False)    
# get one representative motif per cluster
pick_per_cluster = sorted([np.random.choice(i[1][0].split(', ')) for i in clusters.iterrows()])
# write motifs to file
with open(os.path.join(outdir,'marker_sequences.txt'), 'w') as out_file:
     [out_file.write('%s\n'%i) for i in pick_per_cluster]

# write in fasta format
sequence_collection = []
for i in range(len(pick_per_cluster)):
    sequence = pick_per_cluster[i]
    sequence_collection.append(SeqRecord(seq=Seq(sequence, IUPAC.ambiguous_dna), id='marker_sequence_%i'%i, name='marker_sequence_%i'%i,description=''))

SeqIO.write(sequence_collection, os.path.join(outdir,'marker_sequences.fasta'), 'fasta')







