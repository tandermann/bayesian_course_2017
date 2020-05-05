#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 15:36:19 2019

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import numpy as np
import pandas as pd
import gzip
from Bio import SeqIO
import glob
import os
import re
from fuzzysearch import find_near_matches
import sys
sys.path.append("/Users/tobias/GitHub/ssr_extraction/bin/python_functions/")
import ssr_extraction as ssr
from flashtext import KeywordProcessor
from fuzzywuzzy import process


wordsize = 30

# get the ssr markers sequences containing ssr and trailing sequences
ssr_marker_file = '/Users/tobias/GitHub/ssr_extraction/data/processed/ssr_markers/marker_sequences.txt'
ssr_markers = [line.rstrip() for line in open(ssr_marker_file, 'rU').readlines()]

# read the fastq files and iterate through them
read_dir = '/Users/tobias/GitHub/ssr_extraction/data/example_files/'
read_files = glob.glob(os.path.join(read_dir,'*.fq.gz'))

for file in read_files:
    with gzip.open(file, "rt") as handle:
        sequences_sample = [str(record.seq) for record in SeqIO.parse(handle, "fastq")]
        sequences_sample = np.unique(sequences_sample)
        #[[sequence.find(i) for sequence in sequences_sample] for i in ssr_markers]




from sklearn.feature_extraction.text import TfidfVectorizer
from scipy.sparse import csr_matrix
import sparse_dot_topn.sparse_dot_topn as ct


def ngrams(string, n=wordsize):
    ngrams = zip(*[string[i:] for i in range(n)])
    return [''.join(ngram) for ngram in ngrams]


def awesome_cossim_top(A, B, ntop, lower_bound=0):
    # force A and B as a CSR matrix.
    # If they have already been CSR, there is no overhead
    A = A.tocsr()
    B = B.tocsr()
    M, _ = A.shape
    _, N = B.shape
 
    idx_dtype = np.int32
 
    nnz_max = M*ntop
 
    indptr = np.zeros(M+1, dtype=idx_dtype)
    indices = np.zeros(nnz_max, dtype=idx_dtype)
    data = np.zeros(nnz_max, dtype=A.dtype)

    ct.sparse_dot_topn(
        M, N, np.asarray(A.indptr, dtype=idx_dtype),
        np.asarray(A.indices, dtype=idx_dtype),
        A.data,
        np.asarray(B.indptr, dtype=idx_dtype),
        np.asarray(B.indices, dtype=idx_dtype),
        B.data,
        ntop,
        lower_bound,
        indptr, indices, data)

    return csr_matrix((data,indices,indptr),shape=(M,N))

def get_matches_df(sparse_matrix, name_vector, top=100):
    non_zeros = sparse_matrix.nonzero()
    
    sparserows = non_zeros[0]
    sparsecols = non_zeros[1]
    
    if top:
        nr_matches = top
    else:
        nr_matches = sparsecols.size
    
    left_side = np.empty([nr_matches], dtype=object)
    right_side = np.empty([nr_matches], dtype=object)
    similairity = np.zeros(nr_matches)
    
    for index in range(0, nr_matches):
        left_side[index] = name_vector[sparserows[index]]
        right_side[index] = name_vector[sparsecols[index]]
        similairity[index] = sparse_matrix.data[index]
    
    return pd.DataFrame({'left_side': left_side,
                          'right_side': right_side,
                           'similarity': similairity})




sequences = sequences_sample[:1000]

vectorizer = TfidfVectorizer(min_df=3)
seq_tf_idf_matrix = vectorizer.fit_transform(sequences)
ssr_tf_idf_matrix = vectorizer.fit_transform(ssr_markers)
print(ssr_tf_idf_matrix[0])

from sklearn.metrics.pairwise import cosine_similarity

matches = cosine_similarity(seq_tf_idf_matrix)



        


matches_df = get_matches_df(matches, sequences, top=100000)
matches_df = matches_df[matches_df['similairity'] < 0.99999] # Remove all exact matches
matches_df.sample(20)








choices = sequences_sample
process.extract(search_string, choices, limit=2)



sequences_text = ' '.join(sequences_sample+['GACATATTCTCACATTATATAATCTGTTGCATATATATAT'])
search_string = ssr_markers[0]
processor = KeywordProcessor()
processor.add_keyword(search_string)
#processor.add_keywords_from_dict({'batman':['batman','bruce wayne']})
found = processor.extract_keywords(sequences_text, span_info=True)






test_string = 'CATGCTGGCTGTCCTCCACACAGCCCAGCTAGGTCCCCCATTATATATAATGTTGTCAAATTAATACACGACACACACACTCTGATACCACATTGCATGG'


for i in ssr_markers:
    sim = similarity(i,test_string)
    if sim >= 0.55:
        print(i)
    
    matches = find_near_matches(i, sequences_sample, max_deletions=2, max_insertions=2, max_substitutions=2)
    print(matches)
#   break
#



import difflib

def similarity(word, pattern):
    return difflib.SequenceMatcher(a=word.lower(), b=pattern.lower()).ratio()

text = "Somme text with google or gooole or goofle";
lookup = "google";
threshold = 0.8

for word in text.split():
    if similarity(word, lookup) > threshold:
        print(word)






test = sequences_sample[:100]

# 1. search allowing up to 3 substitutions


        
        





    hits = np.array([sequence.find(i) for sequence in sequences_sample])
    if len(hits[hits>-1]) > 1:
        break

np.array(sequences_sample)[hits>-1]



#________________CODE SCRAPS___________________________________________________
#from fuzzywuzzy import process
#from fuzzywuzzy import fuzz      
#str2Match = i
#strOptions = sequences_sample
#Ratios = process.extract(str2Match,strOptions,limit=10)
#process.extractOne(i,sequences_sample,scorer=fuzz.partial_ratio)
#fuzz.partial_ratio(i,sequences_sample[0])


