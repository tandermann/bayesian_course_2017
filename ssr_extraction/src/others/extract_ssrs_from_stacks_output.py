#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 10:04:48 2019

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import numpy as np
import pandas as pd
import itertools
import os
import re
import sys
sys.path.append("/Users/tobias/GitHub/ssr_extraction/bin/python_functions/")
import ssr_extraction as ssr
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO

#sequences_shared_by_n_samples = 40
delta_bp_start_marker = 30
buffer_bp = 5
cluster_merge_similarity_threshold = 0.1

stacks_in = '/Users/tobias/GitHub/ssr_extraction/data/raw/stacks/batch_1.catalog.tags.tsv'
outdir = '/Users/tobias/GitHub/ssr_extraction/data/processed/ssr_markers/'
run=True

if run:
    # define the rep. patterns to search for
    permutations = np.unique(list(map(''.join, itertools.chain(itertools.product(['A','C','G','T'], ['A','C','G','T']), itertools.product(['A','C','G','T'], ['A','C','G','T'])))))
    # remove repeats of same nucleotide, e.g. AA
    permutations = [i for i in permutations if not i[0] == i[1]]
    repetitive_patterns = [''.join([i]*5) for i in permutations]
    
    # read sequences from stacks catalog
    stacks_catalog = pd.read_csv(stacks_in,sep='\t',header=None,skiprows=[0])
    stacks_catalog.columns = ['sql_id','sample_id','locus_id','chromosome','basepair','strand','sequence_type','stack_component','sequence_id','sequence','delivered_flag','blacklisted_flag','lumberjackstack_flag','log_likelihood']
    sequences = stacks_catalog.sequence.values
    
    
    # search each sequence for repetitive patterns
    hits = np.array([[sequence.find(i) for sequence in sequences] for i in repetitive_patterns])
    
    # summarize the matches of repetitive patterns for each sequence
    matches = np.array([len(row[row>-1]) for row in hits.T])
    # extract those reads where we had exactly one match of a repetitive pattern (not multiple)
    n_singular_matches = len(matches[matches==1])
    match_sequences = sequences[matches==1]
    
    # find motifs containing rep-patterns and trailing sequences
    
    
    # the following line ensures that only those sequences are reported that have a sufficient sequence buffer before the motif, yet only the part before the motif is exported, including a =nbp buffer
    target_markers_step1 = [[re.findall('[ACTG]{%i}%s'%(delta_bp_start_marker,j),i) for i in match_sequences if len(re.findall('[ACTG]{%i}%s'%(delta_bp_start_marker,j),i)) > 0] for j in repetitive_patterns]
    #target_markers_step1 = [[re.findall('[ACTG]{%i}%s[ACTG]{%i}'%(trailing_basepairs,j,trailing_basepairs),i) for i in match_sequences if len(re.findall('[ACTG]{%i}%s[ACTG]{%i}'%(trailing_basepairs,j,trailing_basepairs),i)) > 0] for j in repetitive_patterns]
    target_markers_step2 = [item for sublist in target_markers_step1 for item in sublist]
    target_markers = target_markers_step2
    #target_markers = [i[0][:(delta_bp_start_marker-buffer_bp)] for i in target_markers_step2]
    # sort out duplicate motifs
    target_markers = np.unique(target_markers)
    # create similarity matrix and cluster similar markers
    matrix = ssr.create_simmatrix(target_markers,target_markers)
    clusters = ssr.smart_cluster_based_on_matrix(matrix,target_markers,cluster_merge_similarity_threshold)
    # write the identified motif clusters to file
    clusters.to_csv(os.path.join(outdir,'motif_clusters.txt'), sep = '\t', index = False, header=False)    
    # get one representative motif per cluster
    pick_per_cluster = [np.random.choice(i[1][0].split(', ')) for i in clusters.iterrows()]
    # write motifs to file
    with open(os.path.join(outdir,'marker_sequences.txt'), 'w') as out_file:
         [out_file.write('%s\n'%i) for i in pick_per_cluster]
else:
    pick_per_cluster = pd.read_csv(os.path.join(outdir,'marker_sequences.txt'),header=None).values
    pick_per_cluster = [i[0] for i in pick_per_cluster]


# write in fasta format
sequence_collection = []
for i in range(len(pick_per_cluster)):
    sequence = pick_per_cluster[i]
    sequence_collection.append(SeqRecord(seq=Seq(sequence, IUPAC.ambiguous_dna), id='marker_sequence_%i'%i, name='marker_sequence_%i'%i,description=''))

SeqIO.write(sequence_collection, os.path.join(outdir,'marker_sequences.fasta'), 'fasta')





# CODE SCRAPS__________________________________________________________________
## get the sequences that are shared by a minimum number of samples
#number_of_hits = [len(np.unique([j.split('_')[0] for j in i.split(',')])) for i in stacks_catalog.sequence_id.values]
## how many samples represent n % of all samples, set by sequences_shared_by_n_percent_samples flag
#common_sequences = np.array(sequences)[np.array(number_of_hits)>=sequences_shared_by_n_samples]
#common_sequence_ids = stacks_catalog.index[np.array(number_of_hits)>=sequences_shared_by_n_samples].values
#
## search each sequence for repetitive patterns
#hits = np.array([[sequence.find(i) for sequence in common_sequences] for i in repetitive_patterns])
#
## summarize the matches of repetitive patterns for each sequence
#matches = np.array([len(row[row>-1]) for row in hits.T])
## extract those reads where we had exactly one match of a repetitive pattern (not multiple)
#n_singular_matches = len(matches[matches==1])
#match_sequences = common_sequences[matches==1]
