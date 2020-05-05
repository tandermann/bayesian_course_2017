#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 11:26:23 2019

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import os
import numpy as np
import pandas as pd
import glob
from Bio import SeqIO
import matplotlib.pyplot as plt
import itertools
import re
from fuzzywuzzy import process


permutations = np.unique(list(map(''.join, itertools.chain(itertools.product(['A','C','G','T'], ['A','C','G','T']), itertools.product(['A','C','G','T'], ['A','C','G','T'])))))
repetitive_patterns = [''.join([i]*5) for i in permutations]
folder = '/Users/tobias/GitHub/seqcap_processor/data/processed/cleaned_trimmed_reads/1061_clean'
files = glob.glob(os.path.join(folder,'*.fastq'))
filenames = [path.split('/')[-1] for path in files]
for index, file in enumerate(files):
    filename = filenames[index]
    if not 'single' in filename:
        reads_with_microsattelites = []
        for record in SeqIO.parse(file, "fastq"):
            hits = np.array([str(record.seq).find(i) for i in repetitive_patterns])
            if len(hits[hits>-1]) > 0:
                reads_with_microsattelites.append(record)
        len(reads_with_microsattelites)
        
        # write to output
        sequences = reads_with_microsattelites
        with open("/Users/tobias/GitHub/SSR_extraction/data/processed/reads_with_microsattelites/extracted_reads_clean/extracted_reads_%s"%filename, "w") as output_handle:
            SeqIO.write(sequences, output_handle, "fastq")


read_sequences = [str(i.seq) for i in reads_with_microsattelites]

# extract reads that have repetitive pattern and n nucleotides on each side trailing it
all_matches_raw = [[re.findall('[ACTG]{40}%s[ACTG]{40}'%j,i) for i in read_sequences if len(re.findall('[ACTG]{40}%s[ACTG]{40}'%j,i)) > 0] for j in repetitive_patterns]
all_matches = [item[0] for sublist in all_matches_raw for item in sublist]
all_matches_unique = np.unique(all_matches)

file = files[1]
all_reads = [str(record.seq) for record in SeqIO.parse(file, "fastq") if len(record.seq) > 200]

str2Match = all_matches_unique[0]
strOptions = all_reads
Ratios = process.extract(str2Match,strOptions,limit=10)
print(Ratios)
# You can also select the string with the highest matching percentage
highest = process.extractOne(str2Match,strOptions)
print(highest)



# assemble reads (READ1 and READ2 from each sample) into contigs with SECAPR
# BASH: secapr assemble_reads --input data/processed/reads_with_microsattelites/ --output data/processed/assembled_microsattelite_reads/ --disable_stats

contig_file = '/Users/tobias/GitHub/SSR_extraction/data/processed/assembled_microsattelite_reads/extracted_reads.fa'
headers = np.loadtxt(contig_file,dtype=str,delimiter='\t')[::2]
sequences = np.loadtxt(contig_file,dtype=str,delimiter='\t')[1::2]
contig_lengths = np.array([int(i.split(' ')[1]) for i in headers])
output = list(zip(headers[contig_lengths>300],sequences[contig_lengths>300]))

with open('/Users/tobias/GitHub/SSR_extraction/data/processed/assembled_microsattelite_reads/longest_contigs_300.txt', 'w') as out_file:
     [out_file.write('%s\n%s\n'%(i[0],i[1])) for i in output]








