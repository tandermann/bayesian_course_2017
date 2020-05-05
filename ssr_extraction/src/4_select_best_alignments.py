#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 15:43:36 2020

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import numpy as np
np.set_printoptions(suppress=True)
import pandas as pd
import matplotlib.pyplot as plt
np.random.seed(1234)
import os,glob,shutil
from Bio import SeqIO


def get_fasta_header_list(fasta_file_path):
    headers = ['_'.join(i.id.split('_')[:-1]) for i in list(SeqIO.parse(fasta_file_path, "fasta"))]
    return headers[1:]


min_samples = 3
min_read_count_per_sample = 4


al_folder = '/Users/tobias/GitHub/ssr_extraction/data/extracted_reads_aligned_with_markers/alignments/*.fasta'
al_files = np.array(glob.glob(al_folder))

good_files = np.array([i for i in al_files if len(np.unique(get_fasta_header_list(i),return_counts=True)[0])>=min_samples and np.all(np.greater(np.unique(get_fasta_header_list(i),return_counts=True)[1],min_read_count_per_sample-1))])




file_sizes = np.array([[os.path.getsize(file),i] for i,file in enumerate(good_files)]).astype(float)
file_sizes[:,0] = np.log(file_sizes[:,0])
log_file_sizes = file_sizes



mean = np.mean(log_file_sizes[:,0])
std = np.std(log_file_sizes[:,0])
log_file_sizes = log_file_sizes[log_file_sizes[:,0]>(mean-std)]
log_file_sizes = log_file_sizes[log_file_sizes[:,0]<(mean+std)]

plt.hist(log_file_sizes[:,0],50)
indeces = log_file_sizes[:,1].astype(int)
selected_files = good_files[indeces]
destination_path = al_folder.replace('*.fasta','selected_alignments')
os.makedirs(destination_path)


[shutil.copy(src_file,os.path.join(destination_path,os.path.basename(src_file))) for src_file in selected_files]


