#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 15:38:41 2019

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import numpy as np
import pandas as pd
import os
from Bio import SeqIO
import glob
import sys
from Bio.Align.Applications import MafftCommandline

def fix_line_wrap(alignment_file):
    file_in = open(alignment_file)
    final = {}    
    for line in file_in:
        line = line.strip()
        if line.startswith(">"):
            id = line
            final[id] = " "
        else:
            final[id] += line    
    file_out = open(alignment_file, "w")    
    for key, val in final.items():
        file_out.write(key)
        file_out.write("\n")
        file_out.write(val)
        file_out.write("\n")

# get the names of the loci that were recovered in all samples
match_table = '/Users/tobias/GitHub/ssr_extraction/data/processed/reads_matching_markers/match_table.txt'
match_table_df = pd.read_csv(match_table,sep='\t',index_col=0)
all_present_markers = list(match_table_df[np.sum(match_table_df.values,axis=1)==len(match_table_df.columns)].index.values)

# get all the target sequence collection files
alignment_folder = '/Users/tobias/GitHub/ssr_extraction/data/processed/alignments_reads_with_markers'
sample_folders = [os.path.basename(name) for name in os.scandir(alignment_folder) if os.path.isdir(name)]
# make sure we don't get the output folder of this script, in case it has run previously
sample_folders = [name for name in  sample_folders if not 'all_samples_merged' in name]
all_target_files = []
for sample in sample_folders:
    sequence_collections_path = os.path.join(alignment_folder,sample,'paralog_seq_collections')
    files = [os.path.join(sequence_collections_path,'paralog_reads_collection_locus_%s.fasta'%i) for i in all_present_markers]
    all_target_files.append(files)
all_target_files = [item for sublist in all_target_files for item in sublist]

# read and merge seq-collections of all samples for each marker
outdir = '/Users/tobias/GitHub/ssr_extraction/data/processed/alignments_reads_with_markers/all_samples_merged'
if not os.path.exists(outdir):
    os.makedirs(outdir)
outdir_seqcol = os.path.join(outdir,'sequence_collections')
if not os.path.exists(outdir_seqcol):
    os.makedirs(outdir_seqcol)
for marker in all_present_markers:
    target_files = [filename for filename in all_target_files if marker in filename]
    all_seq_export = []
    for i, file in enumerate(target_files):
        seq_data = list(SeqIO.parse(file, "fasta"))
        sample = os.path.dirname(os.path.dirname(file.replace(alignment_folder,''))).replace('/','')
        for seq_obj in seq_data:
            if not 'marker_sequence' in seq_obj.name:
                seq_obj.id = sample+'_'+seq_obj.id
        # make sure to only get the ref seq one time
        if i == 0:
            sequences = seq_data
        else:
            sequences = seq_data[1:]
        all_seq_export.append(sequences)
    all_seq_export = [item for sublist in all_seq_export for item in sublist]
    SeqIO.write(all_seq_export, os.path.join(outdir_seqcol,'paralog_reads_collection_all_samples_locus_%s.fasta'%marker),'fasta')

# align the sequence collections
outdir_align = os.path.join(outdir,'sequence_alignments')
if not os.path.exists(outdir_align):
    os.makedirs(outdir_align)


# align the sequences and print as fasta alignment file
seq_colls = glob.glob(os.path.join(outdir_seqcol,'*.fasta'))
for counter, sequence_collection in enumerate(seq_colls):
    filename = sequence_collection.split('/')[-1].replace('paralog_reads_collection_all_samples_','alignment_')
    cline = MafftCommandline(input=sequence_collection,op=1.53,ep=0.123)
    stdout, stderr = cline()
    alignment_out = os.path.join(outdir_align,filename)
    sys.stdout.write('\rAligning sequence collections %i/%i '%(int(counter+1),len(all_present_markers)))
    sys.stdout.flush()
    with open(alignment_out, "w") as handle:
        handle.write(stdout)
    fix_line_wrap(alignment_out)



