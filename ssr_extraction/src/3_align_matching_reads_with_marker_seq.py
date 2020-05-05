#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 11:59:53 2019

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import numpy as np
import pandas as pd
import shutil
import os
import sys
import glob
from numpy import genfromtxt
from Bio import SeqIO
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



marker_fasta = '/Users/tobias/GitHub/ssr_extraction/data/marker_sequences/marker_sequences.fasta'
extracted_reads_repo = '/Users/tobias/GitHub/ssr_extraction/data/reads_with_markers/*.fasta'
outdir = '/Users/tobias/GitHub/ssr_extraction/data/extracted_reads_aligned_with_markers/'
if not os.path.exists(outdir):
    os.makedirs(outdir)
seqcol_outdir = os.path.join(outdir,'sequence_collections')
if not os.path.exists(seqcol_outdir):
    os.makedirs(seqcol_outdir)
align_outdir = os.path.join(outdir,'alignments')
if not os.path.exists(align_outdir):
    os.makedirs(align_outdir)
marker_fasta_content = list(SeqIO.parse(marker_fasta, "fasta"))
ref_seq_dict = dict([(fasta.name,fasta) for fasta in marker_fasta_content])

sequence_files = glob.glob(extracted_reads_repo)
for sequence_file in sequence_files:
    records = []
    locus_name = os.path.basename(sequence_file).replace('_reads_all_samples.fasta','')
    ref_seq = ref_seq_dict[locus_name]
    records.append(ref_seq)
    read_collection = list(SeqIO.parse(sequence_file, "fasta"))
    records = records + read_collection
    SeqIO.write(records, os.path.join(seqcol_outdir,'sequence_collection_locus_%s.fasta'%locus_name),'fasta')



# align the sequences and print as fasta alignment file
seq_colls = glob.glob(os.path.join(seqcol_outdir,'*.fasta'))
for counter, sequence_collection in enumerate(seq_colls):
    #locus_name = 'marker_sequence_'+sequence_collection.replace('.fasta','').split('_')[-1]
    #ref_seq = ref_seq_dict[locus_name]
    #marker_out = os.path.join(seqcol_outdir,'marker_seq_locus_%s.fasta'%locus_name)
    #SeqIO.write(ref_seq, marker_out,'fasta')
    filename = os.path.basename(sequence_collection).replace('sequence_collection','alignment')
    cline = MafftCommandline(input=sequence_collection,op=1.53,ep=1.0)
    stdout, stderr = cline()
    alignment_out = os.path.join(align_outdir,filename)
    sys.stdout.write('\rAligning sequence collections %i/%i '%(int(counter+1),len(seq_colls)))
    sys.stdout.flush()
    with open(alignment_out, "w") as handle:
        handle.write(stdout)
    fix_line_wrap(alignment_out)













# read_folder = '/Users/tobias/GitHub/ssr_extraction/data/example_files'
# ref_file = '/Users/tobias/GitHub/ssr_extraction/data/processed/ssr_markers/marker_sequences.fasta'
# root_dir = '/Users/tobias/GitHub/ssr_extraction/data/processed/reads_matching_markers'
# outdir = '/Users/tobias/GitHub/ssr_extraction/data/processed/alignments_reads_with_markers'
# run_sequence_extraction = True

# if run_sequence_extraction:
#     if os.path.exists(outdir):
#         shutil.rmtree(outdir)
#     os.makedirs(outdir)



# subdirs = list(os.walk(root_dir))[0][1]
# #subdirs = list(root_dir)
# for sample in subdirs:
#     print('\nProcessing sample %s'%sample)
#     # get the paths for the sample
#     sample_path = os.path.join(root_dir,sample)
#     read_orientation = os.path.join(sample_path,'read_orientation.txt')
#     para_info = os.path.join(sample_path,'info_paralogous_loci.txt')
#     para_data = genfromtxt(para_info, delimiter='\t', dtype=str)
#     read_file = os.path.join(read_folder,'%s.fq'%sample)
#     sample_out_dir = os.path.join(outdir,sample)
#     sample_sequence_outdir = os.path.join(sample_out_dir,'paralog_seq_collections')

#     if run_sequence_extraction:
#         # read the data
#         ref_seqs = list(SeqIO.parse(ref_file, "fasta"))
#         read_seqs = list(SeqIO.parse(read_file, "fastq"))
#         fastq_headers = [i.name for i in read_seqs]
#         read_orientation_df = pd.read_csv(read_orientation,sep='\t')
#         print('%i markers with matching reads found.'%len(para_data))

#         if not os.path.exists(sample_out_dir):
#             os.makedirs(sample_out_dir)
#         if not os.path.exists(sample_sequence_outdir):
#             os.makedirs(sample_sequence_outdir)
#         # print ref and read sequences for each paralogous locus into a separate sequence collection
#         for counter,i in enumerate(para_data):
#             records = []
#             marker = i[0]
#             marker_id = marker
#             ref_seq = [ref for ref in ref_seqs if marker_id==ref.id][0]
#             records.append(ref_seq)
#             read_list_tmp = i[1:]
#             read_list = np.unique(read_list_tmp[read_list_tmp!=''])
#             for read_id in read_list:
#                 read_seq = read_seqs[fastq_headers.index(read_id)]
#                 orientation = read_orientation_df[read_orientation_df.read_id == read_id].orientation.values[0]
#                 if orientation == '+':
#                     pass
#                 else:
#                     read_seq.seq = read_seq.seq.reverse_complement()
#                 records.append(read_seq)
#             sys.stdout.write('\rPrinting sequence collections %i/%i '%(int(counter+1),len(para_data)))
#             sys.stdout.flush()
#             SeqIO.write(records, os.path.join(sample_sequence_outdir,'paralog_reads_collection_locus_%s.fasta'%marker),'fasta')

#     # align the sequences and print as fasta alignment file
#     sample_alignment_outdir = os.path.join(sample_out_dir,'paralog_alignments')
#     if not os.path.exists(sample_alignment_outdir):
#         os.makedirs(sample_alignment_outdir)
#     seq_colls = glob.glob(os.path.join(sample_sequence_outdir,'*.fasta'))
#     for counter, sequence_collection in enumerate(seq_colls):
#         filename = sequence_collection.split('/')[-1].replace('paralog_reads_collection_','alignment_')
#         cline = MafftCommandline(input=sequence_collection,op=1.53,ep=1.0)
#         stdout, stderr = cline()
#         alignment_out = os.path.join(sample_alignment_outdir,filename)
#         sys.stdout.write('\rAligning sequence collections %i/%i '%(int(counter+1),len(para_data)))
#         sys.stdout.flush()
#         with open(alignment_out, "w") as handle:
#             handle.write(stdout)
#         fix_line_wrap(alignment_out)
        
