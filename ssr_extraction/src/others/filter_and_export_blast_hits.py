#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 08:00:10 2019

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import numpy as np
np.set_printoptions(suppress=True)
import pandas as pd
import os,glob
import subprocess
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

dev_mode = False
min_similarity = 0.8
min_length_fraction = 0.9
target_length = 15
seed_length = 10
marker_sequences = '/Users/tobias/GitHub/ssr_extraction/data/marker_sequences/marker_sequences.fasta'
reads_folder = '/Users/tobias/GitHub/ssr_extraction/data/fastq'
out_dir = '/Users/tobias/GitHub/ssr_extraction/data/reads_with_markers'
clean_up = False

read_files = glob.glob(reads_folder+'/*.fq.gz')
if len(read_files) == 0:
    read_files = glob.glob(reads_folder+'/*.fastq.gz')
    
tmpdir = os.path.join(reads_folder,'temp')

if not os.path.exists(tmpdir):
    os.makedirs(tmpdir)
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

marker_match_dict = {}
for fastq_file in read_files:
    filename = fastq_file.split('/')[-1].replace('.fq.gz','')
    print('Processing sample',filename)

    # convert fastq into fasta adn store as temporary file
    print('Converting fastq into fasta (temporary files) ...')
    seqtk_err = os.path.join(tmpdir,'%s_seqtk_screen.txt'%filename)
    seqtk_out = os.path.join(tmpdir,'%s.fasta'%filename)
    seqtk_cmd = ['seqtk', 'seq', '-A', fastq_file]
    if not dev_mode:
        with open(seqtk_out, 'w') as out, open(seqtk_err, 'w') as err:
            seqtk = subprocess.Popen(seqtk_cmd, stderr = err, stdout=out)
            seqtk.wait()

    # create blast database for file
    print('Creating BLAST database (temporary files) ...')
    makeblastdb_cmd = ['makeblastdb', '-in', seqtk_out, '-dbtype', 'nucl']
    makeblastdb_err = os.path.join(tmpdir,'%s+makeblastdb_screen.txt'%filename)
    if not dev_mode:
        with open(makeblastdb_err, 'w') as out:
            makeblastdb = subprocess.Popen(makeblastdb_cmd, stdout=out)
            makeblastdb.wait()

    # run blast between read db and marker sequences
    print('Finding all reads matching markers with BLAST (temporary files) ...')
    blast_cmd = ['blastn', '-task', 'blastn-short', '-db', seqtk_out, '-query', marker_sequences, '-word_size', str(seed_length), '-gapopen', '100', '-strand', 'both', '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand']
    blast_out = os.path.join(tmpdir,'%s_all_blast_hits.txt'%filename)
    blast_err = os.path.join(tmpdir,'%s_blast_screen.txt'%filename)
    #print(blast_cmd)
    with open(blast_out, 'w') as out, open(blast_err, 'w') as err:
        blast = subprocess.Popen(blast_cmd, stderr = err, stdout=out)
        blast.wait()

    # filtering best matches
    print('Filtering best matches, using min_similarity and min_length_fraction values ...')
    selected_matches_out = os.path.join(out_dir,'%s_selected_blast_hits.txt'%filename)
    blast_hits = pd.read_csv(blast_out,sep='\t',header=None)
    blast_hits.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','sstrand']
    selected_matches_sim = blast_hits[blast_hits.pident>=min_similarity*100]
    selected_matches = selected_matches_sim[selected_matches_sim.length>=min_length_fraction*target_length]
    selected_matches.to_csv(selected_matches_out,sep='\t',index=False)
    
    for marker_matches in selected_matches.groupby(['qseqid']):
        marker = marker_matches[0]
        data = marker_matches[1]
        marker_match_dict.setdefault(marker,[])
        marker_match_dict[marker].append([fastq_file,data.sseqid.values,data.sstrand.values])

    # remove temporary directory if clean mode is activated (default):
    if clean_up:
        print('Removing temporary files.')
        files = glob.glob(tmpdir+'/*')
        for f in files:
            os.remove(f)



# extract sequences for matching reads for each locus
for marker in marker_match_dict.keys():
    print('Extracting sequences for',marker)
    with open(os.path.join(out_dir,'%s_reads_all_samples.fasta'%marker), "w") as out_file:
        for sample_data in marker_match_dict[marker]:
            fastq_file = sample_data[0]
            sample_name = fastq_file.split('/')[-1].replace('.fq.gz','')
            list_of_match_headers = list(sample_data[1])
            # export the fastq headers to temporary test file
            header_list_file = os.path.join(tmpdir,'%s_fastq_headers.txt'%marker)
            np.savetxt(header_list_file,list_of_match_headers,fmt='%s')
            # extract fastq sequences with seqtk
            seqtk_fastq_out = os.path.join(tmpdir,'%s_fastq_reads.fastq'%marker)
            cmd = ['seqtk', 'subseq', fastq_file, header_list_file]
            with open(seqtk_fastq_out, 'w') as out:
                run = subprocess.Popen(cmd, stdout=out)
                run.wait()        
            read_seqs = list(SeqIO.parse(seqtk_fastq_out, "fastq"))
            orientations = list(sample_data[2])
            for i,seq in enumerate(read_seqs):
                orientation = orientations[i]
                if orientation == 'minus':
                    sequence = seq.seq.reverse_complement()
                else:
                    sequence = seq.seq
                new_fasta = SeqRecord(sequence, id=marker+'_'+sample_name+'_'+str(i), name='', description='')
                out_file.write(new_fasta.format('fasta-2line'))

       




















