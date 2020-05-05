#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 09:39:15 2019

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import numpy as np
import pandas as pd
import subprocess
import re
import os
import glob
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, AlignIO


def reads_matching_markers(lastz_df):
    # make a dictionary with all read names that match a ssr marker
    marker_read_dict = {}
    read_marker_dict = {}
    read_orientation_dict = {}
    read_multi_marker_dict = {}
    for row in lastz_df.iterrows():
        locus = row[1].name2
        locus_name = re.sub('\_p[0-9]* \|.*', '', locus)
        locus_name = re.sub('^>', '', locus_name)
        read_header = row[1].name1
        #print(read_header)
        #read_name = re.sub('^\>([0-9]*) .*', '\\1', read_header)
        read_name = re.sub('^\>([^\s]*) .*', '\\1', read_header)
        read_name = re.sub('^>', '', read_name)
        #print(read_name)
        marker_read_dict.setdefault(locus_name,[])
        marker_read_dict[locus_name].append(read_name)
        read_marker_dict.setdefault(read_name,[])
        read_marker_dict[read_name].append(locus_name)
        orientation = row[1].strand2
        read_orientation_dict.setdefault(read_name,orientation)
    
    orientation_df = pd.DataFrame(index=np.arange(0,len(read_orientation_dict)), columns=['read_id','orientation'])
    orientation_df['read_id'] = read_orientation_dict.keys()
    orientation_df['orientation'] = read_orientation_dict.values()
    for read in read_marker_dict.keys():
        if len(read_marker_dict[read]) > 1:
            read_multi_marker_dict.setdefault(read,read_marker_dict[read])
    return marker_read_dict, read_marker_dict, read_orientation_dict, read_multi_marker_dict, orientation_df


def find_duplicates(marker_read_dict,read_marker_dict):
    # get markers that have multiple reads matching them
    invalid_marker_loci = []
    markers_with_multiple_hits = []
    for marker in marker_read_dict.keys():
        if len(marker_read_dict[marker]) > 1:
            markers_with_multiple_hits.append(marker)
            invalid_marker_loci.append(marker)
    # get markers that match on multiple reads
    reads_matching_multiple_markers = []
    for read in read_marker_dict.keys():
        if len(read_marker_dict[read]) > 1:
            reads_matching_multiple_markers.append(read)
            for marker in read_marker_dict[read]:
                invalid_marker_loci.append(marker)
    return invalid_marker_loci, markers_with_multiple_hits, reads_matching_multiple_markers


def find_longest_read(read_names,lastz_df):
    read_header_values = np.array([i.split(' ')[0].replace('>','') for i in lastz_df.name1.values if i.split(' ')[0].replace('>','') in read_names]).astype(int)
    read_length_values = np.array([i.split(' ')[1] for i in lastz_df.name1.values if i.split(' ')[0].replace('>','') in read_names]).astype(int)
    longest_read = read_header_values[list(read_length_values).index(np.max(read_length_values))]
    return longest_read


#marker_read_dict,duplicate_loci,markers_with_multiple_hits,reads_matching_multiple_marker_dict,outdir,lastz_df = marker_read_dict,duplicate_loci,possible_paralogous,read_multi_marker_dict,subfolder,lastz_df
def get_list_of_valid_markers_and_reads(marker_read_dict,duplicate_loci,markers_with_multiple_hits,reads_matching_multiple_marker_dict,outdir,lastz_df,keep_reads_matching_multiple_markers=True,keep_markers_with_multiple_read_hits=True):
    # summarize all markers that should be excluded form further processing (duplicates)
    if keep_reads_matching_multiple_markers:
        # then only mark the list markers_with_multiple_hits as bad markers
        invalid_markers_temp = list(set(markers_with_multiple_hits))
        valid_reads_matching_multiple_marker_dict = {}
        for marker in reads_matching_multiple_marker_dict.keys():
            if not marker in invalid_markers_temp:
                valid_reads_matching_multiple_marker_dict.setdefault(marker,reads_matching_multiple_marker_dict[marker])
        dupl_info = pd.DataFrame.from_dict(valid_reads_matching_multiple_marker_dict, orient='index')
        dupl_info.to_csv(os.path.join(outdir,'info_reads_matching_multiple_markers.txt'),header=False,sep="\t")
    if keep_markers_with_multiple_read_hits:
        invalid_markers_unique = []
        invalid_markers_temp = list(set(duplicate_loci))
        paralogous_markers = {}
        for marker in markers_with_multiple_hits:
            paralogous_markers.setdefault(marker,marker_read_dict[marker])
        paralog_info = pd.DataFrame.from_dict(paralogous_markers, orient='index')
        paralog_info.to_csv(os.path.join(outdir,'info_paralogous_loci.txt'),header=False,sep="\t")
        print('Found %i markers with multiple read hits. Selecting one random read per cluster.'%len(invalid_markers_temp))
    else:
        # remove all duplicates
        invalid_markers_unique = list(set(duplicate_loci))
        print(len(invalid_markers_unique), 'possibly paralogous markers detected - excluded from processing')
    # get list of valid read names
    valid_read_names = []
    for marker in marker_read_dict:
        if marker not in invalid_markers_unique:
            read_name = marker_read_dict[marker]
            read_name = np.random.choice(read_name)
            #read_name = find_longest_read(read_name,lastz_df)
            valid_read_names.append(str(read_name).replace('>',''))
    return valid_read_names


def get_list_of_valid_markers_and_reads_including_all_loci_in_paralog_output_file(marker_read_dict,duplicate_loci,markers_with_multiple_hits,reads_matching_multiple_marker_dict,outdir,lastz_df,keep_reads_matching_multiple_markers=True,keep_markers_with_multiple_read_hits=True):
    # summarize all markers that should be excluded form further processing (duplicates)
    if keep_reads_matching_multiple_markers:
        # then only mark the list markers_with_multiple_hits as bad markers
        invalid_markers_temp = list(set(markers_with_multiple_hits))
        valid_reads_matching_multiple_marker_dict = {}
        for marker in reads_matching_multiple_marker_dict.keys():
            if not marker in invalid_markers_temp:
                valid_reads_matching_multiple_marker_dict.setdefault(marker,reads_matching_multiple_marker_dict[marker])
        dupl_info = pd.DataFrame.from_dict(valid_reads_matching_multiple_marker_dict, orient='index')
        dupl_info.to_csv(os.path.join(outdir,'info_reads_matching_multiple_markers.txt'),header=False,sep="\t")
    if keep_markers_with_multiple_read_hits:
        invalid_markers_unique = []
        invalid_markers_temp = list(set(duplicate_loci))
        paralogous_markers = marker_read_dict
        paralog_info = pd.DataFrame.from_dict(paralogous_markers, orient='index')
        paralog_info.to_csv(os.path.join(outdir,'info_paralogous_loci.txt'),header=False,sep="\t")
        print('Found %i markers with multiple read hits. Selecting one random read per cluster.'%len(invalid_markers_temp))
    else:
        # remove all duplicates
        invalid_markers_unique = list(set(duplicate_loci))
        print(len(invalid_markers_unique), 'possibly paralogous markers detected - excluded from processing')
    # get list of valid read names
    valid_read_names = []
    for marker in marker_read_dict:
        if marker not in invalid_markers_unique:
            read_name = marker_read_dict[marker]
            read_name = np.random.choice(read_name)
            #read_name = find_longest_read(read_name,lastz_df)
            valid_read_names.append(str(read_name).replace('>',''))
    return valid_read_names


def extract_target_reads(sample_id,read_sequences,valid_read_names,read_marker_dict,read_orientation_dict,subfolder):
    printed_reads_counter = 0
    # define the output file where extracted reads will be stored
    global_match_output_name = 'extracted_target_reads_all_samples.fasta'
    global_match_output_file = os.path.join('/'.join(subfolder.split('/')[:-1]),global_match_output_name)
    sample_match_output_name = 'extracted_target_reads_%s.fasta'%sample_id
    sample_match_output_file = os.path.join(subfolder,sample_match_output_name)
    # extract valid reads form read file and print to fasta file with marker-names+ sample_id as headers
    with open(global_match_output_file, "a") as out_file:
        with open(sample_match_output_file, "w") as sample_file:
            for fasta in read_sequences:
                if str(fasta.id) in valid_read_names:
                    orientation = read_orientation_dict[fasta.id]
                    if orientation == '-':
                        seq = fasta.seq.reverse_complement()
                    else:
                        seq = fasta.seq
                    # get the corresponding marker locus name from the dictionary
                    if len(read_marker_dict[fasta.id])>1:
                        for matching_locus in read_marker_dict[fasta.id]:
                            header = '%s_%s |%s' %(matching_locus,sample_id,matching_locus)
                            new_fasta = SeqRecord(seq, id=header, name='', description='')
                            out_file.write(new_fasta.format('fasta'))
                            sample_file.write(new_fasta.format('fasta'))
                            printed_reads_counter += 1                         
                    else:        
                        header = '%s_%s |%s' %(read_marker_dict[fasta.id][0],sample_id,read_marker_dict[fasta.id][0])
                        new_fasta = SeqRecord(seq, id=header, name='', description='')
                        out_file.write(new_fasta.format('fasta'))
                        sample_file.write(new_fasta.format('fasta'))
                        printed_reads_counter += 1 
    return printed_reads_counter


marker_library = '/Users/tobias/GitHub/ssr_extraction/data/processed/marker_sequences.fasta'
min_coverage = 80
min_identity = 80 # 4% is equivalent to 1 allowed mismatches on the 25bp long reference sequence
outfolder = '/Users/tobias/GitHub/ssr_extraction/data/processed/reads_matching_markers/'
if not os.path.exists(outfolder):
    os.makedirs(outfolder)

# find all fastq files
read_dir = '/Users/tobias/GitHub/ssr_extraction/data/raw/other'
fastq_files = glob.glob(os.path.join(read_dir,'*.fq'))

# Create a dataframe filled with 0's
marker_list = [i.name for i in list(AlignIO.read(marker_library,format='fasta'))]
zeros_array = np.zeros([len(marker_list),len(fastq_files)]).astype(int)
column_names = [os.path.basename(i).replace('.fq','') for i in fastq_files]
read_match_df = pd.DataFrame(data=zeros_array, index=marker_list, columns=column_names)

for reads_file in fastq_files:
    sample_name = os.path.basename(reads_file).replace('.fq','')
    subfolder = os.path.join(outfolder,sample_name)
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)    
        
    lastz_output = os.path.join(subfolder,'%s_read_matches.lastz'%sample_name)
    # Blast the the reads against the reference file
    with open(lastz_output, 'w') as lastz_out_file:
        lastz_command = [
            'lastz',
            '%s[multiple,nameparse=full]'%marker_library,
            '%s[nameparse=full]'%reads_file,
            '--strand=both',
#            '--seed=half20',
#            '--seed=match15',
#            '--seed=12of19',
#            '--seed=1111111111111111111111',
#            '--word=15',
#            '--transition',
#            '--nochain',
#            '--gapped',
#            '--noentropy',
            '--coverage=%i'%min_coverage,
            '--identity=%i'%min_identity,
            '--ambiguous=iupac',
            '--format=general:score,name1,strand1,zstart1,end1,length1,name2,strand2,zstart2,end2,length2,diff,cigar,identity,continuity'
        ]
        run_lastz = subprocess.Popen(lastz_command, stdout=lastz_out_file, stderr=None)
        run_lastz.communicate()
    # load the lastz matches from the previous command
    lastz_df = pd.read_csv(lastz_output,sep='\t')
    
    marker_read_dict, read_marker_dict, read_orientation_dict, read_multi_marker_dict, orientation_df = reads_matching_markers(lastz_df)
    


    
    orientation_df.to_csv(os.path.join(subfolder,'read_orientation.txt'),index=False,sep='\t')
    # mark duplicate loci
    duplicate_loci, possible_paralogous, reads_covering_several_loci = find_duplicates(marker_read_dict,read_marker_dict)
    # remove duplicate loci from the list of targeted loci and reads
    target_reads = get_list_of_valid_markers_and_reads_including_all_loci_in_paralog_output_file(marker_read_dict,duplicate_loci,possible_paralogous,read_multi_marker_dict,subfolder,lastz_df)
    # load the actual read sequences
    read_sequences = SeqIO.parse(open(reads_file),'fastq')
    # write those reads that match the reference library to the file
    extracted_read_counter = extract_target_reads(sample_name,read_sequences,target_reads,read_marker_dict,read_orientation_dict,subfolder)
    # Fill the extracted target read into the dataframe
    for read in target_reads:
        for marker in read_marker_dict[read]:
            read_match_df.loc[marker,sample_name] = 1        
    print('Extracted %i reads matching reference markers\n' %extracted_read_counter)

read_match_df.to_csv(os.path.join(outfolder,'match_table.txt'),sep='\t',index=True,encoding='utf-8')

