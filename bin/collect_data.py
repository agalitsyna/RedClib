#!/usr/bin/env python3

"""
Parse intermediary data files and produce single hdf5 file with collected information about each read.
All the intermediary txt files can be deleted after this step.

Usage:
python collect_data.py test output.hdf5 {..long list of input files..} [.. optional files with restriction annot..]

Example:
python bin/collect_data.py tmp.hdf5 \
  results/table/test.001.fastq.txt results/table/test.001.ids.unique.txt results/table/test.001.trimtable.txt \
  results/cout/test.001.1.bridge_forward.txt results/table/test.001.GA.txt results/cout/test.001.2.ggg.txt \
  results/filtered_fastq/test.001.dna.info.txt results/filtered_fastq/test.001.rna1.info.txt results/filtered_fastq/test.001.rna2.info.txt \
  results/sam/test.001.dna.extended.sam results/sam/test.001.dna.sam results/sam/test.001.rna1.sam results/sam/test.001.rna2.sam \
  results/bed/test.001.dna_ext.bed results/bed/test.001.rna1.bed results/bed/test.001.rna2.bed \
  results/table/test.001.rna1.mme-.distances.txt results/table/test.001.rna1.mme+.distances.txt \
  results/table/test.001.rna1.nla.distances.txt results/table/test.001.rna2.mme-.distances.txt \
  results/table/test.001.rna2.mme+.distances.txt results/table/test.001.rna2.nla.distances.txt \
  results/table/test.001.dna.nla.distances.txt
"""

import numpy as np
import h5py
import re
from sys import argv
from utils import *

if len(argv)<18:
    raise Exception("Number of inputs is {len(argv)-1}, not 17 as expected. See usage in the file.")

output_filename = argv[1]

table_file, is_dup_file, trim_file, \
    bridge_file, br_err, ggg_file, \
    dna_len_file, rna1_len_file, rna2_len_file, \
    dna_map_file_sam, dna_map_file_sam_nonextended, rna1_map_file_sam, rna2_map_file_sam, \
    dna_map_file, rna1_map_file, rna2_map_file = argv[2:18]

renz_files = []
if len(argv)>=19:
    renz_files = argv[18:]


outfile = h5py.File(output_filename, "w")

raw_file = raw_read_file(table_file,
                         [1, 3, 5, 3, 5],
                         ['id', 'seqR1_len', 'seqR2_len', 'seqR1', 'seqR2'],
                         modifiers=[lambda x: x[1:], lambda x: len(x), lambda x: len(x), str, str], header=0)
update_hdf5(outfile, raw_file, {'id': 'S100', 'seqR1': 'S250', 'seqR2': 'S250', 'seqR1_len': int, 'seqR2_len': int})

bridge_raw_file = raw_read_file(bridge_file,
                                [3, 4, 5, 6],
                                ['has_nobridge', 'bridge_nmm', 'bridge_start', 'bridge_end'],
                                modifiers=[lambda x: int(x) for y in range(4)], header=1, bottom=1)
update_hdf5(outfile, bridge_raw_file, {'has_nobridge':bool, 'bridge_nmm':int, 'bridge_start':int, 'bridge_end':int})
del bridge_raw_file

br_err_raw_file = raw_read_file(br_err,
                                [2],
                                ['has_GA'],
                                modifiers=[lambda x: int(x)], header=0)
update_hdf5(outfile, br_err_raw_file, bool)
del br_err_raw_file

trim_raw_file = raw_read_file(trim_file,
                              [5, 7],
                              ['trimF', 'trimR'],
                              modifiers=[lambda x: int(x) for y in range(2)], header=0)
update_hdf5(outfile, trim_raw_file, int)
del trim_raw_file

ggg_raw_file = raw_read_file(ggg_file,
                             [3, 5, 6],
                             ['has_noggg', 'ggg_start', 'ggg_end'],
                             modifiers=[lambda x: int(x) for y in range(3)], header=1, bottom=1)
update_hdf5(outfile, ggg_raw_file, {'has_noggg':bool, 'ggg_start':int, 'ggg_end':int})
del ggg_raw_file

dna_len_raw_file = raw_read_file(dna_len_file,
                                 [2, 3, 4, 5, 6],
                                 ['dna_R1_start', 'dna_R1_end', 'dna_R1_len_notrim', 'dna_R1_end_trim',
                                  'dna_R1_len_trim'],
                                 modifiers=[lambda x: int(x) for y in range(5)], header=0)
update_hdf5(outfile, dna_len_raw_file, int)
del dna_len_raw_file

rna2_len_raw_file = raw_read_file(rna2_len_file,
                                 [2, 3, 4, 5, 6],
                                 ['rna2_R2_start', 'rna2_R2_end', 'rna2_R2_len_notrim', 'rna2_R2_end_trim',
                                  'rna2_R2_len_trim'],
                                 modifiers=[lambda x: int(x) for y in range(5)], header=0)
update_hdf5(outfile, rna2_len_raw_file, int)
del rna2_len_raw_file

rna1_len_raw_file = raw_read_file(rna1_len_file,
                                  [2, 3, 4, 5, 6],
                                  ['rna1_R1_start', 'rna1_R1_end', 'rna1_R1_len_notrim', 'rna1_R1_end_trim',
                                   'rna1_R1_len_trim'],
                                  modifiers=[(lambda x: int(x)) for y in range(5)], header=0)
update_hdf5(outfile, rna1_len_raw_file, int)
del rna1_len_raw_file

dct_ids = {k: 0 for k in raw_file['id']}
with open(is_dup_file, 'r') as inf:
    l = inf.readline().strip()
    while len(l) > 0:
        l = inf.readline().strip()
        # Skip if fastuniq output is for larger file and contains extra ids
        # (dedup is done on full library, not the chunks)
        try:
            dct_ids[l] = 1
        except Exception as e:
            pass

fastuniq_dct = {'is_notPCRdup': [dct_ids[l] for l in raw_file['id']]}
update_hdf5(outfile, fastuniq_dct, int)
del fastuniq_dct, dct_ids

dna_map_raw_file = raw_read_file(dna_map_file_sam,
                                 [1, 2, -1, 3, 6],
                                 ['id', 'dna_is_mapped', 'dna_is_not_multi', 'dna_chr', 'dna_cigar'],
                                 modifiers=[str, lambda x: 0 if x == '4' else 1, lambda x: 1 if x == 'NH:i:1' else 0,
                                            str, str], header=0, comment="@")

dna_map_raw_file_inferred = \
    reconstruct_by_ids(dna_map_raw_file, 'id', raw_file['id'],
                       default_dct={'dna_is_mapped': 0, 'dna_is_not_multi': 0, 'dna_chr': '-', 'dna_cigar': ''})

del dna_map_raw_file_inferred['id']
update_hdf5(outfile, dna_map_raw_file_inferred,
            {'dna_is_mapped': bool, 'dna_is_not_multi': bool, 'dna_chr': 'S8', 'dna_cigar': 'S20'})
outfile.create_dataset('dna_nlen', data=np.array( \
    [np.sum([int(m[0]) if m[1] == 'N' else 0 for m in re.findall(r'(\d+)([A-Z]{1})', x)]) for x in
     dna_map_raw_file_inferred['dna_cigar']]))
del dna_map_raw_file_inferred, dna_map_raw_file

rna2_map_raw_file = raw_read_file(rna2_map_file_sam,
                                 [1, 2, -1, 3, 6],
                                 ['id', 'rna2_is_mapped', 'rna2_is_not_multi', 'rna2_chr', 'rna2_cigar'],
                                 modifiers=[str, lambda x: 0 if x == '4' else 1, lambda x: 1 if x == 'NH:i:1' else 0,
                                            str, str], header=0, comment="@")

rna2_map_raw_file_inferred = \
    reconstruct_by_ids(rna2_map_raw_file, 'id', raw_file['id'],
                       default_dct={'rna2_is_mapped': 0, 'rna2_is_not_multi': 0, 'rna2_chr': '-', 'rna2_start': -1,
                                    'rna2_cigar': ''})
del rna2_map_raw_file_inferred['id']
update_hdf5(outfile, rna2_map_raw_file_inferred,
            {'rna2_is_mapped': bool, 'rna2_is_not_multi': bool, 'rna2_chr': 'S8', 'rna2_cigar': 'S20'})
outfile.create_dataset('rna2_nlen', data=np.array( \
    [np.sum([int(m[0]) if m[1] == 'N' else 0 for m in re.findall(r'(\d+)([A-Z]{1})', x)]) for x in
     rna2_map_raw_file_inferred['rna2_cigar']]))

del rna2_map_raw_file_inferred, rna2_map_raw_file

rna1_map_raw_file = raw_read_file(rna1_map_file_sam,
                                  [1, 2, -1, 3, 6],
                                  ['id', 'rna1_is_mapped', 'rna1_is_not_multi', 'rna1_chr', 'rna1_cigar'],
                                  modifiers=[str, lambda x: 0 if x == '4' else 1, lambda x: 1 if x == 'NH:i:1' else 0,
                                             str, str], header=0, comment="@")

rna1_map_raw_file_inferred = \
    reconstruct_by_ids(rna1_map_raw_file, 'id', raw_file['id'],
                       default_dct={'rna1_is_mapped': 0, 'rna1_is_not_multi': 0, 'rna1_chr': '-', 'rna1_cigar': ''})
del rna1_map_raw_file_inferred['id']
update_hdf5(outfile, rna1_map_raw_file_inferred,
            {'rna1_is_mapped': bool, 'rna1_is_not_multi': bool, 'rna1_chr': 'S8', 'rna1_cigar': 'S20'})
outfile.create_dataset('rna1_nlen', data=np.array(
    [np.sum([int(m[0]) if m[1] == 'N' else 0 for m in re.findall(r'(\d+)([A-Z]{1})', x)]) for x in
     rna1_map_raw_file_inferred['rna1_cigar']]))
del rna1_map_raw_file_inferred, rna1_map_raw_file

dna_map_raw_file_nonextended = raw_read_file(dna_map_file_sam_nonextended,
                                 [1, 2, -1, 3, 6],
                                 ['id', 'dna_nonextended_is_mapped', 'dna_nonextended_is_not_multi', 'dna_nonextended_chr',
                                  'dna_nonextended_cigar'],
                                 modifiers=[str, lambda x: 0 if x == '4' else 1, lambda x: 1 if x == 'NH:i:1' else 0,
                                            str, str], header=0, comment="@")

dna_map_raw_file_inferred = \
    reconstruct_by_ids(dna_map_raw_file_nonextended, 'id', raw_file['id'],
                       default_dct={'dna_nonextended_is_mapped': 0, 'dna_nonextended_is_not_multi': 0, 'dna_nonextended_chr': '-',
                                    'dna_nonextended_cigar': ''})

del dna_map_raw_file_inferred['id']
update_hdf5(outfile, dna_map_raw_file_inferred,
            {'dna_nonextended_is_mapped': bool, 'dna_nonextended_is_not_multi': bool, 'dna_nonextended_chr': 'S8',
             'dna_nonextended_cigar': 'S20'})
del dna_map_raw_file_inferred, dna_map_raw_file_nonextended

# Reading extended bed file:
dna_map_raw_file = raw_read_file(dna_map_file,
                                 [2, 3, 4, 6],
                                 ['dna_start', 'dna_end', 'id', 'dna_strand'],
                                 modifiers=[int, int, str, lambda x: 1 if x == '+' else 0], header=0, comment="@")

dna_map_raw_file_inferred = \
    reconstruct_by_ids(dna_map_raw_file, 'id', raw_file['id'],
                       default_dct={'dna_start': 0, 'dna_end': 0, 'dna_strand': 0})

del dna_map_raw_file_inferred['id']
update_hdf5(outfile, dna_map_raw_file_inferred, {'dna_start': int, 'dna_end': int, 'dna_strand': bool})
del dna_map_raw_file_inferred, dna_map_raw_file

rna2_map_raw_file = raw_read_file(rna2_map_file,
                                 [2, 3, 4, 6],
                                 ['rna2_start', 'rna2_end', 'id', 'rna2_strand'],
                                 modifiers=[int, int, str, lambda x: 1 if x == '+' else 0], header=0, comment="@")

rna2_map_raw_file_inferred = \
    reconstruct_by_ids(rna2_map_raw_file, 'id', raw_file['id'],
                       default_dct={'rna2_start': 0, 'rna2_end': 0, 'rna2_strand': 0})

del rna2_map_raw_file_inferred['id']
update_hdf5(outfile, rna2_map_raw_file_inferred, {'rna2_start': int, 'rna2_end': int, 'rna2_strand': bool})
del rna2_map_raw_file_inferred, rna2_map_raw_file

rna1_map_raw_file = raw_read_file(rna1_map_file,
                                  [2, 3, 4, 6],
                                  ['rna1_start', 'rna1_end', 'id', 'rna1_strand'],
                                  modifiers=[int, int, str, lambda x: 1 if x == '+' else 0], header=0, comment="@")

rna1_map_raw_file_inferred = \
    reconstruct_by_ids(rna1_map_raw_file, 'id', raw_file['id'],
                       default_dct={'rna1_start': 0, 'rna1_end': 0, 'rna1_strand': 0})

del rna1_map_raw_file_inferred['id']
update_hdf5(outfile, rna1_map_raw_file_inferred, {'rna1_start': int, 'rna1_end': int, 'rna1_strand': bool})
del rna1_map_raw_file_inferred, rna1_map_raw_file

# Adjust the trimmed length by oligos found in read:
dna_len_trim  = np.minimum(outfile['dna_R1_len_trim'][()], outfile['dna_R1_len_notrim'][()] )
rna2_len_trim = np.minimum( outfile['rna2_R2_len_trim'][()], outfile['rna2_R2_len_notrim'][()] )
rna1_len_trim = np.minimum( outfile['rna1_R1_len_trim'][()], outfile['rna1_R1_len_notrim'][()] )
outfile.create_dataset('dna_R1_len_trim_adjusted', data=np.array(dna_len_trim, dtype=int))
outfile.create_dataset('rna2_R2_len_trim_adjusted', data=np.array(rna2_len_trim, dtype=int))
outfile.create_dataset('rna1_R1_len_trim_adjusted', data=np.array(rna1_len_trim, dtype=int))


### Filtering by restiriction enzyme recognition site:
renz_keys = []
for renz_file in renz_files:
    header = open(renz_file, 'r').readline().strip().split()
    renz_keys += header[1:]
    renz_raw_file = raw_read_file(renz_file,
                                      [1, 2, 3, 4, 5],
                                      header,
                                      sep=' ',
                                      modifiers=[str]+[int for x in range(4)],
                                      header=1)

    renz_raw_file_inferred = \
        reconstruct_by_ids(renz_raw_file, 'id', raw_file['id'],
                           default_dct={x: 0 for x in header[1:]})

    del renz_raw_file_inferred['id']
    update_hdf5(outfile, renz_raw_file_inferred, {x: int for x in header[1:]})
    del renz_raw_file_inferred, renz_raw_file

outfile.close()
