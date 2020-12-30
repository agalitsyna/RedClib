#!/usr/bin/env python3

"""
Get distances to the closest sites to the start and end of mapped reads.
Output file format: tsv file with fields:
read id,
distance to the closest site (strand-specific) to the read start from the left,
distance to the closest site to the read start from the right,
distance to the closest site to the read end from the left,
distance to the closest site to the read end from the right.

If the read is located before the first site or after the last site in the chromosome,
then the reported distance will be artificially large (>> 1e9).

Usage:
python get_closest_sites.py input_reads.bed input_rsites.bed strand output.txt

Example:
python get_closest_sites.py ./AG66_2.001.dna.bed ./dm3.DpnI.bed + tmp.txt

"""

import numpy as np
import pandas as pd
from sys import argv

if not len(argv)==5:
    raise Exception("Number of arguments is too small")

filename_bed    = argv[1]
filename_rsites = argv[2]
strand  = argv[3]
outfile = argv[4]

rsites = pd.read_csv(filename_rsites, header=None, sep='\t')
rsites.columns = ['chrom', 'start', 'end', 'name', '_', 'strand']
rsites = rsites.loc[rsites.strand == strand, :].sort_values(['chrom', 'start'])
rsites_grouped = rsites.groupby('chrom')

bed_file = pd.read_csv(filename_bed, sep='\t', header=None)
bed_file.columns = ['chrom', 'start', 'end', 'read_id', 'q', 'strand', 'cigar']

chs = bed_file.chrom.values.astype(str)
ids = bed_file.read_id.values.astype(str)
bgn = bed_file.start.values.astype(int)
end = bed_file.end.values.astype(int)
l = len(chs)

dct = {}
for k in ['start_left', 'start_right', 'end_left', 'end_right']:
    dct[k] = np.full(l, -1).astype(int)

for ch in rsites_grouped.groups.keys():
    rs = np.concatenate( [[-1e10], rsites_grouped.get_group(ch)['start'].values, [1e10]] )
    mask = chs == ch

    idx = np.digitize(bgn[mask], rs)
    bgns = rs[idx - 1]
    ends = rs[idx]
    dct['start_left'][mask] = bgns - bgn[mask]
    dct['start_right'][mask] = ends - bgn[mask]

    idx = np.digitize(end[mask], rs)
    bgns = rs[idx - 1]
    ends = rs[idx]
    dct['end_left'][mask] = bgns - end[mask]
    dct['end_right'][mask] = ends - end[mask]

with open(outfile, 'a') as outf:
    for i in range(l):
        outf.write(
            ids[i] + " " + \
            str(dct['start_left'][i]) + " " + \
            str(dct['start_right'][i]) + " " + \
            str(dct['end_left'][i]) + " " + \
            str(dct['end_right'][i]) + "\n")
