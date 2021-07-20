#!/usr/bin/env python3

"""
Detect the START POSITIONS of recognition sites of restriction enzymes.
This is not the same as restriction sites, see http://biopython.org/DIST/docs/cookbook/Restriction.html#mozTocId447698
Also report the strand.

End of restriction site is formal (same as start).

Usage:
python detect_restriction_sites hg19.fa MmeI output.bed
"""

from sys import argv
import numpy as np
import pandas as pd

import Bio.Restriction as biorst
from Bio.SeqIO import parse

input_genome = argv[1]
enzyme_name= argv[2]
output_file = argv[3]

fasta_records = parse(input_genome, "fasta")
enzyme = getattr(biorst, enzyme_name)

rsites = []
for seq_record in fasta_records:
    chrom = seq_record.id
    enzyme.search(seq_record.seq)
    if enzyme.is_palindromic():
        rsites.append(
            pd.DataFrame({'chrom': chrom, 'site': enzyme.results, 'strand': '+'})
        )
    else:
        minus = enzyme.on_minus
        rsites.append(
            pd.DataFrame({'chrom': chrom, 'site': minus, 'strand': '-'})
        )
        plus = np.setdiff1d(enzyme.results, minus)
        rsites.append(
            pd.DataFrame({'chrom': chrom, 'site': plus, 'strand': '+'})
        )

rsites = pd.concat(rsites).sort_values(['chrom', 'site']).reset_index(drop=True)

def detect_recognition_start(row):
    if row.strand == "+" :
        return row.site - enzyme.fst5 - 1
    else:
        return row.site + enzyme.fst3 - 1

# Retrieving actual position of recognition site
idx_pos = (rsites.strand=="+")
rsites.loc[idx_pos,'start'] = rsites.loc[idx_pos,'site'] - enzyme.fst5 - 1
idx_neg = (rsites.strand=="-")
rsites.loc[idx_neg,'start'] = rsites.loc[idx_neg,'site'] + enzyme.fst3 - 1

rsites.loc[:, 'end'] = rsites.start.astype(int)
rsites.loc[:, 'name'] = [f'{enzyme_name}_{idx+1}' for idx in rsites.index]
rsites.loc[:, 'foo'] = '.'

rsites.to_csv(output_file,
              sep='\t',
              columns=['chrom', 'start', 'end', 'name', 'foo', 'strand'],
              header=False,
              index=False)
