"""
Parse all the intermediary data files and produce single hdf5 file with collected information about each read.
All the intermediary txt files can be deleted after this step.
"""

from RedClib import *

import pandas as pd
import numpy as np
import h5py
import re
from sys import argv

prefix = argv[1]

# Change of the relative path to the absolute path
PATH_ABSOLUTE = os.path.dirname(os.path.realpath(__file__))
PATH_OUT = os.path.join(PATH_ABSOLUTE, "data/out/")

if not os.path.isdir(PATH_OUT):
    os.mkdir(PATH_OUT)

filename = os.path.join(PATH_OUT, f"table_{prefix}.hdf5")

try:
    outfile = h5py.File(filename, "w")
    
except Exception as e:
    import os
    os.remove(filename)
    logging.error(e)
    exit()
    

def update_hdf5(file_hdf5, dct, dtype=int):
    """dct might be either dictionary or a single data type"""

    ks = dct.keys()
    if not isinstance(dtype, dict):
        dtype = {k:dtype for k in ks}

    for k in ks:
        file_hdf5.create_dataset(k, data=np.array(dct[k], dtype=dtype[k]), chunks=True)


table_file = os.path.join(PATH_ABSOLUTE, f"data/tables/{prefix}.fastq.txt")
raw_file = raw_read_file(table_file,
                         [1, 3, 5, 3, 5],
                         ['id', 'seqR1_len', 'seqR2_len', 'seqR1', 'seqR2'],
                         modifiers=[lambda x: x[1:], lambda x: len(x), lambda x: len(x), str, str], header=0)

update_hdf5(outfile, raw_file, {'id':'S100', 'seqR1':'S250', 'seqR2':'S250', 'seqR1_len':int, 'seqR2_len':int})

bridge_file = os.path.join(PATH_ABSOLUTE, f"data/cout/{prefix}_R1.37br_for.txt")
bridge_raw_file = raw_read_file(bridge_file,
                                [3, 4, 5, 6],
                                ['has_nobridge', 'bridge_nmm', 'bridge_bgn', 'bridge_end'],
                                modifiers=[lambda x: int(x) for y in range(4)], header=1, bottom=1)
update_hdf5(outfile, bridge_raw_file, int)
del bridge_raw_file


br_err = os.path.join(PATH_ABSOLUTE, f"data/cout/{prefix}_R1.GA.txt")
br_err_raw_file = raw_read_file(br_err,
                                [2],
                                ['has_GA'],
                                modifiers=[lambda x: int(x)], header=0)
update_hdf5(outfile, br_err_raw_file, int)
del br_err_raw_file

trim_file = os.path.join(PATH_ABSOLUTE, f"data/tables/{prefix}.trimtable.txt")
trim_raw_file = raw_read_file(trim_file,
                              [5, 7],
                              ['trimF', 'trimR'],
                              modifiers=[lambda x: int(x) for y in range(2)], header=0)
update_hdf5(outfile, trim_raw_file, int)
del trim_raw_file

ggg_file = os.path.join(PATH_ABSOLUTE, f"data/cout/{prefix}_R2.ggg.txt")
ggg_raw_file = raw_read_file(ggg_file,
                             [3, 5, 6],
                             ['has_noggg', 'ggg_bgn', 'ggg_end'],
                             modifiers=[lambda x: int(x) for y in range(3)], header=1, bottom=1)
update_hdf5(outfile, ggg_raw_file, int)
del ggg_raw_file

dna_len_file = os.path.join(PATH_ABSOLUTE, f"data/filtered/{prefix}.dna.info.txt")
dna_len_raw_file = raw_read_file(dna_len_file,
                                 [2, 3, 4, 5, 6],
                                 ['dna_R1_bgn', 'dna_R1_end', 'dna_R1_len_notrim', 'dna_R1_end_trim', 'dna_R1_len_trim'],
                                 modifiers=[lambda x: int(x) for y in range(5)], header=0)
update_hdf5(outfile, dna_len_raw_file, int)
del dna_len_raw_file

rna_len_file = os.path.join(PATH_ABSOLUTE, f"data/filtered/{prefix}.rna.info.txt")
rna_len_raw_file = raw_read_file(rna_len_file,
                                 [2, 3, 4, 5, 6],
                                 ['rna_R2_bgn', 'rna_R2_end', 'rna_R2_len_notrim', 'rna_R2_end_trim', 'rna_R2_len_trim'],
                                 modifiers=[lambda x: int(x) for y in range(5)], header=0)
update_hdf5(outfile, rna_len_raw_file, int)
del rna_len_raw_file

rna1_len_file = os.path.join(PATH_ABSOLUTE, f"data/filtered/{prefix}.rna1.info.txt")
rna1_len_raw_file = raw_read_file(rna1_len_file,
                                  [2, 3, 4, 5, 6],
                                  ['rna1_R1_bgn', 'rna1_R1_end', 'rna1_R1_len_notrim', 'rna1_R1_end_trim', 'rna1_R1_len_trim'],
                                  modifiers=[(lambda x: int(x)) for y in range(5)], header=0)
update_hdf5(outfile, rna1_len_raw_file, int)
del rna1_len_raw_file


dct_ids = {k:0 for k in raw_file['id']}
is_dup_file = os.path.join(PATH_ABSOLUTE, f"data/tables/{prefix}.unique_idx.txt")
with open(is_dup_file, 'r') as inf:
    l = inf.readline().strip()
    while len(l)>0:
        l = inf.readline().strip()
        dct_ids[l] = 1

fastuniq_dct = {'is_notPCRdup':[dct_ids[l] for l in raw_file['id']]}
update_hdf5(outfile, fastuniq_dct, int)
del fastuniq_dct, dct_ids

dna_map_file = os.path.join(PATH_ABSOLUTE, f"data/sam/{prefix}.dna.sam")
dna_map_raw_file = raw_read_file(dna_map_file,
  [1, 2, -1, 3, 6],
  ['id', 'dna_is_mapped', 'dna_is_not_multi', 'dna_chr', 'dna_cigar'],
  modifiers=[str, lambda x: 0 if x=='4' else 1, lambda x: 1 if x=='NH:i:1' else 0, str, str], header=0, comment="@")

dna_map_raw_file_inferred = \
  reconstruct_by_ids(dna_map_raw_file, 'id', raw_file['id'],
                     default_dct={'dna_is_mapped':0, 'dna_is_not_multi':0, 'dna_chr': '-','dna_cigar':''})

del dna_map_raw_file_inferred['id']
update_hdf5(outfile, dna_map_raw_file_inferred, {'dna_is_mapped':int, 'dna_is_not_multi':int, 'dna_chr':'S8', 'dna_cigar':'S20'})
outfile.create_dataset('dna_nlen', data=np.array(\
    [ np.sum([int(m[0]) if m[1]=='N' else 0 for m in re.findall(r'(\d+)([A-Z]{1})', x)]) for x in dna_map_raw_file_inferred['dna_cigar'] ]))
del dna_map_raw_file_inferred, dna_map_raw_file


rna_map_file = os.path.join(PATH_ABSOLUTE, f"data/sam/{prefix}.rna.sam")
rna_map_raw_file = raw_read_file(rna_map_file,
  [1, 2, -1, 3, 6],
  ['id', 'rna_is_mapped', 'rna_is_not_multi', 'rna_chr', 'rna_cigar'],
  modifiers=[str, lambda x: 0 if x=='4' else 1, lambda x: 1 if x=='NH:i:1' else 0, str, str], header=0, comment="@")

rna_map_raw_file_inferred = \
  reconstruct_by_ids(rna_map_raw_file, 'id', raw_file['id'],
                     default_dct={'rna_is_mapped':0, 'rna_is_not_multi':0, 'rna_chr': '-', 'rna_bgn': -1, 'rna_cigar':''})
del rna_map_raw_file_inferred['id']
update_hdf5(outfile, rna_map_raw_file_inferred, {'rna_is_mapped':int, 'rna_is_not_multi':int, 'rna_chr':'S8', 'rna_cigar':'S20'})
outfile.create_dataset('rna_nlen', data=np.array(\
    [ np.sum([int(m[0]) if m[1]=='N' else 0 for m in re.findall(r'(\d+)([A-Z]{1})', x)]) for x in rna_map_raw_file_inferred['rna_cigar'] ]))

del rna_map_raw_file_inferred, rna_map_raw_file

rna1_map_file = os.path.join(PATH_ABSOLUTE, f"data/sam/{prefix}.rna1.sam")
rna1_map_raw_file = raw_read_file(rna1_map_file,
  [1, 2, -1, 3, 6],
  ['id', 'rna1_is_mapped', 'rna1_is_not_multi', 'rna1_chr', 'rna1_cigar'],
  modifiers=[str, lambda x: 0 if x=='4' else 1, lambda x: 1 if x=='NH:i:1' else 0, str, str], header=0, comment="@")

rna1_map_raw_file_inferred = \
  reconstruct_by_ids(rna1_map_raw_file, 'id', raw_file['id'],
                     default_dct={'rna1_is_mapped':0, 'rna1_is_not_multi':0, 'rna1_chr': '-', 'rna1_cigar':''})
del rna1_map_raw_file_inferred['id']
update_hdf5(outfile, rna1_map_raw_file_inferred, {'rna1_is_mapped':int, 'rna1_is_not_multi':int, 'rna1_chr':'S8', 'rna1_cigar':'S20'})
outfile.create_dataset('rna1_nlen', data=np.array([ np.sum([int(m[0]) if m[1]=='N' else 0 for m in re.findall(r'(\d+)([A-Z]{1})', x)]) for x in rna1_map_raw_file_inferred['rna1_cigar'] ]))
del rna1_map_raw_file_inferred, rna1_map_raw_file


dna_map_file = os.path.join(PATH_ABSOLUTE, f"data/sam/{prefix}.dna.sam.nonla")
dna_map_raw_file = raw_read_file(dna_map_file,
  [1, 2, -1, 3, 6],
  ['id', 'dna_noNla_is_mapped', 'dna_noNla_is_not_multi', 'dna_noNla_chr', 'dna_noNla_cigar'],
  modifiers=[str, lambda x: 0 if x=='4' else 1, lambda x: 1 if x=='NH:i:1' else 0, str, str], header=0, comment="@")

dna_map_raw_file_inferred = \
  reconstruct_by_ids(dna_map_raw_file, 'id', raw_file['id'],
                     default_dct={'dna_noNla_is_mapped':0, 'dna_noNla_is_not_multi':0, 'dna_noNla_chr': '-','dna_noNla_cigar':''})

del dna_map_raw_file_inferred['id']
update_hdf5(outfile, dna_map_raw_file_inferred, {'dna_noNla_is_mapped':int, 'dna_noNla_is_not_multi':int, 'dna_noNla_chr':'S8', 'dna_noNla_cigar':'S20'})
del dna_map_raw_file_inferred, dna_map_raw_file


dna_map_file = os.path.join(PATH_ABSOLUTE, f"data/bed/{prefix}.dna.bed")
dna_map_raw_file = raw_read_file(dna_map_file,
  [2, 3, 4, 6],
  ['dna_bgn', 'dna_end', 'id', 'dna_strand'],
  modifiers=[int, int, str, lambda x: 1 if x=='+' else 0], header=0, comment="@")

dna_map_raw_file_inferred = \
  reconstruct_by_ids(dna_map_raw_file, 'id', raw_file['id'],
                     default_dct={'dna_bgn':0, 'dna_end':0, 'dna_strand': 0})

del dna_map_raw_file_inferred['id']
update_hdf5(outfile, dna_map_raw_file_inferred, {'dna_bgn':int, 'dna_end':int, 'dna_strand':bool})
del dna_map_raw_file_inferred, dna_map_raw_file


rna_map_file = os.path.join(PATH_ABSOLUTE, f"data/bed/{prefix}.rna.bed")
rna_map_raw_file = raw_read_file(rna_map_file,
  [2, 3, 4, 6],
  ['rna_bgn', 'rna_end', 'id', 'rna_strand'],
  modifiers=[int, int, str, lambda x: 1 if x=='+' else 0], header=0, comment="@")

rna_map_raw_file_inferred = \
  reconstruct_by_ids(rna_map_raw_file, 'id', raw_file['id'],
                     default_dct={'rna_bgn':0, 'rna_end':0, 'rna_strand': 0})

del rna_map_raw_file_inferred['id']
update_hdf5(outfile, rna_map_raw_file_inferred, {'rna_bgn':int, 'rna_end':int, 'rna_strand':bool})
del rna_map_raw_file_inferred, rna_map_raw_file


rna1_map_file = os.path.join(PATH_ABSOLUTE, f"data/bed/{prefix}.rna1.bed")
rna1_map_raw_file = raw_read_file(rna1_map_file,
  [2, 3, 4, 6],
  ['rna1_bgn', 'rna1_end', 'id', 'rna1_strand'],
  modifiers=[int, int, str, lambda x: 1 if x=='+' else 0], header=0, comment="@")

rna1_map_raw_file_inferred = \
  reconstruct_by_ids(rna1_map_raw_file, 'id', raw_file['id'],
                     default_dct={'rna1_bgn':0, 'rna1_end':0, 'rna1_strand': 0})

del rna1_map_raw_file_inferred['id']
update_hdf5(outfile, rna1_map_raw_file_inferred, {'rna1_bgn':int, 'rna1_end':int, 'rna1_strand':bool})
del rna1_map_raw_file_inferred, rna1_map_raw_file


### Reading already prepared dataset with restriction fragments
f_rsites = h5py.File(os.path.join(PATH_ABSOLUTE, "data/genome/rsites.hdf5"), "r")


### Preparation of dataset for checking rfrag positions
dna_bgn = outfile['dna_bgn'][()]
dna_end = outfile['dna_end'][()]
dna_chr = outfile['dna_chr'][()]
rna_bgn = outfile['rna_bgn'][()]
rna_end = outfile['rna_end'][()]
rna_chr = outfile['rna_chr'][()]
rna1_bgn = outfile['rna1_bgn'][()]
rna1_end = outfile['rna1_end'][()]
rna1_chr = outfile['rna1_chr'][()]

l = len(outfile['id'][()])

rng = ['rna', 'rna1', 'dna']
dct = {}
for k in rng:
    for s in ['{}_bgn_dpn_left',  '{}_bgn_dpn_right',  '{}_end_dpn_left',  '{}_end_dpn_right',
              '{}_bgn_nla_left',  '{}_bgn_nla_right',  '{}_end_nla_left',  '{}_end_nla_right',
              '{}_bgn_mmep_left', '{}_bgn_mmep_right', '{}_end_mmep_left', '{}_end_mmep_right',
              '{}_bgn_mmen_left', '{}_bgn_mmen_right', '{}_end_mmen_left', '{}_end_mmen_right']:
        dct[s.format(k)] = np.full(l, np.nan)

rng = [('rna1', rna1_bgn, rna1_end, rna1_chr),
       ('rna', rna_bgn, rna_end, rna_chr),
       ('dna', dna_bgn, dna_end, dna_chr)]

for renz_raw in f_rsites.keys():
    renz = renz_raw.replace('+','p').replace('-','n')
    for ch in f_rsites['nla'].keys():
        if ch=='chrM': continue
        rs = f_rsites[renz_raw][ch][()]
        for k, bgn, end, chs in rng:
            mask = chs==ch.encode()
            idx = np.digitize(bgn[mask], rs)
            bgns = rs[idx-1]
            ends = rs[idx]
            dct['{}_bgn_{}_left'.format(k, renz)][mask] = bgns - bgn[mask]
            dct['{}_bgn_{}_right'.format(k, renz)][mask] = ends - bgn[mask]
            idx = np.digitize(end[mask], rs)
            bgns = rs[idx-1]
            ends = rs[idx]
            dct['{}_end_{}_left'.format(k, renz)][mask] = bgns - end[mask]
            dct['{}_end_{}_right'.format(k, renz)][mask] = ends - end[mask]

for k in dct.keys():
    dct[k][np.isnan(dct[k])] = -1
    dct[k] = np.array(dct[k], dtype=int)

update_hdf5(outfile, dct)
del dct


v = ( outfile['rna_chr'][()]==outfile['rna1_chr'][()] )&( outfile['rna_strand'][()]!=outfile['rna1_strand'][()] )
outfile.create_dataset('RNAsDirectionPassed', data=np.array(v, dtype=bool))

v = ( outfile['dna_chr'][()]==outfile['rna_chr'][()] )&( \
     (((outfile['dna_strand'][()]==1)&(outfile['rna_strand'][()]==0)&(np.abs(outfile['dna_bgn'][()]-outfile['rna_end'][()])<=1))\
     |((outfile['dna_strand'][()]==0)&(outfile['rna_strand'][()]==1)&(np.abs(outfile['dna_end'][()]-outfile['rna_bgn'][()])<=1)) ))
outfile.create_dataset('RNADNASamePos', data=np.array(v, dtype=bool))

chrms = [str(x).encode() for x in pd.read_csv("./data/genome/chromSizes.txt", sep='\t', header=None)[0].values[:-1]]
v = np.in1d(outfile['dna_chr'][()], chrms)
outfile.create_dataset("dna_chr_canonical", data=np.array(v, dtype=bool))

v = np.in1d(outfile['rna_chr'][()], chrms)
outfile.create_dataset("rna_chr_canonical", data=np.array(v, dtype=bool))

v = np.in1d(outfile['rna1_chr'][()], chrms)
outfile.create_dataset("rna1_chr_canonical", data=np.array(v, dtype=bool))


print(outfile['rna_bgn_nla_left'][()][outfile[ 'rna_chr_canonical'][()] ], outfile['rna_strand'][()][outfile[ 'rna_chr_canonical'][()] ])


nla_rna_failed = ((outfile['rna_strand'][()]==1) & \
  (((outfile['rna_bgn_nla_left'][()]<=0) & (outfile['rna_bgn_nla_left'][()]>=-5)) | \
   ((outfile['rna_bgn_nla_right'][()]<=1) & (outfile['rna_bgn_nla_right'][()]>=0)))) | \
((outfile['rna_strand'][()]==0) & \
  (((outfile['rna_end_nla_left'][()]<=0) & (outfile['rna_end_nla_left'][()]>=-5)) | \
   ((outfile['rna_end_nla_right'][()]<=1) & (outfile['rna_end_nla_right'][()]>=0))))

nla_rna1_failed = ((outfile['rna1_strand'][()]==1) & \
  (((outfile['rna1_bgn_nla_left'][()]<=0) & (outfile['rna1_bgn_nla_left'][()]>=-5)) | \
   ((outfile['rna1_bgn_nla_right'][()]<=1) & (outfile['rna1_bgn_nla_right'][()]>=0)))) | \
((outfile['rna1_strand'][()]==0) & \
  (((outfile['rna1_end_nla_left'][()]<=0) & (outfile['rna1_end_nla_left'][()]>=-5)) | \
   ((outfile['rna1_end_nla_right'][()]<=1) & (outfile['rna1_end_nla_right'][()]>=0))))

mme_rna_failed = ((outfile['rna_strand'][()]==1) & \
  ((outfile['rna_bgn_mmep_left'][()] == -26) | \
   (outfile['rna_bgn_mmep_left'][()] == -27) | \
   (outfile['rna_bgn_mmen_right'][()] == 18) | \
   (outfile['rna_bgn_mmen_right'][()] == 19))) | \
 ((outfile['rna_strand'][()]==0) & \
  ((outfile['rna_end_mmep_left'][()] == -25) | \
   (outfile['rna_end_mmep_left'][()] == -26) | \
   (outfile['rna_end_mmen_right'][()] == 21) | \
   (outfile['rna_end_mmen_right'][()] == 20)))


outfile.create_dataset("rna_nla_failed", data=np.array(nla_rna_failed, dtype=bool))
outfile.create_dataset("rna1_nla_failed", data=np.array(nla_rna1_failed, dtype=bool))
outfile.create_dataset("rna_mme_failed", data=np.array(mme_rna_failed, dtype=bool))

# # Filter for + part:
# -5  CATGN|->N
# -4   CATG|->NN
# -3    CAT|->GNN
# -2     CA|->TGNN
# -1      C|->ATGNN
#  0       |->CATGNN
#  1       |->NCATGN
# # Filter for - part:
# -5       |->NCATGN
# -4       |->CATGNN
# -3      C|->ATGNN
# -2     CA|->TGNN
# -1    CAT|->GNN
#  0   CATG|->NN
#  1  CATGN|->N

outfile.close()

# except Exception as e:
#     print("Failed: {}, {}".format(prefix, e))
#     outfile.close()
#     import os
#     os.remove(filename)