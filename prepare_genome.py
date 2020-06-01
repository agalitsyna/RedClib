from mirnylib import genome

import time
from datetime import timedelta

import logging
logging.basicConfig(level=logging.DEBUG)

import h5py
import os
import numpy as np

### Reading genome
logging.info("Preparation of genome...")
cmd_bgn_time = time.time()

PATH_ABSOLUTE = os.path.dirname(os.path.realpath(__file__))
fasta_path   = os.path.join(PATH_ABSOLUTE, "data/genome/hg19.fa")
chrms             = ['10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '1', '20', '21', '22', '2', '3', '4', '5', '6', '7', '8', '9', 'M', 'X', 'Y' ]
genome_db    = genome.Genome(fasta_path, readChrms=chrms)
enzymes = ['NlaIII', 'MmeI', 'DpnII', 'MseI']

# Applying restriction enzymes
genome_db.setEnzyme('MmeI')

strands = [[] for i in range(len(genome_db.seqs))]
true_poss = [[] for i in range(len(genome_db.seqs))]
for i in range(len(genome_db.seqs)):
    s = genome_db.seqs[i]
    for pos in genome_db.rsites[i]:
        print('bad case! {} {}, sequences: {} {}'.format(i, pos, s.seq[pos+16:pos+22], s.seq[pos-28:pos-22]))
        if (str(s.seq[pos+16:pos+18]).upper()=='GT')&(str(s.seq[pos+19:pos+22]).upper()=='GGA'):
            strand = 0
            true_pos = pos+16
        elif (str(s.seq[pos-28:pos-25]).upper()=='TCC')&(str(s.seq[pos-24:pos-22]).upper()=='AC'):
            strand = 1
            true_pos = pos-28
        elif (str(s.seq[pos+16:pos+18]).upper()=='NN')&(str(s.seq[pos+19:pos+22]).upper()=='NNN'):
            strand = -1
            true_pos = -1
        elif (str(s.seq[pos-28:pos-25]).upper()=='NNN')&(str(s.seq[pos-24:pos-22]).upper()=='NN'):
            strand = -1
            true_pos = -1
        elif (len(str(s.seq[pos+16:pos+18]))<2)|(len(str(s.seq[pos+19:pos+22]))<3):
            strand = -1
            true_pos = -1
        else:
            raise Exception('bad case! {} {}, sequences: {} {}'.format(i, pos, s.seq[pos+16:pos+22], s.seq[pos-28:pos-22]))
        strands[i].append(strand)
        true_poss[i].append(true_pos)
    strands[i] = np.array(strands[i])
    true_poss[i] = np.array(true_poss[i])
    
genome_db.setEnzyme('NlaIII')

cmd_end_time = time.time()
cmd_runtime = timedelta(seconds=cmd_end_time - cmd_bgn_time)
logging.info( "Genome preparation finished, timing: {}".format(cmd_runtime) )

# Writing output for the genome
logging.info("Writing output hdf5 file...")
f = h5py.File( os.path.join(PATH_ABSOLUTE, "data/genome/rsites.hdf5"), "w")
n_chrms = len(genome_db.chrmLabels)
dset_nla = f.create_group("nla", (n_chrms,) )
dset_dpn = f.create_group("dpn", (n_chrms,) )
dset_mme_p = f.create_group("mme+", (n_chrms,) )
dset_mme_m = f.create_group("mme-", (n_chrms,) )

for ch in genome_db.label2idx.keys():
    rs = genome_db.rsites[genome_db.label2idx[ch]]-6
    rs = np.append(rs, len(genome_db.seqs[genome_db.label2idx[ch]]))
    dset = dset_nla.create_dataset(ch, data=rs)
    rs = true_poss[genome_db.label2idx[ch]][strands[genome_db.label2idx[ch]]==1]
    rs = np.append(rs, len(genome_db.seqs[genome_db.label2idx[ch]]))
    dset = dset_mme_p.create_dataset(ch, data=rs)
    rs = true_poss[genome_db.label2idx[ch]][strands[genome_db.label2idx[ch]]==0]
    rs = np.append(rs, len(genome_db.seqs[genome_db.label2idx[ch]]))
    dset = dset_mme_m.create_dataset(ch, data=rs)
    
genome_db.setEnzyme('DpnII')
for ch in genome_db.label2idx.keys():
    rs = genome_db.rsites[genome_db.label2idx[ch]]-6
    rs = np.append(rs, len(genome_db.seqs[genome_db.label2idx[ch]]))
    dset = dset_dpn.create_dataset(ch, data=rs)
    
f.close()
logging.info("Genome restriction sites preparation done!")
