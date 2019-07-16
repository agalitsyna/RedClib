from RedClib import *
from sys import argv
import numpy as np
import h5py

prefix = argv[1]
formatting_dct = {'prefix': prefix}

PATH_ABSOLUTE = os.path.dirname(os.path.realpath(__file__))
PATH_OUT = os.path.join(PATH_ABSOLUTE, "data/out/")
PATH_OUT1 = os.path.join(PATH_ABSOLUTE, "data/results/")

if not os.path.isdir(PATH_OUT1):
    os.mkdir(PATH_OUT1)

filename = os.path.join(PATH_OUT, "table_{prefix}.hdf5".format(**formatting_dct))

outfile = h5py.File(filename, "r")

l = len(outfile['id'][()])

# Final table masks

V0 = l

# Filter 1: Read is not PCR duplication
filt1 = outfile['is_notPCRdup'][()].astype(bool)

V1 = np.sum(filt1)

# Filter 2: there is a bridge in the read, with last AG letters, it is not cut by the quality filter
filt2 = (np.logical_not(outfile['has_nobridge'][()]) & outfile['has_GA'][()]).astype(bool)  # & (outfile['trimF'][()]>=outfile['bridge_end'][()])

V2 = np.sum(filt1&filt2)

# Filter 3: Reverse read starts with GGG
filt3 = (outfile['has_noggg'][()]==0).astype(bool)


dna_len_trim  = np.minimum( outfile['dna_R1_len_trim'][()], outfile['dna_R1_len_notrim'][()] )
rna_len_trim  = np.minimum( outfile['rna_R2_len_trim'][()], outfile['rna_R2_len_notrim'][()] )
rna1_len_trim = np.minimum( outfile['rna1_R1_len_trim'][()], outfile['rna1_R1_len_notrim'][()] )
dna_len  = outfile['dna_R1_len_notrim'][()]
rna_len  = outfile['rna_R2_len_notrim'][()]
rna1_len = outfile['rna1_R1_len_notrim'][()]

# Filter 4: DNA length 18-20
filt4 = ((dna_len>=18) & (dna_len<=20)).astype(bool)

# Filter 5: RNA is >= 14 bp before trimming and RNA1 >= 14 bp
filt5 = ( (rna_len>=14) & (rna1_len>=14) ).astype(bool)

V3 = np.sum(filt1&filt2&filt3)
V4 = np.sum(filt1&filt2&filt3&filt4)
V35 = np.sum(filt1&filt2&filt3&filt5)
V5 = np.sum(filt1&filt2&filt3&filt4&filt5)

V24 = np.sum(filt1&filt2&filt4)
V25 = np.sum(filt1&filt2&filt5)

# Filter 6: Trimming has passed
filt6 = ((outfile['trimF'][()]>0)&(outfile['trimR'][()]>0)).astype(bool)

# Filter 6A: Trimmed parts are > 14 bp
filt6a = ( (dna_len_trim>=18) & \
             (dna_len_trim<=20) & \
             (rna_len_trim>=14) & \
             (rna1_len_trim>=14) ).astype(bool)


# Filter 7: All three parts are uniquely mapped to the canonical chrms
filt7    = ((outfile[ 'dna_is_not_multi'][()]) & \
            (outfile[ 'rna_is_not_multi'][()]) & \
            (outfile['rna1_is_not_multi'][()]) & \
            (outfile[ 'dna_is_mapped'][()]) & \
            (outfile[ 'rna_is_mapped'][()]) & \
            (outfile['rna1_is_mapped'][()]) & \
            outfile[ 'dna_chr_canonical'][()] & \
            outfile[ 'rna_chr_canonical'][()] & \
            outfile['rna1_chr_canonical'][()] 
            ).astype(bool)

filt7a = ( (outfile['rna_chr'][()]==outfile['rna1_chr'][()]) ).astype(bool)
filt7b = ( (outfile['rna_strand'][()]!=outfile['rna1_strand'][()]) ).astype(bool)


# Filter 8: Same pos of DNA and RNA. Removal of NlaIII and MmeI contamination
filt8a = (outfile['RNADNASamePos'][()]==0).astype(bool)

filt8b = ((outfile['rna_nla_failed'][()]==0) & \
    (outfile['rna_mme_failed'][()]==0) & \
    (outfile['rna1_nla_failed'][()]==0) ).astype(bool)

# Filter 9: RNA1 and RNA distance is small enough
filt9 = ( np.abs(outfile['rna_bgn'][()]-outfile['rna1_bgn'][()])<1e4 ).astype(bool)


# Sequential filteres applied to initial set of reads:
V6  = np.sum(filt1&filt2&filt3&filt4&filt5&filt6)
V6a = np.sum(filt1&filt2&filt3&filt4&filt5&filt6a)

V7  = np.sum(filt1&filt2&filt3&filt4&filt5&filt6&filt6a&filt7)
V7a = np.sum(filt1&filt2&filt3&filt4&filt5&filt6&filt6a&filt7&filt7a)
V7b = np.sum(filt1&filt2&filt3&filt4&filt5&filt6&filt6a&filt7&filt7a&filt7b)

V8a = np.sum(filt1&filt2&filt3&filt4&filt5&filt6&filt6a&filt7&filt7a&filt7b&filt8a)
V8b = np.sum(filt1&filt2&filt3&filt4&filt5&filt6&filt6a&filt7&filt7a&filt7b&filt8a&filt8b)
V9 = np.sum(filt1&filt2&filt3&filt4&filt5&filt6&filt6a&filt7&filt7a&filt7b&filt8a&filt8b&filt9)


subsets = [V0, V1, V2, V3, V4, V5, V35, V24, V25, V6, V6a, V7, V7a, V7b, V8a, V8b, V9]

print('\t'.join([str(x) for x in [">", "new", prefix]+subsets]))


mask = filt1&filt2&filt3&filt4&filt5&filt6&filt6a&filt7&filt7a&filt7b&filt8a&filt8b&filt9
idx = np.where(mask)[0]

columns = """id rna1_chr rna1_bgn rna1_end rna1_strand rna1_cigar rna_chr rna_bgn rna_end rna_strand rna_cigar
dna_chr dna_bgn dna_end dna_strand dna_cigar""".split()

# Change the names from rna1-rna to 3'-5' terminology:
columns_names = """id 3rna_chr 3rna_bgn 3rna_end 3rna_strand 3rna_cigar 
5rna_chr 5rna_bgn 5rna_end 5rna_strand 5rna_cigar 
dna_chr dna_bgn dna_end dna_strand dna_cigar""".split()
dct = {k:outfile[k][()][idx] for k in columns}

output = os.path.join(PATH_ABSOLUTE, "data/results/{prefix}_passed.tsv".format(**formatting_dct))

def strand(x):
    if x:
        return "+"
    return "-"

with open(output, 'w') as outf:
    line = [""]+[k for k in columns_names]
    outf.write("\t".join(line)+'\n')
    for i in range(len(idx)):
        line = [ (str(dct[k][i]) if not ('chr' in k or 'id' in k or 'cigar' in k) else dct[k][i].decode()) if not 'strand' in k else strand(dct[k][i])  for k in columns ]
        outf.write(" ".join(line)+'\n')