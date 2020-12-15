#!/usr/bin/env python3

"""
Apply filters for the input hdf5 file.

Usage:
python filter_data.py test.hdf5 output.hdf5 chromsizes.file chroms_pattern filters.file

Example run:
python bin/filter_data.py tmp.hdf5 output.hdf5 \
 results/genome/hg19.chromsizes.txt "chr[0-9,X,Y,M]|chr[0-9][0-9]" filters.txt


"""

from sys import argv

import h5py
import numpy as np
import pandas as pd

import re
import ast

input_file = argv[1]
output_file = argv[2]

chrms_file = argv[3]
canonical_chrms_filt = argv[4]

filters_file = argv[5]

infile = h5py.File(input_file, "r")
outfile = h5py.File(output_file, "w")

outfile.create_dataset("id", data=np.array(infile['id'][()], dtype='S100'))

### Select canonical chromosomes:

canonical_chrms_filt = "("+canonical_chrms_filt+")$"
chrms = [str(x).encode() for x in
         pd.read_csv(chrms_file, sep='\t', header=None)[0].values[:-1]
         if re.match(canonical_chrms_filt, x, re.S) is not None]

v = np.in1d(infile['dna_chr'][()], chrms)
outfile.create_dataset("dna_chr_canonical", data=np.array(v, dtype=bool))

v = np.in1d(infile['rna2_chr'][()], chrms)
outfile.create_dataset("rna2_chr_canonical", data=np.array(v, dtype=bool))

v = np.in1d(infile['rna1_chr'][()], chrms)
outfile.create_dataset("rna1_chr_canonical", data=np.array(v, dtype=bool))


### Restriction and additional filters:
with open(filters_file, 'r') as filt_file:
    filters = filt_file.readlines()
    filters = [f.strip().split(':') for f in filters]
    for key, filter in filters:
        print(key, filter)
        # Loading dataset keys as variables:
        loaded_ids = []
        syntax_tree = ast.parse(filter)
        for node in ast.walk(syntax_tree):
            if type(node) is ast.Name:
                # The element is not already loaded and not a builtin name:
                if (node.id not in list(vars().keys())) and (node.id not in dir(__builtins__)):
                    try:
                        vars()[node.id] = infile[node.id][()]
                        loaded_ids.append(str(node.id))
                    except KeyError as e:
                        try:
                            vars()[node.id] = outfile[node.id][()]
                            loaded_ids.append(str(node.id))
                        except KeyError as e:
                            raise Exception(f"Variable {node.id} is not available from input/created hdf5 file. "
                                            f"List of variables that can be loaded:\n{', '.join(list(infile.keys())+list(outfile.keys()))}")
                        except Exception as e:
                            raise e
                    except Exception as e:
                        raise e
        # Evaluate expression:
        result = eval(filter)
        # Save resulting filter to hdf5:
        outfile.create_dataset(key, data=np.array(result, dtype=bool))
        # Clear workspace from heady loaded variables:
        for id in loaded_ids:
            del vars()[id]

outfile.close()
