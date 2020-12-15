#!/usr/bin/env python3

import numpy as np
import h5py

def update_hdf5(file_hdf5, dct, dtype=int):
    """dct might be ceither dictionary or a single data type"""

    ks = dct.keys()
    if not isinstance(dtype, dict):
        dtype = {k: dtype for k in ks}

    for k in ks:
        file_hdf5.create_dataset(k, data=np.array(dct[k], dtype=dtype[k]), chunks=True)


def reconstruct_by_ids(raw_dct, key, reference_list, default_dct):
    """
    Filling missing values in the dict read by raw_read_file.
    :param raw_dct: initial dict
    :param key: key with indexes in the raw_dct
    :param reference_list: reference list of indexes containing all values
    :param default_dct: dict with default values for columns from raw_dct that will be filled for missing indexes
    :return:
    """
    ret_dct = {k: [] for k in raw_dct.keys()}
    ref_dct = {k: i for i, k in enumerate(raw_dct[key])}

    keys = [x for x in raw_dct.keys() if x != key]

    for i in reference_list:
        try:
            idx = ref_dct[i]
            for k in keys:
                ret_dct[k].append(raw_dct[k][idx])
        except Exception as e:
            for k in keys:
                ret_dct[k].append(default_dct[k])
        ret_dct[key].append(i)

    return ret_dct

def raw_read_file(filename, columns, colnames, sep="\t", header=1, nrows=-1, modifiers=None, comment="#", bottom=0):
    """Read input text file with several columns line by line and saving as dictionary.

    :param filename: input file
    :param columns: numbers of columns to be parsed from file (ordered, 1-indexed as in awk)
    :param colnames: names of columns for the order from "columns" parameter
    :param sep: field separator in the input file
    :param header: number of rows to skip as header. Default: 1
    :param nrows: number of rows to read from file (after header). Default: -1 - read all lines
    :param modifiers: the function(s) to be applied to the columns from "columns" parameter. Can be None, one function or list of functions
    :param comment: skip the lines starting with "comment" character
    :param bottom: number of rows to skip from bottom (note that they should be parsed first)
    :return: dict with colnames as keys and list of columns in values.
    """
    assert len(columns) == len(colnames)
    assert header >= 0
    ret = {k: [] for k in colnames}
    with open(filename, "r") as inf:
        for i in range(header + 1):
            l = inf.readline().split(sep)
        nrows_passed = 0
        while len(l) > 1:
            if nrows > 0 and nrows_passed > nrows:
                break

            if len(l[0]) > 1 and l[0][0] == comment:
                l = inf.readline().split(sep)
                continue
            for i_mod, (i, col) in enumerate(zip(columns, colnames)):
                i1 = i - 1 if i > 0 else i
                val = l[i1].strip()
                if not modifiers is None:
                    val = modifiers[i_mod](val)
                ret[col].append(val)
            l = inf.readline().split(sep)
            nrows_passed += 1
    if bottom > 0:
        ret = {k: ret[k][:-bottom] for k in colnames}
    return (ret)