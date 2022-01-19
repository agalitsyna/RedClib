#!/usr/bin/env python3
from sys import argv

if len(argv) < 8:
    print(
        """Convert parquet to fastq file
    The result of evaluation should be a vector of type column_format with the number of entries equal to the input size of array columns.
    Usage: parquet2fastq.py fragment_name key_readid key_seq key_qual ${selection_criteria} output_filename pq1 pq2 ...}"""
    )
    exit()

import pyarrow as pa
import pyarrow.parquet as pq
import numpy as np
import ast

# Read the input
fragment_name = argv[1]
key_readID = argv[2]
key_sequence = argv[3]
key_quality = argv[4]
key_start, key_end = [f"{fragment_name}_start", f"{fragment_name}_end"]

selection_criteria = argv[5]
output_file = argv[6]
input_parquets = argv[7:]

for k in [key_start, key_end, key_readID, key_sequence, key_quality]:
    for input_parquet in input_parquets:
        pq_table = pa.parquet.read_table(input_parquet, memory_map=True)
        if k in pq_table.column_names:
            vars()[k] = pq_table[k].to_numpy()
            break
    else:
        raise Exception(
            f"Column {k} not fount in parquets: {', '.join(input_parquets)}"
        )

# Loading dataset keys as variables:
loaded_ids = []
syntax_tree = ast.parse(selection_criteria)
for node in ast.walk(syntax_tree):
    if type(node) is ast.Name:
        # The element is not already loaded and not a builtin name:
        if (node.id not in list(vars().keys())) and (node.id not in dir(__builtins__)):
            for input_parquet in input_parquets:
                pq_table = pa.parquet.read_table(input_parquet, memory_map=True)
                if node.id in pq_table.column_names:
                    vars()[node.id] = pq_table[node.id].to_numpy()
                    loaded_ids.append(str(node.id))
                    break
            else:
                raise Exception(
                    f"Column {k} not fount in parquets: {', '.join(input_parquets)}"
                )

mask = eval(selection_criteria)
selected = np.where(mask)[0]

with open(output_file, "w") as outf:
    readIDs = vars()[key_readID]
    seqs = vars()[key_sequence]
    quals = vars()[key_quality]
    starts = vars()[key_start]
    ends = vars()[key_end]
    for i in selected:
        outf.write(readIDs[i] + "\n")  # Sequence name
        outf.write(seqs[i][starts[i] : ends[i]] + "\n")  # Sequence
        outf.write("+\n")
        outf.write(quals[i][starts[i] : ends[i]] + "\n")  # Qualities

print(f"Done writing {len(selected)} sequences into {output_file} !")
