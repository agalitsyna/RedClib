#!/usr/bin/env python3

from sys import argv
if len(argv)!=4 and len(argv)!=5:
    print("Usage: tsv2parquet.py csv_file parquet_file chunksize (suffix, optional)")
    exit()

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

csv_file = argv[1]
parquet_file = argv[2]
chunksize = int(argv[3])
if len(argv)==5:
    suffix = argv[4]
else:
    suffix = ''

csv_stream = pd.read_csv(csv_file, sep='\t', chunksize=chunksize, low_memory=True)

for i, chunk in enumerate(csv_stream):
    if i == 0:
        frame = pa.Table.from_pandas(df=chunk.rename(columns={x: x+suffix for x in chunk.columns}))
        parquet_schema = frame.schema
        parquet_writer = pq.ParquetWriter(parquet_file, parquet_schema, compression='snappy')
    table = pa.Table.from_pandas(chunk.rename(columns={x: x+suffix for x in chunk.columns}), schema=parquet_schema)
    parquet_writer.write_table(table)

parquet_writer.close()
