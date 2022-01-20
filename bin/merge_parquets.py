#!/usr/bin/env python3
from sys import argv

if len(argv) < 3:
    print(
        'Merge parquet tables by columns. Suffixes will be added to all column in the dataset after "__". '
        'Usage: merge_parquets.py output.pq file1.pq file2.pq file3.pq file4:suffix ... (any number of files)'
    )
    exit()

import pyarrow as pa
import pyarrow.parquet as pq

output = argv[1]

columns = []
schema = []
for input in argv[2:]:
    filename = input.split(":")[0]
    suffix = ("__" + "".join(input.split(":")[1:])) if (":" in input) else ""
    pq_array = pq.read_table(filename, memory_map=True)
    pq_array = pq_array.rename_columns([x + suffix for x in pq_array.column_names])
    columns += pq_array.columns
    schema.append(pq_array.schema)

parquet_schema = pa.unify_schemas(schema)
pq_merged = pa.Table.from_arrays(columns, schema=parquet_schema)
parquet_writer = pq.ParquetWriter(output, parquet_schema, compression="snappy")
parquet_writer.write_table(pq_merged)
parquet_writer.close()
