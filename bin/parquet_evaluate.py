#!/usr/bin/env python3
# TODO: add check for the presence of the oclumn in the set or input arrays
from sys import argv

if len(argv) < 4:
    print(
        """Add columns to parquet file. Column schema should be tab-separated and contain three columns: columns_name, column_format and column_expression.
    Column expression is one-liner that does not contain lambda expressions, list comprehensions, and can use only input_parquet column names as variables, built-in functions and numpy for their evaluation.
    The result of evaluation should be a vector of type column_format with the number of entries equal to the input size of array columns.
    Usage: parquet_evaluate.py input.pq columns_schema output.pq"""
    )
    exit()

import pyarrow as pa
import pyarrow.parquet as pq
import numpy as np
import ast

input_parquet = argv[1]
columns_schema = argv[2]
output_parquet = argv[3]

prohibited_symbols = [":", ".", "-", "/", "!", "?", "&", "|", "'", "%", "@"]
pq_table = pa.parquet.read_table(input_parquet, memory_map=True)

arrays_dict = {}
n_evaluated = 0
with open(columns_schema, "r") as input_file:
    for line in input_file.readlines():
        column_name, column_format, column_expression = line.split("\t")

        assert np.all(
            [x not in column_name for x in prohibited_symbols]
        ), "Check the column name. It cannot contain " + ",".join(prohibited_symbols)

        # Loading dataset keys as variables:
        loaded_ids = []
        syntax_tree = ast.parse(column_expression)
        for node in ast.walk(syntax_tree):
            if type(node) is ast.Name:
                # The element is not already loaded and not a builtin name:
                if (node.id not in list(vars().keys())) and (
                    node.id not in dir(__builtins__)
                ):
                    try:
                        vars()[node.id] = pq_table[node.id].to_numpy()
                        loaded_ids.append(str(node.id))
                    except KeyError as e:
                        try:
                            vars()[node.id] = arrays_dict[node.id].to_numpy()
                            loaded_ids.append(str(node.id))
                        except KeyError as e:
                            raise ValueError(
                                f"Variable {node.id} is not available from input/created pyarrow file. "
                                f"List of variables that can be loaded:\n{ ', '.join( list(pq_table.column_names) ) }"
                            )
                    except Exception as e:
                        raise ValueError(e)

        # Evaluate expression:
        print(column_expression, column_name)
        result = eval(column_expression)

        if column_format.lower() == "str":
            pyarrow_format = pa.string()
        elif column_format.lower() == "int":
            pyarrow_format = pa.int32()
        elif column_format.lower() == "int8":
            pyarrow_format = pa.int8()
        elif column_format.lower() == "int16":
            pyarrow_format = pa.int16()
        elif column_format.lower() == "int32":
            pyarrow_format = pa.int32()
        elif column_format.lower() == "bool":
            pyarrow_format = pa.bool_()
        else:
            raise ValueError("Supported formats: str, int and bool for now")

        arrays_dict[column_name] = pa.array(result, type=pyarrow_format)
        n_evaluated += 1
        print(column_name, column_format, column_expression, arrays_dict.keys())


if n_evaluated > 0:
    print(
        f"Evaluated {n_evaluated} expressions, including columns: {', '.join(arrays_dict.keys())}"
    )
    pq_table_output = pa.Table.from_pydict(arrays_dict)
    # Write the output to the same input file:
    parquet_schema = pq_table_output.schema
    parquet_writer = pq.ParquetWriter(
        output_parquet, parquet_schema, compression="snappy"
    )
    parquet_writer.write_table(pq_table_output)
    parquet_writer.close()
else:
    print("No evaluated expression. Empty input table?")
