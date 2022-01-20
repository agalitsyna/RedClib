#!/usr/bin/env python3

import numpy as np
import pandas as pd
import click

# Read the arguments:
@click.command()
@click.argument("filename_bed", metavar="BED", type=click.Path(exists=True))
@click.argument("filename_rsites", metavar="RSITES", type=click.Path(exists=True))
@click.option(
    "-o",
    "--out",
    help="Save output restriction sites positions as TSV/PARQUET/HDF5.",
    required=True,
)
@click.option(
    "-s",
    "--strand",
    help="Strand to search for restriction sites: '+' plus strand, '-' minus strand or 'b' - both strands. ",
    default="b",
    show_default=True,
    type=click.Choice(["+", "-", "b"], case_sensitive=False),
)
@click.option(
    "-a",
    "--align-ids",
    help="Column to align IDs to this file. FORMAT::FILENAME::COLNAME. Optional. ",
    default=None,
)
@click.option(
    "-c",
    "--columns",
    help="Name of columns. Should be in the default order (start then end; left-side then right). Optional. ",
    default="id,start_left,start_right,end_left,end_right",
    type=str,
    show_default=True,
)
@click.option(
    "--out-format",
    help="Type of output. Optional.",
    default="TSV",
    type=click.Choice(["TSV", "PARQUET", "HDF5"], case_sensitive=False),
    show_default=True,
)
def get_closest_sites(
    filename_bed, filename_rsites, strand, align_ids, columns, out, out_format
):
    """
    Get distances to the closest sites to the start and end of mapped reads.
    Output file format: tsv file with fields:
    read id,
    distance to the closest site (strand-specific) to the read start from the left,
    distance to the closest site to the read start from the right,
    distance to the closest site to the read end from the left,
    distance to the closest site to the read end from the right.
    If the read is located before the first site or after the last site in the chromosome,
    then the reported distance will be artificially large (>> 1e9).
    """

    # Read the restriction sites:
    try:
        rsites = pd.read_csv(filename_rsites, header=None, sep="\t")
    except pd.errors.EmptyDataError:
        rsites = pd.DataFrame(columns=np.arange(6))

    rsites.columns = ["chrom", "start", "end", "name", "_", "strand"]
    if strand != "b":
        rsites = rsites.loc[rsites.strand == strand, :].sort_values(["chrom", "start"])
    else:
        rsites = rsites.sort_values(["chrom", "start"])
    rsites_grouped = rsites.groupby("chrom")

    try:
        bed_file = pd.read_csv(filename_bed, sep="\t", header=None)
    except pd.errors.EmptyDataError:
        bed_file = pd.DataFrame(columns=np.arange(7))
    bed_file.columns = ["chrom", "start", "end", "read_id", "q", "strand", "cigar"]

    # Convert dataframe to more effective numpy arrays:
    chs = bed_file.chrom.values.astype(str)
    ids = bed_file.read_id.values.astype(str)
    bgn = bed_file.start.values.astype(int)
    end = bed_file.end.values.astype(int)
    l = len(chs)

    # Numpy checks, more effective than pandas:
    dct = {}
    for k in ["start_left", "start_right", "end_left", "end_right"]:
        dct[k] = np.full(l, -1).astype(int)

    for ch in rsites_grouped.groups.keys():
        rs = np.concatenate(
            [[-1e10], rsites_grouped.get_group(ch)["start"].values, [1e10]]
        )
        mask = chs == ch

        idx = np.digitize(bgn[mask], rs)
        bgns = rs[idx - 1]
        ends = rs[idx]
        dct["start_left"][mask] = bgns - bgn[mask]
        dct["start_right"][mask] = ends - bgn[mask]

        idx = np.digitize(end[mask], rs)
        bgns = rs[idx - 1]
        ends = rs[idx]
        dct["end_left"][mask] = bgns - end[mask]
        dct["end_right"][mask] = ends - end[mask]

    if align_ids is not None:
        align_format, align_file, align_col = align_ids.split("::")
        if align_format.upper() == "PARQUET":
            import pyarrow as pa
            import pyarrow.parquet as pq

            pq_table = pq.read_table(align_file, memory_map=True)
            # print(pq_table.column_names)
            aln_ids = pq_table[align_col].to_numpy()
        elif align_format.upper() == "HDF5":
            import h5py

            h = h5py.File(align_file, "r")
            aln_ids = h[align_col].values()
        elif align_format.upper() == "TSV":
            df_ref = pd.read_csv(align_file, sep="\t")
            aln_ids = df_ref[align_col].values.astype(str)
        elif align_format.upper() == "CSV":
            df_ref = pd.read_csv(align_file)
            aln_ids = df_ref[align_col].values.astype(str)
        else:
            raise ValueError(
                f"Format {align_format} is not supported for alignment of ids. Available formats: PARQUET, TSV, CSV, HDF5."
            )

        df = pd.DataFrame({"id": aln_ids}).reset_index().set_index("id")
        ids_order = df.loc[ids, "index"]  # TODO: check and replace
        aln_l = len(aln_ids)
        assert (
            len(ids_order) == l
        ), f"Detected repeated ids, check the inputs. Ordered ids: {len(ids_order)}, input ids: {l}"
        dct_updated = {}
        for k in ["start_left", "start_right", "end_left", "end_right"]:
            dct_updated[k] = np.full(aln_l, -1).astype(int)
            dct_updated[k][ids_order] = dct[k]

        # Replace with the new values:
        ids = aln_ids
        dct = dct_updated
        l = len(ids)

    if out_format.upper() == "TSV":
        SEP = "\t"
        with open(out, "w") as outf:
            outf.write(columns.replace(",", SEP) + "\n")

            for i in range(l):
                outf.write(
                    ids[i]
                    + SEP
                    + str(dct["start_left"][i])
                    + SEP
                    + str(dct["start_right"][i])
                    + SEP
                    + str(dct["end_left"][i])
                    + SEP
                    + str(dct["end_right"][i])
                    + "\n"
                )

    elif out_format.upper() == "PARQUET":

        import pyarrow as pa
        import pyarrow.parquet as pq

        columns = columns.split(",")
        df = pd.DataFrame(
            {
                columns[0]: ids,
                columns[1]: dct["start_left"],
                columns[2]: dct["start_right"],
                columns[3]: dct["end_left"],
                columns[4]: dct["end_right"],
            }
        )

        frame = pa.Table.from_pandas(df=df)
        parquet_schema = frame.schema
        parquet_writer = pq.ParquetWriter(out, parquet_schema, compression="snappy")
        table = pa.Table.from_pandas(df, schema=parquet_schema)
        parquet_writer.write_table(table)

    elif out_format.upper() == "HDF5":
        import h5py

        h = h5py.File(out, "w")
        h.create_dataset(columns[0], data=ids)
        h.create_dataset(columns[1], data=dct["start_left"])
        h.create_dataset(columns[2], data=dct["start_right"])
        h.create_dataset(columns[3], data=dct["end_left"])
        h.create_dataset(columns[4], data=dct["end_right"])
        h.close()
        
    else:
        raise ValueError(f"Not implemented for format: {output_format}")


if __name__ == "__main__":
    get_closest_sites()
