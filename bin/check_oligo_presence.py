#!/usr/bin/env python3

"""
Script that takes mapped oligos positions and input sequencing table and checks that certain positions in the 
reads have the requested sequence. This is needed when initial search of oligos allowed mistakes, but 
you want to make sure some positions are retained with no mismatches. 

This is technical script used in original RedC paper. Bridge adaptor there required to have GA at its end. 

Usage:
check_oligo_presence.py file_table_fastq file_oligo_hits(output of rk) oligo(sequence) orientation(F,R) position_in_oligo(int) file_output
"""

from sys import argv

if len(argv)!= 7:
    print("""Script that takes mapped oligos positions and input sequencing table and checks that certain positions in the 
reads have the requested sequence. This is needed when initial search of oligos allowed mistakes, but 
you want to make sure some positions are retained with no mismatches. 

This is technical script used in original RedC paper. Bridge adaptor there required to have GA at its end. 

Usage:
check_oligo_presence.py file_table_fastq file_oligo_hits(output of rk) oligo(sequence) orientation(F,R) position_in_oligo(int) file_output

Example usage: 
python bin/check_oligo_presence.py output/table/test-download_01.fastq.txt output/oligos/test-download_01.bridge_forward.R1.tsv GA F 35 tmp.txt
""")
    exit()

file_fastq   = argv[1]
file_oligo_hits  = argv[2]
oligo = argv[3]
orientation = argv[4]
position_in_oligo = int(argv[5])
file_output  = argv[6]

with open(file_output, 'w') as outf:
    outf.write(f"#entry_index\toligo_{oligo}_present_at_{position_in_oligo}\n")
    with open(file_fastq, 'r') as in_f:
        with open(file_oligo_hits, 'r') as hits_f:
            fastq_table_line = in_f.readline()
            hits_table_line = hits_f.readline()

            if fastq_table_line.startswith('#'):
                fastq_table_line = in_f.readline()
            if hits_table_line.startswith('#'):
                hits_table_line = hits_f.readline()

            while len(fastq_table_line) > 0:

                if orientation=='F':
                    read = fastq_table_line.split()[2] # Read forward read
                elif orientation=='R':
                    read = fastq_table_line.split()[4] # Read reverse read
                else:
                    raise ValueError(f"Orientation '{orientation}' is not supported. Please, specify 'F' or 'R' instead")

                oligo_position_start = int(hits_table_line.split()[3]) # Start position of the oligo in the read
                idx = hits_table_line.split()[0] # Number of the read in the table with hits
                if oligo_position_start > len(read):
                    ret = 0
                else:
                    if read[oligo_position_start+position_in_oligo : oligo_position_start+position_in_oligo+len(oligo)] == oligo:
                        ret = 1
                    else:
                        ret = 0
                outf.write("{}\t{}\n".format(idx, ret))

                # Continue to the next read:
                hits_table_line = hits_f.readline()
                fastq_table_line = in_f.readline()
