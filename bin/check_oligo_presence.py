#!/usr/bin/env python3

"""
Usage:
check_oligo_presence.py input_fastq input_table output_file
"""

from sys import argv

file_fastq   = argv[1]
file_bridge  = argv[2]
file_output  = argv[3]
oligo = 'GA'

with open(file_output, 'w') as outf:
    with open(file_fastq, 'r') as in_f:
        with open(file_bridge, 'r') as br_f:
            br = br_f.readline()

            t = in_f.readline()
            br = br_f.readline()

            while len(t) > 0:
                readF = t.split()[2]
                end_pos = int(br.split()[4])
                idx = br.split()[0]
                if end_pos > len(readF):
                    ret = 0
                else:
                    if readF[end_pos - 2:end_pos] == oligo:
                        ret = 1
                    else:
                        ret = 0
                outf.write("{}\t{}\n".format(idx, ret))
                t = in_f.readline()
                br = br_f.readline()