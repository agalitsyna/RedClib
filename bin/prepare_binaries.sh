#!/bin/bash

gcc -O3 -o bin/align_universal lib/align_universal.c lib/compare.c
gcc -O3 -o bin/align_pairwise  lib/align_pairwise.c  lib/compare.c
gcc -O3 -o bin/fasta2bin lib/fasta2bin.c lib/compare.c
gcc -O3 -o bin/fastq2bin lib/fastq2bin.c lib/compare.c

wget https://raw.githubusercontent.com/DaehwanKimLab/hisat2/master/hisat2_extract_splice_sites.py -O ./bin/hisat2_extract_splice_sites.py

echo "Binaries compiled!"