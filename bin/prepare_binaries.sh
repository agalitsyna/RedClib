#!/bin/bash

gcc -O3 -o bin/rk_querysearch lib/rk_querysearch.c lib/librk.c
gcc -O3 -o bin/rk_pairwise  lib/rk_pairwise.c  lib/librk.c
gcc -O3 -o bin/fasta2hash lib/fasta2hash.c lib/librk.c
gcc -O3 -o bin/fastq2hash lib/fastq2hash.c lib/librk.c

wget https://raw.githubusercontent.com/DaehwanKimLab/hisat2/master/hisat2_extract_splice_sites.py -O ./bin/hisat2_extract_splice_sites.py

echo "Binaries compiled!"