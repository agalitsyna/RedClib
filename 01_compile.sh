#!/bin/bash
mkdir -p bin
gcc -O3 -o bin/align_universal lib/align_universal.c lib/compare.c
gcc -O3 -o bin/align_pairwise  lib/align_pairwise.c  lib/compare.c
gcc -O3 -o bin/fasta2bin lib/fasta2bin.c lib/compare.c
gcc -O3 -o bin/fastq2bin lib/fastq2bin.c lib/compare.c

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz -O ./data/genome/hg19.fa.gz
gzip -d ./data/genome/hg19.fa.gz
mkdir ./data/genome/hisat2/
hisat2-build ./data/genome/hg19.fa ./data/genome/hisat2/hg19