#!/bin/bash

# Compiling C binaries
mkdir -p bin
gcc -O3 -o bin/align_universal lib/align_universal.c lib/compare.c
gcc -O3 -o bin/align_pairwise  lib/align_pairwise.c  lib/compare.c
gcc -O3 -o bin/fasta2bin lib/fasta2bin.c lib/compare.c
gcc -O3 -o bin/fastq2bin lib/fastq2bin.c lib/compare.c

# Downloading and compiling genome
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz -O ./data/genome/hg19.fa.gz
gzip -d ./data/genome/hg19.fa.gz
mkdir ./data/genome/hisat2/
hisat2-build ./data/genome/hg19.fa ./data/genome/hisat2/hg19

# Downloading trimmomatic
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip -O ./bin/Trimmomatic-0.39.zip
unzip bin/Trimmomatic-0.39.zip -d bin/

# Preparing restriction sites of the genome
python prepare_genome.py