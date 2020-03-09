#!/bin/bash

mkdir -p logs

#bash 01_compile.sh &>> logs/log_compilation.txt

PREF="sample"
#gzip -d ./data/fastq/${PREF}_*
python 02_run_analysis.py $PREF &>> logs/log_${PREF}.txt
python 03_collect_hdf5.py $PREF &>> logs/log_${PREF}.txt
python 04_filtering.py    $PREF &>> logs/log_${PREF}.txt
