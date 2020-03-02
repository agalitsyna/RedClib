#!/bin/bash

mkdir -p logs

#python 00_download_GEO.py &>> logs/log_download.txt
#bash 01_compile.sh &>> logs/log_compilation.txt

for PREF in K562_rep1_add K562_rep1_wo-ligase
do
  echo $PREF
	gzip -d ./data/fastq/${PREF}_*.fastq.gz
	python 02_run_analysis.py $PREF &>> logs/log_${PREF}.txt &
	python 03_collect_hdf5.py $PREF &>> logs/log_${PREF}.txt
	python 04_filtering.py    $PREF &>> logs/log_${PREF}.txt
	gzip ../data/fastq/${PREF}_*.fastq.gz
done
wait
