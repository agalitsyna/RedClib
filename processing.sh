#!/bin/bash

mkdir -p logs

#python 00_download_GEO.py
bash 01_compile.sh

for PREF in K562_rep1_wo-ligase K562_rep1_rnase #K562_rep1 K562_rep1_add K562_rep2 K562_rep1 K562_rep2 HeLa_M HeLa_DRB HeLa_G1 fibr_rep1_add fibr_rep1 fibr_rep2
do
	gzip -d ../data/fastq/${PREF}*
	python 02_run_analysis.py $PREF &>> logs/log_${PREF}.txt
	python 03_collect_hdf5.py $PREF &>> logs/log_${PREF}.txt
	python 04_filtering.py    $PREF &>> logs/log_${PREF}.txt"
done