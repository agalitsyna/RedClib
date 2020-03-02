# Red-C processing library

## Introduction

Welcome to Red-C RNA-DNA interactome mapping library!

Red-C is a novel method for genome-wide detection of RNA-DNA interactions in chromatin.
The manuscript is currently under preparation and after its release you can reproduce the Red-C results with RedClib.

## System requirements

You will need (updated upon release) of free space to download raw datasets from GEO, including Red-C and RNA-Seq datasets.
Keep in mind that in our pipeline we store intermediary processing results, so that you might need up to (updated upon release) of memory
to reproduce the pipeline for all datasets.

The code is adapted for Linux systems and bash shell. We recommend using conda for requirements management.

## Environment setup

1. After conda installation and activation run:

```bash
conda create --name RedC python=3.7 numpy h5py pandas
conda activate RedC
conda install -c anaconda git pip cython 
conda install -c bioconda hisat2 fastuniq pyfaidx
pip install https://bitbucket.org/mirnylab/mirnylib/get/tip.tar.gz
conda install -c gusdunn geoparse    # For data download with GEOparse; skip if not needed
conda install -c bioconda sra-tools  # For data download with GEOparse; skip if not needed
```

This will create RedC conda environment with all the required packages for data download and processing.

2. Download RedClib source code:

```bash
git clone git@github.com:agalitsyna/RedClib.git
cd RedClib
```

## Pipeline run

You can run processing pipeline as a whole:

```bash
bash processing.sh
```

You may run the steps of the script one-by-one.

Note that data download step not working until the public data release, thus make sure you have folder data/fastq with all fastq files in advance:

```bash
ls data/fastq
    fibr_rep1_add_R1.fastq.gz
    fibr_rep1_add_R2.fastq.gz
    fibr_rep1_R1.fastq.gz
    fibr_rep1_R2.fastq.gz
    fibr_rep2_R1.fastq.gz
    fibr_rep2_R2.fastq.gz
    HeLa_DRB_R1.fastq.gz
    HeLa_DRB_R2.fastq.gz
    HeLa_G1_R1.fastq.gz
    HeLa_G1_R2.fastq.gz
    HeLa_M_R1.fastq.gz
    HeLa_M_R2.fastq.gz
    K562_rep1_add_R1.fastq.gz
    K562_rep1_add_R2.fastq.gz
    K562_rep1_R1.fastq.gz
    K562_rep1_R2.fastq.gz
    K562_rep1_RNASeq.fastq.gz
    K562_rep1_rnase_R1.fastq.gz
    K562_rep1_rnase_R2.fastq.gz
    K562_rep1_wo-ligase_R1.fastq.gz
    K562_rep1_wo-ligase_R2.fastq.gz
    K562_rep2_R1.fastq.gz
    K562_rep2_R2.fastq.gz
    K562_rep2_RNASeq.fastq.gz
```

Note that some elements of code use file prefix name to determine the pre-defined read length in fastq files (e.g. 02_run_analysis.py), so make sure to use the mentioned files and names, or you need to modify some elements of scrits. 


For a test run you can use sample dataset provided with the library: 

```bash
ls data/fastq/sample*
    sample_R1.fastq.gz
    sample_R2.fastq.gz
```

```bash
bash processing_test.sh
```


## Processing steps description

### 00_download_GEO.py

You can download GEO Red-C datasets manually or use [GEOparse](https://github.com/guma44/GEOparse) for batch download.
In the latter case make sure you have installed GEOparse and sra-tools in your environment, pass GEO accession of Red-C
and any e-mail to download the datasets automatically.

```bash
python 00_download_GEO.py GSE... any@email.com
```

### 01_compile.sh

RedClib uses custom implementation of [Rabin-Karp algorithm](https://ieeexplore.ieee.org/document/5390135/) for DNA
hashing and matching adapter, bridge and reverse complement sequences. The implementation is in C, thus first you need
to compile the binaries with gcc. 01_compile.sh will do that and install compiled binaries to ./data/bin/ folder.

RedClib uses hisat2 for RNA and DNA parts mapping to hg19 reference genome. 01_compile.sh will download the
 reference genome and run hisat2-build to build hisat2 index (will be located in ./data/genome/hisat2 folder).
 
Trimmomatic can be installed with conda, but we prefer to download the compiled source from the [Trimmomatic website](http://www.usadellab.org), it will be installed to .data/bin/ folder as well. 

### 02_run_analysis.py

You can run reads processing for each dataset separately:

```bash
python 02_run_analysis.py K562_rep1_wo-ligase
```

This command takes FASTQ forward and reverse reads for Red-C and produces a set of intermediary text files.
It prints extended log files with information about success of each processing step.

Mapping is done in parallel (3 threads). Feel free to modify this parameter inside 02_run_analysis.py (nthreads).

### 03_collect_hdf5.py

```bash
python 03_collect_hdf5.py K562_rep1_wo-ligase
```

This command collects output from 02_run_analysis.py and produces a single [hdf5](https://www.hdfgroup.org/solutions/hdf5/)
file with relevant information for each read.

### 04_filtering.py

```bash
python 04_filtering.py K562_rep1_wo-ligase
```

This command applies filters for each read and produces filtering statistics and final TSV file with observed RNA-DNA contacts.

For more details, see Methods section of upcoming manuscript.

## Running analysis on your data

The RedClib code was optimized and tested for a given dataset (including FASTQ files with their specific lengths, oligos with adapters, and bridge).
Before running it for your custom files, please, note several important details.

1. The read length associated with the FASTQ file is pre-defined. It is explicitly hard-coded in 02_run_analysis.py for the given file names.

If you want to introduce your file name, you should add its prefix to the "lengths_dict" dictionary in 02_run_analysis.py. Note that all reads in FASTQ file should be the same length.

2. The read length is passed to the C-program "align_universal.c" that performs a search of oligos in the FASTQ sequences. It is essential for the search algorithm, that can now work with typical Illumina reads lengths: 151, 133, 125, 101, 75, 80. All other read lengths are not currently implemented.

3.d One can use the custom adapter or oligos sequences. Example files are stored as FASTA files in data/oligos/ folder. If you want to use custom FASTA files, note that files should have Unix newline (\n) at the end of each line, including the last one. FASTA file is a text file, and according to [unix specification](https://pubs.opengroup.org/onlinepubs/009695399/basedefs/xbd_chap03.html#tag_03_392), a text file is a single or a series of lines, each of which ends with a newline character (\n).

Check the presence of newline character at the end of your files with:
```bash
od -t c br_37_for.fasta
```
