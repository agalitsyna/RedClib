# Red-C nextflow pipeline (DSL2 version)

## Introduction

Welcome to Red-C RNA-DNA interactome mapping with [nextflow](https://www.nextflow.io/), DSL2 version.

Red-C is a method for genome-wide detection of RNA-DNA interactions in chromatin.
The details are in our [NAR paper](https://doi.org/10.1093/nar/gkaa457) from 09 July, 2020:

```
Gavrilov, A.A., Zharikova, A.A., Galitsyna, A.A., Luzhin, A.V., Rubanova, N.M., Golov, A.K., 
Petrova, N.V., Logacheva, M.D., Kantidze, O.L., Ulianov, S.V., Magnitov, M.D., Mironov, A.A., and Razin S.V. 2020. 
Studying RNA–DNA interactome by Red-C identifies noncoding RNAs associated with various chromatin 
types and reveals transcription dynamics. 
Nucleic Acids Research.
```

The initial version of the pipeline (RedClib) was implemented in Python, and it was hard to scale. 
Next we switched to *nextflow* version of it, which was still lacking modularity and flexibility. 
Here we develop a DSL2 version of the pipeline that will be easy to adapt to custom user scenarios of RNA-DNA interactions. 

Disclaimer: DSL2 version is not complete yet!

## System requirements

1. You need a system supporting nextflow pipelines. For example, you can create conda environment with nextflow, 
   which will contain redc-nf dependencies:

```
conda env create -n redc-nf
conda activate redc-nf
conda install -c bioconda nextflow
bash ./bin/prepare_binaries.sh
```
   
2. Run test example:

```
nextflow run main.nf -profile test,conda,debug
```

This example should take up to 4 minutes. 
you may want to check that the results are the same as we provide in `output/` folder. 
