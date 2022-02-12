Red-C nextflow pipeline (DSL2 version)
=====================================

Introduction
------------

Welcome to Red-C RNA-DNA interactome mapping with `nextflow <https://www.nextflow.io/>`_, DSL2 version.

Red-C is a method for genome-wide detection of RNA-DNA interactions in chromatin.
The details are in our `NAR paper
<https://doi.org/10.1093/nar/gkaa457/>`_ from 09 July, 2020:

    Gavrilov, A.A., Zharikova, A.A., Galitsyna, A.A., Luzhin, A.V., Rubanova, N.M., Golov, A.K.,
    Petrova, N.V., Logacheva, M.D., Kantidze, O.L., Ulianov, S.V., Magnitov, M.D., Mironov, A.A., and Razin S.V. 2020.
    Studying RNA–DNA interactome by Red-C identifies noncoding RNAs associated with various chromatin
    types and reveals transcription dynamics.
    Nucleic Acids Research.

The initial version of the pipeline (RedClib) was implemented in Python, and it was hard to scale.
Next we switched to *nextflow* version of it, which was still lacking modularity and flexibility. 
Here we develop a DSL2 version of the pipeline that will be easy to adapt to custom user scenarios of RNA-DNA interactions. 

Disclaimer: DSL2 version is not complete yet!

System requirements and test run
-------------------

1. You need a system supporting nextflow pipelines. For example, you can create conda environment with nextflow, 
   which will contain redc-nf dependencies: ::


    conda env create -n redc-nf
    conda activate redc-nf
    conda install -c bioconda "nextflow>=21.10"
    bash ./bin/prepare_binaries.sh


RNA-DNA pipelines are powered by nextflow DSL2, thus you need a version that supports DSL2.
It is safe to use nextflow >= 21.10 version. Check the version of your nextflow *and update if it is below 21.10.*: ::

    nextflow -v

2. Run test example: ::

    nextflow run redc.nf -profile test,conda -params-file params-redc.yml

This example should take up to 30 minutes (depends a lot on your conda performance, as there are around five 
different environments constructed along the run).
You may want to check that the results in the `output/` folder. 

In the future, we plan to add a small test that will allow you to validate the installation and run on a test dataset. 

Disable nextflow's conda
------------------------

Sometimes nextflow's conda management is unsatisfactory (e.g. installation of environments takes too long for each run).
In that case, you might install all the requirements locally and disable conda in the config.
For that:

1. Make sure you have installed all the dependencies from the `environment.yml` file (main repository).

2. Change the line in the params file from `true` to `false`: ::

    enable_conda: false



Full run of the pipeline
------------------------

To use RNA-DNA pipelines to your own data, consider modifying or creating your own samplesheet (`samplesheet.csv`),
config file (`conf/test.config`) and parameters file (`params-redc.yml`).


Resume and other nextflow hints
-------------------------------

1. Resume after interrupted run:::

    nextflow run redc.nf -profile test,conda -params-file params-redc.yml -resume -cache true

2. Resume after parameters change looks the same. We do not recommend to resume if only the values of parameters were changed in the config.

3. Don't remove `work/` folder and `.nextflow*` files if you plan to run your pipeline multiple times on the same dataset in future. They serve as caching history.

4. When you are satisfied with the result and do not need to store cache, or if you want a fresh new start for your run, use: ::

    nextflow clean -f

5. Execution tracing. Nextflow allows [tracing and visualization of the execution](https://www.nextflow.io/docs/latest/tracing.html)
of your pipeline: ::

    nextflow run redc.nf -profile test,conda,debug -with-dag flowchart.png -with-report report.html -with-timeline timeline.html

