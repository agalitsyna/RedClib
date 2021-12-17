// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PARQUET_MERGE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "${moduleDir}/environment_pyarrow.yml" : null)
//    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
//        container "https://depot.galaxyproject.org/singularity/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0"
//    } else {
//        container "quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0"
//    }

    input:
    val(meta)
    path(pqs)
    val(suffixes)

    output:
    //tuple val(meta), path("*.pq"), emit: parquet
    path  "*.version.txt"         , emit: version

    tuple val(meta), path("*fastq.gz"), optional:true, emit: fastq

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def chunksize = options.get('chunksize', 1_000_000)
    def filename = table.take(table.lastIndexOf('.'))

    """
    #merge_parquets.py ${table} ${filename}.pq ${chunksize}

    python -c "import pyarrow; print(pyarrow.__version__)" > ${software}.version.txt
    """
}
