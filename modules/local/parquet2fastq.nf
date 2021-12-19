// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '0.0'

process PARQUET2FASTQ {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "${moduleDir}/environment_pyarrow.yml" : null)

    input:
    tuple val(meta), path(parquets)

    output:
    tuple val(meta), path(parquets), path("*.fq.gz"), emit: output
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def parquets_concatenated = parquets.join(' ')
    def fragment_name = meta.fragment
    def side = meta.side
    def selection_criteria = meta.selection_criteria
    """
    convert_parquet2fastq.py ${fragment_name} readID R${side} Q${side} "${selection_criteria}" ${prefix}.${fragment_name}.fq ${parquets_concatenated}

    gzip ${prefix}.${fragment_name}.fq

    echo $VERSION > ${software}.version.txt
    """
}
