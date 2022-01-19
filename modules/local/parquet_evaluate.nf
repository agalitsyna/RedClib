// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '0.0'

process PARQUET_EVALUATE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "${moduleDir}/environment_pyarrow.yml" : null)

    input:
    tuple val(meta), path(parquet)
    val(filters)

    output:
    tuple val(meta), path(parquet), path("*.pq"), emit: parquet
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def format = options.args.format
    def cmd_create_evaluation_scheme = "printf '' > evaluation_schema.txt\n"
    for (filter in filters.keySet()) {
        cmd_create_evaluation_scheme += "printf \"${filter}\\t${format}\\t${filters[filter]}\\n\" >> evaluation_schema.txt\n"
    }
    """
    ${cmd_create_evaluation_scheme}
    parquet_evaluate.py ${parquet} evaluation_schema.txt ${prefix}.pq

    echo $VERSION > ${software}.version.txt
    """
}
