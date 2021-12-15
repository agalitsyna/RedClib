// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '0.0'

process OLIGOS_CHECK {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(table)
    tuple val(meta_oligos), path(aligned)

    output:
    tuple val(meta), path("*.tsv"), emit: hits
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    check_oligo_presence.py ${table} ${aligned} ${options.args.oligo} ${options.args.orientation} ${options.args.position} ${prefix}.${meta_oligos.oligo}.${options.args.oligo}.tsv

    echo $VERSION > ${software}.version.txt
    """

}
