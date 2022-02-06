// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '0.0'

process RKLIB_QUERYSEARCH {
    tag "$meta.id $meta_oligos.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }
//    cache "${params.cache}"

    input:
    tuple val(meta), path(bin_reads)
    tuple val(meta_oligos), path(bin_oligos)

    output:
    tuple val(meta), path("*.tsv"), emit: aligned
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}.${meta_oligos.id}${options.suffix}" : "${meta.id}.${meta_oligos.id}"
//    meta.put('oligo', meta_oligos.id)
//    meta.put('side', meta_oligos.side)

    def left_allowed_shift = meta_oligos.left_allowed_shift
    def right_allowed_shift = (meta_oligos.right_allowed_shift.toInteger() > 0) ? meta_oligos.right_allowed_shift : meta.rlen.toInteger()+meta_oligos.right_allowed_shift.toInteger()

    def read_length = meta.rlen
    def oligo_length = meta_oligos.expected_length

    def right_to_left = meta_oligos['right_to_left'] ? 1: 0

    if (meta.single_end) {
        """
        rk_querysearch ${bin_oligos} ${bin_reads[0]} ${oligo_length} ${read_length} \
            ${meta_oligos.n_oligos} ${left_allowed_shift} ${right_allowed_shift} \\
            ${meta_oligos.mismatches_allowed} ${right_to_left} \\
            > ${prefix}.tsv

        echo $VERSION > ${software}.version.txt
        """
    } else {
        if (meta_oligos.side==1) {
            """
            rk_querysearch ${bin_oligos} ${bin_reads[0]} ${oligo_length} ${read_length} \\
                ${meta_oligos.n_oligos} ${left_allowed_shift} ${right_allowed_shift} \\
                ${meta_oligos.mismatches_allowed} ${right_to_left} \\
                > ${prefix}.R1.tsv

            echo $VERSION > ${software}.version.txt
            """
        } else if (meta_oligos.side==2) {
            """
            rk_querysearch ${bin_oligos} ${bin_reads[1]} ${oligo_length} ${read_length} \\
                ${meta_oligos.n_oligos} ${left_allowed_shift} ${right_allowed_shift} \\
                ${meta_oligos.mismatches_allowed} ${right_to_left} \\
                > ${prefix}.R2.tsv

            echo $VERSION > ${software}.version.txt
            """
        } else { // Apply to both sides
            """
            rk_querysearch ${bin_oligos} ${bin_reads[0]} ${oligo_length} ${read_length} \\
                ${meta_oligos.n_oligos} ${left_allowed_shift} ${right_allowed_shift} \\
                ${meta_oligos.mismatches_allowed} ${right_to_left} \\
                > ${prefix}.R1.tsv

            rk_querysearch ${bin_oligos} ${bin_reads[1]} ${oligo_length} ${read_length} \\
                ${meta_oligos.n_oligos} ${left_allowed_shift} ${right_allowed_shift} \\
                ${meta_oligos.mismatches_allowed} ${right_to_left} \\
                > ${prefix}.R1.tsv

            echo $VERSION > ${software}.version.txt
            """
        }
    }
}
