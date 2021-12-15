// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '0.0'

process OLIGOS_ALIGN {
    tag "$meta_reads.id $meta_oligos.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
            meta_reads.oligo = meta_oligos.id
            meta_reads.side = meta_oligos.side
            saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta_reads, publish_by_meta:['id']) }

    input:
    tuple val(meta_oligos), path(bin_oligos), val(meta_reads), path(bin_reads)

    output:
    tuple val(meta_reads), path("*.tsv"), emit: aligned
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta_reads.id}.${meta_oligos.id}${options.suffix}" : "${meta_reads.id}.${meta_oligos.id}"
    meta_reads.oligo = meta_oligos.id
    meta_reads.side = meta_oligos.side

    if (!lengths.containsKey(meta_reads.rlen.toInteger())) {
        throw new Exception("Read length is not supported for this type of script")
    }

    def left_allowed_shift = meta_oligos.left_allowed_shift
    def right_allowed_shift = (meta_oligos.right_allowed_shift.toInteger() > 0) ? meta_oligos.right_allowed_shift : meta_reads.rlen.toInteger()+meta_oligos.right_allowed_shift.toInteger()

    if (meta_reads.single_end) {
        """
        rk_querysearch ${bin_oligos} ${bin_reads[0]} ${meta_reads.rlen} ${meta_oligos.expected_length} \\
            ${meta_oligos.n_oligos} ${left_allowed_shift} ${right_allowed_shift} \\
            ${meta_oligos.mismatches_allowed} \\
            > ${prefix}.tsv

        echo $VERSION > ${software}.version.txt
        """
    } else {
        if (meta_oligos.side==1) {
            """
            rk_querysearch ${bin_oligos} ${bin_reads[0]} ${meta_reads.rlen} ${meta_oligos.expected_length} \\
                ${meta_oligos.n_oligos} ${left_allowed_shift} ${right_allowed_shift} \\
                ${meta_oligos.mismatches_allowed} \\
                > ${prefix}.R1.tsv

            echo $VERSION > ${software}.version.txt
            """
        } else if (meta_oligos.side==2) {
            """
            rk_querysearch ${bin_oligos} ${bin_reads[1]} ${meta_reads.rlen} ${meta_oligos.expected_length} \\
                ${meta_oligos.n_oligos} ${left_allowed_shift} ${right_allowed_shift} \\
                ${meta_oligos.mismatches_allowed} \\
                > ${prefix}.R2.tsv

            echo $VERSION > ${software}.version.txt
            """
        } else { // Apply to both sides
            """
            rk_querysearch ${bin_oligos} ${bin_reads[0]} ${meta_reads.rlen} ${meta_oligos.expected_length} \\
                ${meta_oligos.n_oligos} ${left_allowed_shift} ${right_allowed_shift} \\
                ${meta_oligos.mismatches_allowed} \\
                > ${prefix}.R1.tsv

            rk_querysearch ${bin_oligos} ${bin_reads[1]} ${meta_reads.rlen} ${meta_oligos.expected_length} \\
                ${meta_oligos.n_oligos} ${left_allowed_shift} ${right_allowed_shift} \\
                ${meta_oligos.mismatches_allowed} \\
                > ${prefix}.R1.tsv

            echo $VERSION > ${software}.version.txt
            """
        }
    }
}
