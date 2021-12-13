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
            saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta_reads, publish_by_meta:['id']) }

    input:
    tuple val(meta_oligos), path(bin_oligos)
    tuple val(meta_reads), path(bin_reads)

    output:
    tuple val(meta_reads), path("*.txt"), emit: aligned
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta_reads.id}.${meta_oligos.id}${options.suffix}" : "${meta_reads.id}.${meta_oligos.id}"
    meta_reads.oligo = meta_oligos.id

    // Encoding the length of the reads for C program:
    def lengths = [101:15, 151:21, 125:18, 80:12, 133:19, 251:34]

    if (!lengths.containsKey(meta_reads.rlen.toInteger())) {
        throw new Exception("Read length is not supported for this type of script")
    }
    
    def rlen_encoded = lengths[ meta_reads.rlen.toInteger() ]
    def left_allowed_shift = meta_oligos.left_allowed_shift
    def right_allowed_shift = (meta_oligos.right_allowed_shift > 0) ? meta_oligos.right_allowed_shift : meta_reads.rlen-meta_oligos.right_allowed_shift

    if (meta_reads.single_end) {
        """
        align_universal ${bin_oligos} ${bin_reads[0]} 1 ${meta_reads.rlen} ${rlen_encoded} \\
            ${meta_oligos.n_oligos} ${left_allowed_shift} ${right_allowed_shift} \\
            ${meta_oligos.mismateches_allowed} 0 ${meta_oligos.expected_length} \\
            > ${prefix}.txt

        echo $VERSION > ${software}.version.txt
        """
    } else {
        """
        align_universal ${bin_oligos} ${bin_reads[0]} 1 ${meta_reads.rlen} ${rlen_encoded} \\
            ${meta_oligos.n_oligos} ${left_allowed_shift} ${right_allowed_shift} \\
            ${meta_oligos.mismateches_allowed} 0 ${meta_oligos.expected_length} \\
            > ${prefix}.R1.txt

        align_universal ${bin_oligos} ${bin_reads[1]} 1 ${meta_reads.rlen} ${rlen_encoded} \\
            ${meta_oligos.n_oligos} ${left_allowed_shift} ${right_allowed_shift} \\
            ${meta_oligos.mismateches_allowed} 0 ${meta_oligos.expected_length} \\
            > ${prefix}.R2.txt

        echo $VERSION > ${software}.version.txt
        """
    }
}
