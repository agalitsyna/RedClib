// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '0.0'

process FASTQ_EXTEND {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path  "*.version.txt"         , emit: version


    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    /* Read input parameters: */
    def replacement_prefix = options.args.get("prefix", "")
    def replacement_suffix = options.args.get("suffix", "")

    /* Fill in the default quality: */
    def replacement_prefix_qual = "~"*replacement_prefix.size()
    def replacement_suffix_qual  = "~"*replacement_suffix.size()

    /* FASTQ reading command */
    def cmd_read_fastq = reads[0].getName().endsWith(".gz") ? "gzip -dc " : "cat "
    def cmd_write_fastq = options.args2.get("gzip", true) ? ' | gzip -c ' : ''
    def format = options.args2.get("gzip", true) ? '.gz' : ''

    if (meta.single_end) {
        def fastq_r1 = reads[0]
        """
        ${cmd_read_fastq} ${fastq_r1} | \\
            sed '2~4s/^\\(.*\\)\$/${replacement_prefix}\\1${replacement_suffix}/' | \\
            sed '4~4s/^\\(.*\\)\$/${replacement_prefix_qual}\\1${replacement_suffix_qual}/' ${cmd_write_fastq} > ${prefix}.extend.fastq${format}

        echo $VERSION > ${software}.version.txt
        """
    } else {
        def fastq_r1 = reads[0]
        def fastq_r2 = reads[1]
        """
        ${cmd_read_fastq} ${fastq_r1} | \\
            sed '2~4s/^\\(.*\\)\$/${replacement_prefix}\\1${replacement_suffix}/' | \\
            sed '4~4s/^\\(.*\\)\$/${replacement_prefix_qual}\\1${replacement_suffix_qual}/' ${cmd_write_fastq} > ${prefix}.extend_R1.fastq${format}

        ${cmd_read_fastq} ${fastq_r2} | \\
            sed '2~4s/^\\(.*\\)\$/${replacement_prefix}\\1${replacement_suffix}/' | \\
            sed '4~4s/^\\(.*\\)\$/${replacement_prefix_qual}\\1${replacement_suffix_qual}/' ${cmd_write_fastq} > ${prefix}.extend_R2.fastq${format}

        echo $VERSION > ${software}.version.txt
        """
    }
}
