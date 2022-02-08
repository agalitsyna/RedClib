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
//    cache "${params.cache}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}*.fq*"), emit: fastq
    path  "*.version.txt"         , emit: version


    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    /* Read input parameters: */
    def replacement_prefix = options.args.getOrDefault("ext_prefix", "")
    def replacement_suffix = options.args.getOrDefault("ext_suffix", "")
    /* Metadata overrides parameters: */
    replacement_prefix = meta.getOrDefault("ext_prefix", replacement_prefix)
    replacement_suffix = meta.getOrDefault("ext_suffix", replacement_suffix)

    /* Fill in the default quality: */
    def replacement_prefix_qual = "~"*replacement_prefix.size()
    def replacement_suffix_qual  = "~"*replacement_suffix.size()

    /* FASTQ reading command */
    def cmd_read_fastq = reads[0].getName().endsWith(".gz") ? "gzip -dc " : "cat "
    def cmd_write_fastq = options.args.getOrDefault("gzip", true) ? ' | gzip -c ' : ''
    def format = options.args.getOrDefault("gzip", true) ? '.gz' : ''

    def fastq_r1 = ""
    def fastq_r2 = ""
    if (meta.single_end) {
        fastq_r1 = reads[0]
//        def filename = fastq_r1.getName().replaceFirst(~/\.[^\.]+$/, '') // Get the filename without extension
        """
        ${cmd_read_fastq} ${fastq_r1} | \\
            sed '2~4s/^\\(.*\\)\$/${replacement_prefix}\\1${replacement_suffix}/' | \\
            sed '4~4s/^\\(.*\\)\$/${replacement_prefix_qual}\\1${replacement_suffix_qual}/' ${cmd_write_fastq} > ${prefix}.fq${format}

        echo $VERSION > ${software}.version.txt
        """
    } else {
        fastq_r1 = reads[0]
        fastq_r2 = reads[1]
//        def filename1 = fastq_r1.getName().replaceFirst(~/\.[^\.]+$/, '') // Get the filename without extension
//        def filename2 = fastq_r2.getName().replaceFirst(~/\.[^\.]+$/, '')
        """
        ${cmd_read_fastq} ${fastq_r1} | \\
            sed '2~4s/^\\(.*\\)\$/${replacement_prefix}\\1${replacement_suffix}/' | \\
            sed '4~4s/^\\(.*\\)\$/${replacement_prefix_qual}\\1${replacement_suffix_qual}/' ${cmd_write_fastq} > ${prefix}_1.fq${format}

        ${cmd_read_fastq} ${fastq_r2} | \\
            sed '2~4s/^\\(.*\\)\$/${replacement_prefix}\\1${replacement_suffix}/' | \\
            sed '4~4s/^\\(.*\\)\$/${replacement_prefix_qual}\\1${replacement_suffix_qual}/' ${cmd_write_fastq} > ${prefix}_2.fq${format}

        echo $VERSION > ${software}.version.txt
        """
    }
}
