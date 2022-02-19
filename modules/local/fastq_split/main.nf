// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

VERSION = '0.0'

/* Useful Groovy methods */
// Return True if the file is gzipped
Boolean isGZ(line) {
    return (line.endsWith(".gz"))
}

process FASTQ_SPLIT {
    tag "$meta.id"
    label 'process_low'
    label 'no_parallel'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

//    cache "$params.cache"
    conda (params.enable_conda ? "bioconda::tabix=1.11 conda-forge::coreutils=8.31" : null)

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.1.fq"), emit: fastq_r1
    tuple val(meta), path("*.2.fq"), emit: fastq_r2, optional: true
    path  "*.version.txt"        , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    // By default we set large chunks that almost guarantee no chunking:
    def chunksize = options.args.getOrDefault('chunksize', 10000000000)
    def suffix_length = options.args.getOrDefault('suffix_length', 4)

    def Cmd = ""

    if (meta.single_end) {

        def input_fq1 = reads
        def readCmd = (isGZ(input_fq1.toString())) ?  "bgzip -dc -@ ${task.cpus}" : "cat"
        Cmd = "${readCmd} ${input_fq1} | split -l ${chunksize} --suffix-length=${suffix_length} --numeric-suffixes=1 --additional-suffix='.1.fq' - ${prefix}."

    } else {

        def input_fq1 = reads[0]
        def input_fq2 = reads[1]
        def readCmd = (isGZ(input_fq1.toString())) ?  "bgzip -dc -@ ${task.cpus}" : "cat"
        Cmd = "${readCmd} ${input_fq1} | split -l ${chunksize} --suffix-length=${suffix_length} --numeric-suffixes=1 --additional-suffix='.1.fq' - ${prefix}.\n"
        Cmd += "${readCmd} ${input_fq2} | split -l ${chunksize} --suffix-length=${suffix_length} --numeric-suffixes=1 --additional-suffix='.2.fq' - ${prefix}."

    }

    """
    ${Cmd}

    echo 'fastq_split ${VERSION}' > ${software}.version.txt
    """
}