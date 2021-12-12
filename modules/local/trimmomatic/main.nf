// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process TRIMMOMATIC {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::trimmomatic=0.39" : null)
//    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
//        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
//    } else {
//        container "quay.io/biocontainers/YOUR-TOOL-HERE"
//    }

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*fastq*"), emit: fastq
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def format = options.args2.get("gzip", true) ? '.gz' : ''

    if (meta.single_end) {
        def fastq_r1 = fastq[0]
        """
        trimmomatic SE \\
            -threads $task.cpus \\
            ${fastq_r1} \\
            ${prefix}.trim.fastq${format} \\
            $options.args

        trimmomatic -version > ${software}.version.txt
        """
    }  else {
        def fastq_r1 = fastq[0]
        def fastq_r2 = fastq[1]
        """
        trimmomatic PE \\
            -threads $task.cpus \\
            ${fastq_r1} \\
            ${fastq_r2} \\
            ${prefix}_R1.fastq${format} \\
            ${prefix}_R1.unpaired.fastq${format} \\
            ${prefix}_R2.fastq${format} \\
            ${prefix}_R2.unpaired.fastq${format} \\
            $options.args

        trimmomatic -version > ${software}.version.txt
        """
    }
}
