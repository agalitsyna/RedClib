// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FASTUNIQ {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::fastuniq=1.1 anaconda::gawk=5.1.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/fastuniq:1.1--0"
    } else {
        container "quay.io/biocontainers/fastuniq:1.1--h470a237_1"
    }

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.nodup.txt"), emit: nodup_readids
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def VERSION = '1.1'

    if (meta.single_end) {
        def fastq_r1 = fastq[0]
        """
        # Create input file for fastuniq
        echo ${fastq_r1} >  filelist.txt

        # Run fastuniq
        fastuniq -i filelist.txt -tq -c 0 -o ${prefix}.1.unique.fq

        # Parse fastuniq output
        awk 'NR%4==1' ${prefix}.1.unique.fq | gawk '{match(\$0, "@([^ ,/]+)", a)} {print a[1]}' \
            > ${prefix}.nodup.txt

        echo $VERSION > ${software}.version.txt
        """
    } else {
        def fastq_r1 = fastq[0]
        def fastq_r2 = fastq[1]
        """
        # Create input file for fastuniq
        echo ${fastq_r1} >  filelist.txt
        echo ${fastq_r2} >> filelist.txt

        # Run fastuniq
        fastuniq -i filelist.txt -tq -c 0 -o ${prefix}.1.unique.fq -p ${prefix}.2.unique.fq

        # Parse fastuniq output
        awk 'NR%4==1' ${prefix}.1.unique.fq | gawk '{match(\$0, "@([^ ,/]+)", a)} {print a[1]}' \
            > ${prefix}.nodup.txt

        echo $VERSION > ${software}.version.txt
        """
    }
}