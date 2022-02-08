// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FASTQ_DOWNLOAD {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

//    cache "${params.cache}"
    conda (params.enable_conda ? "${moduleDir}/environment.yml" : null)

    input:
    tuple val(meta), val(entry)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path  "*.version.txt"              , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def sra_query = entry[0]
    def sra_query2 = entry[1]
    def single_end = meta.get('single_end').toBoolean()

    if (single_end) {
    // Single-end data, look at the first entry in the list only:

        def fastqdumpCmd = ""
        def sra = ( sra_query=~ /SRR\d+/ )[0]
        def start = ( sra_query.contains('start=') ) ? ( sra_query =~ /start=(\d+)/ )[0][1] : 0
        def end   = ( sra_query.contains('end=') ) ? ( sra_query =~ /end=(\d+)/)[0][1] : 0
        if ((start>0) || (end>0)) {
            fastqdumpCmd += "fastq-dump ${sra} -Z --minSpotId ${start} --maxSpotId ${end}"
            }
        else {
            fastqdumpCmd += "fastq-dump ${sra} -Z"
        }

        """
        ${fastqdumpCmd} | bgzip -c -@${task.cpus} > ${prefix}.fastq.gz)

        fastq-dump -V >> ${software}.version.txt
        """

    } else {
    // Paired-end data:
        if (sra_query2=="") {
        // Both forward and reverse reads are provided in a single SRA:

            def fastqdumpCmd = ""
            def sra = ( sra_query=~ /SRR\d+/ )[0]
            def start = ( sra_query.contains('start=') ) ? ( sra_query =~ /start=(\d+)/ )[0][1] : 0
            def end   = ( sra_query.contains('end=') ) ? ( sra_query =~ /end=(\d+)/)[0][1] : 0
            if ((start>0) || (end>0)) {
                fastqdumpCmd += "fastq-dump ${sra} -Z --split-spot --minSpotId ${start} --maxSpotId ${end}"
                }
            else {
                fastqdumpCmd += "fastq-dump ${sra} -Z --split-spot"
            }

            """
            ${fastqdumpCmd} | pyfilesplit --lines 4 \
                             >(bgzip -c -@${task.cpus} > ${prefix}_1.fastq.gz) \
                             >(bgzip -c -@${task.cpus} > ${prefix}_2.fastq.gz) \
                             | cat

            fastq-dump -V >> ${software}.version.txt
            """

        } else {
        // Forward and reverse reads are provided as separate SRAs:

            def fastqdumpCmd1 = ""
            def fastqdumpCmd2 = ""

            def sra1 = ( sra_query=~ /SRR\d+/ )[0]
            def start1 = ( sra_query.contains('start=') ) ? ( sra_query =~ /start=(\d+)/ )[0][1] : 0
            def end1   = ( sra_query.contains('end=') ) ? ( sra_query =~ /end=(\d+)/)[0][1] : 0

            if ((start1>0) || (end1>0)) {
                fastqdumpCmd1 += "fastq-dump ${sra1} -Z --minSpotId ${start1} --maxSpotId ${end1}"
            } else {
                fastqdumpCmd1 += "fastq-dump ${sra1} -Z"
            }

            def sra2 = ( sra_query2=~ /SRR\d+/ )[0]
            def start2 = ( sra_query2.contains('start=') ) ? ( sra_query2 =~ /start=(\d+)/ )[0][1] : 0
            def end2   = ( sra_query2.contains('end=') ) ? ( sra_query2 =~ /end=(\d+)/)[0][1] : 0

            if ((start2>0) || (end2>0)) {
                fastqdumpCmd2 += "fastq-dump ${sra2} -Z --minSpotId ${start2} --maxSpotId ${end2}"
            } else {
                fastqdumpCmd2 += "fastq-dump ${sra2} -Z"
            }

            """
            ${fastqdumpCmd1} | bgzip -c -@${task.cpus} > ${prefix}_1.fastq.gz)

            ${fastqdumpCmd2} | bgzip -c -@${task.cpus} > ${prefix}_2.fastq.gz)

            fastq-dump -V >> ${software}.version.txt
            """

        }
    }

}
