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

    conda (params.enable_conda ? "${moduleDir}/environment.yml" : null)

    input:
    tuple val(meta), val(entry)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path  "*.version.txt"              , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def sra_query = entry[0] // all the rest is ignored

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

    echo 'Python=3.7' > ${software}.version.txt
    echo 'sra-tools=2.11.0' >> ${software}.version.txt
    """

}
