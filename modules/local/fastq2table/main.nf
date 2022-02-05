// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '0.0.0'

/* Useful Groovy methods */
// Return True if the file is gzipped
Boolean isGZ(line) {
    return (line.endsWith(".gz"))
}

process FASTQ2TSV {
    tag "$meta.id"
    label 'process_low'
    label 'no_parallel'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }
//    cache "${params.cache}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.tsv"), emit: table
    path  "*.version.txt"               , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def readCmd = ""

    if (meta.single_end) {
        def input_fq1 = reads[0]
        readCmd = (isGZ(input_fq1.toString())) ?  "bgzip -dc -@ ${task.cpus}" : "cat"
        """
        echo "#readID\tsample\tR\tQ\trlen"  > ${prefix}.fastq.tsv
        paste <(${readCmd} ${input_fq1} | awk '{print \$1}' | sed 'N;N;N;s/\\n/ /g' | \
                awk 'BEGIN{OFS="\\t"}{print substr(\$1,2), "${prefix}", \$2, \$4, length(\$2)}' ) >> ${prefix}.fastq.txt

        echo $VERSION > ${software}.version.txt
        """
    } else {
        def input_fq1 = reads[0]
        def input_fq2 = reads[1]
        readCmd1 = (isGZ(input_fq1.toString())) ?  "bgzip -dc -@ ${task.cpus}" : "cat"
        readCmd2 = (isGZ(input_fq2.toString())) ?  "bgzip -dc -@ ${task.cpus}" : "cat"
        """
        echo "#readID\tsample\tR1\tQ1\trlen1\tR2\tQ2\trlen2"  > ${prefix}.fastq.tsv
        paste <(${readCmd1} ${input_fq1} | awk '{print \$1}' | sed 'N;N;N;s/\\n/ /g' | \
                awk 'BEGIN{OFS="\\t"}{print substr(\$1,2), "${prefix}", \$2, \$4, length(\$2)}' ) \
              <(${readCmd2} ${input_fq2} | awk '{print \$1}' | sed 'N;N;N;s/\\n/ /g' | \
                awk 'BEGIN{OFS="\\t"}{print \$2, \$4, length(\$2)}' ) >> ${prefix}.fastq.tsv

        echo $VERSION > ${software}.version.txt
        """
    }
}
