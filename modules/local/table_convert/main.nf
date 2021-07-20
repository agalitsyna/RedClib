// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '0.0.0'

process TABLE_CONVERT {
    tag "$meta.id"
    label 'process_low'
    label 'no_parallel'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.txt"), emit: table
    path  "*.version.txt"               , emit: version

    tuple val(meta), path("*fastq.gz"), optional:true, emit: fastq

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    if (meta.single_end) {
         def input_fq1 = reads[0]
        """
        paste <(awk '{print \$1}' ${input_fq1} | sed 'N;N;N;s/\\n/ /g' | \
                awk 'BEGIN{OFS="\\t"}{print \$1, "${prefix}", \$2, \$4}' ) > ${prefix}.fastq.txt

        echo $VERSION > ${software}.version.txt
        """
    } else {
        def input_fq1 = reads[0]
        def input_fq2 = reads[1]
        """
        paste <(awk '{print \$1}' ${input_fq1} | sed 'N;N;N;s/\\n/ /g' | \
                awk 'BEGIN{OFS="\\t"}{print \$1, "${prefix}", \$2, \$4}' ) \
              <(awk '{print \$1}' ${input_fq2} | sed 'N;N;N;s/\\n/ /g' | \
                awk 'BEGIN{OFS="\\t"}{print \$2, \$4}' ) > ${prefix}.fastq.txt

        echo $VERSION > ${software}.version.txt
        """
    }
}
