// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '0.0'

process TABLE_TRIM {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(table)
    tuple val(meta_trim), path(trimmed_reads)

    output:
    tuple val(meta), path("*.trimtable.tsv"), emit: table
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def input_fq1 = trimmed_reads[0]
    def input_fq2 = trimmed_reads[1]

    """
    # Truncated read -> reads length:
    paste <(sed -n '1~4p' ${input_fq1} | awk '{print substr(\$1, 2);}') \\
          <(sed -n '2~4p' ${input_fq1} | awk '{print length(\$0);}') \\
          <(sed -n '2~4p' ${input_fq2} | awk '{print length(\$0);}') \\
          > ${prefix}.trim.info

    # Reads might be in peculiar order, we want to guarantee each read has a line in file in appropriate order:
    echo "#readID__trim\tpos_R1__trim\tpos_R2__trim" > ${prefix}.trimtable.tsv
    awk 'NR==FNR {vals[\$1] = \$1 "\\t" \$2 "\\t" \$3 ; next} \\
        !(\$1 in vals) {vals[\$1] = \$1 "\\t" "0\\t0"} \\
        {\$(NF+1) = vals[\$1]; print vals[\$1]}' ${prefix}.trim.info <(tail -n +2 ${table}) \\
        >> ${prefix}.trimtable.tsv

    echo $VERSION > ${software}.version.txt
    """
}
