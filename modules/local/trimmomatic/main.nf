// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '0.39'

process TRIM_TRIMMOMATIC {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::trimmomatic=0.39" : null)

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.trimmed.fq"), emit: trimmed
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def params_trimmomatic = options.args.get('params_trimmomatic', "")
    if (meta.single_end) {
        def input_fq1 = reads[0]
        """
        trimmomatic SE -phred33 -threads ${task.cpus} ${input_fq1} \
            ${prefix}.1.trimmed.fq ${prefix}_1U ${params_trimmomatic}

        echo $VERSION > ${software}.version.txt
        """
    } else {
        def input_fq1 = reads[0]
        def input_fq2 = reads[1]
        """
        trimmomatic PE -phred33 -threads ${task.cpus} ${input_fq1} ${input_fq2} \
            ${prefix}.1.trimmed.fq ${prefix}_1U \
            ${prefix}.2.trimmed.fq ${prefix}_2U ${params_trimmomatic}

        echo $VERSION > ${software}.version.txt
        """
    }
}
