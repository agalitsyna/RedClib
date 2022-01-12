// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BAM2BED {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bedtools=2.30 bioconda::samtools=1.14" : null)

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def cmd = ""
    def filter1 = params.options.get("filter", "")
    def filter2 = params.options.get("filter2", "")
    def filter3 = params.options.get("filter3", "")

    if (filter1) {
        cmd += "${filter1} ${bam} "
    }
    if (filter2) {
        cmd += " | ${filter2} "
    }
    if (filter3) {
        cmd += " | ${filter3} "
    }

    cmd += " | bedtools bamtobed -cigar -i stdin "

    """
    ${cmd} > ${prefix}.bed

    bedtools --version &> ${software}.version.txt
    samtools --version | head -n 2 > ${software}.version.txt
    """
}
