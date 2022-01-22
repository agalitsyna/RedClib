// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process RNADNATOOLS_GENOME_RECSITES {
    tag "$meta.assembly $meta.renzyme"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:'') }

    conda (params.enable_conda ? "${moduleDir}/../environment.yml" : null)

    input:
    val(meta)
    file(genome_fasta)

    output:
    tuple val(meta), path("${meta.assembly}.${meta.renzyme}.bed"), emit: genome_restricted
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def assembly = meta['assembly']
    def renzyme = meta['renzyme']
    """
    rnadnatools genome renzymes-recsites ${genome_fasta} ${renzyme} ${assembly}.${renzyme}.nonsorted.bed
    sort -k1,1 -k2,2n ${assembly}.${renzyme}.nonsorted.bed > ${assembly}.${renzyme}.bed
    rm ${assembly}.${renzyme}.nonsorted.bed

    python -c "import rnadnatools; print('rnadnatools', rnadnatools.__version__)" > ${software}.version.txt
    """
    }
