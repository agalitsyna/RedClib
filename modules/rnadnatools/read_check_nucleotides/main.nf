// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process RNADNATOOLS_READ_CHECK_NUCLEOTIDES {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "${moduleDir}/../environment.yml" : null)

    input:
    tuple val(meta), path(table)
    tuple val(meta_oligos), path(aligned)

    output:
    tuple val(meta), path("*.tsv"), emit: hits
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    rnadnatools read check-nucleotides --oligo ${options.args.oligo} \
        -o ${prefix}.${meta_oligos.oligo}.${meta_oligos.id}.${options.args.oligo}.tsv \
        --readid-colname readID --seq-colname R1 \
        --position-colname start_hit --shift 35 \
        ${table} ${aligned}

    python -c "import rnadnatools; print('rnadnatools', rnadnatools.__version__)" > ${software}.version.txt
    """

}
