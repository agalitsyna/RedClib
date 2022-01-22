// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
params.options.args.input_format = params.options.args.get('input_format', 'parquet')
params.options.args.output_format = params.options.args.get('output_format', 'parquet')
options        = initOptions(params.options)

process RNADNATOOLS_TABLE_MERGE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "${moduleDir}/../environment.yml" : null)

    input:
    val(meta)
    path(pqs)
    val(suffixes)

    output:
    tuple val(meta), path("*.${options.args.output_format}"), emit: table
    path  "*.version.txt"         , emit: version

    tuple val(meta), path("*fastq.gz"), optional:true, emit: fastq

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def chunksize = options.get('chunksize', 1_000_000)
    def filename = table.take(table.lastIndexOf('.'))

    """
    rnadnatools table merge -i ${options.args.input_format} \
        -o ${options.args.output_format} \
        ${filename}.${options.args.output_format} ${table}

    python -c "import rnadnatools; print('rnadnatools', rnadnatools.__version__)" > ${software}.version.txt
    """
}
