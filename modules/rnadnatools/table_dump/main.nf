// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
params.options.args.input_format = params.options.args.getOrDefault('input_format', 'parquet')
params.options.args.output_format = params.options.args.getOrDefault('output_format', 'parquet')
options        = initOptions(params.options)

process RNADNATOOLS_TABLE_DUMP {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }
//    cache "${params.cache}"

    conda (params.enable_conda ? "${moduleDir}/../environment.yml" : null)

    input:
    tuple val(meta), path(input)
    val(filter)
    val(columns)

    output:
    tuple val(meta), path("*.${options.args.output_format}"), emit: table
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def chunksize = options.getOrDefault('chunksize', 1_000_000)

    def format = options.args.format

    def inputs = ""
    if (input instanceof List) {
            inputs = input.join(" ")
        }
    else {
        inputs = input
    }
    """
    rnadnatools table dump -i ${options.args.input_format} \
                           -o ${options.args.output_format} \
                           -f ${filter} \
                           -c ${columns.join(',')} \
                           ${prefix}.${options.args.output_format} ${inputs}

    python -c "import rnadnatools; print('rnadnatools', rnadnatools.__version__)" > ${software}.version.txt
    """
}
