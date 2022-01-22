// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
params.options.args.input_format = params.options.args.get('input_format', 'parquet')
params.options.args.output_format = params.options.args.get('output_format', 'parquet')
options        = initOptions(params.options)

process RNADNATOOLS_TABLE_EVALUATE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "${moduleDir}/../environment.yml" : null)

    input:
    tuple val(meta), path(input)
    val(filters)

    output:
    tuple val(meta), path(input), path("*.${options.args.output_format}"), emit: output
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def format = options.args.format
    def cmd_create_evaluation_scheme = "printf '' > evaluation_schema.txt\n"
    for (filter in filters.keySet()) {
        cmd_create_evaluation_scheme += "printf \"${filter}\\t${format}\\t${filters[filter]}\\n\" >> evaluation_schema.txt\n"
    }

    def inputs = ""
    if (input instanceof List) {
            inputs = input.join(" ")
        }
    else {
        inputs = input
    }
    """
    ${cmd_create_evaluation_scheme}
    rnadnatools table evaluate -i ${options.args.input_format} \
                               -o ${options.args.output_format} \
                               evaluation_schema.txt ${prefix}.${options.args.output_format} ${inputs}

    python -c "import rnadnatools; print('rnadnatools', rnadnatools.__version__)" > ${software}.version.txt
    """
}
