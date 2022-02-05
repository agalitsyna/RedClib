// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
params.options.args.input_format = params.options.args.getOrDefault('input_format', 'parquet')
params.options.args.output_format = params.options.args.getOrDefault('output_format', 'parquet')
options        = initOptions(params.options)

process RNADNATOOLS_TABLE_MERGE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

//    cache "${params.cache}"
    conda (params.enable_conda ? "${moduleDir}/../environment.yml" : null)

    input:
    tuple val(meta), path(tables)
    val(suffixes)

    output:
    tuple val(meta), path("*.${options.args.output_format}"), emit: table
    path  "*.version.txt"         , emit: version

    tuple val(meta), path("*fastq.gz"), optional:true, emit: fastq

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def chunksize = options.getOrDefault('chunksize', 1_000_000)

    def inputs = ""
    if (tables instanceof List) {
            inputs = tables.join(" ")
        }
    else {
        inputs = tables
    }

    def suffixes_full = []
    def colModifier = ""
    if (suffixes instanceof List) {
        if (suffixes.size()==1) {
            suffixes_full = suffixes.collect() * tables.size()
        } else {
            suffixes_full = suffixes.collect()
        }
        assert suffixes_full.size()==tables.size(): sprintf("Number of suffixes not equal to the number of files: %d vs %d", suffixes.size(), tables.size())
        suffixes_full = suffixes_full.collect { "{col_name}_"+it.toString() }
        colModifier = "--col-modifiers " + suffixes_full.join(",")
    }


    """
    rnadnatools table merge -i ${options.args.input_format} \
        -o ${options.args.output_format} \
        ${colModifier} \
        ${prefix}.${options.args.output_format} ${inputs}

    python -c "import rnadnatools; print('rnadnatools', rnadnatools.__version__)" > ${software}.version.txt
    """
}
