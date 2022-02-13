// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
params.options.args.input_format = params.options.args.getOrDefault('input_format', 'tsv')
params.options.args.ref_format = params.options.args.getOrDefault('ref_format', 'tsv')
params.options.args.output_format = params.options.args.getOrDefault('output_format', 'parquet')
options        = initOptions(params.options)

process RNADNATOOLS_TABLE_ALIGN {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

//    cache "${params.cache}"
    conda (params.enable_conda ? "${moduleDir}/../environment.yml" : null)

    input:
    tuple val(meta), path(table)
    tuple val(meta_ref), path(table_ref)

    output:
    tuple val(meta), path("*.${options.args.output_format}"), emit: table
    path  "*.version.txt"         , emit: version

    tuple val(meta), path("*fastq.gz"), optional:true, emit: fastq

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def outfile   = "${prefix}.${options.args.output_format}"

    def chunksize = options.args.getOrDefault('chunksize', 0)
    def chunkOpt = ""
//    if (chunksize) {
//        chunkOpt = "--chunksize ${chunksize}"
//    }
    def chunksize_writer = options.args.getOrDefault('chunksize_writer', 0)
    def chunkWriterOpt = ""
//    if (chunksize_writer) {
//        chunkWriterOpt = "--chunksize ${chunksize}"
//    }

    """
    rnadnatools table align -i ${options.args.input_format} \
                            -r ${options.args.ref_format} \
                            -o ${options.args.output_format} \
                            ${chunkOpt} \
                            ${chunkWriterOpt} \
                            ${options.args.params} \
                            ${table} ${table_ref} ${outfile}

    python -c "import rnadnatools; print('rnadnatools', rnadnatools.__version__)" > ${software}.version.txt
    """
}
