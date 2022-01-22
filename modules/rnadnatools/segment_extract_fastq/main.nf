// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
params.options.args.input_format = params.options.args.get('input_format', 'tsv')
options        = initOptions(params.options)

process RNADNATOOLS_SEGMENT_EXTRACT_FASTQ {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "${moduleDir}/../environment.yml" : null)

    input:
    tuple val(meta), path(inputs)

    output:
    tuple val(meta), path("*.fq.gz"), path(inputs), emit: output
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def inputs_concatenated = inputs.join(' ')
    def fragment_name = meta.fragment
    def side = meta.side
    def selection_criteria = meta.selection_criteria
    """
    rnadnatools segment extract-fastq -i ${options.args.input_format} \
                                      -s "${selection_criteria}" \
                                      --key-start ${fragment_name}_start \
                                      --key-end ${fragment_name}_end \
                                      --key-readid readID \
                                      --key-seq R${side} \
                                      --key-qual Q${side} \
                                      ${prefix}.${fragment_name}.fq ${inputs_concatenated}

    gzip ${prefix}.${fragment_name}.fq

    python -c "import rnadnatools; print('rnadnatools', rnadnatools.__version__)" > ${software}.version.txt
    """
}
