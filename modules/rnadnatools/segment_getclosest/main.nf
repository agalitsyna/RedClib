// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
params.options.args.output_format = params.options.args.getOrDefault('output_format', 'parquet')
options        = initOptions(params.options)

process RNADNATOOLS_SEGMENT_GETCLOSEST {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }
//    cache "${params.cache}"
    
    conda (params.enable_conda ? "${moduleDir}/../environment.yml" : null)

    input:
    tuple val(meta), path(bed)
    tuple val(meta_table), path(table)
    tuple val(meta_restr), path(restriction_sites)

    output:
    tuple val(meta), path("*.${options.args.output_format}"), emit: table
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}" // Ideally should contain: ${segment_name}.${renz_key}${renz_strand}

    def renz_strand = meta_restr['renz_strand']
    def renz_key = meta_restr['renzyme']
    def segment_name = meta['fragment']

    def renz_strand_key = (renz_strand=="+") ? "p" : (renz_strand=="-") ? "n" : ""  // p, n or empty
    def renz_strand_sub = (renz_strand=="+") ? "+" : (renz_strand=="-") ? "-" : "b"
    def renz_shortname = "${renz_key}${renz_strand_key}"
    def columns = [
        "readID_${segment_name}_${renz_shortname}",
        "${segment_name}_start_${renz_shortname}_left",
        "${segment_name}_start_${renz_shortname}_right",
        "${segment_name}_end_${renz_shortname}_left",
        "${segment_name}_end_${renz_shortname}_right"
        ]
    def header = columns.join(",")

    """
    rnadnatools segment get-closest-sites \
        -o ${prefix}.${renz_key}${renz_strand}.distances.${options.args.output_format} \
        --strand ${renz_strand_sub} \
        --out-format ${options.args.output_format} \
        --columns ${header} \
        --align-ids parquet::${table}::readID \
        ${bed} ${restriction_sites}

    python -c "import rnadnatools; print('rnadnatools', rnadnatools.__version__)" > ${software}.version.txt
    """
}
