// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '0.0'

process BED_ANNOTATE_RESTRICTION {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "${moduleDir}/environment_pyarrow_hdf5_cli.yml" : null)

    input:
    tuple val(meta), path(bed)
    tuple val(meta), path(table)
    tuple val(meta_restr), path(restriction_sites)

    output:
    tuple val(meta), path("*.pq"), emit: output
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}" // Ideally should contain: ${segment_name}.${renz_key}${renz_strand}

    def renz_strand = meta_restr['renz_strand']
    def renz_key = meta_restr['renzyme']
    def segment_name = meta['fragment']

    def renz_strand_key = (renz_strand=="+") ? "p" : (renz_strand=="-") ? "n" : ""  // p, n or empty
    def renz_strand_sub = (renz_strand=="+") ? "+" : (renz_strand=="-") ? "-" : "b"
    def columns = [
        "read_id_${segment_name}",
        "${segment_name}_start_${renz_key}${renz_strand_key}_left",
        "${segment_name}_start_${renz_key}${renz_strand_key}_right",
        "${segment_name}_end_${renz_key}${renz_strand_key}_left",
        "${segment_name}_end_${renz_key}${renz_strand_key}_right"
        ]
    def header = columns.join(",")

    """
    get_closest_sites.py -o ${prefix}.${renz_key}${renz_strand}.distances.pq \
        --strand ${renz_strand_sub} \
        --out-format parquet \
        --columns ${header} \
        --align-ids parquet::${table}::readID \
        ${bed} ${restriction_sites}

    echo $VERSION > ${software}.version.txt
    """
}
