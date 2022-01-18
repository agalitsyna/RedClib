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

    conda (params.enable_conda ? "anaconda::numpy anaconda::pandas" : null)

    input:
    tuple val(meta), path(bed)
    tuple val(meta_restr), path(restriction_sites)

    output:
    tuple val(meta), path("*.tsv"), emit: output
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}" // Ideally should contain: ${segment_name}.${renz_key}${renz_strand}

    def renz_strand = meta_restr['renz_strand']
    def renz_key = meta_restr['renzyme']
    def segment_name = meta['fragment']

    def renz_strand_key = (renz_strand=="+") ? "p" : (renz_strand=="-") ? "n" : ""  // p, n or empty
    def renz_strand_sub = (renz_strand=="+") ? "+" : (renz_strand=="-") ? "-" : "+" // + or - (also + is empty)
    def columns = [
        "${segment_name}_start_${renz_key}${renz_strand_key}_left",
        "${segment_name}_start_${renz_key}${renz_strand_key}_right",
        "${segment_name}_end_${renz_key}${renz_strand_key}_left",
        "${segment_name}_end_${renz_key}${renz_strand_key}_right"
        ]
    def header = (["id"]+columns).join(" ")

    """
    echo "${header}" > ${prefix}.${segment_name}.${renz_key}${renz_strand}.distances.txt
    get_closest_sites.py ${bed} ${restriction_sites} ${renz_strand_sub} ${prefix}.${renz_key}${renz_strand}.distances.tsv

    echo $VERSION > ${software}.version.txt
    """
}
