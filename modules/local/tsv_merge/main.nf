// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '0.0'

process TSV_MERGE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

//    conda (params.enable_conda ? "${moduleDir}/environment_pyarrow.yml" : null)
//    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
//        container "https://depot.galaxyproject.org/singularity/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0"
//    } else {
//        container "quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0"
//    }

    input:
    tuple val(meta), path(files)
    val(suffixes)

    output:
    tuple val(meta), path("*.tsv"), emit: output
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def chunksize = options.get('chunksize', 1_000_000)

    header_command = ''
    for (v in [files, suffixes].transpose()) {
        file = v[0]
        suffix = v[1]
        header_command += "head -n 1 ${file} | sed 's/\\t/__${suffix}\\t/g' | sed 's/\$/__${suffix}/' | sed 's/#//g' | tr '\\n' '\\t' >> header.txt\n"
    }

    """
    touch header.txt
    printf '#' > header.txt
    ${header_command}
    sed -i "s/\\t\$/\\n/" header.txt

    cat header.txt > ${prefix}.OligosHits.tsv
    paste ${files} | tail -n +2 >> ${prefix}.OligosHits.tsv

    echo $VERSION > ${software}.version.txt
    """
}
