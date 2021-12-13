// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '0.0'

process BIN_ENCODE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bin"), emit: bin
    path  "*.version.txt"         , emit: version

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def mode   = options.args.get('mode', 'fasta')

    if (mode=='fasta') {
        def software = "fasta2bin"
        if (meta.single_end) {
            """
            fasta2hash ${reads[0]} ${prefix}.bin
            echo $VERSION > ${software}.version.txt
            """
        } else {
            """
            fasta2hash ${reads[0]} ${prefix}.R1.bin
            fasta2hash ${reads[1]} ${prefix}.R2.bin
            echo $VERSION > ${software}.version.txt
            """
        }
    } else {
        def software = "fastq2bin"
        if (meta.single_end) {
            """
            fastq2hash ${reads[0]} ${prefix}.bin
            echo $VERSION > ${software}.version.txt
            """
        } else {
            """
            fastq2hash ${reads[0]} ${prefix}.R1.bin
            fastq2hash ${reads[1]} ${prefix}.R2.bin
            echo $VERSION > ${software}.version.txt
            """
        }
    }
}
