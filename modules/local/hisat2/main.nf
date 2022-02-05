// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '2.2.1'

process HISAT2_ALIGN {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    // Without lenient cache, long list of index files can occasionally produce wrong hash
    cache "lenient"

    conda (params.enable_conda ? "bioconda::hisat2=2.2.1 bioconda::samtools=1.10" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0"
    }

    input:
    tuple val(meta), path(reads)
    path  index
    path  splicesites

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.log"), emit: summary
    path  "*.version.txt"         , emit: version

    //tuple val(meta), path("*fastq.gz"), optional:true, emit: fastq

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def args     = options.args.tokenize() + meta.getOrDefault('mapping_args', '').tokenize()
    def args_str = options.args + ' ' + meta.getOrDefault('mapping_args', '')

    def strandedness = ''
    if (meta.strandedness == 'forward') {
        strandedness = meta.single_end ? '--rna-strandness F' : '--rna-strandness FR'
    } else if (meta.strandedness == 'reverse') {
        strandedness = meta.single_end ? '--rna-strandness R' : '--rna-strandness RF'
    }

    def splicesites_arg = (args_str.contains('--known-splicesite-infile') ) ?
            "--known-splicesite-infile ${splicesites}" : ''
    args.removeIf { it.contains('--known-splicesite-infile') }

    def seq_center_param = params.getOrDefault('seq_center', false)
    def seq_center = seq_center_param ? "--rg-id ${prefix} --rg SM:$prefix --rg CN:${seq_center_param.replaceAll('\\s','_')}" : "--rg-id ${prefix} --rg SM:$prefix"
    def unaligned = ""
    if (meta.single_end) {
        unaligned = params.getOrDefault('save_unaligned', false) ? "--un-gz ${prefix}.unmapped.fastq.gz" : ''
        """
        INDEX=`find -L ./ -name "*.1.ht2" | sed 's/.1.ht2//'`
        hisat2 \\
            -x \$INDEX \\
            -U $reads \\
            $strandedness \\
            $splicesites_arg \\
            --summary-file ${prefix}.hisat2.summary.log \\
            --threads $task.cpus \\
            $seq_center \\
            $unaligned \\
            ${args.join(' ')} \\
            | samtools view -bS -F 4 -F 256 - > ${prefix}.bam

        echo $VERSION > ${software}.version.txt
        """
    } else {
        unaligned = params.getOrDefault('save_unaligned', false) ? "--un-conc-gz ${prefix}.unmapped.fastq.gz" : ''
        """
        INDEX=`find -L ./ -name "*.1.ht2" | sed 's/.1.ht2//'`
        hisat2 \\
            -x \$INDEX \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            $strandedness \\
            $splicesites_arg \\
            --summary-file ${prefix}.hisat2.summary.log \\
            --threads $task.cpus \\
            $seq_center \\
            $unaligned \\
            --no-mixed \\
            --no-discordant \\
            ${args.join(' ')} \\
            | samtools view -bS -F 4 -F 8 -F 256 - > ${prefix}.bam

        echo $VERSION > ${software}.version.txt
        """
    }
}

//        if [ -f ${prefix}.unmapped.fastq.1.gz ]; then
//            mv ${prefix}.unmapped.fastq.1.gz ${prefix}.unmapped_1.fastq.gz
//        fi
//        if [ -f ${prefix}.unmapped.fastq.2.gz ]; then
//            mv ${prefix}.unmapped.fastq.2.gz ${prefix}.unmapped_2.fastq.gz
//        fi
