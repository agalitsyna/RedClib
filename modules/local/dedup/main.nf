// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '2.2.0'

/*
In RedC.nf I take an unusual approach and deduplicate before the mapping.
This decision was made because we need the statistics on bridge/adapters presence for already deduplicated reads.

Detect unique sequences based on N nucleotides at 5'-ends. N is defined by params.run.fastuniq_basepairs (50 by default)
Two-step procedure:
    1) crop first basepairs for both forward and reverse FASTQ file. Usually, errors are less frequent at 5'-end of read
    2) remove duplicated sequences with fastuniq
*/

process DEDUP {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::trimmomatic=0.39 bioconda::fastuniq=1.1 anaconda::gawk=5.1.0" : null)

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.ids.unique.txt"), emit: unique_index
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def cropped_length = options.args.dedup_crop_length

    if (meta.single_end) {
        def input_fq1 = reads[0]
        """
        # Trim N basepairs
        trimmomatic SE -threads ${task.cpus} ${input_fq1} \
                                             ${prefix}_1P ${prefix}_1U \
                                             CROP:${cropped_length}
        
        # Create input file for fastuniq
        echo ${prefix}_1P >  filelist.txt
        
        # Run fastuniq
        fastuniq -i filelist.txt -tq -c 0 -o ${prefix}.1.unique.fq
        
        # Parse fastuniq output
        awk 'NR%4==1' ${prefix}.1.unique.fq | gawk '{match(\$0, "@([^ ,/]+)", a)} {print a[1]}' \
            > ${prefix}.ids.unique.txt
        
        echo $VERSION > ${software}.version.txt
        """
    } else {
        def input_fq1 = reads[0]
        def input_fq2 = reads[1]
        """
        # Trim N basepairs
        trimmomatic PE -threads ${task.cpus} ${input_fq1} ${input_fq2} \
                                             ${prefix}_1P ${prefix}_1U \
                                             ${prefix}_2P ${prefix}_2U CROP:${cropped_length}
        
        # Create input file for fastuniq
        echo ${prefix}_1P >  filelist.txt
        echo ${prefix}_2P >> filelist.txt
        
        # Run fastuniq
        fastuniq -i filelist.txt -tq -c 0 -o ${prefix}.1.unique.fq -p ${prefix}.2.unique.fq
        
        # Parse fastuniq output
        awk 'NR%4==1' ${prefix}.1.unique.fq | gawk '{match(\$0, "@([^ ,/]+)", a)} {print a[1]}' \
            > ${prefix}.ids.unique.txt
        
        echo $VERSION > ${software}.version.txt
        """
    }
}
