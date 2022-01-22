// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '0.0'

process RKLIB_CHECK_COMPLEMENTARY {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(table)
    tuple val(meta_reads), path(bin_reads)
    tuple val(meta_oligos_bridge_forward), path(aligned_bridge_forward)
    tuple val(meta_oligos_ggg), path(aligned_ggg)

    output:
    tuple val(meta), path("*.complementaryHits.tsv"), emit: hits
    path  "*.version.txt"         , emit: version


    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def rna_complementary_length = options.args.rna_complementary_length
    def extraN = "N"*(500+3*rna_complementary_length)

    def bridge_length = options.args.bridge_length
    def read_length = meta.rlen

    """
    # Get the complemetary regions:

    paste <(awk '{print \$1, \$3, \$4}' ${table} | tail -n +2) \\
          <(head -n -1 ${aligned_bridge_forward} | tail -n +2 | awk '{print \$4}') \\
          | awk 'BEGIN{OFS="\\n";} {bgn=1; if (\$4<500) bgn=\$4+${bridge_length}+1; \\
          print \$1, substr(\$2"${extraN}", bgn, ${rna_complementary_length}), \\
          "+", substr(\$3"${extraN}", bgn, ${rna_complementary_length})}' > ${prefix}.rna-end.1.fq

    paste <(awk '{print \$1, \$5, \$6}' ${table} | tail -n +2) \\
          <(head -n -1 ${aligned_ggg} | tail -n +2 | awk '{print \$5+1}') \\
          | awk 'BEGIN{OFS="\\n";} {bgn=1; if (\$4<500) bgn=\$4; \\
          print \$1, substr(\$2"${extraN}", bgn, ${rna_complementary_length}), \\
          "+", substr(\$3"${extraN}", bgn, ${rna_complementary_length})}' > ${prefix}.rna-end.2.fq

    # Convert to reverse complement:

    paste <(sed -n '1~4p' ${prefix}.rna-end.1.fq) \\
          <(sed -n '2~4p' ${prefix}.rna-end.1.fq | rev | tr "ATGC" "TACG") \\
          <(sed -n '3~4p' ${prefix}.rna-end.1.fq) \\
          <(sed -n '4~4p' ${prefix}.rna-end.1.fq | rev) | tr "\\t" "\\n" \\
          > ${prefix}.rna-end.1.revcomp.fq

    paste <(sed -n '1~4p' ${prefix}.rna-end.2.fq) \\
          <(sed -n '2~4p' ${prefix}.rna-end.2.fq | rev | tr "ATGC" "TACG") \\
          <(sed -n '3~4p' ${prefix}.rna-end.2.fq) \\
          <(sed -n '4~4p' ${prefix}.rna-end.2.fq | rev) | tr "\\t" "\\n" \\
          > ${prefix}.rna-end.2.revcomp.fq

    # Convert to binary file

    fastq2hash ${prefix}.rna-end.1.revcomp.fq ${prefix}.rna-end.1.revcomp.bin
    fastq2hash ${prefix}.rna-end.2.revcomp.fq ${prefix}.rna-end.2.revcomp.bin

    # Align complementary regions to each library:
    rk_pairwise ${prefix}.rna-end.2.revcomp.bin ${bin_reads[0]} ${rna_complementary_length} ${read_length} \\
                   0 ${read_length-rna_complementary_length} 0 > ${prefix}.${meta_reads.id}.${meta_oligos_bridge_forward.id}.${meta_oligos_ggg.id}.hits_complementary.R1.tsv

    rk_pairwise ${prefix}.rna-end.1.revcomp.bin ${bin_reads[0]} ${rna_complementary_length} ${read_length} \\
               0 ${read_length-rna_complementary_length} 0 > ${prefix}.${meta_reads.id}.${meta_oligos_bridge_forward.id}.${meta_oligos_ggg.id}.hits_complementary.R2.tsv

    # Merge two tabular outputs:
    # design header
    head -n 1 ${prefix}.${meta_reads.id}.${meta_oligos_bridge_forward.id}.${meta_oligos_ggg.id}.hits_complementary.R1.tsv \\
     | sed 's/\\t/__complementary_R1\\t/g' | sed 's/\$/__complementary_R1/' | sed 's/#//g' | tr '\\n' '\\t' > ${prefix}.complementaryHits.tsv
    head -n 1 ${prefix}.${meta_reads.id}.${meta_oligos_bridge_forward.id}.${meta_oligos_ggg.id}.hits_complementary.R2.tsv \\
     | sed 's/\\t/__complementary_R2\\t/g' | sed 's/\$/__complementary_R2/' | sed 's/#//g' | tr '\\n' '\\t' >> ${prefix}.complementaryHits.tsv
    sed -i "s/\\t\$/\\n/" ${prefix}.complementaryHits.tsv
    # write to the body
    paste <(tail -n +2 ${prefix}.${meta_reads.id}.${meta_oligos_bridge_forward.id}.${meta_oligos_ggg.id}.hits_complementary.R1.tsv) \\
        <(tail -n +2 ${prefix}.${meta_reads.id}.${meta_oligos_bridge_forward.id}.${meta_oligos_ggg.id}.hits_complementary.R2.tsv) >> ${prefix}.complementaryHits.tsv

    echo $VERSION > ${software}.version.txt
    """
}

