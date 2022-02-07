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
//    cache "${params.cache}"

    input:
    tuple val(meta), path(table)
    tuple val(meta_reads), path(bin_reads)
    tuple val(meta_oligos_left), path(aligned_left)
    tuple val(meta_oligos_right), path(aligned_right)
    val(meta_compl)

    output:
    tuple val(meta), path("*.complHits.tsv"), emit: hits
    path  "*.version.txt"         , emit: version


    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def rna_compl_length = meta_compl.length as Integer
    def extraN = "N"*(500+3*rna_compl_length)

    def left_oligo_length = meta_oligos_left.expected_length as Integer
    def right_oligo_length = meta_oligos_right.expected_length as Integer
    def read_length = meta.rlen as Integer
    def max_shift_right = read_length-rna_compl_length

    // Define what columns have information on reads based in the input table:
    def table_cols_left = (meta_compl.left_side==1) ? " \$3, \$4" : "\$6, \$7"
    def table_cols_right = (meta_compl.right_side==1) ? " \$3, \$4" : "\$6, \$7"

    """
    # Take the ends of reference oligos and get potentially complemetary regions:

    paste <(awk '{print \$1, ${table_cols_left}}' ${table} | tail -n +2) \\
          <(head -n -1 ${aligned_left} | tail -n +2 | awk '{print \$5}') \\
          | awk 'BEGIN{OFS="\\n";} {bgn=1; if (\$4<500) bgn=\$4+${left_oligo_length}+1; \\
          print \$1, substr(\$2"${extraN}", bgn, ${rna_compl_length}), \\
          "+", substr(\$3"${extraN}", bgn, ${rna_compl_length})}' > ${prefix}.rna-end.1.fq

    paste <(awk '{print \$1, ${table_cols_right}}' ${table} | tail -n +2) \\
          <(head -n -1 ${aligned_right} | tail -n +2 | awk '{print \$5}') \\
          | awk 'BEGIN{OFS="\\n";} {bgn=1; if (\$4<500) bgn=\$4+${right_oligo_length}+1; \\
          print \$1, substr(\$2"${extraN}", bgn, ${rna_compl_length}), \\
          "+", substr(\$3"${extraN}", bgn, ${rna_compl_length})}' > ${prefix}.rna-end.2.fq

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

    # Align compl regions to each library, 1 mismatch allowed, screen from right to left:

    rk_pairwise ${prefix}.rna-end.2.revcomp.bin ${bin_reads[0]} ${rna_compl_length} ${read_length} \\
               0 ${max_shift_right} 1 0 > ${prefix}.hits_compl.R1.tsv

    rk_pairwise ${prefix}.rna-end.1.revcomp.bin ${bin_reads[1]} ${rna_compl_length} ${read_length} \\
               0 ${max_shift_right} 1 0 > ${prefix}.hits_compl.R2.tsv

    # Merge two tabular outputs:
    # 1. Design header
    head -n 1 ${prefix}.hits_compl.R1.tsv \\
     | sed 's/\\t/_compl_R1\\t/g' | sed 's/\$/_compl_R1/' | sed 's/#//g' | tr '\\n' '\\t' > ${prefix}.complHits.tsv
    head -n 1 ${prefix}.hits_compl.R2.tsv \\
     | sed 's/\\t/_compl_R2\\t/g' | sed 's/\$/_compl_R2/' | sed 's/#//g' | tr '\\n' '\\t' >> ${prefix}.complHits.tsv
    sed -i "s/\\t\$/\\n/" ${prefix}.complHits.tsv
    # 2. Write to the body
    paste <(tail -n +2 ${prefix}.hits_compl.R1.tsv) \\
        <(tail -n +2 ${prefix}.hits_compl.R2.tsv) >> ${prefix}.complHits.tsv

    echo $VERSION > ${software}.version.txt
    """
}

