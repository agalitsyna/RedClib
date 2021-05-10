#!/usr/bin/env nextflow

// Run test project with:
// nextflow redc.nf -params-file project.yml

def helpMessage() {
    log.info"""
    RedC-nextflow is a pipeline to map RNA-DNA interactions in paired-end sequencing data.

    Usage:
    The typical command for launching the pipeline:
      nextflow redc.nf -params-file project.yml

    All the parameters should be listed in the project file.
    The explanation is provided with the default example.
    """.stripIndent()
}

// Show help message
if (params.get('help', 'false').toBoolean()) {
    helpMessage()
    exit 0
}

/* Useful Groovy methods */
// Get the full path of the output directory
String getOutputDir(output_type) {
    new File(params.output.dirs.get(output_type, output_type)).getCanonicalPath()
}
// Return True if the file is URL
Boolean isURL(line) {
    return (line.startsWith("http://") || line.startsWith("https://") || line.startsWith("ftp://"))
}
// Return True if the file is gzipped
Boolean isGZ(line) {
    return (line.endsWith(".gz"))
}
// Return True if the file is sra
Boolean isSRA(line) {
    return (line.startsWith("sra:"))
}
// Parse forward (left) and reverse (right) pieces of the file,
// assert that they correspond each other and return index
String parseChunkPair(left_file, right_file) {
    chunk_left_idx  =  left_file.toString().tokenize('.')[-3]
    chunk_right_idx = right_file.toString().tokenize('.')[-3]
    chunk_left_side  =  left_file.toString().tokenize('.')[-2]
    chunk_right_side = right_file.toString().tokenize('.')[-2]
    assert chunk_left_idx == chunk_right_idx
    assert chunk_left_side != chunk_right_side
    return chunk_left_idx
}

/* Parameters for the run */

def check_restriction = params.run.get('check_restriction', 'false').toBoolean()
def dna_extension = params.run.get('dna_extension', '')
def auto_download_genome = params.genome.get('auto_download_genome', 'true').toBoolean()
def available_renz = []

/////////////////////////////////////
/*   PREPARE FOLDER WITH BINARIES  */
/////////////////////////////////////

def align_universal = new File("bin/align_universal")
def align_pairwise  = new File("bin/align_pairwise")
def fasta2bin = new File("bin/fasta2bin")
def fastq2bin = new File("bin/fastq2bin")
def trimmomatic = new File("bin/Trimmomatic-0.39/trimmomatic-0.39.jar")
def missing_executables = ((!align_universal.exists()) || (!align_pairwise.exists()) ||
    (!fasta2bin.exists()) || (!fastq2bin.exists()) || (!trimmomatic.exists()))

def preprocessingCmd = ""
if (missing_executables) {
    preprocessingCmd = "bash ./bin/prepare_binaries.sh"
}

Channel.from(
    preprocessingCmd.execute().text
).into{PREPROCESSING_TO_RNA; PREPROCESSING}

///////////////////////////
/* Running the processes */
///////////////////////////

////////////////////////////
/*   PREPARE THE GENOME   */
////////////////////////////

GENOME_ASSEMBLY = params.genome.get('assembly_name', 'genome')
process DOWNLOAD_GENOME{
    tag "${assembly}"
    storeDir getOutputDir('genome')

    input:
    val assembly from GENOME_ASSEMBLY

    output:
    file "${assembly}.fa" into GENOME_FASTA
    file "${assembly}.chromsizes.txt" into GENOME_CHROMSIZES
    //file "${assembly}.fa.*" into GENOME_INDEX
    set "${assembly}.fa", file("${assembly}.fa.*") into GENOME_INDEX

    script:
    if (auto_download_genome) {
        """
        wget http://hgdownload.cse.ucsc.edu/goldenPath/${assembly}/bigZips/${assembly}.fa.gz -O ${assembly}.fa.gz
        bgzip -d -@ ${task.cpus} ${assembly}.fa.gz
        faidx ${assembly}.fa -i chromsizes > ${assembly}.chromsizes.txt
        hisat2-build -p ${task.cpus} ${assembly}.fa ${assembly}.fa
        """
    }
    else {
        // Download/copy genomic fasta
        def fasta_file = params.genome.get("fasta", "")
        def suffix = isGZ(fasta_file) ? ".gz" : ""

        // Get FASTA from URL or copy
        def getGenomeCmd = ""
        if (isURL(fasta_file)){
            getGenomeCmd = """
                wget ${fasta_file} -O ${assembly}.fa${suffix}
            """
        }
        else {
            getGenomeCmd = """
                cp ${fasta_file} ${assembly}.fa${suffix}
            """
        }

        // Unpack genomic fasta
        def unpackGenomeCmd = ""
        if (isGZ(fasta_file)){
            unpackGenomeCmd = """
                bgzip -d -@ ${task.cpus} ${assembly}.fa.gz
            """
        }

        // Download/copy/build chromosome sizes
        def chromsizes_file = params.genome.get("chromsizes", "")
        def getChromsizesCmd = ""
        if (isURL(chromsizes_file)){
            getChromsizesCmd = """
                wget ${chromsizes_file} -O ${assembly}.chromsizes.txt
            """
        }
        else if (chromsizes_file.size()>0) {
            getChromsizesCmd = """
                cp ${chromsizes_file} ${assembly}.chromsizes.txt
            """
        }
        else {
            getChromsizesCmd = """
                faidx ${assembly}.fa -i chromsizes > ${assembly}.chromsizes.txt
            """
        }

        // Run index if it does not exist or fasta file was downloaded from URL
        suffix_list = [".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2", ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2"]
        def getIndexCmd = ""
        if ((params.genome.get('index_prefix', '').size()==0) | (isURL(fasta_file))) {
            getIndexCmd = "hisat2-build -p ${task.cpus} ${assembly}.fa ${assembly}.fa"
        } else {
            for (suf in suffix_list) {
                getIndexCmd = getIndexCmd + "cp ${params.genome.index_prefix}${suf} ${assembly}.fa${suf}; "
            }
        }

        """
        ${getGenomeCmd}
        ${unpackGenomeCmd}
        ${getChromsizesCmd}
        ${getIndexCmd}
        """
    }
}

if (check_restriction) {
    LIST_RENZ = Channel.fromList(params.protocol.renz.collect{k, v -> [k, v]})
    process RESTRICT_GENOME{
        tag "${assembly} ${renz}"
        storeDir getOutputDir('genome')

        input:
        val(assembly) from GENOME_ASSEMBLY
        file(genome_fasta) from GENOME_FASTA // "${assembly}.fa.gz"
        set renz_key, renz from LIST_RENZ

        output:
        set renz_key, "${assembly}.${renz}.bed" into GENOME_RENZ

        script:
        """
        detect_restriction_sites.py ${genome_fasta} ${renz} ${assembly}.${renz}.nonsorted.bed
        sort -k1,1 -k2,2n --parallel=${task.cpus} ${assembly}.${renz}.nonsorted.bed > ${assembly}.${renz}.bed
        """
    }
}

////////////////////////////////
/*   PREPARE RNA ANNOTATION   */
////////////////////////////////

Channel.from(params.rna_annotation.get('rna_annotation_name', 'rna'))
       .combine(PREPROCESSING_TO_RNA).set{GENOME_RNA_ANNOT_NAME}

process PREPARE_RNA_ANNOTATION{
    tag "${rna_annot_name}"
    storeDir getOutputDir('genome')

    input:
    set val(rna_annot_name), val(preprocessing_output) from GENOME_RNA_ANNOT_NAME

    output:
    file "${rna_annot_name}.spliced_genes.txt" into GENOME_SPLICESITES
    file "${rna_annot_name}.gtf" into RNA_ANNOT_FILE

    script:
    def suffix = isGZ(params.rna_annotation.genes_gtf) ? ".gz" : ""
    def getRNAAnnot = ""
    if (isURL(params.rna_annotation.genes_gtf)) {
        getRNAAnnot = """
        wget ${params.rna_annotation.genes_gtf} -O ${rna_annot_name}.gtf${suffix}
        """
    }
    else {
        getRNAAnnot = """
        cp ${params.rna_annotation.genes_gtf} > ${rna_annot_name}.gtf${suffix}
        """
    }
    def unpackRNAAnnot = ""
    if (isGZ(params.rna_annotation.genes_gtf)){
        unpackRNAAnnot = """
            bgzip -d -@ ${task.cpus} ${rna_annot_name}.gtf.gz
        """
    }

    """
    ${getRNAAnnot}
    ${unpackRNAAnnot}
    hisat2_extract_splice_sites.py ${rna_annot_name}.gtf > ${rna_annot_name}.spliced_genes.txt
    """
}

////////////////////////////////////////////////
/*      PREPARE FASTQ FILES AND READ TABLE    */
////////////////////////////////////////////////
/*
download fastq from SRA, split into chunks.
*/

Channel.from(params.input.fastq_paths
    .collect{k, v ->  [k, v]}
    )
    .branch{
        for_download: isSRA(it[1][0])
        local:       !isSRA(it[1][0])
    }.set{FASTQ_PATHS}

/*
 * STEP 0: Download SRA
 */

process DOWNLOAD_FASTQ {
    tag "library:${library}"

    storeDir getOutputDir('fastq')

    input:
    tuple val(library), val(name) from FASTQ_PATHS.for_download

    output:
    tuple val(library), "${library}_1.fastq.gz", "${library}_2.fastq.gz" into DOWNLOADED

    script:
    def sra = ( name=~ /SRR\d+/ )[0]
    """
    fastq-dump ${sra} -Z --split-spot \
                   | pyfilesplit --lines 4 \
                     >(bgzip -c -@${task.cpus} > ${library}_1.fastq.gz) \
                     >(bgzip -c -@${task.cpus} > ${library}_2.fastq.gz) \
                     | cat
    """

}

FASTQ_PATHS.local.map { it -> [ it[0], file(it[1][0]), file(it[1][1]) ] }
    .concat(DOWNLOADED).into{ LIB_FASTQ; LIB_VIEW;
                              LIB_FASTQ_TO_FASTUNIQ }

//Channel.from(
//    params.input.fastq_paths.collect{k, v ->  [k, v[0], file(v[0]), v[1], file(v[1])]}
//    ).set{LIB_FASTQ}
//Channel.from(
//    params.input.fastq_paths.collect{k, v ->  [k, file(v[0]), file(v[1])]}
//    ).set{LIB_FASTQ_TO_FASTUNIQ}
//
def chunksize = params.run.chunksize*4

process SPLIT_FASTQ_INTO_CHUNKS{
    tag "library:${library}"

    storeDir getOutputDir('fastq')

    input:
    set val(library), path(input_fq1), path(input_fq2) from LIB_FASTQ

    output:
    set val(library), "${library}.*.1.fq", "${library}.*.2.fq" into LIB_SPLIT_FASTQ_RAW

    script:
    def readCmd = (isGZ(input_fq1.toString())) ?  "bgzip -dc -@ ${task.cpus}" : "cat"

    """
    echo
    ${readCmd} ${input_fq1} | split -l ${chunksize} --numeric-suffixes=1 \
        --additional-suffix=".1.fq" - ${library}.
    ${readCmd} ${input_fq2} | split -l ${chunksize} --numeric-suffixes=1 \
        --additional-suffix=".2.fq" - ${library}.
    """
}

LIB_SPLIT_FASTQ_RAW
    .transpose()
    .map{[it[0],
          parseChunkPair(it[1],it[2]), // index of the chunk (checked for safety)
          it[1],
          file(it[1]),
          it[2],
          file(it[2])]}
    .into{ LIB_SPLIT_FASTQ_TO_TABLE;
           LIB_SPLIT_FASTQ_TO_TRIM;
           LIB_SPLIT_FASTQ_TO_CINDEX }

/*
Merge fastq files into read table.
Each line in the resulting table corresponds to a single read.
This step saves time for further processing.
Next, it will be processed with bash-only approach.
*/

process CREATE_READS_TABLE_CHUNKS{
    tag "library:${library} chunk:${chunk}"

    storeDir getOutputDir('table')

    input:
    set val(library), val(chunk), val(input1), file(input_fq1), val(input2), file(input_fq2) from LIB_SPLIT_FASTQ_TO_TABLE

    output:
    set library, chunk, "${library}.${chunk}.fastq.txt" into LIB_TABLE_FASTQ

    script:
    """
    paste <(awk '{print \$1}' ${input_fq1} | sed 'N;N;N;s/\\n/ /g' | \
            awk 'BEGIN{OFS="\\t"}{print \$1, "${library}.${chunk}", \$2, \$4}' ) \
          <(awk '{print \$1}' ${input_fq2} | sed 'N;N;N;s/\\n/ /g' | \
            awk 'BEGIN{OFS="\\t"}{print \$2, \$4}' ) > ${library}.${chunk}.fastq.txt
    """
}

LIB_TABLE_FASTQ.into{ LIB_TABLE_FASTQ_FOR_TRIM;
                      LIB_TABLE_FASTQ_FOR_GA;
                      LIB_TABLE_FASTQ_FOR_RNACOMP;
                      LIB_TABLE_FASTQ_FOR_SUBSTR;
                      LIB_TABLE_FASTQ_FOR_COLLECT}

////////////////////////////////
/*       DEDUPLICATION        */
////////////////////////////////
/*
In RedC.nf I take an unusual approach and deduplicate before the mapping.
This decision was made because we need the statistics on bridge/adapters presence for already deduplicated reads.

Detect unique sequences based on N nucleotides at 5'-ends. N is defined by params.run.fastuniq_basepairs (50 by default)
Two-step procedure:
    1) crop first basepairs for both forward and reverse FASTQ file. Usually, errors are less frequent at 5'-end of read
    2) remove duplicated sequences with fastuniq
*/

def cropped_length = params.run.fastuniq_crop

process DEDUP{
    tag "library:${library}"

    storeDir getOutputDir('table')

    input:
    set val(library), file(input_fq1), file(input_fq2) from LIB_FASTQ_TO_FASTUNIQ

    output:
    set library, "${library}.ids.unique.txt" into IDS_FASTUNIQ

    script:
    """
    # Trim N basepairs
    trimmomatic PE -threads ${task.cpus} ${input_fq1} ${input_fq2} \
                                         ${library}_1P ${library}_1U \
                                         ${library}_2P ${library}_2U CROP:${cropped_length}

    # Create input file for fastuniq
    echo ${library}_1P >  filelist.txt
    echo ${library}_2P >> filelist.txt

    # Run fastuniq
    fastuniq -i filelist.txt -tq -c 0 -o ${library}.1.unique.fq -p ${library}.2.unique.fq

    # Parse fastuniq output
    awk 'NR%4==1' ${library}.1.unique.fq | gawk '{match(\$0, "@([^ ,/]+)", a)} {print a[1]}' \
        > ${library}.ids.unique.txt

    """
}

////////////////////////////////
/*          TRIMMING          */
////////////////////////////////
/*
Trim reads by quality with trimmomatic.
*/

def params_trimmomatic = params.run.params_trimmomatic

process TRIM_CHUNKS{
    tag "library:${library} chunk:${chunk}"

    storeDir getOutputDir('table')

    input:
    set val(library), val(chunk),
        val(input1), file(input_fq1),
        val(input2), file(input_fq2) from LIB_SPLIT_FASTQ_TO_TRIM

    output:
    set library, chunk, "${library}.${chunk}.1.trimmed.fq", "${library}.${chunk}.2.trimmed.fq" into LIB_TRIMMED

    script:
    """
    # Trim with specified parameters
    trimmomatic PE -phred33 -threads ${task.cpus} ${input_fq1} ${input_fq2} \
                            ${library}.${chunk}.1.trimmed.fq ${library}_1U \
                            ${library}.${chunk}.2.trimmed.fq ${library}_2U ${params_trimmomatic}
    """
}

// Structure of the channel:
// library, chunk, processed_table, processed_fq1, processed_fq2
LIB_TABLE_FASTQ_FOR_TRIM
    .combine(LIB_TRIMMED, by: [0, 1] )
    .map { it ->  [it[0], it[1], file(it[2]), file(it[3]), file(it[4]) ] }
    .set { LIB_FOR_GET_TRIM_OUTPUT }

process CREATE_TRIM_TABLE_CHUNKS{
    tag "library:${library} chunk:${chunk}"

    storeDir getOutputDir('table')

    input:
    set val(library), val(chunk), file(input_table), file(input_fq1), file(input_fq2) from LIB_FOR_GET_TRIM_OUTPUT

    output:
    set library, chunk, "${library}.${chunk}.trimtable.txt" into LIB_TRIMTABLE

    script:
    """
    paste <(sed -n '1~4p' ${input_fq1} | awk 'BEGIN{OFS="\\t"}{print "${library}.${chunk}", "trimmomatic", \$1}') \
          <(sed -n '2~4p' ${input_fq1} | awk 'BEGIN{OFS="\\t"}{print 0, length(\$0);}') \
          <(sed -n '2~4p' ${input_fq2} | awk 'BEGIN{OFS="\\t"}{print 0, length(\$0);}') \
          > ${library}.${chunk}.trim.info

    awk 'NR==FNR {vals[\$3] = \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 "\\t" \$6 "\\t" \$7 ; next} \
        !(\$1 in vals) {vals[\$1] = "${library}.${chunk}" "\\t" "trimmomatic" "\\t" \$1 "\\t" "0\\t0\\t0\\t0"} \
        {\$(NF+1) = vals[\$1]; print vals[\$1]}' ${library}.${chunk}.trim.info ${input_table} \
        > ${library}.${chunk}.trimtable.txt
    """
}

////////////////////////////////
/*      OLIGOS ALIGNMENT      */
////////////////////////////////
/*
For the alignment of oligos, we run the custom mapper based on Karp–Rabin algorithm.
Input fastq files with reads and files with oligos are indexed first, and then processed with custon C aligner.
*/

LIB_OLIGOS_RAW = Channel.from(
    params.input.oligos.collect{k, v ->  [k, v, file(v)]}
    ).combine(PREPROCESSING)

process INDEX_OLIGOS{
    tag "oligo:${oligo}"

    storeDir getOutputDir('cindex')

    input:
    set val(oligo), input_oligo, file(input_oligo_fa), val(preprocessing_output) from LIB_OLIGOS_RAW

    output:
    set oligo, "${oligo}.bin" into LIB_OLIGOS_CINDEX

    script:
    """
    fasta2bin ${input_oligo_fa} ${oligo}.bin
    """
}

process INDEX_CHUNKS{
    tag "library:${library} chunk:${chunk}"

    storeDir getOutputDir('cindex')

    input:
    set val(library), val(chunk),
        val(input1), file(input_fq1),
        val(input2), file(input_fq2) from LIB_SPLIT_FASTQ_TO_CINDEX

    output:
    set library, chunk, "${library}.${chunk}.1.bin", "${library}.${chunk}.2.bin" into LIB_FASTQ_CINDEX

    script:
    """
    fastq2bin ${input_fq1} ${library}.${chunk}.1.bin
    fastq2bin ${input_fq2} ${library}.${chunk}.2.bin
    """
}

def br_length = params.protocol.bridge_length

// Read length for each library
Channel.fromList(params.protocol.read_length.collect{k, v -> [k, v]})
    .into{ LIST_RLENGTHS_COLLECTION1; LIST_RLENGTHS_COLLECTION2}

// def read_length = params.protocol.read_length

// Set required oligos mapping and parameters of mapping calls:
MAPPING_COLLECTION1 = LIST_RLENGTHS_COLLECTION1.flatMap { lib, read_length ->
    [[library: lib, oligo: "adaptor_forward", apply_to:1, right_shift:read_length-14, read_length:read_length,
    n_primers:2, left_shift:-6, mismatch_general:1, report_len:20],
    [library: lib, oligo: "adaptor_reverse", apply_to:1, right_shift: read_length-14, read_length:read_length,
    n_primers:2, left_shift:-6, mismatch_general:1, report_len:20],
    [library: lib, oligo: "adaptor_reverse_short", apply_to:1, right_shift:read_length-14, read_length:read_length,
    n_primers:1, left_shift:0, mismatch_general:1, report_len:16],
    [library: lib, oligo: "bridge_forward", apply_to:1, right_shift:read_length-14, read_length:read_length,
    n_primers:1, left_shift:0, mismatch_general:1, report_len:br_length],
    [library: lib, oligo: "bridge_reverse", apply_to:1, right_shift:read_length-14, read_length:read_length,
    n_primers:1, left_shift:0, mismatch_general:1, report_len:br_length]]
}

MAPPING_COLLECTION2 = LIST_RLENGTHS_COLLECTION2.flatMap { lib, read_length ->
    [[library: lib, oligo: "adaptor_forward", apply_to:2, right_shift:read_length-14, read_length:read_length,
    n_primers:2, left_shift:-6, mismatch_general:1, report_len:20],
    [library: lib, oligo: "adaptor_reverse", apply_to:2, right_shift:read_length-14, read_length:read_length,
    n_primers:2, left_shift:-6, mismatch_general:1, report_len:20],
    [library: lib, oligo: "bridge_forward", apply_to:2, right_shift:read_length-14, read_length:read_length,
    n_primers:1, left_shift:0, mismatch_general:1, report_len:br_length],
    [library: lib, oligo: "bridge_reverse", apply_to:2, right_shift:read_length-14, read_length:read_length,
    n_primers:1, left_shift:0, mismatch_general:1, report_len:br_length],
    [library: lib, oligo: "ggg", apply_to:2, right_shift:3, read_length:read_length,
    n_primers:1, left_shift:0, mismatch_general:0, report_len:3]]
}

// Split channels for forward and reverse sides of read:
LIB_OLIGOS_CINDEX.into { LIB_OLIGOS_CINDEX1; LIB_OLIGOS_CINDEX2 }
LIB_FASTQ_CINDEX.into { LIB_FASTQ_CINDEX1; LIB_FASTQ_CINDEX2;
                        LIB_FASTQ_CINDEX_FOR_RNACOMP;
                        LIB_FASTQ_CINDEX_FOR_SUBSTR}


// 0          1            2       3     4          5          6
// oligo_name olifo_cindex library chunk fq1_cindex fq2_cindex params
LIB_OLIGOS_CINDEX1.combine(LIB_FASTQ_CINDEX1).combine(MAPPING_COLLECTION1)
    .filter { it[0]==it[6]["oligo"] && it[2]==it[6]["library"] }
    .map { it ->  [it[0], file(it[1]), it[2], it[3], file(it[4]), it[6].apply_to, it[6], it[6]["read_length"] ] }
    .set { LIB_FOR_OLIGOS_MAPPING1 }

LIB_OLIGOS_CINDEX2.combine(LIB_FASTQ_CINDEX2).combine(MAPPING_COLLECTION2)
    .filter { it[0]==it[6]["oligo"] && it[2]==it[6]["library"] }
    .map { it ->  [it[0], file(it[1]), it[2], it[3], file(it[5]), it[6].apply_to, it[6], it[6]["read_length"]  ] }
    .set { LIB_FOR_OLIGOS_MAPPING2 }

LIB_FOR_OLIGOS_MAPPING = LIB_FOR_OLIGOS_MAPPING1.concat(LIB_FOR_OLIGOS_MAPPING2)



// Encoding the length of the reads for C program:
def lengths = [101:15, 151:21, 125:18, 80:12, 133:19, 251:34]

process SEARCH_OLIGOS_CHUNKS{
    tag "library:${library} chunk:${chunk} side:${apply_to} oligo:${oligo}"

    storeDir getOutputDir('cout')

    input:
    set val(oligo), file(oligo_cindex), val(library), val(chunk),
        file(fq_cindex), val(apply_to), val(map_params), val(read_length) from LIB_FOR_OLIGOS_MAPPING

    output:
    set library, chunk, oligo, apply_to, "${library}.${chunk}.${apply_to}.${oligo}.txt", read_length into LIB_MAPPED_OLIGOS

    script:
    def seqlen_converted = lengths[ map_params.read_length ]
    """
    align_universal ${oligo_cindex} ${fq_cindex} 1 ${map_params.read_length} ${seqlen_converted} \
      ${map_params.n_primers} ${map_params.left_shift} ${map_params.right_shift} \
      ${map_params.mismatch_general} 0 ${map_params.report_len} > ${library}.${chunk}.${apply_to}.${oligo}.txt
    """
}

LIB_MAPPED_OLIGOS.into{ LIB_MAPPED_OLIGOS_FOR_GA;
                        LIB_MAPPED_OLIGOS_FOR_RNACOMP;
                        LIB_MAPPED_OLIGOS_FOR_SUBSTR;
                        LIB_MAPPED_OLIGOS_FOR_COLLECT}

// Resulting Channel structure: library, chunk, fastq_table, mapped_oligo
LIB_TABLE_FASTQ_FOR_GA
    .combine( LIB_MAPPED_OLIGOS_FOR_GA, by: [0, 1])
    .filter{ it[3]=="bridge_forward" && it[4]==1 }
    .map{ [it[0], it[1], file(it[2]), file(it[5])]  }
    .set{ LIB_FOR_GA }

process CHECK_GA{
    tag "library:${library} chunk:${chunk}"

    storeDir getOutputDir('table')

    input:
    set val(library), val(chunk), file(table_fq), file(cout_br_for) from LIB_FOR_GA

    output:
    set library, chunk, "${library}.${chunk}.GA.txt" into LIB_MAPPED_GA

    script:
    def checkGACmd=""
    if (params.run.get('check_GA', true)) {
        checkGACmd="check_oligo_presence.py ${table_fq} ${cout_br_for} ${library}.${chunk}.GA.txt"
    } else {
        checkGACmd="awk '{print NR-1\"\\t0\"}' ${table_fq} > ${library}.${chunk}.GA.txt"
    }
    """
    ${checkGACmd}
    """
}


////////////////////////////////
/*    CHECK COMPLEMENTARY     */
////////////////////////////////
/*
Extract RNA regions and check for presence at the opposite side of read in a pair.
*/

// Resulting Channel structure: library, chunk, fastq_table, mapped_oligo1, mapped_oligo2
LIB_TABLE_FASTQ_FOR_RNACOMP
    .combine(LIB_MAPPED_OLIGOS_FOR_RNACOMP, by: [0, 1])
    .branch {
        bridge_forward: it[3]=="bridge_forward"
        ggg: it[3]=="ggg"
    }.combine()
    .filter{
        it[0]==it[6] && it[1]==it[7] && it[4]==1 && it[10]==2
    }.map{
        [it[0], it[1], file(it[2]), file(it[5]), file(it[11]), it[12] ]
    }.combine(LIB_FASTQ_CINDEX_FOR_RNACOMP, by:[0,1])
    .set { LIB_FOR_RNACOMP }

process CHECK_COMPLEMENTARY_RNA_CHUNKS{
    tag "library:${library} chunk:${chunk}"

    storeDir getOutputDir('cout')

    input:
    set val(library), val(chunk), file(table_fq),
        file(cout_br_for), file(cout_ggg_rev),
        file(cindex_fq1), file(cindex_fq2), read_length from LIB_FOR_RNACOMP

    output:
    set library, chunk, "${library}.${chunk}.1.rnacomp.txt", "${library}.${chunk}.2.rnacomp.txt" into LIB_COUT_RNACOMP

    script:

    def rna_complementary_length = params.run.rna_complementary_length
    def rna_comp_length_converted = ((rna_complementary_length+3).intdiv(8))+1
    def extraN = "N"*(500+3*rna_complementary_length)

    """
    # Get the complemetary regions
    paste <(awk '{print \$1, \$3, \$4}' ${table_fq}) \
          <(head -n -1 ${cout_br_for} | tail -n +2 | awk '{print \$4}') \
          | awk 'BEGIN{OFS="\\n";} {bgn=1; if (\$4<500) bgn=\$4+${br_length}+1; \
          print \$1, substr(\$2"${extraN}", bgn, ${rna_complementary_length}), \
          "+", substr(\$3"${extraN}", bgn, ${rna_complementary_length})}' > ${library}.${chunk}.rna-end.1.fq

    paste <(awk '{print \$1, \$5, \$6}' ${table_fq}) \
          <(head -n -1 ${cout_ggg_rev} | tail -n +2 | awk '{print \$5+1}') \
          | awk 'BEGIN{OFS="\\n";} {bgn=1; if (\$4<500) bgn=\$4;
          print \$1, substr(\$2"${extraN}", bgn, ${rna_complementary_length}), \
          "+", substr(\$3"${extraN}", bgn, ${rna_complementary_length})}' > ${library}.${chunk}.rna-end.2.fq

    # Convert to reverse complement
    paste <(sed -n '1~4p' ${library}.${chunk}.rna-end.1.fq) \
          <(sed -n '2~4p' ${library}.${chunk}.rna-end.1.fq | rev | tr "ATGC" "TACG") \
          <(sed -n '3~4p' ${library}.${chunk}.rna-end.1.fq) \
          <(sed -n '4~4p' ${library}.${chunk}.rna-end.1.fq | rev) | tr "\\t" "\\n" \
          > ${library}.${chunk}.rna-end.1.revcomp.fq

    paste <(sed -n '1~4p' ${library}.${chunk}.rna-end.2.fq) \
          <(sed -n '2~4p' ${library}.${chunk}.rna-end.2.fq | rev | tr "ATGC" "TACG") \
          <(sed -n '3~4p' ${library}.${chunk}.rna-end.2.fq) \
          <(sed -n '4~4p' ${library}.${chunk}.rna-end.2.fq | rev) | tr "\\t" "\\n" \
          > ${library}.${chunk}.rna-end.2.revcomp.fq

    # Convert to binary file
    fastq2bin ${library}.${chunk}.rna-end.1.revcomp.fq ${library}.${chunk}.rna-end.1.revcomp.bin
    fastq2bin ${library}.${chunk}.rna-end.2.revcomp.fq ${library}.${chunk}.rna-end.2.revcomp.bin

    # Align complementary regions to each library
    align_pairwise ${library}.${chunk}.rna-end.2.revcomp.bin ${cindex_fq1} 1 ${read_length} \
                   ${seqlen_converted} ${rna_comp_length_converted} 0 \
                   ${read_length-rna_complementary_length} \
                   1 0 ${rna_complementary_length} > ${library}.${chunk}.1.rnacomp.txt

    align_pairwise ${library}.${chunk}.rna-end.1.revcomp.bin ${cindex_fq2} 1 ${read_length} \
               ${seqlen_converted} ${rna_comp_length_converted} 0 \
               ${read_length-rna_complementary_length} \
               1 0 ${rna_complementary_length} > ${library}.${chunk}.2.rnacomp.txt

    """
}

////////////////////////////////
/*       GET SUBSTRINGS       */
////////////////////////////////
/*
Extract DNA, and two RNA regions from FASTQ baed on the mappings.
*/

LIB_MAPPED_OLIGOS_FOR_SUBSTR.into {
    LIB_MAPPED_OLIGOS_FOR_SUBSTR_DNA;
    LIB_MAPPED_OLIGOS_FOR_SUBSTR_RNA1;
    LIB_MAPPED_OLIGOS_FOR_SUBSTR_RNA2;
}
LIB_TABLE_FASTQ_FOR_SUBSTR.into {
    LIB_TABLE_FASTQ_FOR_SUBSTR_DNA;
    LIB_TABLE_FASTQ_FOR_SUBSTR_RNA1;
    LIB_TABLE_FASTQ_FOR_SUBSTR_RNA2
}

LIB_TRIMTABLE.into {
    LIB_TRIMTABLE_FOR_SUBSTR_DNA;
    LIB_TRIMTABLE_FOR_SUBSTR_RNA1;
    LIB_TRIMTABLE_FOR_SUBSTR_RNA2;
    LIB_TRIMTABLE_FOR_COLLECT
}

/* Get DNA substrings */

LIB_MAPPED_OLIGOS_FOR_SUBSTR_DNA
    .branch{
        bridge_forward:  it[2]=="bridge_forward"  && it[3]==1
        adaptor_forward: it[2]=="adaptor_forward" && it[3]==1
    }.set{LIB_MAPPED_BRANCHED_FOR_DNA}

LIB_TABLE_FASTQ_FOR_SUBSTR_DNA
    .combine(LIB_MAPPED_BRANCHED_FOR_DNA.bridge_forward, by:[0,1])  // oligo1
    .combine(LIB_MAPPED_BRANCHED_FOR_DNA.adaptor_forward, by:[0,1]) // oligo2
    .combine(LIB_TRIMTABLE_FOR_SUBSTR_DNA, by:[0,1])
    .map{
        library, chunk, table_fastq,
        oligo1, side1, file_oligo1, read_length1,
        oligo2, side2, file_oligo2, read_length2,
        trim_table
        -> [library, chunk, table_fastq, file_oligo1, file_oligo2, trim_table, read_length1]
    }.set{ LIB_FOR_SUBSTR_DNA }

process GET_DNA_FRAGMENTS_CHUNKS{
    tag "library:${library} chunk:${chunk}"

    storeDir getOutputDir('filtered_fastq')

    input:
    set val(library), val(chunk), file(fastq_table),
        file(cout_r1_br_for), file(cout_r1_for),
        file(trim_table), val(read_length) from LIB_FOR_SUBSTR_DNA

    output:
    set library, chunk, "${library}.${chunk}.dna_nonextended.fq" into LIB_SUBSTR_DNA
    set library, chunk, "${library}.${chunk}.dna.fq" optional true into LIB_SUBSTR_DNA_EXT
    set library, chunk, "${library}.${chunk}.dna.info.txt" into LIB_SUBSTR_DNA_INFO

    script:
    def limit = params.run.min_substring_size

    def extended_Cmd = ""
    def qual_extension = ""
    if (dna_extension.size()>0) {
        qual_extension = "~"*dna_extension.size()
        extended_Cmd = """
        paste <(awk '{print \$1, \$3, \$4}' ${fastq_table}) \
              <(head -n -1 ${cout_r1_for} | tail -n +2 | awk '{print \$5+1}') \
              <(awk '{print \$5}' ${trim_table}) \
              <(head -n -1 ${cout_r1_br_for} | tail -n +2 | awk '{print \$4}') \
              | awk 'BEGIN{OFS="\\n";} {bgn=1; if (\$4<500) bgn=\$4; end=\$5; if (\$6<end) end=\$6; \
              if (end-bgn+1>=${limit}) print \$1, substr(\$2, bgn, end-bgn+1)"${dna_extension}", \
              "+", substr(\$3, bgn, end-bgn+1)"${qual_extension}"}' > ${library}.${chunk}.dna.fq
        """
    }
    """
    ${extended_Cmd}
    paste <(awk '{print \$1, \$3, \$4}' ${fastq_table}) \
          <(head -n -1 ${cout_r1_for} | tail -n +2 | awk '{print \$5+1}') \
          <(awk '{print \$5}' ${trim_table}) \
          <(head -n -1 ${cout_r1_br_for} | tail -n +2 | awk '{print \$4}') \
          | awk 'BEGIN{OFS="\\n";} {bgn=1; if (\$4<500) bgn=\$4; end=\$5; if (\$6<end) end=\$6; \
          if (end-bgn+1>=${limit}) print \$1, substr(\$2, bgn, end-bgn+1), \
          "+", substr(\$3, bgn, end-bgn+1)}' > ${library}.${chunk}.dna_nonextended.fq

    # Write info about selected segments:
    paste <(awk '{print \$1, \$3, \$4}' ${fastq_table}) \
          <(head -n -1 ${cout_r1_for} | tail -n +2 | awk '{print \$5+1}') \
          <(awk '{print \$5}' ${trim_table}) \
          <(head -n -1 ${cout_r1_br_for} | tail -n +2 | awk '{print \$4}') \
          | awk 'BEGIN{OFS="\\t";} {bgn=1; if (\$4<500) bgn=\$4; \
          end=${read_length}; if (\$6<end) end=\$6; \
          print \$1, bgn, end, end-bgn+1, \$5, \$5-bgn+1}' > ${library}.${chunk}.dna.info.txt
    """
}

LIB_COUT_RNACOMP.into { LIB_RNACOMP_TO_SUBSTR_RNA1;
                        LIB_RNACOMP_TO_SUBSTR_RNA2 }

/* Get RNA1 substrings */

LIB_MAPPED_OLIGOS_FOR_SUBSTR_RNA1
    .branch{
        bridge_forward:  it[2]=="bridge_forward"  && it[3]==1
        adaptor_reverse_short: it[2]=="adaptor_reverse_short" && it[3]==1
    }.set{LIB_MAPPED_BRANCHED_FOR_RNA1}

LIB_TABLE_FASTQ_FOR_SUBSTR_RNA1
    .combine(LIB_MAPPED_BRANCHED_FOR_RNA1.bridge_forward, by:[0,1])  // oligo1
    .combine(LIB_MAPPED_BRANCHED_FOR_RNA1.adaptor_reverse_short, by:[0,1]) // oligo2
    .combine(LIB_RNACOMP_TO_SUBSTR_RNA1, by:[0,1])
    .combine(LIB_TRIMTABLE_FOR_SUBSTR_RNA1, by:[0,1])
    .map{
        library, chunk, table_fastq,
        oligo1, side1, file_oligo1, read_length1,
        oligo2, side2, file_oligo2, read_length2,
        rnacomp1, rnacomp2, trim_table
        -> [library, chunk, table_fastq, file_oligo1, file_oligo2, rnacomp1, trim_table, read_length1]
    }.set{ LIB_FOR_SUBSTR_RNA1 }

process GET_RNA1_FRAGMENTS_CHUNKS{
    tag "library:${library} chunk:${chunk}"

    storeDir getOutputDir('filtered_fastq')

    input:
    set val(library), val(chunk), file(fastq_table),
        file(cout_r1_br_for), file(cout_r1_rev),
        file(cout_compl_1), file(trim_table), val(read_length) from LIB_FOR_SUBSTR_RNA1

    output:
    set library, chunk, "${library}.${chunk}.rna1.fq" into LIB_SUBSTR_RNA1
    set library, chunk, "${library}.${chunk}.rna1.info.txt" into LIB_SUBSTR_RNA1_INFO

    script:
    def limit = params.run.min_substring_size

    """
    paste <(awk '{print \$1, \$3, \$4}' ${fastq_table}) <(awk '{print \$5}' ${trim_table}) \
        <(head -n -1 ${cout_r1_br_for} | tail -n +2 | awk '{print \$4}') \
        <(head -n -1 ${cout_r1_rev} | tail -n +2 | awk '{print \$4}') \
        <(head -n -1 ${cout_compl_1} | tail -n +2 | awk '{print \$5}') \
        | awk 'BEGIN{OFS="\\n";} {bgn=\$4; if (\$5<500) bgn=\$5+${br_length}+1; \
        end=\$4; if (\$6<end) end=\$6; if (\$7<end) end=\$7; \
        if (end-bgn+1>=${limit}) print \$1, substr(\$2, bgn, end-bgn+1), \
        "+", substr(\$3, bgn, end-bgn+1)}' > ${library}.${chunk}.rna1.fq

    paste <(awk '{print \$1, \$3, \$4}' ${fastq_table}) \
        <(awk '{print \$5}' ${trim_table}) \
        <(head -n -1 ${cout_r1_br_for} | tail -n +2 | awk '{print \$4}') \
        <(head -n -1 ${cout_r1_rev} | tail -n +2 | awk '{print \$4}') \
        <(head -n -1 ${cout_compl_1} | tail -n +2 | awk '{print \$5}') \
        | awk 'BEGIN{OFS="\\t";} {bgn=1; if (\$5<500) bgn=\$5+${br_length}+1; \
        end=${read_length}; if (\$6<end) end=\$6; if (\$7<end) end=\$7; \
        print \$1, bgn, end, end-bgn+1, \$4, \$4-bgn+1}' > ${library}.${chunk}.rna1.info.txt
    """
}


/* Get RNA2 substrings */

LIB_MAPPED_OLIGOS_FOR_SUBSTR_RNA2
    .branch{
        ggg:  it[2]=="ggg"  && it[3]==2
        adaptor_forward: it[2]=="adaptor_forward" && it[3]==2
        bridge_reverse:  it[2]=="bridge_reverse"  && it[3]==2
    }.set{LIB_MAPPED_BRANCHED_FOR_RNA2}

LIB_TABLE_FASTQ_FOR_SUBSTR_RNA2
    .combine(LIB_MAPPED_BRANCHED_FOR_RNA2.ggg, by:[0,1])             // oligo1
    .combine(LIB_MAPPED_BRANCHED_FOR_RNA2.adaptor_forward, by:[0,1]) // oligo2
    .combine(LIB_MAPPED_BRANCHED_FOR_RNA2.bridge_reverse, by:[0,1])  // oligo3
    .combine(LIB_RNACOMP_TO_SUBSTR_RNA2, by:[0,1])
    .combine(LIB_TRIMTABLE_FOR_SUBSTR_RNA2, by:[0,1])
    .map{
        library, chunk, table_fastq,
        oligo1, side1, file_oligo1, read_length1,
        oligo2, side2, file_oligo2, read_length2,
        oligo3, side3, file_oligo3, read_length3,
        rnacomp1, rnacomp2, trim_table
        -> [library, chunk, table_fastq, file_oligo1, file_oligo2, file_oligo3, rnacomp2, trim_table, read_length1]
    }.set{ LIB_FOR_SUBSTR_RNA2 }

process GET_RNA2_FRAGMENTS_CHUNKS{
    tag "library:${library} chunk:${chunk}"

    storeDir getOutputDir('filtered_fastq')

    input:
    set val(library), val(chunk), file(fastq_table),
        file(cout_r2_ggg), file(cout_r2_for), file(cout_r2_br_rev),
        file(cout_compl_2), file(trim_table), val(read_length) from LIB_FOR_SUBSTR_RNA2

    output:
    set library, chunk, "${library}.${chunk}.rna2.fq" into LIB_SUBSTR_RNA2
    set library, chunk, "${library}.${chunk}.rna2.info.txt" into LIB_SUBSTR_RNA2_INFO

    script:
    def limit = params.run.min_substring_size

    """
    paste <(awk '{print \$1, \$5, \$6}' ${fastq_table}) \
        <(head -n -1 ${cout_r2_ggg} | tail -n +2 | awk '{print \$5+1}') \
        <(head -n -1 ${cout_r2_for} | tail -n +2 | awk '{print \$5+1}') \
        <(awk '{print \$7}' ${trim_table}) \
        <(head -n -1 ${cout_r2_br_rev} | tail -n +2 | awk '{print \$4}') \
        <(head -n -1 ${cout_compl_2} | tail -n +2 | awk '{print \$5}') \
        | awk 'BEGIN{OFS="\\n";} {bgn=1; if (\$4<500) bgn=\$4; if (\$5<500 && \$5>bgn) bgn=\$5; \
        end=\$6; if (\$7<end) end=\$7; if (\$8<end) end=\$8; \
        if (end-bgn+1>=${limit}) print \$1, substr(\$2, bgn, end-bgn+1), \
        "+", substr(\$3, bgn, end-bgn+1)}' > ${library}.${chunk}.rna2.fq

    paste <(awk '{print \$1, \$5, \$6}' ${fastq_table}) \
        <(head -n -1 ${cout_r2_ggg} | tail -n +2 | awk '{print \$5+1}') \
        <(head -n -1 ${cout_r2_for} | tail -n +2 | awk '{print \$5+1}') \
        <(awk '{print \$7}' ${trim_table}) \
        <(head -n -1 ${cout_r2_br_rev} | tail -n +2 | awk '{print \$4}') \
        <(head -n -1 ${cout_compl_2} | tail -n +2 | awk '{print \$5}') \
        | awk 'BEGIN{OFS="\\t";} {bgn=1; if (\$4<500) bgn=\$4; if (\$5<500 && \$5>bgn) bgn=\$5; \
        end=${read_length}; if (\$7<end) end=\$7; if (\$8<end) end=\$8; \
        print \$1, bgn, end, end-bgn+1, \$6, \$6-bgn+1}' > ${library}.${chunk}.rna2.info.txt
    """
}

////////////////////////////////
/*           MAPPING          */
////////////////////////////////
/*
Map DNA, RNA1 and RNA2 parts
*/

GENOME_INDEX.into{
    GENOME_INDEX_FOR_DNA;
    GENOME_INDEX_FOR_DNA_EXT;
    GENOME_INDEX_FOR_RNA1;
    GENOME_INDEX_FOR_RNA2
}

/* Map DNA */
LIB_SUBSTR_DNA.combine(GENOME_INDEX_FOR_DNA)
    .set{LIB_FOR_DNA_MAPPING}

process MAP_DNA_NONEXTENDED_CHUNKS{
    tag "library:${library} chunk:${chunk}"

    storeDir getOutputDir('sam')

    input:
    set val(library), val(chunk), file(input_dna), val(index_pref), file(genome_index) from LIB_FOR_DNA_MAPPING

    output:
    set library, chunk, "${library}.${chunk}.dna_nonextended.sam" into LIB_SAM_DNA

    script:
    """
    hisat2 -p ${task.cpus} -x ${index_pref} --no-spliced-alignment -k 100 \
           --no-softclip -U ${input_dna} > ${library}.${chunk}.dna_nonextended.sam
    """
}

if (dna_extension.size()>0) {
    LIB_SUBSTR_DNA_EXT.combine(GENOME_INDEX_FOR_DNA_EXT)
        .set{LIB_FOR_DNA_MAPPING_EXT}

    process MAP_DNA_EXTENDED_CHUNKS{
        tag "library:${library} chunk:${chunk}"

        storeDir getOutputDir('sam')

        input:
        set val(library), val(chunk), file(input_dna), val(index_pref), file(genome_index) from LIB_FOR_DNA_MAPPING_EXT

        output:
        set library, chunk, "${library}.${chunk}.dna.sam" into LIB_SAM_DNA_EXT

        script:
        """
        hisat2 -p ${task.cpus} -x ${index_pref} --no-spliced-alignment -k 100 \
               --no-softclip -U ${input_dna} > ${library}.${chunk}.dna.sam
        """
    }
}

GENOME_SPLICESITES.into{GENOME_SPLICESITES_RNA1; GENOME_SPLICESITES_RNA2}

/* Map RNA1 */
LIB_SUBSTR_RNA1.combine(GENOME_INDEX_FOR_RNA1).combine(GENOME_SPLICESITES_RNA1)
    .set{LIB_FOR_RNA1_MAPPING}

process MAP_RNA1_CHUNKS{
    tag "library:${library} chunk:${chunk}"

    storeDir getOutputDir('sam')

    input:
    set val(library), val(chunk), file(input_rna1), val(index_pref), file(genome_index), file(known_splicesites) from LIB_FOR_RNA1_MAPPING

    output:
    set library, chunk, "${library}.${chunk}.rna1.sam" into LIB_SAM_RNA1

    script:
    """
    hisat2 -p ${task.cpus} -x ${index_pref} -k 100 --no-softclip --known-splicesite-infile ${known_splicesites} \
        --dta-cufflinks --novel-splicesite-outfile ${library}.${chunk}.novel.splicesites.txt \
        -U ${input_rna1} > ${library}.${chunk}.rna1.sam
    """
}

/* Map RNA2 */
LIB_SUBSTR_RNA2.combine(GENOME_INDEX_FOR_RNA2).combine(GENOME_SPLICESITES_RNA2)
    .set{LIB_FOR_RNA2_MAPPING}

process MAP_RNA2_CHUNKS{
    tag "library:${library} chunk:${chunk}"

    storeDir getOutputDir('sam')

    input:
    set val(library), val(chunk), file(input_rna2), val(index_pref), file(genome_index), file(known_splicesites) from LIB_FOR_RNA2_MAPPING

    output:
    set library, chunk, "${library}.${chunk}.rna2.sam" into LIB_SAM_RNA2

    script:
    """
    hisat2 -p ${task.cpus} -x ${index_pref} -k 100 --no-softclip --known-splicesite-infile ${known_splicesites} \
        --dta-cufflinks --novel-splicesite-outfile ${library}.${chunk}.novel.splicesites.txt \
        -U ${input_rna2} > ${library}.${chunk}.rna2.sam
    """
}


////////////////////////////////
/*       PARSE SAM TO BED     */
////////////////////////////////
/*
Parse SAM files to BED files
*/

LIB_SAM_RNA1.into{LIB_SAM_RNA1_FOR_BED; LIB_SAM_RNA1_FOR_COLLECT}
LIB_SAM_RNA2.into{LIB_SAM_RNA2_FOR_BED; LIB_SAM_RNA2_FOR_COLLECT}

if (dna_extension.size()>0){
    LIB_SAM_DNA.into{LIB_SAM_DNA_FOR_BED; LIB_SAM_DNA_FOR_COLLECT}
    LIB_SAM_DNA_EXT.into{LIB_SAM_DNA_EXT_FOR_BED; LIB_SAM_DNA_EXT_FOR_COLLECT}
    LIB_SAM_DNA_FOR_BED.map{ [it[0], it[1], "dna_nonextended", file(it[2]) ] }
        .concat(LIB_SAM_DNA_EXT_FOR_BED.map{ [it[0], it[1], "dna", file(it[2]) ] })
        .concat(LIB_SAM_RNA1_FOR_BED.map{ [it[0], it[1], "rna1", file(it[2]) ] })
        .concat(LIB_SAM_RNA2_FOR_BED.map{ [it[0], it[1], "rna2", file(it[2]) ] }).set{LIB_SAM2BED}
} else {
    LIB_SAM_DNA.set{LIB_SAM_DNA_EXT_FOR_COLLECT}
    LIB_SAM_DNA_FOR_BED.map{ [it[0], it[1], "dna_nonextended", file(it[2]) ] }
        .concat(LIB_SAM_RNA1_FOR_BED.map{ [it[0], it[1], "rna1", file(it[2]) ] })
        .concat(LIB_SAM_RNA2_FOR_BED.map{ [it[0], it[1], "rna2", file(it[2]) ] }).set{LIB_SAM2BED}
}

process SAM2BED_CHUNKS{
    tag "library:${library} chunk:${chunk} ${segment_name}"

    storeDir getOutputDir('bed')

    input:
    set val(library), val(chunk), val(segment_name), file(sam) from LIB_SAM2BED

    output:
    set library, chunk, segment_name, "${library}.${chunk}.${segment_name}.bed" into LIB_BED

    script:
    // 2 mismatches, up to 1 alignment
    """
    samtools view -Sh -F 4 ${sam} | grep -E 'XM:i:[0-2]\\s.*NH:i:1\$|^@' | samtools view -Sbh - \
        | bedtools bamtobed -cigar -i stdin > ${library}.${chunk}.${segment_name}.bed
    """
}

LIB_BED.into{LIB_BED_FOR_RESTR; LIB_BED_FOR_COLLECT}

////////////////////////////////
/*     CHECK RESTRICTION      */
////////////////////////////////
/*
Check presence of restriction enzyme recognition sites at the ends of RNA parts.
*/

if (!check_restriction) {
    Channel.create().into{LIB_DISTANCES; LIB_RESTR_COMBINATIONS}
} else {
    available_renz = params.protocol.renz

    // Channel structure: key of the renz, segment (rna1, rna2 or dna) and strand of the renz (+ or -)
    Channel.fromList(
        params.run.restriction_check.collect{ segment, renz_keys ->
            renz_keys.collect{ renz_key ->
                available_renz.containsKey(renz_key) ?
                    [renz_key, segment, ''] :
                    renz_key.endsWith('p') ?
                        [renz_key[0..-2], segment, '+'] :
                        [renz_key[0..-2], segment, '-']
            }
        }.sum()
    ).into{LIST_FOR_RESTR_RUN; LIB_RESTR_COMBINATIONS}

    LIST_FOR_RESTR_RUN
        .combine(GENOME_RENZ, by:0)
        .combine(LIB_BED_FOR_RESTR)
        .filter{
            renz_key, segment_left, renz_strand, renz_file, library, chunk, segment_right, bed_file ->
            segment_left == segment_right
         }.map{
            renz_key, segment_left, renz_strand, renz_file, library, chunk, segment_right, bed_file ->
            [library, chunk, segment_left, file(bed_file), renz_key, renz_strand, file(renz_file)]
         }.set{ LIB_FOR_RESTR_RUN }

    process ANNOTATE_RENZYMES_CHUNKS{
        tag "library:${library} ${chunk} ${segment_name} ${renz_key}${renz_strand}"

        storeDir getOutputDir('table')

        input:
        set library, chunk, segment_name, file(bed_file), renz_key, renz_strand, file(renz_file) from LIB_FOR_RESTR_RUN

        output:
        set library, chunk, segment_name, renz_key, renz_strand,
            "${library}.${chunk}.${segment_name}.${renz_key}${renz_strand}.distances.txt" into LIB_DISTANCES

        script:
        def renz_strand_key = (renz_strand=="+") ? "p" : (renz_strand=="-") ? "n" : ""
        def renz_strand_sub = (renz_strand=="+") ? "+" : (renz_strand=="-") ? "-" : "+"
        def columns = [
            "${segment_name}_start_${renz_key}${renz_strand_key}_left",
            "${segment_name}_start_${renz_key}${renz_strand_key}_right",
            "${segment_name}_end_${renz_key}${renz_strand_key}_left",
            "${segment_name}_end_${renz_key}${renz_strand_key}_right"
            ]
        def header = (["id"]+columns).join(" ")
        """
        echo "${header}" > ${library}.${chunk}.${segment_name}.${renz_key}${renz_strand}.distances.txt
        get_closest_sites.py ${bed_file} ${renz_file} ${renz_strand_sub} ${library}.${chunk}.${segment_name}.${renz_key}${renz_strand}.distances.txt
        """
    }
}

///////////////////////////////////////////
/*      COLLECT DATASETS AND FILTERS     */
///////////////////////////////////////////
/*
Collect all outputs into two hdf5 files: collected data and filters
*/

LIB_MAPPED_OLIGOS_FOR_COLLECT.branch{
            bridge_forward: it[2]=="bridge_forward" & it[3]==1
            ggg: it[2]=="ggg" & it[3]==2
        }.set{LIB_MAPPED_OLIGOS_FOR_COLLECT_BRANCHED}

dna_mode = (dna_extension.size()>0) ? "dna" : "dna_nonextended"
LIB_BED_FOR_COLLECT.branch{
        dna : it[2]==dna_mode
        rna1 : it[2]=="rna1"
        rna2 : it[2]=="rna2"
    }.set{LIB_BED_FOR_COLLECT_BRANCHED}

LIB_DISTANCES
    .groupTuple(by: [0,1]).map{it -> [it[0], it[1], it[5]]}.set{ LIB_DISTANCES_FOR_COLLECT }

IDS_FASTUNIQ
    .combine(LIB_TABLE_FASTQ_FOR_COLLECT, by:0).map{ v -> [v[0], v[2], v[3], v[1]]}
    .combine(LIB_TRIMTABLE_FOR_COLLECT, by: [0,1])
    .combine(LIB_MAPPED_OLIGOS_FOR_COLLECT_BRANCHED.bridge_forward.map{it -> [it[0], it[1], it[4]]}, by: [0,1])
    .combine(LIB_MAPPED_GA, by: [0,1])
    .combine(LIB_MAPPED_OLIGOS_FOR_COLLECT_BRANCHED.ggg.map{it -> [it[0], it[1], it[4]]}, by: [0,1])
    .combine(LIB_SUBSTR_DNA_INFO, by:[0,1])
    .combine(LIB_SUBSTR_RNA1_INFO, by:[0,1])
    .combine(LIB_SUBSTR_RNA2_INFO, by:[0,1])
    .combine(LIB_SAM_DNA_EXT_FOR_COLLECT, by:[0,1])
    .combine(LIB_SAM_DNA_FOR_COLLECT, by:[0,1])
    .combine(LIB_SAM_RNA1_FOR_COLLECT, by:[0,1])
    .combine(LIB_SAM_RNA2_FOR_COLLECT, by:[0,1])
    .combine(LIB_BED_FOR_COLLECT_BRANCHED.dna.map{it -> [it[0], it[1], it[3]]}, by:[0,1])
    .combine(LIB_BED_FOR_COLLECT_BRANCHED.rna1.map{it -> [it[0], it[1], it[3]]}, by:[0,1])
    .combine(LIB_BED_FOR_COLLECT_BRANCHED.rna2.map{it -> [it[0], it[1], it[3]]}, by:[0,1])
    .combine(LIB_DISTANCES_FOR_COLLECT, by:[0,1])
    .map{ v -> [ v[0], v[1], v[2..-2]+v[-1] ] }.set{ LIB_COLLECT }

process COLLECT_DATA_CHUNKS{
        tag "library:${library} ${chunk}"

        storeDir getOutputDir('hdf5')

        input:
        set val(library), val(chunk), file(input) from LIB_COLLECT

        output:
        set library, chunk, "${library}.${chunk}.data.hdf5" into LIB_COLLECTED

        script:
        """
        collect_data.py ${library}.${chunk}.data.hdf5 ${input}
        """
}

def restriction_patterns = params.filters.restriction
                   .collect{k, v -> k+":"+"("+v.join(") | (")+")"}.join("\\n")
                   .replaceAll(/"\+"/, "1")
                   .replaceAll(/"-"/, "0")
def additional_patterns = params.filters.additional_filters
                   .collect{k, v -> k+":"+v}.join("\\n")
patterns = [restriction_patterns, additional_patterns].join("\\n")

LIB_COLLECTED.combine(GENOME_CHROMSIZES).set{LIB_COLLECTED_FOR_FILTERS}

process COLLECT_FILTERS_CHUNKS{
        tag "library:${library} ${chunk}"

        storeDir getOutputDir('hdf5')

        input:
        set val(library), val(chunk), file(input), file(genome_chromsizes) from LIB_COLLECTED_FOR_FILTERS

        output:
        set library, chunk, "${library}.${chunk}.data.hdf5", "${library}.${chunk}.filters.hdf5" into LIB_FILTERS

        script:
        def chrom_pattern = params.filters.canonical_chromosomes

        """
        printf '${patterns}' > filters.txt
        filter_data.py ${input} ${library}.${chunk}.filters.hdf5 \
            ${genome_chromsizes} "${chrom_pattern}" filters.txt
        """
}

LIB_FILTERS.into{ LIB_FILTERS_STATS; LIB_FILTERS_TABLE }

//////////////////////////////////
/*         WRITE STATS          */
//////////////////////////////////
/*
Write params.report_stats from filters.hdf5 into text file.
Then we merge stats on chunks into a single file.
*/

def stats_list = params.report_stats.collect()

process WRITE_STATS_CHUNKS{
        tag "library:${library} chunk:${chunk}"

        storeDir getOutputDir('stats')

        input:
        set val(library), val(chunk), file(data_hdf5), file(filters_hdf5) from LIB_FILTERS_STATS

        output:
        set library, chunk, "${library}.${chunk}.stats.txt" into FILES_STATS

        script:
        def stats_str = '["'+ stats_list.join('", "')+'"]'
        """
        #!/usr/bin/env python

import h5py
f = h5py.File("${filters_hdf5}", 'r')
with open("${library}.${chunk}.stats.txt", "w") as outfile:
    for filt in ${stats_str}:
        n = int(sum(f[filt][()]))
        outfile.write(f"{filt}\\t{n}\\n")
f.close()
        """
}

FILES_STATS.groupTuple(by:0).set{FILES_STATS_FOR_MERGE}

process MERGE_STATS{
        tag "library:${library}"

        storeDir getOutputDir('stats')

        input:
        set val(library), val(chunk), file(files_stats) from FILES_STATS_FOR_MERGE

        output:
        set library, "${library}.stats.txt" into FILES_STATS_MERGED

        script:
        def files_str = '["'+ files_stats.join('", "')+'"]'
        """
        #!/usr/bin/env python

import pandas as pd
res = []
for f in ${files_str}:
    tmp = pd.read_csv(f, sep='\\t', header=None).set_index(0)
    if len(res)==0:
        res = tmp.copy()
    else:
        res += tmp
res.to_csv("${library}.stats.txt", sep='\t', header=False, index=True)
        """
}

//////////////////////////////////
/*       WRITE FINAL TABLES     */
//////////////////////////////////
/*
Write final tables from filters.hdf5 and data.hdf5 into text file.
Then we merge chunks into a single file.
You can request one or several files with different information in params.final_table.
*/

if (params.final_table.create_final_table){
    Channel.fromList(
        tables_list = params.final_table.tables.collect{ k, v -> [k, v['filter'], v['header']]}
        ).set{LIST_TABLES}

    LIB_FILTERS_TABLE.combine(LIST_TABLES).set{LIB_FOR_WRITING}

    process WRITE_FINAL_TABLE_CHUNKS{
            tag "library:${library} chunk:${chunk} table:${table_name}"

            storeDir getOutputDir('final_table')

            input:
            set val(library), val(chunk), file(data_hdf5), file(filters_hdf5),
                val(table_name), val(filter), val(header) from LIB_FOR_WRITING

            output:
            set library, chunk, table_name, "${library}.${chunk}.${table_name}.tsv" into FILES_TABLE

            script:
            def stats_str = '["'+ stats_list.join('", "')+'"]'
            """
            #!/usr/bin/env python

    import numpy as np

    import h5py
    f_filt = h5py.File("${filters_hdf5}", 'r')
    f_data = h5py.File("${data_hdf5}", 'r')

    header = "${header}".split()
    res = {}
    for col in header+["${filter}"]:
        try:
            res[col] = f_filt[col][()]
        except Exception as e:
            res[col] = f_data[col][()]
    f_filt.close()
    f_data.close()

    def strand(x):
        if x:
            return "+"
        return "-"

    indexes = np.where(res["${filter}"])[0]
    with open("${library}.${chunk}.${table_name}.tsv", "w") as outfile:
        outfile.write( "\\t".join(header)+"\\n" )
        for i in indexes:
            line = [res[col][i] if isinstance(res[col][i], str) else
                    strand(res[col][i]) if 'strand' in col else
                    res[col][i].decode() if isinstance(res[col][i], np.bytes_) else
                    str(res[col][i]) for col in header]
            line = "\\t".join( line )+"\\n"
            outfile.write( line )

            """
    }

    FILES_TABLE.groupTuple(by:[0,2]).set{FILES_TABLE_FOR_MERGE}

    process MERGE_TABLE{
            tag "library:${library} table:${table_name}"

            storeDir getOutputDir('final_table')

            input:
            set val(library), val(chunk), val(table_name), file(files_stats) from FILES_TABLE_FOR_MERGE

            output:
            set library, "${library}.${table_name}.tsv" into FILES_TABLE_MERGED

            script:
            """
            head -n 1 "${files_stats[0]}" > ${library}.${table_name}.tsv
            for FILE in ${files_stats}
            do
                tail -n +2 \$FILE >> ${library}.${table_name}.tsv
            done
            """
    }
}
