Processes of redc-nf
====================

download_genome
---------------
For genome download, redc-nf uses the field ``genome`` of ``project.yml`` file.

- ``genome.assembly_name``, "genome" by default. Sets the name of genome files.
- ``genome.auto_download``, true by default. Sets whether the genome will be downloaded automatically from UCSC.
This will also force the creation of new index for hisat2 (time-consuming step).
- ``genome.fasta``, ``genome.chromsizes`` and ``genome.index_prefix`` can be optionally provided as the paths if ``genome.auto_download`` is false.
Note that ``genome.fasta`` and ``genome.chromsizes`` can be ftp/http/https links instead.
New index will be built if no value is provided or fasta file is a remote file.

restrict_genome
---------------
redc-nf will annotate the restriction recognition sites for all the enzymes provided in ``protocol.renz`` field.
Throughout the pipeline, the keys of renz will be used instead of full names. These keys will be used for the creation of restriction filters.

redc-nf uses BioPython for the recognition of restriction sites, annotates for forward and reverse strands of the genome and
provides the tart of the recognition site as well as its strand. End of restriction site is formal (same as start).
See ``bin/detect_restriction_sites.py`` for more details.

prepare_rna_annot
-----------------
redc-nf requires the annotation of splicing sites for appropriate run of hisat2 for RNA segments.
The RNA annotation can be specified as ``rna_annotation.genes_gtf`` with ``rna_annotation.rna_annotation_name``.
``rna_annotation.genes_gtf`` can be local or remote file.
Note that if file is compressed, its extension should be .gz.

split_fastq
-----------
Input .fastq file will be split into chunks of the size provided in ``run.chunksize`` field.
Each chunk will have its index (starting from 001 and so on) and will be processed independently.
Each chunk will be stored as compressed fastq file in the folder named as specified in ``output.dirs.fastq``.

prepare_fastq
-------------
redc-nf requires information of each read sequence and id in "one line-one entry" format. This step will create table with reads
in the folder names as specified in ``output.dirs.table``

run_fastuniq
------------
Detection of PCR duplicates. Note that his step is done for the complete fastq input file, not the chunks.
Prior to fastuniq, the reads will be cropped to the equal length (5'-end retained), as specified in ``run.fastuniq_crop``.

run_trimmomatic and get_trim_output
-----------------------------------
Trimming of the input fastq files with trimmomatic. Parameters of trimmomatic can be specified in ``run.params_trimmomatic``.

index_oligos and index_fastq
------------
Custom alignment procedure of redc-nf requires indexing of fasta files with input oligos and fastq files with reads.
Very important parameters for these processes are: ``run.read_length``, ``run.bridge_length`` and all the ``input`` field.

map_oligos
----------
Custom alignment of oligos to the reads with C implementation of Karp–Rabin algorithm.
If you see poor alignment of oligos, check whether ``run.read_length`` field corresponds to real read lengths in the file.

check_GA
--------
Additional step of search of GA in the bridge. This check can be switched off by setting ``run.check_GA`` to false.

check_rna_complementary
----------------------
Forward and reverse reads can contain the remnants of the bridge/adaptors at the ends, which are not detected by the algorithm
(because they are too short). Thus we search for complementary regions of RNAs in forward and reverse and
select the positions of RNA1 and RNA2 segments more precise.

You can vary the length of complementary regions by ``run.rna_complementary_length``.


get_dna_substrings, get_rna1_substrings, get_rna2_substrings
------------------------------------------------------------
Cuts DNA, RNA1 and RNA2 segments out of initial reads, and creates fastq files for them.

If you set ``run.dna_extension``, then the extension will be added to DNA. By default, no extension is added.

If the segment length is shorter than ``run.min_substring_size`` (in basepairs), then it will be omitted.
14 bp is selected for human genome as the size when expect at least two mappings in the genome of this size.

map_dna_nonextended, map_dna_extended, map_rna1, map_rna2
---------------------------------------------------------
Mapping of DNA, RNA1 and RNA2 segments with hisat2.

sam2bed
-------
Conversion of sam files into bed.

annotate_renzymes
----------------
Annotation of recognition sites of the restriction enzymes.

collect_data
------------
Collect all the relevant data into a single hdf5 file (for each chunk independently).

collect_filters
---------------
Calculate all the filters that are present in ``filters`` fields
and write hdf5 file with filters (for each chunk independently).

write_stats and merge_stats
-----------
Calculate stats specified in ``report_stats`` field and report for each chunk and merged library.

write_table and merge_table
-----------
Write tables as specified in ``final_table.tables``, for each chunk and merged library.
You may switch off this step by setting ``final_table.create_final_table`` to false.
