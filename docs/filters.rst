Filters of redc-nf
==================

Filters are named criteria that can be used for selection of
the subset of reads or reporting the statistics.
you can specify your filters in `project.yml` file.

There are three types of filters in redc-nf:
- filter for canonical chromosomes (the field called `filters.canonical_chromosomes`)
- restriction filters (`filters.restriction`)
- additional filters (`filters.additional_filters`)

**Filter for canonical chromosomes** is a regular expression to retain only chromosomes of interest.
It should correspond to `python re module syntax <https://docs.python.org/3/library/re.html>`_ for regular expressions.

**Other redc-nf filters** are strings with Python syntax. These strings are  *expressions* that operate with *variables*. Each expression will be evaluated for each read in the dataset. The variables are the columns of redc-nf output hdf5 file.

Functions in expressions of filters
--------------------
You can use any builtin functions of Python or numpy functions (assuming numpy is loaded as `np`).

Variables in expressions of filters
--------------------
**Sources of variables**
As a result of redc-nf run, the data on oligos alignment, RNA and DNA mapping, etc. are collected in a single
hdf5 file, where each column is stored as numpy vector.
You can use any of these columns as variables in the filter's expression.
Note that some of the columns have string data type, some are integers and some are boolean.

Restriction of the DNA/RNA segments is a separated step.
You specify the list of restirction enzymes in `protocol.renz` with short names (or keys) and actual names of these restiriction enzymes.
Then, redc-nf uses these enzyme names for detection of restriction enzymes recognition sites around RNA and DNA segments (with BioPython).
Not all the combinations of restriction enzymes and segments are meaningful, thus you specify only the annotations that you want
in `run.restriction_check`.

Also, if you specify some filter in the beginning of the list, its name can be used in downstream filters.

**List of variables loaded by default, without restriction**

    - `id` numpy S100: read ID
    - `seqR1`, `seqR2` numpy S250: forward and reverse read sequences
    - `seqR1_len`, `seqR2_len` numpy int: length of forward and reverse read
    - `is_notPCRdup` bool: flag if read is not PCR duplication
    - `trimF`, `trimR` int: position of trimming on forward and reverse read
    - `has_nobridge` bool: flag if read has no bridge in it
    - `has_GA` bool: flag if bridge has appropriate GA in it
    - `bridge_nmm` numpy int: number or mismatches in bridge in forward read
    - `bridge_start`,  `bridge_end` numpy int: start and end position of bridge in forward read
    - `has_noggg` bool: flag if reverse read does not start with GGG
    - `ggg_start`, `ggg_end` numpy int: start and end position of GGG in reverse read
    - `dna_R1_start`, `dna_R1_end`, `dna_R1_end_trim`, `dna_R1_len_notrim`, `dna_R1_len_trim` numpy int: 
DNA segment start, end in forward read and length, before and after trimming
    - `rna1_R1_start`, `rna1_R1_end`, `rna1_R1_end_trim`, `rna1_R1_len_notrim`, `rna1_R1_len_trim` numpy int:
rna1 segment start, end in forward read and length, before and after trimming
    - `rna2_R1_start`, `rna2_R1_end`, `rna2_R1_end_trim`, `rna2_R1_len_notrim`, `rna2_R1_len_trim` numpy int:
rna2 segment start, end in forward read and length, before and after trimming
    - `dna_is_mapped`, `dna_nonextended_is_mapped`, `rna1_is_mapped`, `rna2_is_mapped` bool: 
flag if DNA with extension / DNA / RNA1 / RNA2 mapped (to all possible chromosomes in index file)
    - `dna_is_not_multi`, `dna_nonextended_is_not_multi`, `rna1_is_not_multi`, `rna2_is_not_multi` bool: 
flag if DNA with extension / DNA / RNA1 / RNA2 are not multimappers
    - `dna_chr` numpy S8; `dna_start`, `dna_end` int; `dna_strand` bool; `dna_cigar` numpy S10 : 
information on mapping of DNA with extension
    - `dna_nonextended_chr` numpy S8; `dna_nonextended_start`, `dna_nonextended_end` int; `dna_nonextended_strand` bool;
`dna_nonextended_cigar` numpy S10 : information on mapping of DNA without extension
    - `rna1_chr` numpy S8; `rna1_start`, `rna1_end` int; `rna1_strand` bool; `rna1_cigar` numpy S10 : 
information on mapping of RNA1 segment
    - `rna2_chr` numpy S8; `rna2_start`, `rna2_end` int; `rna2_strand` bool; `rna2_cigar` numpy S10 : 
information on mapping of RNA2 segment
    - `dna_chr_canonical`, `rna1_chr_canonical`, `rna2_chr_canonical` - flags if segments are mapped to canonical chromosomes

**List of variables loaded after annotation of restriction**

1. If the restriction enzyme with name nla appeared for rna1 fragment in `run.restriction_check`, and it corresponds
to enzyme with palindromic recognition site (e.g. NlaIII), the following columns will be added (integer type):

.. code-block:: python

    rna1_end_nla_left, rna1_end_nla_right, rna1_start_nla_left, rna1_start_nla_right

2. If the restriction enzyme is not palindromic, then the DNA strand matters, and we hit + (p) and - (n) strands separately.
These columns will be added:

.. code-block:: python

    rna1_start_mmep_left, rna1_start_mmep_right, rna1_start_mmen_left, rna1_start_mmen_right,
    rna1_end_mmep_left, rna1_end_mmep_right, rna1_end_mmen_left, rna1_end_mmen_right

Examples of filters
------------------
In project.yml file we provide the filters used in `original paper on RedC<https://doi.org/10.1093/nar/gkaa457>`_.

If you want to run custom output of redc-nf, you can design your own filters.

**Simple indicator filters**
Read is mapped to the canonical chromosomes:

.. code-block:: Python
    dna_chr_canonical & rna1_chr_canonical & rna2_chr_canonical

Restriction filters passed successfully (thus all individual restriction filters have not failed):

.. code-block:: Python
    ~rna1_nla_failed & ~rna2_nla_failed & ~rna2_mme_failed

**Complex conitions**
Length of DNA segment after trimming is between 14 and 21 basepairs:

.. code-block:: Python
    (dna_R1_len_notrim>=18)&(dna_R1_len_notrim<=20)

Distance between RNA1 and RNA2 segments mapping positions is small enough and can be considered a single molecule:

.. code-block:: Python
    np.abs(rna2_start-rna1_start)<1e5

Reporting the filters
--------------------
You can use filters in the file with final statistics. Specify the filters of interest in `report_stats` field of
`project.yml`.

Also, you can report the outcome of evaluation of each filter for each read.
For that, specify the filter name in the header of table in
`final_table.tables field`.
