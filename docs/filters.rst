Filters of redc-nf
==================

Filters are named criteria that can be used for selection of
the subset of reads or reporting the statistics.
you can specify your filters in `project.yml` file.

There are three types of filters used in redc-nf:
- filter for canonical chromosomes (`filters.canonical_chromosomes`)
- restriction filters (the field called `filters.restriction`)
- additional filters (`filters.additional_filters`)

**Filter for canonical chromosomes** is a regular expression to retain only chromosomes of interest.
It should correspond to python re module syntax for regular expressions.

**Other redc-nf filters** are strings with Python syntax that will be evaluated for each read in the dataset.
You can use a broad list of variables in this syntax.

Functions in filters
--------------------
You can use any builtin functions of Python or numpy functions (assuming numpy is loaded as np).

Variables in filters
--------------------
**Source of variables**
As a result of redc-nf, the data on mapping of oligos, RNA and DNA parts are collected in a single
hdf5 file, where each column is stored as numpy vector of a particular data type.
You can use any of these columns as input variables for evaluation of the filter.

Full list of variables:

.. code-block:: python

    id, is_notPCRdup, seqR1, seqR1_len, seqR2, seqR2_len, trimF, trimR,
    has_GA, has_nobridge, bridge_end, bridge_nmm, bridge_start, has_noggg, ggg_end, ggg_start,
    dna_R1_len_notrim, dna_R1_len_trim, dna_R1_start, dna_R1_end, dna_R1_end_trim,
    dna_chr, dna_start, dna_end, dna_strand, dna_cigar, dna_is_mapped, dna_is_not_multi, dna_nlen,
    dna_noNla_chr, dna_noNla_is_mapped, dna_noNla_is_not_multi, dna_noNla_cigar,
    rna1_R1_start, rna1_R1_end, rna1_R1_end_trim, rna1_R1_len_notrim, rna1_R1_len_trim,
    rna1_chr, rna1_start, rna1_end, rna1_strand, rna1_cigar, rna1_is_mapped, rna1_is_not_multi, rna1_nlen,
    rna2_R2_start, rna2_R2_end, rna2_R2_end_trim, rna2_R2_len_notrim, rna2_R2_len_trim,
    rna2_chr, rna2_start, rna2_end, rna2_strand, rna2_cigar, rna2_is_mapped, rna2_is_not_multi, rna2_nlen


Restriction of the DNA/RNA parts is a separated step.
You specify the list of restirction enzymes in `protocol.renz` with short names (or keys) and actual names of these restiriction enzymes.
Then, these enzymes are used for detection of restriction enzymes recognition sites around RNA and DNA parts.
Not all the combinations of restriction enzymes and segments are meaningful, thus you specify only the annotations that you want
in `run.restriction_check`.

Running annotation of restriction recodnition sites will create additional columns in hdf5 file.

1. If the restriction enzyme with name nla appeared for rna1 fragment in `run.restriction_check`, and it corresponds
to enzyme with palindromic recognition site (e.g. NlaIII), the following columns will be added:

.. code-block:: python
    rna1_end_nla_left, rna1_end_nla_right, rna1_start_nla_left, rna1_start_nla_right

2. If the restriction enzyme is not palindromic, then the DNA strand matters, and we hit + (p) and - (n) strands separately.
These columns will be added:

.. code-block:: python
    rna1_start_mmep_left, rna1_start_mmep_right, rna1_start_mmen_left, rna1_start_mmen_right,
    rna1_end_mmep_left, rna1_end_mmep_right, rna1_end_mmen_left, rna1_end_mmen_right

