lengths_dict = {'HeLa_M': 151,
 'HeLa_DRB': 151,
 'HeLa_G1': 151,
 'fibr_rep1_add': 133,
 'fibr_rep1': 125,
 'fibr_rep2': 125,
 'K562_rep1_add': 101,
 'K562_rep1_wo-ligase': 125,
 'K562_rep1': 75,
 'K562_rep1_rnase': 80,
 'K562_rep2': 75, 
 'sample': 151}

from sys import argv
prefix = argv[1]

formatting_dct = {'prefix': prefix}

from RedClib import *

import os
PATH_ABSOLUTE = os.path.dirname(os.path.realpath(__file__))

logging.debug("Working directory: {}".format(PATH_ABSOLUTE))

PATH_FASTQ = os.path.join(PATH_ABSOLUTE, 'data/fastq/')
PATH_TABLE = os.path.join(PATH_ABSOLUTE, 'data/tables/')
PATH_CBIN  = os.path.join(PATH_ABSOLUTE, 'data/cbin/')
PATH_HITS  = os.path.join(PATH_ABSOLUTE, 'data/cout/')
PATH_SEQS  = os.path.join(PATH_ABSOLUTE, 'data/filtered/')
PATH_SAM   = os.path.join(PATH_ABSOLUTE, 'data/sam/')
PATH_BED   = os.path.join(PATH_ABSOLUTE, 'data/bed/')

for d in [PATH_TABLE, PATH_CBIN, PATH_HITS, PATH_SEQS, PATH_SAM, PATH_BED]:
    if not os.path.isdir(d):
        os.mkdir(d)

nthreads = 3 # For mapping only

a = RedCprocessing()

infile1_fastq = PATH_FASTQ + '{prefix}_R1_001.fastq'.format(**formatting_dct)
infile2_fastq = PATH_FASTQ + '{prefix}_R2_001.fastq'.format(**formatting_dct)
file_fastq    = PATH_TABLE + '{prefix}.fastq.txt'.format(**formatting_dct)

a.prepare_FASTQ(infile1_fastq, infile2_fastq, file_fastq, prefix)

outfile1_fastuniq    = PATH_TABLE + '{prefix}.R1.50bp.fastq'.format(**formatting_dct)
outfile2_fastuniq    = PATH_TABLE + '{prefix}.R2.50bp.fastq'.format(**formatting_dct)
outfile_fastuniq_idx = PATH_TABLE + '{prefix}.unique_idx.txt'.format(**formatting_dct)

a.run_fastuniq(file_fastq, outfile1_fastuniq, outfile2_fastuniq, outfile_fastuniq_idx, 
               length=50, 
               fastuniq_binary="fastuniq") # Running the system version of fastuniq installed by conda

file1_trim = PATH_TABLE + '{prefix}.forward.trimmomatic_output.fastq'.format(**formatting_dct)
file2_trim = PATH_TABLE + '{prefix}.reverse.trimmomatic_output.fastq'.format(**formatting_dct)
file_trim = PATH_TABLE + '{prefix}.trimtable.txt'.format(**formatting_dct)

a.run_trimmomatic(infile1_fastq, infile2_fastq, file1_trim, file2_trim,
                  window=5, 
                  qual_th=26, 
                  minlen=0, 
                  trimmomatic_path=os.path.join(PATH_ABSOLUTE, "./bin/Trimmomatic-0.39/trimmomatic-0.39.jar") # Running version downloaded and installed by 01_compile.sh
                 )
a.prepare_TRIMOUT(file1_trim, file2_trim, file_trim, prefix, 'len15th26w5', file_fastq)

# Adaptors mapping

PATH_OLIGOS    = os.path.join(PATH_ABSOLUTE, 'data/oligos/')
PATH_CBINARIES = os.path.join(PATH_ABSOLUTE, 'bin/')

file_r1_for = PATH_HITS + '{prefix}_R1.for.txt'.format(**formatting_dct)
file_r1_rev = PATH_HITS + '{prefix}_R1.rev.txt'.format(**formatting_dct)
file_r2_for = PATH_HITS + '{prefix}_R2.for.txt'.format(**formatting_dct)
file_r2_rev = PATH_HITS + '{prefix}_R2.rev.txt'.format(**formatting_dct)


a.run_oligos_alignment(infile1_fastq, PATH_CBIN,
                   file_r1_for, file_r1_rev,
                   PATH_OLIGOS + 'for_20.fasta', PATH_OLIGOS + 'rev_20.fasta',
                   PATH_CBINARIES + 'fastq2bin', PATH_CBINARIES + 'fasta2bin', PATH_CBINARIES + 'align_universal',
                   seq_length=lengths_dict[prefix], n_primers=2,
                   left_shift=-6, right_shift=lengths_dict[prefix]-14,
                   mismatch_general=1, mode=1, report_len=20)

a.run_oligos_alignment(infile2_fastq, PATH_CBIN,
                   file_r2_for, file_r2_rev,
                   PATH_OLIGOS + 'for_20.fasta', PATH_OLIGOS + 'rev_20.fasta',
                   PATH_CBINARIES + 'fastq2bin', PATH_CBINARIES + 'fasta2bin', PATH_CBINARIES + 'align_universal',
                   seq_length=lengths_dict[prefix], n_primers=2,
                   left_shift=-6, right_shift=lengths_dict[prefix]-14,
                   mismatch_general=1, mode=1, report_len=20)

file_r1_rev16 = PATH_HITS + '{prefix}_R1.rev16.txt'.format(**formatting_dct)

a.run_oligos_alignment(infile1_fastq, PATH_CBIN,
                   file_r1_rev16, None,
                   PATH_OLIGOS + 'rev_16.fasta', None,
                   PATH_CBINARIES + 'fastq2bin', PATH_CBINARIES + 'fasta2bin', PATH_CBINARIES + 'align_universal',
                   seq_length=lengths_dict[prefix], n_primers=1,
                   left_shift=0, right_shift=lengths_dict[prefix]-14,
                   mismatch_general=1, mode=1, report_len=16)

# Bridge mapping

file_r1_br_for = PATH_HITS + '{prefix}_R1.37br_for.txt'.format(**formatting_dct)
file_r1_br_rev = PATH_HITS + '{prefix}_R1.37br_rev.txt'.format(**formatting_dct)
file_r2_br_for = PATH_HITS + '{prefix}_R2.37br_for.txt'.format(**formatting_dct)
file_r2_br_rev = PATH_HITS + '{prefix}_R2.37br_rev.txt'.format(**formatting_dct)

a.run_oligos_alignment(infile1_fastq, PATH_CBIN,
                   file_r1_br_for, file_r1_br_rev,
                   PATH_OLIGOS + 'br_37_for.fasta', PATH_OLIGOS + 'br_37_rev.fasta',
                   PATH_CBINARIES + 'fastq2bin', PATH_CBINARIES + 'fasta2bin', PATH_CBINARIES + 'align_universal',
                   seq_length=lengths_dict[prefix], n_primers=1,
                   left_shift=0, right_shift=lengths_dict[prefix]-14,
                   mismatch_general=1, mode=1, report_len=37)

a.run_oligos_alignment(infile2_fastq, PATH_CBIN,
                   file_r2_br_for, file_r2_br_rev,
                   PATH_OLIGOS + 'br_37_for.fasta', PATH_OLIGOS + 'br_37_rev.fasta',
                   PATH_CBINARIES + 'fastq2bin', PATH_CBINARIES + 'fasta2bin', PATH_CBINARIES + 'align_universal',
                   seq_length=lengths_dict[prefix], n_primers=1,
                   left_shift=0, right_shift=lengths_dict[prefix]-14,
                   mismatch_general=1, mode=1, report_len=37)

# Checking bridges for the proper last 2 letters

output_br_GA = PATH_HITS + "{prefix}_R1.GA.txt".format(**formatting_dct)

logging.debug("Checking presence of AG at the end of bridges...")

cmd_bgn_time = time.time()
with open(output_br_GA, 'w') as outf:
    with open(file_fastq, 'r') as in_f:
        with open(file_r1_br_for, 'r') as br_f:
            br = br_f.readline()

            t = in_f.readline()
            br = br_f.readline()

            while len(t) > 0:
                readF = t.split()[2]
                end_pos = int(br.split()[4])
                idx = br.split()[0]
                if end_pos > len(readF):
                    ret = 0
                else:
                    if readF[end_pos - 2:end_pos] == 'GA':
                        ret = 1
                    else:
                        ret = 0
                outf.write("{}\t{}\n".format(idx, ret))
                t = in_f.readline()
                br = br_f.readline()

cmd_end_time = time.time()
cmd_runtime = str(timedelta(seconds=cmd_end_time - cmd_bgn_time))
logging.debug("Bridge checking finished, outfile: {}\nRuntime: {}".format(output_br_GA, cmd_runtime))

# GGG mapping

file_r2_ggg = PATH_HITS + '{prefix}_R2.ggg.txt'.format(**formatting_dct)

a.run_oligos_alignment(infile2_fastq, PATH_CBIN,
                   file_r2_ggg, None,
                   PATH_OLIGOS + 'ggg.fasta', None,
                   PATH_CBINARIES + 'fastq2bin', PATH_CBINARIES + 'fasta2bin', PATH_CBINARIES + 'align_universal',
                   seq_length=lengths_dict[prefix], n_primers=1, left_shift=0, right_shift=3,
                   mismatch_general=0, mode=1, report_len=3)


# Checking complementary regions

outfile_R1 = PATH_HITS + '{prefix}.rna1_14bp.fastq'.format(**formatting_dct)
outfile_R2 = PATH_HITS + '{prefix}.rna_14bp.fastq'.format(**formatting_dct)

a.get_substrings_rna_fragments(file_fastq,
                               file_r1_br_for, file_r2_ggg,
                               outfile_R1, outfile_R2,
                               br_len=37, length=14)

a.reverse_complement_fastq(outfile_R1, outfile_R1+'_revcomp')
a.reverse_complement_fastq(outfile_R2, outfile_R2+'_revcomp')

os.remove(outfile_R1)
os.remove(outfile_R2)

file_r1_compl = PATH_HITS + '{prefix}_R1.complement.txt'.format(**formatting_dct)
file_r2_compl = PATH_HITS + '{prefix}_R2.complement.txt'.format(**formatting_dct)

a.run_oligos_alignment_paired(outfile_R2+'_revcomp', infile1_fastq,
                            PATH_CBIN, file_r1_compl,
                            PATH_CBINARIES + 'fastq2bin',
                            None,
                            PATH_CBINARIES + 'align_pairwise',
                            seq_length=lengths_dict[prefix], oligo_length=14,
                            left_shift=0, right_shift=lengths_dict[prefix]-14,
                            mismatch_general=1, mode=1, report_len=14)


a.run_oligos_alignment_paired(outfile_R1+'_revcomp', infile2_fastq,
                            PATH_CBIN, file_r2_compl,
                            PATH_CBINARIES + 'fastq2bin',
                            None,
                            PATH_CBINARIES + 'align_pairwise',
                            seq_length=lengths_dict[prefix], oligo_length=14,
                            left_shift=0, right_shift=lengths_dict[prefix]-14,
                            mismatch_general=1, mode=1, report_len=14)

# Get substring from R1 and R2

outfile1 = PATH_SEQS + '{prefix}.dna.fastq'.format(**formatting_dct)
outfile2 = PATH_SEQS + '{prefix}.rna.fastq'.format(**formatting_dct)
outfile3 = PATH_SEQS + '{prefix}.rna1.fastq'.format(**formatting_dct)

a.get_substrings(file_fastq, file_trim,
                file_r1_for, file_r2_for, file_r1_rev16,
                file_r1_br_for, file_r2_br_rev,
                file_r2_ggg, file_r1_compl, file_r2_compl,
                outfile1, outfile2, outfile3,
                sequence_length=lengths_dict[prefix], limit=14, add_dna="CATG")

a.get_substrings(file_fastq, file_trim,
                file_r1_for, file_r2_for, file_r1_rev16,
                file_r1_br_for, file_r2_br_rev,
                file_r2_ggg, file_r1_compl, file_r2_compl,
                outfile1+'.nonla', None, None,
                sequence_length=lengths_dict[prefix], limit=14, add_dna='')

outfile1_info = PATH_SEQS + '{prefix}.dna.info.txt'.format(**formatting_dct)
outfile2_info = PATH_SEQS + '{prefix}.rna.info.txt'.format(**formatting_dct)
outfile3_info = PATH_SEQS + '{prefix}.rna1.info.txt'.format(**formatting_dct)

a.write_substrings_info(file_fastq, file_trim,
                file_r1_for, file_r2_for, file_r1_rev16,
                file_r1_br_for, file_r2_br_rev,
                file_r2_ggg, file_r1_compl, file_r2_compl,
                outfile1_info, outfile2_info, outfile3_info,
                sequence_length=lengths_dict[prefix])

# Alignment

outfile1_map = PATH_SAM + '{prefix}.dna.sam'.format(**formatting_dct)
outfile2_map = PATH_SAM + '{prefix}.rna.sam'.format(**formatting_dct)
outfile3_map = PATH_SAM + '{prefix}.rna1.sam'.format(**formatting_dct)

a.align(outfile1, outfile1_map, 'dna',
        nthreads=nthreads)

a.align(outfile1+'.nonla', outfile1_map+'.nonla', 'dna',
        nthreads=nthreads)

a.align(outfile2, outfile2_map, 'rna',
        known_splice=os.path.join(PATH_ABSOLUTE, 'data/genome/spliced_genes_hisat2.gencode_v19.txt'),
        novel_splice=outfile2+'.novel_splice',
        nthreads=nthreads)

a.align(outfile3, outfile3_map, 'rna1',
        known_splice=os.path.join(PATH_ABSOLUTE, 'data/genome/spliced_genes_hisat2.gencode_v19.txt'),
        novel_splice=outfile3+'.novel_splice',
        nthreads=nthreads)

### Aligned to bed
outfile1_bed = PATH_BED + '{prefix}.dna.bed'.format(**formatting_dct)
outfile2_bed = PATH_BED + '{prefix}.rna.bed'.format(**formatting_dct)
outfile3_bed = PATH_BED + '{prefix}.rna1.bed'.format(**formatting_dct)

a.aligned_to_bed(outfile1_map+'.nonla', outfile1_bed+'.nonla', mode='dna')
a.aligned_to_bed(outfile1_map, outfile1_bed, mode='dna')
a.aligned_to_bed(outfile2_map, outfile2_bed, mode='rna')
a.aligned_to_bed(outfile3_map, outfile3_bed, mode='rna1')

logging.info("Processing done!")
