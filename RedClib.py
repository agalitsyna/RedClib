"""
rnadnaDB is a library of functions that are used for NGS RNA-DNA contacts data processing.

Processing of NGS data requires manipulations with huge text files were each read corresponds to one line.
bash commands such as awk, sed and sort are extremely effective for the manipulations of this kind. Thus most of
functions are wrappers for the set of low-level bash commands.

For the convenience and structure, some functions are combined into a class, although it is not used and not needed.

"""

import logging
logging.basicConfig(level=logging.DEBUG)

import time
from datetime import timedelta
import subprocess
import os

def raw_read_file(filename, columns, colnames, sep="\t", header=1, nrows=-1, modifiers=None, comment="#", bottom=0):
    """Read input text file with several columns line by line and saving as dictionary.

    :param filename: input file
    :param columns: numbers of columns to be parsed from file (ordered, 1-indexed as in awk)
    :param colnames: names of columns for the order from "columns" parameter
    :param sep: field separator in the input file
    :param header: number of rows to skip as header. Default: 1
    :param nrows: number of rows to read from file (after header). Default: -1 - read all lines
    :param modifiers: the function(s) to be applied to the columns from "columns" parameter. Can be None, one function or list of functions
    :param comment: skip the lines starting with "comment" character
    :param bottom: number of rows to skip from bottom (note that they should be parsed first)
    :return: dict with colnames as keys and list of columns in values.
    """
    assert len(columns) == len(colnames)
    assert header>=0
    ret = {k:[] for k in colnames}
    with open(filename, "r") as inf:
        for i in range(header+1):
            l = inf.readline().split(sep)
        nrows_passed = 0
        while len(l)>1:
            if nrows>0 and nrows_passed>nrows:
                break
        
            if len(l[0])>1 and l[0][0]==comment:
                l = inf.readline().split(sep)
                continue
            for i_mod, (i, col) in enumerate(zip(columns, colnames)):
                i1 = i-1 if i>0 else i
                val = l[i1].strip()
                if not modifiers is None:
                    val = modifiers[i_mod](val)
                ret[col].append(val)
            l = inf.readline().split(sep)
            nrows_passed += 1
    if bottom>0:
        ret = {k:ret[k][:-bottom] for k in colnames}
    return(ret)


def reconstruct_by_ids(raw_dct, key, reference_list, default_dct):
    """
    Filling missing values in the dict read by raw_read_file.
    :param raw_dct: initial dict
    :param key: key with indexes in the raw_dct
    :param reference_list: reference list of indexes containing all values
    :param default_dct: dict with default values for columns from raw_dct that will be filled for missing indexes
    :return:
    """
    ret_dct = {k:[] for k in raw_dct.keys()}
    ref_dct = {k:i for i, k in enumerate(raw_dct[key])}
    
    keys = [x for x in raw_dct.keys() if x!=key]
    
    for i in reference_list:
        try:
            idx = ref_dct[i]
            for k in keys:
                ret_dct[k].append(raw_dct[k][idx])
        except Exception as e:
            for k in keys:
                ret_dct[k].append(default_dct[k])
        ret_dct[key].append(i)
        
    return ret_dct


class RedCprocessing:
    """
    Formal class for organisation of functions for text files data processing.
    """

    def __init__(self):
        """
        Basic constructor for fake RedCprocessing class.
        """
        pass

    @staticmethod
    def call_and_check_errors(command):
        """ Basic function to run bash commands from Python and report errors. """
        logging.debug("Starting command:\n\t{}\n".format(command))
        process = subprocess.Popen(command,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                shell=True,
                                executable="/bin/bash")
        (stdout, stderr) = process.communicate()
        if stderr:
            logging.info("Stderr is not empty for command:\n\t{}".format(command))
            logging.info("Check stderr:\n{}".format(stderr))
            return stderr   # Something might be wrong, error. Some commands might have output to stderr for normal runs.
        else:
            return 0        # Definitely no error, great!

    def run_and_check_command(self, command):

        cmd_bgn_time = time.time()
        is_err = self.call_and_check_errors(command)
        cmd_end_time = time.time()
        cmd_runtime = timedelta(seconds=cmd_end_time - cmd_bgn_time)

        return is_err, cmd_runtime



    def prepare_FASTQ(self, file1, file2, file_tmp, id_exp=None):
        """Write input pair-end FASTQ files into a single table with sequences.
        This step is time-consuming and output files are huge.

        :param file1: File with forward (R1) reads
        :param file2: File with reverse (R2) reads
        :param file_tmp: File for output with all reads, ids and qualities
        :param id_exp: Name of experiment (will be written into the file
        :return: None
        """

        fmt_dct = {"id_exp": id_exp, "file1": file1, "file2": file2, "file_tmp": file_tmp}

        logging.info("Preparing FASTQ files ({id_exp}):\n\t{file1}\n\t{file2}\n\tto:\t{file_tmp}".format(**fmt_dct))

        if "gz" in file1:
            if id_exp:
                command = str(r"""paste <(gzip -dc {file1} | sed 'N;N;N;s/\n/ /g' | awk 'BEGIN{{OFS="\t"}}{{print $1, "{id_exp}", $3, $5}}') <(sed 'N;N;N;s/\n/ /g' {file2} | awk 'BEGIN{{OFS="\t"}}{{print $3, $5}}') > {file_tmp}""") \
                      .format(**fmt_dct)
            else:
                command = str(r"""paste <(gzip -dc {file1} | sed 'N;N;N;s/\n/ /g' | awk 'BEGIN{{OFS="\t"}}{{print $1, $3, $5}}') <(sed 'N;N;N;s/\n/ /g' {file2} | awk 'BEGIN{{OFS="\t"}}{{print $3, $5}}') > {file_tmp}""")\
                      .format(**fmt_dct)
        else:
            if id_exp:
                command = str(r"""paste <(sed 'N;N;N;s/\n/ /g' {file1} | awk 'BEGIN{{OFS="\t"}}{{print $1, "{id_exp}", $3, $5}}') <(sed 'N;N;N;s/\n/ /g' {file2} | awk 'BEGIN{{OFS="\t"}}{{print $3, $5}}') > {file_tmp}""") \
                      .format(**fmt_dct)
            else:
                command = str(r"""paste <(sed 'N;N;N;s/\n/ /g' {file1} | awk 'BEGIN{{OFS="\t"}}{{print $1, $3, $5}}') <(sed 'N;N;N;s/\n/ /g' {file2} | awk 'BEGIN{{OFS="\t"}}{{print $3, $5}}') > {file_tmp}""")\
                      .format(**fmt_dct)

        is_err, cmd_runtime = self.run_and_check_command(command)

        if is_err:
            logging.debug("Preparing FASTQ files. Runtime: {}. Filenames:\n\t{}\n\t{}".format(cmd_runtime, file1, file2))
        else:
            logging.debug("Preparing FASTQ files finished successfully. Runtime: {}".format(cmd_runtime))


    def get_fastq_substr(self, infile_fastq, outfile_R1, outfile_R2, length=50):

        formatting_dct = {"infile_fastq": infile_fastq,
                          "outfile_R1": outfile_R1,
                          "outfile_R2": outfile_R2,
                          "length": length}

        command1 = str(r"""awk 'BEGIN{{OFS="\n"}}{{print $1, substr($3, 0, {length}), "+", substr($4, 0, {length})}}' {infile_fastq} > {outfile_R1}""".format(**formatting_dct))
        command2 = str(r"""awk 'BEGIN{{OFS="\n"}}{{print $1, substr($5, 0, {length}), "+", substr($6, 0, {length})}}' {infile_fastq} > {outfile_R2}""".format(**formatting_dct))

        cmd_bgn_time = time.time()
        is_err1 = self.call_and_check_errors(command1)
        is_err2 = self.call_and_check_errors(command2)
        cmd_end_time = time.time()
        cmd_runtime = str(timedelta(seconds=cmd_end_time - cmd_bgn_time))

        if is_err1 or is_err2:
            logging.debug("Subsequences retrieval failed. Runtime: {}. Filenames:\n\t{}\n\t{}".format(cmd_runtime, infile_fastq, outfile_R1))
        else:
            logging.debug("Subsequences retrieval finished successfully. Runtime: {}".format(cmd_runtime))


    def run_fastuniq(self, infile_fastq, outfile1, outfile2, outfile_idx,
                     length=50, fastuniq_binary="fastuniq", remove_intermediary=True):
        """
        Running fastuniq command for marking PCR duplicates.

        :param infile1:
        :param infile2:
        :param outfile1:
        :param outfile2:
        :param outfile_idx:
        :param length: number of basepairs to be cut from 5'-end of each read. Full-length reads are used if length==0 or smaller
        :param fastuniq_binary: binary for fastuniq
        :remove_intermediary: remove temoprary intermediary files, default: True
        :return: None
        """

        formatting_dct = {"infile_fastq":infile_fastq,
                          "outfile_idx":outfile_idx,
                          "length":length,
                          "fastuniq_binary":fastuniq_binary,
                          "outfile1":outfile1,
                          "outfile2":outfile2,
                          "outfile_idx_tmp":outfile_idx+'.tmp',
                          "outfile_idx_fastq1":outfile1+".fastquniq",
                          "outfile_idx_fastq2":outfile2+".fastquniq"}


        # Write shortened fastq files (not needed if length <=0)
        if length>0:
            self.get_fastq_substr(infile_fastq, outfile1, outfile2, length=length)

        # Write file with paired fastq files for analysis:
        with open('{outfile_idx_tmp}'.format(**formatting_dct), 'w') as tmp:
            tmp.write("{outfile1}\n{outfile2}".format(**formatting_dct))

        # Run fastuniq
        command1 = str(r"""{fastuniq_binary} -i {outfile_idx_tmp} -tq -c 0 -o {outfile_idx_fastq1} -p {outfile_idx_fastq2}""".format(**formatting_dct))
        # Retreive the list of ids of unique reads
        command2 = str(r"""awk 'NR % 4 == 1' {outfile_idx_fastq1} | gawk '{{match($0, "@([^ ,/]+)", a)}} {{print a[1]}}' > {outfile_idx}""".format(**formatting_dct))

        cmd_bgn_time = time.time()
        is_err1 = self.call_and_check_errors(command1)
        is_err2 = self.call_and_check_errors(command2)
        cmd_end_time = time.time()
        cmd_runtime = str(timedelta(seconds=cmd_end_time - cmd_bgn_time))

        if remove_intermediary:
            for k in ["outfile_idx_tmp", "outfile_idx_fastq1", "outfile_idx_fastq2"]:
                os.remove(formatting_dct[k])
            if length>0:
                for k in ["outfile1", "outfile2"]:
                    os.remove(formatting_dct[k])

        if is_err1 or is_err2:
            logging.debug("Fastuniq run with stderr output. Stderr:\n\t{}\n\t{} \n\tRuntime: {}. Filenames:\n\t{}\n\t{}".format(is_err1, is_err2, cmd_runtime, infile_fastq, outfile_idx))
        else:
            logging.debug("Fastuniq run successfully. Runtime: {}. Filenames:\n\t{}\n\t{}".format(cmd_runtime, infile_fastq, outfile_idx))


    def run_trimmomatic(self, file1, file2, outfile1, outfile2,
                        trimmomatic_path="./bin/Trimmomatic-0.39/trimmomatic-0.39.jar",
                        window=5, qual_th=26, minlen=15):
        """
        Run trimmomatic (paired mode) for two files with reads separately

        :param file1:
        :param outfile1:
        :param trimmomatic_path: path to trimmomatic .jar file
        :return:
        """

        logging.info("Running paired trimmomatic ({2}) \n for files: {0} {1} \n with parameters: {3} {4} {5}".format(
            file1, file2, trimmomatic_path, window, qual_th, minlen))

        command = "java -jar {0} PE -phred33 {1} {2} {3} {4} {5} {6} SLIDINGWINDOW:{7}:{8} MINLEN:{9}".format(
            trimmomatic_path, file1, file2,
            outfile1, outfile1+".unpaired",
            outfile2, outfile2+".unpaired",
            int(window), int(qual_th), int(minlen))

        is_err, cmd_runtime = self.run_and_check_command(command)

        if is_err:
            logging.debug(
                "Trimmomatic run with stderr output for file: {} \n Possible error: \n\t{}\n\t"
                "Trimmomatic runtime: {}.".format(file1, is_err, cmd_runtime))
        else:
            logging.debug("Running trimmomatic succeeded for file: {} (outfile: {})"
                "Trimmomatic runtime: {}.".format(file1, outfile1, cmd_runtime))


    #############################################
    #
    #   beginning of the TRIMOUT section
    #
    #############################################

    def prepare_TRIMOUT(self, file1, file2, file_tmp, id_exp, id_trim_meth, fastq_table, remove_intermediary=True):
        # Add: awk 'NR==FNR {vals[$3] = $4 " " $5 " " $6 " " $7 ; next} !($1 in vals) {vals[$1] = "0 0 0 0"} {$(NF+1) = vals[$1]; print}' <(head -n 1000 AG_9.table.txt) <(head -n 1000 AG_9_S70.fastq.txt)

        logging.info("Preparing TRIMOUT files for upload (%s, %s):\n\t%s\n\t%s\n\tto:\t%s", id_exp, id_trim_meth, file1, file2, file_tmp)

        command1 = str(r"""paste <(sed -n '1~4p' {0} | awk 'BEGIN{{OFS="\t"}}{{print "{3}", "{4}", $1}}') """
                      r"""<(sed -n '2~4p' {0} | awk 'BEGIN{{OFS="\t"}}{{print 0, length($0);}}') """
                      r"""<(sed -n '2~4p' {1} | awk 'BEGIN{{OFS="\t"}}{{print 0, length($0);}}') > {2}""")\
                .format(file1, file2, file_tmp+".tmp", id_exp, id_trim_meth)

        command2 = str(r"""awk 'NR==FNR {{vals[$3] = $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 ; next}} """
                    r"""!($1 in vals) {{vals[$1] = "{3}" "\t" "{4}" "\t" $1 "\t" "0\t0\t0\t0"}} """
                    r"""{{$(NF+1) = vals[$1]; print vals[$1]}}' """
                    r"""{0} {1} > {2}"""\
                .format(file_tmp+".tmp", fastq_table, file_tmp, id_exp, id_trim_meth))

        is_err1, cmd_runtime1 = self.run_and_check_command(command1)
        is_err2, cmd_runtime2 = self.run_and_check_command(command2)

        if remove_intermediary:
            try:
                for f in [file_tmp+".tmp", file1, file2, file1+'.unpaired', file2+'.unpaired']:
                    os.remove(f)
            except Exception as e:
                logging.error("Deletion of intermediary file {} failed.".format(f))

        if is_err1 or is_err2:
            logging.error("Preparing FASTQ files for upload run with stderr: \n\t{}\n\t{}\n\t Runtime: {}. Filenames: {} {}".format(is_err1, is_err2, cmd_runtime1+cmd_runtime2, file1, file2))
        else:
            logging.debug("Preparing FASTQ files for upload finished successfully. Runtime: {} {} Total: {}".format(cmd_runtime1 ,cmd_runtime2, cmd_runtime1+cmd_runtime2))

    ###############################
    #
    # Oligos alignment tool
    #
    ###############################

    def run_oligos_alignment(self, infile, bin_tmp_dir,
                          outfile_for, outfile_rev,
                          primer_for, primer_rev,
                          bin_infile2bin, bin_primer2bin, bin_align,
                          seq_length, n_primers,
                          left_shift=0, right_shift=80,
                          mismatch_general=0, mode=1, report_len=20):

        # dct_formatting = {"infile":infile,
        #                   "bin_tmp_dir": bin_tmp_dir,
        #                   "outfile_for": outfile_for,
        #                   "outfile_rev": outfile_rev,
        #                   "primer_for": primer_for,
        #                   "primer_rev": primer_rev,
        #                   "bin_infile2bin": bin_infile2bin,
        #                   "bin_primer2bin": bin_primer2bin,
        #                   "bin_align": bin_align,
        #                   "seq_length": seq_length,
        #                   "n_primers": n_primers,
        #                   "left_shift": left_shift,
        #                   "right_shift": right_shift,
        #                   "mismatch_general": mismatch_general,
        #                   "mismatch_add": 0, # deprecated parameter. Enter whatever interpretable as unsigned int. TODO: fix
        #                   "mode": mode,
        #                   "report_len": report_len}

        mismatch_add = 0 # deprecated parameter of align_universal. Can be whatever interpretable as unsigned int. TODO: fix
        tmp_filename = bin_tmp_dir + infile.split("/")[-1] + ".bin"

        command1 = f"{bin_infile2bin} {infile} {tmp_filename}"

        tmp_for = bin_tmp_dir + primer_for.split("/")[-1] + ".bin"
        if primer_rev:
            tmp_rev = bin_tmp_dir + primer_rev.split("/")[-1] + ".bin"
        else:
            tmp_rev = ""
        command2 = f"{bin_primer2bin} {primer_for} {tmp_for}"
        command3 = f"{bin_primer2bin} {primer_rev} {tmp_rev}"

        if seq_length==101:
            bin_seq_length = 15
        elif seq_length==151:
            bin_seq_length = 21
        elif seq_length==125:
            bin_seq_length = 18
        elif seq_length==80:
            bin_seq_length = 12
        elif seq_length==133:
            bin_seq_length = 19
        elif seq_length==251:
            bin_seq_length = 34
        else:
            raise Exception("Length not found!")

        command4 = f"{bin_align} {tmp_for} {tmp_filename} {mode} {seq_length} {bin_seq_length} {n_primers} \
{left_shift} {right_shift} {mismatch_general} {mismatch_add} {report_len} > {outfile_for}"
        command5 = f"{bin_align} {tmp_rev} {tmp_filename} {mode} {seq_length} {bin_seq_length} {n_primers} \
{left_shift} {right_shift} {mismatch_general} {mismatch_add} {report_len} > {outfile_rev}"

        if not bin_infile2bin:
            command1 = "echo -"
        if not bin_primer2bin:
            command2 = "echo -"
        if not primer_rev:
            command3 = "echo -"
            command5 = "echo -"

        is_err1, cmd_runtime1 = self.run_and_check_command(command1)
        if is_err1:
            logging.error("File to bin failed: {}.".format(is_err1))
            return 1

        is_err2, cmd_runtime2 = self.run_and_check_command(command2)
        if is_err2:
            logging.error("For to bin failed: {}.".format(is_err2))
            return 1

        is_err3, cmd_runtime3 = self.run_and_check_command(command3)
        if is_err3:
            logging.error("Rev to bin failed: {}.".format(is_err3))
            return 1

        is_err4, cmd_runtime4 = self.run_and_check_command(command4)
        if is_err4:
            logging.error("Align forward failed: {}.".format(is_err4))
            return 1

        is_err5, cmd_runtime5 = self.run_and_check_command(command5)
        if is_err5:
            logging.error("Align forward failed: {}.".format(is_err5))
            return 1

        logging.debug(f"Alignment of oligos finished! Runtime of steps: {cmd_runtime1} {cmd_runtime4} {cmd_runtime5}. \
                        Total: {cmd_runtime1+cmd_runtime4+cmd_runtime5}")
        return 0

    def run_oligos_alignment_paired(self, infile1, infile2, bin_tmp_dir,
                          outfile,
                          bin_infile12bin, bin_infile22bin, bin_align,
                          seq_length, oligo_length,
                          left_shift=0, right_shift=80,
                          mismatch_general=0, mismatch_add=1, mode=1, report_len=15):

        tmp_filename1 = bin_tmp_dir + infile1.split("/")[-1] + ".bin"
        command1 = "{} {} {}".format(bin_infile12bin, infile1, tmp_filename1)

        tmp_filename2 = bin_tmp_dir + infile2.split("/")[-1] + ".bin"
        command2 = "{} {} {}".format(bin_infile22bin, infile2, tmp_filename2)

        if seq_length==101:
            bin_seq_length = 15
        elif seq_length==151:
            bin_seq_length = 21
        elif seq_length==125:
            bin_seq_length = 18
        elif seq_length==80:
            bin_seq_length = 12
        elif seq_length==133:
            bin_seq_length = 19
        elif seq_length==251:
            bin_seq_length = 34
        else:
            raise Exception("Length not found!")

        oligo_length = ((oligo_length+3)//8)+1

        command3 = "{} {} {} {} {} {} {} {} {} {} {} {} > {}".format(bin_align,
                                                                     tmp_filename1,
                                                                     tmp_filename2,
                                                                     mode,
                                                                     seq_length, bin_seq_length,
                                                                     oligo_length,
                                                                     left_shift, right_shift, mismatch_general,
                                                                     mismatch_add, report_len,
                                                                     outfile)

        if not bin_infile12bin:
            command1 = "echo -"
        if not bin_infile22bin:
            command2 = "echo -"

        is_err1, cmd_runtime1 = self.run_and_check_command(command1)
        if is_err1:
            logging.error("File to bin failed: {}.".format(is_err1))
            return 1

        is_err2, cmd_runtime2 = self.run_and_check_command(command2)
        if is_err2:
            logging.error("For to bin failed: {}.".format(is_err2))
            return 1

        is_err3, cmd_runtime3 = self.run_and_check_command(command3)
        if is_err3:
            logging.error("Rev to bin failed: {}.".format(is_err3))
            return 1

        logging.debug("Alignment of oligos finished! Runtime of steps: {} {} {}. Total: {}".format(cmd_runtime1,
                                                                           cmd_runtime2,
                                                                           cmd_runtime3,
                                                                           cmd_runtime1+cmd_runtime2+cmd_runtime3))
        return 0


    def get_substrings(self, fastq_table, trim_table,
                       cout_r1_for, cout_r2_for, cout_r1_rev,
                       cout_r1_br_for, cout_r2_br_rev,
                       cout_r2_ggg, cout_compl_1, cout_compl_2,
                       outfile_dna1, outfile_rna2, outfile_rna1,
                       sequence_length=80, br_len=37, limit=14, add_dna="CATG"):

        formatting_dct = {"fastq_table":fastq_table,
                          "trim_table": trim_table,
                          "cout_r1_for": cout_r1_for,
                          "cout_r2_for": cout_r2_for,
                          "cout_r1_rev": cout_r1_rev,
                          "cout_r1_br_for": cout_r1_br_for,
                          "cout_r2_br_rev": cout_r2_br_rev,
                          "cout_r2_ggg": cout_r2_ggg,
                          "cout_compl_1": cout_compl_1,
                          "cout_compl_2": cout_compl_2,
                          "outfile_dna1": outfile_dna1,
                          "outfile_rna2": outfile_rna2,
                          "outfile_rna1": outfile_rna1,
                          "sequence_length": sequence_length,
                          "br_len": br_len,
                          "limit": limit,
                          "add_dna": add_dna,
                          "add_dna_qual": "".join(["~" for i in range(len(add_dna))])}

        logging.debug("Running getting substrings from input files (e.g. {fastq_table}). Writing output to: {outfile_dna1} {outfile_rna2}".format(**formatting_dct))

        if outfile_dna1:
            command1 = str(r"""paste """
                       """<(awk '{{print $1, $3, $4}}' {fastq_table}) """                      # data on fastq sequence: id $1, seq $2, qual $3
                       """<(head -n -1 {cout_r1_for} | tail -n +2 | awk '{{print $5+1}}') """  # end of forward adaptor on r1 ($4)
                       """<(awk '{{print $5}}' {trim_table}) """                               # position for trimming for forward read ($5)
                       """<(head -n -1 {cout_r1_br_for} | tail -n +2 | awk '{{print $4}}') """ # start of forward bridge on r1 ($6)
                        # DNA bgn = the largest of: end of forward adaptor ($4), if it's found; read start (1) otherwise
                        # DNA end = the smallest of: start of forward bridge on r1 ($6), if it's found; trimming position ($5) otherwise
                       """| awk 'BEGIN{{OFS="\\n";}} {{bgn=1; if ($4<500) bgn=$4;"""
                       """end=$5; if ($6<end) end=$6; """
                       """if (end-bgn+1>={limit}) print $1, substr($2, bgn, end-bgn+1)"{add_dna}", "+", substr($3, bgn, end-bgn+1)"{add_dna_qual}"}}' > {outfile_dna1}"""
                       ).format(**formatting_dct)
        else:
            command1 = "echo -"

        if outfile_rna2:
            command2 = str( r"""paste """
                        """<(awk '{{print $1, $5, $6}}' {fastq_table}) """                      # data on fastq sequence: id $1, seq $2, qual $3
                        """<(head -n -1 {cout_r2_ggg} | tail -n +2 | awk '{{print $5+1}}') """  # end of ggg hit on r2 ($4)
                        """<(head -n -1 {cout_r2_for} | tail -n +2 | awk '{{print $5+1}}') """  # end of forward adaptor with GGG hit on r2 ($5)
                        """<(awk '{{print $7}}' {trim_table}) """                               # position for trimming of reverse read ($6)
                        """<(head -n -1 {cout_r2_br_rev} | tail -n +2 | awk '{{print $4}}') """ # start of reverse bridge on r2 ($7)
                        """<(head -n -1 {cout_compl_2} | tail -n +2 | awk '{{print $5}}') """ # end of rna1 complementary region hit on r2 ($8)
                        # RNA bgn = the largest of: end of ggg hit ($4), if it's found, or end of adaptor with GGG at the end, if it's found; read start (1) otherwise
                        # RNA end = the smallest of: start of reverse bridge on r2 ($7), if it's found, or end of rna complementary region hit on r2 ($8); trimming position ($6) otherwise
                        """| awk 'BEGIN{{OFS="\\n";}} {{bgn=1; if ($4<500) bgn=$4; if ($5<500 && $5>bgn) bgn=$5;"""
                        """end=$6; if ($7<end) end=$7; if ($8<end) end=$8;"""
                        """if (end-bgn+1>={limit}) print $1, substr($2, bgn, end-bgn+1), "+", substr($3, bgn, end-bgn+1)}}' > {outfile_rna2}"""
                        ).format(**formatting_dct)
        else:
            command2 = "echo -"

        if outfile_rna1:
            command3 = str(r"""paste """
                       """<(awk '{{print $1, $3, $4}}' {fastq_table}) """                      # data on fastq sequence: id $1, seq $2, qual $3
                       """<(awk '{{print $5}}' {trim_table}) """                               # position for trimming of forward read ($4)
                       """<(head -n -1 {cout_r1_br_for} | tail -n +2 | awk '{{print $4}}') """ # start of forward bridge on r1 ($5)
                       """<(head -n -1 {cout_r1_rev} | tail -n +2 | awk '{{print $4}}') """    # start of reverse adaptor with GGG hit on r1 ($6)
                       """<(head -n -1 {cout_compl_1} | tail -n +2 | awk '{{print $5}}') """ # end of rna complementary region hit on r1 ($7)
                        # RNA1 bgn = end of forward bridge on r1 ($5), if it's found; trimming start otherwise (to guarantee zero length)
                        # RNA1 end = the smallest of: start of reverse adaptor with GGG hit on r1 ($6), end of rna complementary region hit on r1 ($7) or position for trimming of forward read ($4)
                       """| awk 'BEGIN{{OFS="\\n";}} {{bgn=$4; if ($5<500) bgn=$5+{br_len}+1;"""
                       """end=$4; if ($6<end) end=$6; if ($7<end) end=$7;"""
                       """if (end-bgn+1>={limit}) print $1, substr($2, bgn, end-bgn+1), "+", substr($3, bgn, end-bgn+1)}}' > {outfile_rna1}"""
                       ).format(**formatting_dct)
        else:
            command3 = "echo -"

        is_err1, cmd_runtime1 = self.run_and_check_command(command1)

        if is_err1:
            logging.error("Substr selection failed: {}.".format(is_err1))
            return 1

        is_err2, cmd_runtime2 = self.run_and_check_command(command2)

        if is_err2:
            logging.error("Substr selection failed: {}.".format(is_err2))
            return 1

        is_err3, cmd_runtime3 = self.run_and_check_command(command3)

        if is_err3:
            logging.error("Substr selection failed: {}.".format(is_err3))
            return 1

        logging.debug("Substr selection finished! Runtime of steps: {} {} {}. Total: {}".format(cmd_runtime1,
                                                                                                cmd_runtime2,
                                                                                                cmd_runtime3,
                                                                                                cmd_runtime1 + cmd_runtime2 + cmd_runtime3 ))
        return 0

    def get_substrings_rna_fragments(self, fastq_table,
                       cout_r1_br_for,
                       cout_r2_ggg,
                       outfile_rna1, outfile_rna2,
                       br_len=37, length=15):

        formatting_dct = {"fastq_table":fastq_table,
                          "cout_r1_br_for": cout_r1_br_for,
                          "cout_r2_ggg": cout_r2_ggg,
                          "outfile_rna2": outfile_rna2,
                          "outfile_rna1": outfile_rna1,
                          "br_len": br_len,
                          "length": length,
                          "extraN": "N"*(500+3*length)} # Addd extra NNNN

        logging.debug("Running getting rna substrings from input files (e.g. {fastq_table}). "
                      "Writing output to: {outfile_rna1} {outfile_rna2}".format(**formatting_dct))


        if outfile_rna1:
            command1 = str(r"""paste """
                       """<(awk '{{print $1, $3, $4}}' {fastq_table}) """
                       """<(head -n -1 {cout_r1_br_for} | tail -n +2 | awk '{{print $4}}') """ # r1 br for $4
                       """| awk 'BEGIN{{OFS="\\n";}} {{bgn=1; if ($4<500) bgn=$4+{br_len}+1;"""
                       """print $1, substr($2"{extraN}", bgn, {length}), "+", substr($3"{extraN}", bgn, {length})}}' > {outfile_rna1}"""
                       ).format(**formatting_dct)
        else:
            command1 = "echo -"

        if outfile_rna2:
            command2 = str( r"""paste """
                        """<(awk '{{print $1, $5, $6}}' {fastq_table}) """ # fastq
                        """<(head -n -1 {cout_r2_ggg} | tail -n +2 | awk '{{print $5+1}}') """ # r2 ggg $4
                        """| awk 'BEGIN{{OFS="\\n";}} {{bgn=1; if ($4<500) bgn=$4;"""
                        """print $1, substr($2"{extraN}", bgn, {length}), "+", substr($3"{extraN}", bgn, {length})}}' > {outfile_rna2}"""
                        ).format(**formatting_dct)
        else:
            command2 = "echo -"


        is_err1, cmd_runtime1 = self.run_and_check_command(command1)

        if is_err1:
            logging.error("Substr selection failed: {}.".format(is_err1))
            return 1

        is_err2, cmd_runtime2 = self.run_and_check_command(command2)

        if is_err2:
            logging.error("Substr selection failed: {}.".format(is_err2))
            return 1

        logging.debug("Substr selection finished! Runtime of steps: {} {}. Total: {}".format(cmd_runtime1,
                                                                                                cmd_runtime2,
                                                                                                cmd_runtime1 + cmd_runtime2 ))
        return 0


    def write_substrings_info(self, fastq_table, trim_table,
                       cout_r1_for, cout_r2_for, cout_r1_rev,
                       cout_r1_br_for, cout_r2_br_rev,
                       cout_r2_ggg, cout_compl_1, cout_compl_2,
                       outfile_dna1, outfile_rna2, outfile_rna1,
                       sequence_length=80, br_len=37, limit=14):

        formatting_dct = {"fastq_table":fastq_table,
                          "trim_table": trim_table,
                          "cout_r1_for": cout_r1_for,
                          "cout_r2_for": cout_r2_for,
                          "cout_r1_rev": cout_r1_rev,
                          "cout_r1_br_for": cout_r1_br_for,
                          "cout_r2_br_rev": cout_r2_br_rev,
                          "cout_r2_ggg": cout_r2_ggg,
                          "cout_compl_1": cout_compl_1,
                          "cout_compl_2": cout_compl_2,
                          "outfile_dna1": outfile_dna1,
                          "outfile_rna2": outfile_rna2,
                          "outfile_rna1": outfile_rna1,
                          "sequence_length": sequence_length,
                          "br_len": br_len,
                          "limit": limit}

        logging.debug("Running writing substrings metainfo from input files (e.g. {fastq_table}). Writing output to: {outfile_dna1} {outfile_rna2}".format(**formatting_dct))

        if outfile_dna1:
            command1 = str(r"""paste """
                           """<(awk '{{print $1, $3, $4}}' {fastq_table}) """                      # data on fastq sequence: id $1, seq $2, qual $3
                           """<(head -n -1 {cout_r1_for} | tail -n +2 | awk '{{print $5+1}}') """  # end of forward adaptor on r1 ($4)
                           """<(awk '{{print $5}}' {trim_table}) """                               # position for trimming for forward read ($5)
                           """<(head -n -1 {cout_r1_br_for} | tail -n +2 | awk '{{print $4}}') """ # start of forward bridge on r1 ($6)
                           # DNA bgn = end of forward adaptor ($4), if it's found; read start (1) otherwise
                           # DNA end = start of forward bridge on r1 ($6), if it's found; read length otherwise (sequence_length)
                           """| awk 'BEGIN{{OFS="\\t";}} {{bgn=1; if ($4<500) bgn=$4;"""
                           """end={sequence_length}; if ($6<end) end=$6;"""
                           """print $1, bgn, end, end-bgn+1, $5, $5-bgn+1}}' > {outfile_dna1}"""
                           ).format(**formatting_dct)
        else:
            command1 = "echo -"

        if outfile_rna2:
            command2 = str(r"""paste """
                           """<(awk '{{print $1, $5, $6}}' {fastq_table}) """                       # data on fastq sequence: id $1, seq $2, qual $3
                           """<(head -n -1 {cout_r2_ggg} | tail -n +2 | awk '{{print $5+1}}') """   # end of ggg hit on r2 ($4)
                           """<(head -n -1 {cout_r2_for} | tail -n +2 | awk '{{print $5+1}}') """   # end of forward adaptor with GGG hit on r2 ($5)
                           """<(awk '{{print $7}}' {trim_table}) """                                # position for trimming of reverse read ($6)
                           """<(head -n -1 {cout_r2_br_rev} | tail -n +2 | awk '{{print $4}}') """  # start of reverse bridge on r2 ($7)
                           """<(head -n -1 {cout_compl_2} | tail -n +2 | awk '{{print $5}}') """  # end of rna1 complementary region hit on r2 ($8)
                           # RNA bgn = the largest of: end of ggg hit ($4), if it's found, or end of adaptor with GGG at the end, if it's found; read start (1) otherwise
                           # RNA end = start of reverse bridge on r2 ($7), if it's found, or end of rna complementary region hit on r2 ($8); read length (sequence_length) otherwise
                           """| awk 'BEGIN{{OFS="\\t";}} {{bgn=1; if ($4<500) bgn=$4; if ($5<500 && $5>bgn) bgn=$5;"""
                           """end={sequence_length}; if ($7<end) end=$7; if ($8<end) end=$8;"""
                           """print $1, bgn, end, end-bgn+1, $6, $6-bgn+1}}' > {outfile_rna2}"""
                           ).format(**formatting_dct)
        else:
            command2 = "echo -"

        if outfile_rna1:
            command3 = str(r"""paste """
                           """<(awk '{{print $1, $3, $4}}' {fastq_table}) """                       # data on fastq sequence: id $1, seq $2, qual $3
                           """<(awk '{{print $5}}' {trim_table}) """                                # position for trimming of forward read ($4)
                           """<(head -n -1 {cout_r1_br_for} | tail -n +2 | awk '{{print $4}}') """  # start of forward bridge on r1 ($5)
                           """<(head -n -1 {cout_r1_rev} | tail -n +2 | awk '{{print $4}}') """     # start of reverse adaptor with GGG hit on r1 ($6)
                           """<(head -n -1 {cout_compl_1} | tail -n +2 | awk '{{print $5}}') """  # end of rna complementary region hit on r1 ($7)
                           # RNA1 bgn = end of forward bridge on r1 ($5), if it's found; 1 otherwise
                           # RNA1 end = the smallest of: start of reverse adaptor with GGG hit on r1 ($6), end of rna complementary region hit on r1 ($7) or read length (sequence_length)
                           """| awk 'BEGIN{{OFS="\\t";}} {{bgn=1; if ($5<500) bgn=$5+{br_len}+1;"""
                           """end={sequence_length}; if ($6<end) end=$6; if ($7<end) end=$7;"""
                           """print $1, bgn, end, end-bgn+1, $4, $4-bgn+1}}' > {outfile_rna1}"""
                           ).format(**formatting_dct)
        else:
            command3 = "echo -"

        is_err1, cmd_runtime1 = self.run_and_check_command(command1)

        if is_err1:
            logging.error("Substr selection failed: {}.".format(is_err1))
            return 1

        is_err2, cmd_runtime2 = self.run_and_check_command(command2)

        if is_err2:
            logging.error("Substr selection failed: {}.".format(is_err2))
            return 1

        is_err3, cmd_runtime3 = self.run_and_check_command(command3)

        if is_err3:
            logging.error("Substr selection failed: {}.".format(is_err3))
            return 1

        logging.debug("Substr selection finished! Runtime of steps: {} {} {}. Total: {}".format(cmd_runtime1,
                                                                                                cmd_runtime2,
                                                                                                cmd_runtime3,
                                                                                                cmd_runtime1 + cmd_runtime2 + cmd_runtime3 ))
        return 0


    def align(self, infile, outfile, mode,
              genome="./data/genome/hisat2/hg19", # hisat2 index file
              nthreads=15,
              novel_splice="",
              known_splice="./data/genome/spliced_genes_hisat2.gencode_v19.txt",
              bin_hisat2="hisat2",
              custom_command=None):
        
        if isinstance(custom_command, str):
            command=custom_command
        elif mode=="dna":
            command = r"""{} -p {} -x {} --no-spliced-alignment -k 100 --no-softclip """ \
                      r"""-U {} > {}""".format(bin_hisat2, nthreads, genome, infile, outfile)
        elif "rna" in mode:
            command = r"""{} -p {} -x {}  -k 100 --no-softclip --known-splicesite-infile {} --dta-cufflinks """ \
                      r"""--novel-splicesite-outfile {} -U {} > {}""".format(bin_hisat2, nthreads, genome, known_splice,
                                                                             novel_splice, infile, outfile)
        else:
            raise Exception("Mode unknown!")

        is_err, cmd_runtime = self.run_and_check_command(command)

        logging.debug("{} alignment finished! Runtime: {}.".format(mode, cmd_runtime))

        return 0


    def aligned_to_bed(self, infile, outfile, mode):

        ### 2 mismatches, up to 1 alignment
        command = r"""samtools view -Sh -F 4 {} | grep -E 'XM:i:[0-2]\s.*NH:i:1$|^@' | samtools view -Sbh - | bedtools bamtobed -cigar -i stdin > {}""".format(infile, outfile)

        is_err, cmd_runtime = self.run_and_check_command(command)

        logging.debug("{} sam file parser finished! Runtime: {}.".format(mode, cmd_runtime))

        return 0


    def reverse_complement_fastq(self, infile, outfile):

        command = str(r"""paste <(sed -n '1~4p' {0}) <(sed -n '2~4p' {0} | rev | tr "ATGC" "TACG") <(sed -n '3~4p' {0}) <(sed -n '4~4p' {0} | rev) | tr "\t" "\n" > {1}""").format(infile, outfile)

        logging.info("Running reverse complement with command:\n{}".format(command))

        is_err, cmd_runtime = self.run_and_check_command(command)

        logging.debug("Running custom command finished! Runtime: {}.".format(cmd_runtime))

        return 0
