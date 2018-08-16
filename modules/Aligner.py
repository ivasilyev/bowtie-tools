#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import subprocess
import logging
from modules.PathsKeeper import PathsKeeper
from modules.Utilities import Utilities


class Aligner:
    """
    The class handles the main high-load single-queued pipeline:
    Bowtie (or Bowtie2) -> SamTools
    Consumes: SampleData object
    """
    def __init__(self, path_keeper: PathsKeeper, threads_number: int):
        self._pk = path_keeper
        self._threads_number = threads_number
        Utilities.batch_remove(self._pk.mapped_reads_file_name,
                               self._pk.samtools_converted_file_name,
                               self._pk.samtools_sorted_file_name,
                               self._pk.unmapped_reads_file_name,
                               *self._pk.pairwise_unmapped_reads_files_list,
                               self._pk.aligner_log_file_name)

    def _get_cmd(self):
        """
        The function for detection of colorspace reads which should be processed with Bowtie.
        All other formats are normally processed with Bowtie2.
        """
        if self._pk.raw_reads_file_extension == "csfasta":
            cmd = self.__compose_bowtie_cmds_list()
        else:
            cmd = self.__compose_bowtie2_cmds_list()
        return [str(i) for i in cmd]

    def ___is_large_index(self):
        """
        The function assignes '--large-index' argument to Bowtie command line.
        Bowtie (Bowtie2) large indexes have extension 'ebwtl' ('bt2l') while common indexes end with 'ebwt' ('bt2').
        However, Bowtie2 does not require this flag.
        """
        return any(i.strip().endswith(".ebwtl") for i in subprocess.getoutput("ls -d {}*".format(self._pk.bowtie_index_mask)).split("\n"))

    def __compose_bowtie_cmds_list(self):
        if len(self._pk.raw_reads_files_list) > 1:
            logging.warning("More than one colorspace reads files per sample: '{}'. Only the first file would be processed".format("', '".join(self._pk.raw_reads_files_list)))
        raw_reads_file = self._pk.raw_reads_files_list[0]
        index_mask = self._pk.bowtie_index_mask
        cmd = ['bowtie', '-f', '-C', '-t', '-v', '3', '-k', '1', '--threads', self._threads_number]
        if self.___is_large_index():
            cmd.append('--large-index')
        cmd.extend(['--un', self._pk.unmapped_reads_file_name, index_mask, raw_reads_file, '-S', self._pk.mapped_reads_file_name])
        """
        Bowtie manual page: https://github.com/BenLangmead/bowtie/blob/master/MANUAL.markdown
        -f: The query input files (specified either as <m1> and <m2>, or as <s>) are FASTA files (usually having extension .fa, .mfa, .fna or similar). All quality values are assumed to be 40 on the Phred quality scale.
        -C: Align in colorspace. Read characters are interpreted as colors. The index specified must be a colorspace index (i.e. built with bowtie-build -C, or bowtie will print an error message and quit. See Colorspace alignment for more details.
        -S: Print alignments in SAM format. See the SAM output section of the manual for details. To suppress all SAM headers, use --sam-nohead in addition to -S/--sam. To suppress just the @SQ headers (e.g. if the alignment is against a very large number of reference sequences), use --sam-nosq in addition to -S/--sam. bowtie does not write BAM files directly, but SAM output can be converted to BAM on the fly by piping bowtie's output to samtools view.
        -t: Print the amount of wall-clock time taken by each phase.
        -v 3: Report alignments with at most <int> mismatches. -e and -l options are ignored and quality values have no effect on what alignments are valid. -v is mutually exclusive with -n.
        -k 1: Report up to <int> valid alignments per read or pair (default: 1). Validity of alignments is determined by the alignment policy (combined effects of -n, -v, -l, and -e). If more than one valid alignment exists and the --best and --strata options are specified, then only those alignments belonging to the best alignment "stratum" will be reported. Bowtie is designed to be very fast for small -k but bowtie can become significantly slower as -k increases. If you would like to use Bowtie for larger values of -k, consider building an index with a denser suffix-array sample, i.e. specify a smaller -o/--offrate when invoking bowtie-build for the relevant index (see the Performance tuning section for details).
        --large-index: Force usage of a 'large' index (those ending in '.ebwtl'), even if a small one is present. Default: off.
        --threads <>: Launch <int> parallel search threads (default: 1). Threads will run on separate processors/cores and synchronize when parsing reads and outputting alignments. Searching for alignments is highly parallel, and speedup is fairly close to linear. This option is only available if bowtie is linked with the pthreads library (i.e. if BOWTIE_PTHREADS=0 is not specified at build time).
        --un <>: Write all reads that could not be aligned to a file with name <filename>. Written reads will appear as they did in the input, without any of the trimming or translation of quality values that may have taken place within Bowtie. Paired-end reads will be written to two parallel files with _1 and _2 inserted in the filename, e.g., if <filename> is unaligned.fq, the #1 and #2 mates that fail to align will be written to unaligned_1.fq and unaligned_2.fq respectively. Unless --max is also specified, reads with a number of valid alignments exceeding the limit set with the -m option are also written to <filename>.
        """
        return cmd

    def __compose_bowtie2_cmds_list(self):
        cmd = ["bowtie2", "--very-sensitive", "--threads", self._threads_number]
        if len(self._pk.raw_reads_files_list) == 1:
            if any([self._pk.raw_reads_file_extension == i for i in ["gz", "zip"]]):
                cmd.append("--un-gz")
            else:
                cmd.append("--un")
            cmd.extend([self._pk.unmapped_reads_file_name, "-x", self._pk.bowtie2_index_mask, self._pk.raw_reads_files_list[0]])
        else:
            if any([self._pk.raw_reads_file_extension == i for i in ["gz", "zip"]]):
                cmd.append("--un-conc-gz")
            else:
                cmd.append("--un-conc")
            cmd.extend([self._pk.unmapped_reads_file_name, "-x", self._pk.bowtie2_index_mask, "-1", self._pk.raw_reads_files_list[0], "-2", self._pk.raw_reads_files_list[1]])
        cmd.extend(["-S", self._pk.mapped_reads_file_name])
        """
        Bowtie2 manual page: https://github.com/BenLangmead/bowtie2/blob/master/MANUAL.markdown
        --un-conc-gz <>: Write paired-end reads that fail to align concordantly to file(s) at <path>. These reads correspond to the SAM records with the FLAGS 0x4 bit set and either the 0x40 or 0x80 bit set (depending on whether it's mate #1 or #2). .1 and .2 strings are added to the filename to distinguish which file contains mate #1 and mate #2. If a percent symbol, %, is used in <path>, the percent symbol is replaced with 1 or 2 to make the per-mate filenames. Otherwise, .1 or .2 are added before the final dot in <path> to make the per-mate filenames. Reads written in this way will appear exactly as they did in the input files, without any modification (same sequence, same name, same quality string, same quality encoding). Reads will not necessarily appear in the same order as they did in the inputs.
        --un <>: Write unpaired reads that fail to align to file at <path>. These reads correspond to the SAM records with the FLAGS 0x4 bit set and neither the 0x40 nor 0x80 bits set. If --un-gz is specified, output will be gzip compressed. If --un-bz2 or --un-lz4 is specified, output will be bzip2 or lz4 compressed. Reads written in this way will appear exactly as they did in the input file, without any modification (same sequence, same name, same quality string, same quality encoding). Reads will not necessarily appear in the same order as they did in the input.
        -x <>: The basename of the index for the reference genome. The basename is the name of any of the index files up to but not including the final .1.bt2 / .rev.1.bt2 / etc. bowtie2 looks for the specified index first in the current directory, then in the directory specified in the BOWTIE2_INDEXES environment variable.
        -f: Reads (specified with <m1>, <m2>, <s>) are FASTA files. FASTA files usually have extension .fa, .fasta, .mfa, .fna or similar. FASTA files do not have a way of specifying quality values, so when -f is set, the result is as if --ignore-quals is also set.
        -1 <>: Comma-separated list of files containing mate 1s (filename usually includes _1), e.g. -1 flyA_1.fq,flyB_1.fq. Sequences specified with this option must correspond file-for-file and read-for-read with those specified in <m2>. Reads may be a mix of different lengths. If - is specified, bowtie2 will read the mate 1s from the "standard in" or "stdin" filehandle.
        -2 <>: Comma-separated list of files containing mate 2s (filename usually includes _2), e.g. -2 flyA_2.fq,flyB_2.fq. Sequences specified with this option must correspond file-for-file and read-for-read with those specified in <m1>. Reads may be a mix of different lengths. If - is specified, bowtie2 will read the mate 2s from the "standard in" or "stdin" filehandle.
        --threads <>: Launch NTHREADS parallel search threads (default: 1). Threads will run on separate processors/cores and synchronize when parsing reads and outputting alignments. Searching for alignments is highly parallel, and speedup is close to linear. Increasing -p increases Bowtie 2's memory footprint. E.g. when aligning to a human genome index, increasing -p from 1 to 8 increases the memory footprint by a few hundred megabytes. This option is only available if bowtie is linked with the pthreads library (i.e. if BOWTIE_PTHREADS=0 is not specified at build time).
        --very-sensitive: Same as: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
        -S: File to write SAM alignments to. By default, alignments are written to the "standard out" or "stdout" filehandle (i.e. the console).
        """
        return cmd

    def run(self):
        bwt_cmd_string = " ".join(self._get_cmd())
        logging.info("Alignment launch command line: '{}'".format(bwt_cmd_string))
        log = subprocess.getoutput(bwt_cmd_string)
        Utilities.dump_string(string=log, file=self._pk.aligner_log_file_name)
