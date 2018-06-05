#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
cd /data2/bio/sandbox/nbee2
sudo chmod -R 777 /data2/bio/sandbox/nbee2

head -n 5 /data1/bio/projects/dsafina/hp_checkpoints/srr_hp_checkpoints.sampledata > /data2/bio/sandbox/nbee2/test.sampledata
tail -n 5 /data1/bio/projects/dsafina/hp_checkpoints/srr_hp_checkpoints.sampledata >> /data2/bio/sandbox/nbee2/test2.sampledata


docker pull ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it ivasilyev/bwt_filtering_pipeline_worker python3 \
/data2/bio/sandbox/nbee2/nBee2.py -i /data2/bio/sandbox/nbee2/test.sampledata -r /data/reference/TADB/index/tadb_v2.0.refdata -m test -o /data2/bio/sandbox/nbee2/test

"""

import argparse
import pandas as pd
import re
import os
import subprocess
from nBee import CoverageExtractor, Utilities, RefDataParser, RefDataLine


class Initializer:
    def __init__(self):
        self._namespace = self.parse_args()
        # Aligner class parameters mirroring
        self.mapped_file_name = self._namespace.input
        self.refdata_file_name = self._namespace.refdata
        self.chunk_number = self._namespace.chunk
        self.non_zero_bool = self._namespace.non_zero
        self.mapped_reads_directory = Utilities.ends_with_slash(os.path.dirname(os.path.abspath(self.mapped_file_name)))
        self._output_directory = Utilities.ends_with_slash("/".join(os.path.dirname(os.path.abspath(self.mapped_file_name)).split("/")[:-1]))
        self.logs_directory = "{}Logs/".format(self._output_directory)
        self.statistics_directory = "{}Statistics/".format(self._output_directory)
        self.create_dirs()
    @staticmethod
    def parse_args():
        starting_parser = argparse.ArgumentParser(description="""
        This script performs coverage extraction from SAM or BAM (faster) file
        Required software: bowtie, bowtie2, samtools, bedtools
        '../Statistics' directory also would be created based to input file path""")
        starting_parser.add_argument("-i", "--input", required=True,
                                     help="SAM or BAM file")
        starting_parser.add_argument("-r", "--refdata", required=True,
                                     help="Linker file generated by the 'cook_the_reference.py' script. Contains JSON or tab-delimited columns: original FASTA file name, bowtie colorspace index, bowtie2 index, FASTA samtools index, genome lengths file, annotation")
        starting_parser.add_argument("-c", "--chunk", default=0, type=int,
                                     help="(Optional) REFDATA chunk number (zero based), 0 by default")
        starting_parser.add_argument("-n", "--non_zero", action='store_true', default=False,
                                     help="(Optional) If specified, only genomes with non-zero coverage will be reported")
        return starting_parser.parse_args()
    def create_dirs(self):
        tmp = [os.makedirs(i, exist_ok=True) for i in [self.logs_directory, self.statistics_directory]]
        del tmp


class PathsKeeper:
    def __init__(self):
        self._refdata = refDataParser.get_parsed_list()[mainInitializer.chunk_number]
        self.samtools_index_file = self._refdata.samtools_index_file
        self.bedtools_genome_file = self._refdata.bedtools_genome_file
        self.annotation_file = self._refdata.annotation_file
        self.sample_mask = Utilities.filename_only(mainInitializer.mapped_file_name)
        self._mapped_reads_file_mask = ".".join(mainInitializer.mapped_file_name.split(".")[:-1])
        if self._mapped_reads_file_mask.endswith("_sorted"):
            self._mapped_reads_file_mask = re.sub("_sorted$", "", self._mapped_reads_file_mask)
        #
        self.mapped_reads_file_name = "{}.sam".format(self._mapped_reads_file_mask)
        self.samtools_converted_file_name = "{}.bam".format(self._mapped_reads_file_mask)
        self.samtools_sorted_file_name = "{}_sorted.bam".format(self._mapped_reads_file_mask)
        self.samtools_index_file_name = "{}.bai".format(self.samtools_sorted_file_name)
        #
        self._statistics_file_mask = "{a}{b}".format(a=mainInitializer.statistics_directory,
                                                     b=self.sample_mask)
        self.samtools_idxstats_file_name = "{}_idxstats.txt".format(self._statistics_file_mask)
        self.samtools_stats_file_name = "{}_sam_stats.txt".format(self._statistics_file_mask)
        self.bedtools_histogram_file_name = "{}_genomeCoverageBed.txt".format(self._statistics_file_mask)
        self.stacked_coverage_file_name = "{}_pos_bp.txt".format(self._statistics_file_mask)
        self.final_coverage_file_name = "{}_coverage.tsv".format(self._statistics_file_mask)
        #
        self._logs_file_mask = "{a}{b}".format(a=mainInitializer.logs_directory,
                                               b=self.sample_mask)
        #
        self.samtools_converted_log_file_name = "{}_sam2bam.log".format(self._logs_file_mask)
        self.samtools_index_log_file_name = "{}_index_bam.log".format(self.samtools_sorted_file_name)


if __name__ == '__main__':
    mainInitializer = Initializer()
    refDataParser = RefDataParser()
    pathsKeeper = PathsKeeper()
    extractor = CoverageExtractor(pathsKeeper)
    extractor.run()
