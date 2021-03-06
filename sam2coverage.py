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
import re
import os
from modules.Utilities import Utilities
from modules.SampleDataLine import SampleDataLine
from modules.RefDataLine import RefDataLine


class Initializer:
    def __init__(self):
        namespace = self._parse_args()
        self.input_file_name = namespace.input
        self.refdata_file_name = namespace.refdata
        self.chunk_id = namespace.chunk
        self.mapped_reads_directory = Utilities.ends_with_slash(os.path.dirname(os.path.abspath(self.input_file_name)))

        self._output_directory = Utilities.ends_with_slash("/".join(os.path.dirname(os.path.abspath(self.input_file_name)).split("/")[:-1]))
        self.logs_directory = "{}Logs/".format(self._output_directory)
        self.statistics_directory = "{}Statistics/".format(self._output_directory)
        self._create_dirs()

    @staticmethod
    def _parse_args():
        starting_parser = argparse.ArgumentParser(description="""
This script performs coverage extraction from SAM or BAM (faster) file
Required software: bowtie, bowtie2, samtools, bedtools
'../Statistics' directory also would be created based to input file path
""")
        starting_parser.add_argument("-i", "--input", required=True,
                                     help="SAM or BAM file")
        starting_parser.add_argument("-r", "--refdata", required=True,
                                     help="Linker file generated by the 'cook_the_reference.py' script. Contains JSON or tab-delimited columns: original FASTA file name, bowtie colorspace index, bowtie2 index, FASTA samtools index, genome lengths file, annotation")
        starting_parser.add_argument("-c", "--chunk",
                                     help="(Optional) Sequence ID if multiple references are provided")
        starting_parser.add_argument("-n", "--non_zero", action='store_true', default=False,
                                     help="(Optional) If specified, only genomes with non-zero coverage will be reported")
        return starting_parser.parse_args()

    def _create_dirs(self):
        tmp = [os.makedirs(i, exist_ok=True) for i in [self.logs_directory, self.statistics_directory]]
        del tmp


class PathsKeeper:
    def __init__(self, refdata: RefDataLine):
        # Output directories
        output_dir = Utilities.ends_with_slash(os.path.dirname(os.path.realpath(mainInitializer.input_file_name)))
        mapped_reads_directory = "{}Mapped_reads/".format(output_dir)
        statistics_directory = "{}Statistics/".format(output_dir)
        logs_directory = "{}Logs/".format(output_dir)
        for path in [mapped_reads_directory, statistics_directory, logs_directory]:
            os.makedirs(path, exist_ok=True)
        # Reference data
        self.samtools_index_file = refdata.samtools_index_file
        self.bedtools_genome_file = refdata.bedtools_genome_file
        self.annotation_file = refdata.annotation_file
        # Output files
        if mainInitializer.input_file_name.endswith("_sorted.bam"):
            sample_name = re.sub("_sorted.bam$", "", mainInitializer.input_file_name)
        else:
            sample_name = ".".join(mainInitializer.input_file_name.split(".")[:-1])
        mapped_reads_file_mask = "{a}{b}".format(a=mapped_reads_directory, b=sample_name)
        self.mapped_reads_file_name = "{}.sam".format(mapped_reads_file_mask)
        self.samtools_converted_file_name = "{}.bam".format(mapped_reads_file_mask)
        self.samtools_sorted_file_name = "{}_sorted.bam".format(mapped_reads_file_mask)

        self._mapped_output_mask = refdata.db_name
        unmapped_output_mask = "_".join(["no", self._mapped_output_mask])
        # Sample data
        sample_name = sampledata.name
        self.raw_reads_files_list = sampledata.raw_reads_files_list
        self.raw_reads_file_extension = self.raw_reads_files_list[0].split(".")[-1]
        unmapped_reads_file_mask = "{a}{b}_{c}".format(a=unmapped_reads_directory, b=sample_name, c=unmapped_output_mask)
        self.unmapped_reads_file_name = "{a}.{b}".format(a=unmapped_reads_file_mask,
                                                         b=self.raw_reads_file_extension)
        self.pairwise_unmapped_reads_files_list = ["{a}.{i}.{b}".format(a=unmapped_reads_file_mask,
                                                                        i=i,
                                                                        b=self.raw_reads_file_extension) for i in [1, 2]]
        mapped_reads_file_mask = "{a}{b}_{c}".format(a=mapped_output_mask, b=sample_name, c=self._mapped_output_mask)
        self.mapped_reads_file_name = "{}.sam".format(mapped_reads_file_mask)
        self.samtools_converted_file_name = "{}.bam".format(mapped_reads_file_mask)
        self.samtools_sorted_file_name = "{}_sorted.bam".format(mapped_reads_file_mask)
        self.samtools_index_file_name = "{}.bai".format(self.samtools_sorted_file_name)
        statistics_file_mask = "{a}{b}_{c}".format(a=statistics_directory, b=sample_name, c=self._mapped_output_mask)
        self.samtools_idxstats_file_name = "{}_idxstats.txt".format(statistics_file_mask)
        self.samtools_stats_file_name = "{}_sam_stats.txt".format(statistics_file_mask)
        self.bedtools_histogram_file_name = "{}_genomeCoverageBed.txt".format(statistics_file_mask)
        self.stacked_coverage_file_name = "{}_pos_bp.txt".format(statistics_file_mask)
        self.final_coverage_file_name = "{}_coverage.tsv".format(statistics_file_mask)
        logs_file_mask = "{a}{b}_{c}".format(a=logs_directory, b=sample_name, c=self._mapped_output_mask)
        self.aligner_log_file_name = "{}_aligner.log".format(logs_file_mask)
        self.samtools_converted_log_file_name = "{}_sam2bam.log".format(logs_file_mask)
        self.samtools_index_log_file_name = "{}_index_bam.log".format(logs_file_mask)
        


if __name__ == '__main__':
    mainInitializer = Initializer()
    refDataParser = RefDataParser(mainInitializer.refdata_file_name)
    pathsKeeper = PathsKeeper()
    extractor = CoverageExtractor(pathsKeeper)
    extractor.run()
