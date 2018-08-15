#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
from modules.SampleDataLine import SampleDataLine
from modules.RefDataLine import RefDataLine
from modules.Utilities import Utilities


class PathsKeeper:
    def __init__(self, sampledata: SampleDataLine, refdata: RefDataLine, output_dir: str):
        # Output directories
        output_dir = Utilities.ends_with_slash(output_dir)
        unmapped_reads_directory = "{}Unmapped_reads/".format(output_dir)
        mapped_reads_directory = "{}Mapped_reads/".format(output_dir)
        statistics_directory = "{}Statistics/".format(output_dir)
        logs_directory = "{}Logs/".format(output_dir)
        for path in [unmapped_reads_directory, mapped_reads_directory, statistics_directory, logs_directory]:
            os.makedirs(path, exist_ok=True)
        # Reference data
        self.bowtie_index_mask = refdata.bowtie_index_mask
        self.bowtie2_index_mask = refdata.bowtie2_index_mask
        self.samtools_index_file = refdata.samtools_index_file
        self.bedtools_genome_file = refdata.bedtools_genome_file
        self.annotation_file = refdata.annotation_file
        # Output masks
        mapped_output_mask = refdata.db_name
        unmapped_output_mask = "_".join(["no", mapped_output_mask])
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
        mapped_reads_file_mask = "{a}{b}_{c}".format(a=mapped_reads_directory, b=sample_name, c=mapped_output_mask)
        self.mapped_reads_file_name = "{}.sam".format(mapped_reads_file_mask)
        self.samtools_converted_file_name = "{}.bam".format(mapped_reads_file_mask)
        self.samtools_sorted_file_name = "{}_sorted.bam".format(mapped_reads_file_mask)
        self.samtools_index_file_name = "{}.bai".format(self.samtools_sorted_file_name)
        statistics_file_mask = "{a}{b}_{c}".format(a=statistics_directory, b=sample_name, c=mapped_output_mask)
        self.samtools_idxstats_file_name = "{}_idxstats.txt".format(statistics_file_mask)
        self.samtools_stats_file_name = "{}_sam_stats.txt".format(statistics_file_mask)
        self.bedtools_histogram_file_name = "{}_genomeCoverageBed.txt".format(statistics_file_mask)
        self.stacked_coverage_file_name = "{}_pos_bp.txt".format(statistics_file_mask)
        self.final_coverage_file_name = "{}_coverage.tsv".format(statistics_file_mask)
        logs_file_mask = "{a}{b}_{c}".format(a=logs_directory, b=sample_name, c=mapped_output_mask)
        self.aligner_log_file_name = "{}_aligner.log".format(logs_file_mask)
        self.samtools_converted_log_file_name = "{}_sam2bam.log".format(logs_file_mask)
        self.samtools_index_log_file_name = "{}_index_bam.log".format(logs_file_mask)




