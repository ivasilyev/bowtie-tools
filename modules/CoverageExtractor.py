#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import re
import os
import subprocess
import logging
import pandas as pd
from modules.Utilities import Utilities
from modules.PathsKeeper import PathsKeeper


class CoverageExtractor:
    def __init__(self, paths_keeper: PathsKeeper, non_zero_bool: bool = False):
        self._non_zero_bool = non_zero_bool
        self._pk = paths_keeper
        self._samtools_idxstats = ""
        self._samtools_stats = ""
        self._bedtools_histogram_2d_array = []
        self._stacked_coverages_2d_array = []


    def _sam2bam2sorted_bam(self):
        Utilities.batch_remove(self._pk.samtools_sorted_file_name, self._pk.samtools_converted_log_file_name)
        # Avoiding self._pk.samtools_converted_file_name
        s = subprocess.getoutput("samtools view -bu -@ 1 {a} | \
                                  samtools sort - -o -@ 1 {b}".format(a=self._pk.mapped_reads_file_name,
                                                                      b=self._pk.samtools_sorted_file_name))
        Utilities.dump_string(string=s, file=self._pk.samtools_converted_log_file_name)
        logging.info("Sorted SAM file: '{}'".format(self._pk.samtools_sorted_file_name))

    def _index_bam(self):
        Utilities.batch_remove(self._pk.samtools_index_file_name, self._pk.samtools_index_log_file_name)
        s = subprocess.getoutput("samtools index {}".format(self._pk.samtools_sorted_file_name))
        Utilities.dump_string(string=s, file=self._pk.samtools_index_log_file_name)
        logging.info("Indexed BAM file: '{}'".format(self._pk.samtools_index_file_name))

    def _bam2idxstats(self):
        Utilities.batch_remove(self._pk.samtools_idxstats_file_name)
        self._samtools_idxstats = subprocess.getoutput("samtools idxstats {}".format(self._pk.samtools_sorted_file_name))
        Utilities.dump_string(string=self._samtools_idxstats, file=self._pk.samtools_idxstats_file_name)
        logging.info("Mapped reads statistics: '{}'".format(self._pk.samtools_idxstats_file_name))

    def _bam2stats(self):
        Utilities.batch_remove(self._pk.samtools_stats_file_name)
        self._samtools_stats = subprocess.getoutput("samtools stats {}".format(self._pk.samtools_sorted_file_name))
        Utilities.dump_string(string=self._samtools_stats, file=self._pk.samtools_stats_file_name)
        logging.info("Total coverage statistics: '{}'".format(self._pk.samtools_stats_file_name))

    def _bam2histogram(self):
        Utilities.batch_remove(self._pk.bedtools_histogram_file_name)
        s = subprocess.getoutput("genomeCoverageBed -ibam {}".format(self._pk.samtools_sorted_file_name))
        for row_string in Utilities.remove_empty_values([i.split("\n") for i in re.sub("[\r\n]+", "\n", s)]):
            self._bedtools_histogram_2d_array.append(Utilities.remove_empty_values([i.strip() for i in row_string.split("\t")]))
        if len(self._bedtools_histogram_2d_array) == 0:
            Utilities.log_and_raise("Bad alignment: no histogram to save!")
        Utilities.dump_string(string=s, file=self._pk.bedtools_histogram_file_name)
        logging.info("Coverage histogram data: '{}'".format(self._pk.bedtools_histogram_file_name))

    def _stack_coverage(self):
        Utilities.batch_remove(self._pk.stacked_coverage_file_name)
        self._stacked_coverages_2d_array = [["reference_id", "id_bp", "id_coverage_breadth", "id_mapped_bp", "id_coverage_breadth_to_id_bp", "id_maximal_coverage_depth", "id_total_relative_abundance", "id_mapped_relative_abundance", "id_mapped_reads", "sample_total_reads", "sample_mapped_reads", "sample_total_bp", "sample_mapped_bp", "sample_average_total_reads_bp", "sample_average_mapped_reads_bp", "sample_mapped_reads_to_total_reads", "id_mapped_reads_per_million_sample_total_reads", "id_mapped_reads_per_million_sample_mapped_reads"]]
        # genomecov file columns: reference sequence name, depth of coverage, breadth of coverage with that depth, sequence length, coverage ratio
        row_processing_2d_array = []
        counting_id = ""
        for row_list in self._bedtools_histogram_2d_array:
            if len(row_list) != 5:
                logging.warning("Cannot parse coverage histogram row '{a}' from file '{b}'".format(a=row_list, b=self._pk.bedtools_histogram_file_name))
                continue
            reference_id, id_local_coverage_depth, id_local_coverage_breadth, id_bp, id_local_coverage_ratio = row_list
            if reference_id == 'genome' or '*' in reference_id:
                continue
            if reference_id == counting_id and int(id_local_coverage_depth) > 0:
                row_processing_2d_array.append(row_list)
            else:
                if len(row_processing_2d_array) > 0:
                    # output file columns: reference sequence name, maximal depth of coverage, total breadth of coverage, sequence length, coverage ratio, total mapped bases
                    id_maximal_coverage_depth = max([int(i[1]) for i in row_processing_2d_array])
                    id_coverage_breadth = sum([int(i[2]) for i in row_processing_2d_array])
                    id_bp = int(row_processing_2d_array[0][3])
                    id_coverage_breadth_to_id_bp = sum([float(i[4]) for i in row_processing_2d_array])
                    id_mapped_bp = sum([int(i[1]) * int(i[2]) for i in row_processing_2d_array])
                    self._stacked_coverages_2d_array.append([counting_id, id_maximal_coverage_depth, id_coverage_breadth, id_bp, id_coverage_breadth_to_id_bp, id_mapped_bp])
                row_processing_2d_array = []
                counting_id = reference_id
        if len(self._stacked_coverages_2d_array) == 0:
            Utilities.log_and_raise("Bad alignment: no coverage to stack!")
        Utilities.dump_2d_array(array=self._stacked_coverages_2d_array, file=self._pk.stacked_coverage_file_name)
        logging.info("Stacked coverage: '{}'".format(self._pk.stacked_coverage_file_name))

    def __get_base_alignment_stats(self):
        s = subprocess.getoutput("grep ^SN {} | cut -f 2-".format(self._pk.samtools_stats_file_name))
        d = {"total_reads": re.findall('raw total sequences:\t(\d*)', s)[0],
             "mapped_reads": re.findall('reads mapped:\t(\d*)', s)[0],
             "total_bp": re.findall('total length:\t(\d*)', s)[0],
             "mapped_bp": re.findall('bases mapped:\t(\d*)', s)[0]}
        d = {"sample_{}".format(k): int(d[k]) for k in d}
        return d

    def _reference2statistics(self):
        output_file = self._pk.final_coverage_file_name
        pos_bp_file = self._pk.stacked_coverage_file_name
        Utilities.batch_remove(output_file)
        stats_dict = self.__get_base_alignment_stats()
        reference_df = pd.read_table(self._pk.bedtools_genome_file, header='infer', sep='\t', names=['reference_id', 'id_bp'])
        stacked_coverages_df = pd.read_table(pos_bp_file, sep='\t', header='infer', names=["reference_id", "id_maximal_coverage_depth", "id_coverage_breadth", "id_bp", "id_coverage_breadth_to_id_bp", "id_mapped_bp"])
        if len(stacked_coverages_df) == 0:
            logging.critical("Bad alignment: empty coverage file '{}'".format(self._pk.stacked_coverage_file_name))
            return
        genomes_coverages_df = pd.merge(reference_df, stacked_coverages_df.loc[:, [column for column in list(stacked_coverages_df) if column != "id_bp"]], on="reference_id", how='outer')
        genomes_coverages_df = genomes_coverages_df.loc[genomes_coverages_df['reference_id'] != 'genome']
        genomes_coverages_df["id_total_relative_abundance"] = genomes_coverages_df.loc[:, "id_mapped_bp"] / (genomes_coverages_df.loc[:, "id_bp"] * stats_dict["sample_total_bp"])
        genomes_coverages_df["id_mapped_relative_abundance"] = genomes_coverages_df.loc[:, "id_mapped_bp"] / (genomes_coverages_df.loc[:, "id_bp"] * stats_dict["sample_mapped_bp"])
        genomes_coverages_df["sample_total_reads"] = stats_dict["sample_total_reads"]
        genomes_coverages_df["sample_mapped_reads"] = stats_dict["sample_mapped_reads"]
        genomes_coverages_df["sample_total_bp"] = stats_dict["sample_total_bp"]
        genomes_coverages_df["sample_mapped_bp"] = stats_dict["sample_mapped_bp"]
        genomes_coverages_df["sample_average_total_reads_bp"] = float(stats_dict["sample_total_reads"]) / float(stats_dict["sample_total_bp"])
        genomes_coverages_df["sample_average_mapped_reads_bp"] = float(stats_dict["sample_mapped_reads"]) / float(stats_dict["sample_total_bp"])
        genomes_coverages_df["sample_mapped_reads_to_total_reads"] = float(stats_dict["sample_mapped_reads"]) / float(stats_dict["sample_total_reads"])
        idxstats_df = pd.read_table(self._pk.samtools_idxstats_file_name, sep='\t', header='infer', names=["reference_id", "id_bp", "id_mapped_reads", "id_unmapped_reads"])
        genomes_coverages_df = pd.merge(genomes_coverages_df, idxstats_df.loc[:, [column for column in list(idxstats_df) if column != "id_bp"]], on="reference_id", how='outer')
        genomes_coverages_df["id_mapped_reads_per_million_sample_total_reads"] = genomes_coverages_df.loc[:, "id_mapped_reads"] * (10 ** 6) / float(stats_dict["sample_total_reads"])
        genomes_coverages_df["id_mapped_reads_per_million_sample_mapped_reads"] = genomes_coverages_df.loc[:, "id_mapped_reads"] * (10 ** 6) / float(stats_dict["sample_mapped_reads"])
        output_df = genomes_coverages_df.loc[:, ["reference_id", "id_bp", "id_coverage_breadth", "id_mapped_bp", "id_coverage_breadth_to_id_bp", "id_maximal_coverage_depth", "id_total_relative_abundance", "id_mapped_relative_abundance", "id_mapped_reads", "sample_total_reads", "sample_mapped_reads", "sample_total_bp", "sample_mapped_bp", "sample_average_total_reads_bp", "sample_average_mapped_reads_bp", "sample_mapped_reads_to_total_reads", "id_mapped_reads_per_million_sample_total_reads", "id_mapped_reads_per_million_sample_mapped_reads"]]
        if self._non_zero_bool:
            output_df = output_df[output_df.id_coverage_breadth.notnull()]
        else:
            output_df = output_df.fillna(0)
        for int_column in ["id_bp", "id_coverage_breadth", "id_mapped_bp", "id_maximal_coverage_depth", "id_mapped_reads", "sample_total_reads", "sample_mapped_reads", "sample_total_bp", "sample_mapped_bp"]:
            output_df[int_column] = output_df.loc[:, int_column].astype(int)
        output_df = output_df.loc[~output_df['reference_id'].str.contains('*', regex=False)]
        output_df.to_csv(output_file, sep='\t', index=False)
        logging.info("Extracted coverage table: '{}'".format(output_file))

    def run(self):
        functions_list = [self._sam2bam2sorted_bam, self._index_bam, self._bam2idxstats, self._bam2stats, self._bam2histogram, self._stack_coverage, self._reference2statistics]
        if os.path.isfile(self._pk.samtools_sorted_file_name):
            functions_list = functions_list[1:]
        elif os.path.isfile(self._pk.samtools_converted_file_name) or os.path.isfile(self._pk.mapped_reads_file_name):
            pass
        else:
            Utilities.log_and_raise("Not recognized input file: '{}'".format(self._pk.mapped_reads_file_name))
        tmp = [i() for i in functions_list]
        del tmp
