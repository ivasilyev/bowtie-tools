#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import os
import gc
import subprocess
import logging
import pandas as pd
from modules.Utilities import Utilities
from modules.PathsKeeper import PathsKeeper


class CoverageExtractor:
    def __init__(self, paths_keeper: PathsKeeper, non_zero_bool: bool = False):
        self._index_column = "reference_id"
        self._non_zero_bool = non_zero_bool
        self._pk = paths_keeper
        self._samtools_idxstats_df = pd.DataFrame()
        self._samtools_stats_dict = {}
        self._bedtools_histogram_2d_array = []
        self._stacked_coverages_df = pd.DataFrame()

    def _sam2bam2sorted_bam(self):
        subprocess.getoutput("rm -f {}*".format(self._pk.samtools_sorted_file_name))
        Utilities.batch_remove(self._pk.samtools_converted_log_file_name)
        # SamTools details: http://www.htslib.org/doc/samtools.html
        # Avoiding self._pk.samtools_converted_file_name
        s = subprocess.getoutput("samtools view -bu -@ 1 {a} | \
                                  samtools sort - -o -@ 1 {b}".format(a=self._pk.mapped_reads_file_name,
                                                                      b=self._pk.samtools_sorted_file_name))
        Utilities.dump_string(string=s, file=self._pk.samtools_converted_log_file_name)
        logging.info("Sorted SAM file: '{}'".format(self._pk.samtools_sorted_file_name))
        del s

    def _index_bam(self):
        Utilities.batch_remove(self._pk.samtools_index_file_name, self._pk.samtools_index_log_file_name)
        s = subprocess.getoutput("samtools index {}".format(self._pk.samtools_sorted_file_name))
        Utilities.dump_string(string=s, file=self._pk.samtools_index_log_file_name)
        logging.info("Indexed BAM file: '{}'".format(self._pk.samtools_index_file_name))
        del s

    def _bam2idxstats(self):
        Utilities.batch_remove(self._pk.samtools_idxstats_file_name, self._pk.samtools_idxstats_log_file_name)
        s = subprocess.getoutput("samtools idxstats {a} 2> {b}".format(a=self._pk.samtools_sorted_file_name, b=self._pk.samtools_idxstats_log_file_name))
        Utilities.dump_string(string=s, file=self._pk.samtools_idxstats_file_name)
        logging.info("Saved SAMTools mapped reads statistics: '{}'".format(self._pk.samtools_idxstats_file_name))
        self._samtools_idxstats_df = pd.DataFrame(Utilities.string_to_2d_array(s), columns=[self._index_column,
                                                                                            "id_bp",
                                                                                            "id_mapped_reads",
                                                                                            "id_unmapped_reads"])
        del s

    def _bam2stats(self):
        def __get_base_alignment_stats(string: str):
            d = {}
            # SamTools stats file columns: ID, stat, value, comment
            for line_list in Utilities.string_to_2d_array(string):
                if len(line_list) < 3 or line_list[0] != "SN":
                    continue
                d[re.sub(":$", "", line_list[1])] = line_list[2]
            if len(d) == 0:
                logging.critical("Bad alignment: no SAMTools stats to extract!")
                return {}
            try:
                out = {"total_reads": d["raw total sequences"],
                       "mapped_reads": d["reads mapped"],
                       "total_bp": d["total length"],
                       "mapped_bp": d["bases mapped"]}
            except KeyError:
                return {}
            return {"sample_{}".format(k): int(out[k]) for k in out}

        Utilities.batch_remove(self._pk.samtools_stats_file_name, self._pk.samtools_stats_log_file_name)
        s = subprocess.getoutput("samtools stats {a} 2> {b}".format(a=self._pk.samtools_sorted_file_name, b=self._pk.samtools_stats_log_file_name))
        Utilities.dump_string(string=s, file=self._pk.samtools_stats_file_name)
        logging.info("Saved SAMTools total coverage statistics: '{}'".format(self._pk.samtools_stats_file_name))
        self._samtools_stats_dict = __get_base_alignment_stats(s)
        del s

    def _bam2histogram(self):
        Utilities.batch_remove(self._pk.bedtools_histogram_file_name, self._pk.genomeCoverageBed_log_file_name)
        s = subprocess.getoutput("genomeCoverageBed -ibam {a} 2> {b}".format(a=self._pk.samtools_sorted_file_name, b=self._pk.genomeCoverageBed_log_file_name))
        # GenomeCoverageBed details: https://bedtools.readthedocs.io/en/stable/content/tools/genomecov.html
        # Cannot be converted to DataFrame before stacking
        Utilities.dump_string(string=s, file=self._pk.bedtools_histogram_file_name)
        self._bedtools_histogram_2d_array = Utilities.string_to_2d_array(s)
        if len(self._bedtools_histogram_2d_array) == 0:
            logging.critical("Bad alignment: no BEDTools coverage histogram to save!")
        logging.info("Saved BEDTools coverage histogram data: '{}'".format(self._pk.bedtools_histogram_file_name))
        del s

    def _stack_coverage(self):
        Utilities.batch_remove(self._pk.stacked_coverage_file_name)
        # genomecov file columns: reference sequence name, depth of coverage, breadth of coverage with that depth, sequence length, coverage ratio
        stacked_coverages_2d_array = []
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
                    stacked_coverages_2d_array.append([counting_id,
                                                       id_maximal_coverage_depth,
                                                       id_coverage_breadth,
                                                       id_bp,
                                                       id_coverage_breadth_to_id_bp,
                                                       id_mapped_bp])
                row_processing_2d_array = []
                counting_id = reference_id
        if len(stacked_coverages_2d_array) == 0:
            logging.critical("Bad alignment: no coverage to stack!")
            return
        self._stacked_coverages_df = pd.DataFrame(stacked_coverages_2d_array, columns=[self._index_column,
                                                                                       "id_maximal_coverage_depth",
                                                                                       "id_coverage_breadth",
                                                                                       "id_bp",
                                                                                       "id_coverage_breadth_to_id_bp",
                                                                                       "id_mapped_bp"])
        self._stacked_coverages_df.to_csv(self._pk.stacked_coverage_file_name, sep='\t', index=False)
        logging.info("Stacked BEDTools coverage: '{}'".format(self._pk.stacked_coverage_file_name))
        del self._bedtools_histogram_2d_array, stacked_coverages_2d_array
        gc.collect()

    def _reference2statistics(self):
        Utilities.batch_remove(self._pk.final_coverage_file_name)
        stats_dict = self._samtools_stats_dict
        if len(stats_dict) == 0:
            logging.critical("Bad alignment: empty SAMTools stats: '{}'".format(self._pk.samtools_stats_file_name))
            return
        if len(self._stacked_coverages_df) == 0:
            logging.critical("Bad alignment: empty stacked BEDTools coverage: '{}'".format(self._pk.stacked_coverage_file_name))
            return
        chunk_size = 10 ** 6
        reader = pd.read_table(self._pk.bedtools_genome_file, sep='\t', header="infer", names=[self._index_column, "id_bp"], chunksize=chunk_size)
        for chunk_number, reference_df in enumerate(reader):
            genomes_coverages_df = reference_df.merge(self._stacked_coverages_df.loc[:, [self._index_column] + [i for i in list(self._stacked_coverages_df) if i not in list(reference_df)]], on=self._index_column, how="left")
            genomes_coverages_df = genomes_coverages_df[~genomes_coverages_df[self._index_column].isin(["*", "genome"])]
            if self._non_zero_bool:
                genomes_coverages_df = genomes_coverages_df[genomes_coverages_df.id_coverage_breadth.notnull()]
            else:
                genomes_coverages_df = genomes_coverages_df.fillna(0)
            genomes_coverages_df["id_total_relative_abundance"] = (10 ** 12) * genomes_coverages_df["id_mapped_bp"].astype(int) / (genomes_coverages_df["id_bp"].astype(int) * int(stats_dict["sample_total_bp"]))
            genomes_coverages_df["id_mapped_relative_abundance"] = (10 ** 12) * genomes_coverages_df["id_mapped_bp"].astype(int) / (genomes_coverages_df["id_bp"].astype(int) * int(stats_dict["sample_mapped_bp"]))
            # MRA details: http://www.ibmc.msk.ru/content/thesisDocs/TyakhtAV_thesis.pdf (p.63)
            genomes_coverages_df["sample_total_reads"] = stats_dict["sample_total_reads"]
            genomes_coverages_df["sample_mapped_reads"] = stats_dict["sample_mapped_reads"]
            genomes_coverages_df["sample_total_bp"] = stats_dict["sample_total_bp"]
            genomes_coverages_df["sample_mapped_bp"] = stats_dict["sample_mapped_bp"]
            genomes_coverages_df["sample_average_total_reads_bp"] = float(stats_dict["sample_total_reads"]) / float(stats_dict["sample_total_bp"])
            genomes_coverages_df["sample_average_mapped_reads_bp"] = float(stats_dict["sample_mapped_reads"]) / float(stats_dict["sample_total_bp"])
            genomes_coverages_df["sample_mapped_reads_to_total_reads"] = float(stats_dict["sample_mapped_reads"]) / float(stats_dict["sample_total_reads"])
            genomes_coverages_df = genomes_coverages_df.merge(self._samtools_idxstats_df.loc[:, [self._index_column] + [i for i in list(self._samtools_idxstats_df) if i not in list(genomes_coverages_df)]], on=self._index_column, how="left")
            genomes_coverages_df["id_mapped_reads_per_million_sample_total_reads"] = genomes_coverages_df["id_mapped_reads"].astype(int) * (10 ** 6) / int(stats_dict["sample_total_reads"])
            genomes_coverages_df["id_mapped_reads_per_million_sample_mapped_reads"] = genomes_coverages_df["id_mapped_reads"].astype(int) * (10 ** 6) / int(stats_dict["sample_mapped_reads"])
            # RPM details: https://www.biostars.org/p/273537/
            genomes_coverages_df["id_mapped_reads_per_kbp_per_million_sample_total_reads"] = genomes_coverages_df["id_mapped_reads"].astype(int) * (10 ** 9) / (int(stats_dict["sample_total_reads"]) * genomes_coverages_df["id_bp"])
            genomes_coverages_df["id_mapped_reads_per_kbp_per_million_sample_mapped_reads"] = genomes_coverages_df["id_mapped_reads"].astype(int) * (10 ** 9) / (int(stats_dict["sample_mapped_reads"]) * genomes_coverages_df["id_bp"])
            # RPKM details: https://www.biostars.org/p/273537/
            for int_column in ["id_bp", "id_coverage_breadth", "id_mapped_bp", "id_maximal_coverage_depth",
                               "id_mapped_reads", "sample_total_reads", "sample_mapped_reads", "sample_total_bp",
                               "sample_mapped_bp"]:
                genomes_coverages_df[int_column] = genomes_coverages_df[int_column].astype(int)
            if chunk_number == 0:
                genomes_coverages_df.to_csv(self._pk.final_coverage_file_name, sep='\t', header=True, index=False)
            else:
                with open(file=self._pk.final_coverage_file_name, mode="a", encoding="utf-8") as f:
                    genomes_coverages_df.to_csv(f, sep='\t', header=False, index=False)
                logging.info("Processed chunk {} with size of {} lines".format(chunk_number, chunk_size))
        logging.info("Finished processing coverage table: '{}'".format(self._pk.final_coverage_file_name))

    def run(self):
        functions_list = [self._sam2bam2sorted_bam,
                          self._index_bam,
                          self._bam2idxstats,
                          self._bam2stats,
                          self._bam2histogram,
                          self._stack_coverage,
                          self._reference2statistics]
        if os.path.isfile(self._pk.samtools_sorted_file_name):
            functions_list = functions_list[1:]
        elif not (os.path.isfile(self._pk.samtools_converted_file_name) or os.path.isfile(self._pk.mapped_reads_file_name)):
            logging.critical("Not recognized input file: '{}'".format(self._pk.mapped_reads_file_name))
        for f in functions_list:
            f()
