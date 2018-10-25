#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
import argparse
import pandas as pd
from modules.Utilities import Utilities


class Initializer:
    def __init__(self):
        self._namespace = self._parse_args()
        self.sampledata = self._namespace.input
        self.target_length = CoveragesVerifier.get_wc_l(self._namespace.genome) + 1
        self.prefix = self._namespace.prefix
        self.suffix = self._namespace.suffix
        self.debugging_bool = self._namespace.debug
        self.output = self._namespace.output
        if len(self.output) == 0:
            self.output = "{}sampledata/{}.sampledata".format(Utilities.ends_with_slash(os.path.dirname(self.prefix)), Utilities.get_time())

    @staticmethod
    def _parse_args():
        starting_parser = argparse.ArgumentParser(description="The script to detect nBee or refdata2schedule pipeline faults."
                                                              "Search example: <prefix><sample_name><suffix>")
        starting_parser.add_argument("-i", "--input", required=True,
                                     help="Input table containing two tab-delimited columns for colorspace or non-colorspace sequences and three for paired-end sequences: sample name and absolute path(s). May contain a header")
        starting_parser.add_argument("-g", "--genome", required=True,
                                     help="BEDTools GenomeCoverageBed genome index file")
        starting_parser.add_argument("-p", "--prefix", required=True,
                                     help="Output files prefix before sample name, e.g '/path/to/coverage/'.")
        starting_parser.add_argument("-s", "--suffix", required=True,
                                     help="Output files suffix after sample name, e.g 'coverage_mask.extension'")
        starting_parser.add_argument("-d", "--debug", default=False, action='store_true',
                                     help="Create debugging table")
        starting_parser.add_argument("-o", "--output", default="",
                                     help="Output file name, by default <prefix_dir>/sampledata/<current_time>.sampledata")
        return starting_parser.parse_args()


class CoveragesVerifier:
    def __init__(self):
        self._queue = Utilities.load_2d_array(mainInitializer.sampledata)
        self._verified_df = pd.DataFrame()
        self._no_coverages_df = pd.DataFrame()

    @staticmethod
    def get_wc_l(file_name: str):
        out = 0
        if os.path.isfile(file_name):
            try:
                out = int(subprocess.getoutput("wc -l < {}".format(file_name)).split('\n')[0])
            except ValueError:
                print("Failed to count lines for file: '{}'".format(file_name))
        return out

    @staticmethod
    def dict2pd_series(d: dict):
        output = pd.Series()
        for key in d:
            output.at[key] = d[key]
        return output

    def process_sampledata_row(self, sampledata_row: list):
        sample_name = sampledata_row[0]
        raw_reads_files_list = sampledata_row[1:]
        coverage_path = "{a}{b}{c}".format(a=mainInitializer.prefix, b=sample_name, c=mainInitializer.suffix)
        coverage_rows_count = self.get_wc_l(coverage_path)
        d = {"sample_name": sample_name,
             "source_paths": ",".join(raw_reads_files_list),
             "is_source_exists": all([os.path.isfile(i) for i in raw_reads_files_list]),
             "coverage_path": coverage_path,
             "is_coverage_exists": os.path.isfile(coverage_path),
             "is_coverage_valid": coverage_rows_count == mainInitializer.target_length}
        return self.dict2pd_series(d)

    @staticmethod
    def assembly_df_from_series_list(series_list: list, sorting_col_name: str):
        output_df = pd.DataFrame()
        for series in series_list:
            if len(output_df) == 0:
                output_df = pd.DataFrame([series])
            else:
                output_df = output_df.append(series, ignore_index=True)
        return output_df.sort_values(sorting_col_name)

    @staticmethod
    def multi_core_queue(func, queue: list, processes: int = int(subprocess.getoutput("nproc").strip())):
        import multiprocessing
        pool = multiprocessing.Pool(processes=processes)
        output = pool.map(func, queue)
        pool.close()
        pool.join()
        return output

    def combine_dfs(self):
        mp_result = self.multi_core_queue(func=self.process_sampledata_row, queue=self._queue)
        self._verified_df = self.assembly_df_from_series_list(series_list=mp_result, sorting_col_name="sample_name")
        no_raw_reads_df = self._verified_df.loc[self._verified_df["is_source_exists"] == False]
        if len(no_raw_reads_df) > 0:
            print("""Warning! The sources of following samples have not been found: '{a}'.
                  Please check the provided sample data file: '{b}'
                  """.format(a="', '".join(no_raw_reads_df["sample_name"].values.tolist()), b=mainInitializer.sampledata))
        self._no_coverages_df = self._verified_df.loc[(self._verified_df["is_source_exists"] == True) & (self._verified_df["is_coverage_valid"] == False)]

    def export(self):
        sampledatas_2d_array = [i for i in self._queue if i[0] in self._no_coverages_df["sample_name"].values.tolist()]
        os.makedirs(os.path.dirname(mainInitializer.output), exist_ok=True)
        Utilities.dump_2d_array(sampledatas_2d_array, file=mainInitializer.output)
        print("Done. \nFiles to process: {} \nDumped sample data: '{}'".format(len(self._no_coverages_df), mainInitializer.output))
        if mainInitializer.debugging_bool:
            debug_table = "{}_debug.tsv".format(mainInitializer.output)
            self._verified_df.to_csv(debug_table, sep='\t', header=True, index=False)
            print("Dumped debug table: '{}'".format(debug_table))

    def main(self):
        self.combine_dfs()
        self.export()


if __name__ == '__main__':
    mainInitializer = Initializer()
    verifier = CoveragesVerifier()
    verifier.main()
