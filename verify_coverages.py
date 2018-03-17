#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import multiprocessing
import subprocess
import argparse
import pandas as pd
from nBee import ends_with_slash, file_to_list, get_time


def parse_args():
    starting_parser = argparse.ArgumentParser(description="The script to detect nBee or refdata2schedule pipeline faults")
    starting_parser.add_argument("-s", "--sampledata", required=True,
                                 help="Input list containing two tab-delimited columns for colorspace or non-colorspace sequences and three for paired-end sequences: sample name and absolute path(s). May contain a header")
    starting_parser.add_argument("-a", "--annotation", required=True,
                                 help="Text file supplied with a header containing additional information about sequence headers across all sequence chunks")
    starting_parser.add_argument("-m", "--mask", default='',
                                 help="Mask after the '<sample_name>_' and before the '_coverage.txt' substrings")
    starting_parser.add_argument("-d", "--debug", default=False, action='store_true',
                                 help="Create debugging table")
    starting_parser.add_argument("-o", "--output", required=True,
                                 help="Path to look for coverage files. Must contain subdirectories created by 'nBee.py' script")
    return starting_parser.parse_args()


def parse_namespace():
    namespace = parse_args()
    namespace.output = ends_with_slash(namespace.output)
    if not os.path.isdir(namespace.output):
        raise ValueError("Not a directory:", outputDir)
    return namespace.sampledata, namespace.annotation, namespace.mask, namespace.debug, namespace.output


class SampleDataString:
    def __init__(self, single_sampledata_row):
        self._single_sampledata_list = [i.strip() for i in single_sampledata_row.split('\t') if len(i) > 0]
        self.sample_name = self._single_sampledata_list[0]
        self.source_paths_list = self._single_sampledata_list[1:]
        self.coverage_path = outputDir + "Statistics/" + self.sample_name + "_" + outputMask + "_coverage.txt"

    def is_source(self):
        return all([os.path.isfile(i) for i in self.source_paths_list])

    def is_coverage(self):
        return os.path.isfile(self.coverage_path)

    def is_length(self):
        if self.is_coverage():
            return int(subprocess.getoutput("wc -l < " + self.coverage_path).split('\n')[0]) == annotationRowsCount
        return False


def dict2pd_series(dictionary):
    output = pd.Series()
    for key in dictionary:
        output.at[key] = dictionary[key]
    return output


def process_single_sampledata(row):
    sampledata = SampleDataString(row)
    return dict2pd_series({"sample_name": sampledata.sample_name, "source_paths": ":".join(sampledata.source_paths_list), "is_source_exists": sampledata.is_source(), "coverage_paths": sampledata.coverage_path, "is_coverage_exists": sampledata.is_coverage(), "is_coverage_valid": sampledata.is_length()})


def multi_core_queue(function_name, queue):
    pool = multiprocessing.Pool()
    output = pool.map(function_name, queue)
    pool.close()
    pool.join()
    return output


def assembly_df_from_series_list(series_list, sorting_col_name):
    output_df = pd.DataFrame()
    for series in series_list:
        if len(output_df) == 0:
            output_df = pd.DataFrame([series])
        else:
            output_df = output_df.append(series, ignore_index=True)
    return output_df.sort_values(sorting_col_name)


def list_based_dict_export(input_list, input_dict, output_file):
    output_string = ""
    for key in input_list:
        output_string += ('\t'.join([key, *input_dict[key]]) + '\n')
    file_wrapper = open(output_file, 'a+')
    file_wrapper.write(output_string)
    file_wrapper.close()


if __name__ == '__main__':
    inputSampleDataFileName, inputReferenceAnnotation, outputMask, debuggingBool, outputDir = parse_namespace()
    annotationRowsCount = int(subprocess.getoutput("wc -l < " + inputReferenceAnnotation).split('\n')[0])
    inputSampleDataRowsList = file_to_list(inputSampleDataFileName)
    verifiedSamplesDataFrame = assembly_df_from_series_list(multi_core_queue(process_single_sampledata, inputSampleDataRowsList), "sample_name")
    noSampleDataDF = verifiedSamplesDataFrame.loc[verifiedSamplesDataFrame["is_source_exists"] == False]
    toDoDF = verifiedSamplesDataFrame.loc[(verifiedSamplesDataFrame["is_source_exists"] == True) & (verifiedSamplesDataFrame["is_coverage_valid"] == False)]
    if len(noSampleDataDF) > 0:
        print("Warning! The sources of following samples have not been found: '" + "', '".join(noSampleDataDF["sample_name"].values.tolist()) + "'.\nPlease check the provided data: " + inputSampleDataFileName)
    timeString = get_time()
    outputSampleDataFileName = outputDir + timeString + ".sampledata"
    list_based_dict_export(toDoDF["sample_name"].values.tolist(), {i.split('\t')[0].strip(): i.split('\t')[1:] for i in inputSampleDataRowsList if len(i.split('\t')[0].strip()) > 0}, outputSampleDataFileName)
    if debuggingBool:
        verifiedSamplesDataFrame.to_csv(outputDir + timeString + "_debug.txt", sep='\t', index=False)
