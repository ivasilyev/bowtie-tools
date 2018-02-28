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
    starting_parser.add_argument("-o", "--output", required=True,
                                 help="Output file or directory to place auto named file")
    return starting_parser.parse_args()


def parse_namespace():
    namespace = parse_args()
    return namespace.sampledata, namespace.annotation, namespace.mask, namespace.output


def multi_core_queue(function_name, queue):
    pool = multiprocessing.Pool()
    output = pool.map(function_name, queue)
    pool.close()
    pool.join()
    return output


def is_sampledata_exists(single_sampledata_row):
    single_sampledata_list = single_sampledata_row.split('\t')
    if len(single_sampledata_list) > 1:
        return all([os.path.isfile(i.strip()) for i in single_sampledata_list[1:] if len(i.strip()) > 0])
    return False


def is_coverage_exists(single_sampledata_row):
    sample_name = single_sampledata_row.split('\t')[0].strip()
    coverage_file = outputDir + "Statistics/" + sample_name + "_" + outputMask + "_coverage.txt"
    if len(sample_name) > 0:
        if os.path.isfile(coverage_file):
            return sample_name, int(subprocess.getoutput("wc -l < " + coverage_file).split('\n')[0]) == annotationRowsCount
    return sample_name, False


def dict2pd_series(dictionary):
    output = pd.Series()
    for key in dictionary:
        output = output.set_value(key, dictionary[key])
    return output


def mp_verify_sd_and_cov(single_sampledata_row):
    sample_name, coverage_bool = is_coverage_exists(single_sampledata_row)
    return dict2pd_series({"sample_name": sample_name, "is_sampledata_exists": is_sampledata_exists(single_sampledata_row), "is_coverage_exists": coverage_bool})


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
    inputSampleDataFileName, inputReferenceAnnotation, outputMask, outputDir = parse_namespace()
    annotationRowsCount = int(subprocess.getoutput("wc -l < " + inputReferenceAnnotation).split('\n')[0])
    inputSampleDataRowsList = file_to_list(inputSampleDataFileName)
    verifiedSamplesDataFrame = assembly_df_from_series_list(multi_core_queue(mp_verify_sd_and_cov, inputSampleDataRowsList), "sample_name")
    processedDF = verifiedSamplesDataFrame.loc[(verifiedSamplesDataFrame["is_sampledata_exists"] == True) & (verifiedSamplesDataFrame["is_coverage_exists"] == True)]
    noSampleDataDF = verifiedSamplesDataFrame.loc[verifiedSamplesDataFrame["is_sampledata_exists"] == False]
    toDoDF = verifiedSamplesDataFrame.loc[(verifiedSamplesDataFrame["is_sampledata_exists"] == True) & (verifiedSamplesDataFrame["is_coverage_exists"] == False)]
    for dfDescription, dfName in zip(["Successfully processed files with extracted coverage:", "Files without source:", "Files without coverage:"], [processedDF, noSampleDataDF, toDoDF]):
        print(dfDescription, len(dfName))
    if len(noSampleDataDF) > 0:
        print("Warning! The sources of following samples have not been found: '" + "', '".join(noSampleDataDF["sample_name"].values.tolist()) +  "'.\nPlease check the provided data: " + inputSampleDataFileName)
    if len(toDoDF) > 0:
        if os.path.isdir(outputDir):
            outputDir = ends_with_slash(outputDir)
            outputSampleDataFileName = outputDir + get_time() + ".sampledata"
        else:
            outputSampleDataFileName = outputDir
        list_based_dict_export(toDoDF["sample_name"].values.tolist(), {i.split('\t')[0].strip(): i.split('\t')[1:] for i in inputSampleDataRowsList if len(i.split('\t')[0].strip()) > 0}, outputSampleDataFileName)
        print("Files to process have been dumped into:", outputSampleDataFileName)
