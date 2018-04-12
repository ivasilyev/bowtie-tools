#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import subprocess
import re
import os
import logging
import argparse
import pandas as pd

assert sys.version_info >= (2, 7)


def parse_args():
    starting_parser = argparse.ArgumentParser(description="This tool automatically converts, sorts, indexes and extracts coverage from SAM file. \nRequired software: samtools, bedtools")
    starting_parser.add_argument("-i", "--input", required=True,
                                 help="Alignment data file in SAM or BAM format")
    starting_parser.add_argument("-f", "--fai", required=True,
                                 help="FAI file created by the 'samtools faidx' tool")
    starting_parser.add_argument("-g", "--genome", required=True,
                                 help="Genome file without a header containing tab-delimited sequence headers and lengths for single sequence chunk")
    starting_parser.add_argument("-a", "--annotation", default=None,
                                 help="(Optional) Genome file supplied with a header containing tab-delimited sequence headers and lengths across all sequence chunks")
    starting_parser.add_argument("-t", "--threshold", type=float, default=None,
                                 help="(Optional) Minimal part of the mapped reference entry length to the entry length to be reported")
    starting_parser.add_argument("-n", "--non_zero", action='store_true', default=False,
                                 help="(Optional) If specified, only genomes with non-zero coverage will be reported")
    starting_parser.add_argument("-l", "--log", default=None,
                                 help="(Optional) Log file, the 'Statistics/<input sam file name>_coverage.log' by default")
    starting_parser.add_argument("-b", "--buffer", default=10000, type=int,
                                 help="(Optional) Numbers of row to fetch in time, 10000 by default")
    starting_parser.add_argument("-o", "--output", required=True,
                                 help="Output directory. The subdirectories 'Statistics', 'Mapped_reads' also would be created")
    return starting_parser.parse_args()


def is_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        print(exception)


def ends_with_slash(string):
    if string.endswith("/"):
        return string
    else:
        return str(string + "/")


def parse_namespace():
    namespace = parse_args()
    if (not namespace.input.endswith(".sam")) and (not namespace.input.endswith(".bam")):
        print("Wrong file extension! Supported formats: SAM, BAM")
        sys.exit(2)
    namespace.output = create_dirs([namespace.output, namespace.output + "/Statistics", namespace.output + "/Mapped_reads"])[0]
    if not namespace.log:
        namespace.log = namespace.output + "Statistics/" + filename_only(namespace.input) + "_sam2coverage.log"
    return namespace.input, namespace.fai, namespace.genome, namespace.annotation, namespace.threshold, namespace.non_zero, namespace.log, namespace.buffer, namespace.output


def create_dirs(paths_list):
    paths_list_with_slashes = [ends_with_slash(path_with_slash) for path_with_slash in paths_list]
    for path4check in paths_list_with_slashes:
        is_path_exists(path4check)
    return paths_list_with_slashes


def external_route(input_direction_list, output_direction):
    process = subprocess.Popen(input_direction_list, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (output, error) = process.communicate()
    process.wait()
    if error:
        print(error)
    if not output_direction:
        return output.decode("utf-8")
    else:
        file_append(output.decode("utf-8"), output_direction)


def file_to_list(file):
    file_buffer = open(file, 'rU')
    output_list = []
    for file_row in file_buffer:
        file_row = re.sub('[\n\r]', '', file_row)
        if len(file_row) > 0:
            output_list.append(file_row)
    file_buffer.close()
    return output_list


def var_to_file(var_to_write, file_to_write):
    file = open(file_to_write, 'w')
    file.write(var_to_write)
    file.close()


def file_append(string, file_to_append):
    file = open(file_to_append, 'a+')
    file.write(string)
    file.close()


def filename_only(string):
    return str(".".join(string.rsplit("/", 1)[-1].split(".")[:-1]))


def sam2bam(sam_file):
    sample_name = filename_only(sam_file)
    external_input = sam_file
    external_output = outputDir + "Mapped_reads/" + sample_name + ".bam"
    external_log = outputDir + "Statistics/" + sample_name + "_sam2bam.log"
    make_cleanup([external_output, external_log])
    external_route(['samtools', 'import', referenceFai, external_input, external_output], external_log)
    logging.info("Successfully converted: " + external_output)


def sort_bam(sam_file):
    sample_name = filename_only(sam_file)
    external_input = outputDir + "Mapped_reads/" + sample_name + ".bam"
    external_output = outputDir + "Mapped_reads/" + sample_name + ".sorted.bam"
    external_log = outputDir + "Statistics/" + sample_name + "_sort_bam.log"
    make_cleanup([external_output, external_output + ".bai", external_log])
    tmp = subprocess.getoutput("ls -d " + outputDir + "Mapped_reads/" + sample_name + "*.bam | xargs rm -f")
    del tmp
    external_route(['samtools', 'sort', external_input, '-o', external_output], external_log)
    logging.info("Successfully sorted: " + external_input)


def bam2coverage(sam_file):
    sample_name = filename_only(sam_file)
    external_input = outputDir + "Mapped_reads/" + sample_name + ".sorted.bam"
    external_output = outputDir + "Statistics/" + sample_name + "_genomeCoverageBed.txt"
    make_cleanup([external_output])
    external_route(['genomeCoverageBed', '-ibam', external_input], external_output)
    logging.info("Successfully calculated coverage: " + external_input)


def stack_coverage(sam_file):
    sample_name = filename_only(sam_file)
    input_file = outputDir + "Statistics/" + sample_name + "_genomeCoverageBed.txt"
    output_file = outputDir + "Statistics/" + sample_name + "_pos_bp.txt"
    make_cleanup([output_file])
    coverages_list = subprocess.getoutput("sort " + input_file).split('\n')
    output_list = []
    counting_id = ""
    row_processing_buffer = []
    # genomecov file columns: reference sequence name, depth of coverage, breadth of coverage with that depth, sequence length, coverage ratio
    for row in coverages_list:
        row_list = row.split('\t')
        # remove service lines
        if len(row_list[0].strip()) > 0 and row_list[0].strip() != 'genome' and '*' not in row_list[0]:
            if row_list[0] == counting_id and int(row_list[1]) > 0:
                row_processing_buffer.append(row_list)
            else:
                if len(row_processing_buffer) > 0:
                    # output file columns: reference sequence name, maximal depth of coverage, total breadth of coverage, sequence length, coverage ratio, total mapped bases
                    row_output_buffer = [counting_id, 0, 0, int(row_processing_buffer[0][3]), 0, 0]
                    for mapped_row_list in row_processing_buffer:
                        row_output_buffer[1] = max(int(mapped_row_list[1]), int(row_output_buffer[1]))
                        row_output_buffer[2] += int(mapped_row_list[2])
                        row_output_buffer[4] += float(mapped_row_list[4])
                        row_output_buffer[5] = int(mapped_row_list[1]) * int(mapped_row_list[2])
                    output_list.append('\t'.join(str(column) for column in row_output_buffer))
                row_processing_buffer = []
                counting_id = row_list[0]
    var_to_file('\n'.join(output_row for output_row in output_list), output_file)
    logging.info("Successfully stacked: " + input_file)


def bam_index(sam_file):
    sample_name = filename_only(sam_file)
    external_input = outputDir + "Mapped_reads/" + sample_name + ".sorted.bam"
    external_log = outputDir + "Statistics/" + sample_name + "_index_bam.log"
    make_cleanup([external_log])
    external_route(['samtools', 'index', external_input], external_log)
    logging.info("Successfully indexed: " + external_input)


def bam2idxstats(sam_file):
    sample_name = filename_only(sam_file)
    external_input = outputDir + "Mapped_reads/" + sample_name + ".sorted.bam"
    external_output = outputDir + "Statistics/" + sample_name + "_idxstats.txt"
    make_cleanup([external_output])
    external_route(['samtools', 'idxstats', external_input], external_output)
    logging.info("Successfully retrieved mapped reads number: " + external_input)


def sam2stats(sam_file):
    sam_stats = subprocess.getoutput("samtools stats " + sam_file + " | grep ^SN | cut -f 2-")
    output_file = outputDir + "Statistics/" + filename_only(sam_file) + '_sam_stats.txt'
    make_cleanup([output_file])
    sam_total_reads = re.findall('raw total sequences:\t(\d*)', sam_stats)[0]
    sam_mapped_reads = re.findall('reads mapped:\t(\d*)', sam_stats)[0]
    sam_total_bp = re.findall('total length:\t(\d*)', sam_stats)[0]
    sam_mapped_bp = re.findall('bases mapped:\t(\d*)', sam_stats)[0]
    var_to_file(sam_stats, output_file)
    return sam_total_reads, sam_mapped_reads, sam_total_bp, sam_mapped_bp


def bedtools2coverage_sum(bedtools_coverage_file):
    file_parsed = open(bedtools_coverage_file, 'rU')
    coverage_sum = 0
    for raw_line in file_parsed:
        if raw_line is not None:
            processed_line = raw_line.replace('\r', '\n').replace('\n', '').split('\t')
            if not processed_line[0].startswith("genome\t") and '*' not in processed_line[0]:
                if processed_line[1] != 0:
                    coverage_sum += (int(processed_line[1]) * int(processed_line[2]))
    return coverage_sum


def idxstats2reads_sum(idxstats_file):
    file_parsed = open(idxstats_file, 'rU')
    reads_sum = 0
    for raw_line in file_parsed:
        if raw_line is not None:
            processed_line = raw_line.replace('\r', '\n').replace('\n', '').split('\t')
            if not processed_line[0].startswith("genome\t") and '*' not in processed_line[0]:
                reads_sum += int(processed_line[2])
    return reads_sum


def reference2statistics(sam_file):
    sample_name = filename_only(sam_file)
    output_coverage_file = outputDir + "Statistics/" + sample_name + "_coverage.txt"
    make_cleanup([output_coverage_file])
    pos_bp_file_name = outputDir + "Statistics/" + sample_name + "_pos_bp.txt"
    if len(file_to_list(pos_bp_file_name)) == 0:
        logging.critical("Failed to process: " + output_coverage_file)
        return
    sample_total_reads, sample_mapped_reads, sample_total_bp, sample_mapped_bp = [int(var) for var in sam2stats(sam_file)]
    reference_df = pd.read_table(referenceGenomeLengths, header='infer', sep='\t', names=['reference_id', 'id_bp'])
    stacked_coverages_df = pd.read_table(pos_bp_file_name, sep='\t', header='infer', names=["reference_id", "id_maximal_coverage_depth", "id_coverage_breadth", "id_bp", "id_coverage_breadth_to_id_bp", "id_mapped_bp"])
    genomes_coverages_df = pd.merge(reference_df, stacked_coverages_df.loc[:, [column for column in list(stacked_coverages_df) if column != "id_bp"]], on="reference_id", how='outer')
    genomes_coverages_df = genomes_coverages_df.loc[genomes_coverages_df['reference_id'] != 'genome']
    genomes_coverages_df["id_total_relative_abundance"] = genomes_coverages_df.loc[:, "id_mapped_bp"] / (genomes_coverages_df.loc[:, "id_bp"] * sample_total_bp)
    genomes_coverages_df["id_mapped_relative_abundance"] = genomes_coverages_df.loc[:, "id_mapped_bp"] / (genomes_coverages_df.loc[:, "id_bp"] * sample_mapped_bp)
    genomes_coverages_df["sample_total_reads"] = sample_total_reads
    genomes_coverages_df["sample_mapped_reads"] = sample_mapped_reads
    genomes_coverages_df["sample_total_bp"] = sample_total_bp
    genomes_coverages_df["sample_mapped_bp"] = sample_mapped_bp
    genomes_coverages_df["sample_average_total_reads_bp"] = float(sample_total_reads) / float(sample_total_bp)
    genomes_coverages_df["sample_average_mapped_reads_bp"] = float(sample_mapped_reads) / float(sample_total_bp)
    genomes_coverages_df["sample_mapped_reads_to_total_reads"] = float(sample_mapped_reads) / float(sample_total_reads)
    idxstats_df = pd.read_table(outputDir + "Statistics/" + sample_name + "_idxstats.txt", sep='\t', header='infer', names=["reference_id", "id_bp", "id_mapped_reads", "id_unmapped_reads"])
    genomes_coverages_df = pd.merge(genomes_coverages_df, idxstats_df.loc[:, [column for column in list(idxstats_df) if column != "id_bp"]], on="reference_id", how='outer')
    genomes_coverages_df["id_mapped_reads_per_million_sample_total_reads"] = genomes_coverages_df.loc[:, "id_mapped_reads"] * (10 ** 6) / float(sample_total_reads)
    genomes_coverages_df["id_mapped_reads_per_million_sample_mapped_reads"] = genomes_coverages_df.loc[:, "id_mapped_reads"] * (10 ** 6) / float(sample_mapped_reads)
    output_df = genomes_coverages_df.loc[:, ["reference_id", "id_bp", "id_coverage_breadth", "id_mapped_bp", "id_coverage_breadth_to_id_bp", "id_maximal_coverage_depth", "id_total_relative_abundance", "id_mapped_relative_abundance", "id_mapped_reads", "sample_total_reads", "sample_mapped_reads", "sample_total_bp", "sample_mapped_bp", "sample_average_total_reads_bp", "sample_average_mapped_reads_bp", "sample_mapped_reads_to_total_reads", "id_mapped_reads_per_million_sample_total_reads", "id_mapped_reads_per_million_sample_mapped_reads"]]
    if onlyCoveredOoutput:
        output_df = output_df[output_df.id_coverage_breadth.notnull()]
    else:
        output_df = output_df.fillna(0)
    for int_column in ["id_bp", "id_coverage_breadth", "id_mapped_bp", "id_maximal_coverage_depth", "id_mapped_reads", "sample_total_reads", "sample_mapped_reads", "sample_total_bp", "sample_mapped_bp"]:
        output_df[int_column] = output_df.loc[:, int_column].astype(int)
    output_df = output_df.loc[~output_df['reference_id'].str.contains('*', regex=False)]
    output_df.to_csv(output_coverage_file, sep='\t', index=False)
    logging.info("Successfully extracted coverage into: " + output_coverage_file)


def main_workflow(sam_file):
    if sam_file.endswith(".sam"):
        coverage_count_function_list = [sam2bam, sort_bam, bam2coverage, stack_coverage, bam_index, bam2idxstats, reference2statistics]
    else:
        coverage_count_function_list = [sort_bam, bam2coverage, stack_coverage, bam_index, bam2idxstats, reference2statistics]
    for coverage_count_function in coverage_count_function_list:
        coverage_count_function(sam_file)


def make_cleanup(files_list):
    removed_files_list = []
    removed_directories_list = []
    for file_name in files_list:
        try:
            os.remove(file_name)
            removed_files_list.append(file_name)
        except FileNotFoundError:
            continue
        except IsADirectoryError:
            import shutil
            shutil.rmtree(file_name)
            removed_directories_list.append(file_name)
    if len(removed_files_list) > 0:
        print("Removed files: " + ', '.join(removed_file for removed_file in removed_files_list))
    if len(removed_directories_list) > 0:
        print("Removed directories: " + ', '.join(removed_dir for removed_dir in removed_directories_list))


if __name__ == "__main__":
    inputSamFile, referenceFai, referenceGenomeLengths, referenceAnnotation, bedtoolsCoverageThresholdRate, onlyCoveredOoutput, existing_log, exportBuffer, outputDir = parse_namespace()
    logging.basicConfig(format=u'%(levelname)-8s [%(asctime)s] %(message)s', level=logging.DEBUG, filename=existing_log)
    logging.info("Started " + sys.argv[0])
    logging.info("Command: " + ' '.join(str(arg) for arg in sys.argv))
    main_workflow(inputSamFile)
    logging.info("Completed running " + sys.argv[0] + " for " + inputSamFile)
    logging.info("COVERAGE EXTRACTION COMPLETED")
    print("COVERAGE EXTRACTION COMPLETED")
