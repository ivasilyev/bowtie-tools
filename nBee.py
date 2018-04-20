#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import multiprocessing
import subprocess
import sys
import logging
import argparse

assert sys.version_info >= (2, 7)


def parse_args():
    starting_parser = argparse.ArgumentParser(description="This script will perform single alignment & postprocessing for every colorspace or non-colorspace sequence file specified in the input file list. \nRequired software: bowtie, bowtie2, samtools, bedtools. \nIf all coverage extract parameters specified, the script also will look for 'sam2coverage.py' within the same directory")
    starting_parser.add_argument("-i", "--input", required=True,
                                 help="Input list containing two tab-delimited columns for colorspace or non-colorspace sequences and three for paired-end sequences: sample name and absolute path(s). May contain a header")
    starting_parser.add_argument("-r", "--refdata", default=None,
                                 help="Linker file generated by the 'cook_the_reference.py' script. Contains tab-delimited columns: original FASTA file name, bowtie colorspace index, bowtie2 index, FASTA samtools index, genome lengths file etc.")
    starting_parser.add_argument("-b", "--bwti", default=None,
                                 help="Mask of bowtie or bowtie2 indexes")
    starting_parser.add_argument("-f", "--fai", default=None,
                                 help="(Optional) FAI file created by the 'samtools faidx' command")
    starting_parser.add_argument("-g", "--genome", default=None,
                                 help="(Optional) Genome file without a header containing tab-delimited sequence headers and lengths for single sequence chunk")
    starting_parser.add_argument("-a", "--annotation", default=None,
                                 help="(Optional) Text file supplied with a header containing additional information about sequence headers across all sequence chunks")
    starting_parser.add_argument("-m", "--mask", default="",
                                 help="Mask to be added to resulting files containing reference database name, date etc")
    starting_parser.add_argument("-t", "--threads", default=None, type=int,
                                 help="(Optional) Number of CPU cores to use, maximal by default")
    starting_parser.add_argument("-n", "--no_coverage", default=False, action='store_true',
                                 help="(Optional) If selected, cancels coverage extraction")
    starting_parser.add_argument("-o", "--output", required=True,
                                 help="Output directory")
    return starting_parser.parse_args()


def verify_path(path, file_or_dir):
    if path:
        if file_or_dir == "file" and os.path.isfile(path) is True:
            return os.path.abspath(path)
        if file_or_dir == "dir" and os.path.isdir(path) is True:
            return os.path.abspath(ends_with_slash(path))
        else:
            print("Cannot comprehend " + path + " as \"" + file_or_dir + "\"! Exiting...")
            sys.exit(2)


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


def create_dirs(paths_list):
    paths_list_with_slashes = [ends_with_slash(path_with_slash) for path_with_slash in paths_list]
    for path4check in paths_list_with_slashes:
        is_path_exists(path4check)
    return paths_list_with_slashes


def parse_namespace():
    namespace = parse_args()
    namespace.input = verify_path(namespace.input, 'file')
    if namespace.refdata:
        namespace.refdata = verify_path(namespace.refdata, 'file')
    elif not namespace.refdata and all(i for i in [namespace.bwti, namespace.fai, namespace.genome]):
        namespace.bwti, namespace.fai, namespace.genome = [verify_path(i, 'file') for i in [namespace.bwti, namespace.fai, namespace.genome]]
    else:
        for input_flag, input_arg in zip(["bwti", "fai", "genome"], [namespace.bwti, namespace.fai, namespace.genome]):
            if not input_arg:
                raise ValueError("Missing argument: '" + input_flag + "'")
    namespace.fai = verify_path(namespace.fai, 'file')
    namespace.genome = verify_path(namespace.genome, 'file')
    namespace.annotation = verify_path(namespace.annotation, 'file')
    if not namespace.threads:
        namespace.threads = subprocess.getoutput("nproc")
    namespace.output = create_dirs([namespace.output, namespace.output + "/Statistics", namespace.output + "/Mapped_reads", namespace.output + "/Non-mapped_reads"])[0]
    return namespace.input, namespace.refdata, namespace.bwti, namespace.fai, namespace.genome, namespace.annotation, namespace.mask, str(namespace.threads), namespace.no_coverage, namespace.output


def file_to_list(file):
    import re
    file_buffer = open(file, 'rU')
    output_list = [j for j in [re.sub('[\r\n]', '', i) for i in file_buffer] if len(j) > 0]
    file_buffer.close()
    return output_list


def file_append(string, file_to_append):
    file = open(file_to_append, 'a+')
    file.write(string)
    file.close()


def filename_only(string):
    return str(".".join(string.rsplit("/", 1)[-1].split(".")[:-1]))

# Wrappers


def single_core_queue(function_to_do, queue):
    return [function_to_do(queue_row) for queue_row in queue]


def multi_core_queue(function_name, queue):
    pool = multiprocessing.Pool(processes=int(cpuThreadsString))
    pool.map(function_name, queue)
    pool.close()
    pool.join()


def external_route(input_direction_list, output_direction):
    process = subprocess.Popen(input_direction_list, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (output, error) = process.communicate()
    process.wait()
    if error:
        print(error)
    if not output_direction:
        return output.decode("utf-8").replace('\r', '').replace('\n', '')
    else:
        file_append(output.decode("utf-8"), output_direction)


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

# Main workflow


def bowtie_it(single_sampledata):
    def is_path_endswith(mask, string):
        """
        This function is required to assign '--large-index' argument.
        Bowtie (Bowtie2) large indexes have extension 'ebwtl' ('bt2l') while common indexes end with 'ebwt' ('bt2').
        """
        return any(i.endswith(string) for i in subprocess.getoutput("ls -d " + mask + "*").split("\n"))
    sample_name, reads_file_path = [single_sampledata.split('\t')[0], '\t'.join(single_sampledata.split('\t')[1:])]
    reads_file_extension = '.' + reads_file_path.split('.')[-1]
    external_output_0 = outputDir + "Mapped_reads/" + sample_name + "_" + outputMask + ".sam"
    external_output_1 = outputDir + "Non-mapped_reads/" + sample_name + "_no_" + outputMask + reads_file_extension
    unmapped_reads_files_list = [external_output_1]
    if reads_file_path.endswith(".csfasta"):
        if len(list(filter(None, reads_file_path.split('\t')))) != 1:
            logging.fatal("Failed to parse single reads! The sampledata must contain exactly 2 columns!")
            sys.exit(2)
        logging.info("Mapping colorspace reads with bowtie")
        external_log = outputDir + "Statistics/" + sample_name + "_" + outputMask + "_bowtie.log"
        cmd = ['bowtie', '-f', '-C', '-S', '-t', '-v', '3', '-k', '1']
        if inputRefData:
            bwt_index = inputRefDataList[1]
        else:
            bwt_index = referenceBwtMask
        if is_path_endswith(bwt_index, "ebwtl"):
            cmd.append("--large-index")
        cmd.extend(['--threads', cpuThreadsString, '--un', external_output_1, bwt_index, reads_file_path, external_output_0])
    else:
        external_log = outputDir + "Statistics/" + sample_name + "_" + outputMask + "_bowtie2.log"
        paired_reads_list = [i for i in reads_file_path.split('\t') if len(i) > 0]
        if inputRefData:
            bwt_index = inputRefDataList[2]
        else:
            bwt_index = referenceBwtMask
        if any(reads_file_extension == "." + i for i in ["zip", "tar", "gz", "bz2"]):
            cmd = ['bowtie2', '--un-conc-gz', external_output_1]
        else:
            cmd = ['bowtie2', '--un', external_output_1]
        if len(paired_reads_list) > 2 or len(paired_reads_list) == 0:
            logging.warning("Failed to parse paired reads! The sampledata must contain exactly 3 columns!")
            return
        elif len(paired_reads_list) == 1:
            logging.info("Mapping single reads with bowtie2")
            cmd.extend(['-x', bwt_index, reads_file_path, '--threads', cpuThreadsString, '-S', external_output_0])
        else:
            logging.info("Mapping mate-pair reads with bowtie2")
            cmd.extend(['-x', bwt_index, '-1', paired_reads_list[0], '-2', paired_reads_list[1], '--threads', cpuThreadsString, '-S', external_output_0])
            unmapped_reads_files_list = [external_output_1.replace(reads_file_extension, '.1' + reads_file_extension), external_output_1.replace(reads_file_extension, '.2' + reads_file_extension)]
    make_cleanup([external_output_0, external_output_1, external_log])
    external_route(cmd, external_log)
    if all(os.path.isfile(i) for i in unmapped_reads_files_list):
        file_append('\t'.join(j for j in [sample_name] + unmapped_reads_files_list) + '\n', outputDir + "Statistics/_non-mapped_reads_" + outputMask + '_' + currentTime + ".sampledata")
    else:
        logging.warning("Some output files are missing: ", str(unmapped_reads_files_list))
    if os.path.isfile(external_output_0):
        file_append(sample_name + '\t' + external_output_0 + '\n', outputDir + "Statistics/_mapped_reads_" + outputMask + '_' + currentTime + ".sampledata")
        logging.info("Successfully aligned: " + reads_file_path)
    else:
        logging.warning("The output SAM file is missing: " + external_output_0)


def coverage_extract(single_sampledata):
    sample_name = single_sampledata.split('\t')[0]
    external_input = outputDir + "Mapped_reads/" + sample_name + "_" + outputMask + ".sam"
    external_log_0 = outputDir + "Statistics/" + sample_name + "_" + outputMask + "_coverage_extract.log"
    external_log_1 = str(outputDir + "Statistics/" + filename_only(sys.argv[0]) + "_" + outputMask + "_master.log")
    make_cleanup([external_log_0])
    cmd_list = ['python3', ends_with_slash(os.path.dirname(os.path.realpath(sys.argv[0]))) + 'sam2coverage.py', '-i', external_input, '-f', referenceFai, '-g', referenceGenomeLengths, '-l', external_log_1, '-o', outputDir]
    if referenceAnnotation:
        cmd_list = cmd_list[:-4] + ['-a', referenceAnnotation] + cmd_list[-4:]
    external_route(cmd_list, external_log_0)
    file_append(sample_name + '\t' + external_log_0.replace("_coverage_extract.log", "_coverage.txt") + '\n', outputDir + "Statistics/_coverages_" + outputMask + '_' + currentTime + ".sampledata")
    logging.info("Successfully extracted coverage from: " + external_input)


def get_time():
    from datetime import datetime
    now = datetime.now()
    output_list = []
    for time_unit in [now.year, now.month, now.day, now.hour, now.minute, now.second]:
        time_unit = str(time_unit)
        if len(time_unit) < 2:
            time_unit = '0' + time_unit
        output_list.append(time_unit)
    return '-'.join(output_list)


if __name__ == "__main__":
    inputSampleData, inputRefData, referenceBwtMask, referenceFai, referenceGenomeLengths, referenceAnnotation, outputMask, cpuThreadsString, noCoverageExtractionBool, outputDir = parse_namespace()
    scriptDir = ends_with_slash(os.path.dirname(os.path.realpath(sys.argv[0])))
    if inputRefData:
        inputRefDataList = file_to_list(inputRefData)[0].split('\t')
        # refdata columns: fasta, bwti, bwti2, fai, genome, annotation
        referenceFai, referenceGenomeLengths = inputRefDataList[3:5]
    inputSampleDataList = file_to_list(inputSampleData)
    for processed_sampledata_row in inputSampleDataList:
        if processed_sampledata_row.count('\t') < 1:
            print("Cannot parse sample data table: not enough columns!")
            sys.exit(2)
        if False in [os.path.isfile(input_path) for input_path in processed_sampledata_row.split('\t')[1:]]:
            inputSampleDataList.remove(processed_sampledata_row)
    currentTime = get_time()
    logging.basicConfig(format=u'%(levelname)-8s [%(asctime)s] %(message)s', level=logging.DEBUG, filename=outputDir + "Statistics/" + "_".join([filename_only(sys.argv[0]), outputMask, currentTime, "master.log"]))
    logging.info("Command: " + ' '.join(str(i) for i in sys.argv))
    logging.info("Main process ID: " + str(os.getpid()))
    logging.info("Using CPU threads: " + cpuThreadsString)
    single_core_queue(bowtie_it, inputSampleDataList)
    if referenceFai and referenceGenomeLengths and not noCoverageExtractionBool:
        if os.path.isfile(scriptDir + "sam2coverage.py"):
            logging.info("Starting coverage extract...")
            single_core_queue(coverage_extract, inputSampleDataList)
        else:
            logging.fatal("Input parameters were specified but the coverage extractor script has not been found!")
            logging.fatal("Please put 'sam2coverage.py' into the directory: " + scriptDir)
            sys.exit(2)
    logging.info("ALIGNMENT COMPLETED")
    print("ALIGNMENT COMPLETED")
