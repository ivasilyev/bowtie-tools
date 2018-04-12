#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import sys
import argparse
import threading
import multiprocessing
import subprocess


def parse_args():
    starting_parser = argparse.ArgumentParser(description="This tool automatically fixes, cuts and indexes DNA sequences in FASTA format. \nRequired software: bowtie, bowtie2, samtools")
    starting_parser.add_argument("-i", "--input", required=True,
                                 help="Reference DNA sequence file in FASTA format")
    starting_parser.add_argument("-c", "--correct_headers", default=False, action="store_true",
                                 help="(Optional) Enables headers correction instead of the default headers counter")
    starting_parser.add_argument("-p", "--preserve_headers", default=False, action="store_true",
                                 help="(Optional) Do not modify headers")
    starting_parser.add_argument("-n", "--not_large_index", default=False, action="store_true",
                                 help="(Optional) Avoid usage of a 'large' index, enables sequence splitting")
    starting_parser.add_argument("-s", "--size", type=float, default=3.6,
                                 help="(Optional) Sequence chunk size in billion charactersif '--not_large_index' is present, 3.6 by default")
    starting_parser.add_argument("-t", "--threads", type=int, default=int(multiprocessing.cpu_count()),
                                 help="(Optional) Number of CPU cores to speed up the chunks processing, maximal by default")
    starting_parser.add_argument("-o", "--output", required=True,
                                 help="Output directory. Must not exist and has to be created")
    return starting_parser.parse_args()


def is_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        print(exception)
        print("Please set up new folder or remove the existing: " + path)
        sys.exit(2)


def parse_namespace():
    namespace = parse_args()
    namespace.input = str(os.path.abspath(namespace.input))
    is_path_exists(namespace.output)
    namespace.output = ends_with_slash(os.path.abspath(namespace.output))
    for path in [namespace.input, namespace.output]:
        if ' ' in path:
            print("The path must not contain spaces: " + path + "\n Exiting...")
            sys.exit(2)
    if namespace.preserve_headers is True:
        print("Preserving headers")
    return str(os.path.abspath(namespace.input)), namespace.correct_headers, namespace.preserve_headers, int(float(namespace.size) * 10 ** 9), int(namespace.threads), namespace.not_large_index, namespace.output


def file_to_str(file):
    file_parsed = open(file, "rU").read()
    return file_parsed


def file_append(string, file_to_append):
    file = open(file_to_append, "a+")
    file.write(string)
    del string
    file.close()


def var_to_file(var_to_write, file_to_write):
    file = open(file_to_write, 'w')
    file.write(var_to_write)
    del var_to_write
    file.close()


def ends_with_slash(string):
    if string.endswith("/"):
        return string
    else:
        return str(string + "/")


def filename_only(string):
    return str(".".join(string.rsplit("/", 1)[-1].split(".")[:-1]))


def list_to_file(list_to_write, file_to_write):
    file = open(file_to_write, "w")
    file.write("".join(str(i) for i in list_to_write if i if i is not None))
    file.close()


def multi_core_queue(function_name):
    pool = multiprocessing.Pool(threadsNumber)
    pool.map(function_name, chunks)
    pool.close()
    pool.join()


def external_route(input_direction, output_direction):
    cmd = input_direction.split(" ")
    process = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (output, error) = process.communicate()
    process.wait()
    if error:
        print(error)
    if not output_direction:
        return output.decode("utf-8")
    else:
        file_append(output.decode("utf-8"), output_direction)

# Main workflow


def process_header(input_header):
    if preserveHeadersBool:
        return input_header
    elif correctHeadersBool:
        return '>' + re.sub('_+', '_', re.sub('\W+', '_', input_header)).strip('_') + '\n'


def fasta_headers_fix(sequence_file):
    print("Parsing the whole sequence...")
    file_parsed = open(sequence_file, 'rU')
    output_file = outputDir + filename_only(sequence_file) + ".fasta"
    headers_number_string = subprocess.getoutput("grep -c \> " + sequence_file)
    headers_zfill_number = len(headers_number_string)
    # raise ValueError("Cannot count sequences number:", sequence_file)
    try:
        output_buffer = ""
        fixed_headers_dict = {}
        headers_counter = 0
        for string in file_parsed:
            if string.startswith('>'):
                headers_counter += 1
                processed_string = process_header(string)
                if not processed_string:
                    processed_string = ">ID" + str(headers_counter).zfill(headers_zfill_number) + "\n"
                fixed_headers_dict[re.sub('[\n>]', '', processed_string)] = re.sub('[\n>]', '', string).replace('\t', ' ')
                string = processed_string
            output_buffer += string
        var_to_file(output_buffer, output_file)
        del output_buffer
    except MemoryError:
        print("Not enough memory! Using the per-line processing...")
        fixed_headers_dict = {}
        headers_counter = 0
        for string in file_parsed:
            if string.startswith('>'):
                headers_counter += 1
                processed_string = process_header(string)
                if not processed_string:
                    processed_string = ">ID" + str(headers_counter).zfill(headers_zfill_number) + "\n"
                fixed_headers_dict[re.sub('[\n>]', '', processed_string)] = re.sub('[\n>]', '', string).replace('\t', ' ')
                string = processed_string
            file_append(string, output_file)
    file_parsed.close()
    print("Completed FASTA headers fixing for " + sequence_file + ". \n" + str(len(fixed_headers_dict)) + " of " + headers_number_string + " headers were processed!")
    del file_parsed
    return fixed_headers_dict


def sequence_chop(sequence_file, input_chunk_size):
    default_max_chunk_size = int(3.6 * 10 ** 9)  # Bowtie limit
    if input_chunk_size > default_max_chunk_size:
        print("The chunk size must not be larger than 3.6 billion characters!")
        input_chunk_size = default_max_chunk_size
    print("Using the maximal chunk size: " + str(input_chunk_size))
    sequence_buffer = open(sequence_file, 'rU')
    sequences_number = sum(1 for sequence_buffer_line in sequence_buffer)
    sequence_buffer.close()
    sequence_buffer = open(sequence_file, 'rU')
    chunks_files_list = []
    chunk_buffer = ""
    chunk_index = 1
    sequence_index = 0
    single_sequence = ""
    for sequence_buffer_line in sequence_buffer:
        single_sequence += sequence_buffer_line
        sequence_index += 1
        if not sequence_buffer_line.startswith('>'):
            if len(single_sequence) + len(chunk_buffer) > input_chunk_size or (len(single_sequence) + len(chunk_buffer) < input_chunk_size and sequence_index == sequences_number):
                chunk_file = outputDir + filename_only(inputFile) + "_chunk_" + str(chunk_index) + "." + inputFile.split(".")[-1]
                var_to_file(chunk_buffer, chunk_file)
                chunks_files_list.append(chunk_file)
                chunk_buffer = ""
                chunk_index += 1
            chunk_buffer += single_sequence
            single_sequence = ""
    del single_sequence, chunk_buffer
    print(str(len(chunks_files_list)) + " chunks were created!")
    return chunks_files_list


def bowtie_build(chunk):
    external_route("bowtie-build -C " + chunk + " " + outputDir + filename_only(chunk) + "_colorspace",
                   outputDir + filename_only(chunk) + "_bowtie-build.log")
    print("Created bowtie colorspace index for " + chunk)


def bowtie2_build(chunk):
    external_route("bowtie2-build " + chunk + " " + outputDir + filename_only(chunk) + "_bowtie2",
                   outputDir + filename_only(chunk) + "_bowtie2-build.log")
    print("Created bowtie2 index for " + chunk)


def samtools_faidx(chunk):
    external_route("samtools faidx " + chunk,
                   outputDir + filename_only(chunk) + "_samtools_faidx.log")
    os.rename(chunk + ".fai", outputDir + filename_only(chunk) + "_samtools.fai")
    print("Created samtools index for " + chunk)
    fai2genome(chunk)


def string_process(string):
    string = re.sub('[\r\n]', '', string)
    if string:
        columns = []
        for column in [0, 1]:
            if column is not None:
                try:
                    columns.append(string.strip().split("\t")[column])
                except IndexError:
                    print("The table string \"" + string + "\" does not contain the column with index " + str(column) + "!")
        out = "\t".join(str(i) for i in columns) + str("\n")
        return out


def fai2genome(chunk):
    """Process existing fasta index, depends from 'samtools_faidx' function"""
    strings = list(filter(None, file_to_str(outputDir + filename_only(chunk) + "_samtools.fai").replace("\r", "\n").split("\n")))
    strings_processed = []
    for string in strings:
        strings_processed.append(string_process(string))
    output = "".join(str(i) for i in strings_processed if len(i) > 0)
    file_append(output, outputDir + filename_only(chunk) + "_samtools.genome")
    print("Created genome index for " + chunk)
    output_annotation_file = outputDir + filename_only(inputFile) + "_annotation.txt"
    if not os.path.isfile(output_annotation_file):
        var_to_file("Reference_ID\tReference_Length\n", output_annotation_file)
    file_append(output, output_annotation_file)
    del strings, strings_processed, output
    print("Compiled compact annotation for " + chunk)


def asynchronous_chunk_processing(chunk):
    """Process chopped or whole FASTA"""
    # Initializing threads
    t1 = threading.Thread(target=bowtie_build, args=(chunk,))
    t2 = threading.Thread(target=bowtie2_build, args=(chunk,))
    t3 = threading.Thread(target=samtools_faidx, args=(chunk,))
    # Starting threads
    t1.start()
    t2.start()
    t3.start()
    # Joining threads
    t1.join()
    t2.join()
    t3.join()


def add_former_headers(headers_dict, annotation_file):
    annotation_file_buffer = open(annotation_file, 'rU')
    annotation_file_output = ""
    for annotation_file_row in annotation_file_buffer:
        reference_id, reference_length = annotation_file_row.split('\t')
        annotation_file_output += '\t'.join([reference_id, headers_dict[reference_id], reference_length])
    annotation_file_buffer.close()
    annotation_file_output_buffer = open(annotation_file, 'w')
    annotation_file_output_buffer.write(annotation_file_output)
    annotation_file_output_buffer.close()
    del annotation_file, annotation_file_output, annotation_file_buffer, annotation_file_output_buffer
    print("Added former headers to annotation file!")


def make_refdata(chunks_list):
    chunks_iteration = 0
    for chunk in chunks_list:
        if chunks_iteration == 0:
            var_to_file("", outputDir + filename_only(inputFile) + ".refdata")
        file_append(str(chunk) + '\t' +
                    outputDir + filename_only(chunk) + "_colorspace" + '\t' +
                    outputDir + filename_only(chunk) + "_bowtie2" + '\t' +
                    outputDir + filename_only(chunk) + "_samtools.fai" + '\t' +
                    outputDir + filename_only(chunk) + "_samtools.genome" + '\t' +
                    outputDir + filename_only(inputFile) + "_annotation.txt" + '\n',
                    outputDir + filename_only(inputFile) + ".refdata")
        chunks_iteration += 1


def sequence2chunks_list(sequence_file):
    if notLargeIndexes and os.path.getsize(sequence_file) > 2 ** 32:  # Bowtie limit
        print("The reference file is too large! Splitting sequences...")
        chunks_list = sequence_chop(sequence_file, chunkSize)
    else:
        chunks_list = [outputDir + filename_only(sequence_file) + ".fasta"]
    return chunks_list


if __name__ == "__main__":
    inputFile, correctHeadersBool, preserveHeadersBool, chunkSize, threadsNumber, notLargeIndexes, outputDir = parse_namespace()
    fixedHeadersDict = fasta_headers_fix(inputFile)
    fixedHeadersDict["Reference_ID"] = "Former_ID"
    chunks = sequence2chunks_list(outputDir + filename_only(inputFile) + ".fasta")
    multi_core_queue(asynchronous_chunk_processing)
    add_former_headers(fixedHeadersDict, outputDir + filename_only(inputFile) + "_annotation.txt")
    make_refdata(chunks)
    print("Created refdata linker: " + outputDir + filename_only(inputFile) + ".refdata\n\nReference indexing has been completed.")
