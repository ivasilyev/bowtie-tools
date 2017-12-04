#!/usr/bin/env python
# -*- coding: utf-8 -*-

import getopt
import sys
import re


def usage():
    print("\nUsage: " + sys.argv[0] + " -i <file> -n <int> -p <regex>" + "\n\n" +
          "-i/--input <file> \tA text  file;" + "\n" +
          "-n/--number <int> \tNumber of chunks to chop;" + "\n" +
          "-p/--pattern <regex> \tA pattern or regular expression to split supplied by backslash symbol if required;" + "\n" +
          "-o/--output <regex> \tOutput directory." + "\n")
    sys.exit(2)


def main():
    opts = ""
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:n:p:o:", ["help", "input=", "number=", "pattern=", "output="])
    except getopt.GetoptError as arg_err:
        print(str(arg_err))
        usage()
    i, n, p = (None,) * 3
    o = ''
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
        elif opt in ("-i", "--input"):
            i = str(arg)
        elif opt in ("-n", "--number"):
            try:
                n = int(arg)
            except ValueError:
                print("Chunks number must be integer! Exiting...")
                sys.exit(2)
        elif opt in ("-p", "--pattern"):
            p = str(arg)
        elif opt in ("-o", "--output"):
            o = str(arg)
    if not any(var is None for var in [i, n]):
        return i, n, p, o
    print("The parameters are not yet specified!")
    usage()


def is_path_exists(path):
    import os
    try:
        os.makedirs(path)
    except OSError as exception:
        print(exception)


def ends_with_slash(string):
    if string.endswith("/"):
        return string
    else:
        return str(string + "/")


def file_to_str(file):
    file_parsed = open(file, 'rU').read()
    return file_parsed


def regex_split(string, pattern):
    if pattern is None:
        pattern = "\n"
    print("Using pattern: " + repr(pattern))
    list_to_join = list(filter(None, re.split("(" + pattern + ")", string)))
    return list_to_join


def filename_only(string):
    return str(".".join(string.rsplit("/", 1)[-1].split(".")[:-1]))


def list_chop(input_list, number_of_chunks):
    max_chunks = len(input_list)
    if max_chunks < number_of_chunks:
        print("Cannot make " + str(number_of_chunks) + " chunks! \n" +
              "Will make " + str(max_chunks) + " chunks.")
        number_of_chunks = max_chunks
    chunks_length = int(len(input_list) / number_of_chunks)
    if sys.version_info >= (3, 0):
        chunks_list = [input_list[i:i + chunks_length] for i in range(0, (number_of_chunks - 1) * chunks_length, chunks_length)]
    else:
        chunks_list = [input_list[i:i + chunks_length] for i in xrange(0, (number_of_chunks - 1) * chunks_length, chunks_length)]
    chunks_list.append(input_list[(number_of_chunks - 1) * int(chunks_length):])
    return chunks_list


def export_two_dimensional_array(array):
    chunk = 1
    chunks_list = []
    for dim in array:
        output_file_name = outputDir + filename_only(inputFile) + "_chunk_" + str(chunk) + '.' + inputFile.split('.')[-1]
        list_to_file('', dim, output_file_name)
        chunks_list.append(output_file_name + '\n')
        chunk += 1
    chunks_list_file_name = outputDir + filename_only(inputFile) + ".chunkslist"
    list_to_file('', chunks_list, chunks_list_file_name)
    print("The chunks list saved into", chunks_list_file_name)


def list_to_file(header, list_to_write, file_to_write):
    header += ''.join(str(i) for i in list_to_write if i is not None)
    file = open(file_to_write, 'w')
    file.write(header)
    file.close()


########################################################
inputFile, chunksNumber, inputPattern, outputDir = main()
is_path_exists(outputDir)
outputDir = ends_with_slash(outputDir)
rawStrings = file_to_str(inputFile)
rawList = regex_split(rawStrings, inputPattern)
processedArray = list_chop(rawList, chunksNumber)
export_two_dimensional_array(processedArray)
print("File slicing completed")
