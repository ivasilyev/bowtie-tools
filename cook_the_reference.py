#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
from modules.Utilities import Utilities
from modules.RefDataArray import RefDataArray


class Initializer:
    def __init__(self):
        self._namespace = self.parse_args()
        self.input_nfasta = self._namespace.input
        self.preserve_headers_bool = self._namespace.preserve_headers
        self.not_large_index_bool = self._namespace.not_large_index
        self.chunk_length = int(self._namespace.size * 10 ** 9)
        self.output_dir = Utilities.ends_with_slash(self._namespace.output)
        os.makedirs(self.output_dir, exist_ok=True)

    @staticmethod
    def parse_args():
        starting_parser = argparse.ArgumentParser(description="This tool automatically fixes, cuts and indexes DNA sequences in FASTA format. \nRequired software: bowtie, bowtie2, samtools")
        starting_parser.add_argument("-i", "--input", required=True,
                                     help="Reference DNA sequence file in FASTA format")
        starting_parser.add_argument("-p", "--preserve_headers", default=False, action="store_true",
                                     help="(Optional) Do not modify headers")
        starting_parser.add_argument("-n", "--not_large_index", default=False, action="store_true",
                                     help="(Optional) Avoid usage of a 'large' index, enables sequence splitting if sequence is too large")
        starting_parser.add_argument("-s", "--size", type=float, default=3.6,
                                     help="(Optional) Sequence chunk size in billion charactersif '--not_large_index' is present, 3.6 by default")
        starting_parser.add_argument("-o", "--output", required=True,
                                     help="Output directory")
        return starting_parser.parse_args()


if __name__ == '__main__':
    mainInitializer = Initializer()
    RefDataArray.compile(input_file=mainInitializer.input_nfasta,
                         output_dir=mainInitializer.output_dir,
                         preserve_headers=mainInitializer.preserve_headers_bool,
                         chop=mainInitializer.not_large_index_bool,
                         chunk_length=mainInitializer.chunk_length)
