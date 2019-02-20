#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
from modules.Utilities import Utilities


class RefDataLine:
    """
    The class describes a single REFDATA dictionary
    """

    def __init__(self, parsed_dictionary: dict):
        self._nfasta = parsed_dictionary["reference_nfasta"]
        self.db_name = parsed_dictionary.get("alias")
        if not self.db_name:
            self.db_name = Utilities.filename_only(self._nfasta)
        self._reference_mask = Utilities.ends_with_slash(os.path.dirname(os.path.realpath(self._nfasta))) + self.db_name
        self.bowtie_index_mask = parsed_dictionary["ebwt_mask"]
        self.bowtie2_index_mask = parsed_dictionary["bt2_mask"]
        self.samtools_index_file = parsed_dictionary["fai"]
        self.bedtools_genome_file = parsed_dictionary["genome"]
        self.annotation_file = parsed_dictionary["annotation"]

    @staticmethod
    def fill_dict(nfasta_file: str):
        mask = Utilities.ends_with_slash(os.path.dirname(os.path.realpath(nfasta_file))) + Utilities.filename_only(nfasta_file)
        d = {"ebwt_mask": "{}_colorspace".format(mask),
             "bt2_mask": "{}_bowtie2".format(mask),
             "fai": "{}_samtools.fai".format(mask),
             "genome": "{}_samtools.genome".format(mask),
             "annotation": "{}_annotation.tsv".format(mask)}
        return d

    def __bowtie_build(self):
        s = subprocess.getoutput("bowtie-build -C {a} {b}".format(a=self._nfasta, b=self.bowtie_index_mask))
        Utilities.dump_string(string=s, file="{}_bowtie-build.log".format(self._reference_mask))
        print("Created bowtie colorspace index with mask: '{}'".format(self.bowtie_index_mask))

    def __bowtie2_build(self):
        s = subprocess.getoutput("bowtie2-build {a} {b}".format(a=self._nfasta, b=self.bowtie2_index_mask))
        Utilities.dump_string(string=s, file="{}_bowtie2-build.log".format(self._reference_mask))
        print("Created bowtie2 index with mask: '{}'".format(self.bowtie2_index_mask))

    def __samtools_faidx(self):
        s = subprocess.getoutput("samtools faidx {}".format(self._nfasta))
        Utilities.dump_string(string=s, file="{}_samtools_faidx.log".format(self._reference_mask))
        os.rename("{}.fai".format(self._nfasta), self.samtools_index_file)
        print("Created SAMTools FAI file: '{}'".format(self.samtools_index_file))
        self.___fai2genome()

    def ___fai2genome(self):
        """Process existing fasta index, depends from 'samtools_faidx' function"""
        def ____parse_fai_line(split_line: list):
            if len(split_line) >= 2:
                return split_line[:2]
            print("Bad FAI file line: {}".format("\t".join(split_line)))

        fai_2d_array = Utilities.load_2d_array("{}_samtools.fai".format(self._reference_mask))
        genome_2d_array = []
        for line in fai_2d_array:
            genome_2d_array.append(____parse_fai_line(line))
        out = "{}_samtools.genome".format(self._reference_mask)
        Utilities.dump_2d_array(array=Utilities.remove_empty_values(genome_2d_array), file=out)
        print("Created BEDTools genome index: '{}'".format(out))

    def index(self):
        import threading
        """Process chopped or whole FASTA"""
        # Initializing threads
        t1 = threading.Thread(target=self.__bowtie_build)
        t2 = threading.Thread(target=self.__bowtie2_build)
        t3 = threading.Thread(target=self.__samtools_faidx)
        # Starting threads
        t1.start()
        t2.start()
        t3.start()
        # Joining threads
        t1.join()
        t2.join()
        t3.join()
