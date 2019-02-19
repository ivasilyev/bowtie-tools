#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import logging
from modules.Utilities import Utilities


class SampleDataLine(object):
    """
    The class describes single sample data file row
    """
    def __init__(self, single_sampledata_row):
        body_list = Utilities.remove_empty_values(
            [i.strip() for i in re.sub("[\r\n]", "", single_sampledata_row).split("\t")])
        if len(body_list) < 2:
            Utilities.log_and_raise(
                "Failed to parse sample data row (not enough columns): {}".format(single_sampledata_row))
        self.name = body_list[0]
        self.raw_reads_files_list = body_list[1:]
        if len(self.raw_reads_files_list) > 2:
            logging.warning("Only up to two input files are supported for alignment, using first two values. "
                            "Given sample data row: '{}'".format(single_sampledata_row))
            self.raw_reads_files_list = self.raw_reads_files_list[:2]
        for file in self.raw_reads_files_list:
            if not os.path.isfile(file):
                logging.warning("Not found the raw reads file: '{}'".format(file))

    def __len__(self):
        """
        :return:
        Number of existing files
        """
        return sum([os.path.isfile(i) for i in self.raw_reads_files_list])

    def export(self):
        return "\t".join([self.name] + self.raw_reads_files_list)
