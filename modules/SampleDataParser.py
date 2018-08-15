#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from modules.Utilities import Utilities
from modules.SampleDataLine import SampleDataLine


class SampleDataParser:
    """
    The class analyzes text files containing 2 or 3 tab-delimited columns
    The first column contains sample name
    Second and third (if present) columns contain paths to sample files
    """
    def __init__(self, sampledata_file_name):
        if not os.path.isfile(sampledata_file_name):
            Utilities.log_and_raise("Sample data linker file not found: {}".format(sampledata_file_name))
        self._sampledatas_list = []
        with open(sampledata_file_name, "r", encoding="utf-8") as f:
            for r in f:
                r = r.strip()
                if len(r) > 0:
                    try:
                        self._sampledatas_list.append(SampleDataLine(r))
                    except ValueError:
                        continue
        self._sampledatas_list = Utilities.remove_empty_values(self._sampledatas_list)

    def __len__(self):
        return len(self._sampledatas_list)

    def get_parsed_list(self):
        return self._sampledatas_list
