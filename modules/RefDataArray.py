#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import traceback
from modules.Utilities import Utilities
from modules.RefDataLine import RefDataLine


class RefDataArray:
    """
    The class analyzes JSONs or REFDATA, tab-delimited text files containing paths to reference index files.
    REFDATA format:
    <Reference DNA fasta chunk 1> <Bowtie index mask> <Bowtie2 index mask> <SamTools index> <BedTools genome length> <Annotation file>
    <Reference DNA fasta chunk 2> <Bowtie index mask> <Bowtie2 index mask> <SamTools index> <BedTools genome length> <Annotation file>
    ...
    JSON keys (corresponding to above):
    {'sequence_1': {'reference_nfasta': '', 'ebwt_mask': '', 'bt2_mask': '', 'fai': '', 'genome': '', 'annotation': ''},
    'sequence_2': {'reference_nfasta': '', 'ebwt_mask': '', 'bt2_mask': '', 'fai': '', 'genome': '', 'annotation': ''},
    ...}
    """
    _refdata_keys_list = ["reference_nfasta", "ebwt_mask", "bt2_mask", "fai", "genome", "annotation"]

    def __init__(self, input_dict: dict, verify: bool = False):
        self._body_dict = input_dict
        if verify:
            self._verify_json_refdata(self._body_dict)
        self.refdata_lines_dict = {k: RefDataLine(self._body_dict[k]) for k in self._body_dict}

    def _verify_json_refdata(self, d: dict):
        if len(d) == 0:
            raise ValueError("Empty sample data!")
        for k1 in d:
            if len([i for i in list(d) if i == k1]) > 1:
                raise ValueError("Repeating key: {}".format(k1))
            for k_r in [i for i in self._refdata_keys_list if i != "reference_nfasta"]:
                if k_r not in list(d[k1]):
                    raise ValueError("Missing value for the key: {}".format(k_r))
            for k2 in d[k1]:
                f = d[k1][k2]
                if k2 not in ["reference_nfasta", "ebwt_mask", "bt2_mask"] and not os.path.isfile(f):
                    raise ValueError("Not found file: '{}', keys: '{}', '{}'".format(f, k1, k2))

    def get_parsed_list(self):
        return [self.refdata_lines_dict[k] for k in self.refdata_lines_dict]

    def get_refdata_line_by_index(self, idx: int):
        return self.get_parsed_list()[idx]

    def get_refdata_line_by_key(self, key: str):
        return self._body_dict.get(key)

    @staticmethod
    def _parse_json_refdata(file_wrapper):
        import json
        return RefDataArray(json.load(file_wrapper))

    @staticmethod
    def _parse_table_refdata(file_wrapper):
        d = {}
        counter = 1
        for line in file_wrapper:
            if len(line.strip()):
                d["sequence_{}".format(counter)] = {k: v for k, v in zip(RefDataArray._refdata_keys_list, [i.strip() for i in line.split("\t")])}
        return RefDataArray(d)

    @staticmethod
    def read(file: str):
        wrapper = open(file=file, mode="r", encoding="utf-8")
        try:
            if file.endswith(".json"):
                return RefDataArray._parse_json_refdata(wrapper)
            else:
                return RefDataArray._parse_table_refdata(wrapper)
        except ValueError:
            traceback.print_exc()
            Utilities.log_and_raise("Bad reference data file: {}".format(file))

    @staticmethod
    def compile(input_file: str, output_dir: str, preserve_headers: bool = False, chop: bool = False, chunk_length: int = int(3.6 * 10 ** 9)):
        import json
        from modules.FASTAArray import FASTAArray
        from modules.RefDataLine import RefDataLine

        output_dir = Utilities.ends_with_slash(output_dir)
        os.makedirs(output_dir, exist_ok=True)
        refdatas_dict = FASTAArray.prepare_nfasta_for_indexing(input_file=input_file, output_dir=output_dir, preserve_headers=preserve_headers, chop=chop, chunk_length=chunk_length)
        output_dict = {}
        for sequence_id in refdatas_dict:
            annotation_dict = refdatas_dict[sequence_id]
            nfasta_file = annotation_dict.get("reference_nfasta")
            if not nfasta_file:
                continue
            indexing_dict = {"alias": Utilities.filename_only(nfasta_file)}
            indexing_dict.update(RefDataLine.fill_dict(nfasta_file))
            indexing_dict.update(annotation_dict)
            print("Processing nFASTA: '{}'".format(nfasta_file))
            refdata = RefDataLine(indexing_dict)
            refdata.index()
            output_dict[sequence_id] = indexing_dict
        output_file = "{a}{b}_refdata.json".format(a=Utilities.ends_with_slash(output_dir), b=Utilities.filename_only(input_file))
        Utilities.dump_string(string=json.dumps(output_dict, sort_keys=False, indent=4) + "\n", file=output_file)
        print("Created reference data linker: '{}'".format(output_file))
        return output_file
