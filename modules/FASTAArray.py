#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import re
from modules.FASTALine import FASTALine
from modules.Utilities import Utilities


class FASTAArray(object):
    """
    This class keeps FASTALine objects list
    """
    def __init__(self, parsed_fastas_list: list):
        self._parsed_fastas_list = Utilities.remove_empty_values(parsed_fastas_list)
        self._parsed_fastas_list.sort(key=len, reverse=True)
        self._annotations_2d_array = [["reference_id", "id_bp"]]
        for fasta in self._parsed_fastas_list:
            self._annotations_2d_array.append([fasta.header, str(len(fasta))])

    def __len__(self):
        return sum([len(i) for i in self._parsed_fastas_list])

    def get_fastas_list(self):
        return self._parsed_fastas_list

    def get_headers_list(self):
        return [i.header for i in self._parsed_fastas_list]

    def get_total_sequence(self):
        return "\n".join([i.to_str() for i in self._parsed_fastas_list]) + "\n"

    def get_total_length(self):
        return len(self.get_total_sequence())

    def dump_annotation(self, file: str):
        Utilities.dump_2d_array(array=self._annotations_2d_array, file=file)
        print("Created annotation file: '{}'".format(file))

    def dump_fastas(self, file: str):
        Utilities.dump_list(lst=[i.to_str() for i in self._parsed_fastas_list], file=file)
        print("Created FASTA file: '{}'".format(file))

    def _chop_sequences(self, input_chunk_length: int = int(3.6 * 10 ** 9)):
        default_max_chunk_length = int(3.6 * 10 ** 9)  # Bowtie limit
        if input_chunk_length > default_max_chunk_length:
            print("The chunk size must not be larger than 3.6 billion characters!")
            input_chunk_length = default_max_chunk_length
        print("Using chunk size: {}".format(input_chunk_length))
        chunks_counter = 0
        chunk_fastas_list = []
        output_chunk_length = 0
        output_dict = {}
        for fasta in self._parsed_fastas_list:
            current_chunk_length = output_chunk_length + fasta.get_total_length()
            if current_chunk_length < input_chunk_length:
                chunk_fastas_list.append(fasta)
                output_chunk_length = current_chunk_length
            else:
                chunks_counter += 1
                output_dict["chunk_{}".format(chunks_counter)] = FASTAArray(chunk_fastas_list)
                print("Processed chunk {}".format(chunks_counter))
                chunk_fastas_list = [fasta]
                output_chunk_length = fasta.get_total_length()
        print("Total chunks created: {}".format(len(output_dict)))
        return output_dict

    @staticmethod
    def __process_header(header: str):
        return re.sub('_+', '_', re.sub('\W+', '_', header)).strip('_') + '\n'

    def _fix_headers(self):
        headers_zfill_number = len(str(len(self._parsed_fastas_list)))
        print("Managing {a} headers with {b}-digit filling. All changes would be saved into annotation file".format(a=len(self._parsed_fastas_list), b=headers_zfill_number))
        fixed_fastas_list = []
        annotations_2d_array = [["reference_id", "former_id", "id_bp"]]
        headers_counter = 0
        for fasta in self._parsed_fastas_list:
            headers_counter += 1
            old_header = fasta.header
            fasta.set_header("ID" + str(headers_counter).zfill(headers_zfill_number))
            fixed_fastas_list.append(fasta)
            annotations_2d_array.append([fasta.header, old_header, str(len(fasta))])
        print("{} headers have been processed".format(headers_counter))
        self._parsed_fastas_list = fixed_fastas_list
        self._annotations_2d_array = annotations_2d_array

    @staticmethod
    def parse(string: str):
        from modules.Utilities import Utilities
        string = re.sub("[\r\n]+", "\n", string)
        q = [">{}".format(j) for j in Utilities.remove_empty_values([i.strip() for i in re.split("^>", string)])]
        return FASTAArray([FASTALine(i) for i in q])

    @staticmethod
    def prepare_nfasta_for_indexing(input_file: str, output_dir: str, preserve_headers: bool = False, chop: bool = False, chunk_length: int = int(3.6 * 10 ** 9)):
        array = FASTAArray.parse(Utilities.load_string(input_file))
        if not preserve_headers:
            array._fix_headers()
        output_dir = Utilities.ends_with_slash(output_dir)
        os.makedirs(output_dir, exist_ok=True)
        output_file_mask = (output_dir + Utilities.filename_only(input_file))
        annotation_file = "{}_annotation.tsv".format(output_file_mask)
        array.dump_annotation(annotation_file)
        arrays_dict = {"{}.fasta".format(output_file_mask): array}
        if chop and array.get_total_length() >= chunk_length:
            print("Too large reference nFASTA: '{}'. Splitting sequences".format(input_file))
            arrays_dict = array._chop_sequences(chunk_length)
            arrays_dict = {"{a}{i}.fasta".format(a=output_file_mask, i=i): arrays_dict[i] for i in arrays_dict}
        refdatas_dict = {}
        counter = 0
        for chunk_file in arrays_dict:
            counter += 1
            arrays_dict[chunk_file].dump_fastas(chunk_file)
            refdatas_dict["sequence_{}".format(counter)] = {"reference_nfasta": chunk_file, "annotation": annotation_file}
        print("FASTA files created: {}".format(counter))
        return refdatas_dict
