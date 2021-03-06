#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import re


class FASTALine:
    """
    This class is an attempt to apply NCBI standards to single FASTA.
    Consumes one header followed by sequence.
    """
    def __init__(self, single_fasta: str):
        self._body = re.sub("[\r\n]+", "\n", single_fasta).strip()
        if self._body.startswith(">"):
            self.header = re.sub("^>", "", self._body.split("\n")[0].strip())
            self.sequence = re.sub("[^A-Z]", "", "".join([i.strip() for i in self._body.split("\n")[1:]]).upper())
            # Nucleotide sequence has only AT(U)GC letters. However, it may be also protein FASTA.
        else:
            print("Cannot parse the header for sequence: '{}'".format(self._body))
            self.header = ""
            self.sequence = ""

    @staticmethod
    def _slice_str_by_len(string: str, length: int):
        return [string[0 + i:length + i] for i in range(0, len(string), length)]

    def __len__(self) -> int:
        return len(self.sequence)

    def __hash__(self):
        return hash(self.sequence)

    def __eq__(self, other):
        return self.sequence == other.sequence

    def __lt__(self, other):
        return len(self) < len(other)

    def get_total_length(self):
        return len(self.to_str())

    def _export_sequence(self):
        return "\n".join(self._slice_str_by_len(self.sequence, 70))

    def to_dict(self):
        return {self.header: self._export_sequence()}

    def to_str(self):
        return ">{a}\n{b}".format(a=self.header, b=self._export_sequence())

    def set_header(self, header):
        if header:
            if len(header.strip()) > 0:
                self.header = header.strip()
