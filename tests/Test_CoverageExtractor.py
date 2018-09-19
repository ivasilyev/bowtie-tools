#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
echo Deploy worker image
export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it ${IMG} python3

"""

import pandas as pd

stats_dict = {"sample_total_reads": 999, "sample_mapped_reads": 888, "sample_total_bp": 777, "sample_mapped_bp": 666}

self__stacked_coverages_df = pd.read_table("/data1/bio/sandbox/nBee/igc/Statistics/1-14_igc_v2014.03_chunk_1_pos_bp.txt", header=0, sep='\t')
reference_df = pd.read_table("/data/reference/IGC/igc_v2014.03/index/igc_v2014.03_chunk_1_samtools.genome", header='infer', sep='\t', names=['reference_id', 'id_bp'])

genomes_coverages_df = pd.merge(reference_df, self__stacked_coverages_df.loc[:, [i for i in list(self__stacked_coverages_df) if i != "id_bp"]], on="reference_id", how='left')

genomes_coverages_df = genomes_coverages_df.loc[genomes_coverages_df['reference_id'] != 'genome']

genomes_coverages_df["id_total_relative_abundance"] = genomes_coverages_df.loc[:, "id_mapped_bp"] / (genomes_coverages_df.loc[:, "id_bp"] * stats_dict["sample_total_bp"])
genomes_coverages_df["id_mapped_relative_abundance"] = genomes_coverages_df.loc[:, "id_mapped_bp"] / (genomes_coverages_df.loc[:, "id_bp"] * stats_dict["sample_mapped_bp"])
genomes_coverages_df["sample_total_reads"] = stats_dict["sample_total_reads"]
genomes_coverages_df["sample_mapped_reads"] = stats_dict["sample_mapped_reads"]
genomes_coverages_df["sample_total_bp"] = stats_dict["sample_total_bp"]
genomes_coverages_df["sample_mapped_bp"] = stats_dict["sample_mapped_bp"]
genomes_coverages_df["sample_average_total_reads_bp"] = float(stats_dict["sample_total_reads"]) / float(stats_dict["sample_total_bp"])
genomes_coverages_df["sample_average_mapped_reads_bp"] = float(stats_dict["sample_mapped_reads"]) / float(stats_dict["sample_total_bp"])
genomes_coverages_df["sample_mapped_reads_to_total_reads"] = float(stats_dict["sample_mapped_reads"]) / float(stats_dict["sample_total_reads"])

self__samtools_idxstats_df = pd.read_table("/data1/bio/sandbox/nBee/igc/Statistics/1-14_igc_v2014.03_chunk_2_idxstats.txt", header='infer', sep='\t', names=["reference_id", "id_bp", "id_mapped_reads", "id_unmapped_reads"])

genomes_coverages_df = pd.merge(genomes_coverages_df, self__samtools_idxstats_df.loc[:, [i for i in list(self__samtools_idxstats_df) if i != "id_bp"]], on="reference_id", how='left')
