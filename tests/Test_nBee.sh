#!/usr/bin/env bash
# -*- coding: utf-8 -*-

echo Deploy worker image
export IMG=ivasilyev/bwt_filtering_pipeline_worker:latest && \
docker pull ${IMG} && \
docker run --rm -v /data:/data -v /data1:/data1 -v /data2:/data2 -it ${IMG} bash

echo Edit RefDataArray.py
nano /home/docker/scripts/modules/RefDataArray.py

echo Edit Aligner.py
nano /home/docker/scripts/modules/Aligner.py

echo Edit nBee.py
nano /home/docker/scripts/nBee.py

echo Test nBee
python3 /home/docker/scripts/nBee.py -i /data1/bio/projects/tgrigoreva/ecoli_plates_suspension/raw.sampledata -r /data/reference/IGC/igc_v2014.03/index/igc_v2014.03_refdata.json -m no_hg19 -t half -o /data2/bio/Metagenomes/IGC/test

echo Low-level pipeline processing
bowtie2 --un /data2/bio/Metagenomes/IGC/test/unmapped -x /data/reference/IGC/igc_v2014.03/index/igc_v2014.03_chunk_1_bowtie2 --threads 20 --very-sensitive -S -1 /data2/bio/ecoli_komfi/raw_reads/1-39-0-5_S84_R1_001.fastq.gz -2 /data2/bio/ecoli_komfi/raw_reads/1-39-0-5_S84_R2_001.fastq.gz | \
samtools view -bu -@ 20 - | \
samtools sort - -@ 20 -o /data2/bio/Metagenomes/IGC/test/mapped.sam

echo Low-level no-pipeline processing
bowtie2 --threads 20 --very-sensitive --un-conc /data2/bio/Metagenomes/IGC/test/unmapped.fa -x /data/reference/IGC/igc_v2014.03/index/igc_v2014.03_chunk_1_bowtie2 -1 /data2/bio/ecoli_komfi/raw_reads/1-39-0-5_S84_R1_001.fastq.gz -2 /data2/bio/ecoli_komfi/raw_reads/1-39-0-5_S84_R2_001.fastq.gz -S /data2/bio/Metagenomes/IGC/test/mapped.sam
samtools view -bu -@ 1 /data2/bio/Metagenomes/IGC/test/mapped.sam | \
samtools sort - -@ 1 -o /data2/bio/Metagenomes/IGC/test/mapped_sorted.bam
genomeCoverageBed -ibam /data2/bio/Metagenomes/IGC/test/mapped_sorted.bam > /data2/bio/Metagenomes/IGC/test/mapped_gc.tsv

echo Low-level no-pipeline colorspace processing
bowtie -f -C -t -v 3 -k 1 --threads 20 --un /data2/bio/Metagenomes/IGC/test/col_unmapped.csfasta /data/reference/IGC/igc_v2014.03/index/igc_v2014.03_chunk_1_colorspace /data2/bio/Metagenomes/HG19/Non-mapped_reads/107VZK_no_hg19.csfasta -S /data2/bio/Metagenomes/IGC/test/col_mapped.sam
samtools view -bu -@ 1 /data2/bio/Metagenomes/IGC/test/col_mapped.sam | \
samtools sort - -@ 1 -o /data2/bio/Metagenomes/IGC/test/col_mapped_sorted.bam
genomeCoverageBed -ibam /data2/bio/Metagenomes/IGC/test/col_mapped_sorted.bam > /data2/bio/Metagenomes/IGC/test/col_mapped_gc.tsv

echo Low-level no-pipeline processing with STOUT and STDERR split
bowtie2 --threads 20 --very-sensitive --un-conc /data2/bio/Metagenomes/IGC/test/unmapped.fa -x /data/reference/IGC/igc_v2014.03/index/igc_v2014.03_chunk_1_bowtie2 -1 /data2/bio/ecoli_komfi/raw_reads/1-39-0-5_S84_R1_001.fastq.gz -2 /data2/bio/ecoli_komfi/raw_reads/1-39-0-5_S84_R2_001.fastq.gz 1> /data2/bio/Metagenomes/IGC/test/mapped.sam 2> /data2/bio/Metagenomes/IGC/test/mapped.log | \

echo Low-level pipeline processing with STOUT and STDERR split
bowtie2 --threads 20 --very-sensitive --un-conc /data2/bio/Metagenomes/IGC/test/unmapped.fa -x /data/reference/IGC/igc_v2014.03/index/igc_v2014.03_chunk_1_bowtie2 -1 /data2/bio/ecoli_komfi/raw_reads/1-39-0-5_S84_R1_001.fastq.gz -2 /data2/bio/ecoli_komfi/raw_reads/1-39-0-5_S84_R2_001.fastq.gz 2> /data2/bio/Metagenomes/IGC/test/mapped.log | \
samtools view - -bu -@ 20 | \
samtools sort - -@ 20 -o /data2/bio/Metagenomes/IGC/test/mapped_sorted.bam
genomeCoverageBed -ibam /data2/bio/Metagenomes/IGC/test/mapped_sorted.bam > /data2/bio/Metagenomes/IGC/test/mapped_gc.tsv

echo Low-level pipeline colorspace processing
bowtie -f -C -t -v 3 -k 1 --threads 20 --un /data2/bio/Metagenomes/IGC/test/col_unmapped.csfasta /data/reference/IGC/igc_v2014.03/index/igc_v2014.03_chunk_1_colorspace /data2/bio/Metagenomes/HG19/Non-mapped_reads/107VZK_no_hg19.csfasta -S | \
samtools view - -bu -@ 20 | \
samtools sort - -@ 20 -o /data2/bio/Metagenomes/IGC/test/col_mapped_sorted.bam
genomeCoverageBed -ibam /data2/bio/Metagenomes/IGC/test/col_mapped_sorted.bam > /data2/bio/Metagenomes/IGC/test/col_mapped_gc.tsv

echo Low-level pipeline colorspace processing with STOUT and STDERR split
bowtie -f -C -t -v 3 -k 1 --threads 20 --un /data2/bio/Metagenomes/IGC/test/col_unmapped.csfasta /data/reference/IGC/igc_v2014.03/index/igc_v2014.03_chunk_1_colorspace /data2/bio/Metagenomes/HG19/Non-mapped_reads/107VZK_no_hg19.csfasta -S 2> /data2/bio/Metagenomes/IGC/test/col_mapped.log | \
samtools view - -bu -@ 20 | \
samtools sort - -@ 20 -o /data2/bio/Metagenomes/IGC/test/col_mapped_sorted.bam
genomeCoverageBed -ibam /data2/bio/Metagenomes/IGC/test/col_mapped_sorted.bam > /data2/bio/Metagenomes/IGC/test/col_mapped_gc.tsv
