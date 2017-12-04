# bowtie-tools
Automation & scaling scripts for [bowtie](http://bowtie-bio.sourceforge.net)/[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/).

## nBee
Main script to perform single alignment & postprocessing for every colorspace or non-colorspace sequence reads file specified in the input file list.

## sam2coverage
This tool automatically converts, sorts, indexes and extracts coverage from SAM file.

## cook_the_reference
This tool automatically fixes, cuts and indexes DNA sequences in FASTA format. Also makes the "refdata" reference files linker.

## refdata2schedule
This script will align and extract coverage from colorspace or non-colorspace reads mapped on reference. It requires a "refdata" file.

Required software: bowtie, bowtie2, samtools, bedtools.
