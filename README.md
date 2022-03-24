# bowtie-tools
Automation & scaling scripts for [bowtie](http://bowtie-bio.sourceforge.net)/[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/).

## nBee
Main script to perform single alignment & postprocessing for every colorspace or non-colorspace sequence reads file specified in the input file list.
The coverage table contains the following columns:
* `reference_id`: FASTA header of the reference sequence
* `id_bp`: Reference length (base pairs)
* `id_maximal_coverage_depth`: How many reads mapped simultaneously to any reference part
* `id_coverage_breadth`: The coverage breadth, the length of the reference covered (not additive)
* `id_mapped_bp`: The length of the reference covered (additive)
* `id_mapped_reads`: Number of reads mapped with maximal score on the reference
* `id_unmapped_reads`: Number of reads not mapped with any score on the reference
* `sample_total_reads`: The number of the sample reads, regardless of mapping
* `sample_mapped_reads`: The number of the sample reads mapped on the reference
* `sample_total_bp`: The sum of the sample reads, regardless of mapping
* `sample_mapped_bp`: The sum of the sample reads mapped on the reference
* `id_coverage_breadth_to_id_bp`: `id_coverage_breadth` / `id_bp`
* `id_total_relative_abundance`: 10<sup>12</sup> * `id_mapped_bp` / (`id_bp` * `sample_total_bp`)
* `id_mapped_relative_abundance`: 10<sup>12</sup> * `id_mapped_bp` / (`id_bp` * `sample_mapped_bp`)
* `sample_average_total_reads_bp`: `sample_total_reads` / `sample_total_bp`
* `sample_average_mapped_reads_bp`: `sample_mapped_reads` / `sample_total_bp`
* `sample_mapped_reads_to_total_reads`: `sample_mapped_reads` / `sample_total_reads`
* `id_mapped_reads_per_million_sample_total_reads`: 10<sup>6</sup> * `id_mapped_reads` / `sample_total_reads`
* `id_mapped_reads_per_million_sample_mapped_reads`: 10<sup>6</sup> * `id_mapped_reads` / `sample_mapped_reads`
* `id_mapped_reads_per_kbp_per_million_sample_total_reads`: 10<sup>9</sup> * `id_mapped_reads` / (`sample_total_reads` * `id_bp`)
* `id_mapped_reads_per_kbp_per_million_sample_mapped_reads`: 10<sup>9</sup> * `id_mapped_reads` / (`sample_mapped_reads` * `id_bp`)

## sam2coverage
This tool automatically converts, sorts, indexes and extracts coverage from SAM file.

## cook_the_reference
This tool automatically fixes, cuts and indexes DNA sequences in FASTA format. Also makes the "refdata" reference files linker.

Required software: bowtie, bowtie2, samtools, bedtools.
