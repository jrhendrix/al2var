# Al2var

A tool to calculate the error rate of a bacterial assembly using paired-end Illumina reads. 

## Introduction


## Getting Started
### Requirements
* Python
* bowtie2
* bcftools
* samtools

### Installation
Clone from GitHub repository. If accessing code from GitHub, it is reccomended to install dependencies in a conda environment.

## Usage
As input, Al2var takes a genome assembly file in (.fasta format) and paired end reads (.fastq format).
Example:
'''al2var -r reference.fasta -1 illumina_raeds_pair1.fastq.gz -2 illumina_reads_pair2.fastq.gz'''

Example allowing 1 mismatch in aligning seed sequence:

Example using a seed sequence of 18 bases:

Example with cleanup enabled

Extensive Usage:
```al2var -r reference.fasta -1 illumina_raeds_pair1.fastq.gz -2 illumina_reads_pair2.fastq.gz -p path/to/output -o output_name -s report_name'''


## Output Files
