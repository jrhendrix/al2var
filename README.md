# Al2var

al2var is a tool to find variants between a reference and either paired-end reads or another genome sequence. 

## Introduction
al2var aligns query data to a refence using either bowtie2 or minimap2 and assess the mappings to determine the number of variants and the overall variant rate. When aligning one genome to another (minimap2 mode), al2var detects the differences between the genomes. When aligning paired-end Illumina reads to a reference assembly of the same sample, the variants can be interpretted as errors in the assembly and be used to estimate the assembly error rate.


## Getting Started
### Requirements
* Python 3.7+
* bcftools v1.9+
* samtools v1.6+
* bowtie2 v2.2.5+
* minimap2 v2.21+


al2var was specifically tested with
1. python v3.7.16, bcftools v1.9, samtools v1.6, bowtie2 v2.2.5, and minimap2 v2.21
2. python v3.10.2, bcftools v1.14, samtools v1.14, bowtie2 v2.5.1, and minimap2 v2.24


### Installation
Clone from GitHub repository. If accessing code from GitHub, it is reccomended but not required to install dependencies in a conda environment.

## Output
al2var will create a directory to organize all of the generated output in various subdirectories.

* Subdirectory `vcf/` harbors two files. The file ending in `000.vcf` is a sorted `vcf` that contains a record for each positon in the reference sequence regardless of whether the query sequence had the same or a different genotype. If using the minimap2 mode, the genotype column will reflect the genotype of the query sequence at the mapped position. If using the bowtie2 mode, the genotype will reflect the consensus of the reads that mapped to that given position. The second file ending in `var.vcf` contains only the variants between the reference and the query sequence or input reads.
* `report.txt` file reports on three stats: the alignment rate of the query to the reference, the number of variants between the input, and the variant rate. The variant rate is calculated from the reference length and normalized to 100kb. 
* Subdirectory `bam/` contains the alignmet information in bam format. These files can be viewed in an interactive genome browser such as IGV. These files can be automatically deleted during runtime using the `--cleanup` flag.
* The subdirectories `indexes/`, `mpileup/`, `sam/`, and `unconc/` contain various intermediary files that may be of interest to the user. Some of these files can be large (i.e. `.sam` files) and can be automatically deleted during runtime using the `--cleanup` flag. 


## Usage
Regardless of mode, al2var requires a reference sequence in `fasta` format (any extension is acceptable).

### Paired-end Mode
When using paired-end reads, al2var aligns the reads to the reference using bowtie2. The paired-end reads must be in `fastq` format and can be zipped or not.

Basic usage:
```
al2var bowtie2 -r reference.fasta -1 illumina_raeds_pair1.fastq.gz -2 illumina_reads_pair2.fastq.gz
```

By default, the bowtie2 alignment step will use a seed sequence of 22bp and will not allow any mismatches in the seed sequence. The lenth of the seed sequence can be modified with the `-l` or `--length_seed` flag and the number of permitted mismatches can be increased to 1 using the `-n` or `--num_mismatch` flag. Increasing the number of possible mismatches in the seed sequence increases sensitivity but at the expense of runntime and it it is not recommended to allow for more than 1 mismatch.

To use a seed sequence of 30bp and allow 1 mismatch, use the following command:
```
al2var bowtie2 -r reference.fasta -1 illumina_raeds_pair1.fastq.gz -2 illumina_reads_pair2.fastq.gz --length_seed 30 --num_mismatch 1
```

### Query sequence Mode
When using a second genome, al2var aligns the entire sequence to the reference using minimap2. This second genome sequence must be in `fasta` format (with any extension) and can represent a genome assembly with any number of contigs or an excerpt of a genome.

Basic usage:
```
al2var minimap2 -r reference.fasta -q query.fasta
```

### General usage
By default al2var creates an output file named `out_al2var/` in the current working directory to organize the output files. The name of this file can be changed using the `-o` or `--output_directory` flag and the location can be changed using the `-p` or `--output_path` flag. The prefix of the report file can be specified using the `-s` or `--savename` flag. 

Since the intermediary files can be quite large and may not be of use to the user, al2var provides the user with two options for removing intermediary files during runtime. Use of the `-c` or `--cleanup` flag will remove all data within the directories `bam/`, `indexes/`, `mpileup/`, `sam/`, and `unconc/`. This action will leave the `al2var.log` and `report.txt` files along with the `vcf/` subdirectory. Additional use of the `-b` or `--keep_bam` flags will keep the `bam/` subdirectory rather than deleting it.

