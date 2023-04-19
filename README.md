# Al2var

al2var identifies variants and caclulates the variant rate between a bacterial genome sequence and either paired-end reads or another genome sequence. al2var can also estimate the number of errors in the assembly.

## Introduction
The primary purpose of al2var is to align either paired-end reads or a genome assembly to a given sequence. Al2var aligns the input using either Bowtie2 (paired-end reads) or minimap2 (genome sequence) then interprets the mappings to determine the number of variants and the variant rate. 

Al2var can also be used to estimate the number of errors in an assembly. It is recommended that high quality genome assemblies have fewer than 1 error per 100,000 bases (Chain, 2009).

### Reference-based error estimation
In reference-based error analysis, a sample is compared to a trusted representative of the same species to identify variants. Either genome could be used as the reference or the query for al2var, but the error rate will be normalized to the length of the reference. This application is not recommended for species that are diverse or that lack a sufficiently close and trusted representative sequence because true variants may be misinterpreted as errors.

### _de novo_ error estimation
Errors could also be estimated without a species representative using the al2var bowtie2 mode, where paired-end reads are aligned to a given genome assembly and the variants are interpreted as errors. The paired-end approach to estimating errors requires that the genome and paired-end reads originated from the same sample. In addition, this method assumes that the consensus of multiple Illumina reads at a given location are correct. Though highly accurate, Illumina reads are not impervious to errors so this method is meant to provide an error estimate. To reduce the number of false positives, trim low quality ends off of Illumina reads and use higher read coverage. 


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
al2var is available on [PyPI](https://pypi.org/project/roprokka/) and can be installed using pip.

```pip install al2var```

OR clone from GitHub repository. If accessing code from GitHub, it is recommended but not required to install dependencies in a conda environment.

## Output
al2var will create a directory to organize all of the generated output in various subdirectories.

* Subdirectory `vcf/` harbors two files. The file ending in `000.vcf` is a sorted `vcf` that contains a record for each position in the reference sequence regardless of whether the query sequence had the same or a different genotype. If using the minimap2 mode, the genotype column will reflect the genotype of the query sequence at the mapped position. If using the bowtie2 mode, the genotype will reflect the consensus of the reads that mapped to that given position. The second file ending in `var.vcf` contains only the variants between the reference and the query sequence or input reads.
* `report.txt` file reports on three stats: the alignment rate of the query to the reference, the number of variants between the input, and the variant rate. The variant rate is calculated from the reference length and normalized to 100kb. 
* Subdirectory `bam/` contains the alignment information in bam format. These files can be viewed in an interactive genome browser such as IGV. These files can be automatically deleted during runtime using the `--cleanup` flag.
* The subdirectories `indexes/`, `mpileup/`, `sam/`, and `unconc/` contain various intermediary files that may be of interest to the user. Some of these files can be large (i.e. `.sam` files) and can be automatically deleted during runtime using the `--cleanup` flag. 


## Usage
Regardless of mode, al2var requires a reference sequence in `fasta` format (any extension is acceptable).

### Paired-end Mode
When using paired-end reads, al2var aligns the reads to the reference using bowtie2. The paired-end reads must be in `fastq` format and can be zipped or not.

Basic usage:
```
al2var bowtie2 -r reference.fasta -1 illumina_raeds_pair1.fastq.gz -2 illumina_reads_pair2.fastq.gz
```

By default, the bowtie2 alignment step will use a seed sequence of 22bp and will not allow any mismatches in the seed sequence. The lenth of the seed sequence can be modified with the `-l` or `--length_seed` flag and the number of permitted mismatches can be increased to 1 using the `-n` or `--num_mismatch` flag. Increasing the number of possible mismatches in the seed sequence increases sensitivity but at the expense of runntime and it is not recommended to allow for more than 1 mismatch.

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

## References
Chain, P. S. G., et al. "Genome project standards in a new era of sequencing." Science 326.5950 (2009): 236-237.
