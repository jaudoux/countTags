# countTags - counting a set of k-mers into a set of FASTQ files

## Installation

1. Clone this repository : `git clone https://github.com/jaudoux/countTags.git`
2. Compile the software : `make`
3. Place the binary somewhere acessible in your `$PATH`

## Usage

First you need to create a FASTA file containing the tags you want to quantify.
All tags must have the same length, and the maximum authorize lenght is 32bp.
K-mer length must be provided to countTags using the `-k INT` option.

For example :

`countTags -k 31 file.fa file1.fastq.gz file2.fastq.gz`

### Managing Paired-End files

For paired-end files, you need to apply a post treatment to merge counts
produced separatly for each fastq file.

### Managing stranded files

For stranded fastq, you need to set the `--stranded` option. If you are using PE
stranded FASTQ files, you need to manually reverse one of the pair according to
the protocol orientation.
