# HmnRandomRead

[![Github Version](https://img.shields.io/github/v/release/guillaume-gricourt/HmnRandomRead?display_name=tag&sort=semver)](version)  [![Conda Version](https://img.shields.io/conda/vn/bioconda/hmnrandomread.svg)](https://anaconda.org/bioconda/hmnrandomread)  
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black) [![GitHub Super-Linter](https://github.com/guillaume-gricourt/HmnRandomRead/workflows/Tests/badge.svg)](https://github.com/marketplace/actions/super-linter)  
[![DOI](https://zenodo.org/badge/581311021.svg)](https://zenodo.org/badge/latestdoi/581311021)

## Features

* use one or more references
* control adaptaters and the insert size
* adjust the exact number of sequences
* adapt the error model coming from your sequencer
* eventually add SNPs to introduce diversity
* generate around 1000 sequences by second

## Install

```sh
conda install -c bioconda hmnrandomread
```

## Use

```sh
HmnRandomRead \
    -r/--reference <string, required><double, required> \
    -r1/--read-forward <string, required> \
    -r2/--read-reverse <string, required> \

    -lengthReads/--length-reads <int, optional, 150> \
    -meanInsert/--mean-insert-size <int, optional, 500> \
    -stdInsert/--std-insert-size <int, optional, 50> \

    -profileDiversity/--profile-diversity <string, optional> \
    -profileError/--profile-error <string, optional> \
    -profileErrorId/--profile-error-id <string, optional> \
    -s/--seed <int, optional, 0>
```

### Reference

Use one or more FASTA file used as reference sequence.
Indicate also the number of sequence to generate for each reference.

### Output

`read-forward` and `read-reverse` are required, gzip compressed or not.

### Sequencing size

`length-reads`: the size of the library as sequenced by the sequencer

### Library size

`mean-insert-size` and `std-insert-size`: the gaussian parameters to represent the fragment size.

### Profile diversity

`profile-diversity` a CSV file, comma separated, with header:
* identifier: ID of the fasta file
* Mutation Rate: probability to change the sequence
* Indel Fraction: rate of indel compare to single mutation
* Indel Extend: probability to extend the indel at each base added
* Maximum Insertion Size: maximal size of insertion
The header is mandatory.

### Profile error

`profile-error` a CSV file, comma separated, with header:
* identifier: an ID choose by `profile-error-id`
* sequencer: name of the sequencer
* flowcell: kind of flowcell
* version: the kit version.
* strand: `forward` or `reverse`
* cycles total: by strand
* error by cycle: rate of error by cycle, semi-colon separated. Equal to the number of `cycles total`.
The header is mandatory.

## Test

*pytest* is required:
```sh
make test
```

## Built with these main libraries

* [htslib](https://github.com/samtools/htslib) - Write FASTQ files

## Authors

* **Guillaume Gricourt**
