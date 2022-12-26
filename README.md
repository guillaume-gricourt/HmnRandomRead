# HmnRandomRead

[![Github Version](https://img.shields.io/github/v/release/guillaume-gricourt/HmnRandomRead?display_name=tag&sort=semver)](version)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black) [![GitHub Super-Linter](https://github.com/guillaume-gricourt/HmnRandomRead/workflows/Tests/badge.svg)](https://github.com/marketplace/actions/super-linter)
[![License](https://img.shields.io/github/license/guillaume-gricourt/HmnRandomRead)](license)

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
./HmnRandomRead \
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

## Test

*pytest* is required:
```sh
make test
```

## Built with these main libraries

* [htslib](https://github.com/samtools/htslib) - Write FASTQ files

## Authors

* **Guillaume Gricourt**
