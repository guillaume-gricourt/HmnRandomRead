# HmnRandomRead

[![Github Version](https://img.shields.io/github/v/release/guillaume-gricourt/HmnRandomRead?display_name=tag&sort=semver)](version)  [![Conda Version](https://img.shields.io/conda/vn/bioconda/hmnrandomread.svg)](https://anaconda.org/bioconda/hmnrandomread)
[![DOI](https://zenodo.org/badge/581311021.svg)](https://zenodo.org/badge/latestdoi/581311021)

## Features

- use one or more references
- control adaptaters and the insert size
- adjust the exact number of sequences
- adapt the error model coming from your sequencer
- eventually add SNPs to introduce diversity
- generate around 1000 sequences by second

## Install

```sh
conda install -c bioconda hmnrandomread
```

### Build from source

```sh
git clone git@github.com:guillaume-gricourt/HmnRandomRead.git
cd HmnRandomRead
cargo build --release
```

The binary is built at `target/release/HmnRandomRead`.

## Use

```sh
HmnRandomRead simulate \
    --input-reference-fasta <string, required><int, optional><string, optional> \
    --output-forward-fastq <string, required> \
    --output-reverse-fastq <string, required> \

    --parameter-length-reads-int <int, optional, 150> \
    --parameter-mean-insert-int <int, optional, 500> \
    --parameter-std-insert-int <int, optional, 50> \

    --input-profile-diversity-csv <string, optional> \
    --input-profile-sequencer-csv <string, optional> \
    --parameter-profile-sequencer-id-str <string, optional> \
    --parameter-seed-int <int, optional, 0>

HmnRandomRead build-profile-sequencer \
    --parameter-id-str <string, required> \
    --input-forward-fastq <string, optional> \
    --input-reverse-fastq <string, optional> \
    --input-bam <string, optional> \
    --output-profile-sequencer-csv <string, required>

HmnRandomRead version
```

### Reference

Use one or more FASTA file used as reference sequence (`--input-reference-fasta`, may be
repeated: `path[,nb_reads[,id_diversity]]`).
Indicate also the number of sequence to generate for each reference.

### Output

`--output-forward-fastq` and `--output-reverse-fastq` are required, gzip compressed.

### Sequencing size

`--parameter-length-reads-int`: the size of the library as sequenced by the sequencer

### Library size

`--parameter-mean-insert-int` and `--parameter-std-insert-int`: the gaussian parameters to represent the fragment size.

### Profile diversity

`--input-profile-diversity-csv` a CSV file, comma separated, with header:
- identifier: ID of the fasta file
- Mutation Rate: probability to change the sequence
- Indel Fraction: rate of indel compare to single mutation
- Indel Extend: probability to extend the indel at each base added
- Maximum Insertion Size: maximal size of insertion

The header is mandatory.

### Profile sequencer

`--input-profile-sequencer-csv` a CSV file, comma separated, with header:
- identifier: an ID choose by `--parameter-profile-sequencer-id-str`
- sequencer: name of the sequencer
- flowcell: kind of flowcell
- version: the kit version.
- strand: `forward` or `reverse`
- cycles total: by strand
- error by cycle: rate of error by cycle, semi-colon separated. Equal to the number of `cycles total`.

The header is mandatory.

### Build a profile sequencer from real data

`build-profile-sequencer` produces a `-input-profile-sequencer-csv`-compatible CSV from
either a pair of FASTQ files or a BAM, so you don't have to already know your
sequencer's error curve:

```sh
HmnRandomRead build-profile-sequencer \
    --parameter-id-str n1 \
    --input-forward-fastq lane1_r1.fastq.gz lane2_r1.fastq.gz \
    --input-reverse-fastq lane1_r2.fastq.gz lane2_r2.fastq.gz \
    --output-profile-sequencer-csv profile_sequencer.csv

HmnRandomRead build-profile-sequencer \
    --parameter-id-str n1 \
    --input-bam lane1.bam lane2.bam \
    --output-profile-sequencer-csv profile_sequencer.csv
```

- `--parameter-id-str`: the identifier written to the `identifier` column
  (matches `--parameter-profile-sequencer-id-str` for `simulate`).
- Input is either `--input-forward-fastq`/`--input-reverse-fastq` together, or
  `--input-bam` alone.
- Each of `--input-forward-fastq`, `--input-reverse-fastq`, and `--input-bam`
  accepts a space-separated list of files (e.g. one per lane); their reads
  are pooled into a single profile. `--input-forward-fastq` and
  `--input-reverse-fastq` must list the same number of files, in matching
  order.
- `sequencer` is guessed from the FASTQ read headers, or from the BAM's `@RG
  PM` tag (falling back to its reads' names).
- `flowcell` and `version` can't be recovered from the data and are always
  written as `NA`.
- `cycles` is the maximum read length seen with actual quality data: the
  longest sequence in the FASTQ files, or the furthest BAM cycle covered by
  retained (non-hard-clipped) quality. Hard-clipped bases carry no quality
  in the BAM record, so they're excluded rather than padded into the
  output — their length is only used to place the retained quality at its
  correct cycle (important on the reverse strand, where clipped/retained
  regions get reordered).
- `error_by_cycle` is derived from the average base quality at each cycle.

## Test

```sh
cargo test
```

## Built with these main libraries

- [rust-htslib](https://github.com/rust-bio/rust-htslib) - Indexed FASTA access

## Authors

- **Guillaume Gricourt**
