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
- spike a real sample's FASTQ with a synthetic gene fusion / SV breakpoint
  at a controllable allele fraction

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

HmnRandomRead statistics-insert-size \
    --input-bam <string, required> \
    --output-statistics-csv <string, optional>

HmnRandomRead fusion-in-sample \
    --input-reference-fasta <string, required><int, optional><string, optional> \
    --input-forward-fastq <string, required> \
    --input-reverse-fastq <string, required> \
    --input-bam <string, required> \
    --parameter-breakpoint-primary-roi <string, required> \
    --parameter-breakpoint-secondary-roi <string, required> \
    --parameter-rate-float <float, required> \
    --output-forward-fastq <string, required> \
    --output-reverse-fastq <string, required> \

    --parameter-length-reads-int <int, optional, 150> \
    --parameter-mean-insert-int <int, optional, 500> \
    --parameter-std-insert-int <int, optional, 50> \

    --input-profile-diversity-csv <string, optional> \
    --input-profile-sequencer-csv <string, optional> \
    --parameter-profile-sequencer-id-str <string, optional> \
    --parameter-seed-int <int, optional, 0>

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

### Compute the insert size from real data

`statistics-insert-size` reports the mean and standard deviation of the
fragment insert size (for `--parameter-mean-insert-int` /
`--parameter-std-insert-int`) for each of one or more real BAMs:

```sh
HmnRandomRead statistics-insert-size \
    --input-bam sample.bam

HmnRandomRead statistics-insert-size \
    --input-bam lane1.bam lane2.bam \
    --output-statistics-csv insert_size_stats.csv
```

- `--input-bam` accepts a space-separated list of BAMs (e.g. one per
  sample/lane); each is reported as its own row, not pooled together.
- Only primary, mapped, properly-paired alignments are counted, once per
  pair, from the BAM's `TLEN` field.
- `--output-statistics-csv`, if given, writes a CSV with one row per
  `--input-bam`: `file`, `mean_insert_size`, `std_insert_size`.

### Spike a real sample with a fusion breakpoint

`fusion-in-sample` adds synthetic read pairs supporting a gene fusion /
structural-variant breakpoint to a real sample's paired FASTQ, at a depth
proportional to the real coverage already present at the breakpoint —
useful to build fusion-calling validation datasets with a known, controllable
allele fraction. It reuses `simulate`'s reference, diversity, sequencer
error, insert size, and seed options:

```sh
HmnRandomRead fusion-in-sample \
    --input-reference-fasta genome.fa \
    --input-forward-fastq sample_R1.fastq.gz \
    --input-reverse-fastq sample_R2.fastq.gz \
    --input-bam sample.bam \
    --parameter-breakpoint-primary-roi chr9:130854064 \
    --parameter-breakpoint-secondary-roi chr22:23632600 \
    --parameter-rate-float 0.1 \
    --output-forward-fastq spiked_R1.fastq.gz \
    --output-reverse-fastq spiked_R2.fastq.gz
```

- `--input-forward-fastq`/`--input-reverse-fastq` (plain or gzip-compressed)
  and `--input-bam` (indexed, with a `.bai`/`.csi` sidecar) are the real
  sample's data; all three are required. The output FASTQs contain every
  record from the input FASTQs plus the produced fusion reads appended —
  the input files themselves are never modified.
- `--parameter-breakpoint-primary-roi`/`--parameter-breakpoint-secondary-roi`
  (`chrom:pos`, 1-based) are the two breakpoint partners. `pos` is the last
  reference base kept on the 5' side of the junction for that partner; the
  base right after it starts its 3' side.
- Two chimeric junction sequences are built from the reference: one with the
  primary breakpoint's upstream sequence on the left and the secondary's
  downstream sequence on the right, and the reciprocal (secondary left,
  primary right) — the two derivative junctions of the fusion. Fragments are
  drawn from each junction using the same insert-size gaussian as
  `simulate` (`--parameter-mean-insert-int`/`--parameter-std-insert-int`),
  placed with a gaussian jitter around the junction; a fragment that doesn't
  end up straddling the junction doesn't support the fusion and is
  discarded and redrawn.
- The number of fusion read pairs produced is `round(depth × rate)`, split
  evenly across the two junction orientations, where `depth` is the real
  pileup depth of `--input-bam` at the primary breakpoint position and
  `rate` is `--parameter-rate-float` (0.0-0.5).
- `nb_reads` in an `--input-reference-fasta` spec is ignored by this
  command — the read count comes from depth × rate instead; `id_diversity`
  is still applied (from the reference matching the primary breakpoint) if
  `--input-profile-diversity-csv` is given.
- Produced reads are tagged `fusion` in their FASTQ header, with the
  breakpoint pair (e.g. `chr9:130854064>chr22:23632600`) in place of a
  chromosome name.

## Test

```sh
cargo test
```

## Built with these main libraries

- [rust-htslib](https://github.com/rust-bio/rust-htslib) - Indexed FASTA access

## Authors

- **Guillaume Gricourt**
