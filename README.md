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
- build a pool of synthetic gene fusion / SV breakpoint reads, then dilute it
  into a real sample's FASTQ at a known allelic fraction

## Install

```sh
conda install -c bioconda hmnrandomread
```

### Download a binary

Each [GitHub release](https://github.com/guillaume-gricourt/HmnRandomRead/releases)
ships a prebuilt binary for Linux x86_64 (glibc ≥ 2.35, OpenSSL 3), macOS
x86_64 and macOS arm64 (macOS ≥ 11):

```sh
VERSION=0.12.0
TARGET=x86_64-unknown-linux-gnu  # or x86_64-apple-darwin, aarch64-apple-darwin
curl -L -O "https://github.com/guillaume-gricourt/HmnRandomRead/releases/download/$VERSION/HmnRandomRead-$VERSION-$TARGET.tar.gz"
tar -xzf "HmnRandomRead-$VERSION-$TARGET.tar.gz"
"HmnRandomRead-$VERSION-$TARGET/HmnRandomRead" version
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

HmnRandomRead statistics-bam \
    --input-bam <string, required> \
    --output-statistics-csv <string, optional>

HmnRandomRead fusion-simulate \
    --input-reference-fasta <string, required><int, required><string, optional> \
    --parameter-breakpoint-primary-roi <string, required> \
    --parameter-breakpoint-secondary-roi <string, required> \
    --parameter-forward-strand-rate-float <float, optional, 1.0> \
    --parameter-reciprocal-rate-float <float, optional, 0.5> \
    --output-forward-fastq <string, required> \
    --output-reverse-fastq <string, required> \
    --output-truth-tsv <string, required> \

    --parameter-length-reads-int <int, optional, 150> \
    --parameter-mean-insert-int <int, optional, 500> \
    --parameter-std-insert-int <int, optional, 50> \
    --parameter-minimum-anchor-int <int, optional, 20> \

    --input-profile-diversity-csv <string, optional> \
    --input-profile-sequencer-csv <string, optional> \
    --parameter-profile-sequencer-id-str <string, optional> \
    --parameter-seed-int <int, optional, 0>

HmnRandomRead fusion-spike \
    --input-forward-fastq <string, required> \
    --input-reverse-fastq <string, required> \
    --input-bam <string, required> \
    --input-reference-fasta <string, required> \
    --input-fusion-forward-fastq <string, required> \
    --input-fusion-reverse-fastq <string, required> \
    --input-truth-tsv <string, required> \
    --parameter-allelic-fraction-float <float, required> \
    --output-forward-fastq <string, required> \
    --output-reverse-fastq <string, required> \
    --output-truth-tsv <string, required> \
    --parameter-seed-int <int, optional, 0>

HmnRandomRead version
```

### Reference

Use one or more FASTA file used as reference sequence (`--input-reference-fasta`, may be
repeated: `path[,nb_reads[,id_diversity]]`).
Indicate also the number of sequence to generate for each reference.
`nb_reads` may be left empty (`path,`, `path,,id_diversity`) and then defaults
to 0 — handy for a reference brought in only to resolve a contig, such as a
fusion partner's chromosome, which still needs its `id_diversity`.

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

### Compute statistics from real data

`statistics-bam` reports, for each of one or more real BAMs, the mean and
standard deviation of the fragment insert size (for
`--parameter-mean-insert-int` / `--parameter-std-insert-int`) and the
maximum depth reached, with the position it is first reached at:

```sh
HmnRandomRead statistics-bam \
    --input-bam sample.bam

HmnRandomRead statistics-bam \
    --input-bam lane1.bam lane2.bam \
    --output-statistics-csv bam_stats.csv
```

- `--input-bam` accepts a space-separated list of coordinate-sorted BAMs
  (e.g. one per sample/lane); each is reported as its own row, not pooled
  together. An unsorted BAM is refused.
- The insert size counts only primary, mapped, properly-paired alignments,
  once per pair, from the BAM's `TLEN` field.
- The depth is the one `samtools depth` reports by default: unmapped,
  secondary, QC-fail and duplicate alignments are left out, and deletions
  and skipped regions don't count as covered. The maximum depth bounds the
  lowest allelic fraction `fusion-spike` can carry anywhere in the sample
  (it needs `f × depth ≥ 1` fragment at the breakpoint).
- `--output-statistics-csv`, if given, writes a CSV with one row per
  `--input-bam`: `file`, `mean_insert_size`, `std_insert_size`,
  `max_depth`, `max_depth_position` (`chrom:pos`, 1-based).

### Build a fusion read pool

`fusion-simulate` writes a pool of synthetic read pairs supporting a gene
fusion / structural-variant breakpoint, and nothing else — no sample reads are
merged in, so one pool can be diluted into any sample at any fraction
afterwards. It reuses `simulate`'s reference, diversity, sequencer error,
insert size and seed options:

```sh
HmnRandomRead fusion-simulate \
    --input-reference-fasta genome.fa,50000,human \
    --parameter-breakpoint-primary-roi chr9:130854064 \
    --parameter-breakpoint-secondary-roi chr22:23632600 \
    --parameter-reciprocal-rate-float 0.0 \
    --output-forward-fastq pool_R1.fastq.gz \
    --output-reverse-fastq pool_R2.fastq.gz \
    --output-truth-tsv pool.tsv
```

- `--parameter-breakpoint-primary-roi`/`--parameter-breakpoint-secondary-roi`
  (`chrom:pos`, 1-based) are the two breakpoint partners: the primary one is
  the 5' partner of the primary→secondary junction, the secondary one its 3'
  partner. `pos` is the last reference base before the break, which always
  falls between `pos` and `pos + 1`. If either breakpoint's flank contains
  an `N` (assembly gap), the command fails with a clear error rather than
  emitting reads full of `N`.
- `--parameter-forward-strand-rate-float` (0.0-1.0, default 1.0) is the
  probability for each partner of a read pair to be taken on the `+` strand
  rather than the `-` one, drawn independently for the two partners of
  every pair. The strand decides the side of the break a partner
  contributes and its orientation:

  | strand | as the 5' partner                        | as the 3' partner                        |
  |--------|------------------------------------------|------------------------------------------|
  | `+`    | bases up to `pos`                        | bases from `pos + 1` on                  |
  | `-`    | reverse complement of bases from `pos + 1` on | reverse complement of bases up to `pos` |

  Two partners on the same strand give a translocation (e.g. PML-RARA,
  BCR-ABL1); partners on opposite strands an inversion-type junction (e.g.
  EML4-ALK). With a rate `r`, the pool mixes the four geometries `+/+`,
  `+/-`, `-/+` and `-/-` at `r²`, `r(1-r)`, `(1-r)r` and `(1-r)²`; the
  default 1.0 keeps every pair on `+/+`. Each pair's geometry is in its
  truth row's junction label (`chrom:pos(strand)>chrom:pos(strand)`) and
  the rate in `#forward_strand_rate`. A geometry mix splits the support
  between junctions a caller reports separately, like
  `--parameter-reciprocal-rate-float` does: keep 1.0 (or 0.0) when
  quantifying a single junction.
- The number of read pairs to produce is the `nb_reads` field of the
  `--input-reference-fasta` spec holding the primary breakpoint — the same
  spec that supplies `id_diversity`. Generate the pool with a wide margin: it
  costs little, and `fusion-spike` samples from it.
- Two chimeric junction sequences are built: the primary partner's 5' part
  followed by the secondary partner's 3' part, and the reciprocal (the
  secondary partner's 5' part followed by the primary partner's 3' part) —
  the two derivative junctions of the fusion.
  `--parameter-reciprocal-rate-float` (0.0-1.0, default 0.5) is the fraction
  of produced pairs assigned to the reciprocal orientation: 0.5 splits
  evenly (a balanced reciprocal translocation), 0.0 produces only the
  primary→secondary junction, 1.0 only the reciprocal one. **Use 0.0 (or
  1.0) when quantifying a single junction**: a caller reports the two
  derivatives as separate events, so splitting the pool halves the support
  behind each of them.
- Every fragment crossing the junction is a fusion molecule, and all of
  them are kept, as random shearing produces them: the fragment length is
  drawn from the library's insert-size gaussian
  (`--parameter-mean-insert-int`/`--parameter-std-insert-int`) weighted by
  that length — a molecule crosses a given point with probability
  proportional to its length, so the fragments over a breakpoint are longer
  than the library average, exactly like the reference fragments
  `fusion-spike` counts there — and its position is uniform over every start
  leaving at least one base on each side of the junction. Nothing is
  filtered on where the junction ends up; each pair is instead recorded
  with its support type:
  - `split_read`: a read crosses the junction with at least
    `--parameter-minimum-anchor-int` (default 20) bases on each side — a
    chimeric/supplementary alignment (short-read aligners need a minimum
    seed length — e.g. BWA-MEM's default is 19bp — and SV/fusion callers add
    their own minimum overhang, typically 20-25bp);
  - `short_anchor`: a read crosses the junction with fewer bases on one
    side, and will usually align whole to one partner, soft-clipped;
  - `discordant_pair`: the junction falls in the unsequenced insert, each
    mate aligns to a different partner.

  Their counts are in the truth TSV metadata (`#split_read_pairs`,
  `#short_anchor_pairs`, `#discordant_pairs`); their share follows from the
  insert size and read length, as in the real library.
- Both mates of a pair share one whitespace-free, self-identifying QNAME:
  `hmnrr_fusion_<left_chrom>-<left_pos>_<right_chrom>-<right_pos>_<number>`.
  Everything lives in the name rather than in the FASTQ comment, which
  aligners drop unless told otherwise (`bwa mem -C`), so the reads stay
  findable in an aligned BAM with a single name filter:

  ```sh
  samtools view spiked.bam | grep hmnrr_fusion
  ```

- `--output-truth-tsv` records one row per produced pair — QNAME, junction
  (`chrom:pos(strand)>chrom:pos(strand)`), fragment bounds, the offset of
  the junction within the fragment (the clip point an aligner should find)
  and the support type — under a `#key<TAB>value` metadata block. A
  fragment that can't fit in the junction sequence (beyond `mean + 4 std`,
  or near a contig end) is redrawn and, after 100 attempts, skipped: the
  `#produced_pairs` count, not the requested one, is what was written.

### Dilute a pool into a real sample

`fusion-spike` injects part of a pool into a real sample's paired FASTQ at a
target allelic fraction, scattering the injected pairs at random positions
through the output rather than appending them, and removes the reference
fragments the fusion takes the place of:

```sh
HmnRandomRead fusion-spike \
    --input-forward-fastq sample_R1.fastq.gz \
    --input-reverse-fastq sample_R2.fastq.gz \
    --input-bam sample.bam \
    --input-reference-fasta genome.fa \
    --input-fusion-forward-fastq pool_R1.fastq.gz \
    --input-fusion-reverse-fastq pool_R2.fastq.gz \
    --input-truth-tsv pool.tsv \
    --parameter-allelic-fraction-float 0.001 \
    --output-forward-fastq spiked_R1.fastq.gz \
    --output-reverse-fastq spiked_R2.fastq.gz \
    --output-truth-tsv spiked.tsv
```

- `F` is the number of distinct reference **fragments** crossing the
  primary breakpoint in `--input-bam` (primary, mapped, properly-paired,
  non-duplicate, non-supplementary, QC-pass alignments, counted once per
  template molecule), whether one of their reads covers the break or it
  falls in their insert — the same class of molecule as the pool, which
  holds discordant pairs too. Counting fragments rather than pileup reads is
  what makes `f` (`--parameter-allelic-fraction-float`) the fraction of
  molecules carrying the fusion: when the mean insert size is under twice
  the read length the two mates of a fragment overlap, so one molecule
  contributes two alignments to the pileup. Both counts are reported, so
  their ratio stays visible.
- `--input-bam` is expected **duplicate-marked** (Picard MarkDuplicates,
  samtools markdup…), the way the pipeline under evaluation will process
  the spiked output. `F` counts molecules, so pairs flagged as duplicates
  are left out of it: in a BAM that was never marked, PCR copies count as
  molecules, too many pairs are injected, and once the output is
  deduplicated the realized fraction ends up `F_total / F_unique` times too
  high. A warning is logged when not a single pair crossing a breakpoint of
  at least 100 molecules is flagged as a duplicate. (Leave the BAM unmarked
  only if the pipeline under evaluation doesn't deduplicate either, e.g.
  amplicons.)
- The fusion **replaces** reference molecules instead of adding to them: a
  fusion allele is a reference allele broken and rejoined. For each pair of
  the more represented derivative junction, one reference fragment crossing
  the primary breakpoint and one crossing the secondary breakpoint are
  removed from the output (both mates, matched by read name between the BAM
  and the FASTQ — `--input-bam` must be aligned from the input FASTQ),
  **together with every pair flagged as its duplicate**: a copy left behind
  would be elected the molecule's representative when the output is
  realigned and deduplicated, and the molecule would still be counted.
  Duplicates are matched to their molecule the way markers group them:
  same library (read group `LB`), same unclipped 5' ends of both mates
  (the mate's from its `MC` tag when present) and same strands. With
  `--parameter-reciprocal-rate-float 0` in the pool, that is `round(f * F)`
  pairs injected and as many fragments removed at each breakpoint: the
  molecule count at the primary breakpoint stays `F`. With a balanced pool
  (0.5), both derivatives come from the same two broken homologs, so only
  half as many fragments are removed as pairs injected, and the pair count
  is `round(f * F / (1 - f / 2))`. If the secondary breakpoint has fewer
  reference fragments than needed (e.g. outside a capture panel), as many as
  there are are removed, with a warning.
- The injected reads carry the sample's own alleles: a pool is built from
  the bare reference, so that one pool serves any sample, and the variants
  are added at spike time. The sample's SNVs and indels are called from
  `--input-bam` over both flanks of each breakpoint, one observation per
  fragment — where the two mates of a pair overlap, they read the same
  molecule and count once (agreeing, one vote; disagreeing, the higher base
  quality wins, a tie drops the fragment at that position) — with mapping
  quality ≥ 20, base quality ≥ 13, depth ≥ 10 fragments, and an alternate
  allele fraction ≥ 0.8 for homozygous, ≥ 0.2 for heterozygous. Homozygous
  variants are always carried; heterozygous ones are phased onto a single
  haplotype from the read pairs linking them (a random draw where none
  does), since a fusion happens on one homolog. Each injected pair's
  fragment is then rebuilt from `--input-reference-fasta` (the reference
  the pool was generated from, repeatable), with and without the variants:
  the bases the variants change are rewritten, the pool's own sequencing
  errors and qualities are kept, and the fragment length is unchanged, an
  indel only shifting the bases away from the junction. The variants
  carried are listed in the truth TSV (`#haplotype_variants_primary`,
  `#haplotype_variants_secondary`), with the number of pairs they changed
  (`#haplotype_rewritten_pairs`). A pool built with a diversity profile
  carrying indels loses the sample's variants past such an indel in a read.
- The breakpoint to measure that depth at comes from the pool's
  `--input-truth-tsv`, so a pool can never be diluted against the wrong
  locus.
- Pairs are drawn from the pool uniformly without replacement, never as a
  contiguous slice: a pool split between the two derivative junctions holds
  all of one orientation before any of the other. Use a different
  `--parameter-seed-int` per replicate of a dilution series.
- `--output-truth-tsv` lists the pairs actually injected, under a metadata
  block recording `#reference_fragments` (and `_secondary`),
  `#reference_pileup_reads`, `#target_allelic_fraction`,
  `#realized_allelic_fraction` (and `_secondary`), `#injected_pairs` (split
  into `#injected_primary_secondary_pairs`/`#injected_reciprocal_pairs`),
  `#removed_reference_fragments` (and `_primary`/`_secondary`, counted in
  molecules), `#removed_reference_duplicates` (the duplicate pairs removed
  with them), `#reference_duplicate_pairs`,
  `#sample_pairs` and `#output_pairs`. This is the ground truth of the dataset: quantify
  against it rather than against the reads recovered from the alignment, or
  the measurement ends up conditioned on the caller being evaluated.
- The input FASTQs are never modified, and the command fails before writing
  anything if the two mate files disagree on their record count, if they
  don't hold every fragment picked for removal, or if the pool is too small
  for the requested fraction.

## Test

```sh
cargo test
```

## Built with these main libraries

- [rust-htslib](https://github.com/rust-bio/rust-htslib) - Indexed FASTA and BAM access

## Authors

- **Guillaume Gricourt**
