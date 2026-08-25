//! Builds a chimeric junction sequence from two breakpoints and draws
//! fusion-supporting read pairs from it — the core of `fusion-in-sample`,
//! which spikes a real sample's FASTQ with synthetic reads at a depth
//! derived from real coverage at the primary breakpoint.

use std::error::Error;
use std::io;

use rust_htslib::bam::{self, Read as BamRead};

use crate::diversity::ProfileDiversity;
use crate::fastq::{build_read_pair, FastqRecord};
use crate::io::FastaIndexedReader;
use crate::profile_sequencer::ProfileSequencer;
use crate::rng::RandomGenerator;
use crate::sequence::Sequence;

const MAX_FRAGMENT_ATTEMPTS: usize = 10;

fn parse_err(msg: impl Into<String>) -> Box<dyn Error> {
    io::Error::new(io::ErrorKind::InvalidData, msg.into()).into()
}

/// A single genomic breakpoint coordinate (1-based): the last reference base
/// kept on the 5' side of the junction for this partner; the base right
/// after it starts its 3' side.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Breakpoint {
    pub chrom: String,
    pub pos: u64,
}

impl Breakpoint {
    /// Parse a `chrom:pos` spec (1-based `pos`).
    pub fn parse(spec: &str) -> Result<Self, Box<dyn Error>> {
        let (chrom, pos) = spec
            .rsplit_once(':')
            .ok_or_else(|| parse_err(format!("malformed breakpoint '{spec}', expected 'chrom:pos'")))?;
        if chrom.is_empty() {
            return Err(parse_err(format!("malformed breakpoint '{spec}', chrom is empty")));
        }
        let pos: u64 = pos
            .parse()
            .map_err(|_| parse_err(format!("malformed breakpoint '{spec}', pos isn't a positive integer")))?;
        if pos == 0 {
            return Err(parse_err(format!("malformed breakpoint '{spec}', pos must be >= 1")));
        }
        Ok(Breakpoint { chrom: chrom.to_string(), pos })
    }
}

enum Side {
    Left,
    Right,
}

/// Fetch up to `flank_len` bases immediately to one side of `bp` (clamped at
/// the contig's boundary, so a breakpoint near the start/end of a contig
/// yields a shorter-than-requested flank rather than an error).
fn fetch_flank(
    faidx: &FastaIndexedReader,
    bp: &Breakpoint,
    flank_len: u64,
    side: Side,
) -> io::Result<String> {
    let contig_len = faidx.seq_len(&bp.chrom)?;
    match side {
        Side::Left => {
            let end0 = bp.pos - 1; // 0-based index of the base at 1-based `pos`
            let begin0 = end0.saturating_sub(flank_len.saturating_sub(1));
            if end0 >= contig_len {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("breakpoint position {} exceeds '{}' length ({contig_len})", bp.pos, bp.chrom),
                ));
            }
            faidx.fetch(&bp.chrom, begin0, end0)
        }
        Side::Right => {
            let begin0 = bp.pos; // 0-based index right after `pos`
            if begin0 >= contig_len {
                return Ok(String::new());
            }
            let end0 = (begin0 + flank_len - 1).min(contig_len - 1);
            faidx.fetch(&bp.chrom, begin0, end0)
        }
    }
}

/// Unlike `simulate`, which only ever draws fragments from non-N
/// [`crate::scaffold::Scaffolds`] runs, a breakpoint here is an exact
/// coordinate the user asked for — there's no other position to silently
/// fall back to. Rather than let an assembly gap near the breakpoint
/// produce synthetic reads full of literal `N` bases (which no real
/// sequencer emits, and which downstream aligners would mishandle), flag it
/// as a clear error instead.
fn contains_n(seq: &str) -> bool {
    seq.bytes().any(|b| b.eq_ignore_ascii_case(&b'N'))
}

/// Build a chimeric junction sequence: `left_bp`'s upstream flank followed by
/// `right_bp`'s downstream flank. Returns the concatenated sequence and the
/// 0-based index of its first base belonging to the right partner (i.e. the
/// junction point), which may be less than `flank_len` if the left flank was
/// clamped near a contig boundary.
pub fn build_junction(
    left_faidx: &FastaIndexedReader,
    left_bp: &Breakpoint,
    right_faidx: &FastaIndexedReader,
    right_bp: &Breakpoint,
    flank_len: u64,
) -> io::Result<(String, usize)> {
    let left = fetch_flank(left_faidx, left_bp, flank_len, Side::Left)?;
    if contains_n(&left) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "the {flank_len}bp flank upstream of breakpoint {}:{} contains N bases \
                 (assembly gap?); pick a breakpoint away from reference gaps",
                left_bp.chrom, left_bp.pos
            ),
        ));
    }
    let right = fetch_flank(right_faidx, right_bp, flank_len, Side::Right)?;
    if contains_n(&right) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "the {flank_len}bp flank downstream of breakpoint {}:{} contains N bases \
                 (assembly gap?); pick a breakpoint away from reference gaps",
                right_bp.chrom, right_bp.pos
            ),
        ));
    }
    let junction_index = left.len();
    Ok((left + &right, junction_index))
}

/// Whether an alignment counts toward [`depth_at`]'s depth — matches
/// `samtools depth`'s default exclude filter (unmapped, secondary, QC-fail,
/// duplicate), since htslib's pileup engine itself does not filter anything
/// on its own: without this, `depth_at` would count every alignment
/// (including secondary/duplicate) and systematically overstate real
/// coverage relative to what a user comparing against `samtools depth`
/// would expect.
fn is_countable(record: &bam::Record) -> bool {
    !record.is_unmapped()
        && !record.is_secondary()
        && !record.is_quality_check_failed()
        && !record.is_duplicate()
}

/// Pileup depth at a single 1-based position of an indexed BAM (0 if the
/// position has no coverage), counting only alignments [`is_countable`].
pub fn depth_at(bam_path: &str, chrom: &str, pos_1based: u64) -> Result<u32, Box<dyn Error>> {
    let mut reader = bam::IndexedReader::from_path(bam_path)?;
    let tid = reader
        .header()
        .tid(chrom.as_bytes())
        .ok_or_else(|| parse_err(format!("chromosome '{chrom}' not found in BAM header")))?;
    let pos0 = pos_1based - 1;
    reader.fetch((tid, pos0 as i64, pos0 as i64 + 1))?;

    for p in reader.pileup() {
        let p = p?;
        if p.tid() == tid && p.pos() as u64 == pos0 {
            let depth = p.alignments().filter(|a| is_countable(&a.record())).count() as u32;
            return Ok(depth);
        }
    }
    Ok(0)
}

/// Try up to [`MAX_FRAGMENT_ATTEMPTS`] times to draw a fragment spanning
/// `junction_index`. The fragment's length is drawn from the same
/// insert-size gaussian as `simulate` (`Normal(mean_insert, std_insert)`),
/// but — unlike an earlier version of this function — its *position* is
/// then drawn uniformly across the exact range that keeps it spanning the
/// junction (`[junction_index - size_insert + 1, junction_index]`, clamped
/// to what fits in `[0, seq_len)`), rather than from a second gaussian
/// jitter reusing `std_insert`: that conflated two different sources of
/// variance (fragment-length variability vs. positional variability around
/// a fixed point) and, for a small `std_insert`, systematically
/// under-placed the junction near either end of the fragment. A uniform
/// draw over the valid spanning range matches what a truly random shearing
/// process would produce, conditioned on the fragment spanning the
/// junction at all.
///
/// A fragment merely *containing* the junction still isn't enough on its
/// own: only the first/last `length_reads` bases of a fragment are ever
/// sequenced (see [`crate::fastq::build_read_pair`]), so a junction that
/// falls in the unsequenced middle of a long fragment would silently
/// produce a read pair with no trace of the fusion in its actual sequence
/// (a "supporting" read that isn't represented at the breakpoint at all).
/// Only fragments where the junction lands inside the head read's window
/// (`[start, start+length_reads)`) or the tail read's window
/// (`[stop-length_reads, stop)`) are accepted. Returns `None` if every
/// attempt failed.
fn pick_fragment_near_junction(
    seq_len: usize,
    junction_index: usize,
    length_reads: usize,
    rng: &mut RandomGenerator,
    mean_insert: f64,
    std_insert: f64,
) -> Option<(usize, usize)> {
    for _ in 0..MAX_FRAGMENT_ATTEMPTS {
        let size_insert = rng.normal(mean_insert, std_insert).round().max(1.0) as i64;

        // The range of `start` values for which [start, start+size_insert)
        // still contains `junction_index`, clamped to what fits in the
        // junction sequence.
        let lo = (junction_index as i64 - size_insert + 1).max(0);
        let hi = (junction_index as i64).min(seq_len as i64 - size_insert);
        if hi < lo {
            continue;
        }
        let start = rng.range(lo, hi);
        let stop = start + size_insert;
        let (start, stop) = (start as usize, stop as usize);

        let take = length_reads.min(stop - start);
        let head_end = start + take;
        let tail_start = stop - take;
        let junction_in_head = start <= junction_index && junction_index < head_end;
        let junction_in_tail = tail_start <= junction_index && junction_index < stop;
        if junction_in_head || junction_in_tail {
            return Some((start, stop));
        }
    }
    None
}

pub struct FusionConfig {
    pub length_reads: usize,
    pub mean_insert_size: f64,
    pub std_insert_size: f64,
    pub profile_diversity: Option<ProfileDiversity>,
    pub id_diversity: Option<String>,
    pub profile_sequencer: Option<ProfileSequencer>,
}

pub struct FusionGenerator {
    config: FusionConfig,
}

impl FusionGenerator {
    pub fn new(config: FusionConfig) -> Self {
        FusionGenerator { config }
    }

    /// Generate `n` fusion-supporting read pairs from `junction`, numbered
    /// `start_number..start_number + n`. Read pairs whose fragment can't be
    /// placed so the junction actually lands in a sequenced read after
    /// [`MAX_FRAGMENT_ATTEMPTS`] attempts are skipped (logged), matching
    /// `simulate`'s behavior for an unfulfillable fragment.
    pub fn generate(
        &self,
        junction: &str,
        junction_index: usize,
        n: u64,
        label: &str,
        rng: &mut RandomGenerator,
        start_number: u64,
    ) -> (Vec<FastqRecord>, Vec<FastqRecord>) {
        let mut forward = Vec::with_capacity(n as usize);
        let mut reverse = Vec::with_capacity(n as usize);
        let mut skipped = 0u64;
        let diversity = self
            .config
            .profile_diversity
            .as_ref()
            .zip(self.config.id_diversity.as_ref())
            .and_then(|(profile, id)| profile.get(id));

        for i in 0..n {
            let number = start_number + i;
            let Some((start, stop)) = pick_fragment_near_junction(
                junction.len(),
                junction_index,
                self.config.length_reads,
                rng,
                self.config.mean_insert_size,
                self.config.std_insert_size,
            ) else {
                // Logged at debug rather than warn: with the split-read
                // requirement above, a low success rate is expected for
                // some parameter combinations (e.g. mean insert size much
                // larger than 2*length_reads) and would otherwise flood the
                // log with one line per skipped read; see the aggregated
                // warning below instead.
                log::debug!(
                    "fusion read {number}: no valid fragment landing the junction in a \
                     sequenced read found after {MAX_FRAGMENT_ATTEMPTS} attempts, skipping"
                );
                skipped += 1;
                continue;
            };

            let sequence = Sequence::new(junction[start..stop].to_string());
            let (r1, r2) = build_read_pair(
                sequence,
                rng,
                diversity,
                self.config.length_reads,
                self.config.profile_sequencer.as_ref(),
                number,
                "fusion".to_string(),
                label.to_string(),
                start as u64,
                stop as u64,
            );
            forward.push(r1);
            reverse.push(r2);
        }

        if skipped > 0 {
            log::warn!("{label}: {skipped} fusion read pair(s) out of {n} requested were skipped");
        }
        (forward, reverse)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_fasta(dir: &std::path::Path, contents: &str) -> String {
        let path = dir.join("ref.fa");
        std::fs::File::create(&path).unwrap().write_all(contents.as_bytes()).unwrap();
        path.to_str().unwrap().to_string()
    }

    #[test]
    fn breakpoint_parse_valid() {
        let bp = Breakpoint::parse("chr1:100").unwrap();
        assert_eq!(bp.chrom, "chr1");
        assert_eq!(bp.pos, 100);
    }

    #[test]
    fn breakpoint_parse_rejects_missing_colon() {
        assert!(Breakpoint::parse("chr1-100").is_err());
    }

    #[test]
    fn breakpoint_parse_rejects_zero_pos() {
        assert!(Breakpoint::parse("chr1:0").is_err());
    }

    #[test]
    fn breakpoint_parse_rejects_non_numeric_pos() {
        assert!(Breakpoint::parse("chr1:abc").is_err());
    }

    #[test]
    fn build_junction_concatenates_flanks_at_expected_index() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), &format!(">chr1\n{}\n", "A".repeat(20)));
        let faidx = FastaIndexedReader::open(&path).unwrap();

        // Breakpoint at 1-based pos 10: left flank is bases [1..=10] (10
        // bases if flank_len=10), right flank is bases [11..=15] (5 bases).
        let left_bp = Breakpoint { chrom: "chr1".into(), pos: 10 };
        let right_bp = Breakpoint { chrom: "chr1".into(), pos: 10 };
        let (junction, idx) = build_junction(&faidx, &left_bp, &faidx, &right_bp, 10).unwrap();
        assert_eq!(idx, 10);
        assert_eq!(junction.len(), 20); // left=10 (bases 1..=10), right=10 (bases 11..=20)
    }

    #[test]
    fn build_junction_rejects_n_bases_in_the_left_flank() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), ">chr1\nACGTNACGTACGTACGTACGT\n");
        let faidx = FastaIndexedReader::open(&path).unwrap();
        // pos 5 (1-based) is the 'N': the left flank [1..=5] contains it.
        let bp = Breakpoint { chrom: "chr1".into(), pos: 5 };
        assert!(build_junction(&faidx, &bp, &faidx, &bp, 5).is_err());
    }

    #[test]
    fn build_junction_rejects_n_bases_in_the_right_flank() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), ">chr1\nACGTACGTNACGTACGTACGT\n");
        let faidx = FastaIndexedReader::open(&path).unwrap();
        // pos 8 (1-based): the right flank [9..=13] contains the 'N' at pos 9.
        let bp = Breakpoint { chrom: "chr1".into(), pos: 8 };
        assert!(build_junction(&faidx, &bp, &faidx, &bp, 5).is_err());
    }

    #[test]
    fn fetch_flank_clamps_at_contig_start() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), ">chr1\nACGTACGTAC\n");
        let faidx = FastaIndexedReader::open(&path).unwrap();
        let bp = Breakpoint { chrom: "chr1".into(), pos: 3 };
        // Requesting 10 bases to the left of pos 3 can only yield 3 (bases 1..=3).
        let left = fetch_flank(&faidx, &bp, 10, Side::Left).unwrap();
        assert_eq!(left, "ACG");
    }

    #[test]
    fn fetch_flank_clamps_at_contig_end() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), ">chr1\nACGTACGTAC\n");
        let faidx = FastaIndexedReader::open(&path).unwrap();
        let bp = Breakpoint { chrom: "chr1".into(), pos: 8 };
        // Requesting 10 bases to the right of pos 8 can only yield 2 (bases 9..=10).
        let right = fetch_flank(&faidx, &bp, 10, Side::Right).unwrap();
        assert_eq!(right, "AC");
    }

    #[test]
    fn fetch_flank_right_is_empty_past_contig_end() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), ">chr1\nACGTACGTAC\n");
        let faidx = FastaIndexedReader::open(&path).unwrap();
        let bp = Breakpoint { chrom: "chr1".into(), pos: 10 };
        let right = fetch_flank(&faidx, &bp, 10, Side::Right).unwrap();
        assert_eq!(right, "");
    }

    #[test]
    fn pick_fragment_always_lands_junction_in_a_sequenced_read() {
        // mean_insert (300) is well over 2*length_reads (100), so plenty of
        // straddling fragments would leave the junction stranded in the
        // unsequenced middle; only the ones near enough to an edge (via the
        // gaussian jitter) should be accepted.
        let mut rng = RandomGenerator::new(42);
        let seq_len = 800;
        let junction_index = 400;
        let length_reads = 50;
        let mut found_any = false;
        for _ in 0..200 {
            if let Some((start, stop)) = pick_fragment_near_junction(
                seq_len,
                junction_index,
                length_reads,
                &mut rng,
                300.0,
                50.0,
            ) {
                found_any = true;
                assert!(stop <= seq_len);
                let take = length_reads.min(stop - start);
                let junction_in_head = start <= junction_index && junction_index < start + take;
                let junction_in_tail = stop - take <= junction_index && junction_index < stop;
                assert!(
                    junction_in_head || junction_in_tail,
                    "junction {junction_index} not in a sequenced window of [{start},{stop})"
                );
            }
        }
        assert!(found_any);
    }

    #[test]
    fn pick_fragment_rejects_junction_stuck_in_the_unsequenced_middle() {
        // A fragment fixed to [0, 300) with a junction at 150 and
        // length_reads=50 never puts the junction in a sequenced read
        // (head covers [0,50), tail covers [250,300)): every attempt must
        // be rejected.
        let mut rng = RandomGenerator::new(1);
        // std tuned so size_insert is deterministically ~300: the only
        // start that keeps a 300bp fragment inside a 300bp sequence while
        // spanning position 150 is 0, which strands the junction in the
        // unsequenced middle every time.
        let result = pick_fragment_near_junction(300, 150, 50, &mut rng, 300.0, 0.001);
        assert!(result.is_none());
    }

    #[test]
    fn generator_produces_requested_pair_count() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), &format!(">chr1\n{}\n", "ACGT".repeat(100)));
        let faidx = FastaIndexedReader::open(&path).unwrap();
        let left_bp = Breakpoint { chrom: "chr1".into(), pos: 200 };
        let right_bp = Breakpoint { chrom: "chr1".into(), pos: 200 };
        let (junction, idx) = build_junction(&faidx, &left_bp, &faidx, &right_bp, 150).unwrap();

        let config = FusionConfig {
            length_reads: 50,
            mean_insert_size: 150.0,
            std_insert_size: 20.0,
            profile_diversity: None,
            id_diversity: None,
            profile_sequencer: None,
        };
        let mut rng = RandomGenerator::new(7);
        let generator = FusionGenerator::new(config);
        let (fwd, rev) = generator.generate(&junction, idx, 10, "chr1:200>chr1:200", &mut rng, 0);
        assert_eq!(fwd.len(), rev.len());
        assert!(fwd.len() <= 10);
        assert!(!fwd.is_empty());
    }

    #[test]
    fn depth_at_counts_overlapping_reads() {
        let dir = tempfile::tempdir().unwrap();
        let sam_path = dir.path().join("reads.sam");
        let sam = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n\
            r1\t0\tchr1\t1\t60\t50M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
            r2\t0\tchr1\t1\t60\t50M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
        std::fs::write(&sam_path, sam).unwrap();

        let bam_path = dir.path().join("reads.bam");
        {
            let sam_reader = bam::Reader::from_path(&sam_path).unwrap();
            let header = bam::Header::from_template(sam_reader.header());
            let mut writer =
                bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
            let mut sam_reader = bam::Reader::from_path(&sam_path).unwrap();
            let mut record = bam::Record::new();
            while let Some(result) = sam_reader.read(&mut record) {
                result.unwrap();
                writer.write(&record).unwrap();
            }
        }
        bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();

        let depth = depth_at(bam_path.to_str().unwrap(), "chr1", 10).unwrap();
        assert_eq!(depth, 2);
    }

    #[test]
    fn depth_at_excludes_secondary_supplementary_and_duplicate_reads() {
        let dir = tempfile::tempdir().unwrap();
        let sam_path = dir.path().join("reads.sam");
        // r1: primary, countable. r2: secondary (flag 256). r3: duplicate
        // (flag 1024). r4: QC-fail (flag 512). Only r1 should be counted.
        let sam = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n\
            r1\t0\tchr1\t1\t60\t50M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
            r2\t256\tchr1\t1\t60\t50M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
            r3\t1024\tchr1\t1\t60\t50M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
            r4\t512\tchr1\t1\t60\t50M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
        std::fs::write(&sam_path, sam).unwrap();

        let bam_path = dir.path().join("reads.bam");
        {
            let sam_reader = bam::Reader::from_path(&sam_path).unwrap();
            let header = bam::Header::from_template(sam_reader.header());
            let mut writer =
                bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
            let mut sam_reader = bam::Reader::from_path(&sam_path).unwrap();
            let mut record = bam::Record::new();
            while let Some(result) = sam_reader.read(&mut record) {
                result.unwrap();
                writer.write(&record).unwrap();
            }
        }
        bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();

        let depth = depth_at(bam_path.to_str().unwrap(), "chr1", 10).unwrap();
        assert_eq!(depth, 1);
    }

    #[test]
    fn depth_at_is_zero_without_coverage() {
        let dir = tempfile::tempdir().unwrap();
        let sam_path = dir.path().join("reads.sam");
        let sam = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n\
            r1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n";
        std::fs::write(&sam_path, sam).unwrap();

        let bam_path = dir.path().join("reads.bam");
        {
            let sam_reader = bam::Reader::from_path(&sam_path).unwrap();
            let header = bam::Header::from_template(sam_reader.header());
            let mut writer =
                bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
            let mut sam_reader = bam::Reader::from_path(&sam_path).unwrap();
            let mut record = bam::Record::new();
            while let Some(result) = sam_reader.read(&mut record) {
                result.unwrap();
                writer.write(&record).unwrap();
            }
        }
        bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();

        let depth = depth_at(bam_path.to_str().unwrap(), "chr1", 500).unwrap();
        assert_eq!(depth, 0);
    }
}
