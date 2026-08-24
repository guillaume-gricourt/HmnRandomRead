//! Builds a chimeric junction sequence from two breakpoints and draws
//! fusion-supporting read pairs from it — the core of `fusion-in-sample`,
//! which spikes a real sample's FASTQ with synthetic reads at a depth
//! derived from real coverage at the primary breakpoint.

use std::error::Error;
use std::io;

use rust_htslib::bam::{self, Read as BamRead};

use crate::diversity::ProfileDiversity;
use crate::fastq::FastqRecord;
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
    let right = fetch_flank(right_faidx, right_bp, flank_len, Side::Right)?;
    let junction_index = left.len();
    Ok((left + &right, junction_index))
}

/// Pileup depth at a single 1-based position of an indexed BAM (0 if the
/// position has no coverage).
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
            return Ok(p.depth());
        }
    }
    Ok(0)
}

/// Try up to [`MAX_FRAGMENT_ATTEMPTS`] times to draw a fragment (insert-size
/// gaussian, as in `simulate`) placed by a gaussian jitter around
/// `junction_index`, keeping only fragments that actually straddle the
/// junction (`start <= junction_index < stop`) — the read pairs the user
/// wants discarded are exactly the ones that don't. Returns `None` if every
/// attempt failed.
fn pick_fragment_near_junction(
    seq_len: usize,
    junction_index: usize,
    rng: &mut RandomGenerator,
    mean_insert: f64,
    std_insert: f64,
) -> Option<(usize, usize)> {
    for _ in 0..MAX_FRAGMENT_ATTEMPTS {
        let size_insert = rng.normal(mean_insert, std_insert).round().max(1.0) as i64;
        let offset = rng.normal(0.0, std_insert).round() as i64;
        let center = junction_index as i64 + offset;
        let start = center - size_insert / 2;
        let stop = start + size_insert;
        if start < 0 || stop <= start || stop as usize > seq_len {
            continue;
        }
        let (start, stop) = (start as usize, stop as usize);
        if start <= junction_index && junction_index < stop {
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
    /// placed to straddle the junction after [`MAX_FRAGMENT_ATTEMPTS`]
    /// attempts are skipped (logged), matching `simulate`'s behavior for an
    /// unfulfillable fragment.
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
                rng,
                self.config.mean_insert_size,
                self.config.std_insert_size,
            ) else {
                log::warn!(
                    "fusion read {number}: no valid fragment straddling the junction found after \
                     {MAX_FRAGMENT_ATTEMPTS} attempts, skipping"
                );
                continue;
            };

            let mut sequence = Sequence::new(junction[start..stop].to_string());
            if let Some(diversity) = diversity {
                sequence.make_mutation(rng, diversity);
            }

            let head_seq = sequence.sub_read(self.config.length_reads, true, false);
            let tail_seq = sequence.sub_read(self.config.length_reads, false, true);
            let head_len = head_seq.len() as u64;
            let tail_len = tail_seq.len() as u64;

            let mut head = FastqRecord::new(
                head_seq,
                33,
                number,
                true,
                "fusion".to_string(),
                label.to_string(),
                start as u64,
                start as u64 + head_len,
            );
            let mut tail = FastqRecord::new(
                tail_seq,
                33,
                number,
                false,
                "fusion".to_string(),
                label.to_string(),
                (stop as u64).saturating_sub(tail_len),
                stop as u64,
            );

            head.init_qual(rng);
            tail.init_qual(rng);
            if let Some(profile_sequencer) = &self.config.profile_sequencer {
                head.make_errors(rng, profile_sequencer);
                tail.make_errors(rng, profile_sequencer);
            }

            let (r1, r2) = if rng.unit() >= 0.5 { (head, tail) } else { (tail, head) };
            forward.push(r1);
            reverse.push(r2);
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
    fn pick_fragment_always_straddles_junction() {
        let mut rng = RandomGenerator::new(42);
        let seq_len = 400;
        let junction_index = 200;
        for _ in 0..200 {
            if let Some((start, stop)) =
                pick_fragment_near_junction(seq_len, junction_index, &mut rng, 100.0, 10.0)
            {
                assert!(start <= junction_index && junction_index < stop);
                assert!(stop <= seq_len);
            }
        }
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
