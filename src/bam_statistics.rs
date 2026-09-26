//! Statistics of real paired-end BAM data: the fragment insert-size mean and
//! standard deviation — which feed `simulate`'s `--parameter-mean-insert-int`
//! / `--parameter-std-insert-int` — and the maximum depth reached, which
//! bounds how low an allelic fraction `fusion-spike` can carry anywhere in
//! the sample.

use std::collections::BTreeMap;
use std::error::Error;
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

use rust_htslib::bam::{self, ext::BamRecordExtensions, Read as BamRead};

/// Insert size and depth statistics of one or more BAMs.
pub struct BamStatistics {
    /// Read pairs the insert size is computed from.
    pub n: u64,
    /// Mean and standard deviation of the observed template length (BAM
    /// `TLEN`), pooled across every usable read pair found.
    pub mean: f64,
    pub std: f64,
    /// Largest per-position depth, and the first position (`chrom:pos`,
    /// 1-based) it is reached at; `None` without any countable alignment.
    pub max_depth: u64,
    pub max_depth_position: Option<String>,
}

/// Per-position depth over a coordinate-sorted stream of alignments, kept
/// as +1/-1 events at the start/end of each aligned block: events before
/// the current alignment's start can no longer change, so they are folded
/// into the running depth as soon as the stream moves past them, and memory
/// stays proportional to the alignments overlapping the current position.
struct DepthSweep {
    events: BTreeMap<i64, i64>,
    depth: i64,
    max_depth: u64,
    max_depth_position: Option<(i32, i64)>,
    last: Option<(i32, i64)>,
}

impl DepthSweep {
    fn new() -> Self {
        DepthSweep {
            events: BTreeMap::new(),
            depth: 0,
            max_depth: 0,
            max_depth_position: None,
            last: None,
        }
    }

    /// Fold every event before 0-based `until` on contig `tid` into the
    /// running depth, recording each new maximum.
    fn flush(&mut self, tid: i32, until: i64) {
        while let Some((&pos, _)) = self.events.first_key_value() {
            if pos >= until {
                break;
            }
            let (pos, delta) = self.events.pop_first().unwrap();
            self.depth += delta;
            if self.depth as u64 > self.max_depth {
                self.max_depth = self.depth as u64;
                self.max_depth_position = Some((tid, pos));
            }
        }
    }

    fn add(&mut self, record: &bam::Record) -> Result<(), Box<dyn Error>> {
        let here = (record.tid(), record.pos());
        match self.last {
            Some(last) if here < last => {
                return Err(parse_err(
                    "BAM is not coordinate-sorted, which the maximum depth needs: sort it \
                     with `samtools sort`",
                ));
            }
            Some((tid, _)) if tid != here.0 => {
                self.flush(tid, i64::MAX);
                self.depth = 0;
            }
            _ => {}
        }
        self.last = Some(here);
        self.flush(here.0, here.1);
        for [start, end] in record.aligned_blocks() {
            *self.events.entry(start).or_default() += 1;
            *self.events.entry(end).or_default() -= 1;
        }
        Ok(())
    }

    fn finish(mut self) -> (u64, Option<(i32, i64)>) {
        if let Some((tid, _)) = self.last {
            self.flush(tid, i64::MAX);
        }
        (self.max_depth, self.max_depth_position)
    }
}

impl BamStatistics {
    /// Build from one or more coordinate-sorted BAMs. The insert size is
    /// pooled into a single result: only primary, mapped, properly-paired
    /// alignments are counted, and only once per pair (via
    /// `is_first_in_template`) since a proper pair's `TLEN` has the same
    /// magnitude on both mates. The depth is `samtools depth`'s default one
    /// — unmapped, secondary, QC-fail and duplicate alignments left out,
    /// deletions and skipped regions not counted as covered — and the
    /// maximum is taken over the files, never summed across them.
    pub fn from_bam<S: AsRef<str>>(bam_paths: &[S]) -> Result<Self, Box<dyn Error>> {
        // Welford's online algorithm, so pooling many reads never needs to
        // hold every insert size in memory at once.
        let mut n: u64 = 0;
        let mut mean = 0.0f64;
        let mut m2 = 0.0f64;
        let mut max_depth = 0u64;
        let mut max_depth_position = None;

        for bam_path in bam_paths {
            let mut reader = bam::Reader::from_path(bam_path.as_ref())?;
            let header = reader.header().clone();
            let mut sweep = DepthSweep::new();
            let mut record = bam::Record::new();
            while let Some(result) = reader.read(&mut record) {
                result?;
                if is_depth_counted(&record) {
                    sweep.add(&record)?;
                }
                if !is_usable(&record) {
                    continue;
                }

                n += 1;
                let size = record.insert_size().unsigned_abs() as f64;
                let delta = size - mean;
                mean += delta / n as f64;
                m2 += delta * (size - mean);
            }
            let (depth, position) = sweep.finish();
            if depth > max_depth {
                max_depth = depth;
                max_depth_position = position.map(|(tid, pos)| {
                    format!(
                        "{}:{}",
                        String::from_utf8_lossy(header.tid2name(tid as u32)),
                        pos + 1
                    )
                });
            }
        }

        if n == 0 {
            return Err(parse_err(
                "no usable (primary, mapped, properly-paired) read pairs found in any \
                 --input-bam file",
            ));
        }

        let variance = if n > 1 { m2 / (n - 1) as f64 } else { 0.0 };
        Ok(BamStatistics {
            n,
            mean,
            std: variance.sqrt(),
            max_depth,
            max_depth_position,
        })
    }
}

/// An alignment `samtools depth` counts by default.
fn is_depth_counted(record: &bam::Record) -> bool {
    !record.is_unmapped()
        && !record.is_secondary()
        && !record.is_quality_check_failed()
        && !record.is_duplicate()
}

/// A read counted at most once per pair (first-in-template only), and only
/// when its `TLEN` is actually meaningful.
fn is_usable(record: &bam::Record) -> bool {
    !record.is_secondary()
        && !record.is_supplementary()
        && !record.is_unmapped()
        && !record.is_mate_unmapped()
        && record.is_paired()
        && record.is_proper_pair()
        && record.is_first_in_template()
        && record.insert_size() != 0
}

fn parse_err(msg: impl Into<String>) -> Box<dyn Error> {
    io::Error::new(io::ErrorKind::InvalidData, msg.into()).into()
}

/// Writes one row per `(file, stats)` entry: `file`, `mean_insert_size`,
/// `std_insert_size`, `max_depth`, `max_depth_position`, comma-separated
/// with a header row.
pub fn write_csv<P: AsRef<Path>>(entries: &[(String, BamStatistics)], path: P) -> io::Result<()> {
    let mut file = File::create(path)?;
    writeln!(
        file,
        "file,mean_insert_size,std_insert_size,max_depth,max_depth_position"
    )?;
    for (name, stats) in entries {
        writeln!(
            file,
            "{name},{:.2},{:.2},{},{}",
            stats.mean,
            stats.std,
            stats.max_depth,
            stats.max_depth_position.as_deref().unwrap_or("")
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::path::Path;

    fn write_file(dir: &Path, name: &str, contents: &str) -> String {
        let path = dir.join(name);
        File::create(&path).unwrap().write_all(contents.as_bytes()).unwrap();
        path.to_str().unwrap().to_string()
    }

    /// Two properly-paired reads (SAM, which htslib reads the same as BAM):
    /// one pair with TLEN 200, one pair with TLEN 300.
    const PAIRED_SAM: &str = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10000\n\
        r1\t99\tchr1\t101\t60\t50M\t=\t251\t200\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
        r1\t147\tchr1\t251\t60\t50M\t=\t101\t-200\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
        r2\t99\tchr1\t5001\t60\t50M\t=\t5251\t300\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
        r2\t147\tchr1\t5251\t60\t50M\t=\t5001\t-300\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";

    #[test]
    fn computes_mean_and_std_across_pairs() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_file(dir.path(), "reads.sam", PAIRED_SAM);

        let stats = BamStatistics::from_bam(&[path.as_str()]).unwrap();
        assert_eq!(stats.n, 2);
        assert_eq!(stats.mean, 250.0);
        // Sample std of [200, 300]: variance = ((200-250)^2 + (300-250)^2) / (2-1) = 5000.
        assert!((stats.std - 5000f64.sqrt()).abs() < 1e-9);
    }

    #[test]
    fn max_depth_is_the_deepest_position() {
        // 1-based: r1 covers [1,50]; r2 [5,19] then, past a 10-base
        // deletion, [30,64]; r3 [25,74]; d1 is a duplicate. Depth 2 over
        // [5,19], 1 over [20,24] (the deletion isn't covered), 2 over
        // [25,29], 3 over [30,50]: maximum 3, first reached at 30.
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
        let qual = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
        let sam = format!(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n\
             r1\t0\tchr1\t1\t60\t50M\t*\t0\t0\t{seq}\t{qual}\n\
             r2\t0\tchr1\t5\t60\t15M10D35M\t*\t0\t0\t{seq}\t{qual}\n\
             d1\t1024\tchr1\t10\t60\t50M\t*\t0\t0\t{seq}\t{qual}\n\
             r3\t0\tchr1\t25\t60\t50M\t*\t0\t0\t{seq}\t{qual}\n\
             c1\t0\tchr2\t1\t60\t50M\t*\t0\t0\t{seq}\t{qual}\n\
             c2\t0\tchr2\t1\t60\t50M\t*\t0\t0\t{seq}\t{qual}\n"
        );
        let dir = tempfile::tempdir().unwrap();
        let path = write_file(dir.path(), "reads.sam", &sam);
        // No proper pair: the insert size can't be computed...
        assert!(BamStatistics::from_bam(&[path.as_str()]).is_err());

        // ...so check the sweep on its own.
        let mut reader = bam::Reader::from_path(&path).unwrap();
        let mut sweep = DepthSweep::new();
        let mut record = bam::Record::new();
        while let Some(result) = reader.read(&mut record) {
            result.unwrap();
            if is_depth_counted(&record) {
                sweep.add(&record).unwrap();
            }
        }
        // chr1: 3 over [30,50] (the duplicate d1 left out); chr2: 2. The
        // depth on chr1 doesn't leak into chr2.
        assert_eq!(sweep.finish(), (3, Some((0, 29))));
    }

    #[test]
    fn max_depth_requires_a_coordinate_sorted_bam() {
        let seq = "ACGTACGTAC";
        let qual = "IIIIIIIIII";
        let sam = format!(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n\
             a\t99\tchr1\t500\t60\t10M\t=\t600\t110\t{seq}\t{qual}\n\
             b\t99\tchr1\t100\t60\t10M\t=\t200\t110\t{seq}\t{qual}\n"
        );
        let dir = tempfile::tempdir().unwrap();
        let path = write_file(dir.path(), "reads.sam", &sam);
        assert!(BamStatistics::from_bam(&[path.as_str()]).is_err());
    }

    #[test]
    fn pools_reads_across_multiple_bam_files() {
        let dir = tempfile::tempdir().unwrap();
        let bam1 = write_file(dir.path(), "a.sam", PAIRED_SAM);
        let bam2 = write_file(dir.path(), "b.sam", PAIRED_SAM);

        let pooled = BamStatistics::from_bam(&[bam1.as_str(), bam2.as_str()]).unwrap();
        assert_eq!(pooled.n, 4);
        assert_eq!(pooled.mean, 250.0);
    }

    #[test]
    fn empty_bam_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_file(
            dir.path(),
            "empty.sam",
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10000\n",
        );
        assert!(BamStatistics::from_bam(&[path.as_str()]).is_err());
    }

    #[test]
    fn write_csv_writes_one_row_per_entry() {
        let dir = tempfile::tempdir().unwrap();
        let bam_path = write_file(dir.path(), "reads.sam", PAIRED_SAM);
        let stats = BamStatistics::from_bam(&[bam_path.as_str()]).unwrap();

        let out_path = dir.path().join("stats.csv");
        write_csv(&[("reads.sam".to_string(), stats)], &out_path).unwrap();

        let contents = std::fs::read_to_string(&out_path).unwrap();
        let mut lines = contents.lines();
        assert_eq!(
            lines.next(),
            Some("file,mean_insert_size,std_insert_size,max_depth,max_depth_position")
        );
        assert_eq!(lines.next(), Some("reads.sam,250.00,70.71,1,chr1:101"));
        assert_eq!(lines.next(), None);
    }
}
