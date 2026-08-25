//! Computes the fragment insert-size mean and standard deviation from real
//! paired-end BAM data — feeds `simulate`'s `--parameter-mean-insert-int` /
//! `--parameter-std-insert-int`.

use std::error::Error;
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

use rust_htslib::bam::{self, Read as BamRead};

/// Mean and standard deviation of the observed template length (BAM `TLEN`),
/// pooled across every usable read pair found.
pub struct InsertSizeStats {
    pub n: u64,
    pub mean: f64,
    pub std: f64,
}

impl InsertSizeStats {
    /// Build from one or more coordinate- or name-sorted BAMs, pooled into a
    /// single result. Only primary, mapped, properly-paired alignments are
    /// counted, and only once per pair (via `is_first_in_template`) since a
    /// proper pair's `TLEN` has the same magnitude on both mates.
    pub fn from_bam<S: AsRef<str>>(bam_paths: &[S]) -> Result<Self, Box<dyn Error>> {
        // Welford's online algorithm, so pooling many reads never needs to
        // hold every insert size in memory at once.
        let mut n: u64 = 0;
        let mut mean = 0.0f64;
        let mut m2 = 0.0f64;

        for bam_path in bam_paths {
            let mut reader = bam::Reader::from_path(bam_path.as_ref())?;
            let mut record = bam::Record::new();
            while let Some(result) = reader.read(&mut record) {
                result?;
                if !is_usable(&record) {
                    continue;
                }

                n += 1;
                let size = record.insert_size().unsigned_abs() as f64;
                let delta = size - mean;
                mean += delta / n as f64;
                m2 += delta * (size - mean);
            }
        }

        if n == 0 {
            return Err(parse_err(
                "no usable (primary, mapped, properly-paired) read pairs found in any \
                 --input-bam file",
            ));
        }

        let variance = if n > 1 { m2 / (n - 1) as f64 } else { 0.0 };
        Ok(InsertSizeStats {
            n,
            mean,
            std: variance.sqrt(),
        })
    }
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
/// `std_insert_size`, comma-separated with a header row.
pub fn write_csv<P: AsRef<Path>>(entries: &[(String, InsertSizeStats)], path: P) -> io::Result<()> {
    let mut file = File::create(path)?;
    writeln!(file, "file,mean_insert_size,std_insert_size")?;
    for (name, stats) in entries {
        writeln!(file, "{name},{:.2},{:.2}", stats.mean, stats.std)?;
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

        let stats = InsertSizeStats::from_bam(&[path.as_str()]).unwrap();
        assert_eq!(stats.n, 2);
        assert_eq!(stats.mean, 250.0);
        // Sample std of [200, 300]: variance = ((200-250)^2 + (300-250)^2) / (2-1) = 5000.
        assert!((stats.std - 5000f64.sqrt()).abs() < 1e-9);
    }

    #[test]
    fn pools_reads_across_multiple_bam_files() {
        let dir = tempfile::tempdir().unwrap();
        let bam1 = write_file(dir.path(), "a.sam", PAIRED_SAM);
        let bam2 = write_file(dir.path(), "b.sam", PAIRED_SAM);

        let pooled = InsertSizeStats::from_bam(&[bam1.as_str(), bam2.as_str()]).unwrap();
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
        assert!(InsertSizeStats::from_bam(&[path.as_str()]).is_err());
    }

    #[test]
    fn write_csv_writes_one_row_per_entry() {
        let dir = tempfile::tempdir().unwrap();
        let bam_path = write_file(dir.path(), "reads.sam", PAIRED_SAM);
        let stats = InsertSizeStats::from_bam(&[bam_path.as_str()]).unwrap();

        let out_path = dir.path().join("stats.csv");
        write_csv(&[("reads.sam".to_string(), stats)], &out_path).unwrap();

        let contents = std::fs::read_to_string(&out_path).unwrap();
        let mut lines = contents.lines();
        assert_eq!(lines.next(), Some("file,mean_insert_size,std_insert_size"));
        assert_eq!(lines.next(), Some("reads.sam,250.00,70.71"));
        assert_eq!(lines.next(), None);
    }
}
