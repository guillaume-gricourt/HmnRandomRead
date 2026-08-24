//! Computes the fragment insert-size mean and standard deviation from real
//! paired-end BAM data, optionally restricted to regions in a BED file —
//! feeds `simulate`'s `--parameter-mean-insert-int` / `--parameter-std-insert-int`.

use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{self, Read as BamRead};

/// Per-chromosome list of `(start, end)` BED intervals (0-based, half-open).
type Regions = HashMap<String, Vec<(i64, i64)>>;

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
    ///
    /// `bed_path`, if given, restricts the tally to read pairs whose
    /// first-in-template alignment overlaps at least one region listed in
    /// the (0-based, half-open) BED file.
    pub fn from_bam<S: AsRef<str>>(
        bam_paths: &[S],
        bed_path: Option<&str>,
    ) -> Result<Self, Box<dyn Error>> {
        let regions = bed_path.map(parse_bed).transpose()?;

        // Welford's online algorithm, so pooling many reads never needs to
        // hold every insert size in memory at once.
        let mut n: u64 = 0;
        let mut mean = 0.0f64;
        let mut m2 = 0.0f64;

        for bam_path in bam_paths {
            let mut reader = bam::Reader::from_path(bam_path.as_ref())?;
            let target_names: Vec<String> = reader
                .header()
                .target_names()
                .iter()
                .map(|name| String::from_utf8_lossy(name).into_owned())
                .collect();

            let mut record = bam::Record::new();
            while let Some(result) = reader.read(&mut record) {
                result?;
                if !is_usable(&record) {
                    continue;
                }
                if let Some(regions) = &regions {
                    let chrom = &target_names[record.tid() as usize];
                    if !overlaps(regions, chrom, record.pos(), record.reference_end()) {
                        continue;
                    }
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
                 --input-bam file, within --input-bed if given",
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

fn overlaps(regions: &Regions, chrom: &str, start: i64, end: i64) -> bool {
    regions
        .get(chrom)
        .is_some_and(|intervals| intervals.iter().any(|&(s, e)| start < e && end > s))
}

/// Minimal BED parser: `chrom`, `start`, `end` (0-based, half-open) columns,
/// tab- or space-separated; any further columns are ignored. Comment/header
/// (`#`, `track`, `browser`) and blank lines are skipped, matching the BED
/// format's own conventions.
fn parse_bed(path: &str) -> Result<Regions, Box<dyn Error>> {
    let file = File::open(path)?;
    let mut regions: Regions = HashMap::new();
    for line in BufReader::new(file).lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("track")
            || line.starts_with("browser")
        {
            continue;
        }
        let mut fields = line.split_whitespace();
        let chrom = fields
            .next()
            .ok_or_else(|| parse_err("BED line is missing its chrom column"))?;
        let start: i64 = fields
            .next()
            .and_then(|s| s.parse().ok())
            .ok_or_else(|| parse_err("BED line is missing/has an invalid start column"))?;
        let end: i64 = fields
            .next()
            .and_then(|s| s.parse().ok())
            .ok_or_else(|| parse_err("BED line is missing/has an invalid end column"))?;
        regions.entry(chrom.to_string()).or_default().push((start, end));
    }
    if regions.is_empty() {
        return Err(parse_err("no regions found in BED file"));
    }
    Ok(regions)
}

fn parse_err(msg: impl Into<String>) -> Box<dyn Error> {
    io::Error::new(io::ErrorKind::InvalidData, msg.into()).into()
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
    /// one pair with TLEN 200 on chr1 around pos 100, one pair with TLEN 300
    /// on chr1 around pos 5000 (outside the BED region below).
    const PAIRED_SAM: &str = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10000\n\
        r1\t99\tchr1\t101\t60\t50M\t=\t251\t200\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
        r1\t147\tchr1\t251\t60\t50M\t=\t101\t-200\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
        r2\t99\tchr1\t5001\t60\t50M\t=\t5251\t300\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
        r2\t147\tchr1\t5251\t60\t50M\t=\t5001\t-300\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";

    #[test]
    fn computes_mean_and_std_across_pairs() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_file(dir.path(), "reads.sam", PAIRED_SAM);

        let stats = InsertSizeStats::from_bam(&[path.as_str()], None).unwrap();
        assert_eq!(stats.n, 2);
        assert_eq!(stats.mean, 250.0);
        // Sample std of [200, 300] is 100/sqrt(2) * sqrt(2) = ... compute via
        // formula: variance = ((200-250)^2 + (300-250)^2) / (2-1) = 5000.
        assert!((stats.std - 5000f64.sqrt()).abs() < 1e-9);
    }

    #[test]
    fn restricts_to_bed_regions() {
        let dir = tempfile::tempdir().unwrap();
        let bam_path = write_file(dir.path(), "reads.sam", PAIRED_SAM);
        let bed_path = write_file(dir.path(), "regions.bed", "chr1\t0\t1000\n");

        let stats = InsertSizeStats::from_bam(&[bam_path.as_str()], Some(bed_path.as_str())).unwrap();
        assert_eq!(stats.n, 1);
        assert_eq!(stats.mean, 200.0);
        assert_eq!(stats.std, 0.0);
    }

    #[test]
    fn pools_reads_across_multiple_bam_files() {
        let dir = tempfile::tempdir().unwrap();
        let bam1 = write_file(dir.path(), "a.sam", PAIRED_SAM);
        let bam2 = write_file(dir.path(), "b.sam", PAIRED_SAM);

        let pooled = InsertSizeStats::from_bam(&[bam1.as_str(), bam2.as_str()], None).unwrap();
        assert_eq!(pooled.n, 4);
        assert_eq!(pooled.mean, 250.0);
    }

    #[test]
    fn bed_region_matching_no_reads_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let bam_path = write_file(dir.path(), "reads.sam", PAIRED_SAM);
        let bed_path = write_file(dir.path(), "regions.bed", "chr2\t0\t1000\n");

        assert!(InsertSizeStats::from_bam(&[bam_path.as_str()], Some(bed_path.as_str())).is_err());
    }

    #[test]
    fn empty_bam_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_file(
            dir.path(),
            "empty.sam",
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10000\n",
        );
        assert!(InsertSizeStats::from_bam(&[path.as_str()], None).is_err());
    }

    #[test]
    fn malformed_bed_line_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let bed_path = write_file(dir.path(), "regions.bed", "chr1\tnotanumber\t1000\n");
        assert!(parse_bed(&bed_path).is_err());
    }

    #[test]
    fn bed_comments_and_blank_lines_are_skipped() {
        let dir = tempfile::tempdir().unwrap();
        let bed_path = write_file(
            dir.path(),
            "regions.bed",
            "# comment\ntrack name=x\n\nchr1\t0\t1000\n",
        );
        let regions = parse_bed(&bed_path).unwrap();
        assert_eq!(regions["chr1"], vec![(0, 1000)]);
    }
}
