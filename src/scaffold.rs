//! Non-N intervals within reference sequences, used to pick fragment
//! locations that don't land on stretches of unknown bases.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;

use crate::io::FastaIndexedReader;

/// A contiguous non-N interval `[start, stop)` within one contig.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Scaffold {
    pub name: String,
    pub start: u64,
    pub stop: u64,
}

impl Scaffold {
    pub fn length(&self) -> u64 {
        self.stop - self.start
    }
}

/// All non-N intervals of at least a minimum length across a reference's
/// contigs, plus a running total used to pick one at random weighted by
/// length.
#[derive(Clone, Debug, Default)]
pub struct Scaffolds {
    pub list: Vec<Scaffold>,
    pub intervals: Vec<u64>,
    pub total_length: u64,
}

impl Scaffolds {
    /// Build (or load a cached `<fasta>.scaff`) the set of non-N intervals at
    /// least `min_len` bases long across every contig of `fasta`.
    pub fn build(fasta: &FastaIndexedReader, min_len: u64) -> io::Result<Self> {
        let cache_path = format!("{}.scaff", fasta.path().display());
        if Path::new(&cache_path).is_file() {
            return Self::load(&cache_path);
        }

        let mut list = Vec::new();
        for name in fasta.seq_names()? {
            let seq = fasta.fetch_all(&name)?;
            list.extend(find_non_n_runs(&name, seq.as_bytes(), min_len));
        }
        let scaffolds = Self::from_list(list);
        scaffolds.save(&cache_path)?;
        Ok(scaffolds)
    }

    fn from_list(list: Vec<Scaffold>) -> Self {
        let mut intervals = Vec::with_capacity(list.len());
        let mut total_length = 0u64;
        for s in &list {
            total_length += s.length();
            intervals.push(total_length);
        }
        Scaffolds {
            list,
            intervals,
            total_length,
        }
    }

    fn save(&self, path: &str) -> io::Result<()> {
        let mut f = File::create(path)?;
        for s in &self.list {
            writeln!(f, "{}\t{}\t{}", s.name, s.start, s.stop)?;
        }
        Ok(())
    }

    fn load(path: &str) -> io::Result<Self> {
        let file = File::open(path)?;
        let mut list = Vec::new();
        for line in BufReader::new(file).lines() {
            let line = line?;
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() != 3 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("malformed scaffold cache line: '{line}'"),
                ));
            }
            let parse = |s: &str| {
                s.parse::<u64>().map_err(|e| {
                    io::Error::new(io::ErrorKind::InvalidData, e.to_string())
                })
            };
            list.push(Scaffold {
                name: fields[0].to_string(),
                start: parse(fields[1])?,
                stop: parse(fields[2])?,
            });
        }
        Ok(Self::from_list(list))
    }

    /// The scaffold covering a length-weighted random offset `[0, total_length)`.
    pub fn pick(&self, offset: u64) -> Option<&Scaffold> {
        let idx = self.intervals.partition_point(|&cum| cum <= offset);
        self.list.get(idx)
    }
}

/// Find every run of non-`N` bases at least `min_len` long in `seq`, as
/// scaffolds named `name`.
fn find_non_n_runs(name: &str, seq: &[u8], min_len: u64) -> Vec<Scaffold> {
    let mut runs = Vec::new();
    let mut run_start: Option<usize> = None;

    for (i, &b) in seq.iter().enumerate() {
        let is_n = b == b'N' || b == b'n';
        match (is_n, run_start) {
            (false, None) => run_start = Some(i),
            (true, Some(start)) => {
                push_run(&mut runs, name, start, i, min_len);
                run_start = None;
            }
            _ => {}
        }
    }
    if let Some(start) = run_start {
        push_run(&mut runs, name, start, seq.len(), min_len);
    }
    runs
}

fn push_run(runs: &mut Vec<Scaffold>, name: &str, start: usize, stop: usize, min_len: u64) {
    let len = (stop - start) as u64;
    if len >= min_len {
        runs.push(Scaffold {
            name: name.to_string(),
            start: start as u64,
            stop: stop as u64,
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_non_n_runs_above_min_length() {
        // Runs: [3,11)="ACGTACGT" (len 8, kept), [15,17)="AC" (len 2,
        // filtered by min_len=4), [20,30)="TTTTTTTTTT" (len 10, kept).
        let seq = b"NNNACGTACGTNNNNACNNNTTTTTTTTTTNN";
        let runs = find_non_n_runs("chr1", seq, 4);
        assert_eq!(
            runs,
            vec![
                Scaffold { name: "chr1".into(), start: 3, stop: 11 },
                Scaffold { name: "chr1".into(), start: 20, stop: 30 },
            ]
        );
    }

    #[test]
    fn whole_sequence_is_one_run_when_no_n() {
        let runs = find_non_n_runs("chr1", b"ACGTACGT", 1);
        assert_eq!(runs, vec![Scaffold { name: "chr1".into(), start: 0, stop: 8 }]);
    }

    #[test]
    fn pick_finds_scaffold_covering_offset() {
        let scaffolds = Scaffolds::from_list(vec![
            Scaffold { name: "a".into(), start: 0, stop: 10 },
            Scaffold { name: "b".into(), start: 0, stop: 5 },
        ]);
        assert_eq!(scaffolds.total_length, 15);
        assert_eq!(scaffolds.pick(0).unwrap().name, "a");
        assert_eq!(scaffolds.pick(9).unwrap().name, "a");
        assert_eq!(scaffolds.pick(10).unwrap().name, "b");
        assert_eq!(scaffolds.pick(14).unwrap().name, "b");
        assert!(scaffolds.pick(15).is_none());
    }

    #[test]
    fn save_and_load_roundtrip() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("x.scaff");
        let scaffolds = Scaffolds::from_list(vec![Scaffold {
            name: "chr1".into(),
            start: 5,
            stop: 20,
        }]);
        scaffolds.save(path.to_str().unwrap()).unwrap();
        let loaded = Scaffolds::load(path.to_str().unwrap()).unwrap();
        assert_eq!(loaded.list, scaffolds.list);
        assert_eq!(loaded.total_length, scaffolds.total_length);
    }
}
