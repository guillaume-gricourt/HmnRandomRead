//! Per-cycle sequencer error model and its CSV profile file.

use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

fn parse_err(msg: impl Into<String>) -> Box<dyn Error> {
    io::Error::new(io::ErrorKind::InvalidData, msg.into()).into()
}

/// `0` = reverse strand, `1` = forward strand — matches [`crate::fastq::FastqRecord::mate`].
pub const REVERSE: u8 = 0;
pub const FORWARD: u8 = 1;

/// A parsed `-profileError` CSV file: per-cycle substitution error rate for
/// one or both strands of a single sequencer/id.
#[derive(Clone, Debug, Default)]
pub struct ProfileError {
    pub version: f32,
    errors: HashMap<u8, Vec<f32>>,
}

impl ProfileError {
    /// Parse a profile error CSV for the row(s) matching `id`:
    /// ```text
    /// ##Version: 1.0
    /// #identifiant,sequencer,flowcell,version,strand,cycles total,error by cycle
    /// n1,nextseq 550,high-output,NA,forward,150,2.34;0.23;...
    /// ```
    /// If `is_paired` and only one strand was found for `id`, the same error
    /// curve is reused for the other strand (with a warning), matching the
    /// original tool's fallback for single-strand profiles applied to
    /// paired-end output.
    pub fn parse_csv<P: AsRef<Path>>(
        path: P,
        id: &str,
        is_paired: bool,
    ) -> Result<Self, Box<dyn Error>> {
        let file = File::open(path)?;
        let mut profile = ProfileError::default();
        let mut last_strand = FORWARD;

        for (i, line) in BufReader::new(file).lines().enumerate() {
            let line = line?;
            match i {
                0 => {
                    let fields: Vec<&str> = line.split(' ').collect();
                    if fields.len() != 2 || !fields[0].starts_with("##") {
                        return Err(parse_err("profile error version line malformatted"));
                    }
                    profile.version = fields[1].parse()?;
                }
                1 => {
                    let fields: Vec<&str> = line.split(',').collect();
                    if fields.len() != 7 || !fields[0].starts_with('#') {
                        return Err(parse_err("profile error header malformatted"));
                    }
                }
                _ => {
                    let fields: Vec<&str> = line.split(',').collect();
                    if fields.len() != 7 || fields[0] != id {
                        continue;
                    }
                    let strand = match fields[4] {
                        "reverse" => REVERSE,
                        "forward" => FORWARD,
                        other => {
                            return Err(parse_err(format!("strand unknown: '{other}'")));
                        }
                    };
                    let cycle_total: usize = fields[5].parse()?;
                    let cycle_error: Vec<f32> = fields[6]
                        .split(';')
                        .map(|s| s.parse::<f32>().map(|v| v / 100.0))
                        .collect::<Result<_, _>>()?;
                    if cycle_total != cycle_error.len() {
                        return Err(parse_err(
                            "cycle total discordant with cycles indicated",
                        ));
                    }
                    last_strand = strand;
                    profile.errors.entry(strand).or_insert(cycle_error);
                }
            }
        }

        if is_paired {
            if profile.errors.len() == 1 {
                eprintln!(
                    "Output is paired but only one strand is found in profile file error, \
                     use the same error for the 2 strand"
                );
                let missing = 1 - last_strand;
                let existing = profile.errors.get(&last_strand).unwrap().clone();
                profile.errors.insert(missing, existing);
            } else if profile.errors.len() != 2 {
                return Err(parse_err(format!(
                    "too many identifiers were found for '{id}'"
                )));
            }
        } else if profile.errors.len() == 1 && last_strand == FORWARD {
            return Err(parse_err(
                "only strand reverse error profile was found, only forward profile error \
                 for forward fastq to produce must be indicated",
            ));
        } else if profile.errors.len() > 1 {
            return Err(parse_err(format!(
                "too many identifiers were found for '{id}'"
            )));
        }

        Ok(profile)
    }

    /// Error rate for `strand` at 0-based `cycle`, if the profile covers it.
    pub fn rate(&self, strand: u8, cycle: usize) -> Option<f32> {
        self.errors.get(&strand)?.get(cycle).copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_profile(dir: &Path, contents: &str) -> std::path::PathBuf {
        let path = dir.join("profile_error.csv");
        let mut f = File::create(&path).unwrap();
        f.write_all(contents.as_bytes()).unwrap();
        path
    }

    const HEADER: &str = "##Version: 1.0\n#identifiant,sequencer,flowcell,version,strand,cycles total,error by cycle\n";

    #[test]
    fn parses_paired_profile_with_both_strands() {
        let dir = tempfile::tempdir().unwrap();
        let contents = format!(
            "{HEADER}n1,seq,flow,NA,forward,3,1;2;3\nn1,seq,flow,NA,reverse,3,4;5;6\n"
        );
        let path = write_profile(dir.path(), &contents);
        let profile = ProfileError::parse_csv(&path, "n1", true).unwrap();
        assert_eq!(profile.rate(FORWARD, 0).unwrap(), 0.01);
        assert_eq!(profile.rate(REVERSE, 2).unwrap(), 0.06);
    }

    #[test]
    fn single_strand_is_duplicated_for_paired_output() {
        let dir = tempfile::tempdir().unwrap();
        let contents = format!("{HEADER}n1,seq,flow,NA,reverse,2,1;2\n");
        let path = write_profile(dir.path(), &contents);
        let profile = ProfileError::parse_csv(&path, "n1", true).unwrap();
        assert_eq!(profile.rate(REVERSE, 0), profile.rate(FORWARD, 0));
    }

    #[test]
    fn single_forward_only_profile_rejected_for_single_end_output() {
        let dir = tempfile::tempdir().unwrap();
        let contents = format!("{HEADER}n1,seq,flow,NA,forward,2,1;2\n");
        let path = write_profile(dir.path(), &contents);
        assert!(ProfileError::parse_csv(&path, "n1", false).is_err());
    }

    #[test]
    fn cycle_count_mismatch_is_rejected() {
        let dir = tempfile::tempdir().unwrap();
        let contents = format!("{HEADER}n1,seq,flow,NA,forward,5,1;2\n");
        let path = write_profile(dir.path(), &contents);
        assert!(ProfileError::parse_csv(&path, "n1", false).is_err());
    }
}
