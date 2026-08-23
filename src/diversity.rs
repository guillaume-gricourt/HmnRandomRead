//! SNP/indel diversity model and its CSV profile file.

use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

/// Per-reference mutation rates applied by [`crate::sequence::Sequence::make_mutation`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Diversity {
    pub mut_rate: f64,
    pub indel_frac: f64,
    pub indel_extend: f64,
    pub max_ins_size: usize,
}

fn parse_err(msg: impl Into<String>) -> Box<dyn Error> {
    io::Error::new(io::ErrorKind::InvalidData, msg.into()).into()
}

/// A parsed `-profileDiversity` CSV file: one [`Diversity`] per identifier.
#[derive(Clone, Debug, Default)]
pub struct ProfileDiversity {
    pub version: f32,
    data: HashMap<String, Diversity>,
}

impl ProfileDiversity {
    /// Parse a profile diversity CSV:
    /// ```text
    /// ##Version: 1.0
    /// #identifiant,Mutation Rate,Indel Fraction,Indel Extend,Maximum Insertion Size
    /// human,0.001,0.15,0.3,15
    /// ```
    pub fn parse_csv<P: AsRef<Path>>(path: P) -> Result<Self, Box<dyn Error>> {
        let file = File::open(path)?;
        let mut profile = ProfileDiversity::default();

        for (i, line) in BufReader::new(file).lines().enumerate() {
            let line = line?;
            match i {
                0 => {
                    let fields: Vec<&str> = line.split(' ').collect();
                    if fields.len() != 2 || !fields[0].starts_with("##") {
                        return Err(parse_err("profile diversity version line malformatted"));
                    }
                    profile.version = fields[1].parse()?;
                }
                1 => {
                    let fields: Vec<&str> = line.split(',').collect();
                    if fields.len() != 5 || !fields[0].starts_with('#') {
                        return Err(parse_err("profile diversity header malformatted"));
                    }
                }
                _ => {
                    let fields: Vec<&str> = line.split(',').collect();
                    if fields.len() != 5 {
                        return Err(parse_err(format!(
                            "profile diversity line {} malformatted",
                            i + 1
                        )));
                    }
                    let id = fields[0].to_string();
                    if profile.data.contains_key(&id) {
                        return Err(parse_err(format!("identifier '{id}' isn't unique")));
                    }
                    let diversity = Diversity {
                        mut_rate: fields[1].parse()?,
                        indel_frac: fields[2].parse()?,
                        indel_extend: fields[3].parse()?,
                        max_ins_size: fields[4].parse()?,
                    };
                    profile.data.insert(id, diversity);
                }
            }
        }
        Ok(profile)
    }

    /// Look up the diversity model for a reference's `id_diversity`. Returns
    /// `None` for an unknown id instead of ever fabricating a default entry
    /// (the C++ version used `map::operator[]`, which silently inserted a
    /// null entry for missing ids and crashed the first time it was used).
    pub fn get(&self, id: &str) -> Option<&Diversity> {
        self.data.get(id)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_profile(dir: &Path, contents: &str) -> std::path::PathBuf {
        let path = dir.join("profile_diversity.csv");
        let mut f = File::create(&path).unwrap();
        f.write_all(contents.as_bytes()).unwrap();
        path
    }

    const VALID: &str = "##Version: 1.0\n#identifiant,Mutation Rate,Indel Fraction,Indel Extend,Maximum Insertion Size\nhuman,0.001,0.15,0.3,15\nbacteria,0.4,0.15,0.45,15\n";

    #[test]
    fn parses_valid_profile() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_profile(dir.path(), VALID);
        let profile = ProfileDiversity::parse_csv(&path).unwrap();
        assert_eq!(profile.version, 1.0);
        let human = profile.get("human").unwrap();
        assert_eq!(human.mut_rate, 0.001);
        assert_eq!(human.max_ins_size, 15);
    }

    #[test]
    fn unknown_id_returns_none_not_a_crash() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_profile(dir.path(), VALID);
        let profile = ProfileDiversity::parse_csv(&path).unwrap();
        assert!(profile.get("does-not-exist").is_none());
    }

    #[test]
    fn duplicate_identifier_is_rejected() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_profile(
            dir.path(),
            "##Version: 1.0\n#identifiant,Mutation Rate,Indel Fraction,Indel Extend,Maximum Insertion Size\nhuman,0.001,0.15,0.3,15\nhuman,0.002,0.1,0.2,10\n",
        );
        assert!(ProfileDiversity::parse_csv(&path).is_err());
    }

    #[test]
    fn malformed_version_line_is_rejected() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_profile(dir.path(), "not-a-version-line\n");
        assert!(ProfileDiversity::parse_csv(&path).is_err());
    }
}
