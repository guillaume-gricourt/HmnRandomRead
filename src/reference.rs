//! A single `-r/--reference` entry: a FASTA file, how many reads to draw
//! from it, and (optionally) which diversity profile to mutate them with.

use std::error::Error;
use std::path::Path;

use crate::io::FastaIndexedReader;
use crate::scaffold::Scaffolds;

pub struct Reference {
    pub path: String,
    pub faidx: FastaIndexedReader,
    pub scaffolds: Scaffolds,
    pub nb_reads: u64,
    pub nb_reads_remaining: u64,
    pub id_diversity: Option<String>,
}

impl Reference {
    /// Open `path`'s FASTA index and build its non-N scaffold set (contigs
    /// with no run of at least `min_scaffold_len` non-N bases contribute no
    /// scaffolds and are simply never picked).
    pub fn open(
        path: String,
        nb_reads: u64,
        id_diversity: Option<String>,
        min_scaffold_len: u64,
    ) -> Result<Self, Box<dyn Error>> {
        let faidx = FastaIndexedReader::open(&path)?;
        let scaffolds = Scaffolds::build(&faidx, min_scaffold_len)?;
        Ok(Reference {
            path,
            faidx,
            scaffolds,
            nb_reads,
            nb_reads_remaining: nb_reads,
            id_diversity,
        })
    }

    /// A reference's file base name, stripped of a `.fa`/`.fasta`/`.fna`
    /// extension, used as the "reference" field of a read's name.
    pub fn name(&self) -> String {
        let base = Path::new(&self.path)
            .file_name()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| self.path.clone());
        for ext in [".fasta", ".fna", ".fa"] {
            if let Some(stripped) = base.strip_suffix(ext) {
                return stripped.to_string();
            }
        }
        base
    }

    /// Parse a `-r/--reference` value: `path`, `path,nb_reads`, or
    /// `path,nb_reads,id_diversity`. `nb_reads` may be left empty (`path,`,
    /// `path,,id_diversity`) and then defaults to 0, so a reference brought
    /// in only to resolve a contig — a fusion partner's chromosome, say —
    /// can still carry an `id_diversity` without claiming a read count.
    pub fn parse_spec(spec: &str) -> Result<(String, u64, Option<String>), Box<dyn Error>> {
        let fields: Vec<&str> = spec.split(',').collect();
        match fields.as_slice() {
            [path] => Ok((path.to_string(), 0, None)),
            [path, nb_reads] => Ok((path.to_string(), Self::parse_nb_reads(nb_reads)?, None)),
            [path, nb_reads, id_diversity] => Ok((
                path.to_string(),
                Self::parse_nb_reads(nb_reads)?,
                Some(id_diversity.to_string()),
            )),
            _ => Err(format!("malformed --reference value: '{spec}'").into()),
        }
    }

    /// A spec's `nb_reads` field: an empty field means "unspecified", 0.
    fn parse_nb_reads(nb_reads: &str) -> Result<u64, Box<dyn Error>> {
        if nb_reads.is_empty() {
            return Ok(0);
        }
        Ok(nb_reads.parse()?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_spec_path_only() {
        let (path, nb_reads, id) = Reference::parse_spec("genome.fa").unwrap();
        assert_eq!(path, "genome.fa");
        assert_eq!(nb_reads, 0);
        assert_eq!(id, None);
    }

    #[test]
    fn parse_spec_path_and_count() {
        let (path, nb_reads, id) = Reference::parse_spec("genome.fa,100").unwrap();
        assert_eq!(path, "genome.fa");
        assert_eq!(nb_reads, 100);
        assert_eq!(id, None);
    }

    #[test]
    fn parse_spec_path_count_and_diversity() {
        let (path, nb_reads, id) = Reference::parse_spec("genome.fa,100,human").unwrap();
        assert_eq!(path, "genome.fa");
        assert_eq!(nb_reads, 100);
        assert_eq!(id.as_deref(), Some("human"));
    }

    #[test]
    fn parse_spec_path_and_empty_count() {
        let (path, nb_reads, id) = Reference::parse_spec("genome.fa,").unwrap();
        assert_eq!(path, "genome.fa");
        assert_eq!(nb_reads, 0);
        assert_eq!(id, None);
    }

    #[test]
    fn parse_spec_empty_count_and_diversity() {
        let (path, nb_reads, id) = Reference::parse_spec("genome.fa,,human").unwrap();
        assert_eq!(path, "genome.fa");
        assert_eq!(nb_reads, 0);
        assert_eq!(id.as_deref(), Some("human"));
    }

    #[test]
    fn parse_spec_rejects_too_many_fields() {
        assert!(Reference::parse_spec("genome.fa,100,human,extra").is_err());
    }

    #[test]
    fn parse_spec_rejects_non_numeric_count() {
        assert!(Reference::parse_spec("genome.fa,abc").is_err());
    }
}
