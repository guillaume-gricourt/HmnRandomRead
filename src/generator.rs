//! Orchestrates paired-end FASTQ read generation across one or more
//! references.

use std::error::Error;

use crate::diversity::ProfileDiversity;
use crate::fastq::FastqRecord;
use crate::io::FastqWriter;
use crate::profile_error::ProfileError;
use crate::reference::Reference;
use crate::rng::RandomGenerator;
use crate::sequence::Sequence;

const MAX_FRAGMENT_ATTEMPTS: usize = 10;

pub struct Config {
    pub length_reads: usize,
    pub mean_insert_size: f64,
    pub std_insert_size: f64,
    pub seed: u64,
    pub profile_diversity: Option<ProfileDiversity>,
    pub profile_error: Option<ProfileError>,
}

pub struct Generator {
    references: Vec<Reference>,
    config: Config,
}

impl Generator {
    /// Validates that every reference's `id_diversity` (if any) actually
    /// exists in `config.profile_diversity` up front, before any output is
    /// written — this is the fix for the original tool's segfault: it tried
    /// to run this same check, but `count(id) > 1` can never be true for a
    /// map's `count`, so an unknown id slipped through and only crashed the
    /// first time a read from that reference was mutated.
    pub fn new(references: Vec<Reference>, config: Config) -> Result<Self, Box<dyn Error>> {
        if let Some(profile) = &config.profile_diversity {
            for reference in &references {
                if let Some(id) = &reference.id_diversity {
                    if profile.get(id).is_none() {
                        return Err(format!(
                            "diversity id '{id}' for reference '{}' was not found in the \
                             profile diversity file",
                            reference.path
                        )
                        .into());
                    }
                }
            }
        }
        Ok(Generator { references, config })
    }

    pub fn run(mut self, out_forward: &str, out_reverse: &str) -> Result<(), Box<dyn Error>> {
        let mut writer1 = FastqWriter::create(out_forward)?;
        let mut writer2 = FastqWriter::create(out_reverse)?;
        let mut rng = RandomGenerator::new(self.config.seed);

        let total_reads: u64 = self.references.iter().map(|r| r.nb_reads).sum();
        let mut skipped = 0u64;

        for i in 0..total_reads {
            if self.references.is_empty() {
                log::warn!("no reference left to draw the remaining requested reads from");
                break;
            }
            let ref_idx = if self.references.len() > 1 {
                rng.range(0, self.references.len() - 1)
            } else {
                0
            };

            let fragment = pick_fragment(
                &self.references[ref_idx],
                &mut rng,
                self.config.mean_insert_size,
                self.config.std_insert_size,
            )?;
            let Some((chrom, location_start, location_stop, mut sequence)) = fragment else {
                skipped += 1;
                log::warn!("read {i}: no valid fragment found after {MAX_FRAGMENT_ATTEMPTS} attempts, skipping");
                continue;
            };

            let reference = &self.references[ref_idx];
            if let (Some(profile), Some(id)) =
                (&self.config.profile_diversity, &reference.id_diversity)
            {
                let diversity = profile
                    .get(id)
                    .expect("diversity ids are validated in Generator::new");
                sequence.make_mutation(&mut rng, diversity);
            }
            let ref_name = reference.name();

            // `head` always starts at `location_start` reading forward (no
            // reverse-complement); `tail` always ends at `location_stop`
            // reading backward (reverse-complemented). Which one lands in
            // the forward-strand output file vs the reverse one is
            // randomized below, but content and coordinates always stay
            // paired correctly by construction.
            let head_seq = sequence.sub_read(self.config.length_reads, true, false);
            let tail_seq = sequence.sub_read(self.config.length_reads, false, true);
            let head_len = head_seq.len() as u64;
            let tail_len = tail_seq.len() as u64;

            let mut head = FastqRecord::new(
                head_seq,
                33,
                i,
                true,
                ref_name.clone(),
                chrom.clone(),
                location_start,
                location_start + head_len,
            );
            let mut tail = FastqRecord::new(
                tail_seq,
                33,
                i,
                false,
                ref_name,
                chrom,
                location_stop.saturating_sub(tail_len),
                location_stop,
            );

            head.init_qual(&mut rng);
            tail.init_qual(&mut rng);
            if let Some(profile_error) = &self.config.profile_error {
                head.make_errors(&mut rng, profile_error);
                tail.make_errors(&mut rng, profile_error);
            }

            let (r1, r2) = if rng.unit() >= 0.5 {
                (head, tail)
            } else {
                (tail, head)
            };
            writer1.write_record(&r1)?;
            writer2.write_record(&r2)?;

            self.references[ref_idx].nb_reads_remaining -= 1;
            if self.references[ref_idx].nb_reads_remaining == 0 {
                self.references.remove(ref_idx);
            }
        }

        writer1.finish()?;
        writer2.finish()?;
        if skipped > 0 {
            log::warn!("{skipped} read(s) out of {total_reads} requested were skipped");
        }
        Ok(())
    }
}

type Fragment = (String, u64, u64, Sequence);

/// Try up to [`MAX_FRAGMENT_ATTEMPTS`] times to pick a length-weighted random
/// scaffold, draw an insert size, and fetch the resulting genomic interval.
/// Returns `Ok(None)` if every attempt failed (a degenerate insert size, an
/// unlucky scaffold boundary, or a fetch coming back empty).
fn pick_fragment(
    reference: &Reference,
    rng: &mut RandomGenerator,
    mean_insert_size: f64,
    std_insert_size: f64,
) -> std::io::Result<Option<Fragment>> {
    if reference.scaffolds.total_length == 0 {
        return Ok(None);
    }

    for _ in 0..MAX_FRAGMENT_ATTEMPTS {
        let offset = rng.range(0, reference.scaffolds.total_length - 1);
        let Some(scaffold) = reference.scaffolds.pick(offset) else {
            continue;
        };

        let size_insert = rng.normal(mean_insert_size, std_insert_size).round().max(0.0) as u64;
        let (location_start, location_stop) = if size_insert > scaffold.length() {
            (scaffold.start, scaffold.stop)
        } else {
            let start = rng.range(scaffold.start, scaffold.stop - size_insert);
            (start, start + size_insert)
        };
        // A size-0 insert (or any other degenerate draw) yields an empty
        // range; skip it rather than asking the FASTA reader to fetch an
        // inverted/empty region (undefined behavior in the underlying
        // htslib binding for out-of-range regions).
        if location_stop <= location_start {
            continue;
        }

        match reference
            .faidx
            .fetch(&scaffold.name, location_start, location_stop - 1)
        {
            Ok(seq) if !seq.is_empty() => {
                return Ok(Some((
                    scaffold.name.clone(),
                    location_start,
                    location_stop,
                    Sequence::new(seq),
                )))
            }
            _ => continue,
        }
    }
    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Read, Write};

    fn write_fasta(dir: &std::path::Path, name: &str, contents: &str) -> String {
        let path = dir.join(name);
        std::fs::File::create(&path)
            .unwrap()
            .write_all(contents.as_bytes())
            .unwrap();
        path.to_str().unwrap().to_string()
    }

    fn read_gz(path: &str) -> String {
        let file = std::fs::File::open(path).unwrap();
        let mut decoder = flate2::read::GzDecoder::new(file);
        let mut out = String::new();
        decoder.read_to_string(&mut out).unwrap();
        out
    }

    fn toy_fasta_body() -> String {
        format!(">chr1\n{}\n", "ACGT".repeat(100))
    }

    #[test]
    fn generates_the_requested_number_of_read_pairs() {
        let dir = tempfile::tempdir().unwrap();
        let fasta_path = write_fasta(dir.path(), "ref.fa", &toy_fasta_body());
        let reference = Reference::open(fasta_path, 20, None, 10).unwrap();

        let config = Config {
            length_reads: 50,
            mean_insert_size: 150.0,
            std_insert_size: 20.0,
            seed: 42,
            profile_diversity: None,
            profile_error: None,
        };
        let out1 = dir.path().join("r1.fastq.gz");
        let out2 = dir.path().join("r2.fastq.gz");
        Generator::new(vec![reference], config)
            .unwrap()
            .run(out1.to_str().unwrap(), out2.to_str().unwrap())
            .unwrap();

        let r1 = read_gz(out1.to_str().unwrap());
        let r2 = read_gz(out2.to_str().unwrap());
        assert_eq!(r1.lines().count(), 20 * 4);
        assert_eq!(r2.lines().count(), 20 * 4);
        for line in r1.lines().step_by(4) {
            assert!(line.starts_with('@'));
        }
    }

    #[test]
    fn same_seed_is_deterministic() {
        let dir = tempfile::tempdir().unwrap();
        let fasta_path = write_fasta(dir.path(), "ref.fa", &toy_fasta_body());

        let run = |seed: u64| {
            let reference = Reference::open(fasta_path.clone(), 10, None, 10).unwrap();
            let config = Config {
                length_reads: 30,
                mean_insert_size: 100.0,
                std_insert_size: 10.0,
                seed,
                profile_diversity: None,
                profile_error: None,
            };
            let out1 = dir.path().join(format!("a{seed}.fastq.gz"));
            let out2 = dir.path().join(format!("b{seed}.fastq.gz"));
            Generator::new(vec![reference], config)
                .unwrap()
                .run(out1.to_str().unwrap(), out2.to_str().unwrap())
                .unwrap();
            (read_gz(out1.to_str().unwrap()), read_gz(out2.to_str().unwrap()))
        };

        let first = run(7);
        let second = run(7);
        assert_eq!(first, second);
    }

    #[test]
    fn unknown_diversity_id_is_a_clean_error_not_a_crash() {
        let dir = tempfile::tempdir().unwrap();
        let fasta_path = write_fasta(dir.path(), "ref.fa", &toy_fasta_body());
        let reference =
            Reference::open(fasta_path, 5, Some("does-not-exist".to_string()), 10).unwrap();

        let profile_path = dir.path().join("diversity.csv");
        std::fs::File::create(&profile_path)
            .unwrap()
            .write_all(
                b"##Version: 1.0\n#identifiant,Mutation Rate,Indel Fraction,Indel Extend,Maximum Insertion Size\nhuman,0.001,0.15,0.3,15\n",
            )
            .unwrap();
        let profile = ProfileDiversity::parse_csv(&profile_path).unwrap();

        let config = Config {
            length_reads: 30,
            mean_insert_size: 100.0,
            std_insert_size: 10.0,
            seed: 1,
            profile_diversity: Some(profile),
            profile_error: None,
        };
        assert!(Generator::new(vec![reference], config).is_err());
    }
}
