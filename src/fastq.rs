//! A single simulated FASTQ record: sequence, quality, and read-pair metadata.

use std::fmt;

use crate::profile_sequencer::ProfileSequencer;
use crate::rng::RandomGenerator;
use crate::sequence::choose_base;

/// One mate of a simulated read pair.
#[derive(Clone, Debug)]
pub struct FastqRecord {
    pub sequence: String,
    pub quality: Vec<u8>,
    pub phred_offset: u8,
    pub number: u64,
    /// `true` for the read taken from the forward strand of the fragment.
    pub forward: bool,
    pub reference: String,
    pub chromosome: String,
    pub start: u64,
    pub end: u64,
}

impl FastqRecord {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        sequence: String,
        phred_offset: u8,
        number: u64,
        forward: bool,
        reference: String,
        chromosome: String,
        start: u64,
        end: u64,
    ) -> Self {
        let quality = vec![0u8; sequence.len()];
        FastqRecord {
            sequence,
            quality,
            phred_offset,
            number,
            forward,
            reference,
            chromosome,
            start,
            end,
        }
    }

    /// `0` for the forward-strand mate, `1` for the other — matches
    /// [`crate::profile_sequencer::FORWARD`]/[`crate::profile_sequencer::REVERSE`].
    pub fn mate(&self) -> u8 {
        if self.forward {
            0
        } else {
            1
        }
    }

    fn name(&self) -> String {
        format!(
            "{}/{} {}_{}_{}_{}",
            self.number,
            self.mate(),
            self.reference,
            self.chromosome,
            self.start,
            self.end
        )
    }

    /// Draw a baseline Phred quality score (29-36) for every base.
    pub fn init_qual(&mut self, rng: &mut RandomGenerator) {
        for q in self.quality.iter_mut() {
            *q = rng.range(29u8, 36u8);
        }
    }

    /// Introduce sequencing errors per-cycle according to `profile`, mutating
    /// both the base and its quality score wherever an error is drawn.
    pub fn make_errors(&mut self, rng: &mut RandomGenerator, profile: &ProfileSequencer) {
        let strand = self.mate();
        let mut bytes = std::mem::take(&mut self.sequence).into_bytes();
        for (i, base) in bytes.iter_mut().enumerate() {
            let Some(error_rate) = profile.rate(strand, i) else {
                continue;
            };
            if rng.unit() < error_rate as f64 {
                *base = choose_base(rng, Some(*base));
                self.quality[i] = error_rate_to_qscore(error_rate);
            }
        }
        self.sequence = String::from_utf8(bytes).expect("sequence bytes must remain valid UTF-8");
    }
}

fn error_rate_to_qscore(error_rate: f32) -> u8 {
    (-10.0 * error_rate.log10() + 0.5) as u8
}

impl fmt::Display for FastqRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "@{}", self.name())?;
        writeln!(f, "{}", self.sequence)?;
        writeln!(f, "+")?;
        for q in &self.quality {
            write!(f, "{}", (q + self.phred_offset) as char)?;
        }
        writeln!(f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn to_string_renders_four_fastq_lines() {
        let mut record = FastqRecord::new(
            "ACGT".to_string(),
            33,
            7,
            true,
            "genome".to_string(),
            "chr1".to_string(),
            100,
            104,
        );
        let mut rng = RandomGenerator::new(1);
        record.init_qual(&mut rng);

        let text = record.to_string();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 4);
        assert_eq!(lines[0], "@7/0 genome_chr1_100_104");
        assert_eq!(lines[1], "ACGT");
        assert_eq!(lines[2], "+");
        assert_eq!(lines[3].len(), 4);
    }

    #[test]
    fn mate_matches_strand() {
        let fwd = FastqRecord::new(
            "A".into(),
            33,
            0,
            true,
            "r".into(),
            "c".into(),
            0,
            1,
        );
        let rev = FastqRecord::new(
            "A".into(),
            33,
            0,
            false,
            "r".into(),
            "c".into(),
            0,
            1,
        );
        assert_eq!(fwd.mate(), 0);
        assert_eq!(rev.mate(), 1);
    }
}
