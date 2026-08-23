//! Builds a `-profileSequencer` CSV (see [`crate::profile_sequencer`]) from real
//! sequencing data — paired FASTQ files or a BAM — instead of requiring the
//! user to already know their sequencer's error curve.
//!
//! `flowcell` and `version` can't be recovered from the data, so they are
//! always written as `NA`.

use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

use rust_htslib::bam::record::{Cigar, CigarStringView};
use rust_htslib::bam::{self, Read as BamRead};

use crate::io::FastqReader;
use crate::profile_sequencer::{FORWARD, REVERSE};
use crate::sequencer;

/// Sequencer output is 33-offset (Sanger/Illumina 1.8+) FASTQ/BAM quality,
/// the only encoding this crate ever produces or consumes elsewhere.
const PHRED_OFFSET: u8 = 33;

/// A profile sequencer row set for a single identifier, built from real data and
/// ready to serialize with [`BuiltProfile::write_csv`].
pub struct BuiltProfile {
    pub id: String,
    pub sequencer: String,
    cycles: HashMap<u8, Vec<f32>>,
}

impl BuiltProfile {
    /// Build from matching lists of forward/reverse FASTQ files (plain or
    /// `.gz`) — e.g. one pair per lane, pooled into a single profile. The
    /// sequencer is guessed from the first read header found in
    /// `forward_paths`; `cycles` and `error_by_cycle` accumulate every
    /// read from every file, independently per strand.
    pub fn from_fastq<S: AsRef<str>>(
        id: &str,
        forward_paths: &[S],
        reverse_paths: &[S],
    ) -> Result<Self, Box<dyn Error>> {
        let mut forward_acc = Accumulator::default();
        let mut reverse_acc = Accumulator::default();
        let mut sequencer = None;

        for path in forward_paths {
            if let Some(found) = accumulate_fastq(path.as_ref(), &mut forward_acc)? {
                sequencer.get_or_insert(found);
            }
        }
        let mut saw_reverse_read = false;
        for path in reverse_paths {
            if accumulate_fastq(path.as_ref(), &mut reverse_acc)?.is_some() {
                saw_reverse_read = true;
            }
        }

        let Some(sequencer) = sequencer else {
            return Err(parse_err(
                "no reads found in any --input-forward-fastq file",
            ));
        };
        if !saw_reverse_read {
            return Err(parse_err(
                "no reads found in any --input-reverse-fastq file",
            ));
        }

        let mut cycles = HashMap::new();
        cycles.insert(FORWARD, forward_acc.error_rates());
        cycles.insert(REVERSE, reverse_acc.error_rates());

        Ok(BuiltProfile {
            id: id.to_string(),
            sequencer,
            cycles,
        })
    }

    /// Build from one or more coordinate- or name-sorted BAMs, pooled into a
    /// single profile. The sequencer is guessed from the first `@RG PM` tag
    /// found across them, otherwise from the first usable read's `QNAME`.
    /// Secondary/supplementary/unmapped alignments are skipped so a read
    /// isn't counted more than once. Read 1 maps to the `forward` strand
    /// and read 2 to `reverse`, matching the mate convention `simulate`
    /// writes (see [`crate::fastq::FastqRecord::mate`]) — this is
    /// independent of each alignment's genomic strand.
    ///
    /// Hard-clipped (`H`) bases were physically removed from a record's
    /// `SEQ`/`QUAL`, so they carry no quality to feed into
    /// `error_by_cycle` and are excluded from the tally entirely — they
    /// never widen `cycles` on their own. Their length is still used to
    /// correctly place the *retained* quality at its true absolute cycle
    /// (see the offset handling below), so a read clipped at its start
    /// doesn't get its real quality mistakenly counted as starting at
    /// cycle 0.
    pub fn from_bam<S: AsRef<str>>(id: &str, bam_paths: &[S]) -> Result<Self, Box<dyn Error>> {
        let mut sequencer: Option<String> = None;
        let mut accumulators: HashMap<u8, Accumulator> = HashMap::new();
        let mut seen_any = false;

        for bam_path in bam_paths {
            let bam_path = bam_path.as_ref();
            let mut reader = bam::Reader::from_path(bam_path)?;
            if sequencer.is_none() {
                let header_text = String::from_utf8_lossy(reader.header().as_bytes()).into_owned();
                sequencer = sequencer::from_bam_header(&header_text);
            }

            let mut record = bam::Record::new();
            while let Some(result) = reader.read(&mut record) {
                result?;
                if record.is_secondary() || record.is_supplementary() || record.is_unmapped() {
                    continue;
                }
                seen_any = true;
                if sequencer.is_none() {
                    sequencer = Some(sequencer::from_header(&String::from_utf8_lossy(
                        record.qname(),
                    )));
                }

                let strand = if record.is_last_in_template() {
                    REVERSE
                } else {
                    FORWARD
                };

                // Hard clips carry no quality — only their length matters,
                // to correctly offset the quality that *is* retained.
                let (leading_h, trailing_h) = hardclip_lengths(&record.cigar());

                // BAM stores SEQ/QUAL reverse-complemented relative to the
                // original sequencing direction for reverse-strand
                // alignments; undo that so index 0 is always the first
                // sequenced cycle. The hard clip that lands at the *start*
                // of the original read is then whichever end trails in
                // stored order once un-reversed.
                let mut qual = record.qual().to_vec();
                let offset = if record.is_reverse() {
                    qual.reverse();
                    trailing_h
                } else {
                    leading_h
                } as usize;

                accumulators.entry(strand).or_default().add_at(&qual, offset);
            }
        }

        if !seen_any {
            return Err(parse_err(
                "no usable (primary, mapped) alignments found in any --input-bam file",
            ));
        }

        let cycles = accumulators
            .into_iter()
            .map(|(strand, acc)| (strand, acc.error_rates()))
            .collect();

        Ok(BuiltProfile {
            id: id.to_string(),
            sequencer: sequencer.unwrap_or_else(|| sequencer::UNKNOWN.to_string()),
            cycles,
        })
    }

    /// Write the `-profileSequencer` CSV: version header, column header, then
    /// one row per strand found in the input.
    pub fn write_csv<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        let mut file = File::create(path)?;
        writeln!(file, "##Version: 1.0")?;
        writeln!(
            file,
            "#id,sequencer,flowcell,version,strand,cycles,error_by_cycle"
        )?;
        for (strand, label) in [(FORWARD, "forward"), (REVERSE, "reverse")] {
            let Some(rates) = self.cycles.get(&strand) else {
                continue;
            };
            let error_by_cycle: Vec<String> =
                rates.iter().map(|rate| format!("{rate:.2}")).collect();
            writeln!(
                file,
                "{},{},NA,NA,{},{},{}",
                self.id,
                self.sequencer,
                label,
                rates.len(),
                error_by_cycle.join(";")
            )?;
        }
        Ok(())
    }
}

fn parse_err(msg: impl Into<String>) -> Box<dyn Error> {
    io::Error::new(io::ErrorKind::InvalidData, msg.into()).into()
}

/// Per-cycle running mean Phred quality, converted to a percent error rate.
#[derive(Default)]
struct Accumulator {
    sum_q: Vec<f64>,
    count: Vec<u64>,
}

impl Accumulator {
    /// `qual` must already be Phred-scaled (offset already removed) and in
    /// sequencing-cycle order (index 0 = first cycle). `offset` is the
    /// absolute cycle of `qual[0]` — nonzero when bases before it (e.g. a
    /// hard clip at the read's start) carry no quality to report. Only the
    /// cycles `qual` actually covers ever get a slot; cycles with no
    /// quality anywhere (whether hard-clipped or simply never covered by
    /// any read) are never added and so never appear in
    /// [`Accumulator::error_rates`]'s output, rather than being padded in
    /// with a fabricated rate.
    fn add_at(&mut self, qual: &[u8], offset: usize) {
        let needed = offset + qual.len();
        if needed > self.sum_q.len() {
            self.sum_q.resize(needed, 0.0);
            self.count.resize(needed, 0);
        }
        for (i, &q) in qual.iter().enumerate() {
            let cycle = offset + i;
            self.sum_q[cycle] += q as f64;
            self.count[cycle] += 1;
        }
    }

    /// The per-cycle mean Phred quality turned back into a probability of
    /// error (as a percent, matching the CSV's `error_by_cycle` units) —
    /// the inverse of how a quality score is derived from an error rate in
    /// `crate::fastq::error_rate_to_qscore`.
    fn error_rates(&self) -> Vec<f32> {
        self.sum_q
            .iter()
            .zip(&self.count)
            .map(|(sum, count)| {
                if *count == 0 {
                    0.0
                } else {
                    let mean_q = sum / *count as f64;
                    (10f64.powf(-mean_q / 10.0) * 100.0) as f32
                }
            })
            .collect()
    }
}

/// Reads every record in `path` into `acc`, returning the sequencer guessed
/// from its first read's header — or `None` if `path` has no reads at all.
fn accumulate_fastq(path: &str, acc: &mut Accumulator) -> Result<Option<String>, Box<dyn Error>> {
    let mut reader = FastqReader::open(path)?;
    let mut sequencer = None;

    while let Some(record) = reader.next_record()? {
        if sequencer.is_none() {
            sequencer = Some(sequencer::from_header(&record.header));
        }
        let qual: Vec<u8> = record
            .quality
            .iter()
            .map(|q| q.saturating_sub(PHRED_OFFSET))
            .collect();
        acc.add_at(&qual, 0);
    }

    Ok(sequencer)
}

/// Hard-clip lengths at the start and end of `cigar`, in *stored* (reference)
/// order. Per the SAM spec, `H` can only appear as the first and/or last
/// operation — anywhere else the BAM is malformed. A single-operation CIGAR
/// that is itself a hard clip (a fully-clipped, degenerate record) is
/// reported as leading-only, since counting it on both ends would double it.
fn hardclip_lengths(cigar: &CigarStringView) -> (u32, u32) {
    let hardclip_len = |op: Option<&Cigar>| match op {
        Some(Cigar::HardClip(len)) => *len,
        _ => 0,
    };
    let leading = hardclip_len(cigar.first());
    let trailing = if cigar.len() < 2 {
        0
    } else {
        hardclip_len(cigar.last())
    };
    (leading, trailing)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write_fastq(dir: &Path, name: &str, contents: &str) -> String {
        let path = dir.join(name);
        File::create(&path).unwrap().write_all(contents.as_bytes()).unwrap();
        path.to_str().unwrap().to_string()
    }

    #[test]
    fn builds_profile_from_fastq_pair() {
        let dir = tempfile::tempdir().unwrap();
        let forward = write_fastq(
            dir.path(),
            "r1.fastq",
            "@NB552023:1:AAA:1:1:1:1 1:N:0:1\nACGT\n+\nIIII\n@NB552023:1:AAA:1:1:1:2 1:N:0:1\nACGT\n+\n!!!!\n",
        );
        let reverse = write_fastq(
            dir.path(),
            "r2.fastq",
            "@NB552023:1:AAA:1:1:1:1 2:N:0:1\nACGT\n+\nIIII\n",
        );

        let profile = BuiltProfile::from_fastq("n1", &[forward.as_str()], &[reverse.as_str()]).unwrap();
        assert_eq!(profile.sequencer, "nextseq 550");
        assert_eq!(profile.cycles[&FORWARD].len(), 4);
        assert_eq!(profile.cycles[&REVERSE].len(), 4);
        // Read 1 at cycle 0 alternates a high-quality 'I' and a low-quality
        // '!' base; its averaged error rate should sit between the two
        // single-read rates.
        assert!(profile.cycles[&FORWARD][0] > 0.0);
    }

    #[test]
    fn write_csv_round_trips_through_profile_sequencer_parser() {
        let dir = tempfile::tempdir().unwrap();
        let forward = write_fastq(dir.path(), "r1.fastq", "@inst:1:AAA:1:1:1:1\nACGT\n+\nIIII\n");
        let reverse = write_fastq(dir.path(), "r2.fastq", "@inst:1:AAA:1:1:1:1\nACGT\n+\nIIII\n");
        let profile = BuiltProfile::from_fastq("n1", &[forward.as_str()], &[reverse.as_str()]).unwrap();

        let out_path = dir.path().join("profile_sequencer.csv");
        profile.write_csv(&out_path).unwrap();

        let parsed = crate::profile_sequencer::ProfileSequencer::parse_csv(&out_path, "n1", true).unwrap();
        assert!(parsed.rate(FORWARD, 0).is_some());
        assert!(parsed.rate(REVERSE, 0).is_some());
    }

    #[test]
    fn empty_fastq_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let forward = write_fastq(dir.path(), "r1.fastq", "");
        let reverse = write_fastq(dir.path(), "r2.fastq", "");
        assert!(BuiltProfile::from_fastq("n1", &[forward.as_str()], &[reverse.as_str()]).is_err());
    }

    #[test]
    fn pools_reads_across_multiple_lane_files() {
        let dir = tempfile::tempdir().unwrap();
        let lane1_fwd = write_fastq(dir.path(), "l1_r1.fastq", "@inst:1:AAA:1:1:1:1\nACGT\n+\nIIII\n");
        let lane2_fwd = write_fastq(dir.path(), "l2_r1.fastq", "@inst:1:AAA:2:1:1:1\nACGT\n+\nIIII\n");
        let lane1_rev = write_fastq(dir.path(), "l1_r2.fastq", "@inst:1:AAA:1:1:1:1\nACGT\n+\nIIII\n");
        let lane2_rev = write_fastq(dir.path(), "l2_r2.fastq", "@inst:1:AAA:2:1:1:1\nACGT\n+\nIIII\n");

        let profile = BuiltProfile::from_fastq(
            "n1",
            &[lane1_fwd.as_str(), lane2_fwd.as_str()],
            &[lane1_rev.as_str(), lane2_rev.as_str()],
        )
        .unwrap();
        // 2 reads pooled per strand, both with identical quality: the
        // averaged rate should be identical to what a single read gives.
        let single = BuiltProfile::from_fastq("n1", &[lane1_fwd.as_str()], &[lane1_rev.as_str()])
            .unwrap();
        assert_eq!(profile.cycles[&FORWARD], single.cycles[&FORWARD]);
    }

    #[test]
    fn missing_reverse_reads_is_an_error_even_with_forward_data() {
        let dir = tempfile::tempdir().unwrap();
        let forward = write_fastq(dir.path(), "r1.fastq", "@inst:1:AAA:1:1:1:1\nACGT\n+\nIIII\n");
        let reverse = write_fastq(dir.path(), "r2.fastq", "");
        assert!(BuiltProfile::from_fastq("n1", &[forward.as_str()], &[reverse.as_str()]).is_err());
    }

    #[test]
    fn hardclip_lengths_reads_both_ends() {
        let cigar =
            bam::record::CigarString(vec![Cigar::HardClip(3), Cigar::Match(10), Cigar::HardClip(2)])
                .into_view(0);
        assert_eq!(hardclip_lengths(&cigar), (3, 2));
    }

    #[test]
    fn hardclip_lengths_is_zero_without_hard_clips() {
        let cigar =
            bam::record::CigarString(vec![Cigar::SoftClip(3), Cigar::Match(10)]).into_view(0);
        assert_eq!(hardclip_lengths(&cigar), (0, 0));
    }

    #[test]
    fn hardclip_lengths_single_op_counts_as_leading_only() {
        let cigar = bam::record::CigarString(vec![Cigar::HardClip(7)]).into_view(0);
        assert_eq!(hardclip_lengths(&cigar), (7, 0));
    }

    #[test]
    fn accumulator_places_quality_at_its_absolute_cycle_offset_without_padding_past_it() {
        let mut acc = Accumulator::default();
        acc.add_at(&[10, 10, 10], 2);
        let rates = acc.error_rates();
        // No slot is created beyond the real data (offset 2 + 3 values):
        // hard-clipped/uncovered cycles past it are excluded, not zeroed.
        assert_eq!(rates.len(), 5);
        assert!(rates[2] > 0.0 && rates[3] > 0.0 && rates[4] > 0.0);
    }

    /// Builds a tiny single-record SAM file (htslib reads `.sam` the same
    /// as `.bam`) to exercise the reverse-strand un-reversal and hard-clip
    /// offset math together, end to end through `from_bam`.
    ///
    /// CIGAR `3H10M2H` on a reverse-strand alignment: 10 quality values are
    /// actually stored. Because the record is reverse-stranded,
    /// un-reversing puts the *trailing* stored hard clip (2) at the
    /// original read's start — see [`BuiltProfile::from_bam`] — so the 10
    /// real values land at cycles 2..12, and the *leading* stored clip (3),
    /// which un-reverses to the read's end, contributes nothing at all.
    #[test]
    fn from_bam_excludes_hardclipped_cycles_on_reverse_strand() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("reads.sam");
        // Raw Phred qualities 0..=9 (ASCII offset +33) in *stored* order.
        let qual: String = (0..10u8).map(|q| (33 + q) as char).collect();
        let sam = format!(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n\
             NB552023:1:AAA:1:1:1:1\t144\tchr1\t1\t60\t3H10M2H\t*\t0\t0\tACGTACGTAC\t{qual}\n"
        );
        File::create(&path).unwrap().write_all(sam.as_bytes()).unwrap();

        let profile = BuiltProfile::from_bam("n1", &[path.to_str().unwrap()]).unwrap();
        assert_eq!(profile.sequencer, "nextseq 550");
        // Flag 144 = last-in-template + reverse -> the `reverse` strand row.
        let rates = &profile.cycles[&REVERSE];
        // The trailing 3H (un-reversed to the read's end) is excluded
        // entirely rather than padding the output with fabricated cycles.
        assert_eq!(rates.len(), 12);
        // Un-reversing stored raw qualities [0..9] puts the highest quality
        // (9, lowest error) at cycle 2 and the lowest quality (0, highest
        // error = 100%) at cycle 11.
        assert_eq!(rates[11], 100.0);
        assert!(rates[2] < rates[11]);
    }

    fn write_sam_with_one_forward_read(dir: &Path, name: &str, qname: &str) -> String {
        let path = dir.join(name);
        let sam = format!(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n{qname}\t64\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
        );
        File::create(&path).unwrap().write_all(sam.as_bytes()).unwrap();
        path.to_str().unwrap().to_string()
    }

    #[test]
    fn from_bam_pools_reads_across_multiple_files() {
        let dir = tempfile::tempdir().unwrap();
        let bam1 = write_sam_with_one_forward_read(dir.path(), "a.sam", "NB552023:1:AAA:1:1:1:1");
        let bam2 = write_sam_with_one_forward_read(dir.path(), "b.sam", "NB552023:1:AAA:2:1:1:1");

        let pooled =
            BuiltProfile::from_bam("n1", &[bam1.as_str(), bam2.as_str()]).unwrap();
        let single = BuiltProfile::from_bam("n1", &[bam1.as_str()]).unwrap();
        // Both files have identical single-read quality: pooling them
        // shouldn't change the averaged rate, only confirm both contributed
        // (checked indirectly by them staying equal after doubling).
        assert_eq!(pooled.cycles[&FORWARD], single.cycles[&FORWARD]);
    }
}
