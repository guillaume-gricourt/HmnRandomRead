//! The sample's own variants around a breakpoint, called from its BAM and
//! phased onto one haplotype, so the fusion reads carry the alleles of the
//! homolog the fusion happened on instead of the bare reference.
//!
//! Without this, every heterozygous SNP/indel the sample carries next to the
//! breakpoint is absent from the fusion reads — a signature no real fusion
//! has, which lets a caller (or a phasing-aware filter) tell the injected
//! reads apart from the sample's own.
//!
//! The calling is deliberately minimal — a pileup genotyper with fixed
//! thresholds, not a variant caller: it only has to reproduce the clear-cut
//! germline alleles a real fusion read would carry, over a region of a few
//! hundred bases.

use std::collections::HashMap;
use std::error::Error;
use std::fmt;
use std::io;

use rust_htslib::bam::{self, pileup::Indel, Read as BamRead};

use crate::rng::RandomGenerator;

/// Minimum mapping quality of an alignment counted in the pileup.
const MIN_MAPQ: u8 = 20;
/// Minimum base quality of a base counted in the pileup.
const MIN_BASEQ: u8 = 13;
/// Minimum number of counted fragments (read pairs, however many of their
/// mates cover it) for a position to be genotyped.
const MIN_DEPTH: u32 = 10;
/// Alternate allele fraction from which a position is called heterozygous.
const MIN_HET_FRACTION: f64 = 0.2;
/// Alternate allele fraction from which a position is called homozygous.
const MIN_HOM_FRACTION: f64 = 0.8;

/// What one alignment shows at one reference position: the base it carries
/// there, plus the insertion that follows it or the number of reference
/// bases it deletes right after it — the pileup's own anchoring of indels.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Allele {
    pub base: u8,
    pub insertion: Vec<u8>,
    pub deletion: u32,
}

impl Allele {
    fn reference(base: u8) -> Self {
        Allele {
            base,
            insertion: Vec::new(),
            deletion: 0,
        }
    }
}

/// One variant applied to the fusion sequence: `allele` replaces the
/// reference base at 0-based `pos0`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Variant {
    pub chrom: String,
    pub pos0: u64,
    /// Reference bases the variant spans (the anchor base plus any deleted
    /// ones), for display only.
    pub reference: String,
    pub allele: Allele,
    pub homozygous: bool,
}

impl fmt::Display for Variant {
    /// VCF-like `chrom:pos:REF>ALT:hom|het`, 1-based.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut alt = vec![self.allele.base];
        alt.extend_from_slice(&self.allele.insertion);
        write!(
            f,
            "{}:{}:{}>{}:{}",
            self.chrom,
            self.pos0 + 1,
            self.reference,
            String::from_utf8_lossy(&alt),
            if self.homozygous { "hom" } else { "het" }
        )
    }
}

/// A heterozygous candidate, with which allele each read pair carries (true
/// for the alternate one), for phasing.
struct Heterozygous {
    variant: Variant,
    carriers: HashMap<Vec<u8>, bool>,
}

fn is_counted(record: &bam::Record) -> bool {
    !record.is_unmapped()
        && !record.is_secondary()
        && !record.is_supplementary()
        && !record.is_quality_check_failed()
        && !record.is_duplicate()
        && record.mapq() >= MIN_MAPQ
}

/// Call the sample's variants over `[begin0, end0]` (0-based, inclusive) of
/// `chrom` from `bam_path`, and keep the ones on a single haplotype.
///
/// `reference` is the reference sequence of exactly that region. Homozygous
/// variants are always kept. Heterozygous ones are phased greedily from the
/// reads: the first one's haplotype is drawn at random, each next one takes
/// the allele that the read pairs linking it to the already phased sites
/// vote for, and falls back to a random draw when no read pair links it.
/// That keeps nearby alleles in the sample's own phase, which is all a read
/// pair can ever show.
pub fn sample_haplotype(
    bam_path: &str,
    chrom: &str,
    begin0: u64,
    end0: u64,
    reference: &str,
    rng: &mut RandomGenerator,
) -> Result<Vec<Variant>, Box<dyn Error>> {
    let reference = reference.as_bytes();
    if reference.len() as u64 != end0 - begin0 + 1 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "reference sequence does not match the region to call variants over",
        )
        .into());
    }

    let mut reader = bam::IndexedReader::from_path(bam_path)?;
    let tid = reader.header().tid(chrom.as_bytes()).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("chromosome '{chrom}' not found in the header of '{bam_path}'"),
        )
    })?;
    reader.fetch((tid, begin0 as i64, end0 as i64 + 1))?;

    let mut homozygous = Vec::new();
    let mut heterozygous = Vec::new();
    let mut pileups = reader.pileup();
    pileups.set_max_depth(i32::MAX as u32);
    for pileup in pileups {
        let pileup = pileup?;
        let pos0 = pileup.pos() as u64;
        if pileup.tid() != tid || pos0 < begin0 || pos0 > end0 {
            continue;
        }
        let reference_base = reference[(pos0 - begin0) as usize].to_ascii_uppercase();
        if reference_base == b'N' {
            continue;
        }

        // One observation per fragment, not per alignment: where the two
        // mates of a pair overlap they read the same molecule, and counting
        // both would weigh that molecule twice (and make a pair whose mates
        // disagree look like a heterozygous site). The mates' alleles are
        // merged: agreeing, one vote; disagreeing, the higher base quality
        // wins, and a tie drops the fragment at this position.
        let mut fragments: HashMap<Vec<u8>, Option<(Allele, u8)>> = HashMap::new();
        for alignment in pileup.alignments() {
            let record = alignment.record();
            if !is_counted(&record) || alignment.is_del() || alignment.is_refskip() {
                continue;
            }
            let Some(qpos) = alignment.qpos() else {
                continue;
            };
            let quality = record.qual()[qpos];
            if quality < MIN_BASEQ {
                continue;
            }
            let seq = record.seq();
            let base = seq[qpos].to_ascii_uppercase();
            let allele = match alignment.indel() {
                Indel::Ins(len) => Allele {
                    base,
                    insertion: (qpos + 1..qpos + 1 + len as usize)
                        .map(|i| seq[i].to_ascii_uppercase())
                        .collect(),
                    deletion: 0,
                },
                Indel::Del(len) => Allele {
                    base,
                    insertion: Vec::new(),
                    deletion: len,
                },
                Indel::None => Allele::reference(base),
            };
            fragments
                .entry(record.qname().to_vec())
                .and_modify(|seen| {
                    *seen = match seen.take() {
                        Some((kept, kept_quality)) if kept == allele => {
                            Some((kept, kept_quality.max(quality)))
                        }
                        Some((kept, kept_quality)) => match kept_quality.cmp(&quality) {
                            std::cmp::Ordering::Greater => Some((kept, kept_quality)),
                            std::cmp::Ordering::Less => Some((allele.clone(), quality)),
                            std::cmp::Ordering::Equal => None,
                        },
                        // Already dropped: a third alignment (e.g. a
                        // supplementary one, filtered above) can't settle it.
                        None => None,
                    };
                })
                .or_insert_with(|| Some((allele.clone(), quality)));
        }
        let observations: Vec<(Vec<u8>, Allele)> = fragments
            .into_iter()
            .filter_map(|(qname, seen)| seen.map(|(allele, _)| (qname, allele)))
            .collect();
        let mut counts: HashMap<Allele, u32> = HashMap::new();
        for (_, allele) in &observations {
            *counts.entry(allele.clone()).or_default() += 1;
        }

        let depth: u32 = counts.values().sum();
        if depth < MIN_DEPTH {
            continue;
        }
        let reference_allele = Allele::reference(reference_base);
        // Ties broken on the allele itself so the call doesn't depend on
        // HashMap iteration order.
        let Some((alternate, count)) = counts
            .iter()
            .filter(|(allele, _)| **allele != reference_allele)
            .max_by(|a, b| a.1.cmp(b.1).then_with(|| b.0.base.cmp(&a.0.base)))
            .map(|(allele, count)| (allele.clone(), *count))
        else {
            continue;
        };
        let fraction = count as f64 / depth as f64;
        if fraction < MIN_HET_FRACTION {
            continue;
        }

        let span_end = (pos0 + alternate.deletion as u64).min(end0);
        let variant = Variant {
            chrom: chrom.to_string(),
            pos0,
            reference: String::from_utf8_lossy(
                &reference[(pos0 - begin0) as usize..=(span_end - begin0) as usize],
            )
            .to_ascii_uppercase(),
            allele: alternate.clone(),
            homozygous: fraction >= MIN_HOM_FRACTION,
        };
        if variant.homozygous {
            homozygous.push(variant);
            continue;
        }

        // Mates already merged above: one allele per fragment.
        let carriers = observations
            .into_iter()
            .filter_map(|(qname, allele)| {
                if allele == alternate {
                    Some((qname, true))
                } else if allele == reference_allele {
                    Some((qname, false))
                } else {
                    None
                }
            })
            .collect();
        heterozygous.push(Heterozygous { variant, carriers });
    }

    let mut kept = homozygous;
    let mut phased: Vec<bool> = Vec::with_capacity(heterozygous.len());
    for k in 0..heterozygous.len() {
        let mut vote = 0i64;
        for j in 0..k {
            for (qname, &carries_here) in &heterozygous[k].carriers {
                if let Some(&carries_there) = heterozygous[j].carriers.get(qname) {
                    // On the kept haplotype at site j, the pair's allele here
                    // is the kept one; off it, the other one.
                    let on_kept = carries_there == phased[j];
                    vote += if carries_here == on_kept { 1 } else { -1 };
                }
            }
        }
        let alternate_kept = match vote.cmp(&0) {
            std::cmp::Ordering::Greater => true,
            std::cmp::Ordering::Less => false,
            std::cmp::Ordering::Equal => rng.unit() < 0.5,
        };
        phased.push(alternate_kept);
        if alternate_kept {
            kept.push(heterozygous[k].variant.clone());
        }
    }
    kept.sort_by_key(|v| v.pos0);
    Ok(kept)
}

/// Apply `variants` to `sequence`, the reference over the region starting at
/// 0-based `begin0`. Variants outside the region are ignored, and a
/// deletion running past its end is cut there.
pub fn apply_variants(sequence: &str, begin0: u64, variants: &[Variant]) -> String {
    let by_position: HashMap<u64, &Variant> = variants.iter().map(|v| (v.pos0, v)).collect();
    let mut out = String::with_capacity(sequence.len());
    let mut skip = 0u32;
    for (i, base) in sequence.chars().enumerate() {
        if skip > 0 {
            skip -= 1;
            continue;
        }
        match by_position.get(&(begin0 + i as u64)) {
            Some(variant) => {
                out.push(variant.allele.base as char);
                out.extend(variant.allele.insertion.iter().map(|&b| b as char));
                skip = variant.allele.deletion;
            }
            None => out.push(base),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn variant(pos0: u64, base: u8, insertion: &[u8], deletion: u32) -> Variant {
        Variant {
            chrom: "chr1".into(),
            pos0,
            reference: "N".into(),
            allele: Allele {
                base,
                insertion: insertion.to_vec(),
                deletion,
            },
            homozygous: true,
        }
    }

    #[test]
    fn apply_variants_handles_snv_insertion_and_deletion() {
        // Region starts at 0-based 100: "ACGTACGT" covers 100..=107.
        let seq = "ACGTACGT";
        let variants = [
            variant(101, b'T', b"", 0),   // C>T
            variant(103, b'T', b"GG", 0), // T>TGG
            variant(105, b'C', b"", 2),   // CGT>C
        ];
        assert_eq!(apply_variants(seq, 100, &variants), "ATGTGGAC");
    }

    #[test]
    fn apply_variants_ignores_variants_outside_the_region_and_cuts_deletions() {
        let seq = "ACGT";
        let variants = [variant(10, b'G', b"", 0), variant(2, b'G', b"", 5)];
        assert_eq!(apply_variants(seq, 0, &variants), "ACG");
    }

    fn write_bam(dir: &std::path::Path, sam: &str) -> String {
        let sam_path = dir.join("reads.sam");
        std::fs::write(&sam_path, sam).unwrap();
        let bam_path = dir.join("reads.bam");
        {
            let template = bam::Reader::from_path(&sam_path).unwrap();
            let header = bam::Header::from_template(template.header());
            let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
            let mut reader = bam::Reader::from_path(&sam_path).unwrap();
            let mut record = bam::Record::new();
            while let Some(result) = reader.read(&mut record) {
                result.unwrap();
                writer.write(&record).unwrap();
            }
        }
        bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();
        bam_path.to_str().unwrap().to_string()
    }

    #[test]
    fn keeps_a_homozygous_snv_and_phases_linked_heterozygous_ones_together() {
        // Reference "AAAAAAAAAAAAAAAAAAAA" (20bp). 20 reads over [0,20):
        // every one carries G at 0-based 2 (homozygous); half carry C at 5
        // *and* T at 15, the other half neither (two linked hets).
        let reference = "A".repeat(20);
        let mut sam = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:20\n");
        for i in 0..20 {
            let mut seq = reference.clone().into_bytes();
            seq[2] = b'G';
            if i % 2 == 0 {
                seq[5] = b'C';
                seq[15] = b'T';
            }
            sam.push_str(&format!(
                "r{i}\t0\tchr1\t1\t60\t20M\t*\t0\t0\t{}\t{}\n",
                String::from_utf8(seq).unwrap(),
                "I".repeat(20)
            ));
        }
        let dir = tempfile::tempdir().unwrap();
        let bam = write_bam(dir.path(), &sam);

        for seed in 0..8 {
            let mut rng = RandomGenerator::new(seed);
            let variants = sample_haplotype(&bam, "chr1", 0, 19, &reference, &mut rng).unwrap();
            let applied = apply_variants(&reference, 0, &variants);
            let bytes = applied.as_bytes();
            assert_eq!(bytes[2], b'G', "homozygous SNV always applied");
            // Either both hets or neither — never one without the other.
            assert_eq!(
                bytes[5] == b'C',
                bytes[15] == b'T',
                "seed {seed}: {applied}"
            );
        }
    }

    #[test]
    fn counts_overlapping_mates_once() {
        // 12 fragments whose two mates both cover [0,20). Mate 1 always
        // reads the reference; mate 2 reads a C at 0-based 5 with a lower
        // base quality ('5' = 20 < 'I' = 40): counted per alignment that
        // would be a 50% het, per fragment it is no variant at all. At
        // 0-based 10 the mates disagree at equal quality: fragment dropped.
        let reference = "A".repeat(20);
        let mut sam = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:20\n");
        let mut qual2 = "I".repeat(20).into_bytes();
        qual2[5] = b'5';
        let qual2 = String::from_utf8(qual2).unwrap();
        let mut seq2 = reference.clone().into_bytes();
        seq2[5] = b'C';
        seq2[10] = b'G';
        let seq2 = String::from_utf8(seq2).unwrap();
        for i in 0..12 {
            sam.push_str(&format!(
                "f{i}\t99\tchr1\t1\t60\t20M\t=\t1\t20\t{reference}\t{}\n",
                "I".repeat(20)
            ));
            sam.push_str(&format!(
                "f{i}\t147\tchr1\t1\t60\t20M\t=\t1\t-20\t{seq2}\t{qual2}\n"
            ));
        }
        let dir = tempfile::tempdir().unwrap();
        let bam = write_bam(dir.path(), &sam);
        let mut rng = RandomGenerator::new(0);
        assert!(sample_haplotype(&bam, "chr1", 0, 19, &reference, &mut rng)
            .unwrap()
            .is_empty());
    }

    #[test]
    fn ignores_positions_below_the_minimum_depth() {
        let reference = "A".repeat(20);
        let mut sam = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:20\n");
        for i in 0..(MIN_DEPTH - 1) {
            sam.push_str(&format!(
                "r{i}\t0\tchr1\t1\t60\t20M\t*\t0\t0\tAAGAAAAAAAAAAAAAAAAA\t{}\n",
                "I".repeat(20)
            ));
        }
        let dir = tempfile::tempdir().unwrap();
        let bam = write_bam(dir.path(), &sam);
        let mut rng = RandomGenerator::new(0);
        assert!(sample_haplotype(&bam, "chr1", 0, 19, &reference, &mut rng)
            .unwrap()
            .is_empty());
    }

    #[test]
    fn calls_an_insertion_and_a_deletion() {
        // 12 reads: 1-base insertion (G) after 0-based 4, 2-base deletion
        // of 0-based 11..=12 (anchored on the G at 10).
        let reference = "ACGTACGTACGTACGTACGT";
        let read = "ACGTAGCGTACGCGTACGT"; // 5M 1I 6M 2D 7M, 19 query bases
        let mut sam = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:20\n");
        for i in 0..12 {
            sam.push_str(&format!(
                "r{i}\t0\tchr1\t1\t60\t5M1I6M2D7M\t*\t0\t0\t{read}\t{}\n",
                "I".repeat(19)
            ));
        }
        let dir = tempfile::tempdir().unwrap();
        let bam = write_bam(dir.path(), &sam);
        let mut rng = RandomGenerator::new(0);
        let variants = sample_haplotype(&bam, "chr1", 0, 19, reference, &mut rng).unwrap();
        assert_eq!(variants.len(), 2, "{variants:?}");
        assert_eq!(variants[0].to_string(), "chr1:5:A>AG:hom");
        assert_eq!(variants[1].to_string(), "chr1:11:GTA>G:hom");
        assert_eq!(
            apply_variants(reference, 0, &variants),
            "ACGTAGCGTACGCGTACGT"
        );
    }
}
