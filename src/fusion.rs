//! Builds a chimeric junction sequence from two breakpoints and draws
//! fusion-supporting read pairs from it — the core of `fusion-simulate`,
//! which writes a pool of such pairs plus the truth describing them, and of
//! `fusion-spike`, which dilutes that pool into a real sample at a chosen
//! allelic fraction.

use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::io;

use rust_htslib::bam::{
    self,
    record::{Cigar, CigarString},
    Read as BamRead,
};

use crate::diversity::ProfileDiversity;
use crate::fastq::{build_read_pair, FastqRecord};
use crate::haplotype::{apply_variants, Variant};
use crate::io::FastaIndexedReader;
use crate::profile_sequencer::ProfileSequencer;
use crate::rng::RandomGenerator;
use crate::sequence::{ReverseComplement, Sequence};

const MAX_FRAGMENT_ATTEMPTS: usize = 100;

fn parse_err(msg: impl Into<String>) -> Box<dyn Error> {
    io::Error::new(io::ErrorKind::InvalidData, msg.into()).into()
}

/// Strand of a fusion partner's gene on the reference, which decides which
/// side of its breakpoint it contributes and in which orientation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

impl Strand {
    /// Parse `+` or `-`.
    pub fn parse(spec: &str) -> Result<Self, Box<dyn Error>> {
        match spec {
            "+" => Ok(Strand::Forward),
            "-" => Ok(Strand::Reverse),
            _ => Err(parse_err(format!("malformed strand '{spec}', expected '+' or '-'"))),
        }
    }

    pub fn symbol(self) -> char {
        match self {
            Strand::Forward => '+',
            Strand::Reverse => '-',
        }
    }
}

/// A fusion partner's breakpoint: 1-based `pos` is the last reference base
/// before the break, which always falls between `pos` and `pos + 1` on the
/// reference, whatever the strand. `strand` is the strand the partner is
/// taken on (drawn per read pair by `fusion-simulate`):
/// a `+` gene keeps the bases up to `pos` as its 5' part and from `pos + 1`
/// on as its 3' part; a `-` gene the reverse complement of the bases from
/// `pos + 1` on as its 5' part and of the bases up to `pos` as its 3' part.
/// Two `+` partners (or two `-`) give a translocation; one of each, an
/// inversion-type junction such as EML4-ALK.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Breakpoint {
    pub chrom: String,
    pub pos: u64,
    pub strand: Strand,
}

impl Breakpoint {
    /// Parse a `chrom:pos` spec (1-based `pos`), on the `+` strand.
    pub fn parse(spec: &str) -> Result<Self, Box<dyn Error>> {
        let (chrom, pos) = spec
            .rsplit_once(':')
            .ok_or_else(|| parse_err(format!("malformed breakpoint '{spec}', expected 'chrom:pos'")))?;
        if chrom.is_empty() {
            return Err(parse_err(format!("malformed breakpoint '{spec}', chrom is empty")));
        }
        let pos: u64 = pos
            .parse()
            .map_err(|_| parse_err(format!("malformed breakpoint '{spec}', pos isn't a positive integer")))?;
        if pos == 0 {
            return Err(parse_err(format!("malformed breakpoint '{spec}', pos must be >= 1")));
        }
        Ok(Breakpoint { chrom: chrom.to_string(), pos, strand: Strand::Forward })
    }

    /// `chrom:pos`, as parsed by [`Breakpoint::parse`] (the strand is not
    /// part of it).
    pub fn position(&self) -> String {
        format!("{}:{}", self.chrom, self.pos)
    }
}

/// Prefix carried by every produced read's QNAME, so a spiked BAM can be
/// reduced to the synthetic reads with a single name filter — whatever
/// position the aligner ended up placing them at.
pub const QNAME_PREFIX: &str = "hmnrr_fusion";

/// The QNAME shared by both mates of fusion read pair `number`:
/// `hmnrr_fusion_<left_chrom>-<left_pos>_<right_chrom>-<right_pos>_<number>`.
///
/// Unlike `simulate`'s header, which carries its metadata in the FASTQ
/// *comment* (everything after the first space), everything here lives in
/// the name itself: aligners drop the comment unless explicitly told to keep
/// it (`bwa mem -C`), so a comment-borne tag never reaches the BAM and the
/// reads become impossible to tell apart from the sample's own once aligned.
///
/// Both mates get the byte-identical name, with no `/1`,`/2` suffix, exactly
/// as Illumina's own FASTQ does. That keeps every aligner pairing them
/// correctly and makes the name in the BAM identical to the one in the truth
/// TSV, with no trimming rule to account for — BWA strips a trailing
/// `/<digit>` from a QNAME, other aligners don't.
pub fn build_qname(left_bp: &Breakpoint, right_bp: &Breakpoint, number: u64) -> String {
    format!(
        "{QNAME_PREFIX}_{}-{}_{}-{}_{number:010}",
        left_bp.chrom, left_bp.pos, right_bp.chrom, right_bp.pos
    )
}

/// Human-readable junction identifier,
/// `left_chrom:left_pos(strand)>right_chrom:right_pos(strand)`.
pub fn junction_label(left_bp: &Breakpoint, right_bp: &Breakpoint) -> String {
    format!(
        "{}:{}({})>{}:{}({})",
        left_bp.chrom,
        left_bp.pos,
        left_bp.strand.symbol(),
        right_bp.chrom,
        right_bp.pos,
        right_bp.strand.symbol()
    )
}

/// The 5' partner's `chrom:pos` in a [`junction_label`], strand left out.
pub fn five_prime_position(label: &str) -> &str {
    let left = label.split_once('>').map_or(label, |(left, _)| left);
    left.rsplit_once('(').map_or(left, |(position, _)| position)
}

/// The 5' and 3' partners, with their strands, a [`junction_label`] was
/// built from.
pub fn parse_junction_label(label: &str) -> Result<(Breakpoint, Breakpoint), Box<dyn Error>> {
    let partner = |spec: &str| -> Result<Breakpoint, Box<dyn Error>> {
        let (position, strand) = spec
            .strip_suffix(')')
            .and_then(|s| s.rsplit_once('('))
            .ok_or_else(|| parse_err(format!("malformed junction partner '{spec}' in '{label}'")))?;
        let mut bp = Breakpoint::parse(position)?;
        bp.strand = Strand::parse(strand)?;
        Ok(bp)
    };
    let (five, three) = label
        .split_once('>')
        .ok_or_else(|| parse_err(format!("malformed junction label '{label}'")))?;
    Ok((partner(five)?, partner(three)?))
}

/// The two reads sequenced from a fragment of `junction` carrying
/// `junction_offset` bases of its 5' partner and `fragment_len` bases in
/// all, each `read_len` long: the head read on the fragment's strand, the
/// tail read reverse-complemented, as [`crate::fastq::build_read_pair`]
/// makes them. `None` if the fragment doesn't fit in `junction`.
pub fn fragment_reads(
    junction: &str,
    junction_index: usize,
    junction_offset: usize,
    fragment_len: usize,
    read_len: usize,
) -> Option<(String, String)> {
    let start = junction_index.checked_sub(junction_offset)?;
    let fragment = junction.get(start..start + fragment_len)?;
    let take = read_len.min(fragment_len);
    Some((
        fragment[..take].to_string(),
        fragment[fragment_len - take..].to_string().reverse_complement(),
    ))
}

/// `read`, a pool read sequenced from the reference allele `clean`, made
/// to carry the sample's alleles `carried` instead. Wherever the pool read
/// differs from `clean` — a sequencing error or diversity mutation drawn
/// when the pool was built — its own base is kept, so only the bases the
/// sample's variants change are rewritten. Case follows the pool read.
pub fn carry_variants(read: &str, clean: &str, carried: &str) -> String {
    read.bytes()
        .zip(clean.bytes().chain(std::iter::repeat(b'N')))
        .zip(carried.bytes().chain(std::iter::repeat(b'N')))
        .map(|((r, c), v)| {
            if !r.eq_ignore_ascii_case(&c) {
                r
            } else if r.is_ascii_lowercase() {
                v.to_ascii_lowercase()
            } else {
                v.to_ascii_uppercase()
            }
        })
        .map(char::from)
        .collect()
}

/// Number of positions where `a` and `b` differ, case aside.
pub fn mismatches(a: &str, b: &str) -> usize {
    a.bytes().zip(b.bytes()).filter(|(x, y)| !x.eq_ignore_ascii_case(y)).count()
        + a.len().abs_diff(b.len())
}

/// How a fusion read pair shows the junction once sequenced — decided by
/// where the junction falls in its fragment, not chosen: every fragment
/// crossing the junction is a fusion molecule, and a real library holds all
/// three kinds.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Support {
    /// A read crosses the junction with at least `min_anchor` bases on each
    /// side of it: a split (chimeric/supplementary) alignment.
    SplitRead,
    /// A read crosses the junction, but with fewer than `min_anchor` bases
    /// on one side: typically aligned whole to one partner, soft-clipped.
    ShortAnchor,
    /// The junction falls between the two reads: each mate aligns whole to
    /// a different partner, a discordant pair.
    DiscordantPair,
}

impl Support {
    pub fn as_str(self) -> &'static str {
        match self {
            Support::SplitRead => "split_read",
            Support::ShortAnchor => "short_anchor",
            Support::DiscordantPair => "discordant_pair",
        }
    }

    pub fn parse(spec: &str) -> Result<Self, Box<dyn Error>> {
        match spec {
            "split_read" => Ok(Support::SplitRead),
            "short_anchor" => Ok(Support::ShortAnchor),
            "discordant_pair" => Ok(Support::DiscordantPair),
            _ => Err(parse_err(format!("unknown support type '{spec}'"))),
        }
    }
}

/// One produced fusion read pair, as written to the truth TSV: the ground
/// truth of what was injected, so a quantification never has to be
/// reconstructed from the very alignment it is meant to evaluate.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TruthRow {
    /// QNAME shared by both mates, byte-identical to the one in the BAM.
    pub qname: String,
    /// `left_chrom:left_pos>right_chrom:right_pos`.
    pub junction: String,
    /// Fragment bounds within the junction sequence.
    pub fragment_start: usize,
    pub fragment_stop: usize,
    /// Number of bases the fragment carries from the left partner, i.e. where
    /// the junction falls inside it — the clip point an aligner should find.
    pub junction_offset: usize,
    /// How the pair shows the junction once sequenced.
    pub support: Support,
}

enum Side {
    Left,
    Right,
}

/// Fetch up to `flank_len` bases immediately to one side of `bp`, on the
/// forward strand, with the 0-based position of the first one (clamped at
/// the contig's boundary, so a breakpoint near the start/end of a contig
/// yields a shorter-than-requested flank rather than an error).
fn fetch_flank(
    faidx: &FastaIndexedReader,
    bp: &Breakpoint,
    flank_len: u64,
    side: Side,
) -> io::Result<(u64, String)> {
    let contig_len = faidx.seq_len(&bp.chrom)?;
    if bp.pos > contig_len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("breakpoint position {} exceeds '{}' length ({contig_len})", bp.pos, bp.chrom),
        ));
    }
    match side {
        Side::Left => {
            let end0 = bp.pos - 1; // 0-based index of the base at 1-based `pos`
            let begin0 = end0.saturating_sub(flank_len.saturating_sub(1));
            Ok((begin0, faidx.fetch(&bp.chrom, begin0, end0)?))
        }
        Side::Right => {
            let begin0 = bp.pos; // 0-based index right after `pos`
            if begin0 >= contig_len {
                return Ok((begin0, String::new()));
            }
            let end0 = (begin0 + flank_len - 1).min(contig_len - 1);
            Ok((begin0, faidx.fetch(&bp.chrom, begin0, end0)?))
        }
    }
}

/// The 0-based, inclusive region around `bp` that its two flanks of
/// `flank_len` bases can cover — what the sample's variants have to be
/// called over for either junction to carry them.
pub fn flank_region(
    faidx: &FastaIndexedReader,
    bp: &Breakpoint,
    flank_len: u64,
) -> io::Result<(u64, u64)> {
    let contig_len = faidx.seq_len(&bp.chrom)?;
    if bp.pos > contig_len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("breakpoint position {} exceeds '{}' length ({contig_len})", bp.pos, bp.chrom),
        ));
    }
    let begin0 = bp.pos.saturating_sub(flank_len);
    let end0 = (bp.pos + flank_len - 1).min(contig_len - 1);
    Ok((begin0, end0))
}

/// Unlike `simulate`, which only ever draws fragments from non-N
/// [`crate::scaffold::Scaffolds`] runs, a breakpoint here is an exact
/// coordinate the user asked for — there's no other position to silently
/// fall back to. Rather than let an assembly gap near the breakpoint
/// produce synthetic reads full of literal `N` bases (which no real
/// sequencer emits, and which downstream aligners would mishandle), flag it
/// as a clear error instead.
fn contains_n(seq: &str) -> bool {
    seq.bytes().any(|b| b.eq_ignore_ascii_case(&b'N'))
}

/// One side of a junction: a breakpoint, the reference it lies on, and the
/// sample's variants to carry over its flank (empty for the bare
/// reference).
pub struct JunctionPartner<'a> {
    pub faidx: &'a FastaIndexedReader,
    pub bp: &'a Breakpoint,
    pub variants: &'a [Variant],
}

enum Role {
    FivePrime,
    ThreePrime,
}

/// The flank `partner` contributes to a junction in `role`, oriented 5'->3'
/// along the fusion: which side of the breakpoint it is and whether it is
/// reverse-complemented both follow from the partner's strand (see
/// [`Breakpoint`]). The sample's variants are applied on the forward strand,
/// before any reverse complement.
fn partner_flank(partner: &JunctionPartner, role: Role, flank_len: u64) -> io::Result<String> {
    let bp = partner.bp;
    let side = match (&role, bp.strand) {
        (Role::FivePrime, Strand::Forward) | (Role::ThreePrime, Strand::Reverse) => Side::Left,
        (Role::FivePrime, Strand::Reverse) | (Role::ThreePrime, Strand::Forward) => Side::Right,
    };
    let where_ = match side {
        Side::Left => "up to",
        Side::Right => "after",
    };
    let (begin0, flank) = fetch_flank(partner.faidx, bp, flank_len, side)?;
    if flank.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "no reference sequence {where_} breakpoint {}:{}: a {} strand partner needs \
                 bases on that side of it",
                bp.chrom,
                bp.pos,
                bp.strand.symbol()
            ),
        ));
    }
    if contains_n(&flank) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "the {flank_len}bp flank {where_} breakpoint {}:{} contains N bases \
                 (assembly gap?); pick a breakpoint away from reference gaps",
                bp.chrom, bp.pos
            ),
        ));
    }
    let flank = apply_variants(&flank, begin0, partner.variants);
    Ok(match bp.strand {
        Strand::Forward => flank,
        Strand::Reverse => flank.reverse_complement(),
    })
}

/// Build a chimeric junction sequence: `five_prime`'s 5' part followed by
/// `three_prime`'s 3' part, each taken from the side of its breakpoint and
/// in the orientation its strand calls for. Returns the concatenated
/// sequence and the 0-based index of its first base belonging to the 3'
/// partner (i.e. the junction point), which may be less than `flank_len` if
/// the 5' flank was clamped near a contig boundary or shortened by a
/// deletion.
pub fn build_junction(
    five_prime: &JunctionPartner,
    three_prime: &JunctionPartner,
    flank_len: u64,
) -> io::Result<(String, usize)> {
    let left = partner_flank(five_prime, Role::FivePrime, flank_len)?;
    let right = partner_flank(three_prime, Role::ThreePrime, flank_len)?;
    let junction_index = left.len();
    Ok((left + &right, junction_index))
}

/// Whether an alignment counts toward [`depth_at`]'s depth — matches
/// `samtools depth`'s default exclude filter (unmapped, secondary, QC-fail,
/// duplicate), since htslib's pileup engine itself does not filter anything
/// on its own: without this, `depth_at` would count every alignment
/// (including secondary/duplicate) and systematically overstate real
/// coverage relative to what a user comparing against `samtools depth`
/// would expect.
fn is_countable(record: &bam::Record) -> bool {
    !record.is_unmapped()
        && !record.is_secondary()
        && !record.is_quality_check_failed()
        && !record.is_duplicate()
}

/// Pileup depth at a single 1-based position of an indexed BAM (0 if the
/// position has no coverage), counting only alignments [`is_countable`].
pub fn depth_at(bam_path: &str, chrom: &str, pos_1based: u64) -> Result<u32, Box<dyn Error>> {
    let mut reader = bam::IndexedReader::from_path(bam_path)?;
    let tid = reader
        .header()
        .tid(chrom.as_bytes())
        .ok_or_else(|| parse_err(format!("chromosome '{chrom}' not found in BAM header")))?;
    let pos0 = pos_1based - 1;
    reader.fetch((tid, pos0 as i64, pos0 as i64 + 1))?;

    for p in reader.pileup() {
        let p = p?;
        if p.tid() == tid && p.pos() as u64 == pos0 {
            let depth = p.alignments().filter(|a| is_countable(&a.record())).count() as u32;
            return Ok(depth);
        }
    }
    Ok(0)
}

/// The reference fragments crossing a breakpoint, as
/// [`reference_fragments_at`] collects them.
#[derive(Debug, Default)]
pub struct ReferenceFragments {
    /// QNAME of each molecule's representative: the pair not flagged as a
    /// duplicate.
    pub fragments: HashSet<String>,
    /// QNAMEs of the pairs flagged as PCR/optical duplicates of each
    /// representative, found by the duplicate key a marker groups them
    /// under (see [`duplicate_key`]).
    pub duplicates: HashMap<String, Vec<String>>,
    /// Pairs flagged as duplicates crossing the break, whether or not
    /// their representative was found. Zero at any real depth usually
    /// means the BAM was never duplicate-marked.
    pub duplicate_pairs: u64,
}

impl ReferenceFragments {
    /// Number of molecules (representatives).
    pub fn len(&self) -> usize {
        self.fragments.len()
    }

    pub fn is_empty(&self) -> bool {
        self.fragments.is_empty()
    }

    /// `qname` and every duplicate of it: all the read pairs that carry
    /// the same molecule, and have to go together for it to be gone once
    /// the output is realigned and duplicate-marked again.
    pub fn with_duplicates<'a>(&'a self, qname: &'a str) -> impl Iterator<Item = &'a str> {
        std::iter::once(qname).chain(
            self.duplicates
                .get(qname)
                .into_iter()
                .flatten()
                .map(String::as_str),
        )
    }
}

/// A read pair's 5' ends, strands and library — what Picard MarkDuplicates
/// and samtools markdup group duplicates by: two pairs from the same library
/// whose mates start at the same unclipped 5' positions on the same strands
/// are copies of one molecule. Taken from the pair's leftmost mate; the
/// other mate's unclipped 5' end comes from its `MC` tag, or its aligned
/// start when the aligner didn't write one.
type DuplicateKey = (String, i64, bool, i32, i64, bool);

/// Unclipped 5' end (0-based) of an alignment starting at `pos` with
/// `cigar`: soft/hard clips are added back, since two copies of a molecule
/// may be clipped differently by base-calling errors.
fn unclipped_five_prime(pos: i64, reverse: bool, cigar: &[Cigar]) -> i64 {
    let clip = |c: &Cigar| match c {
        Cigar::SoftClip(n) | Cigar::HardClip(n) => Some(*n as i64),
        _ => None,
    };
    if !reverse {
        pos - cigar.iter().map_while(clip).sum::<i64>()
    } else {
        let reference_len: i64 = cigar
            .iter()
            .map(|c| match c {
                Cigar::Match(n) | Cigar::Del(n) | Cigar::RefSkip(n) | Cigar::Equal(n)
                | Cigar::Diff(n) => *n as i64,
                _ => 0,
            })
            .sum();
        pos + reference_len - 1 + cigar.iter().rev().map_while(clip).sum::<i64>()
    }
}

fn duplicate_key(record: &bam::Record, libraries: &HashMap<String, String>) -> DuplicateKey {
    let library = match record.aux(b"RG") {
        Ok(bam::record::Aux::String(rg)) => libraries.get(rg).cloned().unwrap_or_else(|| rg.to_string()),
        _ => String::new(),
    };
    let own = unclipped_five_prime(record.pos(), record.is_reverse(), &record.cigar());
    let mate_cigar = match record.aux(b"MC") {
        Ok(bam::record::Aux::String(mc)) => CigarString::try_from(mc).ok(),
        _ => None,
    };
    let mate = match mate_cigar {
        Some(cigar) => unclipped_five_prime(record.mpos(), record.is_mate_reverse(), &cigar.0),
        None => record.mpos(),
    };
    (library, own, record.is_reverse(), record.mtid(), mate, record.is_mate_reverse())
}

/// Read group ID -> library (`LB`), for [`duplicate_key`]: duplicates are
/// only ever marked within a library, whatever read groups it spans.
fn read_group_libraries(header: &bam::HeaderView) -> HashMap<String, String> {
    bam::Header::from_template(header)
        .to_hashmap()
        .get("RG")
        .into_iter()
        .flatten()
        .filter_map(|rg| Some((rg.get("ID")?.clone(), rg.get("LB")?.clone())))
        .collect()
}

/// The reference *fragments* (template molecules, counted once however
/// many of their mates overlap) that cross the break between 1-based
/// `pos_1based` and `pos_1based + 1` — whether one of their reads covers it
/// or it falls in the unsequenced insert between the two — with the pairs
/// flagged as duplicates of each.
///
/// This — not [`depth_at`] — is the denominator an allelic fraction is
/// defined against, and it counts the same class of molecule the fusion
/// pool is drawn from: every fragment crossing the junction, split read,
/// short anchor or discordant pair alike. Counting only reads covering the
/// position would leave out the reference counterpart of the discordant
/// pairs, and pileup reads would count a molecule twice whenever its mates
/// overlap. Duplicate-flagged pairs are not molecules of their own and are
/// left out of the count, which therefore expects a duplicate-marked BAM;
/// they are still collected, keyed to their representative, so removing a
/// molecule removes all its copies. Both [`depth_at`] and this are reported
/// by `fusion-spike` so their ratio stays visible.
///
/// A fragment's extent is read off its leftmost mate's position and TLEN.
/// Those mates can start up to one insert size before the break, so the
/// reads covering the break are scanned first for the largest TLEN among
/// them, which bounds the window scanned for leftmost mates; the window is
/// widened for as long as the reads it holds show a longer TLEN still.
pub fn reference_fragments_at(
    bam_path: &str,
    chrom: &str,
    pos_1based: u64,
) -> Result<ReferenceFragments, Box<dyn Error>> {
    let mut reader = bam::IndexedReader::from_path(bam_path)?;
    let tid = reader
        .header()
        .tid(chrom.as_bytes())
        .ok_or_else(|| parse_err(format!("chromosome '{chrom}' not found in BAM header")))?;
    let libraries = read_group_libraries(reader.header());
    let last0 = pos_1based - 1; // last base before the break, 0-based
    reader.fetch((tid, last0 as i64, last0 as i64 + 2))?;
    let mut window = 0i64;
    for record in reader.records() {
        let record = record?;
        if is_reference_template(&record) {
            window = window.max(record.insert_size().abs());
        }
    }

    loop {
        reader.fetch((tid, (last0 as i64 - window).max(0), last0 as i64 + 1))?;
        // Several unflagged pairs can share a key (a BAM never marked, or
        // marked with other criteria): each is a molecule of its own.
        let mut representatives: HashMap<DuplicateKey, Vec<String>> = HashMap::new();
        let mut copies: Vec<(DuplicateKey, String)> = Vec::new();
        let mut widest = window;
        for record in reader.records() {
            let record = record?;
            if !is_reference_template(&record) {
                continue;
            }
            let insert = record.insert_size();
            widest = widest.max(insert.abs());
            // [start, start + insert) must hold both last0 and last0 + 1.
            let start = record.pos();
            if insert > 0 && start <= last0 as i64 && start + insert > last0 as i64 + 1 {
                let qname = String::from_utf8_lossy(record.qname()).into_owned();
                let key = duplicate_key(&record, &libraries);
                if record.is_duplicate() {
                    copies.push((key, qname));
                } else {
                    representatives.entry(key).or_default().push(qname);
                }
            }
        }
        if widest > window {
            window = widest;
            continue;
        }

        let mut result = ReferenceFragments {
            fragments: representatives.values().flatten().cloned().collect(),
            duplicate_pairs: copies.len() as u64,
            ..Default::default()
        };
        // A marker leaves one unflagged pair per key; should there be more,
        // the copies go with the first by name, so the result doesn't
        // depend on the BAM's order.
        for (key, qname) in copies {
            if let Some(representative) = representatives.get(&key).and_then(|r| r.iter().min()) {
                result.duplicates.entry(representative.clone()).or_default().push(qname);
            }
        }
        return Ok(result);
    }
}

/// Whether an alignment is a mate of a reference template for
/// [`reference_fragments_at`]. Duplicates are kept here (and set apart
/// there), and it is stricter than [`is_countable`] on two points, because
/// it collects *reference molecules* rather than reproducing `samtools
/// depth`: a supplementary alignment at a breakpoint is evidence *of* a
/// junction, not of the reference allele, and a pair that isn't properly
/// paired doesn't traverse the position contiguously either.
fn is_reference_template(record: &bam::Record) -> bool {
    !record.is_unmapped()
        && !record.is_secondary()
        && !record.is_quality_check_failed()
        && !record.is_supplementary()
        && record.is_proper_pair()
}

/// Number of [`reference_fragments_at`] the break after `pos_1based`.
pub fn fragment_depth_at(
    bam_path: &str,
    chrom: &str,
    pos_1based: u64,
) -> Result<u32, Box<dyn Error>> {
    Ok(reference_fragments_at(bam_path, chrom, pos_1based)?.len() as u32)
}

/// Draw a fragment crossing `junction_index`, as random shearing of a
/// fusion molecule would produce it, or `None` if [`MAX_FRAGMENT_ATTEMPTS`]
/// draws all failed.
///
/// The fragment's length comes from the library's insert-size gaussian
/// (`Normal(mean_insert, std_insert)`), weighted by that length: a molecule
/// crosses a given point with probability proportional to its length, so
/// the fragments over a breakpoint are longer on average than the library
/// as a whole — the reference fragments `fusion-spike` counts at the
/// breakpoint carry the same bias. Its position is then uniform over every
/// start that keeps at least one base on each side of the junction.
///
/// Nothing else is filtered: whether the junction ends up inside a read
/// (split read, possibly with a short anchor) or between the two (discordant
/// pair) is left to chance, as in a real library, and recorded by
/// [`classify_support`]. A draw only fails when the fragment doesn't fit in
/// the junction sequence, i.e. beyond `mean + 4 std` or near a contig end.
fn pick_fragment_spanning_junction(
    seq_len: usize,
    junction_index: usize,
    rng: &mut RandomGenerator,
    mean_insert: f64,
    std_insert: f64,
) -> Option<(usize, usize)> {
    let longest = (mean_insert + 4.0 * std_insert).max(2.0);
    for _ in 0..MAX_FRAGMENT_ATTEMPTS {
        let size_insert = rng.normal(mean_insert, std_insert).round();
        // One base on each side of the junction at the very least.
        if size_insert < 2.0 {
            continue;
        }
        // Length bias, by rejection: keep a length L with probability
        // L / longest.
        if rng.unit() * longest > size_insert {
            continue;
        }
        let size_insert = size_insert as i64;

        // The range of `start` values for which [start, start+size_insert)
        // keeps at least one base on each side of `junction_index`, clamped
        // to what fits in the junction sequence.
        let lo = (junction_index as i64 - size_insert + 1).max(0);
        let hi = (junction_index as i64 - 1).min(seq_len as i64 - size_insert);
        if hi < lo {
            continue;
        }
        let start = rng.range(lo, hi);
        return Some((start as usize, (start + size_insert) as usize));
    }
    None
}

/// How the pair sequenced from fragment `[start, stop)` shows the junction
/// at `junction_index`: only the first and last `length_reads` bases of a
/// fragment are ever read (see [`crate::fastq::build_read_pair`]).
///
/// `min_anchor` is the shortest segment on either side of the junction a
/// short-read aligner can still place as a split alignment (short-read
/// aligners need a minimum seed length — e.g. BWA-MEM's default is 19bp —
/// and SV/fusion callers layer their own minimum overhang on top of that,
/// typically 20-25bp); below it the read usually aligns whole to its
/// majority partner, soft-clipped.
pub fn classify_support(
    start: usize,
    stop: usize,
    junction_index: usize,
    length_reads: usize,
    min_anchor: usize,
) -> Support {
    let take = length_reads.min(stop - start);
    // Shorter of the two segments a read [read_start, read_stop) splits
    // into at the junction, if it crosses it.
    let anchor = |read_start: usize, read_stop: usize| {
        (read_start < junction_index && junction_index < read_stop)
            .then(|| (junction_index - read_start).min(read_stop - junction_index))
    };
    match anchor(start, start + take).max(anchor(stop - take, stop)) {
        Some(a) if a >= min_anchor => Support::SplitRead,
        Some(_) => Support::ShortAnchor,
        None => Support::DiscordantPair,
    }
}

pub struct FusionConfig {
    pub length_reads: usize,
    pub mean_insert_size: f64,
    pub std_insert_size: f64,
    pub min_anchor: usize,
    pub profile_diversity: Option<ProfileDiversity>,
    pub id_diversity: Option<String>,
    pub profile_sequencer: Option<ProfileSequencer>,
}

pub struct FusionGenerator {
    config: FusionConfig,
}

impl FusionGenerator {
    pub fn new(config: FusionConfig) -> Self {
        FusionGenerator { config }
    }

    /// Generate `n` fusion read pairs from `junction`, numbered
    /// `start_number..start_number + n`, alongside the [`TruthRow`] of each
    /// pair actually produced. A read pair whose fragment can't be fitted
    /// in the junction sequence after [`MAX_FRAGMENT_ATTEMPTS`] draws is
    /// skipped (logged), matching `simulate`'s behavior for an
    /// unfulfillable fragment — which is why the truth rows, not `n`, are
    /// what a quantification must be based on.
    #[allow(clippy::too_many_arguments)]
    pub fn generate(
        &self,
        junction: &str,
        junction_index: usize,
        n: u64,
        left_bp: &Breakpoint,
        right_bp: &Breakpoint,
        rng: &mut RandomGenerator,
        start_number: u64,
    ) -> (Vec<FastqRecord>, Vec<FastqRecord>, Vec<TruthRow>) {
        let label = junction_label(left_bp, right_bp);
        let mut forward = Vec::with_capacity(n as usize);
        let mut reverse = Vec::with_capacity(n as usize);
        let mut truth = Vec::with_capacity(n as usize);
        let mut skipped = 0u64;
        let diversity = self
            .config
            .profile_diversity
            .as_ref()
            .zip(self.config.id_diversity.as_ref())
            .and_then(|(profile, id)| profile.get(id));

        for i in 0..n {
            let number = start_number + i;
            let Some((start, stop)) = pick_fragment_spanning_junction(
                junction.len(),
                junction_index,
                rng,
                self.config.mean_insert_size,
                self.config.std_insert_size,
            ) else {
                // Logged at debug rather than warn so a pathological
                // parameter set doesn't flood the log with one line per
                // skipped read; see the aggregated warning below instead.
                log::debug!(
                    "fusion read {number}: no fragment fitting the junction sequence found \
                     after {MAX_FRAGMENT_ATTEMPTS} attempts, skipping"
                );
                skipped += 1;
                continue;
            };

            let sequence = Sequence::new(junction[start..stop].to_string());
            let (mut r1, mut r2) = build_read_pair(
                sequence,
                rng,
                diversity,
                self.config.length_reads,
                self.config.profile_sequencer.as_ref(),
                number,
                "fusion".to_string(),
                label.clone(),
                start as u64,
                stop as u64,
            );
            let qname = build_qname(left_bp, right_bp, number);
            r1.qname = Some(qname.clone());
            r2.qname = Some(qname.clone());
            forward.push(r1);
            reverse.push(r2);
            truth.push(TruthRow {
                qname,
                junction: label.clone(),
                fragment_start: start,
                fragment_stop: stop,
                junction_offset: junction_index - start,
                support: classify_support(
                    start,
                    stop,
                    junction_index,
                    self.config.length_reads,
                    self.config.min_anchor,
                ),
            });
        }

        if skipped > 0 {
            log::warn!("{label}: {skipped} fusion read pair(s) out of {n} requested were skipped");
        }
        (forward, reverse, truth)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn bp(chrom: &str, pos: u64) -> Breakpoint {
        Breakpoint { chrom: chrom.into(), pos, strand: Strand::Forward }
    }

    fn partner<'a>(faidx: &'a FastaIndexedReader, bp: &'a Breakpoint) -> JunctionPartner<'a> {
        JunctionPartner { faidx, bp, variants: &[] }
    }

    fn write_fasta(dir: &std::path::Path, contents: &str) -> String {
        let path = dir.join("ref.fa");
        std::fs::File::create(&path).unwrap().write_all(contents.as_bytes()).unwrap();
        path.to_str().unwrap().to_string()
    }

    #[test]
    fn breakpoint_parse_valid() {
        let bp = Breakpoint::parse("chr1:100").unwrap();
        assert_eq!(bp.chrom, "chr1");
        assert_eq!(bp.pos, 100);
    }

    #[test]
    fn breakpoint_parse_rejects_missing_colon() {
        assert!(Breakpoint::parse("chr1-100").is_err());
    }

    #[test]
    fn breakpoint_parse_rejects_zero_pos() {
        assert!(Breakpoint::parse("chr1:0").is_err());
    }

    #[test]
    fn breakpoint_parse_rejects_non_numeric_pos() {
        assert!(Breakpoint::parse("chr1:abc").is_err());
    }

    #[test]
    fn build_junction_concatenates_flanks_at_expected_index() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), &format!(">chr1\n{}\n", "A".repeat(20)));
        let faidx = FastaIndexedReader::open(&path).unwrap();

        // Breakpoint at 1-based pos 10: left flank is bases [1..=10] (10
        // bases if flank_len=10), right flank is bases [11..=15] (5 bases).
        let left_bp = bp("chr1", 10);
        let right_bp = bp("chr1", 10);
        let (junction, idx) = build_junction(&partner(&faidx, &left_bp), &partner(&faidx, &right_bp), 10).unwrap();
        assert_eq!(idx, 10);
        assert_eq!(junction.len(), 20); // left=10 (bases 1..=10), right=10 (bases 11..=20)
    }

    #[test]
    fn build_junction_orients_each_partner_by_its_strand() {
        let dir = tempfile::tempdir().unwrap();
        // 1-based:        1234567890
        let path = write_fasta(dir.path(), ">chr1\nAAACCCGGGT\n>chr2\nTTTGGGCCCA\n");
        let faidx = FastaIndexedReader::open(&path).unwrap();
        let mut five = bp("chr1", 5);
        let mut three = bp("chr2", 5);

        // +/+: chr1[1..=5] then chr2[6..=10].
        let (junction, idx) =
            build_junction(&partner(&faidx, &five), &partner(&faidx, &three), 10).unwrap();
        assert_eq!((junction.as_str(), idx), ("AAACCGCCCA", 5));

        // -/+: the 5' part of a '-' gene is the reverse complement of the
        // bases after the break, chr1[6..=10] = CGGGT -> ACCCG.
        five.strand = Strand::Reverse;
        let (junction, idx) =
            build_junction(&partner(&faidx, &five), &partner(&faidx, &three), 10).unwrap();
        assert_eq!((junction.as_str(), idx), ("ACCCGGCCCA", 5));

        // -/-: the 3' part of a '-' gene is the reverse complement of the
        // bases up to the break, chr2[1..=5] = TTTGG -> CCAAA.
        three.strand = Strand::Reverse;
        let (junction, idx) =
            build_junction(&partner(&faidx, &five), &partner(&faidx, &three), 10).unwrap();
        assert_eq!((junction.as_str(), idx), ("ACCCGCCAAA", 5));
    }

    #[test]
    fn build_junction_carries_the_sample_variants_before_orienting() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), ">chr1\nAAACCCGGGT\n>chr2\nTTTGGGCCCA\n");
        let faidx = FastaIndexedReader::open(&path).unwrap();
        let mut five = bp("chr1", 5);
        five.strand = Strand::Reverse;
        let three = bp("chr2", 5);
        // chr1 0-based 6 (1-based 7, a G) -> T on the forward strand, which
        // reads A once the '-' partner is reverse-complemented.
        let variants = [Variant {
            chrom: "chr1".into(),
            pos0: 6,
            reference: "G".into(),
            allele: crate::haplotype::Allele { base: b'T', insertion: vec![], deletion: 0 },
            homozygous: true,
        }];
        let five_partner = JunctionPartner { faidx: &faidx, bp: &five, variants: &variants };
        let (junction, _) = build_junction(&five_partner, &partner(&faidx, &three), 10).unwrap();
        assert_eq!(junction, "ACCAGGCCCA");
    }

    #[test]
    fn build_junction_rejects_a_partner_with_nothing_on_its_side() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), ">chr1\nAAACCCGGGT\n");
        let faidx = FastaIndexedReader::open(&path).unwrap();
        // A '-' 5' partner breaking at the contig's last base has no bases
        // after the break to contribute.
        let mut five = bp("chr1", 10);
        five.strand = Strand::Reverse;
        let three = bp("chr1", 5);
        assert!(build_junction(&partner(&faidx, &five), &partner(&faidx, &three), 10).is_err());
    }

    #[test]
    fn parse_junction_label_round_trips() {
        let mut five = bp("chr2", 42522656);
        five.strand = Strand::Reverse;
        let three = bp("chr2", 29446394);
        let (a, b) = parse_junction_label(&junction_label(&five, &three)).unwrap();
        assert_eq!((a, b), (five, three));
        assert!(parse_junction_label("chr1:10>chr2:20").is_err());
    }

    #[test]
    fn fragment_reads_match_the_pool_generator() {
        let junction = "AAAACCCCGGGGTTTT";
        // Fragment [2, 12): 6 bases of 5' partner (junction at 8).
        let (head, tail) = fragment_reads(junction, 8, 6, 10, 4).unwrap();
        assert_eq!(head, "AACC");
        assert_eq!(tail, "CCCC"); // revcomp of GGGG
        assert!(fragment_reads(junction, 8, 9, 10, 4).is_none());
    }

    #[test]
    fn carry_variants_keeps_the_pool_errors() {
        // Position 1: sample SNV (C>T). Position 3: pool sequencing error.
        assert_eq!(carry_variants("ACGA", "ACGT", "ATGT"), "ATGA");
        assert_eq!(carry_variants("acgt", "ACGT", "ATGT"), "atgt");
    }

    #[test]
    fn five_prime_position_drops_the_strand() {
        assert_eq!(five_prime_position("chr2:42522656(-)>chr2:29446394(+)"), "chr2:42522656");
    }

    #[test]
    fn junction_label_carries_both_strands() {
        let mut left = bp("chr2", 42522656);
        left.strand = Strand::Reverse;
        let right = bp("chr2", 29446394);
        assert_eq!(junction_label(&left, &right), "chr2:42522656(-)>chr2:29446394(+)");
    }

    #[test]
    fn build_junction_rejects_n_bases_in_the_left_flank() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), ">chr1\nACGTNACGTACGTACGTACGT\n");
        let faidx = FastaIndexedReader::open(&path).unwrap();
        // pos 5 (1-based) is the 'N': the left flank [1..=5] contains it.
        let bp = bp("chr1", 5);
        assert!(build_junction(&partner(&faidx, &bp), &partner(&faidx, &bp), 5).is_err());
    }

    #[test]
    fn build_junction_rejects_n_bases_in_the_right_flank() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), ">chr1\nACGTACGTNACGTACGTACGT\n");
        let faidx = FastaIndexedReader::open(&path).unwrap();
        // pos 8 (1-based): the right flank [9..=13] contains the 'N' at pos 9.
        let bp = bp("chr1", 8);
        assert!(build_junction(&partner(&faidx, &bp), &partner(&faidx, &bp), 5).is_err());
    }

    #[test]
    fn fetch_flank_clamps_at_contig_start() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), ">chr1\nACGTACGTAC\n");
        let faidx = FastaIndexedReader::open(&path).unwrap();
        let bp = bp("chr1", 3);
        // Requesting 10 bases to the left of pos 3 can only yield 3 (bases 1..=3).
        let left = fetch_flank(&faidx, &bp, 10, Side::Left).unwrap().1;
        assert_eq!(left, "ACG");
    }

    #[test]
    fn fetch_flank_clamps_at_contig_end() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), ">chr1\nACGTACGTAC\n");
        let faidx = FastaIndexedReader::open(&path).unwrap();
        let bp = bp("chr1", 8);
        // Requesting 10 bases to the right of pos 8 can only yield 2 (bases 9..=10).
        let right = fetch_flank(&faidx, &bp, 10, Side::Right).unwrap().1;
        assert_eq!(right, "AC");
    }

    #[test]
    fn fetch_flank_right_is_empty_past_contig_end() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), ">chr1\nACGTACGTAC\n");
        let faidx = FastaIndexedReader::open(&path).unwrap();
        let bp = bp("chr1", 10);
        let right = fetch_flank(&faidx, &bp, 10, Side::Right).unwrap().1;
        assert_eq!(right, "");
    }

    #[test]
    fn pick_fragment_always_crosses_the_junction() {
        let mut rng = RandomGenerator::new(42);
        let (seq_len, junction_index) = (800, 400);
        for _ in 0..500 {
            let (start, stop) =
                pick_fragment_spanning_junction(seq_len, junction_index, &mut rng, 300.0, 50.0)
                    .unwrap();
            assert!(stop <= seq_len);
            assert!(start < junction_index && junction_index < stop, "[{start},{stop})");
        }
    }

    #[test]
    fn pick_fragment_does_not_filter_on_where_the_junction_lands() {
        // Long fragments with short reads: most of them strand the junction
        // in the unsequenced insert, and those must be kept — they are the
        // discordant pairs a real library holds.
        let mut rng = RandomGenerator::new(5);
        let mut counts = [0usize; 3];
        for _ in 0..2000 {
            let (start, stop) =
                pick_fragment_spanning_junction(800, 400, &mut rng, 300.0, 30.0).unwrap();
            match classify_support(start, stop, 400, 50, 20) {
                Support::SplitRead => counts[0] += 1,
                Support::ShortAnchor => counts[1] += 1,
                Support::DiscordantPair => counts[2] += 1,
            }
        }
        assert!(counts.iter().all(|&c| c > 0), "{counts:?}");
        assert!(counts[2] > counts[0], "{counts:?}");
    }

    #[test]
    fn pick_fragment_is_length_biased() {
        // A molecule crosses a point with probability proportional to its
        // length: the drawn sizes average above the library mean, by
        // std^2 / mean for a gaussian (here 100^2 / 300 ~ 33).
        let mut rng = RandomGenerator::new(9);
        let n = 20000;
        let mut total = 0usize;
        for _ in 0..n {
            let (start, stop) =
                pick_fragment_spanning_junction(2000, 1000, &mut rng, 300.0, 100.0).unwrap();
            total += stop - start;
        }
        let mean = total as f64 / n as f64;
        assert!((mean - 333.0).abs() < 8.0, "mean {mean}");
    }

    #[test]
    fn pick_fragment_gives_up_when_no_fragment_fits() {
        // Fragments of ~300 bases can't fit in a 100-base junction sequence.
        let mut rng = RandomGenerator::new(1);
        assert!(pick_fragment_spanning_junction(100, 50, &mut rng, 300.0, 0.001).is_none());
    }

    #[test]
    fn classify_support_tells_split_short_anchor_and_discordant_apart() {
        // Reads of 50 over fragments starting at 0, junction at 400.
        // Junction 30 bases into the head read: split.
        assert_eq!(classify_support(370, 670, 400, 50, 20), Support::SplitRead);
        // 10 bases into the head read: short anchor.
        assert_eq!(classify_support(390, 690, 400, 50, 20), Support::ShortAnchor);
        // 30 bases before the end of the tail read: split.
        assert_eq!(classify_support(130, 430, 400, 50, 20), Support::SplitRead);
        // In the insert between the reads: discordant.
        assert_eq!(classify_support(250, 550, 400, 50, 20), Support::DiscordantPair);
        // Exactly on the head read's end: its reads each sit on one side.
        assert_eq!(classify_support(350, 650, 400, 50, 20), Support::DiscordantPair);
        // A fragment shorter than a read: both reads cover all of it.
        assert_eq!(classify_support(380, 420, 400, 50, 20), Support::SplitRead);
    }

    #[test]
    fn generator_produces_requested_pair_count() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), &format!(">chr1\n{}\n", "ACGT".repeat(100)));
        let faidx = FastaIndexedReader::open(&path).unwrap();
        let left_bp = bp("chr1", 200);
        let right_bp = bp("chr1", 200);
        let (junction, idx) = build_junction(&partner(&faidx, &left_bp), &partner(&faidx, &right_bp), 150).unwrap();

        let config = FusionConfig {
            length_reads: 50,
            mean_insert_size: 150.0,
            std_insert_size: 20.0,
            min_anchor: 20,
            profile_diversity: None,
            id_diversity: None,
            profile_sequencer: None,
        };
        let mut rng = RandomGenerator::new(7);
        let generator = FusionGenerator::new(config);
        let (fwd, rev, truth) =
            generator.generate(&junction, idx, 10, &left_bp, &right_bp, &mut rng, 0);
        assert_eq!(fwd.len(), rev.len());
        assert_eq!(fwd.len(), truth.len());
        assert_eq!(fwd.len(), 10);
        // Both mates carry the same self-identifying, whitespace-free name.
        for row in &truth {
            assert!(row.qname.starts_with(QNAME_PREFIX));
            assert!(!row.qname.contains(char::is_whitespace));
            assert_eq!(row.junction, "chr1:200(+)>chr1:200(+)");
            assert_eq!(row.junction_offset, idx - row.fragment_start);
            assert_eq!(
                row.support,
                classify_support(row.fragment_start, row.fragment_stop, idx, 50, 20)
            );
        }
    }

    #[test]
    fn generated_mates_share_one_whitespace_free_qname() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), &format!(">chr1\n{}\n", "ACGT".repeat(100)));
        let faidx = FastaIndexedReader::open(&path).unwrap();
        let bp = bp("chr1", 200);
        let (junction, idx) = build_junction(&partner(&faidx, &bp), &partner(&faidx, &bp), 150).unwrap();

        let config = FusionConfig {
            length_reads: 50,
            mean_insert_size: 150.0,
            std_insert_size: 20.0,
            min_anchor: 20,
            profile_diversity: None,
            id_diversity: None,
            profile_sequencer: None,
        };
        let mut rng = RandomGenerator::new(7);
        let generator = FusionGenerator::new(config);
        let (fwd, rev, truth) = generator.generate(&junction, idx, 5, &bp, &bp, &mut rng, 0);

        for ((r1, r2), row) in fwd.iter().zip(rev.iter()).zip(truth.iter()) {
            let header = format!("{r1}");
            let header = header.lines().next().unwrap();
            assert_eq!(header, format!("@{}", row.qname));
            assert_eq!(r1.qname, r2.qname);
            assert_eq!(format!("{r2}").lines().next().unwrap(), header);
        }
    }

    #[test]
    fn qname_is_unique_per_pair_and_encodes_both_breakpoints() {
        let left = bp("chr9", 130854064);
        let right = bp("chr22", 23632600);
        let a = build_qname(&left, &right, 0);
        let b = build_qname(&left, &right, 1);
        assert_eq!(a, "hmnrr_fusion_chr9-130854064_chr22-23632600_0000000000");
        assert_ne!(a, b);
        // The reciprocal junction is a different name, not the same one.
        assert_ne!(a, build_qname(&right, &left, 0));
    }

    #[test]
    fn depth_at_counts_overlapping_reads() {
        let dir = tempfile::tempdir().unwrap();
        let sam_path = dir.path().join("reads.sam");
        let sam = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n\
            r1\t0\tchr1\t1\t60\t50M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
            r2\t0\tchr1\t1\t60\t50M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
        std::fs::write(&sam_path, sam).unwrap();

        let bam_path = dir.path().join("reads.bam");
        {
            let sam_reader = bam::Reader::from_path(&sam_path).unwrap();
            let header = bam::Header::from_template(sam_reader.header());
            let mut writer =
                bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
            let mut sam_reader = bam::Reader::from_path(&sam_path).unwrap();
            let mut record = bam::Record::new();
            while let Some(result) = sam_reader.read(&mut record) {
                result.unwrap();
                writer.write(&record).unwrap();
            }
        }
        bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();

        let depth = depth_at(bam_path.to_str().unwrap(), "chr1", 10).unwrap();
        assert_eq!(depth, 2);
    }

    #[test]
    fn depth_at_excludes_secondary_supplementary_and_duplicate_reads() {
        let dir = tempfile::tempdir().unwrap();
        let sam_path = dir.path().join("reads.sam");
        // r1: primary, countable. r2: secondary (flag 256). r3: duplicate
        // (flag 1024). r4: QC-fail (flag 512). Only r1 should be counted.
        let sam = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n\
            r1\t0\tchr1\t1\t60\t50M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
            r2\t256\tchr1\t1\t60\t50M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
            r3\t1024\tchr1\t1\t60\t50M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
            r4\t512\tchr1\t1\t60\t50M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
        std::fs::write(&sam_path, sam).unwrap();

        let bam_path = dir.path().join("reads.bam");
        {
            let sam_reader = bam::Reader::from_path(&sam_path).unwrap();
            let header = bam::Header::from_template(sam_reader.header());
            let mut writer =
                bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
            let mut sam_reader = bam::Reader::from_path(&sam_path).unwrap();
            let mut record = bam::Record::new();
            while let Some(result) = sam_reader.read(&mut record) {
                result.unwrap();
                writer.write(&record).unwrap();
            }
        }
        bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();

        let depth = depth_at(bam_path.to_str().unwrap(), "chr1", 10).unwrap();
        assert_eq!(depth, 1);
    }

    #[test]
    fn fragment_depth_counts_a_molecule_once_however_many_mates_overlap() {
        let dir = tempfile::tempdir().unwrap();
        let sam_path = dir.path().join("reads.sam");
        // f1's two mates both cover position 10, f2's second mate doesn't,
        // and s1 isn't properly paired: 2 fragments, 3 pileup reads.
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
        let qual = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
        let sam = format!(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n\
             f1\t99\tchr1\t1\t60\t50M\t=\t5\t54\t{seq}\t{qual}\n\
             f2\t99\tchr1\t1\t60\t50M\t=\t500\t549\t{seq}\t{qual}\n\
             s1\t65\tchr1\t1\t60\t50M\t=\t900\t949\t{seq}\t{qual}\n\
             f1\t147\tchr1\t5\t60\t50M\t=\t1\t-54\t{seq}\t{qual}\n\
             f2\t147\tchr1\t500\t60\t50M\t=\t1\t-549\t{seq}\t{qual}\n"
        );
        std::fs::write(&sam_path, sam).unwrap();

        let bam_path = dir.path().join("reads.bam");
        {
            let sam_reader = bam::Reader::from_path(&sam_path).unwrap();
            let header = bam::Header::from_template(sam_reader.header());
            let mut writer =
                bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
            let mut sam_reader = bam::Reader::from_path(&sam_path).unwrap();
            let mut record = bam::Record::new();
            while let Some(result) = sam_reader.read(&mut record) {
                result.unwrap();
                writer.write(&record).unwrap();
            }
        }
        bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();

        let path = bam_path.to_str().unwrap();
        // f1 (twice over) + f2 + the improperly paired s1.
        assert_eq!(depth_at(path, "chr1", 10).unwrap(), 4);
        // f1 counted once, f2 once, s1 not at all.
        assert_eq!(fragment_depth_at(path, "chr1", 10).unwrap(), 2);
    }

    #[test]
    fn reference_fragments_include_those_whose_insert_crosses_the_break() {
        let dir = tempfile::tempdir().unwrap();
        let sam_path = dir.path().join("reads.sam");
        // Break after 1-based 200. i1's reads sit at 1-50 and 451-500:
        // neither covers the break, its insert does. c1 covers it with a
        // read. e1's template ends at 200, the last base before the break,
        // so it doesn't cross it.
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
        let qual = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
        let sam = format!(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n\
             i1\t99\tchr1\t1\t60\t50M\t=\t451\t500\t{seq}\t{qual}\n\
             e1\t99\tchr1\t101\t60\t50M\t=\t151\t100\t{seq}\t{qual}\n\
             e1\t147\tchr1\t151\t60\t50M\t=\t101\t-100\t{seq}\t{qual}\n\
             c1\t99\tchr1\t181\t60\t50M\t=\t631\t500\t{seq}\t{qual}\n\
             i1\t147\tchr1\t451\t60\t50M\t=\t1\t-500\t{seq}\t{qual}\n\
             c1\t147\tchr1\t631\t60\t50M\t=\t181\t-500\t{seq}\t{qual}\n"
        );
        std::fs::write(&sam_path, sam).unwrap();

        let bam_path = dir.path().join("reads.bam");
        {
            let sam_reader = bam::Reader::from_path(&sam_path).unwrap();
            let header = bam::Header::from_template(sam_reader.header());
            let mut writer =
                bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
            let mut sam_reader = bam::Reader::from_path(&sam_path).unwrap();
            let mut record = bam::Record::new();
            while let Some(result) = sam_reader.read(&mut record) {
                result.unwrap();
                writer.write(&record).unwrap();
            }
        }
        bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();

        let fragments = reference_fragments_at(bam_path.to_str().unwrap(), "chr1", 200).unwrap();
        let mut fragments: Vec<_> = fragments.fragments.into_iter().collect();
        fragments.sort();
        assert_eq!(fragments, vec!["c1".to_string(), "i1".to_string()]);
    }

    #[test]
    fn unclipped_five_prime_adds_the_clips_back() {
        let cigar = CigarString::try_from("5S40M2D5M3S").unwrap();
        // Forward: 5' end is the start, minus the leading soft clip.
        assert_eq!(unclipped_five_prime(100, false, &cigar.0), 95);
        // Reverse: 5' end is the last aligned base (100 + 47 - 1), plus
        // the trailing soft clip.
        assert_eq!(unclipped_five_prime(100, true, &cigar.0), 149);
    }

    #[test]
    fn reference_fragments_group_duplicates_under_their_representative() {
        let dir = tempfile::tempdir().unwrap();
        let sam_path = dir.path().join("reads.sam");
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
        let qual = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
        // f1 aligns from 101 with 5 soft-clipped bases, d1 from 96 without:
        // both start at 96 unclipped, a duplicate of f1 (flag 1024). d2 has
        // the same ends but comes from another library: not a copy of f1.
        // u1 is flagged duplicate but matches no representative.
        let sam = format!(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n\
             @RG\tID:a\tLB:libA\n@RG\tID:a2\tLB:libA\n@RG\tID:b\tLB:libB\n\
             d1\t1123\tchr1\t96\t60\t50M\t=\t301\t255\t{seq}\t{qual}\tRG:Z:a2\n\
             d2\t1123\tchr1\t96\t60\t50M\t=\t301\t255\t{seq}\t{qual}\tRG:Z:b\n\
             u1\t1123\tchr1\t98\t60\t50M\t=\t301\t253\t{seq}\t{qual}\tRG:Z:a\n\
             f1\t99\tchr1\t101\t60\t5S45M\t=\t301\t250\t{seq}\t{qual}\tRG:Z:a\n\
             d1\t1171\tchr1\t301\t60\t50M\t=\t96\t-255\t{seq}\t{qual}\tRG:Z:a2\n\
             d2\t1171\tchr1\t301\t60\t50M\t=\t96\t-255\t{seq}\t{qual}\tRG:Z:b\n\
             f1\t147\tchr1\t301\t60\t50M\t=\t101\t-250\t{seq}\t{qual}\tRG:Z:a\n\
             u1\t1171\tchr1\t301\t60\t50M\t=\t98\t-253\t{seq}\t{qual}\tRG:Z:a\n"
        );
        std::fs::write(&sam_path, sam).unwrap();

        let bam_path = dir.path().join("reads.bam");
        {
            let sam_reader = bam::Reader::from_path(&sam_path).unwrap();
            let header = bam::Header::from_template(sam_reader.header());
            let mut writer =
                bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
            let mut sam_reader = bam::Reader::from_path(&sam_path).unwrap();
            let mut record = bam::Record::new();
            while let Some(result) = sam_reader.read(&mut record) {
                result.unwrap();
                writer.write(&record).unwrap();
            }
        }
        bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();

        let reference = reference_fragments_at(bam_path.to_str().unwrap(), "chr1", 120).unwrap();
        assert_eq!(reference.len(), 1);
        assert_eq!(reference.duplicate_pairs, 3);
        assert_eq!(reference.with_duplicates("f1").collect::<Vec<_>>(), vec!["f1", "d1"]);
    }

    #[test]
    fn fragment_depth_is_zero_without_coverage() {
        let dir = tempfile::tempdir().unwrap();
        let sam_path = dir.path().join("reads.sam");
        let sam = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n\
            f1\t99\tchr1\t1\t60\t10M\t=\t20\t29\tACGTACGTAC\tIIIIIIIIII\n";
        std::fs::write(&sam_path, sam).unwrap();

        let bam_path = dir.path().join("reads.bam");
        {
            let sam_reader = bam::Reader::from_path(&sam_path).unwrap();
            let header = bam::Header::from_template(sam_reader.header());
            let mut writer =
                bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
            let mut sam_reader = bam::Reader::from_path(&sam_path).unwrap();
            let mut record = bam::Record::new();
            while let Some(result) = sam_reader.read(&mut record) {
                result.unwrap();
                writer.write(&record).unwrap();
            }
        }
        bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();

        assert_eq!(
            fragment_depth_at(bam_path.to_str().unwrap(), "chr1", 500).unwrap(),
            0
        );
    }

    #[test]
    fn depth_at_is_zero_without_coverage() {
        let dir = tempfile::tempdir().unwrap();
        let sam_path = dir.path().join("reads.sam");
        let sam = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n\
            r1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n";
        std::fs::write(&sam_path, sam).unwrap();

        let bam_path = dir.path().join("reads.bam");
        {
            let sam_reader = bam::Reader::from_path(&sam_path).unwrap();
            let header = bam::Header::from_template(sam_reader.header());
            let mut writer =
                bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
            let mut sam_reader = bam::Reader::from_path(&sam_path).unwrap();
            let mut record = bam::Record::new();
            while let Some(result) = sam_reader.read(&mut record) {
                result.unwrap();
                writer.write(&record).unwrap();
            }
        }
        bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();

        let depth = depth_at(bam_path.to_str().unwrap(), "chr1", 500).unwrap();
        assert_eq!(depth, 0);
    }
}
