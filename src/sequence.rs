//! Read extraction, reverse-complement, and diversity mutation.

use crate::diversity::Diversity;
use crate::rng::RandomGenerator;

const NUCLEOTIDES: [u8; 4] = *b"ATCG";

/// Reverse-complement a nucleotide string, preserving case and passing
/// through anything outside plain IUPAC codes as `N`/`n`.
pub trait ReverseComplement {
    fn reverse_complement(&self) -> Self;
}

impl ReverseComplement for String {
    fn reverse_complement(&self) -> Self {
        self.bytes()
            .rev()
            .map(complement_base)
            .map(char::from)
            .collect()
    }
}

fn complement_base(b: u8) -> u8 {
    let comp = match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' | b'U' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'S' => b'S',
        b'W' => b'W',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        _ => b'N',
    };
    if b.is_ascii_lowercase() {
        comp.to_ascii_lowercase()
    } else {
        comp
    }
}

/// Pick a random nucleotide from A/C/G/T, optionally excluding the base
/// currently at a position (case of the exclusion is preserved in the pick,
/// mirroring how the original tool avoided re-picking the same base in
/// lower-case-masked regions).
pub fn choose_base(rng: &mut RandomGenerator, exclude: Option<u8>) -> u8 {
    let mut choices: Vec<u8> = NUCLEOTIDES.to_vec();
    if let Some(e) = exclude {
        let target = e.to_ascii_uppercase();
        choices.retain(|&b| b != target);
    }
    let base = choices[rng.range(0, choices.len() - 1)];
    match exclude {
        Some(e) if e.is_ascii_lowercase() => base.to_ascii_lowercase(),
        _ => base,
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Mutation {
    None,
    Substitute,
    Insert,
    Delete,
}

/// A nucleotide sequence fetched from a reference, before it's split into
/// forward/reverse reads.
pub struct Sequence {
    bases: String,
}

impl Sequence {
    pub fn new(bases: String) -> Self {
        Sequence { bases }
    }

    pub fn len(&self) -> usize {
        self.bases.len()
    }

    pub fn is_empty(&self) -> bool {
        self.bases.is_empty()
    }

    pub fn as_str(&self) -> &str {
        &self.bases
    }

    /// Extract a `length`-base read from the start or end of the sequence,
    /// optionally reverse-complemented.
    pub fn sub_read(&self, length: usize, from_start: bool, reverse_complement: bool) -> String {
        let bytes = self.bases.as_bytes();
        let sub = if from_start {
            &bytes[..length.min(bytes.len())]
        } else {
            &bytes[bytes.len().saturating_sub(length)..]
        };
        let s = String::from_utf8_lossy(sub).into_owned();
        if reverse_complement {
            s.reverse_complement()
        } else {
            s
        }
    }

    /// Apply SNP/indel mutations per `diversity`'s rates, in place. Ports the
    /// original algorithm 1:1, including its quirk that an insertion's
    /// length is drawn as a single coin-flip extended by a `while` loop that
    /// never redraws its probability (so it lands on either 1 base or
    /// `max_ins_size + 1` bases, never in between) — preserved here rather
    /// than silently changed, since it affects the simulated diversity
    /// characteristics.
    pub fn make_mutation(&mut self, rng: &mut RandomGenerator, diversity: &Diversity) {
        let len_seq = self.bases.len();
        let mut mutations = vec![Mutation::None; len_seq];
        let mut insertions: Vec<Vec<u8>> = Vec::new();
        let mut deleting = false;

        for i in 0..len_seq {
            if deleting {
                if rng.unit() < diversity.indel_extend {
                    mutations[i] = Mutation::Delete;
                    continue;
                }
                deleting = false;
            }
            if rng.unit() < diversity.mut_rate {
                if rng.unit() >= diversity.indel_frac {
                    mutations[i] = Mutation::Substitute;
                } else if rng.unit() < 0.5 && i > 1 {
                    mutations[i] = Mutation::Delete;
                    deleting = true;
                } else {
                    let mut len_ins = 1usize;
                    let p = rng.unit();
                    while len_ins <= diversity.max_ins_size && p < diversity.indel_extend {
                        len_ins += 1;
                    }
                    mutations[i] = Mutation::Insert;
                    insertions.push((0..len_ins).map(|_| choose_base(rng, None)).collect());
                }
            }
        }

        let mut bytes = std::mem::take(&mut self.bases).into_bytes();
        let mut insertions = insertions.into_iter();
        let mut pos = 0usize;
        for m in mutations {
            match m {
                Mutation::Insert => {
                    let insertion = insertions.next().unwrap_or_default();
                    let ins_len = insertion.len();
                    bytes.splice(pos..pos, insertion);
                    // Matches the original's `pos_seq += len - 1` immediately
                    // followed by the unconditional `pos_seq++` below.
                    pos += ins_len.saturating_sub(1);
                }
                Mutation::Substitute => {
                    bytes[pos] = choose_base(rng, Some(bytes[pos]));
                }
                Mutation::Delete => {
                    bytes.remove(pos);
                    // Wrapping, not checked, subtraction: the unconditional
                    // `pos += 1` right below always restores a valid index
                    // before `pos` is read again, matching the transient
                    // negative value the original's signed counter took on.
                    pos = pos.wrapping_sub(1);
                }
                Mutation::None => {}
            }
            pos = pos.wrapping_add(1);
        }
        self.bases = String::from_utf8(bytes).expect("sequence bytes must remain valid UTF-8");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reverse_complement_preserves_case_and_handles_n() {
        assert_eq!("ACGTacgtN".to_string().reverse_complement(), "NacgtACGT");
    }

    #[test]
    fn sub_read_from_start() {
        let seq = Sequence::new("ACGTACGT".to_string());
        assert_eq!(seq.sub_read(4, true, false), "ACGT");
    }

    #[test]
    fn sub_read_from_end_reverse_complemented() {
        let seq = Sequence::new("AAAATTTT".to_string());
        assert_eq!(seq.sub_read(4, false, true), "AAAA");
    }

    #[test]
    fn sub_read_longer_than_sequence_returns_whole_sequence() {
        let seq = Sequence::new("ACGT".to_string());
        assert_eq!(seq.sub_read(100, true, false), "ACGT");
        assert_eq!(seq.sub_read(100, false, false), "ACGT");
    }

    #[test]
    fn make_mutation_preserves_alphabet_and_never_panics() {
        let diversity = Diversity {
            mut_rate: 0.4,
            indel_frac: 0.5,
            indel_extend: 0.5,
            max_ins_size: 5,
        };
        let mut rng = RandomGenerator::new(1);
        for _ in 0..50 {
            let mut seq = Sequence::new("ACGTACGTACGTACGTACGTACGTACGT".to_string());
            seq.make_mutation(&mut rng, &diversity);
            assert!(seq.as_str().bytes().all(|b| NUCLEOTIDES.contains(&b)));
        }
    }
}
