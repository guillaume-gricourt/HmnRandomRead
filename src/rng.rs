//! Seeded random number source for read generation.
//!
//! Determinism here means "the same seed produces the same Rust output", not
//! bit-compatibility with the old C++ tool: `std::normal_distribution` and
//! `std::uniform_real_distribution` algorithms are implementation-defined in
//! C++ (they already differ between libstdc++ and libc++), so no seeded PRNG
//! choice on the Rust side could reproduce the old byte-exact output anyway.

use rand::distributions::uniform::SampleUniform;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rand_distr::{Distribution, Normal};

pub struct RandomGenerator {
    rng: StdRng,
}

impl RandomGenerator {
    pub fn new(seed: u64) -> Self {
        RandomGenerator {
            rng: StdRng::seed_from_u64(seed),
        }
    }

    /// Draw from a normal distribution with the given mean/standard deviation.
    pub fn normal(&mut self, mean: f64, std: f64) -> f64 {
        Normal::new(mean, std)
            .expect("standard deviation must be finite and non-negative")
            .sample(&mut self.rng)
    }

    /// Draw a uniform value in `[from, thru]` (inclusive on both ends, matching
    /// the original tool's semantics).
    pub fn range<T: SampleUniform + PartialOrd>(&mut self, from: T, thru: T) -> T {
        self.rng.gen_range(from..=thru)
    }

    /// Draw a uniform `f64` in `[0.0, 1.0]`.
    pub fn unit(&mut self) -> f64 {
        self.range(0.0, 1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn same_seed_gives_same_sequence() {
        let mut a = RandomGenerator::new(42);
        let mut b = RandomGenerator::new(42);
        for _ in 0..10 {
            assert_eq!(a.range(0, 1000), b.range(0, 1000));
            assert_eq!(a.normal(500.0, 50.0), b.normal(500.0, 50.0));
        }
    }

    #[test]
    fn range_is_within_bounds() {
        let mut rng = RandomGenerator::new(7);
        for _ in 0..1000 {
            let v = rng.range(3, 9);
            assert!((3..=9).contains(&v));
        }
    }
}
