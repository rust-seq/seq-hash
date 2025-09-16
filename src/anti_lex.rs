//! A fast implementation of 'anti lexicographic' hashing:
//! A kmer's hash found by reading it's characters right to left, and by inverting the last (most significant) character.
//! When k > 16, only the last 16 characters are used.

use std::cmp::min;

use packed_seq::Delay;

use crate::{KmerHasher, S};

pub struct AntiLexHasher<const RC: bool> {
    k: usize,
    /// Number of bits of each character.
    b: usize,
    /// Number of bits to shift each new character up to make it the most significant one.
    shift: u32,
    /// Mask to flip the bits of the most significant character.
    anti: u32,
    /// Mask to keep only the lowest k*b bits.
    mask: u32,
}

impl<const RC: bool> AntiLexHasher<RC> {
    pub fn new(k: usize) -> Self {
        let b = 2;
        let shift = if b * k <= 32 { b * (k - 1) } else { 32 - b } as u32;
        let anti = ((1 << b) - 1) << shift;
        let mask = if b * k < 32 {
            (1 << (b * k)) - 1
        } else {
            u32::MAX
        };
        Self {
            k,
            b,
            shift,
            anti,
            mask,
        }
    }
}

impl KmerHasher for AntiLexHasher<false> {
    const RC: bool = false;

    fn k(&self) -> usize {
        self.k
    }

    fn mapper(&self) -> impl FnMut(u8) -> u32 {
        let mut fw: u32 = 0;
        move |a| {
            fw = (fw >> self.b) ^ ((a as u32) << (32 - self.b));
            fw ^ self.anti
        }
    }

    fn in_out_mapper_scalar(&self) -> impl FnMut((u8, u8)) -> u32 {
        let mut fw: u32 = 0;
        move |(a, _r)| {
            fw = (fw >> self.b) ^ ((a as u32) << self.shift);
            fw ^ self.anti
        }
    }

    fn in_out_mapper_simd(&self) -> impl FnMut((S, S)) -> S {
        let mut fw: S = S::splat(0);
        move |(a, _r)| {
            fw = (fw >> self.b as u32) ^ (a << self.shift);
            fw ^ S::splat(self.anti)
        }
    }
}

impl KmerHasher for AntiLexHasher<true> {
    const RC: bool = true;

    fn k(&self) -> usize {
        self.k
    }

    fn delay(&self) -> Delay {
        Delay(self.k.saturating_sub(32 / self.b))
    }

    fn mapper(&self) -> impl FnMut(u8) -> u32 {
        unimplemented!();
        |_| unreachable!()
    }

    fn in_out_mapper_scalar(&self) -> impl FnMut((u8, u8)) -> u32 {
        let mut fw: u32 = 0;
        let mut rc: u32 = 0;
        move |(a, r)| {
            fw = (fw >> self.b) ^ ((a as u32) << self.shift);
            // ^2 for complement.
            rc = ((rc << self.b) & self.mask) ^ (r as u32 ^ 2);
            min(fw ^ self.anti, rc ^ self.anti)
        }
    }

    fn in_out_mapper_simd(&self) -> impl FnMut((S, S)) -> S {
        let mut fw: S = S::splat(0);
        let mut rc: S = S::splat(0);
        move |(a, r)| {
            fw = (fw >> self.b as u32) ^ (a << self.shift);
            rc = ((rc << self.b as u32) & S::splat(self.mask)) ^ (r ^ S::splat(2));
            (fw ^ S::splat(self.anti)).min(rc ^ S::splat(self.anti))
        }
    }
}
