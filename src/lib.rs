#![allow(unreachable_code)]

mod intrinsics;
mod nthash;

use std::iter::{repeat, zip};

pub use nthash::{MulHasher, NtHasher};
use packed_seq::{ChunkIt, Delay, PaddedIt, Seq};
type S = wide::u32x8;

pub trait KmerHasher {
    /// True when the hash function is invariant under reverse-complement.
    const RC: bool;

    /// Add one character to the hash.
    fn mapper(&self) -> impl FnMut(u8) -> u32;
    /// Hash k-mers given the new character and if needed the one k-1 behind to remove.
    fn in_out_mapper_scalar(&self, k: usize) -> impl FnMut((u8, u8)) -> u32;
    /// Hash k-mers given the new character and if needed the one k-1 behind to remove, for each of 4 lanes.
    fn in_out_mapper_simd(&self, k: usize) -> impl FnMut((S, S)) -> S;

    /// Hash the given sequence/kmer.
    fn hash_seq<'s>(&self, seq: impl Seq<'s>) -> u32 {
        seq.iter_bp().map(self.mapper()).last().unwrap_or(0)
    }
    /// Hash all non-empty prefixes of the given sequence.
    fn hash_prefixes<'s>(&self, seq: impl Seq<'s>) -> impl ExactSizeIterator<Item = u32> {
        seq.iter_bp().map(self.mapper())
    }

    /// Hash all k-mers in the given sequence.
    fn hash_kmers_scalar<'s>(&self, k: usize, seq: impl Seq<'s>) -> impl ExactSizeIterator<Item = u32> {
        let mut add = seq.iter_bp();
        let remove = seq.iter_bp();
        let mut mapper = self.in_out_mapper_scalar(k);
        zip(add.by_ref().take(k - 1), repeat(0)).for_each(|a| {
            mapper(a);
        });
        zip(add, remove).map(mapper)
    }

    /// Hash all k-mers in the given sequence, using 4 lanes in parallel.
    fn hash_kmers_simd<'s>(
        &self,
        seq: impl Seq<'s>,
        k: usize,
        context: usize,
    ) -> PaddedIt<impl ChunkIt<S>> {
        seq.par_iter_bp_delayed(context + k - 1, Delay(k - 1))
            .map(self.in_out_mapper_simd(k))
            .advance(k - 1)
    }
}
