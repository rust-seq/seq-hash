#![allow(unreachable_code)]

mod anti_lex;
mod intrinsics;
mod nthash;

pub use anti_lex::AntiLexHasher;
pub use nthash::{MulHasher, NtHasher};
use packed_seq::{ChunkIt, Delay, PaddedIt, Seq};
use std::iter::{repeat, zip};

type S = wide::u32x8;

pub trait KmerHasher {
    /// True when the hash function is invariant under reverse-complement.
    const RC: bool;

    fn is_canonical(&self) -> bool {
        Self::RC
    }

    /// The hasher is already initialized with the value of k.
    fn k(&self) -> usize;

    /// The delay of the 'out' character when `MAPPER_NEEDS_OUT` is true.
    /// Defaults to `k-1`.
    fn delay(&self) -> Delay {
        Delay(self.k() - 1)
    }

    /// Add one character to the hash.
    fn mapper<'s>(&self, seq: impl Seq<'s>) -> impl FnMut(u8) -> u32;
    /// Hash k-mers given the new character and if needed the one k-1 behind to remove.
    fn in_out_mapper_scalar<'s>(&self, seq: impl Seq<'s>) -> impl FnMut((u8, u8)) -> u32;
    /// Hash k-mers given the new character and if needed the one k-1 behind to remove, for each of 4 lanes.
    fn in_out_mapper_simd<'s>(&self, seq: impl Seq<'s>) -> impl FnMut((S, S)) -> S;

    /// Hash the given sequence/kmer.
    fn hash_kmer<'s>(&self, seq: impl Seq<'s>) -> u32 {
        seq.iter_bp().map(self.mapper(seq)).last().unwrap_or(0)
    }
    /// Hash all non-empty prefixes of the given sequence.
    fn hash_prefixes<'s>(&self, seq: impl Seq<'s>) -> impl ExactSizeIterator<Item = u32> {
        seq.iter_bp().map(self.mapper(seq))
    }

    /// Hash all k-mers in the given sequence.
    fn hash_kmers_scalar<'s>(&self, seq: impl Seq<'s>) -> impl ExactSizeIterator<Item = u32> {
        let k = self.k();
        let delay = self.delay();
        let mut add = seq.iter_bp();
        let mut remove = seq.iter_bp();
        let mut mapper = self.in_out_mapper_scalar(seq);
        zip(add.by_ref().take(delay.0), repeat(0)).for_each(|a| {
            mapper(a);
        });
        zip(add.by_ref(), remove.by_ref())
            .take(k - 1 - delay.0)
            .for_each(|a| {
                mapper(a);
            });
        zip(add, remove).map(mapper)
    }

    /// Hash all k-mers in the given sequence, using 4 lanes in parallel.
    fn hash_kmers_simd<'s>(&self, seq: impl Seq<'s>, context: usize) -> PaddedIt<impl ChunkIt<S>> {
        let k = self.k();
        let delay = self.delay();
        seq.par_iter_bp_delayed(context + k - 1, delay)
            .map(self.in_out_mapper_simd(seq))
            .advance(k - 1)
    }
}
