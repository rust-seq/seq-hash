# seq-hash

[![crates.io](https://img.shields.io/crates/v/seq-hash)](https://crates.io/crates/seq-hash)
[![docs](https://img.shields.io/docsrs/seq-hash)](https://docs.rs/seq-hash)

A SIMD-accelerated library for iterating over k-mer hashes of DNA sequences, building on
[`packed_seq`](https://github.com/rust-seq/packed-seq).
Building block for [`simd-minimizers`](https://github.com/rust-seq/simd-minimizers).

**Paper:**
Please cite the
[`simd-minimizers`](https://github.com/rust-seq/simd-minimizers) paper, for which this
crate was developed:

- SimdMinimizers: Computing Random Minimizers, fast.  
  Ragnar Groot Koerkamp, Igor Martayan.
  SEA 2025 [https://doi.org/10.4230/LIPIcs.SEA.2025.20](doi.org/10.4230/LIPIcs.SEA.2025.20)

## Requirements

This library requires AVX2 or NEON instruction sets, which, on x64, requires
either `target-cpu=native` or `target-cpu=x86-64-v3`.
See [this README](https://github.com/ragnargrootkoerkamp/ensure_simd) for details and [this
blog](https://curiouscoding.nl/posts/distributing-rust-simd-binaries/) for background.
The same restrictions apply when using seq-hash in a larger project.

``` sh
RUSTFLAGS="-C target-cpu=native" cargo run --release
```

## Usage example

Full documentation can be found on [docs.rs](https://docs.rs/seq-hash).

```rust
use packed_seq::{AsciiSeqVec, PackedSeqVec, SeqVec};
use seq_hash::{KmerHasher, NtHasher};

let seq = b"ACGGCAGCGCATATGTAGT";
let packed_seq = PackedSeqVec::from_ascii(seq);

let k = 3;
// Default `NtHasher` is canonical.
let hasher = <NtHasher>::new(k);

// Consider a 'context' of a single kmer.
let hashes: Vec<_> = hasher.hash_kmers_simd(packed_seq.as_slice(), 1).collect();
assert_eq!(hashes.len(), seq.len() - (k-1)
```
