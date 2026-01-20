# Changelog

<!-- next-header -->

## git
- perf: slightly improve codegen for `hash_valid_kmers_simd` by avoiding `zip`.
  (Instead just call `.next()` manually.)
- doc: update readme to link to `ensure_simd` docs.
- misc: cargo update

## 0.1.1
- perf: precompute initial hash values.
- perf: more `inline(always)`
- perf: Avoid slow `u32x8::splat` on NEON.
- perf: bump to packed-seq 4.2.0

## 0.1.0
- feat: Use `PackedNSeq` wrapper type introduced by `packed-seq` 4.1.0.
- feat: Add `scalar` feature that disabled the `packed-seq` error when native
  compilation is disabled and no SIMD intrinsics are found.
- perf: Add `inline(always)` to `table_lookup` intrinsic for better perf with thin-lto.

## 0.0.1
- Copy from `simd-minimizers`
- Add `KmerHasher` trait; make `k` a parameter of the hasher itself.
- Add generic to `NtHash` to rotate more than 1.
