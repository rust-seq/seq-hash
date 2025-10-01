# Changelog

## 0.1.0
- feat: Use `PackedNSeq` wrapper type introduced by `packed-seq` 4.1.0.
- feat: Add `scalar` feature that disabled the `packed-seq` error when native
  compilation is disabled and no SIMD intrinsics are found.
- perf: Add `inline(always)` to `table_lookup` intrinsic for better perf with thin-lto.

## 0.0.1
- Copy from `simd-minimizers`
- Add `KmerHasher` trait; make `k` a parameter of the hasher itself.
- Add generic to `NtHash` to rotate more than 1.
