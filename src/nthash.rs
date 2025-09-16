//! NtHash the kmers in a sequence.
use std::array::from_fn;
use std::hash::BuildHasher;
use std::hash::BuildHasherDefault;
use std::hash::DefaultHasher;

use super::intrinsics;
use crate::KmerHasher;
use crate::S;
use packed_seq::complement_base;
use packed_seq::Seq;
use wide::u32x8;

type SeedHasher = BuildHasherDefault<DefaultHasher>;

/// Original ntHash seed values.
// TODO: Update to guarantee unique hash values for k<=16?
const HASHES_F: [u32; 4] = [
    0x3c8b_fbb3_95c6_0474u64 as u32,
    0x3193_c185_62a0_2b4cu64 as u32,
    0x2032_3ed0_8257_2324u64 as u32,
    0x2955_49f5_4be2_4456u64 as u32,
];

pub trait CharHasher: Clone {
    const RC: bool;
    const BITS_PER_CHAR: usize;
    fn new(k: usize) -> Self {
        Self::new_with_seed(k, None)
    }
    fn new_with_seed(k: usize, seed: Option<u32>) -> Self;
    fn k(&self) -> usize;
    fn f(&self, b: u8) -> u32;
    fn c(&self, b: u8) -> u32;
    fn f_rot(&self, b: u8) -> u32;
    fn c_rot(&self, b: u8) -> u32;
    fn simd_f(&self, b: u32x8) -> u32x8;
    fn simd_c(&self, b: u32x8) -> u32x8;
    fn simd_f_rot(&self, b: u32x8) -> u32x8;
    fn simd_c_rot(&self, b: u32x8) -> u32x8;

    fn fw_init(&self) -> u32 {
        let mut fw = 0u32;
        for _ in 0..self.k() - 1 {
            fw = fw.rotate_left(R) ^ self.f(0);
        }
        fw
    }
    fn rc_init(&self) -> u32 {
        let mut rc = 0u32;
        for _ in 0..self.k() - 1 {
            rc = rc.rotate_right(R) ^ self.c_rot(0);
        }
        rc
    }
}

const R: u32 = 7;

#[derive(Clone)]
pub struct NtHasher<const RC: bool> {
    k: usize,
    f: [u32; 4],
    c: [u32; 4],
    f_rot: [u32; 4],
    c_rot: [u32; 4],
    simd_f: u32x8,
    simd_c: u32x8,
    simd_f_rot: u32x8,
    simd_c_rot: u32x8,
}

impl<const RC: bool> NtHasher<RC> {
    pub fn new(k: usize) -> Self {
        CharHasher::new(k)
    }
}

impl<const RC: bool> CharHasher for NtHasher<RC> {
    const RC: bool = RC;
    const BITS_PER_CHAR: usize = 2;

    fn new_with_seed(k: usize, seed: Option<u32>) -> Self {
        // FIXME: Assert that the corresponding sequence type has 2 bits per character.

        let rot = k as u32 - 1;
        let hasher = SeedHasher::new();
        let f = match seed {
            None => HASHES_F,
            Some(seed) => from_fn(|i| hasher.hash_one(HASHES_F[i] ^ seed) as u32),
        };
        let c = from_fn(|i| f[complement_base(i as u8) as usize]);
        let f_rot = f.map(|h| h.rotate_left(rot * R));
        let c_rot = c.map(|h| h.rotate_left(rot * R));
        let idx = [0, 1, 2, 3, 0, 1, 2, 3];
        let simd_f = idx.map(|i| f[i]).into();
        let simd_c = idx.map(|i| c[i]).into();
        let simd_f_rot = idx.map(|i| f_rot[i]).into();
        let simd_c_rot = idx.map(|i| c_rot[i]).into();

        Self {
            k,
            f,
            c,
            f_rot,
            c_rot,
            simd_f,
            simd_c,
            simd_f_rot,
            simd_c_rot,
        }
    }

    fn k(&self) -> usize {
        self.k
    }

    fn f(&self, b: u8) -> u32 {
        unsafe { *self.f.get_unchecked(b as usize) }
    }
    fn c(&self, b: u8) -> u32 {
        unsafe { *self.c.get_unchecked(b as usize) }
    }
    fn f_rot(&self, b: u8) -> u32 {
        unsafe { *self.f_rot.get_unchecked(b as usize) }
    }
    fn c_rot(&self, b: u8) -> u32 {
        unsafe { *self.c_rot.get_unchecked(b as usize) }
    }

    fn simd_f(&self, b: u32x8) -> u32x8 {
        intrinsics::table_lookup(self.simd_f, b)
    }
    fn simd_c(&self, b: u32x8) -> u32x8 {
        intrinsics::table_lookup(self.simd_c, b)
    }
    fn simd_f_rot(&self, b: u32x8) -> u32x8 {
        intrinsics::table_lookup(self.simd_f_rot, b)
    }
    fn simd_c_rot(&self, b: u32x8) -> u32x8 {
        intrinsics::table_lookup(self.simd_c_rot, b)
    }
}

#[derive(Clone)]
pub struct MulHasher<const RC: bool> {
    k: usize,
    rot: u32,
    mul: u32,
}

impl<const RC: bool> MulHasher<RC> {
    pub fn new(k: usize) -> Self {
        CharHasher::new(k)
    }
}

// Mixing constant.
const C: u32 = 0x517cc1b727220a95u64 as u32;

impl<const RC: bool> CharHasher for MulHasher<RC> {
    const RC: bool = RC;
    const BITS_PER_CHAR: usize = 8;

    fn new_with_seed(k: usize, seed: Option<u32>) -> Self {
        Self {
            k,
            rot: (k as u32 - 1) % 32,
            mul: C ^ match seed {
                None => 0,
                // don't change parity,
                Some(seed) => (SeedHasher::new().hash_one(seed) as u32) << 1,
            },
        }
    }

    fn k(&self) -> usize {
        self.k
    }

    fn f(&self, b: u8) -> u32 {
        (b as u32).wrapping_mul(self.mul)
    }
    fn c(&self, b: u8) -> u32 {
        (complement_base(b) as u32).wrapping_mul(self.mul)
    }
    fn f_rot(&self, b: u8) -> u32 {
        (b as u32).wrapping_mul(self.mul).rotate_left(self.rot * R)
    }
    fn c_rot(&self, b: u8) -> u32 {
        (complement_base(b) as u32)
            .wrapping_mul(self.mul)
            .rotate_left(self.rot * R)
    }

    fn simd_f(&self, b: u32x8) -> u32x8 {
        b * self.mul.into()
    }
    fn simd_c(&self, b: u32x8) -> u32x8 {
        packed_seq::complement_base_simd(b) * self.mul.into()
    }
    fn simd_f_rot(&self, b: u32x8) -> u32x8 {
        let r = b * self.mul.into();
        let rot = self.rot * R % 32;
        (r << rot) | (r >> (32 - rot))
    }
    fn simd_c_rot(&self, b: u32x8) -> u32x8 {
        let r = packed_seq::complement_base_simd(b) * self.mul.into();
        let rot = self.rot * R % 32;
        (r << rot) | (r >> (32 - rot))
    }
}

impl<CH: CharHasher> KmerHasher for CH {
    const RC: bool = CH::RC;

    fn k(&self) -> usize {
        self.k()
    }

    #[inline(always)]
    fn mapper<'s>(&self, seq: impl Seq<'s>) -> impl FnMut(u8) -> u32 {
        assert!(seq.bits_per_char() <= CH::BITS_PER_CHAR);

        let mut fw = 0u32;
        let mut rc = 0u32;
        move |a| {
            fw = fw.rotate_left(R) ^ self.f(a);
            if Self::RC {
                rc = rc.rotate_right(R) ^ self.c_rot(a);
                fw.wrapping_add(rc)
            } else {
                fw
            }
        }
    }

    #[inline(always)]
    fn in_out_mapper_scalar<'s>(&self, seq: impl Seq<'s>) -> impl FnMut((u8, u8)) -> u32 {
        assert!(seq.bits_per_char() <= CH::BITS_PER_CHAR);

        let mut fw = self.fw_init();
        let mut rc = self.rc_init();

        move |(a, r)| {
            let fw_out = fw.rotate_left(R) ^ self.f(a);
            fw = fw_out ^ self.f_rot(r);
            if Self::RC {
                let rc_out = rc.rotate_right(R) ^ self.c_rot(a);
                rc = rc_out ^ self.c(r);
                fw_out.wrapping_add(rc_out)
            } else {
                fw_out
            }
        }
    }

    #[inline(always)]
    fn in_out_mapper_simd<'s>(&self, seq: impl Seq<'s>) -> impl FnMut((S, S)) -> S {
        assert!(seq.bits_per_char() <= CH::BITS_PER_CHAR);
        let mut fw = S::splat(self.fw_init());
        let mut rc = S::splat(self.rc_init());

        move |(a, r)| {
            let fw_out = ((fw << R) | (fw >> (32 - R))) ^ self.simd_f(a);
            fw = fw_out ^ self.simd_f_rot(r);
            if Self::RC {
                let rc_out = ((rc >> R) | (rc << (32 - R))) ^ self.simd_c_rot(a);
                rc = rc_out ^ self.simd_c(r);
                // Wrapping SIMD add
                fw_out + rc_out
            } else {
                fw_out
            }
        }
    }
}
