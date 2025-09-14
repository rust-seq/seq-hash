//! NtHash the kmers in a sequence.
use std::array::from_fn;

use super::intrinsics;
use crate::KmerHasher;
use crate::S;
use packed_seq::complement_base;
use wide::u32x8;

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
    fn new(k: usize) -> Self;
    fn f(&self, b: u8) -> u32;
    fn c(&self, b: u8) -> u32;
    fn f_rot(&self, b: u8) -> u32;
    fn c_rot(&self, b: u8) -> u32;
    fn simd_f(&self, b: u32x8) -> u32x8;
    fn simd_c(&self, b: u32x8) -> u32x8;
    fn simd_f_rot(&self, b: u32x8) -> u32x8;
    fn simd_c_rot(&self, b: u32x8) -> u32x8;
}

const R: u32 = 7;

#[derive(Clone)]
pub struct NtHasher<const RC: bool> {
    f: [u32; 4],
    c: [u32; 4],
    f_rot: [u32; 4],
    c_rot: [u32; 4],
    simd_f: u32x8,
    simd_c: u32x8,
    simd_f_rot: u32x8,
    simd_c_rot: u32x8,
    fw_init: u32,
    rc_init: u32,
}

impl<const RC: bool> NtHasher<RC> {
    pub fn new(k: usize) -> Self {
        CharHasher::new(k)
    }
}

impl<const RC: bool> CharHasher for NtHasher<RC> {
    const RC: bool = RC;
    fn new(k: usize) -> Self {
        // FIXME: Assert that the corresponding sequence type has 2 bits per character.

        let rot = k as u32 - 1;
        let f = HASHES_F;
        let c = from_fn(|i| HASHES_F[complement_base(i as u8) as usize]);
        let f_rot = f.map(|h| h.rotate_left(rot));
        let c_rot = c.map(|h| h.rotate_left(rot));
        let idx = [0, 1, 2, 3, 0, 1, 2, 3];
        let simd_f = idx.map(|i| f[i]).into();
        let simd_c = idx.map(|i| c[i]).into();
        let simd_f_rot = idx.map(|i| f_rot[i]).into();
        let simd_c_rot = idx.map(|i| c_rot[i]).into();

        let mut this = Self {
            f,
            c,
            f_rot,
            c_rot,
            simd_f,
            simd_c,
            simd_f_rot,
            simd_c_rot,
            fw_init: 0,
            rc_init: 0,
        };
        for _ in 0..k - 1 {
            this.fw_init = this.fw_init.rotate_left(R) ^ this.f(0);
            this.rc_init = this.rc_init.rotate_right(R) ^ this.c_rot(0);
        }

        this
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
    rot: u32,
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
    fn new(k: usize) -> Self {
        Self {
            rot: (k as u32 - 1) % 32,
        }
    }

    fn f(&self, b: u8) -> u32 {
        (b as u32).wrapping_mul(C)
    }
    fn c(&self, b: u8) -> u32 {
        (complement_base(b) as u32).wrapping_mul(C)
    }
    fn f_rot(&self, b: u8) -> u32 {
        (b as u32).wrapping_mul(C).rotate_left(self.rot * R)
    }
    fn c_rot(&self, b: u8) -> u32 {
        (complement_base(b) as u32)
            .wrapping_mul(C)
            .rotate_left(self.rot * R)
    }

    fn simd_f(&self, b: u32x8) -> u32x8 {
        b * C.into()
    }
    fn simd_c(&self, b: u32x8) -> u32x8 {
        packed_seq::complement_base_simd(b) * C.into()
    }
    fn simd_f_rot(&self, b: u32x8) -> u32x8 {
        let r = b * C.into();
        let rot = self.rot * R % 32;
        (r << rot) | (r >> (32 - rot))
    }
    fn simd_c_rot(&self, b: u32x8) -> u32x8 {
        let r = packed_seq::complement_base_simd(b) * C.into();
        let rot = self.rot * R % 32;
        (r << rot) | (r >> (32 - rot))
    }
}

impl<CH: CharHasher> KmerHasher for CH {
    const RC: bool = CH::RC;

    fn mapper(&self) -> impl FnMut(u8) -> u32 {
        let mut fw = 0u32;
        let mut rc = 0u32;
        move |a| {
            fw = fw.rotate_left(R) ^ self.f(a);
            if !Self::RC {
                rc = rc.rotate_right(R) ^ self.c_rot(a);
                fw.wrapping_add(rc)
            } else {
                fw
            }
        }
    }

    fn in_out_mapper_scalar(&self, k: usize) -> impl FnMut((u8, u8)) -> u32 {
        let mut fw = 0u32;
        let mut rc = 0u32;
        for _ in 0..k - 1 {
            fw = fw.rotate_left(R) ^ self.f(0);
            if Self::RC {
                rc = rc.rotate_right(R) ^ self.c_rot(0);
            }
        }

        move |(a, r)| {
            let fw_out = fw.rotate_left(R) ^ self.f(a);
            fw = fw_out ^ self.f_rot(r);
            if !Self::RC {
                let rc_out = rc.rotate_right(R) ^ self.c_rot(a);
                rc = rc_out ^ self.c(r);
                fw_out.wrapping_add(rc_out)
            } else {
                fw_out
            }
        }
    }

    fn in_out_mapper_simd(&self, k: usize) -> impl FnMut((S, S)) -> S {
        let mut fw = 0u32;
        let mut rc = 0u32;
        for _ in 0..k - 1 {
            fw = fw.rotate_left(R) ^ self.f(0);
            if Self::RC {
                rc = rc.rotate_right(R) ^ self.c_rot(0);
            }
        }

        let mut fw = S::splat(fw);
        let mut rc = S::splat(rc);

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
