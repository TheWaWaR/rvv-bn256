#![no_std]
#![allow(clippy::many_single_char_names)]

extern crate alloc;
use alloc::vec::Vec;

#[link(name = "gfp_rvv")]
extern "C" {
    fn gfp_carry();
    fn gfp_neg();
    fn gfp_add();
    fn gfp_sub();
    fn gfp_mul();
}

#[derive(Clone, Eq, PartialEq, Default, Debug)]
struct Gfp([u64; 4]);

const MASK16: u64 = 0x0000ffff;
const MASK32: u64 = 0xffffffff;
const NUM_BYTES: usize = 256 / 8;

#[derive(Debug)]
pub enum Error {
    UnmarshalGfpFailed,
    UnmarshalG1Failed,
}

impl Gfp {
    // r2 is R^2 where R = 2^256 mod p.
    pub const fn zero() -> Gfp {
        Gfp([0, 0, 0, 0])
    }
    // p2 is p, represented as little-endian 64-bit words.
    pub const fn p2() -> Gfp {
        Gfp([
            0x3c208c16d87cfd47,
            0x97816a916871ca8d,
            0xb85045b68181585d,
            0x30644e72e131a029,
        ])
    }
    // np is the negative inverse of p, mod 2^256.
    pub const fn np() -> Gfp {
        Gfp([
            0x87d20782e4866389,
            0x9ede7d651eca6ac9,
            0xd8afcbd01833da80,
            0xf57a22b791888c6b,
        ])
    }
    // rN1 is R^-1 where R = 2^256 mod p.
    pub const fn rN1() -> Gfp {
        Gfp([
            0xed84884a014afa37,
            0xeb2022850278edf8,
            0xcf63e9cfb74492d9,
            0x2e67157159e5c639,
        ])
    }
    // r2 is R^2 where R = 2^256 mod p.
    pub const fn r2() -> Gfp {
        Gfp([
            0xf32cfc5b538afa89,
            0xb5e71911d44501fb,
            0x47ab1eff0a417ff6,
            0x06d89f71cab8351f,
        ])
    }
    // r3 is R^3 where R = 2^256 mod p.
    pub const fn r3() -> Gfp {
        Gfp([
            0xb1cd6dafda1530df,
            0x62f210e6a7283db6,
            0xef7f0b0c0ada0afb,
            0x20fd6e902d592544,
        ])
    }
    // xiToPSquaredMinus1Over3 is ξ^((p²-1)/3) where ξ = i+9.
    pub const fn xiToPSquaredMinus1Over3() -> Gfp {
        Gfp([
            0x3350c88e13e80b9c,
            0x7dce557cdb5e56b9,
            0x6001b4b8b615564a,
            0x2682e617020217e0,
        ])
    }
    // xiTo2PSquaredMinus2Over3 is ξ^((2p²-2)/3) where ξ = i+9 (a cubic root of unity, mod p).
    pub const fn xiTo2PSquaredMinus2Over3() -> Gfp {
        Gfp([
            0x71930c11d782e155,
            0xa6bb947cffbe3323,
            0xaa303344d4741444,
            0x2c3b3f0d26594943,
        ])
    }
    // xiToPSquaredMinus1Over6 is ξ^((1p²-1)/6) where ξ = i+9 (a cubic root of -1, mod p).
    pub const fn xiToPSquaredMinus1Over6() -> Gfp {
        Gfp([
            0xca8d800500fa1bf2,
            0xf0c5d61468b39769,
            0x0e201271ad0d4418,
            0x04290f65bad856e6,
        ])
    }

    pub fn new(x: i64) -> Gfp {
        let mut out = if x > 0 {
            Gfp([x as u64, 0, 0, 0])
        } else {
            let mut out = Gfp::default();
            out.neg(&Gfp([-x as u64, 0, 0, 0]));
            out
        };
        out.mont_encode(&out.clone());
        out
    }

    // Ops
    pub fn carry(&mut self, head: u64) {
        let a = self;
        let mut b = Gfp::default();

        let mut carry: u64 = 0;
        for i in 0..4 {
            let pi = Gfp::p2().0[i];
            let ai = a.0[i];
            let bi = ai.wrapping_sub(pi).wrapping_sub(carry);
            b.0[i] = bi;
            carry = ((pi & (!ai)) | (pi | (!ai)) & bi) >> 63;
        }
        carry &= !head;

        // If b is negative, then return a.
        // Else return b.
        if carry == 0 {
            // println!("  changed: a: {:?}, b: {:?}, head: {}", a, b, head);
            *a = b;
        } else {
            // println!("unchanged: a: {:?}, head: {}", a, head);
        }

        /*
        carry = 0u64.wrapping_sub(carry);
        let ncarry = !carry;
        for i in 0..4 {
            a.0[i] = (a.0[i] & carry) | (b.0[i] & ncarry);
        }
        */
    }
    pub fn neg(&mut self, a: &Gfp) {
        // unsafe {
        //     gfp_neg();
        // }
        let c = self;
        let mut carry: u64 = 0;
        for i in 0..4 {
            let pi = Gfp::p2().0[i];
            let ai = a.0[i];
            let ci = pi.wrapping_sub(ai).wrapping_sub(carry);
            c.0[i] = ci;
            carry = (ai & (!pi) | (ai | (!pi)) & ci) >> 63;
        }
        c.carry(0);
    }
    pub fn add(&mut self, a: &Gfp, b: &Gfp) {
        let c = self;
        let mut carry: u64 = 0;
        for i in 0..4 {
            let ai = a.0[i];
            let bi = b.0[i];
            let ci = ai.wrapping_add(bi).wrapping_add(carry);
            c.0[i] = ci;
            carry = (ai & bi | (ai | bi) & (!ci)) >> 63;
        }
        c.carry(carry);
    }
    pub fn sub(&mut self, a: &Gfp, b: &Gfp) {
        let c = self;
        let mut t = Gfp::default();
        let mut carry: u64 = 0;
        for i in 0..4 {
            let pi = Gfp::p2().0[i];
            let bi = b.0[i];
            let ti = pi.wrapping_sub(bi).wrapping_sub(carry);
            t.0[i] = ti;
            carry = (bi & (!pi) | (bi | (!pi)) & ti) >> 63;
        }

        carry = 0;
        for i in 0..4 {
            let ai = a.0[i];
            let ti = t.0[i];
            let ci = ai.wrapping_add(ti).wrapping_add(carry);
            c.0[i] = ci;
            carry = (ai & ti | (ai | ti) & (!ci)) >> 63;
        }
        c.carry(carry);
    }
    pub fn mul(&mut self, a: &Gfp, b: &Gfp) {
        let c = self;
        let mut tt = Gfp::_mul(a, b);
        let m = Gfp::_half_mul(&Gfp([tt[0], tt[1], tt[2], tt[3]]), &Gfp::np());
        let t = Gfp::_mul(&Gfp([m[0], m[1], m[2], m[3]]), &Gfp::p2());

        let mut carry: u64 = 0;
        for i in 0..8 {
            let tti = tt[i];
            let ti = t[i];
            let zi = tti.wrapping_add(ti).wrapping_add(carry);
            tt[i] = zi;
            carry = (tti & ti | (tti | ti) & (!zi)) >> 63;
        }

        c.0[0] = tt[4];
        c.0[1] = tt[5];
        c.0[2] = tt[6];
        c.0[3] = tt[7];
        c.carry(carry);
    }
    pub fn _half_mul(a: &Gfp, b: &Gfp) -> [u64; 4] {
        let mut buff = [0u64; 18];
        for i in 0..4 {
            let ai = a.0[i];
            let a0 = ai & MASK16;
            let a1 = (ai >> 16) & MASK16;
            let a2 = (ai >> 32) & MASK16;
            let a3 = ai >> 48;

            for j in 0..4 {
                if i + j > 3 {
                    break;
                }
                let bj = b.0[j];
                let b0 = bj & MASK32;
                let b2 = bj >> 32;

                let off = 4 * (i + j);
                buff[off]     += a0 * b0 + 0;
                buff[off + 1] += a1 * b0 + 0;
                buff[off + 2] += a2 * b0 + a0 * b2;
                buff[off + 3] += a3 * b0 + a1 * b2;
                buff[off + 4] += 0       + a2 * b2;
                buff[off + 5] += 0       + a3 * b2;
            }
        }

        for i in 1..4 {
            let shift = 16 * i;
            let mut head: u64 = 0;
            let mut carry: u64 = 0;
            for j in 0..4 {
                let block = 4 * j;
                let xi = buff[block];
                let yi = (buff[block + i] << shift) + head;
                let zi = xi.wrapping_add(yi).wrapping_add(carry);
                buff[block] = zi;
                carry = (xi & yi | (xi | yi) & (!zi)) >> 63;
                head = buff[block + i] >> (64 - shift);
            }
        }

        [buff[0], buff[4], buff[8], buff[12]]
    }
    pub fn _mul(a: &Gfp, b: &Gfp) -> [u64; 8] {
        let mut buff = [0u64; 32];
        for i in 0..4 {
            let ai = a.0[i];
            let a0 = ai & MASK16;
            let a1 = (ai >> 16) & MASK16;
            let a2 = (ai >> 32) & MASK16;
            let a3 = ai >> 48;

            for j in 0..4 {
                let bj = b.0[j];
                let b0 = bj & MASK32;
                let b2 = bj >> 32;

                let off = 4 * (i + j);
                buff[off]     += a0 * b0 + 0  * 0;
                buff[off + 1] += a1 * b0 + 0  * 0;
                buff[off + 2] += a2 * b0 + a0 * b2;
                buff[off + 3] += a3 * b0 + a1 * b2;
                buff[off + 4] += 0  * 0  + a2 * b2;
                buff[off + 5] += 0  * 0  + a3 * b2;
            }
        }

        for i in 1..4 {
            let shift = 16 * i;
            let mut head: u64 = 0;
            let mut carry: u64 = 0;
            for j in 0..8 {
                let block = 4 * j;
                let xi = buff[block];
                let yi = (buff[block + i] << shift) + head;
                let zi = xi.wrapping_add(yi).wrapping_add(carry);
                buff[block] = zi;
                carry = (xi & yi | (xi | yi) & (!zi)) >> 63;
                head = buff[block + i] >> (64 - shift);
            }
        }

        [
            buff[0], buff[4], buff[8], buff[12], buff[16], buff[20], buff[24], buff[28],
        ]
    }

    pub fn mont_encode(&mut self, a: &Gfp) {
        self.mul(a, &Gfp::r2());
    }
    pub fn mont_decode(&mut self, a: &Gfp) {
        self.mul(a, &Gfp([1, 0, 0, 0]));
    }

    pub fn set(&mut self, other: &Gfp) {
        self.0 = other.0;
    }
    pub fn invert(&mut self, f: Gfp) {
        const BITS: [u64; 4] = [
            0x3c208c16d87cfd45,
            0x97816a916871ca8d,
            0xb85045b68181585d,
            0x30644e72e131a029,
        ];
        let mut sum = Gfp::rN1();
        let mut power = f;

        for word in 0..4 {
            for bit in 0..64 {
                if (BITS[word] >> bit) & 1 == 1 {
                    sum.mul(&sum.clone(), &power);
                }
                let power_clone = power.clone();
                power.mul(&power_clone, &power_clone);
            }
        }
        sum.mul(&sum.clone(), &Gfp::r3());
        self.set(&sum);
    }
    pub fn marshal(&self) -> Vec<u8> {
        let mut out = alloc::vec![0u8; 32];
        for w in 0..4 {
            for b in 0..8 {
                out[8 * w + b] = (self.0[3 - w] >> (56 - 8 * b)) as u8;
            }
        }
        out
    }
    pub fn unmarshal(data: &[u8]) -> Result<Gfp, Error> {
        // Unmarshal the bytes into little endian form
        let mut e = [0u64; 4];
        for w in 0..4 {
            for b in 0..8 {
                e[3 - w] += u64::from(data[8 * w + b]) << (56 - 8 * b);
            }
        }
        // Ensure the point respects the curve modulus
        let mut i = 3;
        loop {
            if e[i] < Gfp::p2().0[i] {
                return Ok(Gfp(e));
            }
            if e[i] > Gfp::p2().0[i] {
                return Err(Error::UnmarshalGfpFailed);
            }

            if i == 0 {
                break;
            } else {
                i -= 1;
            }
        }
        Err(Error::UnmarshalGfpFailed)
    }
}

// gfP2 implements a field of size p² as a quadratic extension of the base field
// where i²=-1.
//   value is xi+y.
struct Gfp2 {
    x: Gfp,
    y: Gfp,
}

impl Gfp2 {
    // xiToPMinus1Over6 is ξ^((p-1)/6) where ξ = i+9.
    pub const fn xiToPMinus1Over6() -> Gfp2 {
        Gfp2 {
            x: Gfp([
                0xa222ae234c492d72,
                0xd00f02a4565de15b,
                0xdc2ff3a253dfc926,
                0x10a75716b3899551,
            ]),
            y: Gfp([
                0xaf9ba69633144907,
                0xca6b1d7387afb78a,
                0x11bded5ef08a2087,
                0x02f34d751a1f3a7c,
            ]),
        }
    }
    // xiToPMinus1Over3 is ξ^((p-1)/3) where ξ = i+9.
    pub const fn xiToPMinus1Over3() -> Gfp2 {
        Gfp2 {
            x: Gfp([
                0x6e849f1ea0aa4757,
                0xaa1c7b6d89f89141,
                0xb6e713cdfae0ca3a,
                0x26694fbb4e82ebc3,
            ]),
            y: Gfp([
                0xb5773b104563ab30,
                0x347f91c8a9aa6454,
                0x7a007127242e0991,
                0x1956bcd8118214ec,
            ]),
        }
    }
    // xiToPMinus1Over2 is ξ^((p-1)/2) where ξ = i+9.
    pub const fn xiToPMinus1Over2() -> Gfp2 {
        Gfp2 {
            x: Gfp([
                0xa1d77ce45ffe77c7,
                0x07affd117826d1db,
                0x6d16bd27bb7edc6b,
                0x2c87200285defecc,
            ]),
            y: Gfp([
                0xe4bbdd0c2936b629,
                0xbb30f162e133bacb,
                0x31a9d1b6f9645366,
                0x253570bea500f8dd,
            ]),
        }
    }
    // xiTo2PMinus2Over3 is ξ^((2p-2)/3) where ξ = i+9.
    pub const fn xiTo2PMinus2Over3() -> Gfp2 {
        Gfp2 {
            x: Gfp([
                0x5dddfd154bd8c949,
                0x62cb29a5a4445b60,
                0x37bc870a0c7dd2b9,
                0x24830a9d3171f0fd,
            ]),
            y: Gfp([
                0x7361d77f843abe92,
                0xa5bb2bd3273411fb,
                0x9c941f314b3e2399,
                0x15df9cddbb9fd3ec,
            ]),
        }
    }
}

// gfP6 implements the field of size p⁶ as a cubic extension of gfP2 where τ³=ξ
// and ξ=i+9.
//   value is xτ² + yτ + z
struct Gfp6 {
    x: Gfp2,
    y: Gfp2,
    z: Gfp2,
}

// gfP12 implements the field of size p¹² as a quadratic extension of gfP6
// where ω²=τ.
//   value is xω + y
struct Gfp12 {
    x: Gfp6,
    y: Gfp6,
}

// curvePoint implements the elliptic curve y²=x³+3. Points are kept in Jacobian
// form and t=z² when valid. G₁ is the set of points of this curve on GF(p).
#[derive(Clone, Eq, PartialEq, Default, Debug)]
struct CurvePoint {
    x: Gfp,
    y: Gfp,
    z: Gfp,
    t: Gfp,
}

impl CurvePoint {
    pub fn add(&mut self, a: &Self, b: &Self) {
        // FIXME: implement use rvv
        if a.is_infinity() {
            self.set(b);
            return;
        }
        if b.is_infinity() {
            self.set(a);
            return;
        }

        // See http://hyperelliptic.org/EFD/g1p/auto-code/shortw/jacobian-0/addition/add-2007-bl.op3

        // Normalize the points by replacing a = [x1:y1:z1] and b = [x2:y2:z2]
        // by [u1:s1:z1·z2] and [u2:s2:z1·z2]
        // where u1 = x1·z2², s1 = y1·z2³ and u1 = x2·z1², s2 = y2·z1³
        let mut z12 = Gfp::default();
        let mut z22 = Gfp::default();
        z12.mul(&a.z, &a.z);
        z22.mul(&b.z, &b.z);

        let mut u1 = Gfp::default();
        let mut u2 = Gfp::default();
        u1.mul(&a.x, &z22);
        u2.mul(&b.x, &z12);

        let mut t = Gfp::default();
        let mut s1 = Gfp::default();
        t.mul(&b.z, &z22);
        s1.mul(&a.y, &t);

        let mut s2 = Gfp::default();
        t.mul(&a.z, &z12);
        s2.mul(&b.y, &t);

        // Compute x = (2h)²(s²-u1-u2)
        // where s = (s2-s1)/(u2-u1) is the slope of the line through
        // (u1,s1) and (u2,s2). The extra factor 2h = 2(u2-u1) comes from the value of z below.
        // This is also:
        // 4(s2-s1)² - 4h²(u1+u2) = 4(s2-s1)² - 4h³ - 4h²(2u1)
        //                        = r² - j - 2v
        // with the notations below.
        let mut h = Gfp::default();
        h.sub(&u2, &u1);
        let x_equal = h == Gfp::zero();

        t.add(&h, &h);
        // i = 4h²
        let mut i = Gfp::default();
        i.mul(&t, &t);
        // j = 4h³
        let mut j = Gfp::default();
        j.mul(&h, &i);

        t.sub(&s2, &s1);
        let y_equal = t == Gfp::zero();
        let c = self;
        if x_equal && y_equal {
            c.double(&a);
            return;
        }
        let mut r = Gfp::default();
        r.add(&t, &t);

        let mut v = Gfp::default();
        v.mul(&u1, &i);

        // t4 = 4(s2-s1)²
        let mut t4 = Gfp::default();
        let mut t6 = Gfp::default();
        t4.mul(&r, &r);
        t.add(&v, &v);
        t6.sub(&t4, &j);

        c.x.sub(&t6, &t);

        // Set y = -(2h)³(s1 + s*(x/4h²-u1))
        // This is also
        // y = - 2·s1·j - (s2-s1)(2x - 2i·u1) = r(v-x) - 2·s1·j
        t.sub(&v, &c.x); // t7
        t4.mul(&s1, &j); // t8
        t6.add(&t4, &t4); // t9
        t4.mul(&r, &t); // t10
        c.y.sub(&t4, &t6);

        // Set z = 2(u2-u1)·z1·z2 = 2h·z1·z2
        t.add(&a.z, &b.z); // t11
        t4.mul(&t, &t); // t12
        t.sub(&t4, &z12); // t13
        t4.sub(&t, &z22); // t14
        c.z.mul(&t4, &h);
    }
    pub fn neg(&mut self, a: &Self) {
        self.x.set(&a.x);
        self.y.neg(&a.y);
        self.z.set(&a.z);
        self.t = Gfp::zero();
    }
    pub fn mul(&mut self, a: &Self, scalar: &Gfp) {
        // FIXME: todo
        unimplemented!()
    }
    pub fn double(&mut self, a: &Self) {
        // See http://hyperelliptic.org/EFD/g1p/auto-code/shortw/jacobian-0/doubling/dbl-2009-l.op3
        let mut A = Gfp::default();
        let mut B = Gfp::default();
        let mut C = Gfp::default();
        A.mul(&a.x, &a.x);
        B.mul(&a.y, &a.y);
        C.mul(&B, &B);

        let mut t = Gfp::default();
        let mut t2 = Gfp::default();
        t.add(&a.x, &B);
        t2.mul(&t, &t);
        t.sub(&t2, &A);
        t2.sub(&t, &C);

        let mut d = Gfp::default();
        let mut e = Gfp::default();
        let mut f = Gfp::default();
        d.add(&t2, &t2);
        t.add(&A, &A);
        e.add(&t, &A);
        f.mul(&e, &e);

        t.add(&d, &d);
        self.x.sub(&f, &t);

        t.add(&C, &C);
        t2.add(&t, &t);
        t.add(&t2, &t2);
        self.y.sub(&d, &self.x);
        t2.mul(&e, &self.y);
        self.y.sub(&t2, &t);

        t.mul(&a.y, &a.z);
        self.z.add(&t, &t);
    }
    pub fn set(&mut self, a: &Self) {
        self.x.set(&a.x);
        self.y.set(&a.y);
        self.z.set(&a.z);
        self.t.set(&a.t);
    }
    pub fn make_affine(&mut self) {
        if self.z == Gfp::new(1) {
            return;
        } else if self.z == Gfp::new(0) {
            self.x = Gfp::zero();
            self.y = Gfp::new(1);
            self.t = Gfp::zero();
            return;
        }

        let mut zInv = Gfp::default();
        zInv.invert(self.z.clone());

        let mut t = Gfp::default();
        let mut zInv2 = Gfp::default();
        t.mul(&self.y, &zInv);
        zInv2.mul(&zInv, &zInv);

        self.x.mul(&self.x.clone(), &zInv2);
        self.y.mul(&t, &zInv2);

        self.z = Gfp::new(1);
        self.t = Gfp::new(1);
    }
    pub fn set_infinity(&mut self) {
        self.x = Gfp::zero();
        self.y = Gfp::new(1);
        self.z = Gfp::zero();
        self.t = Gfp::zero();
    }

    pub fn is_on_curve(&mut self) -> bool {
        self.make_affine();
        if self.is_infinity() {
            return true;
        }

        let mut y2 = Gfp::default();
        let mut x3 = Gfp::default();
        y2.mul(&self.y, &self.y);
        x3.mul(&self.x, &self.x);
        x3.mul(&x3.clone(), &self.x);
        x3.add(&x3.clone(), &Gfp::new(3));

        y2 == x3
    }
    pub fn is_infinity(&self) -> bool {
        self.z == Gfp::zero()
    }
}

// twistPoint implements the elliptic curve y²=x³+3/ξ over GF(p²). Points are
// kept in Jacobian form and t=z² when valid. The group G₂ is the set of
// n-torsion points of this curve over GF(p²) (where n = Order)
struct TwistPoint {
    x: Gfp2,
    y: Gfp2,
    z: Gfp2,
    t: Gfp2,
}

// G1 is an abstract cyclic group. The zero value is suitable for use as the
// output of an operation, but cannot be used as an input.
#[derive(Default, Debug)]
struct G1 {
    p: CurvePoint,
}

impl G1 {
    pub fn add(&mut self, a: &G1, b: &G1) {
        self.p.add(&a.p, &b.p);
    }
    pub fn neg(&mut self, a: &G1) {
        self.p.neg(&a.p);
    }
    pub fn set(&mut self, a: &G1) {
        self.p.set(&a.p);
    }
    pub fn marshal(&mut self) -> Vec<u8> {
        self.p.make_affine();
        let mut temp = Gfp::default();
        temp.mont_decode(&self.p.x);
        let mut ret = temp.marshal();
        temp.mont_decode(&self.p.y);
        ret.append(&mut temp.marshal());
        ret
    }
    pub fn unmarshal(m: &[u8]) -> Result<G1, Error> {
        if m.len() < 2 * NUM_BYTES {
            return Err(Error::UnmarshalG1Failed);
        }

        let mut p = CurvePoint {
            x: Gfp::unmarshal(m)?,
            y: Gfp::unmarshal(&m[NUM_BYTES..])?,
            z: Gfp::default(),
            t: Gfp::default(),
        };
        // Encode into Montgomery form and ensure it's on the curve
        p.x.mont_encode(&p.x.clone());
        p.y.mont_encode(&p.y.clone());

        let zero = Gfp::zero();
        if p.x == zero && p.y == zero {
            // This is the point at infinity.
            p.y = Gfp::new(1);
            p.z = Gfp::zero();
            p.t = Gfp::zero();
        } else {
            p.z = Gfp::new(1);
            p.t = Gfp::new(1);
            if !p.is_on_curve() {
                return Err(Error::UnmarshalG1Failed);
            }
        }
        Ok(G1 { p })
    }
}

// G2 is an abstract cyclic group. The zero value is suitable for use as the
// output of an operation, but cannot be used as an input.
struct G2 {
    p: TwistPoint,
}

// GT is an abstract cyclic group. The zero value is suitable for use as the
// output of an operation, but cannot be used as an input.
struct GT {
    p: Gfp12,
}

#[cfg(test)]
mod test {
    use super::*;
    fn run_bn256_add(input: &[u8]) -> Vec<u8> {
        let x = G1::unmarshal(input).unwrap();
        let y = G1::unmarshal(&input[64..]).unwrap();
        let mut res = G1::default();
        res.add(&x, &y);
        res.marshal()
    }

    fn run_test(input_hex: &str, expected_hex: &str) {
        let input = hex::decode(input_hex).unwrap();
        let expected = hex::decode(expected_hex).unwrap();
        assert_eq!(run_bn256_add(&input), expected);
    }

    // == test case from: go-ethereum/core/vm/testdata/precompiles/bn256Add.json

    #[test]
    fn test_chfast1() {
        run_test(
            "18b18acfb4c2c30276db5411368e7185b311dd124691610c5d3b74034e093dc9063c909c4720840cb5134cb9f59fa749755796819658d32efc0d288198f3726607c2b7f58a84bd6145f00c9c2bc0bb1a187f20ff2c92963a88019e7c6a014eed06614e20c147e940f2d70da3f74c9a17df361706a4485c742bd6788478fa17d7",
            "2243525c5efd4b9c3d3c45ac0ca3fe4dd85e830a4ce6b65fa1eeaee202839703301d1d33be6da8e509df21cc35964723180eed7532537db9ae5e7d48f195c915"
        );
    }

    #[test]
    fn test_chfast2() {
        run_test(
            "2243525c5efd4b9c3d3c45ac0ca3fe4dd85e830a4ce6b65fa1eeaee202839703301d1d33be6da8e509df21cc35964723180eed7532537db9ae5e7d48f195c91518b18acfb4c2c30276db5411368e7185b311dd124691610c5d3b74034e093dc9063c909c4720840cb5134cb9f59fa749755796819658d32efc0d288198f37266",
            "2bd3e6d0f3b142924f5ca7b49ce5b9d54c4703d7ae5648e61d02268b1a0a9fb721611ce0a6af85915e2f1d70300909ce2e49dfad4a4619c8390cae66cefdb204"
        );
    }

    #[test]
    fn test_cdetrio1() {
        // FIXME:: bug
        run_test(
            "0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",
            "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",
        );
    }
}
