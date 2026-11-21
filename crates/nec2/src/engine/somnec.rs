//! Sommerfeld integrals: interpolation tables for the field of a current
//! element over a lossy ground.
//!
//! The tables hold four field functions on three grids in distance and
//! elevation angle, for image distances up to one wavelength. They are
//! computed by numerical integration in the complex plane with Shanks'
//! convergence acceleration.
//!
//! Ported from `somnec.c` of nec2c 1.3.1 (itself a translation of the SOMNEC
//! program). Two places follow the original Fortran instead, because nec2c
//! 1.3.1 departs from it (see `docs/nec2c-deviations.md`): the choice of the
//! integration path in `evlua` and the convergence test of `gshank`.

use num_complex::Complex64;

use super::kernel::test;
use super::{C0, C1, all_finite, pass_bound};
use crate::consts::{
    ACCS, C2, C3, CONST1, CONST4, CRIT, CVEL, GAMMA, MAXH, NM, NTS, P10, P11, P20, P21, PI, PI10, POF, PTP, Q10, Q11,
    Q20, Q21, SP, TOSP, TP, cmplx,
};
use crate::error::{Error, Result};

/// (1 - 2 gamma) / (2 pi), the constant of the small-argument series of Y0'.
const HANKEL_C1: f64 = -0.0245785095;

/// Row lengths and layer sizes of the three grids.
pub(crate) const NDA: [i32; 3] = [11, 17, 9];
pub(crate) const NDPA: [i32; 3] = [110, 85, 72];

const J: Complex64 = Complex64::new(0.0, 1.0);

/// The three interpolation grids (COMMON /GGRID/).
#[derive(Debug, Clone)]
pub(crate) struct Grid {
    pub nxa: [i32; 3],
    pub nya: [i32; 3],
    pub dxa: [f64; 3],
    pub dya: [f64; 3],
    pub xsa: [f64; 3],
    pub ysa: [f64; 3],
    ar: [Vec<Complex64>; 3],
}

impl Default for Grid {
    fn default() -> Self {
        Grid {
            nxa: [11, 17, 9],
            nya: [10, 5, 8],
            dxa: [0.02, 0.05, 0.1],
            dya: [0.1745329252, 0.0872664626, 0.1745329252],
            xsa: [0.0, 0.2, 0.2],
            ysa: [0.0, 0.0, 0.3490658504],
            ar: [vec![C0; 11 * 10 * 4], vec![C0; 17 * 5 * 4], vec![C0; 9 * 8 * 4]],
        }
    }
}

impl Grid {
    /// Values of grid `k`: index `ir + ith * nx + function * nx * ny`.
    pub(crate) fn table(&self, k: usize) -> &[Complex64] {
        &self.ar[k]
    }
}

/// State of the table interpolation `intrp`: the grid region and polynomial
/// coefficients used by the previous call.
#[derive(Debug, Clone)]
pub(crate) struct Interp {
    pub ix: i32,
    pub iy: i32,
    pub ixs: i32,
    pub iys: i32,
    pub igrs: i32,
    pub ixeg: i32,
    pub iyeg: i32,
    pub nxm2: i32,
    pub nym2: i32,
    pub nxms: i32,
    pub nyms: i32,
    pub nd: i32,
    pub ndp: i32,
    pub dx: f64,
    pub dy: f64,
    pub xs: f64,
    pub ys: f64,
    pub xz: f64,
    pub yz: f64,
    pub a: [[Complex64; 4]; 4],
    pub b: [[Complex64; 4]; 4],
    pub c: [[Complex64; 4]; 4],
    pub d: [[Complex64; 4]; 4],
}

impl Default for Interp {
    fn default() -> Self {
        Interp {
            ix: 0,
            iy: 0,
            ixs: -10,
            iys: -10,
            igrs: -10,
            ixeg: 0,
            iyeg: 0,
            nxm2: 0,
            nym2: 0,
            nxms: 0,
            nyms: 0,
            nd: 0,
            ndp: 0,
            dx: 1.0,
            dy: 1.0,
            xs: 0.0,
            ys: 0.0,
            xz: 0.0,
            yz: 0.0,
            a: [[C0; 4]; 4],
            b: [[C0; 4]; 4],
            c: [[C0; 4]; 4],
            d: [[C0; 4]; 4],
        }
    }
}

/// Sommerfeld interpolation tables and their interpolation state.
#[derive(Debug, Clone, Default)]
pub(crate) struct Tables {
    pub grid: Grid,
    pub interp: Interp,
    /// Complex relative permittivity the tables were computed for.
    pub epscf: Complex64,
}

/// Sommerfeld attenuation function for numerical distance `p`.
pub(crate) fn fbar(p: Complex64) -> Complex64 {
    let mut z = J * p.sqrt();
    if z.norm() <= 3.0 {
        // Series expansion.
        let zs = z * z;
        let mut sum = z;
        let mut pow = z;
        for i in 1..=100 {
            pow = -pow * zs / f64::from(i);
            let term = pow / (2.0 * f64::from(i) + 1.0);
            sum += term;
            if term.norm_sqr() / sum.norm_sqr() < ACCS {
                break;
            }
        }
        return 1.0 - (1.0 - sum * TOSP) * z * zs.exp() * SP;
    }

    // Asymptotic expansion.
    let minus = z.re < 0.0;
    if minus {
        z = -z;
    }
    let zs = 0.5 / (z * z);
    let mut sum = C0;
    let mut term = C1;
    for i in 1..=6 {
        term = -term * (2.0 * f64::from(i) - 1.0) * zs;
        sum += term;
    }
    if minus {
        sum -= 2.0 * SP * z * (z * z).exp();
    }
    -sum
}

/// Evaluation of the Sommerfeld integrals for one ground (COMMON /EVLCOM/
/// and /CNTOUR/).
struct Evaluator {
    /// 0 for the Bessel-function form of the integrals, 1 for the Hankel form.
    jh: i32,
    ck2: f64,
    ck2sq: f64,
    tkmag: f64,
    tsmag: f64,
    ck1r: f64,
    zph: f64,
    rho: f64,
    ct1: Complex64,
    ct2: Complex64,
    ct3: Complex64,
    ck1: Complex64,
    ck1sq: Complex64,
    cksm: Complex64,
    /// Ends of the current integration interval in the complex lambda plane.
    a: Complex64,
    b: Complex64,
    /// Series coefficients and term counts of the Bessel and Hankel functions.
    a1: [f64; 25],
    a2: [f64; 25],
    a3: [f64; 25],
    a4: [f64; 25],
    bessel_terms: [usize; 101],
    hankel_terms: [usize; 101],
    /// Shanks' algorithm did not converge at some point.
    not_converged: bool,
    /// An integral is not finite.
    fault: bool,
}

impl Evaluator {
    fn new(epscf: Complex64) -> Self {
        let ck2 = TP;
        let ck2sq = ck2 * ck2;
        // The integrals use exp(-jwt) while NEC uses exp(+jwt), hence the
        // conjugate of the permittivity; the fields are conjugated in `evlua`.
        let ck1sq = ck2sq * epscf.conj();
        let ck1 = ck1sq.sqrt();
        let mut erv = ck1sq * ck1sq;
        let mut ezv = ck2sq * ck2sq;
        let ct2 = 0.125 * (erv - ezv);
        erv *= ck1sq;
        ezv *= ck2sq;
        let ct3 = 0.0625 * (erv - ezv);

        let mut a1 = [0.0; 25];
        let mut a2 = [0.0; 25];
        let mut a3 = [0.0; 25];
        let mut a4 = [0.0; 25];
        let mut psi = -GAMMA;
        for k in 1..=25 {
            let kf = f64::from(k);
            let i = (k - 1) as usize;
            a1[i] = -0.25 / (kf * kf);
            a2[i] = 1.0 / (kf + 1.0);
            psi += 1.0 / kf;
            a3[i] = psi + psi;
            a4[i] = (psi + psi + 1.0 / (kf + 1.0)) / (kf + 1.0);
        }
        let mut bessel_terms = [0; 101];
        let mut hankel_terms = [0; 101];
        for i in 1..=101 {
            let mut tst = 1.0;
            let mut last = 0;
            for (k, a1k) in a1.iter().enumerate().take(24) {
                last = k;
                tst *= -f64::from(i) * a1k;
                if tst < 1.0e-6 {
                    break;
                }
            }
            bessel_terms[(i - 1) as usize] = last + 1;

            let mut tst = 1.0;
            let mut last = 0;
            for k in 0..24 {
                last = k;
                tst *= -f64::from(i) * a1[k];
                if tst * a3[k] < 1.0e-6 {
                    break;
                }
            }
            hankel_terms[(i - 1) as usize] = last + 1;
        }

        Evaluator {
            jh: 0,
            ck2,
            ck2sq,
            tkmag: 100.0 * ck1.norm(),
            tsmag: 100.0 * ck1.norm_sqr(),
            ck1r: ck1.re,
            zph: 0.0,
            rho: 0.0,
            ct1: 0.5 * (ck1sq - ck2sq),
            ct2,
            ct3,
            ck1,
            ck1sq,
            cksm: ck2sq / (ck1sq + ck2sq),
            a: C0,
            b: C0,
            a1,
            a2,
            a3,
            a4,
            bessel_terms,
            hankel_terms,
            not_converged: false,
            fault: false,
        }
    }

    /// Bessel function J0 and its derivative for complex argument.
    fn bessel(&self, z: Complex64) -> (Complex64, Complex64) {
        let zms = z.norm_sqr();
        if zms <= 1.0e-12 {
            return (C1, -0.5 * z);
        }

        let mut blend = false;
        let (mut j0x, mut j0px) = (C0, C0);
        if zms <= 37.21 {
            blend = zms > 36.0;
            // Series expansion.
            let miz = self.bessel_terms[zms as usize];
            let mut j0 = C1;
            let mut j0p = C1;
            let mut zk = C1;
            let zi = z * z;
            for k in 0..miz {
                zk *= self.a1[k] * zi;
                j0 += zk;
                j0p += self.a2[k] * zk;
            }
            j0p *= -0.5 * z;
            if !blend {
                return (j0, j0p);
            }
            j0x = j0;
            j0px = j0p;
        }

        // Asymptotic expansion.
        let zi = 1.0 / z;
        let mut zi2 = zi * zi;
        let p0z = 1.0 + (P20 * zi2 - P10) * zi2;
        let p1z = 1.0 + (P11 - P21 * zi2) * zi2;
        let q0z = (Q20 * zi2 - Q10) * zi;
        let q1z = (Q11 - Q21 * zi2) * zi;
        let mut zk = (J * (z - POF)).exp();
        zi2 = 1.0 / zk;
        let cz = 0.5 * (zk + zi2);
        let sz = J * 0.5 * (zi2 - zk);
        zk = C3 * zi.sqrt();
        let j0 = zk * (p0z * cz - q0z * sz);
        let j0p = -zk * (p1z * sz + q1z * cz);
        if !blend {
            return (j0, j0p);
        }

        // Blend the two expansions where they overlap.
        let w = ((zms.sqrt() - 6.0) * PI10).cos();
        (0.5 * (j0x * (1.0 + w) + j0 * (1.0 - w)), 0.5 * (j0px * (1.0 + w) + j0p * (1.0 - w)))
    }

    /// Hankel function of the first kind and order zero, and its derivative,
    /// for a non-zero complex argument.
    fn hankel(&self, z: Complex64) -> (Complex64, Complex64) {
        let zms = z.norm_sqr();
        debug_assert!(zms != 0.0, "the Hankel function is not defined at zero");

        let mut blend = false;
        let (mut h0x, mut h0px) = (C0, C0);
        if zms <= 16.81 {
            blend = zms > 16.0;
            // Series expansion.
            let miz = self.hankel_terms[zms as usize];
            let mut j0 = C1;
            let mut j0p = C1;
            let mut y0 = C0;
            let mut y0p = C0;
            let mut zk = C1;
            let zi = z * z;
            for k in 0..miz {
                zk *= self.a1[k] * zi;
                j0 += zk;
                j0p += self.a2[k] * zk;
                y0 += self.a3[k] * zk;
                y0p += self.a4[k] * zk;
            }
            j0p *= -0.5 * z;
            let clogz = (0.5 * z).ln();
            y0 = (2.0 * j0 * clogz - y0) / PI + C2;
            y0p = (2.0 / z + 2.0 * j0p * clogz + 0.5 * y0p * z) / PI + HANKEL_C1 * z;
            let h0 = j0 + J * y0;
            let h0p = j0p + J * y0p;
            if !blend {
                return (h0, h0p);
            }
            h0x = h0;
            h0px = h0p;
        }

        // Asymptotic expansion.
        let zi = 1.0 / z;
        let zi2 = zi * zi;
        let p0z = 1.0 + (P20 * zi2 - P10) * zi2;
        let p1z = 1.0 + (P11 - P21 * zi2) * zi2;
        let q0z = (Q20 * zi2 - Q10) * zi;
        let q1z = (Q11 - Q21 * zi2) * zi;
        let zk = (J * (z - POF)).exp() * zi.sqrt() * C3;
        let h0 = zk * (p0z + J * q0z);
        let h0p = J * zk * (p1z + J * q1z);
        if !blend {
            return (h0, h0p);
        }

        let w = ((zms.sqrt() - 4.0) * 31.41592654).cos();
        (0.5 * (h0x * (1.0 + w) + h0 * (1.0 - w)), 0.5 * (h0px * (1.0 + w) + h0p * (1.0 - w)))
    }

    /// Integrands of the six Sommerfeld integrals at parameter `t` along the
    /// interval from `a` to `b`, for source and observer above the ground.
    fn saoa(&self, t: f64) -> [Complex64; 6] {
        let dxl = self.b - self.a;
        let xl = self.a + dxl * t;
        let (ck1, ck1sq, ck2, ck2sq) = (self.ck1, self.ck1sq, self.ck2, self.ck2sq);
        let (mut b0, mut b0p, mut cgam1, mut cgam2);

        if self.jh == 0 {
            // Bessel function form.
            (b0, b0p) = self.bessel(xl * self.rho);
            b0 *= 2.0;
            b0p *= 2.0;
            cgam1 = (xl * xl - ck1sq).sqrt();
            cgam2 = (xl * xl - ck2sq).sqrt();
            if cgam1.re == 0.0 {
                cgam1 = cmplx(0.0, -cgam1.im.abs());
            }
            if cgam2.re == 0.0 {
                cgam2 = cmplx(0.0, -cgam2.im.abs());
            }
        } else {
            // Hankel function form.
            (b0, b0p) = self.hankel(xl * self.rho);
            let com = xl - ck1;
            cgam1 = (xl + ck1).sqrt() * com.sqrt();
            if com.re < 0.0 && com.im >= 0.0 {
                cgam1 = -cgam1;
            }
            let com = xl - ck2;
            cgam2 = (xl + ck2).sqrt() * com.sqrt();
            if com.re < 0.0 && com.im >= 0.0 {
                cgam2 = -cgam2;
            }
        }

        // cgam2 - cgam1, from a series when the two nearly cancel.
        let series = |sign: f64| {
            let dgam = 1.0 / (xl * xl);
            sign * ((self.ct3 * dgam + self.ct2) * dgam + self.ct1) / xl
        };
        let dgam = if xl.norm_sqr() < self.tsmag {
            cgam2 - cgam1
        } else if xl.im < 0.0 {
            series(1.0)
        } else if xl.re < ck2 {
            series(-1.0)
        } else if xl.re > self.ck1r {
            series(1.0)
        } else {
            cgam2 - cgam1
        };

        let den2 = self.cksm * dgam / (cgam2 * (ck1sq * cgam2 + ck2sq * cgam1));
        let den1 = 1.0 / (cgam1 + cgam2) - self.cksm / cgam2;
        let mut com = dxl * xl * (-cgam2 * self.zph).exp();
        let mut ans = [C0; 6];
        ans[5] = com * b0 * den1 / ck1;
        com *= den2;

        if self.rho != 0.0 {
            b0p /= self.rho;
            ans[0] = -com * xl * (b0p + b0 * xl);
            ans[3] = com * xl * b0p;
        } else {
            ans[0] = -com * xl * xl * 0.5;
            ans[3] = ans[0];
        }
        ans[1] = com * cgam2 * cgam2 * b0;
        ans[2] = -ans[3] * cgam2 * self.rho;
        ans[4] = com * b0;
        ans
    }

    /// Integrates the six integrals from `a` to `b` by variable-interval
    /// Romberg integration, starting with `nx` subintervals. A result that
    /// is not finite sets `self.fault`.
    fn rom1(&mut self, nx: i32) -> [Complex64; 6] {
        const N: usize = 6;
        let s = 1.0;
        let ze = 1.0;
        let ep = s / (1.0e4 * f64::from(NM));
        let zend = ze - ep;
        let mut z = 0.0;
        let mut sum = [C0; N];
        let mut ns = nx;
        let mut nt = 0;
        let mut g1 = self.saoa(z);
        let (mut g3, mut g5) = ([C0; N], [C0; N]);
        let (mut t01, mut t10, mut t20) = ([C0; N], [C0; N], [C0; N]);
        let mut dz = 0.0;
        let mut dzot = 0.0;
        let mut refine = false;

        for _ in 0..pass_bound(NM) {
            if !refine {
                dz = s / f64::from(ns);
                if z + dz > ze {
                    dz = ze - z;
                    if dz <= ep {
                        return self.checked(sum);
                    }
                }
                dzot = dz * 0.5;
                g3 = self.saoa(z + dzot);
                g5 = self.saoa(z + dz);
            }

            // Three-point Romberg result.
            let mut nogo = false;
            for i in 0..N {
                let t00 = (g1[i] + g5[i]) * dzot;
                t01[i] = (t00 + dz * g3[i]) * 0.5;
                t10[i] = (4.0 * t01[i] - t00) / 3.0;
                let (tr, ti) = test(t01[i].re, t10[i].re, t01[i].im, t10[i].im, 0.0);
                if tr > CRIT || ti > CRIT {
                    nogo = true;
                }
            }

            if nogo {
                // Five-point Romberg result.
                let g2 = self.saoa(z + dz * 0.25);
                let g4 = self.saoa(z + dz * 0.75);
                nogo = false;
                for i in 0..N {
                    let t02 = (t01[i] + dzot * (g2[i] + g4[i])) * 0.5;
                    let t11 = (4.0 * t02 - t01[i]) / 3.0;
                    t20[i] = (16.0 * t11 - t10[i]) / 15.0;
                    let (tr, ti) = test(t11.re, t20[i].re, t11.im, t20[i].im, 0.0);
                    if tr > CRIT || ti > CRIT {
                        nogo = true;
                    }
                }

                if nogo {
                    nt = 0;
                    if ns < NM {
                        // Halve the step and reuse the points already computed.
                        ns *= 2;
                        dz = s / f64::from(ns);
                        dzot = dz * 0.5;
                        g5 = g3;
                        g3 = g2;
                        refine = true;
                        continue;
                    }
                    // Step size limited: accept the five-point result.
                }
                for i in 0..N {
                    sum[i] += t20[i];
                }
                nt += 1;
            } else {
                for i in 0..N {
                    sum[i] += t10[i];
                }
                nt += 2;
            }

            z += dz;
            if z > zend {
                return self.checked(sum);
            }
            g1 = g5;
            if nt >= NTS && ns > nx {
                ns /= 2;
                nt = 1;
            }
            refine = false;
        }

        // The integration does not advance.
        self.fault = true;
        [C0; N]
    }

    /// Marks the evaluation as failed if `sum` holds a value that is not finite.
    fn checked(&mut self, sum: [Complex64; 6]) -> [Complex64; 6] {
        if !all_finite(&sum) {
            self.fault = true;
        }
        sum
    }

    /// Integrates the six integrals from `start` to infinity in steps of
    /// `dela`, accelerating the convergence with Shanks' algorithm. `seed`
    /// is the value of the integrals up to `start`. With `break_point` the
    /// step changes to `delb` when the path reaches the real part of `bk`.
    fn gshank(
        &mut self,
        start: Complex64,
        dela: Complex64,
        seed: &[Complex64; 6],
        break_point: Option<(Complex64, Complex64)>,
    ) -> [Complex64; 6] {
        const N: usize = 6;
        let mut q1 = [[C0; MAXH]; N];
        let mut q2 = [[C0; MAXH]; N];
        let mut ans1 = [C0; N];
        let mut ans2 = *seed;
        let mut del = dela;
        let mut pending_break = break_point;
        self.b = start;

        let mut int = 1;
        while int <= MAXH {
            self.a = self.b;
            self.b += del;
            if let Some((bk, delb)) = pending_break {
                if self.b.re >= bk.re {
                    // Hit the break point: reset the seed and start over.
                    pending_break = None;
                    self.b = bk;
                    del = delb;
                    let sum = self.rom1(2);
                    for i in 0..N {
                        ans2[i] += sum[i];
                    }
                    int = 1;
                    continue;
                }
            }
            let sum = self.rom1(2);
            for i in 0..N {
                ans1[i] = ans2[i] + sum[i];
            }

            self.a = self.b;
            self.b += del;
            if let Some((bk, delb)) = pending_break {
                if self.b.re >= bk.re {
                    pending_break = None;
                    self.b = bk;
                    del = delb;
                    let sum = self.rom1(2);
                    for i in 0..N {
                        ans2[i] = ans1[i] + sum[i];
                    }
                    int = 1;
                    continue;
                }
            }
            let sum = self.rom1(2);
            for i in 0..N {
                ans2[i] = ans1[i] + sum[i];
            }

            // Shanks' transformation of the sequence of partial integrals.
            let mut den: f64 = 0.0;
            for i in 0..N {
                let mut as1 = ans1[i];
                let mut as2 = ans2[i];
                for jm in 0..int - 1 {
                    let aa = q2[i][jm];
                    let mut a1 = q1[i][jm] + as1 - 2.0 * aa;
                    if a1.re != 0.0 || a1.im != 0.0 {
                        let a2 = aa - q1[i][jm];
                        a1 = q1[i][jm] - a2 * a2 / a1;
                    } else {
                        a1 = q1[i][jm];
                    }
                    let mut a2 = aa + as2 - 2.0 * as1;
                    if a2.re != 0.0 || a2.im != 0.0 {
                        a2 = aa - (as1 - aa) * (as1 - aa) / a2;
                    } else {
                        a2 = aa;
                    }
                    q1[i][jm] = as1;
                    q2[i][jm] = as2;
                    as1 = a1;
                    as2 = a2;
                }
                q1[i][int - 1] = as1;
                q2[i][int - 1] = as2;
                den = den.max(as2.re.abs() + as2.im.abs());
            }

            // Converged when the last diagonals of the two tables agree.
            let denm = 1.0e-3 * den * CRIT;
            let first = int.saturating_sub(3).max(1) - 1;
            let converged = (first..int).all(|j| {
                (0..N).all(|i| {
                    let a1 = q2[i][j];
                    let tolerance = ((a1.re.abs() + a1.im.abs()) * CRIT).max(denm);
                    let diff = q1[i][j] - a1;
                    diff.re.abs() + diff.im.abs() <= tolerance
                })
            });
            if converged || int == MAXH {
                if !converged {
                    self.not_converged = true;
                }
                let mut sum = [C0; N];
                for i in 0..N {
                    sum[i] = 0.5 * (q1[i][int - 1] + q2[i][int - 1]);
                }
                return sum;
            }
            int += 1;
        }
        unreachable!("the loop returns at the last iteration")
    }

    /// Evaluates the Sommerfeld integrals for the point set in `rho` and
    /// `zph`, choosing the integration contour in the complex lambda plane.
    /// Returns `(erv, ezv, erh, eph)`.
    fn evlua(&mut self) -> (Complex64, Complex64, Complex64, Complex64) {
        let (zph, rho, ck1, ck2) = (self.zph, self.rho, self.ck1, self.ck2);
        let mut del = zph.max(rho);
        let mut ans;

        if zph >= 2.0 * rho {
            // Bessel function form of the Sommerfeld integrals.
            self.jh = 0;
            self.a = C0;
            del = 1.0 / del;
            let mut sum;
            if del > self.tkmag {
                self.b = cmplx(0.1 * self.tkmag, -0.1 * self.tkmag);
                sum = self.rom1(2);
                self.a = self.b;
                self.b = cmplx(del, -del);
                let more = self.rom1(2);
                for i in 0..6 {
                    sum[i] += more[i];
                }
            } else {
                self.b = cmplx(del, -del);
                sum = self.rom1(2);
            }
            let delta = cmplx(PTP * del, 0.0);
            ans = self.gshank(self.b, delta, &sum, None);
        } else {
            // Hankel function form of the Sommerfeld integrals.
            self.jh = 1;
            let cp1 = cmplx(0.0, 0.4 * ck2);
            let cp2 = cmplx(0.6 * ck2, -0.2 * ck2);
            let cp3 = cmplx(1.02 * ck2, -0.2 * ck2);
            self.a = cp1;
            self.b = cp2;
            let mut sum = self.rom1(2);
            self.a = cp2;
            self.b = cp3;
            let more = self.rom1(2);
            for i in 0..6 {
                sum[i] = -(sum[i] + more[i]);
            }

            // Path from the imaginary axis to -infinity.
            let slope = if zph > 0.001 * rho { rho / zph } else { 1000.0 };
            del = PTP / del;
            let delta = cmplx(-1.0, slope) * del / (1.0 + slope * slope).sqrt();
            let delta2 = -delta.conj();
            ans = self.gshank(cp1, delta, &sum, None);

            let mut rmis = rho * (ck1.re - ck2);
            let mut between_branch_cuts = rmis >= 2.0 * ck2 && rho >= 1.0e-10;
            if between_branch_cuts && zph >= 1.0e-10 {
                let bk = cmplx(-zph, rho) * (ck1 - cp3);
                rmis = -bk.re / bk.im.abs();
                if rmis > 4.0 * rho / zph {
                    between_branch_cuts = false;
                }
            }

            if between_branch_cuts {
                // Integrate up between the branch cuts, then to +infinity.
                let cp1 = ck1 - cmplx(0.1, 0.2);
                let cp2 = cp1 + 0.2;
                let bk = cmplx(0.0, del);
                let sum = self.gshank(cp1, bk, &ans, None);
                self.a = cp1;
                self.b = cp2;
                ans = self.rom1(1);
                for i in 0..6 {
                    ans[i] -= sum[i];
                }
                let sum = self.gshank(cp3, bk, &ans, None);
                ans = self.gshank(cp2, delta2, &sum, None);
            } else {
                // Integrate below the branch points, then to +infinity.
                let mut sum = [C0; 6];
                for i in 0..6 {
                    sum[i] = -ans[i];
                }
                let rmis = (ck1.re * 1.01).max(ck2 + 1.0);
                let bk = cmplx(rmis, 0.99 * ck1.im);
                let mut delta = bk - cp3;
                delta *= del / delta.norm();
                ans = self.gshank(cp3, delta, &sum, Some((bk, delta2)));
            }
        }

        ans[5] *= ck1;
        // Conjugate, since NEC uses exp(+jwt).
        let erv = (self.ck1sq * ans[2]).conj();
        let ezv = (self.ck1sq * (ans[1] + self.ck2sq * ans[4])).conj();
        let erh = (self.ck2sq * (ans[0] + ans[5])).conj();
        let eph = -(self.ck2sq * (ans[3] + ans[5])).conj();
        (erv, ezv, erh, eph)
    }
}

impl Tables {
    /// Computes the interpolation tables for a ground of relative
    /// permittivity `epr` and conductivity `sig` (S/m) at `fmhz` MHz. A
    /// negative `sig` is minus the imaginary part of the complex permittivity.
    /// Returns `false` if some integral did not converge, and an error if one
    /// is not finite.
    pub(crate) fn generate(&mut self, epr: f64, sig: f64, fmhz: f64) -> Result<bool> {
        self.epscf = if sig >= 0.0 {
            let wlam = CVEL / fmhz;
            cmplx(epr, -sig * wlam * 59.96)
        } else {
            cmplx(epr, sig)
        };
        let epscf = self.epscf;
        let mut ev = Evaluator::new(epscf);
        // Coefficients interpolated from the previous tables must not be reused.
        self.interp = Interp::default();
        let grid = &mut self.grid;

        for k in 0..3 {
            let nr = grid.nxa[k] as usize;
            let nth = grid.nya[k] as usize;
            let dr = grid.dxa[k];
            let dth = grid.dya[k];
            // Grid 1 starts at r = 0, which is filled analytically below.
            let (mut r, first) = if k == 0 { (grid.xsa[k], 1) } else { (grid.xsa[k] - dr, 0) };

            for ir in first..nr {
                r += dr;
                let mut thet = grid.ysa[k] - dth;
                for ith in 0..nth {
                    thet += dth;
                    ev.rho = r * thet.cos();
                    ev.zph = r * thet.sin();
                    if ev.rho < 1.0e-7 {
                        ev.rho = 1.0e-8;
                    }
                    if ev.zph < 1.0e-7 {
                        ev.zph = 0.0;
                    }

                    let (erv, ezv, erh, eph) = ev.evlua();
                    if ev.fault {
                        return Err(Error::Integration { routine: "ROM1", source: None, observation: None });
                    }
                    let rk = ev.ck2 * r;
                    let con = -CONST1 * r / cmplx(rk.cos(), -rk.sin());
                    let at = ir + ith * nr;
                    let layer = nr * nth;
                    grid.ar[k][at] = erv * con;
                    grid.ar[k][at + layer] = ezv * con;
                    grid.ar[k][at + 2 * layer] = erh * con;
                    grid.ar[k][at + 3 * layer] = eph * con;
                }
            }
        }

        // Grid 1 for r equal to zero.
        let cl2 = -CONST4 * (epscf - 1.0) / (epscf + 1.0);
        let cl1 = cl2 / (epscf + 1.0);
        let ezv = epscf * cl1;
        let dth = grid.dya[0];
        let nth = grid.nya[0] as usize;
        let nr = grid.nxa[0] as usize;
        let layer = nr * nth;
        let mut thet = -dth;
        for ith in 0..nth {
            thet += dth;
            let (erv, erh, eph);
            if ith + 1 != nth {
                let mut tfac2 = thet.cos();
                let tfac1 = (1.0 - thet.sin()) / tfac2;
                tfac2 = tfac1 / tfac2;
                erv = epscf * cl1 * tfac1;
                erh = cl1 * (tfac2 - 1.0) + cl2;
                eph = cl1 * tfac2 - cl2;
            } else {
                erv = C0;
                erh = cl2 - 0.5 * cl1;
                eph = -erh;
            }
            let at = ith * nr;
            grid.ar[0][at] = erv;
            grid.ar[0][at + layer] = ezv;
            grid.ar[0][at + 2 * layer] = erh;
            grid.ar[0][at + 3 * layer] = eph;
        }

        Ok(!ev.not_converged)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rom1_reports_values_that_are_not_finite() {
        let mut ev = Evaluator::new(cmplx(13.0, -0.3));
        ev.rho = 0.1;
        ev.zph = 0.1;
        ev.jh = 0;
        ev.a = C0;
        ev.b = cmplx(1.0, -1.0);
        let sum = ev.rom1(2);
        assert!(!ev.fault && all_finite(&sum));

        ev.a = cmplx(f64::NAN, 0.0);
        ev.rom1(2);
        assert!(ev.fault);
    }
}
