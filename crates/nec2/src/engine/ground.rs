//! Ground: parameters, the Norton ground-wave formulas and the field of a
//! segment over a Sommerfeld ground.
//!
//! Ported from `ground.c` (`rom2`, `sflds`), `fields.c` (`gwave`),
//! `calculations.c` (`intrp`) and the ground part of `main.c` of nec2c 1.3.1.

use num_complex::Complex64;

use super::kernel::test;
use super::somnec::{Interp, NDA, NDPA, fbar};
use super::{C0, C1, Engine, all_finite, pass_bound};
use crate::consts::{CONST1, CONST4, PI, POT, TP, TPJ, cmplx};
use crate::error::{Error, Result};

/// How the ground is modeled.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GroundType {
    FreeSpace,
    /// Perfectly conducting ground.
    Perfect,
    /// Finite ground by the reflection-coefficient approximation.
    ReflectionCoefficient,
    /// Finite ground by the Sommerfeld-Norton method.
    Sommerfeld,
}

/// Ground parameters of a GN card.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GroundParams {
    /// 0 reflection coefficients, 1 perfect ground, 2 Sommerfeld-Norton.
    pub iperf: i32,
    /// Number of radial wires of the ground-screen approximation, or 0.
    pub nradl: i32,
    /// Relative permittivity and conductivity (S/m) of the ground. A negative
    /// conductivity is minus the imaginary part of the complex permittivity,
    /// at every frequency.
    pub epsr: f64,
    pub sig: f64,
    /// Length and wire radius of the radials, metres.
    pub scrwlt: f64,
    pub scrwrt: f64,
}

/// The ground as the solution uses it at the current frequency.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GroundInfo {
    pub kind: GroundType,
    /// Number, length and wire radius of the radials, if a screen is modeled.
    pub radials: Option<(i32, f64, f64)>,
    pub epsr: f64,
    /// Conductivity at the current frequency, S/m.
    pub sig: f64,
    /// Complex relative permittivity.
    pub epsc: Complex64,
}

impl Engine {
    /// Sets the ground under the antenna; `None` is free space.
    pub fn set_ground(&mut self, params: Option<GroundParams>) -> Result<()> {
        if params.is_some_and(|p| p.nradl != 0 && p.iperf == 2) {
            return Err(Error::RadialsWithSommerfeld);
        }
        // Nothing of the screen of an earlier ground survives.
        self.gnd.nradl = 0;
        self.gnd.scrwl = 0.0;
        self.gnd.scrwr = 0.0;
        self.gnd.t1 = C0;
        self.gnd.t2 = 0.0;
        self.scrwlt = 0.0;
        self.scrwrt = 0.0;
        let Some(p) = params else {
            self.gnd.ksymp = 1;
            self.gnd.iperf = 0;
            return Ok(());
        };
        self.gnd.iperf = p.iperf;
        self.gnd.ksymp = 2;
        self.epsr = p.epsr;
        self.sig = p.sig;
        // NEC-2 has no screen over a perfect ground.
        if p.nradl > 0 && p.iperf != 1 {
            self.gnd.nradl = p.nradl;
            self.scrwlt = p.scrwlt;
            self.scrwrt = p.scrwrt;
        }
        Ok(())
    }

    /// Whether a ground plane is present.
    pub fn has_ground(&self) -> bool {
        self.gnd.ksymp != 1
    }

    /// Evaluates the ground parameters at the current frequency. For the
    /// Sommerfeld ground this computes the interpolation tables.
    pub fn prepare_ground(&mut self) -> Result<GroundInfo> {
        let mut info =
            GroundInfo { kind: GroundType::FreeSpace, radials: None, epsr: self.epsr, sig: self.sig, epsc: C1 };
        if self.gnd.ksymp == 1 {
            return Ok(info);
        }
        self.gnd.frati = C1;
        if self.gnd.iperf == 1 {
            info.kind = GroundType::Perfect;
            return Ok(info);
        }

        let wlam = self.wlam;
        // A negative conductivity stands for the imaginary part of the
        // permittivity. The value of the card is kept, so that this holds at
        // every frequency (nec2c and the Fortran convert it once).
        let sig = if self.sig < 0.0 { -self.sig / (59.96 * wlam) } else { self.sig };
        let epsc = cmplx(self.epsr, -sig * wlam * 59.96);
        self.gnd.zrati = 1.0 / epsc.sqrt();
        self.gwav.u = self.gnd.zrati;
        self.gwav.u2 = self.gwav.u * self.gwav.u;

        if self.gnd.nradl > 0 {
            self.gnd.scrwl = self.scrwlt / wlam;
            self.gnd.scrwr = self.scrwrt / wlam;
            self.gnd.t1 = cmplx(0.0, 1.0) * 2367.067 / f64::from(self.gnd.nradl);
            self.gnd.t2 = self.gnd.scrwr * f64::from(self.gnd.nradl);
            info.radials = Some((self.gnd.nradl, self.scrwlt, self.scrwrt));
        }

        if self.gnd.iperf != 2 {
            info.kind = GroundType::ReflectionCoefficient;
        } else {
            self.sommerfeld.generate(self.epsr, sig, self.fmhz)?;
            self.gnd.frati = (epsc - 1.0) / (epsc + 1.0);
            if ((self.sommerfeld.epscf - epsc) / epsc).norm() >= 1.0e-3 {
                return Err(Error::GroundParameters);
            }
            info.kind = GroundType::Sommerfeld;
        }
        info.sig = sig;
        info.epsc = epsc;
        Ok(info)
    }

    /// Integrates the field due to the ground over the source segment from
    /// `a` to `b` by variable-interval Romberg integration. The nine results
    /// are the x, y, z field components of the constant, sine and cosine
    /// current terms. If the integrand is not finite the integration stops
    /// and `self.diag.fault` names this routine.
    pub(crate) fn rom2(&mut self, a: f64, b: f64, dmin: f64) -> [Complex64; 9] {
        const NTS: i32 = 4;
        const NX: i32 = 1;
        const RX: f64 = 1.0e-4;
        const N: usize = 9;

        let mut sum = [C0; N];
        let mut z = a;
        let ze = b;
        let s = b - a;
        // Limit of the halving, 65536 as in the NEC-2 Fortran (nec2c uses
        // the number of segments here).
        const NM: i32 = 65536;
        let ep = s / (1.0e4 * f64::from(NM));
        let zend = ze - ep;
        let mut ns = NX;
        let mut nt = 0;
        let mut g1 = self.sflds(z);
        let (mut g3, mut g5) = ([C0; N], [C0; N]);
        let (mut t01, mut t10, mut t20) = ([C0; N], [C0; N], [C0; N]);
        let mut dz = 0.0;
        let mut dzot = 0.0;
        let mut flag = true;

        for _ in 0..pass_bound(NM) {
            if flag {
                dz = s / f64::from(ns);
                if z + dz > ze {
                    dz = ze - z;
                    if dz <= ep {
                        return sum;
                    }
                }
                dzot = dz * 0.5;
                g3 = self.sflds(z + dzot);
                g5 = self.sflds(z + dz);
            }

            // Three-point Romberg result; convergence is judged on the
            // constant-current terms.
            let mut tmag1 = 0.0;
            let mut tmag2 = 0.0;
            for i in 0..N {
                let t00 = (g1[i] + g5[i]) * dzot;
                t01[i] = (t00 + dz * g3[i]) * 0.5;
                t10[i] = (4.0 * t01[i] - t00) / 3.0;
                if i > 2 {
                    continue;
                }
                tmag1 += t01[i].norm_sqr();
                tmag2 += t10[i].norm_sqr();
            }
            let (tr, _) = test(tmag1.sqrt(), tmag2.sqrt(), 0.0, 0.0, dmin);

            if tr <= RX {
                for i in 0..N {
                    sum[i] += t10[i];
                }
                nt += 2;
                z += dz;
                if z > zend {
                    return sum;
                }
                g1 = g5;
                if nt >= NTS && ns > NX {
                    ns /= 2;
                    nt = 1;
                }
                flag = true;
                continue;
            }

            let g2 = self.sflds(z + dz * 0.25);
            let g4 = self.sflds(z + dz * 0.75);
            let mut tmag1 = 0.0;
            let mut tmag2 = 0.0;

            // Five-point Romberg result.
            for i in 0..N {
                let t02 = (t01[i] + dzot * (g2[i] + g4[i])) * 0.5;
                let t11 = (4.0 * t02 - t01[i]) / 3.0;
                t20[i] = (16.0 * t11 - t10[i]) / 15.0;
                if i > 2 {
                    continue;
                }
                tmag1 += t11.norm_sqr();
                tmag2 += t20[i].norm_sqr();
            }
            let (tr, _) = test(tmag1.sqrt(), tmag2.sqrt(), 0.0, 0.0, dmin);

            if tr > RX {
                nt = 0;
                if ns < NM {
                    ns *= 2;
                    dz = s / f64::from(ns);
                    dzot = dz * 0.5;
                    g5 = g3;
                    g3 = g2;
                    flag = false;
                    continue;
                }
                // `%12.5E` of C: the exponent has a sign and two digits.
                let text = format!("{z:.5E}");
                let (mantissa, exponent) = text.split_once('E').expect("the E format has an exponent");
                let exponent: i32 = exponent.parse().expect("the exponent is an integer");
                let sign = if exponent < 0 { '-' } else { '+' };
                let z = format!("{mantissa}E{sign}{:02}", exponent.abs());
                self.diag.warnings.push(format!("ROM2 -- STEP SIZE LIMITED AT Z = {z:>12}"));
            }

            // A value that is not finite fails both convergence tests and
            // arrives here.
            if !all_finite(&t20) {
                break;
            }
            for i in 0..N {
                sum[i] += t20[i];
            }
            nt += 1;
            z += dz;
            if z > zend {
                return sum;
            }
            g1 = g5;
            flag = true;
            if nt < NTS || ns <= NX {
                continue;
            }
            ns /= 2;
            nt = 1;
        }

        // The integrand is not finite, or the integration does not advance.
        self.diag.fault = Some("ROM2");
        [C0; N]
    }

    /// Field due to the ground of a current element on the source segment at
    /// distance `t` from its centre, at the observation point in `self.incom`.
    pub(crate) fn sflds(&mut self, t: f64) -> [Complex64; 9] {
        let (xj, yj, zj) = (self.dj.xj, self.dj.yj, self.dj.zj);
        let (cabj, sabj, salpj, s) = (self.dj.cabj, self.dj.sabj, self.dj.salpj, self.dj.s);
        let (sn, xsn, ysn) = (self.incom.sn, self.incom.xsn, self.incom.ysn);
        let mut e = [C0; 9];

        let xt = xj + t * cabj;
        let yt = yj + t * sabj;
        let zt = zj + t * salpj;
        let mut rhx = self.incom.xo - xt;
        let mut rhy = self.incom.yo - yt;
        let rhs = rhx * rhx + rhy * rhy;
        let rho = rhs.sqrt();
        let (phx, phy);
        if rho <= 0.0 {
            rhx = 1.0;
            rhy = 0.0;
            phx = 0.0;
            phy = 1.0;
        } else {
            rhx /= rho;
            rhy /= rho;
            phx = -rhy;
            phy = rhx;
        }

        let mut cph = rhx * xsn + rhy * ysn;
        let mut sph = rhy * xsn - rhx * ysn;
        if cph.abs() < 1.0e-10 {
            cph = 0.0;
        }
        if sph.abs() < 1.0e-10 {
            sph = 0.0;
        }

        self.gwav.zph = self.incom.zo + zt;
        let zph = self.gwav.zph;
        let zphs = zph * zph;
        let r2s = rhs + zphs;
        self.gwav.r2 = r2s.sqrt();
        let r2 = self.gwav.r2;
        let rk = r2 * TP;
        self.gwav.xx2 = cmplx(rk.cos(), -rk.sin());

        // Norton approximation: the current is lumped at the segment centre
        // with the moment of a constant, sine or cosine distribution.
        if self.incom.isnor != 1 {
            self.gwav.zmh = 1.0;
            self.gwav.r1 = 1.0;
            self.gwav.xx1 = C0;
            let (mut erv, mut ezv, mut erh, mut ezh, mut eph) = self.gwave();

            let mut et = -CONST1 * self.gnd.frati * self.gwav.xx2 / (r2s * r2);
            let er = 2.0 * et * cmplx(1.0, rk);
            et *= cmplx(1.0 - rk * rk, rk);
            let hrv = (er + et) * rho * zph / r2s;
            let hzv = (zphs * er - rhs * et) / r2s;
            let hrh = (rhs * er - zphs * et) / r2s;
            erv -= hrv;
            ezv -= hzv;
            erh += hrh;
            ezh += hrv;
            eph += et;
            erv *= salpj;
            ezv *= salpj;
            erh = erh * sn * cph;
            ezh = ezh * sn * cph;
            eph = eph * sn * sph;
            erh = erv + erh;
            e[0] = (erh * rhx + eph * phx) * s;
            e[1] = (erh * rhy + eph * phy) * s;
            e[2] = (ezv + ezh) * s;
            let mut sfac = PI * s;
            sfac = sfac.sin() / sfac;
            e[6] = e[0] * sfac;
            e[7] = e[1] * sfac;
            e[8] = e[2] * sfac;
            return e;
        }

        // Interpolate in the Sommerfeld field tables.
        let thet = if rho >= 1.0e-12 { (zph / rho).atan() } else { POT };
        let (erv, ezv, erh, eph) = self.intrp(r2, thet);

        // Combine vertical and horizontal components, convert to x, y, z
        // and multiply by exp(-jkr)/r.
        self.gwav.xx2 /= r2;
        let xx2 = self.gwav.xx2;
        let sfac = sn * cph;
        let erh2 = xx2 * (salpj * erv + sfac * erh);
        let ezh = xx2 * (salpj * ezv - sfac * erv);
        let eph2 = sn * sph * xx2 * eph;
        // Constant current.
        e[0] = erh2 * rhx + eph2 * phx;
        e[1] = erh2 * rhy + eph2 * phy;
        e[2] = ezh;
        // Sine current.
        let rk = TP * t;
        let sfac = rk.sin();
        e[3] = e[0] * sfac;
        e[4] = e[1] * sfac;
        e[5] = e[2] * sfac;
        // Cosine current.
        let sfac = rk.cos();
        e[6] = e[0] * sfac;
        e[7] = e[1] * sfac;
        e[8] = e[2] * sfac;
        e
    }

    /// Electric field of a current element over ground, including the ground
    /// wave, by the formulas of K. A. Norton (Proc. IRE, Sept. 1937).
    /// Returns `(erv, ezv, erh, ezh, eph)`.
    pub(crate) fn gwave(&self) -> (Complex64, Complex64, Complex64, Complex64, Complex64) {
        let g = &self.gwav;
        let (u, u2) = (g.u, g.u2);
        let sppp = g.zmh / g.r1;
        let sppp2 = sppp * sppp;
        let cppp2 = (1.0 - sppp2).max(1.0e-20);
        let cppp = cppp2.sqrt();
        let spp = g.zph / g.r2;
        let spp2 = spp * spp;
        let cpp2 = (1.0 - spp2).max(1.0e-20);
        let cpp = cpp2.sqrt();
        let rk1 = -TPJ * g.r1;
        let rk2 = -TPJ * g.r2;
        let t1 = 1.0 - u2 * cpp2;
        let t2 = t1.sqrt();
        let t3 = (1.0 - 1.0 / rk1) / rk1;
        let t4 = (1.0 - 1.0 / rk2) / rk2;
        let p1 = rk2 * u2 * t1 / (2.0 * cpp2);
        let rv = (spp - u * t2) / (spp + u * t2);
        let omr = 1.0 - rv;
        let mut w = 1.0 / omr;
        w = 4.0 * p1 * w * w;
        let f = fbar(w);
        let q1 = rk2 * t1 / (2.0 * u2 * cpp2);
        let rh = (t2 - u * spp) / (t2 + u * spp);
        let mut v = 1.0 / (1.0 + rh);
        v = 4.0 * q1 * v * v;
        let gg = fbar(v);
        let xr1 = g.xx1 / g.r1;
        let xr2 = g.xx2 / g.r2;

        let x1 = cppp2 * xr1;
        let x2 = rv * cpp2 * xr2;
        let x3 = omr * cpp2 * f * xr2;
        let x4 = u * t2 * spp * 2.0 * xr2 / rk2;
        let x5 = xr1 * t3 * (1.0 - 3.0 * sppp2);
        let x6 = xr2 * t4 * (1.0 - 3.0 * spp2);
        let ezv = (x1 + x2 + x3 - x4 - x5 - x6) * (-CONST4);

        let x1 = sppp * cppp * xr1;
        let x2 = rv * spp * cpp * xr2;
        let x3 = cpp * omr * u * t2 * f * xr2;
        let x4 = spp * cpp * omr * xr2 / rk2;
        let x5 = 3.0 * sppp * cppp * t3 * xr1;
        let x6 = cpp * u * t2 * omr * xr2 / rk2 * 0.5;
        let x7 = 3.0 * spp * cpp * t4 * xr2;
        let erv = -(x1 + x2 - x3 + x4 - x5 + x6 - x7) * (-CONST4);
        let ezh = -(x1 - x2 + x3 - x4 - x5 - x6 + x7) * (-CONST4);

        let x1 = sppp2 * xr1;
        let x2 = rv * spp2 * xr2;
        let x4 = u2 * t1 * omr * f * xr2;
        let x5 = t3 * (1.0 - 3.0 * cppp2) * xr1;
        let x6 = t4 * (1.0 - 3.0 * cpp2) * (1.0 - u2 * (1.0 + rv) - u2 * omr * f) * xr2;
        let x7 = u2 * cpp2 * omr * (1.0 - 1.0 / rk2) * (f * (u2 * t1 - spp2 - 1.0 / rk2) + 1.0 / rk2) * xr2;
        let erh = (x1 - x2 - x4 - x5 + x6 + x7) * (-CONST4);

        let x1 = xr1;
        let x2 = rh * xr2;
        let x3 = (rh + 1.0) * gg * xr2;
        let x4 = t3 * xr1;
        let x5 = t4 * (1.0 - u2 * (1.0 + rv) - u2 * omr * f) * xr2;
        let x6 = 0.5 * u2 * omr * (f * (u2 * t1 - spp2 - 1.0 / rk2) + 1.0 / rk2) * xr2 / rk2;
        let eph = -(x1 - x2 + x3 - x4 + x5 + x6) * (-CONST4);

        (erv, ezv, erh, ezh, eph)
    }

    /// Bivariate cubic interpolation in the Sommerfeld tables: the four field
    /// functions at distance `x` (wavelengths) and elevation angle `y`
    /// (radians). The polynomial coefficients of the last 4 by 4 point region
    /// are kept and reused while the point stays in it.
    pub(crate) fn intrp(&mut self, x: f64, y: f64) -> (Complex64, Complex64, Complex64, Complex64) {
        let grid = &self.sommerfeld.grid;
        let st: &mut Interp = &mut self.sommerfeld.interp;

        let mut jump = false;
        if x < st.xs || y < st.ys {
            jump = true;
        } else {
            st.ix = ((x - st.xs) / st.dx) as i32 + 1;
            st.iy = ((y - st.ys) / st.dy) as i32 + 1;
        }

        if st.ix < st.ixeg || st.iy < st.iyeg || (st.ix - st.ixs).abs() >= 2 || (st.iy - st.iys).abs() >= 2 || jump {
            // Determine the grid and the grid region.
            let igr = if x <= grid.xsa[1] {
                0
            } else if y > grid.ysa[2] {
                2
            } else {
                1
            };

            if igr as i32 != st.igrs {
                st.igrs = igr as i32;
                st.dx = grid.dxa[igr];
                st.dy = grid.dya[igr];
                st.xs = grid.xsa[igr];
                st.ys = grid.ysa[igr];
                st.nxm2 = grid.nxa[igr] - 2;
                st.nym2 = grid.nya[igr] - 2;
                st.nxms = ((st.nxm2 + 1) / 3) * 3 + 1;
                st.nyms = ((st.nym2 + 1) / 3) * 3 + 1;
                st.nd = NDA[igr];
                st.ndp = NDPA[igr];
                st.ix = ((x - st.xs) / st.dx) as i32 + 1;
                st.iy = ((y - st.ys) / st.dy) as i32 + 1;
            }

            st.ixs = ((st.ix - 1) / 3) * 3 + 2;
            if st.ixs < 2 {
                st.ixs = 2;
            }
            st.ixeg = -10000;
            if st.ixs > st.nxm2 {
                st.ixs = st.nxm2;
                st.ixeg = st.nxms;
            }

            st.iys = ((st.iy - 1) / 3) * 3 + 2;
            if st.iys < 2 {
                st.iys = 2;
            }
            st.iyeg = -10000;
            if st.iys > st.nym2 {
                st.iys = st.nym2;
                st.iyeg = st.nyms;
            }

            // Coefficients of the 4 cubic polynomials in x for the 4 grid
            // values of y, for each of the 4 functions.
            let ar = grid.table(st.igrs as usize);
            let mut iadz = st.ixs + (st.iys - 3) * st.nd - st.ndp;
            for k in 0..4 {
                iadz += st.ndp;
                let mut iadd = iadz;
                for i in 0..4 {
                    iadd += st.nd;
                    let at = iadd as usize;
                    let p1 = ar[at - 2];
                    let p2 = ar[at - 1];
                    let p3 = ar[at];
                    let p4 = ar[at + 1];
                    st.a[i][k] = (p4 - p1 + 3.0 * (p2 - p3)) * 0.1666666667;
                    st.b[i][k] = (p1 - 2.0 * p2 + p3) * 0.5;
                    st.c[i][k] = p3 - (2.0 * p1 + 3.0 * p2 + p4) * 0.1666666667;
                    st.d[i][k] = p2;
                }
            }

            st.xz = f64::from(st.ixs - 1) * st.dx + st.xs;
            st.yz = f64::from(st.iys - 1) * st.dy + st.ys;
        }

        // Evaluate the polynomials in x and interpolate cubically in y.
        let xx = (x - st.xz) / st.dx;
        let yy = (y - st.yz) / st.dy;
        let eval = |k: usize| {
            let fx = |i: usize| ((st.a[i][k] * xx + st.b[i][k]) * xx + st.c[i][k]) * xx + st.d[i][k];
            let (fx1, fx2, fx3, fx4) = (fx(0), fx(1), fx(2), fx(3));
            let p1 = fx4 - fx1 + 3.0 * (fx2 - fx3);
            let p2 = 3.0 * (fx1 - 2.0 * fx2 + fx3);
            let p3 = 6.0 * fx3 - 2.0 * fx1 - 3.0 * fx2 - fx4;
            ((p1 * yy + p2) * yy + p3) * yy * 0.1666666667 + fx2
        };
        (eval(0), eval(1), eval(2), eval(3))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Geometry;

    fn monopole_over_sommerfeld_ground(epsr: f64) -> (Engine, Result<GroundInfo>) {
        let mut geometry = Geometry::new();
        geometry.wire(1, 3, [0.0, 0.0, 0.1], [0.0, 0.0, 0.4], 0.001, 1.0, 1.0).unwrap();
        let mut engine = Engine::new(&geometry.connect(1).unwrap());
        engine.set_frequency(299.8);
        let ground = GroundParams { iperf: 2, nradl: 0, epsr, sig: 0.005, scrwlt: 0.0, scrwrt: 0.0 };
        engine.set_ground(Some(ground)).unwrap();
        let info = engine.prepare_ground();
        (engine, info)
    }

    #[test]
    fn rom2_stops_on_values_that_are_not_finite() {
        let (mut engine, info) = monopole_over_sommerfeld_ground(13.0);
        info.unwrap();
        engine.set_source_segment(0);
        engine.incom.zo = 0.3;
        for (a, b) in [(f64::NAN, 0.05), (-0.05, f64::NAN)] {
            engine.rom2(a, b, 0.0);
            assert_eq!(engine.diag.fault.take(), Some("ROM2"), "limits {a} {b}");
        }
        let sum = engine.rom2(-0.05, 0.05, 0.0);
        assert!(engine.diag.fault.is_none() && all_finite(&sum));
    }

    #[test]
    fn sommerfeld_tables_that_are_not_finite_are_an_error() {
        let (_, info) = monopole_over_sommerfeld_ground(f64::NAN);
        assert_eq!(info, Err(Error::Integration { routine: "ROM1", source: None, observation: None }));
    }

    #[test]
    fn a_ground_without_a_screen_forgets_the_screen_before_it() {
        let (mut engine, _) = monopole_over_sommerfeld_ground(13.0);
        let screen = GroundParams { iperf: 0, nradl: 12, epsr: 13.0, sig: 0.005, scrwlt: 5.0, scrwrt: 0.001 };
        engine.set_ground(Some(screen)).unwrap();
        assert!(engine.prepare_ground().unwrap().radials.is_some());
        assert!(engine.gnd.scrwl > 0.0 && engine.gnd.t2 > 0.0);

        engine.set_ground(Some(GroundParams { iperf: 1, nradl: 0, ..screen })).unwrap();
        engine.prepare_ground().unwrap();
        assert_eq!((engine.gnd.nradl, engine.gnd.scrwl, engine.gnd.scrwr, engine.gnd.t2), (0, 0.0, 0.0, 0.0));
        assert_eq!((engine.gnd.t1, engine.scrwlt, engine.scrwrt), (C0, 0.0, 0.0));
    }
}
