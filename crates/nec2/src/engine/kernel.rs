//! Electric field of a current-carrying segment: the kernel of the
//! electric-field integral equation.
//!
//! `efld` returns the field of the constant, sine and cosine current terms of
//! the source segment held in `self.dj`, including the image in the ground.
//! Close to the source the field is integrated along the segment (`eksc`, or
//! `ekscx` for the extended thin-wire kernel); far away the segment is treated
//! as a point source.
//!
//! Ported from `fields.c` and `calculations.c` of nec2c 1.3.1.

use num_complex::Complex64;

use super::{C0, C1, Diag, Engine, SourceSeg, pass_bound};
use crate::consts::{CONST1, ETA, PCHCON, PI, TP, cmplx};

/// Parameters of the integrand `gf` (COMMON /TMI/).
#[derive(Debug, Clone, Copy)]
pub(crate) struct Tmi {
    pub ij: i32,
    pub zpk: f64,
    pub rkb2: f64,
}

/// Integrand exp(-jkr)/(kr) for `intx`; with `ij == 0` the 1/(kr)
/// singularity of the self term is taken out.
#[inline]
fn gf(tmi: &Tmi, zk: f64) -> (f64, f64) {
    let zdk = zk - tmi.zpk;
    let rk = (tmi.rkb2 + zdk * zdk).sqrt();
    let si = rk.sin() / rk;
    if tmi.ij != 0 {
        return (rk.cos() / rk, si);
    }
    if rk >= 0.2 {
        return ((rk.cos() - 1.0) / rk, si);
    }
    let rks = rk * rk;
    (((-1.38888889e-3 * rks + 4.16666667e-2) * rks - 0.5) * rk, si)
}

/// Relative difference of two estimates of a complex integral, used as the
/// convergence test of the Romberg integrations.
#[inline]
pub(crate) fn test(f1r: f64, f2r: f64, f1i: f64, f2i: f64, dmin: f64) -> (f64, f64) {
    let mut den = f2r.abs();
    let tr = f2i.abs();
    if den < tr {
        den = tr;
    }
    if den < dmin {
        den = dmin;
    }
    if den < 1.0e-37 {
        return (0.0, 0.0);
    }
    (((f1r - f2r) / den).abs(), ((f1i - f2i) / den).abs())
}

/// Integrates exp(-jkr)/(kr) over the source segment from `el1` to `el2`
/// by Romberg integration with a variable interval width. If the integrand
/// is not finite the integration stops and `diag.fault` names this routine.
pub(crate) fn intx(tmi: &Tmi, el1: f64, el2: f64, b: f64, ij: i32, diag: &mut Diag) -> (f64, f64) {
    const NX: i32 = 1;
    const NMA: i32 = 65536;
    const NTS: i32 = 4;
    const RX: f64 = 1.0e-4;

    let mut z = el1;
    let mut ze = el2;
    if ij == 0 {
        ze = 0.0;
    }
    let s = ze - z;
    let fnm = f64::from(NMA);
    let ep = s / (10.0 * fnm);
    let zend = ze - ep;
    let mut sgr = 0.0;
    let mut sgi = 0.0;
    let mut ns = NX;
    let mut nt = 0;
    let (mut g1r, mut g1i) = gf(tmi, z);
    let (mut g3r, mut g3i, mut g5r, mut g5i) = (0.0, 0.0, 0.0, 0.0);
    let mut dz = 0.0;
    let mut dzot = 0.0;
    let mut flag = true;

    // The self term is integrated over half the segment; the singular part
    // is added analytically.
    let finish = |sgr: f64, sgi: f64| {
        if ij == 0 { (2.0 * (sgr + (((b * b + s * s).sqrt() + s) / b).ln()), 2.0 * sgi) } else { (sgr, sgi) }
    };

    for _ in 0..pass_bound(NMA) {
        if flag {
            let fns = f64::from(ns);
            dz = s / fns;
            let zp = z + dz;
            if zp > ze {
                dz = ze - z;
                if dz.abs() <= ep {
                    return finish(sgr, sgi);
                }
            }
            dzot = dz * 0.5;
            (g3r, g3i) = gf(tmi, z + dzot);
            (g5r, g5i) = gf(tmi, z + dz);
        }

        let t00r = (g1r + g5r) * dzot;
        let t00i = (g1i + g5i) * dzot;
        let t01r = (t00r + dz * g3r) * 0.5;
        let t01i = (t00i + dz * g3i) * 0.5;
        let t10r = (4.0 * t01r - t00r) / 3.0;
        let t10i = (4.0 * t01i - t00i) / 3.0;

        // Three-point Romberg result.
        let (te1r, te1i) = test(t01r, t10r, t01i, t10i, 0.0);
        if te1i <= RX && te1r <= RX {
            sgr += t10r;
            sgi += t10i;
            nt += 2;
            z += dz;
            if z >= zend {
                return finish(sgr, sgi);
            }
            g1r = g5r;
            g1i = g5i;
            if nt >= NTS && ns > NX {
                ns /= 2;
                nt = 1;
            }
            flag = true;
            continue;
        }

        let (g2r, g2i) = gf(tmi, z + dz * 0.25);
        let (g4r, g4i) = gf(tmi, z + dz * 0.75);
        let t02r = (t01r + dzot * (g2r + g4r)) * 0.5;
        let t02i = (t01i + dzot * (g2i + g4i)) * 0.5;
        let t11r = (4.0 * t02r - t01r) / 3.0;
        let t11i = (4.0 * t02i - t01i) / 3.0;
        let t20r = (16.0 * t11r - t10r) / 15.0;
        let t20i = (16.0 * t11i - t10i) / 15.0;

        // Five-point Romberg result.
        let (te2r, te2i) = test(t11r, t20r, t11i, t20i, 0.0);
        if te2i > RX || te2r > RX {
            nt = 0;
            if ns >= NMA {
                diag.warnings.push(format!("STEP SIZE LIMITED AT Z= {z:10.5}"));
            } else {
                ns *= 2;
                let fns = f64::from(ns);
                dz = s / fns;
                dzot = dz * 0.5;
                g5r = g3r;
                g5i = g3i;
                g3r = g2r;
                g3i = g2i;
                flag = false;
                continue;
            }
        }

        // A value that is not finite fails both convergence tests and
        // arrives here.
        if !(t20r.is_finite() && t20i.is_finite()) {
            break;
        }
        sgr += t20r;
        sgi += t20i;
        nt += 1;
        z += dz;
        if z >= zend {
            return finish(sgr, sgi);
        }
        g1r = g5r;
        g1i = g5i;
        if nt >= NTS && ns > NX {
            ns /= 2;
            nt = 1;
        }
        flag = true;
    }

    // The integrand is not finite, or the integration does not advance.
    diag.fault = Some("INTX");
    (0.0, 0.0)
}

/// Segment-end contributions for the thin-wire kernel.
#[inline]
fn gx(zz: f64, rh: f64, xk: f64) -> (Complex64, Complex64) {
    let r2 = zz * zz + rh * rh;
    let r = r2.sqrt();
    let rkz = xk * r;
    let gz = cmplx(rkz.cos(), -rkz.sin()) / r;
    let gzp = -cmplx(1.0, rkz) * gz / r2;
    (gz, gzp)
}

/// Segment-end contributions for the extended thin-wire kernel.
/// Returns `(g1, g1p, g2, g2p, g3, gzp)`.
fn gxx(
    zz: f64,
    rh: f64,
    a: f64,
    a2: f64,
    xk: f64,
    ira: i32,
) -> (Complex64, Complex64, Complex64, Complex64, Complex64, Complex64) {
    let r2 = zz * zz + rh * rh;
    let r = r2.sqrt();
    let r4 = r2 * r2;
    let rk = xk * r;
    let rk2 = rk * rk;
    let rh2 = rh * rh;
    let t1 = 0.25 * a2 * rh2 / r4;
    let mut t2 = 0.5 * a2 / r2;
    let c1 = cmplx(1.0, rk);
    let c2 = 3.0 * c1 - rk2;
    let c3 = cmplx(6.0, rk) * rk2 - 15.0 * c1;
    let mut gz = cmplx(rk.cos(), -rk.sin()) / r;
    let mut g2 = gz * (1.0 + t1 * c2);
    let g1 = g2 - t2 * c1 * gz;
    gz /= r2;
    let mut g2p = gz * (t1 * c3 - c1);
    let mut gzp = t2 * c2 * gz;
    let mut g3 = g2p + gzp;
    let g1p = g3 * zz;

    if ira != 1 {
        g3 = (g3 + gzp) * rh;
        gzp = -zz * c1 * gz;
        if rh <= 1.0e-10 {
            return (g1, g1p, C0, C0, g3, gzp);
        }
        g2 /= rh;
        g2p = g2p * zz / rh;
        return (g1, g1p, g2, g2p, g3, gzp);
    }

    t2 = 0.5 * a;
    g2 = -t2 * c1 * gz;
    g2p = t2 * gz * c2 / r2;
    g3 = rh2 * g2p - a * gz * c1;
    g2p *= zz;
    gzp = -zz * c1 * gz;
    (g1, g1p, g2, g2p, g3, gzp)
}

/// Field components along (z) and across (r) the source segment for the sine,
/// cosine and constant current terms: `(ezs, ers, ezc, erc, ezk, erk)`.
type FieldTerms = (Complex64, Complex64, Complex64, Complex64, Complex64, Complex64);

/// Field of a current filament by the thin-wire approximation.
fn eksc(s: f64, z: f64, rh: f64, xk: f64, ij: i32, diag: &mut Diag) -> FieldTerms {
    let rhk = xk * rh;
    let tmi = Tmi { ij, zpk: xk * z, rkb2: rhk * rhk };
    let sh = 0.5 * s;
    let shk = xk * sh;
    let ss = shk.sin();
    let cs = shk.cos();
    let z2a = sh - z;
    let z1a = -(sh + z);
    let (gz1, gp1) = gx(z1a, rh, xk);
    let (gz2, gp2) = gx(z2a, rh, xk);
    let mut gzp1 = gp1 * z1a;
    let mut gzp2 = gp2 * z2a;
    let ezs = CONST1 * ((gz2 - gz1) * cs * xk - (gzp2 + gzp1) * ss);
    let ezc = -CONST1 * ((gz2 + gz1) * ss * xk + (gzp2 - gzp1) * cs);
    let erk = CONST1 * (gp2 - gp1) * rh;
    let (cint, sint) = intx(&tmi, -shk, shk, rhk, ij, diag);
    let ezk = -CONST1 * (gzp2 - gzp1 + xk * xk * cmplx(cint, -sint));
    gzp1 *= z1a;
    gzp2 *= z2a;

    if rh >= 1.0e-10 {
        let ers = -CONST1 * ((gzp2 + gzp1 + gz2 + gz1) * ss - (z2a * gz2 - z1a * gz1) * cs * xk) / rh;
        let erc = -CONST1 * ((gzp2 - gzp1 + gz2 - gz1) * cs + (z2a * gz2 + z1a * gz1) * ss * xk) / rh;
        return (ezs, ers, ezc, erc, ezk, erk);
    }
    (ezs, C0, ezc, C0, ezk, erk)
}

/// Field of a current filament by the extended thin-wire approximation.
#[allow(clippy::too_many_arguments)]
fn ekscx(bx: f64, s: f64, z: f64, rhx: f64, xk: f64, ij: i32, inx1: i32, inx2: i32, diag: &mut Diag) -> FieldTerms {
    let (rh, b, ira) = if rhx >= bx { (rhx, bx, 0) } else { (bx, rhx, 1) };
    let sh = 0.5 * s;
    let rhk = xk * rh;
    let tmi = Tmi { ij, zpk: xk * z, rkb2: rhk * rhk };
    let shk = xk * sh;
    let ss = shk.sin();
    let cs = shk.cos();
    let z2a = sh - z;
    let z1a = -(sh + z);
    let a2 = b * b;

    let end = |za: f64, inx: i32| {
        if inx != 2 {
            gxx(za, rh, b, a2, xk, ira)
        } else {
            let (gz, grk) = gx(za, rhx, xk);
            let gzp = grk * za;
            (gz, gzp, gz / rhx, gzp / rhx, grk * rhx, C0)
        }
    };
    let (gz1, gzp1, gr1, grp1, grk1, gzz1) = end(z1a, inx1);
    let (gz2, gzp2, gr2, grp2, grk2, gzz2) = end(z2a, inx2);

    let ezs = CONST1 * ((gz2 - gz1) * cs * xk - (gzp2 + gzp1) * ss);
    let ezc = -CONST1 * ((gz2 + gz1) * ss * xk + (gzp2 - gzp1) * cs);
    let ers = -CONST1 * ((z2a * grp2 + z1a * grp1 + gr2 + gr1) * ss - (z2a * gr2 - z1a * gr1) * cs * xk);
    let erc = -CONST1 * ((z2a * grp2 - z1a * grp1 + gr2 - gr1) * cs + (z2a * gr2 + z1a * gr1) * ss * xk);
    let erk = CONST1 * (grk2 - grk1);
    let (cint, sint) = intx(&tmi, -shk, shk, rhk, ij, diag);
    let bk = b * xk;
    let bk2 = bk * bk * 0.25;
    let ezk = -CONST1 * (gzp2 - gzp1 + xk * xk * (1.0 - bk2) * cmplx(cint, -sint) - bk2 * (gzz2 - gzz1));
    (ezs, ers, ezc, erc, ezk, erk)
}

impl Engine {
    /// Loads segment `j` (0-based) as the source segment of the field routines.
    pub(crate) fn set_source_segment(&mut self, j: usize) {
        self.dj.s = self.si[j];
        self.dj.b = self.bi[j];
        self.dj.xj = self.x[j];
        self.dj.yj = self.y[j];
        self.dj.zj = self.z[j];
        self.dj.cabj = self.cab[j];
        self.dj.sabj = self.sab[j];
        self.dj.salpj = self.salp[j];
    }

    /// Decides for both ends of segment `j` (0-based) whether the extended
    /// thin-wire kernel may be used: 0 yes, 1 the end is free, 2 no (the end
    /// meets a bend, a change of radius or a junction of several wires).
    ///
    /// An end on a patch gets 2, except end 1 in the matrix fill, where
    /// NEC-2 uses 0: `patch_at_end1` is that value.
    pub(crate) fn set_extended_kernel_ends(&mut self, j: usize, patch_at_end1: i32) {
        let jx = (j + 1) as i32;
        let (cabj, sabj, salpj, b) = (self.dj.cabj, self.dj.sabj, self.dj.salpj, self.dj.b);
        let collinear_same_radius = |k: usize| {
            let xi = (cabj * self.cab[k] + sabj * self.sab[k] + salpj * self.salp[k]).abs();
            !(xi < 0.999999 || (self.bi[k] / b - 1.0).abs() > 1.0e-6)
        };
        // A vertical segment that meets its own image continues straight.
        let own_image = if cabj * cabj + sabj * sabj > 1.0e-8 { 2 } else { 0 };

        let ipr = self.icon1[j];
        self.dj.ind1 = if ipr > PCHCON {
            patch_at_end1
        } else if ipr < 0 {
            let k = (-ipr - 1) as usize;
            if -self.icon1[k] != jx || !collinear_same_radius(k) { 2 } else { 0 }
        } else if ipr == 0 {
            1
        } else if ipr != jx {
            let k = (ipr - 1) as usize;
            if self.icon2[k] != jx || !collinear_same_radius(k) { 2 } else { 0 }
        } else {
            own_image
        };

        let ipr = self.icon2[j];
        self.dj.ind2 = if ipr > PCHCON {
            2
        } else if ipr < 0 {
            let k = (-ipr - 1) as usize;
            if -self.icon2[k] != jx || !collinear_same_radius(k) { 2 } else { 0 }
        } else if ipr == 0 {
            1
        } else if ipr != jx {
            let k = (ipr - 1) as usize;
            if self.icon1[k] != jx || !collinear_same_radius(k) { 2 } else { 0 }
        } else {
            own_image
        };
    }

    /// Computes the near electric field of the source segment at the point
    /// `(xi, yi, zi)`, for its constant, sine and cosine current terms. The
    /// result is left in `self.dj`. `ai` is the radius of the observation
    /// segment and `ij` is zero when the observation point is on the source
    /// segment itself.
    pub(crate) fn efld(&mut self, xi: f64, yi: f64, zi: f64, ai: f64, ij: i32) {
        let SourceSeg { s, b, xj, yj, zj, cabj, sabj, salpj, rkh, iexk, ind1, ind2, .. } = self.dj;
        let xij = xi - xj;
        let yij = yi - yj;
        let mut ijx = ij;
        let mut rfl = -1.0;

        for ip in 0..self.gnd.ksymp {
            if ip == 1 {
                ijx = 1;
            }
            rfl = -rfl;
            let salpr = salpj * rfl;
            let zij = zi - rfl * zj;
            let zp = xij * cabj + yij * sabj + zij * salpr;
            let mut rhox = xij - cabj * zp;
            let mut rhoy = yij - sabj * zp;
            let mut rhoz = zij - salpr * zp;

            let rh = (rhox * rhox + rhoy * rhoy + rhoz * rhoz + ai * ai).sqrt();
            if rh <= 1.0e-10 {
                rhox = 0.0;
                rhoy = 0.0;
                rhoz = 0.0;
            } else {
                rhox /= rh;
                rhoy /= rh;
                rhoz /= rh;
            }

            let r = (zp * zp + rh * rh).sqrt();
            let (tezc, terc, tezk, terk);
            let (mut txs, mut tys, mut tzs);
            if r >= rkh {
                // Lumped current element approximation for large separations.
                let mut rmag = TP * r;
                let cth = zp / r;
                let px = rh / r;
                let txk = cmplx(rmag.cos(), -rmag.sin());
                let py = TP * r * r;
                let tyk = ETA * cth * txk * cmplx(1.0, -1.0 / rmag) / py;
                let tzk = ETA * px * txk * cmplx(1.0, rmag - 1.0 / rmag) / (2.0 * py);
                let ezk = tyk * cth - tzk * px;
                let erk = tyk * px + tzk * cth;
                rmag = (PI * s).sin() / PI;
                tezc = ezk * rmag;
                terc = erk * rmag;
                tezk = ezk * s;
                terk = erk * s;
                txs = C0;
                tys = C0;
                tzs = C0;
            } else {
                let (tezs, ters);
                (tezs, ters, tezc, terc, tezk, terk) = if iexk != 1 {
                    eksc(s, zp, rh, TP, ijx, &mut self.diag)
                } else {
                    ekscx(b, s, zp, rh, TP, ijx, ind1, ind2, &mut self.diag)
                };
                txs = tezs * cabj + ters * rhox;
                tys = tezs * sabj + ters * rhoy;
                tzs = tezs * salpr + ters * rhoz;
            }

            let mut txk = tezk * cabj + terk * rhox;
            let mut tyk = tezk * sabj + terk * rhoy;
            let mut tzk = tezk * salpr + terk * rhoz;
            let mut txc = tezc * cabj + terc * rhox;
            let mut tyc = tezc * sabj + terc * rhoy;
            let mut tzc = tezc * salpr + terc * rhoz;

            if ip == 1 {
                if self.gnd.iperf <= 0 {
                    // Reflection coefficients of the finite ground.
                    let mut zratx = self.gnd.zrati;
                    let rmag = r;
                    let xymag = (xij * xij + yij * yij).sqrt();

                    if self.gnd.nradl > 0 {
                        // Radial-wire ground screen.
                        let xspec = (xi * zj + zi * xj) / (zi + zj);
                        let yspec = (yi * zj + zi * yj) / (zi + zj);
                        let rhospc = (xspec * xspec + yspec * yspec + self.gnd.t2 * self.gnd.t2).sqrt();
                        if rhospc <= self.gnd.scrwl {
                            let zscrn = self.gnd.t1 * rhospc * (rhospc / self.gnd.t2).ln();
                            zratx = (zscrn * self.gnd.zrati) / (ETA * self.gnd.zrati + zscrn);
                        }
                    }

                    let (px, py, cth, zrsin);
                    if xymag <= 1.0e-6 {
                        px = 0.0;
                        py = 0.0;
                        cth = 1.0;
                        zrsin = C1;
                    } else {
                        px = -yij / xymag;
                        py = xij / xymag;
                        cth = zij / rmag;
                        zrsin = (1.0 - zratx * zratx * (1.0 - cth * cth)).sqrt();
                    }

                    let refs = (cth - zratx * zrsin) / (cth + zratx * zrsin);
                    let mut refps = -(zratx * cth - zrsin) / (zratx * cth + zrsin);
                    refps -= refs;
                    let reflect = |tx: &mut Complex64, ty: &mut Complex64, tz: &mut Complex64| {
                        let mut epy = px * *tx + py * *ty;
                        let epx = px * epy;
                        epy = py * epy;
                        *tx = refs * *tx + refps * epx;
                        *ty = refs * *ty + refps * epy;
                        *tz = refs * *tz;
                    };
                    reflect(&mut txk, &mut tyk, &mut tzk);
                    reflect(&mut txs, &mut tys, &mut tzs);
                    reflect(&mut txc, &mut tyc, &mut tzc);
                }

                let frati = self.gnd.frati;
                self.dj.exk -= txk * frati;
                self.dj.eyk -= tyk * frati;
                self.dj.ezk -= tzk * frati;
                self.dj.exs -= txs * frati;
                self.dj.eys -= tys * frati;
                self.dj.ezs -= tzs * frati;
                self.dj.exc -= txc * frati;
                self.dj.eyc -= tyc * frati;
                self.dj.ezc -= tzc * frati;
                continue;
            }

            self.dj.exk = txk;
            self.dj.eyk = tyk;
            self.dj.ezk = tzk;
            self.dj.exs = txs;
            self.dj.eys = tys;
            self.dj.ezs = tzs;
            self.dj.exc = txc;
            self.dj.eyc = tyc;
            self.dj.ezc = tzc;
        }

        if self.gnd.iperf != 2 {
            return;
        }

        // Field due to the ground by the Sommerfeld-Norton method.
        self.incom.sn = (cabj * cabj + sabj * sabj).sqrt();
        if self.incom.sn >= 1.0e-5 {
            self.incom.xsn = cabj / self.incom.sn;
            self.incom.ysn = sabj / self.incom.sn;
        } else {
            self.incom.sn = 0.0;
            self.incom.xsn = 1.0;
            self.incom.ysn = 0.0;
        }

        // Displace the observation point for the thin-wire approximation.
        let zij = zi + zj;
        let salpr = -salpj;
        let rhox = sabj * zij - salpr * yij;
        let rhoy = salpr * xij - cabj * zij;
        let rhoz = cabj * yij - sabj * xij;
        let mut rh = rhox * rhox + rhoy * rhoy + rhoz * rhoz;
        if rh <= 1.0e-10 {
            self.incom.xo = xi - ai * self.incom.ysn;
            self.incom.yo = yi + ai * self.incom.xsn;
            self.incom.zo = zi;
        } else {
            rh = ai / rh.sqrt();
            if rhoz < 0.0 {
                rh = -rh;
            }
            self.incom.xo = xi + rh * rhox;
            self.incom.yo = yi + rh * rhoy;
            self.incom.zo = zi + rh * rhoz;
        }

        // Ground field of the constant, sine and cosine terms: x, y, z of each.
        let mut egnd;
        let r = xij * xij + yij * yij + zij * zij;
        if r <= 0.95 {
            // The field from the interpolation tables is integrated over the segment.
            self.incom.isnor = 1;
            let dmin = 0.01 * (self.dj.exk.norm_sqr() + self.dj.eyk.norm_sqr() + self.dj.ezk.norm_sqr()).sqrt();
            let shaf = 0.5 * s;
            egnd = self.rom2(-shaf, shaf, dmin);
        } else {
            // Norton field equations and lumped current element approximation.
            self.incom.isnor = 2;
            egnd = self.sflds(0.0);

            let zp = xij * cabj + yij * sabj + zij * salpr;
            let rh = r - zp * zp;
            let dmin = if rh <= 1.0e-10 { 0.0 } else { (rh / (rh + ai * ai)).sqrt() };
            if dmin <= 0.95 {
                let px = 1.0 - dmin;
                for term in egnd.chunks_exact_mut(3) {
                    let along = (term[0] * cabj + term[1] * sabj + term[2] * salpr) * px;
                    term[0] = dmin * term[0] + along * cabj;
                    term[1] = dmin * term[1] + along * sabj;
                    term[2] = dmin * term[2] + along * salpr;
                }
            }
        }

        self.dj.exk += egnd[0];
        self.dj.eyk += egnd[1];
        self.dj.ezk += egnd[2];
        self.dj.exs += egnd[3];
        self.dj.eys += egnd[4];
        self.dj.ezs += egnd[5];
        self.dj.exc += egnd[6];
        self.dj.eyc += egnd[7];
        self.dj.ezc += egnd[8];
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn intx_stops_on_values_that_are_not_finite() {
        let tmi = Tmi { ij: 1, zpk: 0.3, rkb2: 0.01 };
        for (el1, el2) in [(f64::NAN, 0.1), (-0.1, f64::NAN), (-0.1, f64::INFINITY)] {
            let mut diag = Diag::default();
            intx(&tmi, el1, el2, 0.1, 1, &mut diag);
            assert_eq!(diag.fault, Some("INTX"), "limits {el1} {el2}");
        }
        // The integrand itself is not a number.
        let mut diag = Diag::default();
        intx(&Tmi { zpk: f64::NAN, ..tmi }, -0.1, 0.1, 0.1, 1, &mut diag);
        assert_eq!(diag.fault, Some("INTX"));

        // A sound integration reports nothing.
        let mut diag = Diag::default();
        let (re, im) = intx(&tmi, -0.1, 0.1, 0.1, 1, &mut diag);
        assert!(diag.fault.is_none() && diag.warnings.is_empty());
        assert!(re.is_finite() && im.is_finite() && re > 0.0);
    }
}
