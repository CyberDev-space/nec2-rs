//! Near electric and magnetic fields of the solved structure.
//!
//! Ported from `nefld`, `nhfld`, `hsfld`, `hsflx`, `hfk` and `gh`
//! (`fields.c`) of nec2c 1.3.1.

use num_complex::Complex64;

use super::kernel::test;
use super::{C0, C1, Diag, Engine, pass_bound};
use crate::consts::{ETA, PI8, TP, TPJ, cmplx};
use crate::error::Result;

/// Integrand of the H field of a uniform current filament.
#[inline]
fn gh(zk: f64, zpka: f64, rhks: f64) -> (f64, f64) {
    let rs = zk - zpka;
    let rs = rhks + rs * rs;
    let r = rs.sqrt();
    let ckr = r.cos();
    let skr = r.sin();
    let rr2 = 1.0 / rs;
    let rr3 = rr2 / r;
    (skr * rr2 + ckr * rr3, ckr * rr2 - skr * rr3)
}

/// H field of a uniform current filament by variable-interval Romberg
/// integration from `el1` to `el2` (all arguments in radians of phase). If
/// the integrand is not finite the integration stops and `diag.fault` names
/// this routine.
fn hfk(el1: f64, el2: f64, rhk: f64, zpkx: f64, diag: &mut Diag) -> (f64, f64) {
    const NX: i32 = 1;
    const NMA: i32 = 65536;
    const NTS: i32 = 4;
    const RX: f64 = 1.0e-4;

    let rhks = rhk * rhk;
    let mut z = el1;
    let ze = el2;
    let s = ze - z;
    let ep = s / (10.0 * f64::from(NMA));
    let zend = ze - ep;
    let mut sgr = 0.0;
    let mut sgi = 0.0;
    let mut ns = NX;
    let mut nt = 0;
    let (mut g1r, mut g1i) = gh(z, zpkx, rhks);
    let (mut g3r, mut g3i, mut g5r, mut g5i) = (0.0, 0.0, 0.0, 0.0);
    let mut dz = 0.0;
    let mut dzot = 0.0;
    let mut flag = true;
    let finish = |sgr: f64, sgi: f64| (sgr * rhk * 0.5, sgi * rhk * 0.5);

    for _ in 0..pass_bound(NMA) {
        if flag {
            dz = s / f64::from(ns);
            if z + dz > ze {
                dz = ze - z;
                if dz.abs() <= ep {
                    return finish(sgr, sgi);
                }
            }
            dzot = dz * 0.5;
            (g3r, g3i) = gh(z + dzot, zpkx, rhks);
            (g5r, g5i) = gh(z + dz, zpkx, rhks);
        }

        let t00r = (g1r + g5r) * dzot;
        let t00i = (g1i + g5i) * dzot;
        let t01r = (t00r + dz * g3r) * 0.5;
        let t01i = (t00i + dz * g3i) * 0.5;
        let t10r = (4.0 * t01r - t00r) / 3.0;
        let t10i = (4.0 * t01i - t00i) / 3.0;

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

        let (g2r, g2i) = gh(z + dz * 0.25, zpkx, rhks);
        let (g4r, g4i) = gh(z + dz * 0.75, zpkx, rhks);
        let t02r = (t01r + dzot * (g2r + g4r)) * 0.5;
        let t02i = (t01i + dzot * (g2i + g4i)) * 0.5;
        let t11r = (4.0 * t02r - t01r) / 3.0;
        let t11i = (4.0 * t02i - t01i) / 3.0;
        let t20r = (16.0 * t11r - t10r) / 15.0;
        let t20i = (16.0 * t11i - t10i) / 15.0;

        let (te2r, te2i) = test(t11r, t20r, t11i, t20i, 0.0);
        if te2i > RX || te2r > RX {
            nt = 0;
            if ns >= NMA {
                diag.warnings.push(format!("STEP SIZE LIMITED AT Z= {z:10.5}"));
            } else {
                ns *= 2;
                dz = s / f64::from(ns);
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
    diag.fault = Some("HFK");
    (0.0, 0.0)
}

/// H field of the constant, sine and cosine current terms of a segment of
/// length `s` at distance `rh` from its axis and `zpx` along it:
/// `(hpk, hps, hpc)`.
fn hsflx(s: f64, rh: f64, zpx: f64, diag: &mut Diag) -> (Complex64, Complex64, Complex64) {
    if rh < 1.0e-10 {
        return (C0, C0, C0);
    }
    let fjk = -TPJ;
    let j = cmplx(0.0, 1.0);
    let (zp, hss) = if zpx >= 0.0 { (zpx, 1.0) } else { (-zpx, -1.0) };
    let dh = 0.5 * s;
    let z1 = zp + dh;
    let z2a = zp - dh;
    let rhz = if z2a >= 1.0e-7 { rh / z2a } else { 1.0 };

    let dk = TP * dh;
    let cdk = dk.cos();
    let sdk = dk.sin();
    let (hkr, hki) = hfk(-dk, dk, rh * TP, zp * TP, diag);
    let hpk = cmplx(hkr, hki);

    if rhz >= 1.0e-3 {
        let rh2 = rh * rh;
        let r1 = (rh2 + z1 * z1).sqrt();
        let r2 = (rh2 + z2a * z2a).sqrt();
        let ekr1 = (fjk * r1).exp();
        let ekr2 = (fjk * r2).exp();
        let t1 = z1 * ekr1 / r1;
        let t2 = z2a * ekr2 / r2;
        let hps = (cdk * (ekr2 - ekr1) - j * sdk * (t2 + t1)) * hss;
        let hpc = -sdk * (ekr2 + ekr1) - j * cdk * (t2 - t1);
        let cons = -j / (2.0 * TP * rh);
        return (hpk, cons * hps, cons * hpc);
    }

    // Far along the axis of the segment.
    let ekr1 = cmplx(cdk, sdk) / (z2a * z2a);
    let ekr2 = cmplx(cdk, -sdk) / (z1 * z1);
    let t1 = TP * (1.0 / z1 - 1.0 / z2a);
    let t2 = (fjk * zp).exp() * rh / PI8;
    let hps = t2 * (t1 + (ekr1 + ekr2) * sdk) * hss;
    let hpc = t2 * (-j * t1 + (ekr1 - ekr2) * cdk);
    (hpk, hps, hpc)
}

impl Engine {
    /// Radius to keep the observation point away from the axis of a segment
    /// that it lies inside of, or zero.
    fn observation_radius(&self, xob: f64, yob: f64, zob: f64) -> f64 {
        for i in 0..self.n {
            let xj = xob - self.x[i];
            let yj = yob - self.y[i];
            let zj = zob - self.z[i];
            let zp = self.cab[i] * xj + self.sab[i] * yj + self.salp[i] * zj;
            if zp.abs() > 0.5001 * self.si[i] {
                continue;
            }
            let rho2 = xj * xj + yj * yj + zj * zj - zp * zp;
            let radius = self.bi[i];
            if rho2 > 0.9 * radius * radius {
                continue;
            }
            return radius;
        }
        0.0
    }

    /// The current on patch `i` (0-based) along its two tangent vectors.
    fn patch_current(&self, i: usize) -> (Complex64, Complex64) {
        let p = &self.patches[i];
        let jc = self.n + 3 * i;
        let (cx, cy, cz) = (self.cur[jc], self.cur[jc + 1], self.cur[jc + 2]);
        (p.t1x * cx + p.t1y * cy + p.t1z * cz, p.t2x * cx + p.t2y * cy + p.t2z * cz)
    }

    /// Near electric field at a point given in wavelengths, in volts per
    /// metre, for the last solution (`nefld`).
    pub fn near_electric_field(&mut self, xob: f64, yob: f64, zob: f64) -> Result<[Complex64; 3]> {
        let ax = self.observation_radius(xob, yob, zob);
        let (mut ex, mut ey, mut ez) = (C0, C0, C0);
        for i in 0..self.n {
            self.set_source_segment(i);
            if self.dj.iexk != 0 {
                self.set_extended_kernel_ends(i, 2);
            }
            self.efld(xob, yob, zob, ax, 1);
            self.check_integration(i + 1, None)?;
            let acx = cmplx(self.air[i], self.aii[i]);
            let bcx = cmplx(self.bir[i], self.bii[i]);
            let ccx = cmplx(self.cir[i], self.cii[i]);
            ex += self.dj.exk * acx + self.dj.exs * bcx + self.dj.exc * ccx;
            ey += self.dj.eyk * acx + self.dj.eys * bcx + self.dj.eyc * ccx;
            ez += self.dj.ezk * acx + self.dj.ezs * bcx + self.dj.ezc * ccx;
        }

        // The patches and their images in the ground.
        for i in 0..self.m {
            self.set_source_patch(i);
            let (acx, bcx) = self.patch_current(i);
            for ip in 1..=self.gnd.ksymp as i32 {
                self.dj.ipgnd = ip;
                self.unere(xob, yob, zob);
                ex += acx * self.dj.exk + bcx * self.dj.exs;
                ey += acx * self.dj.eyk + bcx * self.dj.eys;
                ez += acx * self.dj.ezk + bcx * self.dj.ezs;
            }
        }
        Ok([ex, ey, ez])
    }

    /// Near magnetic field at a point given in wavelengths, in amperes per
    /// metre, for the last solution (`nhfld`).
    pub fn near_magnetic_field(&mut self, xob: f64, yob: f64, zob: f64) -> Result<[Complex64; 3]> {
        let ax = self.observation_radius(xob, yob, zob);
        let (mut hx, mut hy, mut hz) = (C0, C0, C0);
        for i in 0..self.n {
            self.set_source_segment(i);
            self.hsfld(xob, yob, zob, ax);
            self.check_integration(i + 1, None)?;
            let acx = cmplx(self.air[i], self.aii[i]);
            let bcx = cmplx(self.bir[i], self.bii[i]);
            let ccx = cmplx(self.cir[i], self.cii[i]);
            hx += self.dj.exk * acx + self.dj.exs * bcx + self.dj.exc * ccx;
            hy += self.dj.eyk * acx + self.dj.eys * bcx + self.dj.eyc * ccx;
            hz += self.dj.ezk * acx + self.dj.ezs * bcx + self.dj.ezc * ccx;
        }

        for i in 0..self.m {
            self.set_source_patch(i);
            self.hintg(xob, yob, zob);
            let (acx, bcx) = self.patch_current(i);
            hx += acx * self.dj.exk + bcx * self.dj.exs;
            hy += acx * self.dj.eyk + bcx * self.dj.eys;
            hz += acx * self.dj.ezk + bcx * self.dj.ezs;
        }
        Ok([hx, hy, hz])
    }

    /// H field of the constant, sine and cosine current terms of the source
    /// segment at the point `(xi, yi, zi)`, including the image in the
    /// ground. The result is left in `self.dj`.
    pub(crate) fn hsfld(&mut self, xi: f64, yi: f64, zi: f64, ai: f64) {
        let (s, xj, yj, zj) = (self.dj.s, self.dj.xj, self.dj.yj, self.dj.zj);
        let (cabj, sabj, salpj) = (self.dj.cabj, self.dj.sabj, self.dj.salpj);
        let xij = xi - xj;
        let yij = yi - yj;
        let mut rfl = -1.0;

        for ip in 0..self.gnd.ksymp {
            rfl = -rfl;
            let salpr = salpj * rfl;
            let zij = zi - rfl * zj;
            let zp = xij * cabj + yij * sabj + zij * salpr;
            let mut rhox = xij - cabj * zp;
            let mut rhoy = yij - sabj * zp;
            let mut rhoz = zij - salpr * zp;
            let rh = (rhox * rhox + rhoy * rhoy + rhoz * rhoz + ai * ai).sqrt();

            if rh <= 1.0e-10 {
                self.dj.exk = C0;
                self.dj.eyk = C0;
                self.dj.ezk = C0;
                self.dj.exs = C0;
                self.dj.eys = C0;
                self.dj.ezs = C0;
                self.dj.exc = C0;
                self.dj.eyc = C0;
                self.dj.ezc = C0;
                continue;
            }

            rhox /= rh;
            rhoy /= rh;
            rhoz /= rh;
            let phx = sabj * rhoz - salpr * rhoy;
            let phy = salpr * rhox - cabj * rhoz;
            let phz = cabj * rhoy - sabj * rhox;
            let (hpk, hps, hpc) = hsflx(s, rh, zp, &mut self.diag);

            if ip == 1 {
                // Direction of the image field; for a finite ground it is
                // weighted by the reflection coefficients.
                let (qx, qy, qz);
                if self.gnd.iperf != 1 {
                    let mut zratx = self.gnd.zrati;
                    let rmag = (zp * zp + rh * rh).sqrt();
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

                    let (px, py, cth, mut rrv);
                    if xymag <= 1.0e-6 {
                        px = 0.0;
                        py = 0.0;
                        cth = 1.0;
                        rrv = C1;
                    } else {
                        px = -yij / xymag;
                        py = xij / xymag;
                        cth = zij / rmag;
                        rrv = (1.0 - zratx * zratx * (1.0 - cth * cth)).sqrt();
                    }
                    let mut rrh = zratx * cth;
                    rrh = -(rrh - rrv) / (rrh + rrv);
                    rrv = zratx * rrv;
                    rrv = (cth - rrv) / (cth + rrv);
                    let q = (phx * px + phy * py) * (rrv - rrh);
                    qx = q * px + phx * rrh;
                    qy = q * py + phy * rrh;
                    qz = phz * rrh;
                } else {
                    qx = cmplx(phx, 0.0);
                    qy = cmplx(phy, 0.0);
                    qz = cmplx(phz, 0.0);
                }
                self.dj.exk -= hpk * qx;
                self.dj.eyk -= hpk * qy;
                self.dj.ezk -= hpk * qz;
                self.dj.exs -= hps * qx;
                self.dj.eys -= hps * qy;
                self.dj.ezs -= hps * qz;
                self.dj.exc -= hpc * qx;
                self.dj.eyc -= hpc * qy;
                self.dj.ezc -= hpc * qz;
                continue;
            }

            self.dj.exk = hpk * phx;
            self.dj.eyk = hpk * phy;
            self.dj.ezk = hpk * phz;
            self.dj.exs = hps * phx;
            self.dj.eys = hps * phy;
            self.dj.ezs = hps * phz;
            self.dj.exc = hpc * phx;
            self.dj.eyc = hpc * phy;
            self.dj.ezc = hpc * phz;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hfk_stops_on_values_that_are_not_finite() {
        for (el1, el2) in [(f64::NAN, 0.1), (-0.1, f64::NAN), (-0.1, f64::INFINITY)] {
            let mut diag = Diag::default();
            hfk(el1, el2, 0.1, 0.3, &mut diag);
            assert_eq!(diag.fault, Some("HFK"), "limits {el1} {el2}");
        }
        let mut diag = Diag::default();
        let (re, im) = hfk(-0.1, 0.1, 0.1, 0.3, &mut diag);
        assert!(diag.fault.is_none() && re.is_finite() && im.is_finite());
    }
}
