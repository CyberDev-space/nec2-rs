//! Excitation: the right-hand side of the matrix equation.
//!
//! Ported from `etmns` (`matrix.c`) and `qdsrc` (`input.c`) of nec2c 1.3.1.

use num_complex::Complex64;

use super::{C0, Engine};
use crate::consts::{CCJ, RETA, TP, cmplx};
use crate::error::Result;

/// What drives the structure.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Excitation {
    /// The voltage sources given to [`Engine::set_sources`] (EX types 0 and 5).
    Sources,
    /// A plane wave arriving from direction (`theta`, `phi`), radians.
    ///
    /// `eta` is the polarization angle of the electric field measured from
    /// the theta unit vector. `polarization` is 1 for linear, 2 for
    /// right-hand elliptic and 3 for left-hand elliptic polarization, with
    /// `axial_ratio` = minor axis / major axis (EX types 1 to 3).
    PlaneWave { polarization: i32, theta: f64, phi: f64, eta: f64, axial_ratio: f64 },
    /// An elementary current source at `position` (wavelengths), oriented by
    /// the angles `alpha` (elevation above the XY plane) and `beta` (azimuth),
    /// radians. `moment` is the current moment divided by the wavelength
    /// squared (EX type 4).
    CurrentElement { position: [f64; 3], alpha: f64, beta: f64, moment: f64 },
}

impl Engine {
    /// Sets the voltage sources: applied-field sources (EX 0) and
    /// current-slope-discontinuity sources (EX 5). Segment numbers are 1-based.
    pub fn set_sources(&mut self, applied_field: &[(usize, Complex64)], current_slope: &[(usize, Complex64)]) {
        self.isant = applied_field.iter().map(|s| s.0).collect();
        self.vsant = applied_field.iter().map(|s| s.1).collect();
        self.ivqd = current_slope.iter().map(|s| s.0).collect();
        self.vqd = current_slope.iter().map(|s| s.1).collect();
    }

    /// Fills the right-hand side with minus the electric field that the
    /// excitation applies along every segment, and with its magnetic field
    /// across the tangent vectors of every patch (`etmns`).
    pub fn excite(&mut self, excitation: &Excitation) -> Result<()> {
        let neq = self.neq();
        let mut e = std::mem::take(&mut self.cur);
        self.iqds.clear();
        self.vqds.clear();
        let result = self.etmns(excitation, &mut e[..neq]);
        self.cur = e;
        result
    }

    fn etmns(&mut self, excitation: &Excitation, e: &mut [Complex64]) -> Result<()> {
        let n = self.n;
        match *excitation {
            Excitation::Sources => {
                e.fill(C0);
                for i in 0..self.isant.len() {
                    let is = self.isant[i] - 1;
                    e[is] = -self.vsant[i] / (self.si[is] * self.wlam);
                }
                for i in 0..self.ivqd.len() {
                    self.qdsrc(self.ivqd[i], self.vqd[i], e)?;
                }
            }

            Excitation::PlaneWave { polarization, theta, phi, eta, axial_ratio } => {
                let (cth, sth) = (theta.cos(), theta.sin());
                let (cph, sph) = (phi.cos(), phi.sin());
                let (cet, set) = (eta.cos(), eta.sin());
                let pxl = cth * cph * cet - sph * set;
                let pyl = cth * sph * cet + cph * set;
                let pzl = -sth * cet;
                let wx = -sth * cph;
                let wy = -sth * sph;
                let wz = -cth;
                let qx = wy * pzl - wz * pyl;
                let qy = wz * pxl - wx * pzl;
                let qz = wx * pyl - wy * pxl;

                // Reflection coefficients of the ground for the reflected wave.
                let (mut rrv, mut rrh) = (C0, C0);
                let ground = self.gnd.ksymp != 1;
                if ground {
                    if self.gnd.iperf != 1 {
                        let zrati = self.gnd.zrati;
                        rrv = (1.0 - zrati * zrati * sth * sth).sqrt();
                        rrh = zrati * cth;
                        rrh = (rrh - rrv) / (rrh + rrv);
                        rrv = zrati * rrv;
                        rrv = -(cth - rrv) / (cth + rrv);
                    } else {
                        rrv = cmplx(-1.0, 0.0);
                        rrh = cmplx(-1.0, 0.0);
                    }
                }

                // Polarization vector: real for linear polarization, complex
                // for elliptic.
                let (mut cx, mut cy, mut cz) = (cmplx(pxl, 0.0), cmplx(pyl, 0.0), cmplx(pzl, 0.0));
                if polarization != 1 {
                    let mut tt1 = -cmplx(0.0, 1.0) * axial_ratio;
                    if polarization == 3 {
                        tt1 = -tt1;
                    }
                    cx = pxl + tt1 * qx;
                    cy = pyl + tt1 * qy;
                    cz = pzl + tt1 * qz;
                }

                for i in 0..n {
                    let arg = -TP * (wx * self.x[i] + wy * self.y[i] + wz * self.z[i]);
                    e[i] = -(cx * self.cab[i] + cy * self.sab[i] + cz * self.salp[i]) * cmplx(arg.cos(), arg.sin());
                }

                if ground {
                    let tt = (cy * cph - cx * sph) * (rrh - rrv);
                    let rx = rrv * cx - tt * sph;
                    let ry = rrv * cy + tt * cph;
                    let rz = -rrv * cz;
                    for i in 0..n {
                        let arg = -TP * (wx * self.x[i] + wy * self.y[i] - wz * self.z[i]);
                        e[i] -= (rx * self.cab[i] + ry * self.sab[i] + rz * self.salp[i]) * cmplx(arg.cos(), arg.sin());
                    }
                }

                // Patches: the magnetic field of the wave, then that of the
                // wave reflected by the ground.
                if self.m != 0 {
                    let (mut cx, mut cy, mut cz) = (cmplx(qx, 0.0), cmplx(qy, 0.0), cmplx(qz, 0.0));
                    if polarization != 1 {
                        let mut tt1 = -cmplx(0.0, 1.0) * axial_ratio;
                        if polarization == 3 {
                            tt1 = -tt1;
                        }
                        cx = qx - tt1 * pxl;
                        cy = qy - tt1 * pyl;
                        cz = qz - tt1 * pzl;
                    }
                    for (is, p) in self.patches.iter().enumerate() {
                        let i1 = n + 2 * is;
                        let arg = -TP * (wx * p.x + wy * p.y + wz * p.z);
                        let tt = cmplx(arg.cos(), arg.sin()) * p.salp * RETA;
                        e[i1 + 1] = (cx * p.t1x + cy * p.t1y + cz * p.t1z) * tt;
                        e[i1] = (cx * p.t2x + cy * p.t2y + cz * p.t2z) * tt;
                    }
                    if ground {
                        let tt1 = (cy * cph - cx * sph) * (rrv - rrh);
                        let rx = -(rrh * cx - tt1 * sph);
                        let ry = -(rrh * cy + tt1 * cph);
                        let rz = rrh * cz;
                        for (is, p) in self.patches.iter().enumerate() {
                            let i1 = n + 2 * is;
                            let arg = -TP * (wx * p.x + wy * p.y - wz * p.z);
                            let tt = cmplx(arg.cos(), arg.sin()) * p.salp * RETA;
                            e[i1 + 1] += (rx * p.t1x + ry * p.t1y + rz * p.t1z) * tt;
                            e[i1] += (rx * p.t2x + ry * p.t2y + rz * p.t2z) * tt;
                        }
                    }
                }
            }

            Excitation::CurrentElement { position, alpha, beta, moment } => {
                let [p1, p2, p3] = position;
                let mut wz = alpha.cos();
                let wx = wz * beta.cos();
                let wy = wz * beta.sin();
                wz = alpha.sin();
                let ds = moment * 59.958;

                for i in 0..n {
                    let mut pxl = self.x[i] - p1;
                    let mut pyl = self.y[i] - p2;
                    let mut pzl = self.z[i] - p3;
                    let rs = pxl * pxl + pyl * pyl + pzl * pzl;
                    if rs < 1.0e-30 {
                        // The source sits on the segment centre: no field is
                        // computed there, and no earlier value may survive.
                        e[i] = C0;
                        continue;
                    }
                    let r = rs.sqrt();
                    pxl /= r;
                    pyl /= r;
                    pzl /= r;
                    let cth = pxl * wx + pyl * wy + pzl * wz;
                    let sth = (1.0 - cth * cth).sqrt();
                    let mut qx = pxl - wx * cth;
                    let mut qy = pyl - wy * cth;
                    let mut qz = pzl - wz * cth;
                    let arg = (qx * qx + qy * qy + qz * qz).sqrt();
                    if arg >= 1.0e-30 {
                        qx /= arg;
                        qy /= arg;
                        qz /= arg;
                    } else {
                        qx = 1.0;
                        qy = 0.0;
                        qz = 0.0;
                    }

                    let arg = -TP * r;
                    let tt1 = cmplx(arg.cos(), arg.sin());
                    let tt2 = cmplx(1.0, -1.0 / (r * TP)) / rs;
                    let er = ds * tt1 * tt2 * cth;
                    let et = 0.5 * ds * tt1 * (cmplx(0.0, 1.0) * TP / r + tt2) * sth;
                    let ezh = er * cth - et * sth;
                    let erh = er * sth + et * cth;
                    let cx = ezh * wx + erh * qx;
                    let cy = ezh * wy + erh * qy;
                    let cz = ezh * wz + erh * qz;
                    e[i] = -(cx * self.cab[i] + cy * self.sab[i] + cz * self.salp[i]);
                }

                // Patches: the magnetic field of the element.
                let dsh = moment / (2.0 * TP);
                for (is, p) in self.patches.iter().enumerate() {
                    let i1 = n + 2 * is;
                    let mut pxl = p.x - p1;
                    let mut pyl = p.y - p2;
                    let mut pzl = p.z - p3;
                    let rs = pxl * pxl + pyl * pyl + pzl * pzl;
                    if rs < 1.0e-30 {
                        // The source sits on the patch centre. (nec2c also
                        // loses count of the patches here.)
                        e[i1] = C0;
                        e[i1 + 1] = C0;
                        continue;
                    }
                    let r = rs.sqrt();
                    pxl /= r;
                    pyl /= r;
                    pzl /= r;
                    let cth = pxl * wx + pyl * wy + pzl * wz;
                    let sth = (1.0 - cth * cth).sqrt();
                    let mut qx = pxl - wx * cth;
                    let mut qy = pyl - wy * cth;
                    let mut qz = pzl - wz * cth;
                    let arg = (qx * qx + qy * qy + qz * qz).sqrt();
                    if arg >= 1.0e-30 {
                        qx /= arg;
                        qy /= arg;
                        qz /= arg;
                    } else {
                        qx = 1.0;
                        qy = 0.0;
                        qz = 0.0;
                    }
                    let arg = -TP * r;
                    let tt1 = cmplx(arg.cos(), arg.sin());
                    let hx = wy * qz - wz * qy;
                    let hy = wz * qx - wx * qz;
                    let hz = wx * qy - wy * qx;
                    let tt2 = dsh * tt1 * cmplx(1.0 / r, TP) / r * sth * p.salp;
                    let (cx, cy, cz) = (tt2 * hx, tt2 * hy, tt2 * hz);
                    e[i1 + 1] = cx * p.t1x + cy * p.t1y + cz * p.t1z;
                    e[i1] = cx * p.t2x + cy * p.t2y + cz * p.t2z;
                }
            }
        }
        Ok(())
    }

    /// Adds to `e` the field of a current-slope-discontinuity source of
    /// voltage `v` on segment `is` (1-based).
    fn qdsrc(&mut self, is: usize, v: Complex64, e: &mut [Complex64]) -> Result<()> {
        let isx = is - 1;
        let saved = self.icon1[isx];
        self.icon1[isx] = 0;
        let built = self.tbf(is as i32, 0);
        self.icon1[isx] = saved;
        built?;

        let sh = self.si[isx] * 0.5;
        let last = self.segj.jco.len() - 1;
        let curd = CCJ * v
            / (((2.0 * sh / self.bi[isx]).ln() - 1.0)
                * (self.segj.bx[last] * (TP * sh).cos() + self.segj.cx[last] * (TP * sh).sin())
                * self.wlam);
        self.vqds.push(v);
        self.iqds.push(is);

        for jx in 0..self.segj.jco.len() {
            let j = (self.segj.jco[jx] - 1) as usize;
            let (ax, bx, cx) = (self.segj.ax[jx], self.segj.bx[jx], self.segj.cx[jx]);
            self.set_source_segment(j);
            if self.dj.iexk != 0 {
                self.set_extended_kernel_ends(j, 2);
            }

            for i in 0..self.n {
                let ij = i as i32 - j as i32;
                self.efld(self.x[i], self.y[i], self.z[i], self.bi[i], ij);
                self.check_integration(j + 1, Some(i + 1))?;
                let (cabi, sabi, salpi) = (self.cab[i], self.sab[i], self.salp[i]);
                let etk = self.dj.exk * cabi + self.dj.eyk * sabi + self.dj.ezk * salpi;
                let ets = self.dj.exs * cabi + self.dj.eys * sabi + self.dj.ezs * salpi;
                let etc = self.dj.exc * cabi + self.dj.eyc * sabi + self.dj.ezc * salpi;
                e[i] -= (etk * ax + ets * bx + etc * cx) * curd;
            }

            // The magnetic field of the source at the patches.
            for is in 0..self.m {
                let p = self.patches[is];
                self.hsfld(p.x, p.y, p.z, 0.0);
                self.check_integration(j + 1, None)?;
                let dj = &self.dj;
                for (k, (tx, ty, tz)) in [(p.t2x, p.t2y, p.t2z), (p.t1x, p.t1y, p.t1z)].into_iter().enumerate() {
                    let etk = dj.exk * tx + dj.eyk * ty + dj.ezk * tz;
                    let ets = dj.exs * tx + dj.eys * ty + dj.ezs * tz;
                    let etc = dj.exc * tx + dj.eyc * ty + dj.ezc * tz;
                    e[self.n + 2 * is + k] += (etk * ax + ets * bx + etc * cx) * curd * p.salp;
                }
            }

            if self.nload > 0 {
                e[j] += self.zarray[j] * curd * (ax + cx);
            }
        }
        Ok(())
    }
}
