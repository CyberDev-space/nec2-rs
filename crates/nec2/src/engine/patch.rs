//! Fields of the current on a surface patch.
//!
//! The source patch is the one loaded into `self.dj` by
//! `Engine::set_source_patch`. The results are left in `self.dj`: the field
//! of a unit current along the first tangent vector in `exk`, `eyk`, `ezk`
//! and along the second in `exs`, `eys`, `ezs`.
//!
//! Ported from `hintg`, `unere` and `pcint` (`fields.c`) of nec2c 1.3.1.

use num_complex::Complex64;

use super::{C0, C1, Engine};
use crate::consts::{CONST2, FPI, TP, cmplx};

impl Engine {
    /// Magnetic field at the point `(xi, yi, zi)` of the current on the
    /// source patch, with its image in the ground (`hintg`).
    pub(crate) fn hintg(&mut self, xi: f64, yi: f64, zi: f64) {
        let dj = &mut self.dj;
        let rx = xi - dj.xj;
        let ry = yi - dj.yj;
        let mut rfl = -1.0;
        dj.exk = C0;
        dj.eyk = C0;
        dj.ezk = C0;
        dj.exs = C0;
        dj.eys = C0;
        dj.ezs = C0;

        for ip in 1..=self.gnd.ksymp {
            rfl = -rfl;
            let rz = zi - dj.zj * rfl;
            let rsq = rx * rx + ry * ry + rz * rz;
            if rsq < 1.0e-20 {
                continue;
            }
            let r = rsq.sqrt();
            let rk = TP * r;
            let cr = rk.cos();
            let sr = rk.sin();
            let mut gam = -(cmplx(cr, -sr) + rk * cmplx(sr, cr)) / (FPI * rsq * r) * dj.s;
            dj.exc = gam * rx;
            dj.eyc = gam * ry;
            dj.ezc = gam * rz;
            let t1zr = dj.t1zj * rfl;
            let t2zr = dj.t2zj * rfl;
            let mut f1x = dj.eyc * t1zr - dj.ezc * dj.t1yj;
            let mut f1y = dj.ezc * dj.t1xj - dj.exc * t1zr;
            let mut f1z = dj.exc * dj.t1yj - dj.eyc * dj.t1xj;
            let mut f2x = dj.eyc * t2zr - dj.ezc * dj.t2yj;
            let mut f2y = dj.ezc * dj.t2xj - dj.exc * t2zr;
            let mut f2z = dj.exc * dj.t2yj - dj.eyc * dj.t2xj;

            if ip != 1 {
                if self.gnd.iperf == 1 {
                    f1x = -f1x;
                    f1y = -f1y;
                    f1z = -f1z;
                    f2x = -f2x;
                    f2y = -f2y;
                    f2z = -f2z;
                } else {
                    // Image weighted by the reflection coefficients.
                    let zrati = self.gnd.zrati;
                    let xymag = (rx * rx + ry * ry).sqrt();
                    let (pxx, pyy, cth, mut rrv);
                    if xymag <= 1.0e-6 {
                        pxx = 0.0;
                        pyy = 0.0;
                        cth = 1.0;
                        rrv = C1;
                    } else {
                        pxx = -ry / xymag;
                        pyy = rx / xymag;
                        cth = rz / r;
                        rrv = (1.0 - zrati * zrati * (1.0 - cth * cth)).sqrt();
                    }
                    let mut rrh = zrati * cth;
                    rrh = (rrh - rrv) / (rrh + rrv);
                    rrv = zrati * rrv;
                    rrv = -(cth - rrv) / (cth + rrv);
                    gam = (f1x * pxx + f1y * pyy) * (rrv - rrh);
                    f1x = f1x * rrh + gam * pxx;
                    f1y = f1y * rrh + gam * pyy;
                    f1z *= rrh;
                    gam = (f2x * pxx + f2y * pyy) * (rrv - rrh);
                    f2x = f2x * rrh + gam * pxx;
                    f2y = f2y * rrh + gam * pyy;
                    f2z *= rrh;
                }
            }

            dj.exk += f1x;
            dj.eyk += f1y;
            dj.ezk += f1z;
            dj.exs += f2x;
            dj.eys += f2y;
            dj.ezs += f2z;
        }
    }

    /// Electric field at the point `(xob, yob, zob)` of a unit current on the
    /// source patch, or on its image in the ground if `self.dj.ipgnd` is 2
    /// (`unere`).
    pub(crate) fn unere(&mut self, xob: f64, yob: f64, zob: f64) {
        let dj = &mut self.dj;
        let mut zr = dj.zj;
        let mut t1zr = dj.t1zj;
        let mut t2zr = dj.t2zj;
        if dj.ipgnd == 2 {
            zr = -zr;
            t1zr = -t1zr;
            t2zr = -t2zr;
        }

        let rx = xob - dj.xj;
        let ry = yob - dj.yj;
        let rz = zob - zr;
        let r2 = rx * rx + ry * ry + rz * rz;
        if r2 <= 1.0e-20 {
            dj.exk = C0;
            dj.eyk = C0;
            dj.ezk = C0;
            dj.exs = C0;
            dj.eys = C0;
            dj.ezs = C0;
            return;
        }

        let r = r2.sqrt();
        let tt1 = -TP * r;
        let tt2 = tt1 * tt1;
        let rt = r2 * r;
        let mut er = cmplx(tt1.sin(), -tt1.cos()) * (CONST2 * dj.s);
        let q1 = cmplx(tt2 - 1.0, tt1) * er / rt;
        let q2 = cmplx(3.0 - tt2, -3.0 * tt1) * er / (rt * r2);
        er = q2 * (dj.t1xj * rx + dj.t1yj * ry + t1zr * rz);
        dj.exk = q1 * dj.t1xj + er * rx;
        dj.eyk = q1 * dj.t1yj + er * ry;
        dj.ezk = q1 * t1zr + er * rz;
        er = q2 * (dj.t2xj * rx + dj.t2yj * ry + t2zr * rz);
        dj.exs = q1 * dj.t2xj + er * rx;
        dj.eys = q1 * dj.t2yj + er * ry;
        dj.ezs = q1 * t2zr + er * rz;

        if dj.ipgnd == 1 {
            return;
        }
        if self.gnd.iperf == 1 {
            dj.exk = -dj.exk;
            dj.eyk = -dj.eyk;
            dj.ezk = -dj.ezk;
            dj.exs = -dj.exs;
            dj.eys = -dj.eys;
            dj.ezs = -dj.ezs;
            return;
        }

        let zrati = self.gnd.zrati;
        let xymag = (rx * rx + ry * ry).sqrt();
        let (px, py, cth, mut rrv);
        if xymag <= 1.0e-6 {
            px = 0.0;
            py = 0.0;
            cth = 1.0;
            rrv = C1;
        } else {
            px = -ry / xymag;
            py = rx / xymag;
            cth = rz / (xymag * xymag + rz * rz).sqrt();
            rrv = (1.0 - zrati * zrati * (1.0 - cth * cth)).sqrt();
        }
        let mut rrh = zrati * cth;
        rrh = (rrh - rrv) / (rrh + rrv);
        rrv = zrati * rrv;
        rrv = -(cth - rrv) / (cth + rrv);
        let mut edp = (dj.exk * px + dj.eyk * py) * (rrh - rrv);
        dj.exk = dj.exk * rrv + edp * px;
        dj.eyk = dj.eyk * rrv + edp * py;
        dj.ezk *= rrv;
        edp = (dj.exs * px + dj.eys * py) * (rrh - rrv);
        dj.exs = dj.exs * rrv + edp * px;
        dj.eys = dj.eys * rrv + edp * py;
        dj.ezs *= rrv;
    }

    /// Integrates over the four patches around a wire end the field along a
    /// segment at `(xi, yi, zi)` with the direction cosines `cabi`, `sabi`,
    /// `salpi` (`pcint`). The source patch in `self.dj` is the undivided
    /// patch: its centre is the wire end.
    ///
    /// The first eight results are the matrix elements of the four patches
    /// for the two tangential currents, with the current interpolated between
    /// the patch centres. The ninth belongs to the current that flows out of
    /// the wire end over the surface.
    pub(crate) fn pcint(&mut self, xi: f64, yi: f64, zi: f64, cabi: f64, sabi: f64, salpi: f64) -> [Complex64; 9] {
        const NINT: usize = 10;
        let d = self.dj.s.sqrt() * 0.5;
        let ds = 4.0 * d / NINT as f64;
        let da = ds * ds;
        let gcon = 1.0 / self.dj.s;
        let fcon = 1.0 / (2.0 * TP * d);
        let (xxj, xyj, xzj, xs) = (self.dj.xj, self.dj.yj, self.dj.zj, self.dj.s);
        let (t1xj, t1yj, t1zj) = (self.dj.t1xj, self.dj.t1yj, self.dj.t1zj);
        let (t2xj, t2yj, t2zj) = (self.dj.t2xj, self.dj.t2yj, self.dj.t2zj);
        self.dj.s = da;
        let mut s1 = d + ds * 0.5;
        let mut xss = xxj + s1 * (t1xj + t2xj);
        let mut yss = xyj + s1 * (t1yj + t2yj);
        let mut zss = xzj + s1 * (t1zj + t2zj);
        s1 += d;
        let s2x = s1;
        let mut e = [C0; 9];

        for _ in 0..NINT {
            s1 -= ds;
            let mut s2 = s2x;
            xss -= ds * t1xj;
            yss -= ds * t1yj;
            zss -= ds * t1zj;
            self.dj.xj = xss;
            self.dj.yj = yss;
            self.dj.zj = zss;

            for _ in 0..NINT {
                s2 -= ds;
                self.dj.xj -= ds * t2xj;
                self.dj.yj -= ds * t2yj;
                self.dj.zj -= ds * t2zj;
                self.unere(xi, yi, zi);
                let exk = self.dj.exk * cabi + self.dj.eyk * sabi + self.dj.ezk * salpi;
                let exs = self.dj.exs * cabi + self.dj.eys * sabi + self.dj.ezs * salpi;
                let g1 = (d + s1) * (d + s2) * gcon;
                let g2 = (d - s1) * (d + s2) * gcon;
                let g3 = (d - s1) * (d - s2) * gcon;
                let g4 = (d + s1) * (d - s2) * gcon;
                let mut f2 = (s1 * s1 + s2 * s2) * TP;
                let f1 = s1 / f2 - (g1 - g2 - g3 + g4) * fcon;
                f2 = s2 / f2 - (g1 + g2 - g3 - g4) * fcon;
                e[0] += exk * g1;
                e[1] += exk * g2;
                e[2] += exk * g3;
                e[3] += exk * g4;
                e[4] += exs * g1;
                e[5] += exs * g2;
                e[6] += exs * g3;
                e[7] += exs * g4;
                e[8] += exk * f1 + exs * f2;
            }
        }

        self.dj.xj = xxj;
        self.dj.yj = xyj;
        self.dj.zj = xzj;
        self.dj.s = xs;
        e
    }
}
