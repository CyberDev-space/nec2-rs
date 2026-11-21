//! Current expansion: the NEC-2 basis functions.
//!
//! The current on every segment is a constant plus a sine and a cosine term.
//! A basis function is centred on one segment and extends over the segments
//! connected to its two ends; `tbf` builds one basis function, `sbf` returns
//! the part of a basis function that lies on a given segment, `trio` collects
//! the parts of all basis functions on one segment and `cabc` turns the
//! solution vector into the three coefficients of every segment.
//!
//! Ported from `calculations.c` of nec2c 1.3.1.

use num_complex::Complex64;

use super::Engine;
use crate::consts::{CCJ, PCHCON, PI, TP, cmplx};
use crate::error::{Error, Result};

/// 1 - cos(2 d), with a series for small arguments.
#[inline]
fn one_minus_cos_2d(d: f64, cdh: f64, sdh: f64) -> f64 {
    if d <= 0.015 {
        let omc = 4.0 * d * d;
        ((1.3888889e-3 * omc - 4.1666666667e-2) * omc + 0.5) * omc
    } else {
        1.0 - cdh * cdh + sdh * sdh
    }
}

/// End-cap term of a free wire end of radius `b` (in wavelengths).
#[inline]
fn end_cap(b: f64) -> f64 {
    let qp = PI * b;
    let xxi = qp * qp;
    qp * (1.0 - 0.5 * xxi) / (1.0 - xxi)
}

impl Engine {
    /// Computes basis function `i` (1-based segment number) into `self.segj`.
    ///
    /// `icap == 0` leaves out the end-cap term at free ends.
    pub(crate) fn tbf(&mut self, i: i32, icap: i32) -> Result<()> {
        self.segj.clear();
        let mut pp = 0.0;
        let ix = (i - 1) as usize;
        // An end on a patch is treated like an end on its own image.
        let mut jcox = if self.icon1[ix] > PCHCON { i } else { self.icon1[ix] };
        let mut jend = -1;
        let mut iend = -1;
        let mut sig = -1.0;
        let mut pm = 0.0;
        let mut njun1 = 0;

        // First the segments connected to end 1, then those connected to end 2.
        loop {
            if jcox != 0 {
                if jcox < 0 {
                    jcox = -jcox;
                } else {
                    sig = -sig;
                    jend = -jend;
                }
                let jcoxx = (jcox - 1) as usize;
                let d = PI * self.si[jcoxx];
                let sdh = d.sin();
                let cdh = d.cos();
                let sd = 2.0 * sdh * cdh;
                let omc = one_minus_cos_2d(d, cdh, sdh);
                let aj = 1.0 / ((1.0 / (PI * self.bi[jcoxx])).ln() - 0.577215664);
                pp -= omc / sd * aj;
                self.segj.push(jcox, aj / sd * sig, aj / (2.0 * cdh), -aj / (2.0 * sdh) * sig);

                if jcox != i {
                    jcox = if jend == 1 { self.icon2[jcoxx] } else { self.icon1[jcoxx] };
                    if jcox.abs() != i {
                        if jcox != 0 && jcox <= PCHCON {
                            continue;
                        }
                        return Err(Error::Connection { routine: "TBF", segment: i as usize });
                    }
                } else {
                    // The segment is connected to its own image in the ground.
                    let last = self.segj.bx.len() - 1;
                    self.segj.bx[last] = -self.segj.bx[last];
                }
                if iend == 1 {
                    break;
                }
            }

            pm = -pp;
            pp = 0.0;
            njun1 = self.segj.jco.len();
            jcox = if self.icon2[ix] > PCHCON { i } else { self.icon2[ix] };
            jend = 1;
            iend = 1;
            sig = -1.0;
            if jcox == 0 {
                break;
            }
        }

        let jsno = self.segj.jco.len();
        let njun2 = jsno - njun1;
        let d = PI * self.si[ix];
        let sdh = d.sin();
        let cdh = d.cos();
        let sd = 2.0 * sdh * cdh;
        let cd = cdh * cdh - sdh * sdh;
        let omc = if d <= 0.015 {
            let omc = 4.0 * d * d;
            ((1.3888889e-3 * omc - 4.1666666667e-2) * omc + 0.5) * omc
        } else {
            1.0 - cd
        };
        let ap = 1.0 / ((1.0 / (PI * self.bi[ix])).ln() - 0.577215664);
        let aj = ap;
        let xxi = if icap == 0 { 0.0 } else { end_cap(self.bi[ix]) };
        let segj = &mut self.segj;

        if njun1 == 0 {
            if njun2 == 0 {
                segj.push(i, -1.0, 0.0, 1.0 / (cdh - xxi * sdh));
                return Ok(());
            }
            let qp = -(omc + xxi * sd) / (sd * (ap + xxi * pp) + cd * (xxi * ap - pp));
            let d = cd - xxi * sd;
            let bxp = (sdh + ap * qp * (cdh - xxi * sdh)) / d;
            let cxp = (cdh + ap * qp * (sdh + xxi * cdh)) / d;
            for k in 0..njun2 {
                segj.ax[k] = -segj.ax[k] * qp;
                segj.bx[k] *= qp;
                segj.cx[k] = -segj.cx[k] * qp;
            }
            segj.push(i, -1.0, bxp, cxp);
            return Ok(());
        }

        if njun2 == 0 {
            let qm = (omc + xxi * sd) / (sd * (aj - xxi * pm) + cd * (pm + xxi * aj));
            let d = cd - xxi * sd;
            let bxp = (aj * qm * (cdh - xxi * sdh) - sdh) / d;
            let cxp = (cdh - aj * qm * (sdh + xxi * cdh)) / d;
            for k in 0..njun1 {
                segj.ax[k] *= qm;
                segj.bx[k] *= qm;
                segj.cx[k] *= qm;
            }
            segj.push(i, -1.0, bxp, cxp);
            return Ok(());
        }

        let mut qp = sd * (pm * pp + aj * ap) + cd * (pm * ap - pp * aj);
        let qm = (ap * omc - pp * sd) / qp;
        qp = -(aj * omc + pm * sd) / qp;
        let bxp = (aj * qm + ap * qp) * sdh / sd;
        let cxp = (aj * qm - ap * qp) * cdh / sd;
        for k in 0..njun1 {
            segj.ax[k] *= qm;
            segj.bx[k] *= qm;
            segj.cx[k] *= qm;
        }
        for k in njun1..jsno {
            segj.ax[k] = -segj.ax[k] * qp;
            segj.bx[k] *= qp;
            segj.cx[k] = -segj.cx[k] * qp;
        }
        segj.push(i, -1.0, bxp, cxp);
        Ok(())
    }

    /// Returns the constant, sine and cosine coefficients of basis function
    /// `i` on segment `is` (both 1-based).
    pub(crate) fn sbf(&self, i: i32, is: i32) -> Result<(f64, f64, f64)> {
        let (mut aa, mut bb, mut cc) = (0.0, 0.0, 0.0);
        let mut june = 0;
        let mut jsno = 0;
        let mut pp = 0.0;
        let ix = (i - 1) as usize;
        let mut jcox = if self.icon1[ix] > PCHCON { i } else { self.icon1[ix] };
        let mut jend = -1;
        let mut iend = -1;
        let mut sig = -1.0;
        let mut pm = 0.0;
        let mut njun1 = 0;

        loop {
            if jcox != 0 {
                if jcox < 0 {
                    jcox = -jcox;
                } else {
                    sig = -sig;
                    jend = -jend;
                }
                let jcoxx = (jcox - 1) as usize;
                jsno += 1;
                let d = PI * self.si[jcoxx];
                let sdh = d.sin();
                let cdh = d.cos();
                let sd = 2.0 * sdh * cdh;
                let omc = one_minus_cos_2d(d, cdh, sdh);
                let aj = 1.0 / ((1.0 / (PI * self.bi[jcoxx])).ln() - 0.577215664);
                pp -= omc / sd * aj;

                if jcox == is {
                    aa = aj / sd * sig;
                    bb = aj / (2.0 * cdh);
                    cc = -aj / (2.0 * sdh) * sig;
                    june = iend;
                }

                if jcox != i {
                    jcox = if jend != 1 { self.icon1[jcoxx] } else { self.icon2[jcoxx] };
                    if jcox.abs() != i {
                        if jcox == 0 || jcox > PCHCON {
                            return Err(Error::Connection { routine: "SBF", segment: i as usize });
                        }
                        continue;
                    }
                } else if jcox == is {
                    bb = -bb;
                }
                if iend == 1 {
                    break;
                }
            }

            pm = -pp;
            pp = 0.0;
            njun1 = jsno;
            jcox = if self.icon2[ix] > PCHCON { i } else { self.icon2[ix] };
            jend = 1;
            iend = 1;
            sig = -1.0;
            if jcox == 0 {
                break;
            }
        }

        let njun2 = jsno - njun1;
        let d = PI * self.si[ix];
        let sdh = d.sin();
        let cdh = d.cos();
        let sd = 2.0 * sdh * cdh;
        let cd = cdh * cdh - sdh * sdh;
        let omc = if d <= 0.015 {
            let omc = 4.0 * d * d;
            ((1.3888889e-3 * omc - 4.1666666667e-2) * omc + 0.5) * omc
        } else {
            1.0 - cd
        };
        let ap = 1.0 / ((1.0 / (PI * self.bi[ix])).ln() - 0.577215664);
        let aj = ap;

        if njun1 == 0 {
            let xxi = end_cap(self.bi[ix]);
            if njun2 == 0 {
                aa = -1.0;
                cc = 1.0 / (cdh - xxi * sdh);
                return Ok((aa, bb, cc));
            }
            let qp = -(omc + xxi * sd) / (sd * (ap + xxi * pp) + cd * (xxi * ap - pp));
            if june == 1 {
                aa = -aa * qp;
                bb *= qp;
                cc = -cc * qp;
                if i != is {
                    return Ok((aa, bb, cc));
                }
            }
            aa -= 1.0;
            let d = cd - xxi * sd;
            bb += (sdh + ap * qp * (cdh - xxi * sdh)) / d;
            cc += (cdh + ap * qp * (sdh + xxi * cdh)) / d;
            return Ok((aa, bb, cc));
        }

        if njun2 == 0 {
            let xxi = end_cap(self.bi[ix]);
            let qm = (omc + xxi * sd) / (sd * (aj - xxi * pm) + cd * (pm + xxi * aj));
            if june == -1 {
                aa *= qm;
                bb *= qm;
                cc *= qm;
                if i != is {
                    return Ok((aa, bb, cc));
                }
            }
            aa -= 1.0;
            let d = cd - xxi * sd;
            bb += (aj * qm * (cdh - xxi * sdh) - sdh) / d;
            cc += (cdh - aj * qm * (sdh + xxi * cdh)) / d;
            return Ok((aa, bb, cc));
        }

        let mut qp = sd * (pm * pp + aj * ap) + cd * (pm * ap - pp * aj);
        let qm = (ap * omc - pp * sd) / qp;
        qp = -(aj * omc + pm * sd) / qp;

        if june != 0 {
            if june < 0 {
                aa *= qm;
                bb *= qm;
                cc *= qm;
            } else {
                aa = -aa * qp;
                bb *= qp;
                cc = -cc * qp;
            }
            if i != is {
                return Ok((aa, bb, cc));
            }
        }
        aa -= 1.0;
        bb += (aj * qm + ap * qp) * sdh / sd;
        cc += (aj * qm - ap * qp) * cdh / sd;
        Ok((aa, bb, cc))
    }

    /// Collects into `self.segj` the components on segment `j` (1-based) of
    /// all basis functions that extend onto it.
    pub(crate) fn trio(&mut self, j: i32) -> Result<()> {
        self.segj.clear();
        let jx = (j - 1) as usize;
        let mut jcox = self.icon1[jx];
        let mut jend = -1;
        let mut iend = -1;

        // A free end and an end on a patch have no neighbours.
        if jcox == 0 || jcox > PCHCON {
            jcox = self.icon2[jx];
            jend = 1;
            iend = 1;
        }

        if jcox != 0 && jcox <= PCHCON {
            loop {
                if jcox < 0 {
                    jcox = -jcox;
                } else {
                    jend = -jend;
                }
                let jcoxx = (jcox - 1) as usize;

                if jcox != j {
                    let (a, b, c) = self.sbf(jcox, j)?;
                    self.segj.push(jcox, a, b, c);
                    jcox = if jend != 1 { self.icon1[jcoxx] } else { self.icon2[jcoxx] };
                    if jcox == 0 || jcox > PCHCON {
                        return Err(Error::Connection { routine: "TRIO", segment: j as usize });
                    }
                    continue;
                }

                if iend == 1 {
                    break;
                }
                jcox = self.icon2[jx];
                jend = 1;
                iend = 1;
                if jcox == 0 || jcox > PCHCON {
                    break;
                }
            }
        }

        let (a, b, c) = self.sbf(j, j)?;
        self.segj.push(j, a, b, c);
        Ok(())
    }

    /// Converts the solution vector `curx` (one amplitude per basis function
    /// and two per patch) into the constant, sine and cosine coefficients of
    /// every segment. `curx` then holds the currents at the segment centres,
    /// followed by the x, y and z components of the current on every patch,
    /// so it must have room for three values per patch.
    pub(crate) fn cabc(&mut self, curx: &mut [Complex64]) -> Result<()> {
        let n = self.n;
        for i in 0..n {
            self.air[i] = 0.0;
            self.aii[i] = 0.0;
            self.bir[i] = 0.0;
            self.bii[i] = 0.0;
            self.cir[i] = 0.0;
            self.cii[i] = 0.0;
        }

        for i in 0..n {
            let ar = curx[i].re;
            let ai = curx[i].im;
            self.tbf((i + 1) as i32, 1)?;
            self.add_basis(ar, ai);
        }

        // Current-slope-discontinuity sources add a basis function that is not
        // one of the unknowns.
        for is in 0..self.iqds.len() {
            let i = self.iqds[is] - 1;
            let saved = self.icon1[i];
            self.icon1[i] = 0;
            let built = self.tbf((i + 1) as i32, 0);
            self.icon1[i] = saved;
            built?;
            let sh = self.si[i] * 0.5;
            let last = self.segj.jco.len() - 1;
            let curd = CCJ * self.vqds[is]
                / (((2.0 * sh / self.bi[i]).ln() - 1.0)
                    * (self.segj.bx[last] * (TP * sh).cos() + self.segj.cx[last] * (TP * sh).sin())
                    * self.wlam);
            self.add_basis(curd.re, curd.im);
        }

        for i in 0..n {
            curx[i] = cmplx(self.air[i] + self.cir[i], self.aii[i] + self.cii[i]);
        }

        // Surface currents: from the two tangential components of every
        // patch to x, y and z components, from the last patch backwards
        // because the vector grows in place.
        let m = self.m;
        let mut jco1 = n + 2 * m;
        let mut jco2 = jco1 + m;
        for i in (0..m).rev() {
            jco1 -= 2;
            jco2 -= 3;
            let p = &self.patches[i];
            let cs1 = curx[jco1];
            let cs2 = curx[jco1 + 1];
            curx[jco2] = cs1 * p.t1x + cs2 * p.t2x;
            curx[jco2 + 1] = cs1 * p.t1y + cs2 * p.t2y;
            curx[jco2 + 2] = cs1 * p.t1z + cs2 * p.t2z;
        }
        Ok(())
    }

    /// Adds the basis function held in `self.segj`, scaled by `ar + j ai`, to
    /// the segment coefficients.
    fn add_basis(&mut self, ar: f64, ai: f64) {
        for jx in 0..self.segj.jco.len() {
            let j = (self.segj.jco[jx] - 1) as usize;
            self.air[j] += self.segj.ax[jx] * ar;
            self.aii[j] += self.segj.ax[jx] * ai;
            self.bir[j] += self.segj.bx[jx] * ar;
            self.bii[j] += self.segj.bx[jx] * ai;
            self.cir[j] += self.segj.cx[jx] * ar;
            self.cii[j] += self.segj.cx[jx] * ai;
        }
    }
}
