//! Interaction matrix: fill, factorization and solution.
//!
//! The matrix is stored the way NEC-2 stores it: element `cm[j + i * nrow]`
//! is the field at observation point `i` due to unknown `j`, so the array
//! holds the transpose of the usual impedance matrix. For a structure with
//! symmetry only the observation points of the first symmetric cell are
//! filled and the matrix is split into one block per symmetry mode.
//!
//! The unknowns of one symmetric cell are the basis functions of its
//! segments followed by two surface-current components for each of its
//! patches. The equations are the electric field along the segments and the
//! magnetic field in the two tangential directions at the patch centres.
//!
//! Ported from `matrix.c` of nec2c 1.3.1.

use num_complex::Complex64;

use super::{C0, C1, Engine, ILL_CONDITIONED};
use crate::consts::{PCHCON, PI, TP, cmplx};
use crate::error::{Error, Result};

/// Factors the `n` by `n` matrix held transposed in `a` (leading dimension
/// `ndim`) into a unit lower and an upper triangular matrix by the
/// Gauss-Doolittle method with row pivoting. `ip` receives the pivot rows
/// (1-based).
///
/// Returns the ratio of the smallest to the largest magnitude on the
/// diagonal of the upper factor, an indicator of the conditioning. A pivot
/// that is zero or not finite stops the factorization with the number of
/// its equation (1-based) as the error.
pub(crate) fn factr(
    n: usize,
    a: &mut [Complex64],
    ip: &mut [i32],
    ndim: usize,
    warnings: &mut Vec<String>,
) -> std::result::Result<f64, usize> {
    // Un-transpose the matrix for Gauss elimination.
    for i in 1..n {
        for j in 0..i {
            a.swap(i + j * ndim, j + i * ndim);
        }
    }

    // Smallest and largest squared magnitude of a pivot.
    let (mut smallest, mut largest) = (f64::INFINITY, 0.0f64);
    let mut scm = vec![C0; n];
    for r in 0..n {
        for k in 0..n {
            scm[k] = a[k + r * ndim];
        }

        // Steps 2 and 3.
        for j in 0..r {
            let pj = (ip[j] - 1) as usize;
            let arj = scm[pj];
            a[j + r * ndim] = arj;
            scm[pj] = scm[j];
            for i in (j + 1)..n {
                scm[i] -= a[i + j * ndim] * arj;
            }
        }

        // Step 4: find the pivot.
        let mut dmax = scm[r].norm_sqr();
        let rp1 = r + 1;
        ip[r] = rp1 as i32;
        for i in rp1..n {
            let elmag = scm[i].norm_sqr();
            if elmag >= dmax {
                dmax = elmag;
                ip[r] = (i + 1) as i32;
            }
        }
        if dmax == 0.0 || !dmax.is_finite() {
            return Err(rp1);
        }
        if dmax < 1.0e-10 {
            warnings.push(format!("PIVOT({r})= {dmax:16.8E}"));
        }
        smallest = smallest.min(dmax);
        largest = largest.max(dmax);

        let pr = (ip[r] - 1) as usize;
        a[r + r * ndim] = scm[pr];
        scm[pr] = scm[r];

        // Step 5.
        if rp1 < n {
            let arj = 1.0 / a[r + r * ndim];
            for i in rp1..n {
                a[i + r * ndim] = scm[i] * arj;
            }
        }
    }
    Ok(if n == 0 { 1.0 } else { (smallest / largest).sqrt() })
}

/// Solves `L U x = b` for the factors produced by [`factr`]; the solution
/// replaces `b`.
pub(crate) fn solve(n: usize, a: &[Complex64], ip: &[i32], b: &mut [Complex64], ndim: usize) {
    let mut scm = vec![C0; n];

    // Forward substitution.
    for i in 0..n {
        let pia = (ip[i] - 1) as usize;
        scm[i] = b[pia];
        b[pia] = b[i];
        for j in (i + 1)..n {
            b[j] -= a[j + i * ndim] * scm[i];
        }
    }

    // Backward substitution.
    for k in 0..n {
        let i = n - k - 1;
        let mut sum = C0;
        for j in (i + 1)..n {
            sum += a[i + j * ndim] * b[j];
        }
        b[i] = (scm[i] - sum) / a[i + i * ndim];
    }
}

impl Engine {
    /// Sets up the transform between the symmetric sections of the structure
    /// and its symmetry modes (`fblock`).
    pub(crate) fn fblock(&mut self) -> Result<()> {
        let nrow = self.npeq();
        let ncol = self.neq();
        self.nop = ncol / nrow;
        if nrow == ncol {
            self.ssx.clear();
            return Ok(());
        }
        let nop = self.nop;
        if nop * nrow != ncol {
            return Err(Error::Symmetry { cell: nrow, total: ncol });
        }
        self.ssx = vec![C0; nop * nop];

        // Rotational symmetry.
        if self.ipsym <= 0 {
            let phaz = TP / nop as f64;
            for i in 1..nop {
                for j in i..nop {
                    let arg = phaz * i as f64 * j as f64;
                    self.ssx[i + j * nop] = cmplx(arg.cos(), arg.sin());
                    self.ssx[j + i * nop] = self.ssx[i + j * nop];
                }
            }
            return Ok(());
        }

        // Plane symmetry: the section count is a power of two.
        if !nop.is_power_of_two() {
            return Err(Error::Symmetry { cell: nrow, total: ncol });
        }
        let mut kk = 1;
        self.ssx[0] = C1;
        for _ in 0..nop.trailing_zeros() {
            for i in 0..kk {
                for j in 0..kk {
                    let deter = self.ssx[i + j * nop];
                    self.ssx[i + (j + kk) * nop] = deter;
                    self.ssx[i + kk + (j + kk) * nop] = -deter;
                    self.ssx[i + kk + j * nop] = deter;
                }
            }
            kk *= 2;
        }
        Ok(())
    }

    /// Fills the interaction matrix (`cmset`). `rkhx` is the separation in
    /// wavelengths beyond which segments interact as point sources and
    /// `iexkx` selects the extended thin-wire kernel.
    pub(crate) fn cmset(&mut self, rkhx: f64, iexkx: i32) -> Result<()> {
        let neq = self.neq();
        let npeq = self.npeq();
        let nrow = neq;
        let (n, np, mp) = (self.n, self.np, self.mp);
        let mp2 = 2 * mp;
        self.dj.rkh = rkhx;
        self.dj.iexk = iexkx;
        // A structure too large for the memory is an error; without this the
        // allocation would end the whole process.
        let too_large = Error::TooLarge { equations: neq };
        let size = nrow.checked_mul(npeq).ok_or(too_large.clone())?;
        self.cm.clear();
        self.cm.try_reserve_exact(size).map_err(|_| too_large)?;
        self.cm.resize(size, C0);

        // Wire sources.
        for j in 1..=n {
            self.trio(j as i32)?;
            // The place of a basis function among the unknowns: the patch
            // unknowns of the cells before it come first.
            for jco in self.segj.jco.iter_mut() {
                *jco += ((*jco - 1) / np as i32) * mp2 as i32;
            }
            self.cmww(j, 1, np)?;
            if mp2 != 0 {
                self.cmws(j, mp2, np * nrow)?;
            }

            // Matrix elements modified by loading.
            if self.nload == 0 || j > np {
                continue;
            }
            let zaj = self.zarray[j - 1];
            for i in 0..self.segj.jco.len() {
                let jss = self.segj.jco[i] as usize;
                self.cm[(jss - 1) + (j - 1) * nrow] -= (self.segj.ax[i] + self.segj.cx[i]) * zaj;
            }
        }

        // Patch sources, one symmetric cell after the other.
        if self.m != 0 {
            for cell in 0..self.nop {
                let (jm1, jm2) = (cell * mp + 1, (cell + 1) * mp);
                let jst = cell * npeq + np;
                if np != 0 {
                    self.cmsw(jm1, jm2, np, jst)?;
                }
                self.cmss(jm1, jm2, mp, jst + np * nrow);
            }
        }

        if self.nop == 1 {
            return Ok(());
        }

        // Combine the elements of the symmetric sections into symmetry modes.
        let nop = self.nop;
        let mut scm = vec![C0; nop];
        for i in 0..npeq {
            for j in 0..npeq {
                for k in 0..nop {
                    scm[k] = self.cm[j + k * npeq + i * nrow];
                }
                let mut deter = scm[0];
                for value in &scm[1..] {
                    deter += value;
                }
                self.cm[j + i * nrow] = deter;

                for k in 1..nop {
                    let mut deter = scm[0];
                    for kk in 1..nop {
                        deter += scm[kk] * self.ssx[k + kk * nop];
                    }
                    self.cm[j + k * npeq + i * nrow] = deter;
                }
            }
        }
        Ok(())
    }

    /// Adds to the matrix the interaction of source segment `j` (1-based)
    /// with the observation segments `i1..=i2` (1-based), using the basis
    /// function components in `self.segj`.
    fn cmww(&mut self, j: usize, i1: usize, i2: usize) -> Result<()> {
        let jx = j - 1;
        let nr = self.neq();
        self.set_source_segment(jx);
        if self.dj.iexk != 0 {
            self.set_extended_kernel_ends(jx, 0);
        }

        for (ipr, i) in ((i1 - 1)..i2).enumerate() {
            let ij = i as i32 - jx as i32;
            self.efld(self.x[i], self.y[i], self.z[i], self.bi[i], ij);
            self.check_integration(j, Some(i + 1))?;
            let (cabi, sabi, salpi) = (self.cab[i], self.sab[i], self.salp[i]);
            let etk = self.dj.exk * cabi + self.dj.eyk * sabi + self.dj.ezk * salpi;
            let ets = self.dj.exs * cabi + self.dj.eys * sabi + self.dj.ezs * salpi;
            let etc = self.dj.exc * cabi + self.dj.eyc * sabi + self.dj.ezc * salpi;

            // The matrix is filled transposed.
            for k in 0..self.segj.jco.len() {
                let jc = (self.segj.jco[k] - 1) as usize;
                self.cm[jc + ipr * nr] += etk * self.segj.ax[k] + ets * self.segj.bx[k] + etc * self.segj.cx[k];
            }
        }
        Ok(())
    }

    /// Factors the matrix of every symmetry mode (`factrs`) and records the
    /// conditioning indicator.
    pub(crate) fn factrs(&mut self) -> Result<()> {
        let np = self.npeq();
        let nrow = self.neq();
        let mut ratio = 1.0f64;
        for kk in 0..self.nop {
            let ka = kk * np;
            match factr(np, &mut self.cm[ka..], &mut self.ip[ka..], nrow, &mut self.diag.warnings) {
                Ok(block) => ratio = ratio.min(block),
                Err(row) => {
                    return Err(Error::SingularMatrix {
                        matrix: "INTERACTION",
                        equation: ka + row,
                        frequency: self.fmhz,
                    });
                }
            }
        }
        self.pivot_ratio = ratio;
        if ratio < ILL_CONDITIONED {
            self.diag.warnings.push(format!(
                "WARNING: THE INTERACTION MATRIX IS ILL-CONDITIONED, ABOUT {:.0} DIGITS OF THE RESULT ARE LOST",
                -ratio.log10()
            ));
        }
        Ok(())
    }

    /// Solves the factored matrix equation for the right-hand side `b`
    /// (`solves`): the excitation is transformed to the symmetry modes, each
    /// mode is solved and the solutions are transformed back.
    pub(crate) fn solves(&self, b: &mut [Complex64]) {
        let npeq = self.npeq();
        let neq = self.neq();
        let nop = self.nop;
        let fnorm = 1.0 / nop as f64;
        let mut scm = vec![C0; nop];
        // The right-hand side lists all segments and then all patches. The
        // matrix has them cell by cell.
        let mixed = nop != 1 && self.n != 0 && self.m != 0;

        if nop != 1 {
            if mixed {
                let old = b[..neq].to_vec();
                for (new, &from) in b[..neq].iter_mut().zip(&self.cell_order()) {
                    *new = old[from];
                }
            }
            for i in 0..npeq {
                for k in 0..nop {
                    scm[k] = b[i + k * npeq];
                }
                let mut sum = scm[0];
                for value in &scm[1..] {
                    sum += value;
                }
                b[i] = sum * fnorm;

                for k in 1..nop {
                    let mut sum = scm[0];
                    for j in 1..nop {
                        sum += scm[j] * self.ssx[k + j * nop].conj();
                    }
                    b[i + k * npeq] = sum * fnorm;
                }
            }
        }

        for kk in 0..nop {
            let ia = kk * npeq;
            solve(npeq, &self.cm[ia..], &self.ip[ia..], &mut b[ia..], neq);
        }

        if nop == 1 {
            return;
        }

        for i in 0..npeq {
            for k in 0..nop {
                scm[k] = b[i + k * npeq];
            }
            let mut sum = scm[0];
            for value in &scm[1..] {
                sum += value;
            }
            b[i] = sum;

            for k in 1..nop {
                let mut sum = scm[0];
                for j in 1..nop {
                    sum += scm[j] * self.ssx[k + j * nop];
                }
                b[i + k * npeq] = sum;
            }
        }

        if mixed {
            let old = b[..neq].to_vec();
            for (&value, &to) in old.iter().zip(&self.cell_order()) {
                b[to] = value;
            }
        }
    }

    /// For every unknown in the order of the matrix (cell by cell, segments
    /// before patches), its place in the order of the right-hand side (all
    /// segments, then all patches).
    fn cell_order(&self) -> Vec<usize> {
        let (n, np, mp2) = (self.n, self.np, 2 * self.mp);
        let mut order = Vec::with_capacity(self.neq());
        for cell in 0..self.nop {
            order.extend(cell * np..(cell + 1) * np);
            order.extend(n + cell * mp2..n + (cell + 1) * mp2);
        }
        order
    }

    /// Adds to the matrix the magnetic field at the patches of the first
    /// cell due to source segment `j` (1-based), in the `nobs` patch
    /// equations that start at `offset` in the matrix (`cmws`).
    fn cmws(&mut self, j: usize, nobs: usize, offset: usize) -> Result<()> {
        let nr = self.neq();
        self.set_source_segment(j - 1);

        // Two equations for every patch: the field across the second tangent
        // vector, then across the first.
        for ipr in 0..nobs {
            let js = ipr / 2;
            let p = self.patches[js];
            let (tx, ty, tz) = if ipr % 2 == 0 {
                self.hsfld(p.x, p.y, p.z, 0.0);
                self.check_integration(j, None)?;
                (p.t2x, p.t2y, p.t2z)
            } else {
                (p.t1x, p.t1y, p.t1z)
            };
            let dj = &self.dj;
            let etk = -(dj.exk * tx + dj.eyk * ty + dj.ezk * tz) * p.salp;
            let ets = -(dj.exs * tx + dj.eys * ty + dj.ezs * tz) * p.salp;
            let etc = -(dj.exc * tx + dj.eyc * ty + dj.ezc * tz) * p.salp;
            for k in 0..self.segj.jco.len() {
                let jx = (self.segj.jco[k] - 1) as usize;
                self.cm[offset + jx + ipr * nr] +=
                    etk * self.segj.ax[k] + ets * self.segj.bx[k] + etc * self.segj.cx[k];
            }
        }
        Ok(())
    }

    /// Fills the matrix elements for the electric field along the segments
    /// `1..=i2` due to the currents on the patches `j1..=j2` (1-based),
    /// whose unknowns start at column `column` (`cmsw`).
    fn cmsw(&mut self, j1: usize, j2: usize, i2: usize, column: usize) -> Result<()> {
        let nrow = self.neq();
        // Position among the four patches around a wire end, and the end of
        // the wire that is on them.
        let mut icgo = 0;
        let mut fsign = 1.0;

        for i in 0..i2 {
            let (xi, yi, zi) = (self.x[i], self.y[i], self.z[i]);
            let (cabi, sabi, salpi) = (self.cab[i], self.sab[i], self.salp[i]);
            let mut ipch = 0;
            if self.icon1[i] >= PCHCON {
                ipch = (self.icon1[i] - PCHCON) as usize;
                fsign = -1.0;
            }
            if self.icon2[i] >= PCHCON {
                ipch = (self.icon2[i] - PCHCON) as usize;
                fsign = 1.0;
            }

            let mut emel = [C0; 9];
            for (jl, j) in (j1..=j2).enumerate() {
                let (c1, c2) = (column + 2 * jl + i * nrow, column + 2 * jl + 1 + i * nrow);
                self.set_source_patch(j - 1);

                // The patch, then its image in the ground.
                for ip in 1..=self.gnd.ksymp as i32 {
                    self.dj.ipgnd = ip;
                    if (ipch == j || icgo != 0) && ip != 2 {
                        if icgo == 0 {
                            // The four patches that the wire end divides: the
                            // current flows out of the wire over them.
                            emel = self.pcint(xi, yi, zi, cabi, sabi, salpi);
                            let pyl = PI * self.si[i] * fsign;
                            let pxl = pyl.sin();
                            let pyl = pyl.cos();
                            let exc = emel[8] * fsign;
                            self.trio((i + 1) as i32)?;
                            // The last component is the basis function of
                            // the segment itself. (nec2c takes the index
                            // from an earlier call of trio.)
                            let own = self.segj.jco.len() - 1;
                            self.cm[i + i * nrow] +=
                                exc * (self.segj.ax[own] + self.segj.bx[own] * pxl + self.segj.cx[own] * pyl);
                        }
                        self.cm[c1] = emel[icgo];
                        self.cm[c2] = emel[icgo + 4];
                        icgo += 1;
                        if icgo == 4 {
                            icgo = 0;
                        }
                        continue;
                    }

                    self.unere(xi, yi, zi);
                    let dj = &self.dj;
                    self.cm[c1] += dj.exk * cabi + dj.eyk * sabi + dj.ezk * salpi;
                    self.cm[c2] += dj.exs * cabi + dj.eys * sabi + dj.ezs * salpi;
                }
            }
        }
        Ok(())
    }

    /// Fills the matrix elements for the magnetic field at the `nobs`
    /// patches of the first cell due to the currents on the patches
    /// `j1..=j2` (1-based). `offset` is the place of the first of these
    /// elements in the matrix (`cmss`).
    fn cmss(&mut self, j1: usize, j2: usize, nobs: usize, offset: usize) {
        let nrow = self.neq();
        for il in 0..nobs {
            let pi = self.patches[il];
            let (t1xi, t1yi, t1zi) = (pi.t1x * pi.salp, pi.t1y * pi.salp, pi.t1z * pi.salp);
            let (t2xi, t2yi, t2zi) = (pi.t2x * pi.salp, pi.t2y * pi.salp, pi.t2z * pi.salp);
            let (ii1, ii2) = (2 * il, 2 * il + 1);

            for (jl, j) in (j1..=j2).enumerate() {
                let (jj1, jj2) = (2 * jl, 2 * jl + 1);
                self.set_source_patch(j - 1);
                self.hintg(pi.x, pi.y, pi.z);
                let dj = &self.dj;
                let mut g11 = -(t2xi * dj.exk + t2yi * dj.eyk + t2zi * dj.ezk);
                let g12 = -(t2xi * dj.exs + t2yi * dj.eys + t2zi * dj.ezs);
                let g21 = -(t1xi * dj.exk + t1yi * dj.eyk + t1zi * dj.ezk);
                let mut g22 = -(t1xi * dj.exs + t1yi * dj.eys + t1zi * dj.ezs);
                if il + 1 == j {
                    g11 -= 0.5;
                    g22 += 0.5;
                }
                self.cm[offset + jj1 + ii1 * nrow] = g11;
                self.cm[offset + jj2 + ii1 * nrow] = g12;
                self.cm[offset + jj1 + ii2 * nrow] = g21;
                self.cm[offset + jj2 + ii2 * nrow] = g22;
            }
        }
    }

    /// Fills and factors the interaction matrix for the current frequency,
    /// loads and ground.
    ///
    /// `rkh` is the separation in wavelengths beyond which the lumped-source
    /// approximation is used (1.0 by default in NEC-2) and `extended_kernel`
    /// selects the extended thin-wire kernel.
    pub fn fill_and_factor(&mut self, rkh: f64, extended_kernel: bool) -> Result<()> {
        self.fblock()?;
        self.cmset(rkh, i32::from(extended_kernel))?;
        self.factrs()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Row `i` of the matrix is stored at `a[3 * i..]`.
    fn factor(rows: [[f64; 3]; 3]) -> std::result::Result<f64, usize> {
        let mut a: Vec<Complex64> = rows.iter().flatten().map(|&v| cmplx(v, 0.0)).collect();
        factr(3, &mut a, &mut [0; 3], 3, &mut Vec::new())
    }

    #[test]
    fn two_identical_rows_are_a_singular_matrix() {
        assert_eq!(factor([[2.0, 1.0, 1.0], [4.0, 3.0, 3.0], [4.0, 3.0, 3.0]]), Err(3));
        assert_eq!(factor([[0.0; 3]; 3]), Err(1));
        assert_eq!(factor([[1.0, 0.0, 0.0], [0.0, f64::NAN, 0.0], [0.0, 0.0, 1.0]]), Err(2));
    }

    #[test]
    fn pivot_ratio_is_smallest_over_largest_pivot() {
        assert_eq!(factor([[4.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 0.5]]), Ok(0.125));
    }
}
