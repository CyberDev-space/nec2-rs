//! Impedance loading of segments (LD cards).
//!
//! Ported from `load` and `zint` (`calculations.c`) of nec2c 1.3.1. `zint`
//! follows the original NEC-2 Fortran: nec2c 1.3.1 lost the branches for thin
//! or poorly conducting wires and a constant has the wrong exponent there
//! (see `docs/nec2c-deviations.md`).

use num_complex::Complex64;

use super::{C0, Engine};
use crate::consts::{PI, POT, TP, cmplx};
use crate::error::{Error, Result};

/// One loading card.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Load {
    /// 0 series RLC, 1 parallel RLC, 2 series RLC per metre, 3 parallel RLC
    /// per metre, 4 fixed impedance, 5 wire conductivity.
    pub ldtyp: i32,
    /// Tag of the loaded segments; 0 selects by absolute segment number.
    pub ldtag: i32,
    /// First and last loaded segment within the tag (or absolute numbers);
    /// both zero means all segments.
    pub ldtagf: i32,
    pub ldtagt: i32,
    /// Resistance (ohms), or ohms per metre, or the real part of the
    /// impedance, or the conductivity in S/m, depending on `ldtyp`.
    pub zlr: f64,
    /// Inductance (henry), or henry per metre, or the imaginary part of the
    /// impedance.
    pub zli: f64,
    /// Capacitance (farad), or farad per metre.
    pub zlc: f64,
}

/// One line of the loading table.
#[derive(Debug, Clone, PartialEq)]
pub struct LoadLine {
    pub tag: i32,
    pub from: i32,
    pub thru: i32,
    /// Resistance, inductance, capacitance, impedance (real, imaginary) and
    /// conductivity columns; zero means the column is blank.
    pub values: [f64; 6],
    /// Circuit type as NEC-2 prints it.
    pub kind: &'static str,
}

/// What [`Engine::set_loads`] did.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct LoadingReport {
    pub lines: Vec<LoadLine>,
    /// Some segments were named by more than one card; their impedances add.
    pub loaded_twice: bool,
}

/// Internal impedance per unit length of a round wire, times the wavelength.
/// `sigl` is the conductivity times the wavelength and `rolam` the radius in
/// wavelengths.
fn zint(sigl: f64, rolam: f64) -> Complex64 {
    const TPCMU: f64 = 2.368705e+3;
    const CMOTP: f64 = 60.0;
    const CC: [Complex64; 14] = [
        Complex64::new(6.0e-7, 1.9e-6),
        Complex64::new(-3.4e-6, 5.1e-6),
        Complex64::new(-2.52e-5, 0.0),
        Complex64::new(-9.06e-5, -9.01e-5),
        Complex64::new(0.0, -9.765e-4),
        Complex64::new(0.0110486, -0.0110485),
        Complex64::new(0.0, -0.3926991),
        Complex64::new(1.6e-6, -3.2e-6),
        Complex64::new(1.17e-5, -2.4e-6),
        Complex64::new(3.46e-5, 3.38e-5),
        Complex64::new(5.0e-7, 2.452e-4),
        Complex64::new(-1.3813e-3, 1.3811e-3),
        Complex64::new(-6.25001e-2, -1.0e-7),
        Complex64::new(0.7071068, 0.7071068),
    ];
    const CN: Complex64 = Complex64::new(0.70710678, 0.70710678);

    let x = (TPCMU * sigl).sqrt() * rolam;
    let br1 = if x > 110.0 {
        // Skin depth much smaller than the radius.
        cmplx(0.70710678, -0.70710678)
    } else if x > 8.0 {
        // Asymptotic expansions of the Kelvin functions.
        let th = |d: f64| (((((CC[0] * d + CC[1]) * d + CC[2]) * d + CC[3]) * d + CC[4]) * d + CC[5]) * d + CC[6];
        let ph = |d: f64| (((((CC[7] * d + CC[8]) * d + CC[9]) * d + CC[10]) * d + CC[11]) * d + CC[12]) * d + CC[13];
        let f = (POT / x).sqrt() * (-CN * x + th(-8.0 / x)).exp();
        let g = (CN * x + th(8.0 / x)).exp() / (TP * x).sqrt();
        let br2 = cmplx(0.0, 1.0) * f / PI;
        let br1 = g + br2;
        let br2 = g * ph(8.0 / x) - br2 * ph(-8.0 / x);
        br1 / br2
    } else {
        // Polynomial approximations of ber, bei and their derivatives.
        let mut y = x / 8.0;
        y *= y;
        let s = y * y;
        let ber =
            ((((((-9.01e-6 * s + 1.22552e-3) * s - 0.08349609) * s + 2.6419140) * s - 32.363456) * s + 113.77778) * s
                - 64.0)
                * s
                + 1.0;
        let bei =
            ((((((1.1346e-4 * s - 0.01103667) * s + 0.52185615) * s - 10.567658) * s + 72.817777) * s - 113.77778) * s
                + 16.0)
                * y;
        let br1 = cmplx(ber, bei);
        let ber = (((((((-3.94e-6 * s + 4.5957e-4) * s - 0.02609253) * s + 0.66047849) * s - 6.0681481) * s
            + 14.222222)
            * s
            - 4.0)
            * y)
            * x;
        let bei =
            ((((((4.609e-5 * s - 3.79386e-3) * s + 0.14677204) * s - 2.3116751) * s + 11.377778) * s - 10.666667) * s
                + 0.5)
                * x;
        br1 / cmplx(ber, bei)
    };
    cmplx(0.0, 1.0) * (CMOTP / sigl).sqrt() * br1 / rolam
}

impl Engine {
    /// Removes all loads.
    pub fn clear_loads(&mut self) {
        self.nload = 0;
        self.zarray.clear();
    }

    /// Computes the impedance per unit length of every loaded segment at the
    /// current frequency (`load`).
    pub fn set_loads(&mut self, loads: &[Load]) -> Result<LoadingReport> {
        let tpcj = cmplx(0.0, 1.883698955e+9);
        let n = self.n;
        let wlam = self.wlam;
        self.nload = loads.len();
        self.zarray.clear();
        self.zarray.resize(n, C0);
        let mut report = LoadingReport::default();

        for card in loads {
            let Load { ldtyp, ldtag: ldtags, ldtagf, ldtagt, zlr, zli, zlc } = *card;
            if !(0..=5).contains(&ldtyp) {
                return Err(Error::ImproperLoadType(ldtyp));
            }

            // Segments named by the card: all, a range of absolute numbers,
            // or a range within a tag. A range that reaches outside the
            // structure loads the segments that exist (NEC-2 does not check).
            let (mut l1, mut l2) = (1, n as i32);
            if ldtags == 0 && (ldtagf != 0 || ldtagt != 0) {
                l1 = ldtagf.max(1);
                l2 = ldtagt.min(n as i32);
            }

            let mut ichk = 0;
            for i in (l1 - 1) as usize..l2.max(0) as usize {
                if ldtags != 0 {
                    if ldtags != self.itag[i] {
                        continue;
                    }
                    if ldtagf != 0 {
                        ichk += 1;
                        if ichk < ldtagf || ichk > ldtagt {
                            continue;
                        }
                    } else {
                        ichk = 1;
                    }
                } else {
                    ichk = 1;
                }

                // Impedance per unit length times the wavelength.
                let si = self.si[i];
                let zt = match ldtyp {
                    0 => {
                        let mut zt = zlr / si + tpcj * zli / (si * wlam);
                        if zlc.abs() > 1.0e-20 {
                            zt += wlam / (tpcj * si * zlc);
                        }
                        zt
                    }
                    1 => {
                        let mut zt = tpcj * si * zlc / wlam;
                        if zli.abs() > 1.0e-20 {
                            zt += si * wlam / (tpcj * zli);
                        }
                        if zlr.abs() > 1.0e-20 {
                            zt += si / zlr;
                        }
                        1.0 / zt
                    }
                    2 => {
                        let mut zt = zlr * wlam + tpcj * zli;
                        if zlc.abs() > 1.0e-20 {
                            zt += 1.0 / (tpcj * si * si * zlc);
                        }
                        zt
                    }
                    3 => {
                        let mut zt = tpcj * si * si * zlc;
                        if zli.abs() > 1.0e-20 {
                            zt += 1.0 / (tpcj * zli);
                        }
                        if zlr.abs() > 1.0e-20 {
                            zt += 1.0 / (zlr * wlam);
                        }
                        1.0 / zt
                    }
                    4 => cmplx(zlr, zli) / si,
                    _ => zint(zlr * wlam, self.bi[i]),
                };

                if self.zarray[i].re.abs() + self.zarray[i].im.abs() > 1.0e-20 {
                    report.loaded_twice = true;
                }
                self.zarray[i] += zt;
            }

            if ichk == 0 && ldtags != 0 {
                return Err(Error::LoadTagNotFound(ldtags));
            }

            let (values, kind) = match ldtyp {
                0 => ([zlr, zli, zlc, 0.0, 0.0, 0.0], " SERIES "),
                1 => ([zlr, zli, zlc, 0.0, 0.0, 0.0], "PARALLEL"),
                2 => ([zlr, zli, zlc, 0.0, 0.0, 0.0], "SERIES (PER METER)"),
                3 => ([zlr, zli, zlc, 0.0, 0.0, 0.0], "PARALLEL (PER METER)"),
                4 => ([0.0, 0.0, 0.0, zlr, zli, 0.0], "FIXED IMPEDANCE "),
                _ => ([0.0, 0.0, 0.0, 0.0, 0.0, zlr], "  WIRE  "),
            };
            report.lines.push(LoadLine { tag: ldtags, from: ldtagf, thru: ldtagt, values, kind });
        }

        // A symmetric structure must be loaded symmetrically: the loads of
        // the first cell are copied to the others.
        let nop = n.checked_div(self.np).unwrap_or(1);
        if nop != 1 {
            for i in 0..self.np {
                let zt = self.zarray[i];
                for l2 in 1..nop {
                    self.zarray[i + l2 * self.np] = zt;
                }
            }
        }
        Ok(report)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// DC limit: the internal impedance of a thin wire is its DC resistance,
    /// 1 / (sigma pi a^2) per metre.
    #[test]
    fn zint_tends_to_dc_resistance_for_thin_wires() {
        let sigma = 5.8e7;
        let radius = 1.0e-5;
        let wavelength = 3000.0;
        let z = zint(sigma * wavelength, radius / wavelength) / wavelength;
        let dc = 1.0 / (sigma * std::f64::consts::PI * radius * radius);
        assert!((z.re / dc - 1.0).abs() < 1e-3, "{} vs {}", z.re, dc);
    }

    /// The three branches of `zint` join smoothly.
    #[test]
    fn zint_is_continuous_across_its_branches() {
        for x0 in [8.0, 110.0] {
            let sigl = 5.8e7;
            let rolam = |x: f64| x / (2.368705e+3_f64 * sigl).sqrt();
            let below = zint(sigl, rolam(x0 * 0.9999));
            let above = zint(sigl, rolam(x0 * 1.0001));
            assert!(((below - above) / above).norm() < 1e-2, "jump at x = {x0}: {below} vs {above}");
        }
    }
}
