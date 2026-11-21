//! Numerical constants of NEC-2.
//!
//! The values are the rounded ones used by NEC-2 and nec2c (for example
//! `PI = 3.141592654`, `CVEL = 299.8`). They are kept on purpose: results are
//! compared against other NEC-2 programs digit for digit, and a more precise
//! value of the speed of light or of the free-space impedance would shift
//! every result in the fifth or sixth significant digit.

use num_complex::Complex64;

pub const PI: f64 = 3.141592654;
/// pi / 2
pub const POT: f64 = 1.570796327;
/// 2 pi
pub const TP: f64 = 6.283185308;
/// pi / 5
pub const PTP: f64 = 0.6283185308;
/// j 2 pi
pub const TPJ: Complex64 = Complex64::new(0.0, 6.283185308);
/// 8 pi
pub const PI8: f64 = 25.13274123;
/// 10 pi
pub const PI10: f64 = 31.41592654;
/// Degrees to radians.
pub const TA: f64 = 1.745329252E-02;
/// Radians to degrees.
pub const TD: f64 = 57.29577951;
/// Impedance of free space, ohms.
pub const ETA: f64 = 376.73;
/// Speed of light in units of 1e6 m/s, so that wavelength = CVEL / MHz.
pub const CVEL: f64 = 299.8;
/// 1 / ETA
pub const RETA: f64 = 2.654420938E-3;
/// 2 / sqrt(pi)
pub const TOSP: f64 = 1.128379167;
pub const ACCS: f64 = 1.0E-12;
/// sqrt(pi)
pub const SP: f64 = 1.772453851;
/// 4 pi
pub const FPI: f64 = 12.56637062;
/// -j / 60
pub const CCJ: Complex64 = Complex64::new(0.0, -0.01666666667);
/// j ETA / (8 pi^2)
pub const CONST1: Complex64 = Complex64::new(0.0, 4.771341189);
pub const CONST2: f64 = 4.771341188;
/// -j ETA / (4 pi)
pub const CONST3: Complex64 = Complex64::new(0.0, -29.97922085);
/// j ETA / 2
pub const CONST4: Complex64 = Complex64::new(0.0, 188.365);
pub const GAMMA: f64 = 0.5772156649;
pub const C1: f64 = -0.02457850915;
pub const C2: f64 = 0.3674669052;
pub const C3: f64 = 0.7978845608;
pub const P10: f64 = 0.0703125;
pub const P20: f64 = 0.1121520996;
pub const Q10: f64 = 0.125;
pub const Q20: f64 = 0.0732421875;
pub const P11: f64 = 0.1171875;
pub const P21: f64 = 0.1441955566;
pub const Q11: f64 = 0.375;
pub const Q21: f64 = 0.1025390625;
/// pi / 4
pub const POF: f64 = 0.7853981635;
pub const MAXH: usize = 20;
pub const CRIT: f64 = 1.0E-4;
pub const NM: i32 = 131072;
pub const NTS: i32 = 4;
/// Segment ends closer than `SMIN` times the segment length are connected.
pub const SMIN: f64 = 1.0e-3;

/// `(re, im)` shorthand used throughout the engine.
#[inline]
pub(crate) fn cmplx(re: f64, im: f64) -> Complex64 {
    Complex64::new(re, im)
}

/// Phase angle of `z` in degrees.
#[inline]
pub fn cang(z: Complex64) -> f64 {
    z.im.atan2(z.re) * TD
}

/// Connection numbers above this value name a patch: `PCHCON + k` is patch
/// `k` (from 1), the first of the four that a wire end divides a patch into.
pub const PCHCON: i32 = 100_000;

/// Decibels of a power ratio; values below 1e-20 give -999.99.
#[inline]
pub fn db10(x: f64) -> f64 {
    if x < 1.0e-20 { -999.99 } else { 10.0 * x.log10() }
}

/// Decibels of a field ratio; values below 1e-20 give -999.99.
#[inline]
pub fn db20(x: f64) -> f64 {
    if x < 1.0e-20 { -999.99 } else { 20.0 * x.log10() }
}
