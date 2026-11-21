//! Number formatting in the style of C `printf`, as NEC reports need it.

/// `%{width}.{prec}E`: scientific notation with a signed exponent of at
/// least two digits, right-aligned in `width` columns.
pub(crate) fn e(x: f64, width: usize, prec: usize) -> String {
    format!("{:>width$}", exp_notation(x, prec, false))
}

/// `%+{width}.{prec}E`: like [`e`] with an explicit plus sign.
pub(crate) fn e_plus(x: f64, width: usize, prec: usize) -> String {
    format!("{:>width$}", exp_notation(x, prec, true))
}

fn exp_notation(x: f64, prec: usize, plus: bool) -> String {
    let sign = if x.is_sign_negative() {
        "-"
    } else if plus {
        "+"
    } else {
        ""
    };
    if x.is_nan() {
        return format!("{sign}NAN");
    }
    if x.is_infinite() {
        return format!("{sign}INF");
    }
    let formatted = format!("{:.prec$E}", x.abs());
    let (mantissa, exponent) = formatted.split_once('E').expect("the E format always has an exponent");
    let exponent: i32 = exponent.parse().expect("the exponent is an integer");
    let exponent_sign = if exponent < 0 { '-' } else { '+' };
    format!("{sign}{mantissa}E{exponent_sign}{:02}", exponent.abs())
}

/// `%{width}.{prec}f`.
pub(crate) fn f(x: f64, width: usize, prec: usize) -> String {
    format!("{x:>width$.prec$}")
}

/// `%+{width}.{prec}f`.
pub(crate) fn f_plus(x: f64, width: usize, prec: usize) -> String {
    format!("{x:>+width$.prec$}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scientific_notation_matches_printf() {
        assert_eq!(e(84.823, 11, 4), " 8.4823E+01");
        assert_eq!(e(-0.0050549, 11, 4), "-5.0549E-03");
        assert_eq!(e(0.0, 11, 4), " 0.0000E+00");
        assert_eq!(e(1.0e100, 12, 5), "1.00000E+100");
        assert_eq!(e_plus(0.3, 11, 4), "+3.0000E-01");
        assert_eq!(e(299.8, 12, 5), " 2.99800E+02");
    }

    #[test]
    fn fixed_notation_matches_printf() {
        assert_eq!(f(-15.236, 8, 2), "  -15.24");
        assert_eq!(f(-0.0, 9, 4), "  -0.0000");
        assert_eq!(f_plus(4.0, 7, 4), "+4.0000");
    }
}
