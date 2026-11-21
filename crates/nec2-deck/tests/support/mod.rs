//! Helpers shared by the integration tests: comparison of a report with a
//! reference report of nec2c, and access to the tables of a report.
//!
//! The two texts are compared line by line and token by token. Words must be
//! equal. A number must be printed in the same way, and it may differ from
//! the reference by a few units of the last digit that the reference prints,
//! which absorbs differences in rounding.
//! Quantities that are numerical noise in both programs (the phase of a zero
//! field, a gain below -150 dB, a current component a billion times smaller
//! than the current on its patch) are not compared.

#![allow(dead_code)]

use std::fmt::Write as _;

/// How many units of the last printed digit two numbers may differ by.
const DEFAULT_UNITS: f64 = 1.0;

/// One mismatch between the reports.
#[derive(Debug)]
pub struct Mismatch {
    pub line: usize,
    pub reference: String,
    pub actual: String,
    pub detail: String,
}

/// Lines that are never compared: the banner with the program name and the
/// lines with run times.
fn skipped(line: &str) -> bool {
    line.contains('|') || line.contains("________") || line.contains("msec")
}

/// Value of one unit in the last printed digit of a number token, e.g. 0.01
/// for "-15.24" and 1e-6 for "8.9520E-02".
fn last_digit_unit(token: &str) -> f64 {
    let (mantissa, exponent) = match token.split_once(['E', 'e']) {
        Some((m, x)) => (m, x.parse::<i32>().unwrap_or(0)),
        None => (token, 0),
    };
    let decimals = mantissa.split_once('.').map_or(0, |(_, fraction)| fraction.len() as i32);
    10f64.powi(exponent - decimals)
}

/// How a number is printed: the digits after the decimal point and the
/// digits of the exponent, for example (5, 2) for "-1.25000E-02".
fn layout(token: &str) -> (usize, usize) {
    let (mantissa, exponent) = token.split_once(['E', 'e']).unwrap_or((token, ""));
    let decimals = mantissa.split_once('.').map_or(0, |(_, fraction)| fraction.len());
    (decimals, exponent.trim_start_matches(['+', '-']).len())
}

fn number(token: &str) -> Option<f64> {
    let first = token.as_bytes().first()?;
    if first.is_ascii_digit() || matches!(first, b'-' | b'+' | b'.') { token.parse().ok() } else { None }
}

/// The table a line belongs to, found from the title above it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Section {
    Other,
    Currents,
    PatchCurrents,
    Pattern,
    GroundWave,
    NearField,
    PowerBudget,
}

fn section_of(title: &str) -> Option<Section> {
    if !title.contains("---") {
        return None;
    }
    Some(if title.contains("CURRENTS AND LOCATION") || title.contains("CHARGE DENSITIES") {
        Section::Currents
    } else if title.contains("SURFACE PATCH CURRENTS") {
        Section::PatchCurrents
    } else if title.contains("RADIATION PATTERNS") {
        Section::Pattern
    } else if title.contains("RADIATED FIELDS NEAR GROUND") {
        Section::GroundWave
    } else if title.contains("NEAR ELECTRIC FIELDS") || title.contains("NEAR MAGNETIC FIELDS") {
        Section::NearField
    } else if title.contains("POWER BUDGET") {
        Section::PowerBudget
    } else if title.contains("ANGLES")
        || title.contains("LOCATION")
        || title.contains("CURRENT (AMPS)")
        || title.contains("SURFACE COMPONENTS")
        || title.contains("PATCH CENTER")
    {
        // Column headings inside a table do not end it.
        return None;
    } else {
        Section::Other
    })
}

/// Difference of two angles in degrees, ignoring whole turns of `period`.
fn angle_difference(a: f64, b: f64, period: f64) -> f64 {
    let d = (a - b).rem_euclid(period);
    d.min(period - d)
}

/// Compares `actual` with `reference`. `units` scales the tolerance for
/// decks whose results are sensitive to rounding. With `strict` every number
/// is compared to one unit of its last digit, without the noise rules.
pub fn compare(reference: &str, actual: &str, units: Option<f64>, strict: bool) -> Vec<Mismatch> {
    let units = units.unwrap_or(DEFAULT_UNITS);
    let mut mismatches = Vec::new();
    let reference: Vec<&str> = reference.lines().collect();
    let actual: Vec<&str> = actual.lines().collect();
    let mut section = Section::Other;
    let mut input_power = 0.0;

    for index in 0..reference.len().max(actual.len()) {
        let r = reference.get(index).copied().unwrap_or("<missing line>");
        let a = actual.get(index).copied().unwrap_or("<missing line>");
        let mut fail = |detail: String| {
            mismatches.push(Mismatch { line: index + 1, reference: r.to_string(), actual: a.to_string(), detail });
        };
        if skipped(r) && skipped(a) {
            continue;
        }
        if let Some(new_section) = section_of(r) {
            section = new_section;
        }

        let rt: Vec<&str> = r.split_whitespace().collect();
        let at: Vec<&str> = a.split_whitespace().collect();
        if rt.len() != at.len() {
            fail(format!("{} tokens instead of {}", at.len(), rt.len()));
            continue;
        }

        let value = |tokens: &[&str], k: usize| tokens.get(k).and_then(|t| number(t));
        if section == Section::PowerBudget && r.contains("INPUT POWER") {
            input_power = rt.iter().rev().find_map(|t| number(t)).unwrap_or(0.0);
        }
        for k in 0..rt.len() {
            let (Some(x), Some(y)) = (number(rt[k]), number(at[k])) else {
                if rt[k] != at[k] {
                    // The sense of rotation is undefined for a negligible field.
                    if !strict && section == Section::Pattern && noise_gain(&rt) {
                        continue;
                    }
                    fail(format!("token {}: {:?} instead of {:?}", k + 1, at[k], rt[k]));
                }
                continue;
            };
            if layout(rt[k]) != layout(at[k]) {
                fail(format!("token {}: {} is not printed like {}", k + 1, at[k], rt[k]));
                continue;
            }
            let unit = last_digit_unit(rt[k]);
            let mut tolerance = units * unit * 1.0001;
            let mut difference = (x - y).abs();

            match section {
                _ if strict => {}
                Section::PowerBudget => {
                    // Losses that are zero come out as rounding noise.
                    tolerance = tolerance.max(1.0e-9 * input_power.abs());
                }
                Section::GroundWave if rt.len() == 9 => {
                    // rho phi z |Et| phase |Ep| phase |Er| phase
                    let largest = [3, 5, 7].iter().filter_map(|&c| value(&rt, c)).fold(0.0, f64::max);
                    if k >= 3 {
                        let magnitude = value(&rt, k - (k + 1) % 2).unwrap_or(0.0);
                        if magnitude < 1.0e-9 * largest {
                            continue;
                        }
                        if k % 2 == 0 {
                            difference = angle_difference(x, y, 360.0);
                        }
                    }
                }
                Section::Currents if rt.len() == 10 => {
                    // seg tag x y z length real imag magnitude phase
                    let magnitude = value(&rt, 8).unwrap_or(0.0);
                    if k == 9 {
                        if magnitude < 1.0e-12 {
                            continue;
                        }
                        difference = angle_difference(x, y, 360.0);
                    }
                    if (6..=8).contains(&k) {
                        // Components much smaller than the magnitude carry
                        // the rounding error of the magnitude.
                        tolerance = tolerance.max(units * 1.0e-5 * magnitude);
                    }
                }
                Section::PatchCurrents if rt.len() == 14 => {
                    // patch x y z |J1| phase |J2| phase xr xi yr yi zr zi
                    let largest = [4, 6].iter().filter_map(|&c| value(&rt, c)).fold(0.0, f64::max);
                    if k == 5 || k == 7 {
                        if value(&rt, k - 1).unwrap_or(0.0) < 1.0e-9 * largest {
                            continue;
                        }
                        difference = angle_difference(x, y, 360.0);
                    } else if k >= 4 {
                        tolerance = tolerance.max(units * 1.0e-9 * largest);
                    }
                }
                Section::Pattern if rt.len() >= 11 => {
                    // theta phi gain gain total axial tilt [sense] |Et| phase |Ep| phase
                    let n = rt.len();
                    let negligible = noise_gain(&rt);
                    if (2..=4).contains(&k) {
                        // Gains below the numerical noise floor.
                        if x < -150.0 && y < -150.0 {
                            continue;
                        }
                        // A component far below the total is a small
                        // difference of large numbers.
                        let total = value(&rt, 4).unwrap_or(0.0);
                        if total - x > 40.0 {
                            tolerance = tolerance.max(units * 0.01 * 10f64.powf((total - x - 40.0) / 20.0));
                        }
                    } else if k == 5 || k == 6 {
                        if negligible {
                            continue;
                        }
                        if k == 6 {
                            difference = angle_difference(x, y, 180.0);
                            // The tilt of a nearly circular ellipse is undefined.
                            let axial = value(&rt, 5).unwrap_or(0.0);
                            if axial > 0.99 {
                                continue;
                            }
                        }
                    } else if k == n - 3 || k == n - 1 {
                        // Phases: compared when the magnitude is not noise.
                        let magnitude = value(&rt, k - 1).unwrap_or(0.0);
                        let other = value(&rt, if k == n - 1 { n - 4 } else { n - 2 }).unwrap_or(0.0);
                        if magnitude < 1.0e-9 || magnitude < 1.0e-7 * other {
                            continue;
                        }
                        difference = angle_difference(x, y, 360.0);
                    } else if k == n - 4 || k == n - 2 {
                        let other = value(&rt, if k == n - 2 { n - 4 } else { n - 2 }).unwrap_or(0.0);
                        tolerance = tolerance.max(units * 1.0e-9 * other).max(1.0e-11);
                    }
                }
                _ => {}
            }

            if difference > tolerance {
                fail(format!(
                    "token {}: {} instead of {} (off by {:.1} units of the last digit)",
                    k + 1,
                    at[k],
                    rt[k],
                    difference / unit
                ));
            }
        }
    }
    mismatches
}

/// A pattern row whose total gain is numerical noise.
fn noise_gain(tokens: &[&str]) -> bool {
    tokens.get(4).and_then(|t| number(t)).is_some_and(|gain| gain < -150.0)
}

/// The report without the warnings that nec2-rs adds. They are whole lines
/// that begin with `WARNING:`; nec2c prints no such line, and without them
/// the report has the layout of nec2c.
pub fn without_warnings(report: &str) -> String {
    report.lines().filter(|line| !line.trim_start().starts_with("WARNING:")).collect::<Vec<_>>().join("\n")
}

/// The warnings that nec2-rs added to a report.
pub fn warnings(report: &str) -> Vec<&str> {
    report.lines().map(str::trim_start).filter(|line| line.starts_with("WARNING:")).collect()
}

/// Formats the first mismatches for a test failure message.
pub fn describe(name: &str, mismatches: &[Mismatch]) -> String {
    let mut text = format!("{name}: {} lines differ from the nec2c reference", mismatches.len());
    for m in mismatches.iter().take(8) {
        let _ = write!(
            text,
            "\n  line {}: {}\n    ref: {}\n    got: {}",
            m.line,
            m.detail,
            m.reference.trim_end(),
            m.actual.trim_end()
        );
    }
    text
}

/// A token of a report as a number.
pub fn num(token: &str) -> f64 {
    number(token).unwrap_or_else(|| panic!("{token:?} is not a number"))
}

/// The rows of every table whose title line contains `title`: one list of
/// rows per occurrence, each row split into tokens. The rows of a table are
/// the consecutive lines that begin with a number, after at most ten lines
/// of headings.
pub fn tables<'a>(report: &'a str, title: &str) -> Vec<Vec<Vec<&'a str>>> {
    let lines: Vec<&str> = report.lines().collect();
    let is_row = |line: &str| line.split_whitespace().next().is_some_and(|token| number(token).is_some());
    let mut found = Vec::new();
    for (index, _) in lines.iter().enumerate().filter(|(_, line)| line.contains(title)) {
        let headings = lines[index + 1..].iter().take(10).take_while(|line| !is_row(line)).count();
        let rows = lines[index + 1 + headings..]
            .iter()
            .take_while(|line| is_row(line))
            .map(|line| line.split_whitespace().collect())
            .collect();
        found.push(rows);
    }
    found
}

/// The value printed after `label` and an equals sign or colon, for every
/// occurrence, for example `value(report, "NETWORK LOSS")`.
pub fn values(report: &str, label: &str) -> Vec<f64> {
    report
        .lines()
        .filter_map(|line| line.split_once(label))
        .filter_map(|(_, rest)| rest.split_whitespace().find_map(number))
        .collect()
}
