//! Reading the lines and fields of a NEC deck.
//!
//! Fields are free-format, separated by blanks, tabs or commas, and missing
//! fields are zero, as in nec2c 1.3.1 (`readmn`, `readgm`). Unlike nec2c,
//! lines may be indented and of any length, a line that ends at the end of
//! the file is complete, an apostrophe starts a comment like `#` does, text
//! after `!` or an apostrophe on a card is ignored, an integer field may be
//! written as a real number with an integral value, and `D` is accepted as
//! exponent letter.

/// Mnemonics of the geometry section.
pub(crate) const GEOMETRY_CARDS: [&str; 13] =
    ["GW", "GC", "GX", "GR", "GS", "GE", "GM", "SP", "SM", "GA", "SC", "GH", "GF"];
/// Mnemonics of the program-control section.
pub(crate) const CONTROL_CARDS: [&str; 20] = [
    "FR", "LD", "GN", "EX", "NT", "TL", "XQ", "GD", "RP", "NX", "PT", "KH", "NE", "NH", "PQ", "EK", "CP", "PL", "EN",
    "WG",
];

/// One card: its line number in the deck (from 1) and its text without
/// leading blanks, the mnemonic in upper case.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Line {
    pub number: usize,
    pub text: String,
}

/// Line reader over the text of a deck.
pub(crate) struct Lines<'a> {
    rest: std::iter::Enumerate<std::str::Lines<'a>>,
    /// The card returned last, for error messages; `None` at the end of the deck.
    pub last: Option<Line>,
}

impl<'a> Lines<'a> {
    pub(crate) fn new(text: &'a str) -> Self {
        Lines { rest: text.lines().enumerate(), last: None }
    }

    /// Returns the next card, or `None` at the end of the deck. Empty lines
    /// and comment lines are skipped.
    pub(crate) fn next_line(&mut self) -> Option<Line> {
        loop {
            let Some((index, line)) = self.rest.next() else {
                self.last = None;
                return None;
            };
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') || line.starts_with('\'') {
                continue;
            }
            let mut text = line.to_string();
            let head = text.char_indices().nth(2).map_or(text.len(), |(i, _)| i);
            text[..head].make_ascii_uppercase();
            self.last = Some(Line { number: index + 1, text });
            return self.last.clone();
        }
    }
}

/// The two-letter mnemonic of a card.
pub(crate) fn mnemonic(line: &str) -> &str {
    line.get(..2).unwrap_or(line)
}

/// A card split into its integer and floating-point fields.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Card {
    pub mnemonic: String,
    pub ints: Vec<i32>,
    pub floats: Vec<f64>,
}

/// C `atof` for a string that starts at a sign, digit or decimal point: the
/// longest prefix that reads as a decimal number, or zero.
fn atof(s: &[u8]) -> f64 {
    let mut end = 0;
    let mut i = 0;
    if matches!(s.first(), Some(b'+' | b'-')) {
        i = 1;
    }
    let digits_from = i;
    while s.get(i).is_some_and(u8::is_ascii_digit) {
        i += 1;
    }
    let mut mantissa_digits = i - digits_from;
    if s.get(i) == Some(&b'.') {
        i += 1;
        let fraction_from = i;
        while s.get(i).is_some_and(u8::is_ascii_digit) {
            i += 1;
        }
        mantissa_digits += i - fraction_from;
    }
    if mantissa_digits > 0 {
        end = i;
        if matches!(s.get(i), Some(b'e' | b'E' | b'd' | b'D')) {
            let mut k = i + 1;
            if matches!(s.get(k), Some(b'+' | b'-')) {
                k += 1;
            }
            let exponent_from = k;
            while s.get(k).is_some_and(u8::is_ascii_digit) {
                k += 1;
            }
            if k > exponent_from {
                end = k;
            }
        }
    }
    let text = String::from_utf8_lossy(&s[..end]).replace(['d', 'D'], "E");
    text.parse().unwrap_or(0.0)
}

/// Splits `line` into `nint` integers followed by `nflt` floating-point
/// numbers. `kind` names the card family in error messages ("COMMAND" or
/// "GEOMETRY").
pub(crate) fn parse_card(line: &str, nint: usize, nflt: usize, kind: &str) -> Result<Card, String> {
    // Text after an exclamation mark or an apostrophe is a comment.
    let line = line.split(['!', '\'']).next().unwrap_or("");
    let bytes = line.as_bytes();
    if bytes.len() < 2 {
        return Err(format!("\n  {kind} DATA CARD ERROR:\n  CARD'S MNEMONIC CODE TOO SHORT OR MISSING."));
    }
    let mut card = Card { mnemonic: mnemonic(line).to_string(), ints: vec![0; nint], floats: vec![0.0; nflt] };
    let at = |i: usize| bytes.get(i).copied().unwrap_or(0);
    let starts_number = |c: u8| c.is_ascii_digit() || c == b'+' || c == b'-' || c == b'.';
    let is_separator = |c: u8| c == b' ' || c == b'\t' || c == b',' || c == 0;

    let mut idx = 1;
    for field in 0..nint + nflt {
        let float = field >= nint;
        // Find the first character of the number.
        loop {
            idx += 1;
            let c = at(idx);
            if starts_number(c) {
                break;
            }
            if c == 0 {
                return Ok(card);
            }
        }
        let start = idx;
        while !is_separator(at(idx)) {
            let c = at(idx);
            if !(starts_number(c) || matches!(c, b'E' | b'e' | b'D' | b'd')) {
                let field_kind = if float { "FLOAT" } else { "INTEGER" };
                return Err(format!(
                    "\n  {kind} DATA CARD \"{}\" ERROR:\n  NON-NUMERICAL CHARACTER '{}' IN {field_kind} FIELD AT CHAR. {}\n",
                    card.mnemonic,
                    c as char,
                    idx + 1
                ));
            }
            idx += 1;
        }
        let value = atof(&bytes[start..idx]);
        if float {
            // nec2c reads a number too large for the machine as infinity and
            // goes on to print values that are not numbers.
            if !value.is_finite() {
                return Err(format!(
                    "\n  {kind} DATA CARD \"{}\" ERROR:\n  VALUE {} IN FIELD {} IS TOO LARGE FOR THE MACHINE\n",
                    card.mnemonic,
                    &line[start..idx],
                    field + 1
                ));
            }
            card.floats[field - nint] = value;
        } else {
            // An integer may be written as a real number, but it is never truncated.
            if value.fract() != 0.0 || value.abs() > f64::from(i32::MAX) {
                return Err(format!(
                    "\n  {kind} DATA CARD \"{}\" ERROR:\n  VALUE {} IN INTEGER FIELD {} IS NOT AN INTEGER\n",
                    card.mnemonic,
                    &line[start..idx],
                    field + 1
                ));
            }
            card.ints[field] = value as i32;
        }
        if at(idx) == 0 {
            return Ok(card);
        }
    }
    Ok(card)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn free_format_fields() {
        let card = parse_card("GW 1 21 0 0 -0.25 0 0 .25 1e-3", 2, 7, "GEOMETRY").unwrap();
        assert_eq!(card.mnemonic, "GW");
        assert_eq!(card.ints, vec![1, 21]);
        assert_eq!(card.floats, vec![0.0, 0.0, -0.25, 0.0, 0.0, 0.25, 0.001]);
    }

    #[test]
    fn commas_and_missing_fields() {
        let card = parse_card("EX,0,1,11,0,1.5", 4, 6, "COMMAND").unwrap();
        assert_eq!(card.ints, vec![0, 1, 11, 0]);
        assert_eq!(card.floats, vec![1.5, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let card = parse_card("XQ", 4, 6, "COMMAND").unwrap();
        assert_eq!(card.ints, vec![0; 4]);
    }

    #[test]
    fn letters_and_fractions_in_integers_are_errors() {
        assert!(parse_card("FR 0 1 0 0 14.2MHz", 4, 6, "COMMAND").is_err());
        assert!(parse_card("GW 1 2.5 0 0 0 0 0 1 0.001", 2, 7, "GEOMETRY").is_err());
        assert!(parse_card("KH 0.3", 4, 6, "COMMAND").unwrap_err().contains("INTEGER FIELD 1"));
    }

    #[test]
    fn integral_reals_and_d_exponents() {
        let card = parse_card("GW 1. 9.0 0 0 -0.25 0 0 0.25 1.0D-3", 2, 7, "GEOMETRY").unwrap();
        assert_eq!(card.ints, vec![1, 9]);
        assert_eq!(card.floats[6], 0.001);
        assert_eq!(parse_card("RP 0 1E1", 4, 6, "COMMAND").unwrap().ints[1], 10);
    }

    #[test]
    fn numbers_beyond_the_range_of_the_machine() {
        let error = parse_card("RP 0 4 2 1000 0 1e400 30 90", 4, 6, "COMMAND").unwrap_err();
        assert!(error.contains("VALUE 1e400 IN FIELD 6 IS TOO LARGE FOR THE MACHINE"), "{error}");
        assert!(parse_card("LD 5 0 0 0 -1E400", 4, 6, "COMMAND").is_err());
        assert!(parse_card("GW 1 1e400 0 0 0 0 0 1 0.001", 2, 7, "GEOMETRY").is_err());
        assert_eq!(parse_card("FR 0 1 0 0 1e-400 0", 4, 6, "COMMAND").unwrap().floats[0], 0.0);
    }

    #[test]
    fn inline_comments_are_dropped() {
        let card = parse_card("EX 0 1 5 ! feed 2", 4, 6, "COMMAND").unwrap();
        assert_eq!(card.ints, vec![0, 1, 5, 0]);
    }

    #[test]
    fn comment_blank_and_indented_lines() {
        let mut lines = Lines::new("# note\n\n' other note\n  gw 1 2\r\nEN");
        assert_eq!(lines.next_line(), Some(Line { number: 4, text: "GW 1 2".to_string() }));
        assert_eq!(lines.next_line().map(|line| line.text), Some("EN".to_string()));
        assert_eq!(lines.next_line(), None);
    }

    #[test]
    fn atof_reads_the_longest_number() {
        assert_eq!(atof(b"1.5E"), 1.5);
        assert_eq!(atof(b"-.5,"), -0.5);
        assert_eq!(atof(b"3.72e7 0"), 3.72e7);
        assert_eq!(atof(b"-"), 0.0);
    }
}
