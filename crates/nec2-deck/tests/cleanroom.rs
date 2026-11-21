//! Tests for the items of the clean-room specification (CR-nnn): behaviour
//! that nec2c 1.3.1 does not have. The decks and the expected numbers come
//! from the specification, where they were measured on running programs.

mod support;

use num_complex::Complex64;
use support::{num, tables, values, warnings};

/// Relative tolerance of the specification for impedances and currents.
const REL: f64 = 2.0e-4;

fn run(deck: &str) -> nec2_deck::Output {
    nec2_deck::run(deck)
}

fn run_ok(deck: &str) -> String {
    let output = run(deck);
    assert_eq!(output.error, None, "the run stopped with an error");
    output.report
}

/// The error of a deck that must be rejected before anything is solved.
fn rejected(deck: &str) -> String {
    let output = run(deck);
    let error = output.error.expect("the deck must be rejected");
    assert!(!output.report.contains("ANTENNA INPUT PARAMETERS"), "rejected only after a solution: {error}");
    error
}

/// Input impedances of all sources, in the order of the report.
fn impedances(report: &str) -> Vec<Complex64> {
    tables(report, "ANTENNA INPUT PARAMETERS")
        .iter()
        .flatten()
        .map(|row| Complex64::new(num(row[6]), num(row[7])))
        .collect()
}

/// Segment currents of the `occurrence`-th current table, by segment number.
fn currents(report: &str, occurrence: usize) -> Vec<Complex64> {
    tables(report, "CURRENTS AND LOCATION")[occurrence]
        .iter()
        .map(|row| Complex64::new(num(row[6]), num(row[7])))
        .collect()
}

/// Total gain in dB at the given angles in the first radiation pattern.
fn total_gain(report: &str, theta: f64, phi: f64) -> f64 {
    tables(report, "RADIATION PATTERNS")[0]
        .iter()
        .find(|row| num(row[0]) == theta && num(row[1]) == phi)
        .map(|row| num(row[4]))
        .expect("pattern row")
}

#[track_caller]
fn assert_close(actual: Complex64, re: f64, im: f64, rel: f64) {
    let expected = Complex64::new(re, im);
    assert!((actual - expected).norm() <= rel * expected.norm(), "{actual} instead of {expected}");
}

#[track_caller]
fn assert_near(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual} instead of {expected}");
}

const DIPOLE_9: &str = "\
CM Nine-segment dipole
CE
GW 1 9 0 0 -0.25 0 0 0.25 0.001
GE 0
EX 0 1 5 0 1 0
FR 0 1 0 0 299.8 0
";

// --- 1. Input parsing ---

#[test]
fn cr001_last_line_without_newline() {
    let report = run_ok(&format!("{DIPOLE_9}XQ 0"));
    assert_close(impedances(&report)[0], 83.252, 46.785, 1.0e-4);

    let report = run_ok(&format!("{DIPOLE_9}RP 0 3 1 1000 0 0 45 0"));
    assert_near(total_gain(&report, 45.0, 0.0), -1.96, 0.01);
    assert_near(total_gain(&report, 90.0, 0.0), 2.16, 0.01);
}

#[test]
fn cr002_end_of_file_ends_the_job() {
    let report = run_ok(&format!("{DIPOLE_9}XQ 0\n"));
    assert_close(impedances(&report)[0], 83.252, 46.785, 1.0e-4);

    let error = rejected("CM No GE card\nCE\nGW 1 9 0 0 -0.25 0 0 0.25 0.001\n");
    assert!(error.contains("GEOMETRY"), "{error}");
}

#[test]
fn cr003_indented_cards_are_read() {
    let report = run_ok(
        "\
CM Probe: second wire card is indented by one blank
CE
GW 1 9 0 0 -0.25 0 0 0.25 0.001
 GW 2 9 0.1 0 -0.25 0.1 0 0.25 0.001
GE 0
EX 0 1 5 0 1 0
FR 0 1 0 0 299.8 0
XQ 0
EN
",
    );
    assert!(report.contains("TOTAL SEGMENTS USED: 18"));
    assert_close(impedances(&report)[0], 25.458, 68.574, 1.0e-4);
}

#[test]
fn cr007_segment_references_are_checked() {
    let with = |card: &str| {
        format!(
            "CM Probe\nCE\nGW 1 9 0 0 -0.25 0 0 0.25 0.001\nGE 0\nEX 0 1 5 0 1 0\n{card}\nFR 0 1 0 0 299.8 0\nXQ 0\nEN\n"
        )
    };
    // A source or a line on a segment that does not exist: NEC-2 prints
    // numbers that are not finite.
    let error = rejected(
        "CM Probe\nCE\nGW 1 9 0 0 -0.25 0 0 0.25 0.001\nGE 0\nEX 0 0 15 0 1 0\nFR 0 1 0 0 299.8 0\nXQ 0\nEN\n",
    );
    assert!(error.contains("EX 0 0 15"), "{error}");
    let error = rejected(&with("TL 1 5 0 30 50 0 0 0 0 0"));
    assert!(error.contains("TL 1 5 0 30"), "{error}");

    // A coupling segment that does not exist: NEC-2 stops when it computes
    // the coupling, after the solution. Here a warning comes with the card.
    let output = run(&with("CP 1 5 0 30"));
    assert!(warnings(&output.report).iter().any(|line| line.contains("NO SEGMENT 30 OF TAG 0")));
    assert!(output.error.expect("the coupling cannot be computed").contains("SEGMENT NUMBER 30"));

    // A range that selects no segment: NEC-2 goes on as if the card were
    // not there, without a word. The result is the same, with a warning.
    let alone = impedances(&run_ok(&with("")));
    for card in ["LD 0 0 12 14 50 0 0", "LD 0 1 12 14 50 0 0", "PT 0 1 12 14", "PQ 0 0 12 14"] {
        let report = run_ok(&with(card));
        let expected = "SEGMENTS 12 TO 14 ARE NOT AMONG THE 9 SEGMENTS";
        assert!(warnings(&report).iter().any(|line| line.contains(expected)), "{card}: {:?}", warnings(&report));
        assert_eq!(impedances(&report), alone, "{card}");
    }

    // A range that reaches beyond the last segment selects those that exist.
    let report = run_ok(&with("LD 4 1 5 99 50 0"));
    assert!(warnings(&report).is_empty(), "{:?}", warnings(&report));
    assert_eq!(impedances(&report), impedances(&run_ok(&with("LD 4 1 5 9 50 0"))));
    assert_ne!(impedances(&report), alone);
}

#[test]
fn cr009_blank_incidence_counts_mean_one() {
    let deck = |counts: &str| {
        format!(
            "CM Probe: plane wave and normalized receiving pattern print\nCE\n\
             GW 1 9 0 0 -0.25 0 0 0.25 0.001\nGE 0\nEX 1 {counts} 0 45 0 0\n\
             FR 0 1 0 0 299.8 0\nPT 2 1 5 5\nXQ 0\nEN\n"
        )
    };
    let blank = run_ok(&deck("0 0"));
    let row = &tables(&blank, "RECEIVING PATTERN PARAMETERS")[0][0];
    assert_near(num(row[2]), 2.2223e-3, REL * 2.2223e-3);
    assert_near(num(row[3]), 147.61, 0.01);
    let normalized = &tables(&blank, "NORMALIZED RECEIVING PATTERN")[0];
    assert_eq!(normalized.len(), 1);
    assert_near(num(normalized[0][2]), 0.0, 0.005);
    assert_eq!(blank, run_ok(&deck("1 1")).replace("EX   1     1     1", "EX   1     0     0"));
}

// --- 2. Geometry ---

#[test]
fn cr020_segment_counts_and_tags() {
    let deck = |card: &str| {
        format!(
            "CM Probe\nCE\nGW 1 9 0 0 -0.25 0 0 0.25 0.001\n{card}\nGE 0\nEX 0 1 5 0 1 0\nFR 0 1 0 0 299.8 0\nXQ 0\nEN\n"
        )
    };
    // Counts with which NEC-2 damages the structure or fails.
    for card in [
        "GW 2 -4 0.1 0 -0.25 0.1 0 0.25 0.001",
        "GH 2 -3 0.1 0.2 0.05 0.05 0.05 0.05 0.001",
        "GR 1 -2",
        "GM 0 -1 0 0 0 0.3 0 0 0",
        "GM 0 1 0 0 0 0.3 0 0 -2",
    ] {
        let error = rejected(&deck(card));
        assert!(error.contains(card), "{card}: {error}");
    }

    // NEC-2 ignores a wire without segments. So does this program, with a warning.
    let alone = impedances(&run_ok(&deck("")));
    for card in ["GW 2 0 0.1 0 -0.25 0.1 0 0.25 0.001", "GA 2 0 0.2 0 90 0.001"] {
        let report = run_ok(&deck(card));
        assert!(warnings(&report).iter().any(|line| line.contains("HAS NO SEGMENTS AND IS IGNORED")), "{card}");
        assert_eq!(impedances(&report), alone, "{card}");
        assert!(report.contains("TOTAL SEGMENTS USED: 9 "), "{card}");
    }

    // A negative tag is a label like any other.
    let report = run_ok(&deck("GW -2 4 0.1 0 -0.25 0.1 0 0.25 0.001"));
    assert!(warnings(&report).is_empty(), "{:?}", warnings(&report));
    assert!(report.contains("TOTAL SEGMENTS USED: 13 "));
}

#[test]
fn cr021_a_deck_without_geometry_is_rejected() {
    let error = rejected("CM Probe: no geometry\nCE\nGE 0\nEX 0 1 1 0 1 0\nXQ 0\nEN\n");
    assert!(error.contains("NO SEGMENTS"), "{error}");
}

#[test]
fn cr023_junction_search_terminates() {
    // Very short wire with tiny segments ending on the end of a wire with long segments.
    let at_wire_end = "\
CM Probe: tiny segments at a wire end
CE
GW 1 10 0 0 0 0 0 1 0.001
GW 2 20 0.0005 0 0 0 0 0 0.001
GE 0
EX 0 1 5 0 1 0
FR 0 1 0 0 30 0
XQ 0
EN
";
    // The same at a junction of two long wires.
    let at_junction = "\
CM Probe: tiny segments at a junction
CE
GW 1 10 0 0 0 0 0 1 0.001
GW 2 10 0 0 0 0 1 0 0.001
GW 3 20 0.0005 0 0 0 0 0 0.001
GE 0
EX 0 1 5 0 1 0
FR 0 1 0 0 30 0
XQ 0
EN
";
    for deck in [at_wire_end, at_junction] {
        let started = std::time::Instant::now();
        let output = run(deck);
        assert!(started.elapsed().as_secs_f64() < 1.0);
        // Either a geometry error, or a solution with a consistent junction table.
        if let Some(error) = output.error {
            assert!(error.contains("CONNECTION TOLERANCE"), "{error}");
        } else {
            assert!(output.report.contains("ANTENNA INPUT PARAMETERS"));
        }
    }
}

#[test]
fn cr025_helix_with_zero_length_is_a_flat_spiral() {
    let report = run_ok(
        "\
CM Probe: flat spiral, 16 segments, 8 segments per turn, radius 0.02 to 0.10
CE
GH 1 16 8 0 0.02 0.02 0.10 0.10 0.001
GE 0
EX 0 1 1 0 1 0
FR 0 1 0 0 299.8 0
XQ 0
EN
",
    );
    let segments = &tables(&report, "SEGMENTATION DATA")[0];
    let centre = |row: &[&str]| [num(row[1]), num(row[2]), num(row[3])];
    assert_eq!(centre(&segments[0]), [0.0188, 0.0088, 0.0]);
    assert_eq!(centre(&segments[15]), [0.0836, -0.0336, 0.0]);
    // Measured with the constants of another program, hence the wider tolerance.
    assert_close(impedances(&report)[0], 4.117, -2991.2, 1.0e-3);

    let error = rejected(
        "CM Probe: zero turn spacing\nCE\nGH 1 16 0 0.3 0.05 0.05 0.05 0.05 0.001\nGE 0\nEX 0 1 1 0 1 0\nXQ 0\nEN\n",
    );
    assert!(error.contains("GH 1 16 0 0.3"), "{error}");
}

// --- 4. Ground ---

#[test]
fn cr050_sommerfeld_contour_is_completed() {
    let report = run_ok(
        "\
CM Probe: 1.25-wavelength horizontal wire, 0.3 wavelength high, Sommerfeld ground, 3.5 MHz
CE
GW 1 41 -53.5 0 25.7 53.5 0 25.7 0.001
GE 0
GN 2 0 0 0 13 0.005
EX 0 1 21 0 1 0
FR 0 1 0 0 3.5 0
RP 0 10 2 1000 0 0 10 90
EN
",
    );
    assert_close(impedances(&report)[0], 213.42, -999.36, REL);

    let report = run_ok(
        "\
CM Probe: two parallel half-wave dipoles 0.6 wavelength apart, 0.3 wavelength high, Sommerfeld ground, 3.5 MHz
CE
GW 1 15 -21 0 25.7 21 0 25.7 0.001
GW 2 15 -21 51.4 25.7 21 51.4 25.7 0.001
GE 0
GN 2 0 0 0 13 0.005
EX 0 1 8 0 1 0
FR 0 1 0 0 3.5 0
RP 0 10 2 1000 0 0 10 90
EN
",
    );
    assert_close(impedances(&report)[0], 96.744, 4.740, REL);
    assert_near(total_gain(&report, 20.0, 90.0), 3.86, 0.02);

    let report = run_ok(
        "\
CM Probe: horizontal half-wave dipole 0.05 wavelength above sea water, Sommerfeld ground, 7 MHz
CE
GW 1 21 -10.4 0 2.14 10.4 0 2.14 0.001
GE 0
GN 2 0 0 0 80 5
EX 0 1 11 0 1 0
FR 0 1 0 0 7 0
RP 0 10 2 1000 0 0 10 90
EN
",
    );
    assert_close(impedances(&report)[0], 7.7637, -5.8325, REL);
    assert_near(total_gain(&report, 0.0, 0.0), 7.67, 0.02);
}

#[test]
fn cr051_ground_integration_limit_does_not_depend_on_the_model_size() {
    let report = run_ok(
        "\
CM Probe: 2-segment wire 0.2 wavelength long, 0.0005 wavelength above Sommerfeld ground
CE
GW 1 2 -8.566 0 0.0428 8.566 0 0.0428 0.001
GE 0
GN 2 0 0 0 13 0.005
EX 0 1 1 0 1 0
FR 0 1 0 0 3.5 0
XQ 0
EN
",
    );
    assert!(!report.contains("STEP SIZE LIMITED"));
    assert_close(impedances(&report)[0], 127.96, -1332.4, REL);
}

// --- 5. Loads and networks ---

#[test]
fn cr060_wire_conductivity_uses_the_formula_for_its_argument() {
    let report = run_ok(
        "\
CM Probe: half-wave dipole of stainless steel wire radius 0.2 mm, LD 5, 7 MHz
CE
GW 1 21 0 0 -10.4 0 0 10.4 0.0002
GE 0
LD 5 1 0 0 1.4e6
EX 0 1 11 0 1 0
FR 0 1 0 0 7 0
XQ 0
EN
",
    );
    assert_close(impedances(&report)[0], 135.19, 6.883, REL);
    assert_near(values(&report, "EFFICIENCY")[0], 53.06, 0.01);

    let report = run_ok(
        "\
CM Probe: 80 m half-wave dipole of thin copper wire (radius 0.1 mm), LD 5, 3.5 MHz
CE
GW 1 21 0 0 -21.2 0 0 21.2 0.0001
GE 0
LD 5 1 0 0 5.8e7
EX 0 1 11 0 1 0
FR 0 1 0 0 3.5 0
XQ 0
EN
",
    );
    assert_close(impedances(&report)[0], 95.278, 37.174, REL);
    assert_near(values(&report, "EFFICIENCY")[0], 78.78, 0.01);
}

#[test]
fn cr061_source_outside_the_network_is_not_network_loss() {
    let report = run_ok(
        "\
CM Probe: network with one source on a network port and one source off the network, full-sphere average gain
CE
GW 1 9 0 0 -0.25 0 0 0.25 0.001
GW 2 9 0.3 0 -0.25 0.3 0 0.25 0.001
GW 3 1 0.6 0 -0.01 0.6 0 0.01 0.001
GE 0
NT 1 5 3 1 0.02 -0.01 0.001 0.03 0.015 0.002
NT 2 5 3 1 0.01 0 0 0.02 0.01 0
EX 0 3 1 0 1 0
EX 0 2 2 0 0.5 0.5
FR 0 1 0 0 299.8 0
RP 0 37 73 1002 0 0 5 5
EN
",
    );
    assert_near(values(&report, "NETWORK LOSS")[0], 2.9054e-2, REL * 2.9054e-2);
    assert_near(values(&report, "RADIATED POWER")[0], 8.0054e-3, REL * 8.0054e-3);
    assert_near(values(&report, "EFFICIENCY")[0], 21.60, 0.01);
    let z = impedances(&report);
    assert_close(z[0], 12.991, -1.2609, REL);
    assert_close(z[1], -160.45, 108.25, REL);
    assert_near(values(&report, "AVERAGE POWER GAIN")[0], 0.2151, 0.002);
}

// --- 6. Excitation ---

#[test]
fn cr070_current_source_on_a_segment_centre() {
    let deck = |first_voltage: &str| {
        format!(
            "CM Probe: voltage-source solution first, then EX 4 source exactly at the centre of segment 2\nCE\n\
             GW 1 3 0 0 -0.15 0 0 0.15 0.001\nGW 2 3 0.2 0 -0.15 0.2 0 0.15 0.001\nGE 0\n\
             EX 0 1 2 0 {first_voltage} 0\nFR 0 1 0 0 299.8 0\nXQ 0\n\
             EX 4 0 0 0 0 0 0 90 0 0.001\nXQ 0\nEN\n"
        )
    };
    let report = run_ok(&deck("1000"));
    let second = currents(&report, 1);
    assert_close(second[0], 1.7504e-3, -2.2178e-4, REL);
    assert_close(second[1], 2.5435e-3, -3.5571e-4, REL);
    assert_close(second[4], -3.1888e-4, -3.7971e-4, REL);
    // The result does not depend on the excitation before it.
    assert_eq!(second, currents(&run_ok(&deck("1")), 1));
}

#[test]
fn cr072_extended_kernel_flags_of_end_2() {
    let report = run_ok(
        "\
CM Probe: two fat wires meeting head to head (end 2 to end 2) at a right angle, extended kernel, near E field
CE
GW 1 5 0 0 0 0 0 0.25 0.005
GW 2 5 0.2 0 0.25 0 0 0.25 0.005
GE 0
EK 0
EX 0 1 3 0 1 0
FR 0 1 0 0 299.8 0
NE 0 5 1 5 0.02 0.01 0.20 0.04 0 0.025
EN
",
    );
    let field = &tables(&report, "NEAR ELECTRIC FIELDS")[0];
    let at = |x: f64, z: f64| field.iter().find(|row| num(row[0]) == x && num(row[2]) == z).expect("field point");
    assert_near(num(at(0.06, 0.25)[3]), 0.30738, REL * 0.30738);
    assert_near(num(at(0.02, 0.25)[3]), 0.7575, 2.0 * REL * 0.7575);
    assert_near(num(at(0.02, 0.25)[5]), 30.25, 2.0 * REL * 30.25);

    let report = run_ok(
        "\
CM Probe: two wires head to head, EX 5 current-slope source next to the junction, extended kernel
CE
GW 1 5 0 0 0 0 0 0.25 0.005
GW 2 5 0.2 0 0.25 0 0 0.25 0.005
GE 0
EK 0
EX 5 1 5 0 1 0
FR 0 1 0 0 299.8 0
XQ 0
EN
",
    );
    assert_close(impedances(&report)[0], 39.153, -45.055, REL);
}

#[test]
fn cr073_repeated_execution_after_an_incident_wave_loop() {
    let report = run_ok(
        "\
CM Probe: second and third XQ after an incident-wave loop without a new EX card
CE
GW 1 9 0 0 -0.25 0 0 0.25 0.001
GW 2 9 0.25 0 -0.25 0.25 0 0.25 0.001
GE 0
LD 4 1 5 5 50 0
EX 1 4 3 0 0 0 0 30 45 0
FR 0 1 0 0 299.8 0
PT 1 1 5 5
XQ 0
PT 2 1 5 5
XQ 0
PT 3 1 5 5
XQ 0
EN
",
    );
    // PT 1 and PT 2 list the 12 currents, PT 3 prints the normalized table only.
    let received = tables(&report, "RECEIVING PATTERN PARAMETERS");
    assert_eq!(received.len(), 2);
    assert_eq!(received[0].len(), 12);
    assert_eq!(received[0], received[1]);
    let current = |theta: f64, phi: f64| {
        let row = received[1].iter().find(|row| num(row[0]) == theta && num(row[1]) == phi).expect("angle");
        num(row[2])
    };
    assert_near(current(30.0, 0.0), 4.5297e-4, REL * 4.5297e-4);
    assert_near(current(90.0, 90.0), 2.0033e-3, REL * 2.0033e-3);

    let normalized = tables(&report, "NORMALIZED RECEIVING PATTERN");
    assert_eq!(normalized.len(), 2);
    assert_eq!(normalized[0], normalized[1]);
    assert_eq!(normalized[0].len(), 12);
    assert_eq!(values(&report, "NORMALIZATION FACTOR"), vec![2.0033e-3, 2.0033e-3]);
    assert_near(num(normalized[0][1][2]), -12.91, 0.01);
}

// --- Robustness items ---

fn dipole_with(cards: &str) -> String {
    format!("{DIPOLE_9}{cards}\nEN\n")
}

#[test]
fn cr004_long_lines_are_one_card() {
    let padding = " ".repeat(300);
    let comment = "x".repeat(400);
    let report = run_ok(&format!(
        "CM {comment}\nCE\nGW 1 9 0 0 -0.25 0 0 0.25 0.001{padding}\nGE 0\nEX 0 1 5 0 1 0\nFR 0 1 0 0 299.8 0\nXQ 0\nEN\n"
    ));
    assert_close(impedances(&report)[0], 83.252, 46.785, 1.0e-4);
}

#[test]
fn cr005_comment_styles() {
    let report = run_ok(
        "\
CM Probe: apostrophe comment line, inline exclamation comment, lower case mnemonics
CE
gw 1 9 0 0 -0.25 0 0 0.25 0.001
' a 4nec2 style comment line
ge 0
ex 0 1 5 0 1 0 ! feed in the middle
fr 0 1 0 0 299.8 0
xq 0
en
",
    );
    assert_close(impedances(&report)[0], 83.252, 46.785, 1.0e-4);
}

#[test]
fn cr006_numeric_field_syntax() {
    let report = run_ok(
        "CM Probe\nCE\nGW 1. 9.0 0 0 -0.25 0 0 0.25 1.0D-3\nGE 0\nEX 0 1 5 0 1 0\nFR 0 1 0 0 2.998D2 0\nXQ 0\nEN\n",
    );
    assert_close(impedances(&report)[0], 83.252, 46.785, 1.0e-4);
    let error = rejected(&dipole_with("KH 0.3\nXQ 0"));
    assert!(error.contains("KH 0.3") && error.contains("INTEGER FIELD 1"), "{error}");
}

#[test]
fn cr008_blank_near_field_counts_mean_one() {
    let blank = run_ok(&dipole_with("NE 0 0 0 0 0.1 0 0"));
    let rows = &tables(&blank, "NEAR ELECTRIC FIELDS")[0];
    assert_eq!(rows.len(), 1);
    assert_near(num(rows[0][3]), 5.8080e-8, 1.0e-10);
    assert_near(num(rows[0][7]), 2.6688, REL * 2.6688);
    assert_near(num(rows[0][8]), 146.54, 0.01);
}

#[test]
fn cr010_to_cr013_bad_cards_are_rejected() {
    for card in [
        "FR 0 1 0 0 0 0",
        "FR 0 5 0 0 100 -30",
        "FR 1 3 0 0 100 0",
        "TL 1 5 1 3 0 0",
        "GN 0",
        "GN 0 -4 0 0 13 0.005 1 0.001",
        "GW 2 5 0.1 0 -0.25 0.1 0 0.25 0.001",
        "QQ 1",
    ] {
        let error = rejected(&dipole_with(&format!("{card}\nXQ 0")));
        assert!(error.contains(card), "{card}: {error}");
    }
    // A lossless finite ground is valid.
    run_ok(&dipole_with("GN 0 0 0 0 4 0\nXQ 0").replace("-0.25 0 0 0.25", "0.25 0 0 0.75"));
    // A ground type that NEC-2 does not define is solved the way NEC-2
    // solves it, with a warning.
    let report = run_ok(&dipole_with("GN 3 0 0 0 13 0.005\nXQ 0").replace("-0.25 0 0 0.25", "0.25 0 0 0.75"));
    assert!(warnings(&report).iter().any(|line| line.contains("GROUND TYPE 3 IS NOT DEFINED")));
}

#[test]
fn cr012_ground_flag_and_wires_below_ground() {
    let deck = |wire: &str, ge: &str, gn: &str| {
        format!("CM Probe\nCE\n{wire}\n{ge}\n{gn}EX 0 1 3 0 1 0\nFR 0 1 0 0 299.8 0\nXQ 0\nEN\n")
    };
    let standing = "GW 1 5 0 0 0 0 0 0.25 0.001";
    let sunk = "GW 1 5 0 0 -0.1 0 0 0.25 0.001";
    let lying = "GW 1 5 0 0 0 0.25 0 0 0.001";
    #[track_caller]
    fn warned(deck: &str, text: &str) {
        let report = run_ok(deck);
        assert!(warnings(&report).iter().any(|line| line.contains(text)), "{:?}", warnings(&report));
    }
    // NEC-2 solves these structures, so they are solved, but not without a word.
    warned(&deck(standing, "GE 1", ""), "SEGMENT ENDS JOINED TO THEIR IMAGE BY THE GE CARD: 1,");
    warned(&deck(sunk, "GE -1", "GN 1\n"), "REACH BELOW THE GROUND PLANE OR LIE IN IT: 2,");
    warned(&deck(sunk, "GE 0", "GN 1\n"), "REACH BELOW THE GROUND PLANE OR LIE IN IT: 2,");
    warned(&deck(lying, "GE 0", "GN 1\n"), "REACH BELOW THE GROUND PLANE OR LIE IN IT: 5,");
    // With GE 1 NEC-2 refuses a wire below the ground.
    assert!(rejected(&deck(sunk, "GE 1", "GN 1\n")).contains("BELOW GROUND"));
    // A wire that ends on the ground, and a flag that changes nothing.
    assert!(warnings(&run_ok(&deck(standing, "GE -1", "GN 1\n"))).is_empty());
    assert!(warnings(&run_ok(&deck("GW 1 5 0 0 1 0 0 1.25 0.001", "GE 1", ""))).is_empty());
}

#[test]
fn cr022_rotation_by_one_sector_changes_nothing() {
    let deck = |card: &str| {
        format!("CM Probe\nCE\nGW 1 9 0.1 0 -0.25 0.1 0 0.25 0.001\n{card}GE 0\nEX 0 1 5 0 1 0\nXQ 0\nEN\n")
    };
    let plain = run_ok(&deck(""));
    let rotated = run_ok(&deck("GR 0 1\n"));
    assert_eq!(impedances(&plain), impedances(&rotated));
    assert!(rotated.contains("SYMMETRY FLAG: 0"));
    assert!(rejected(&deck("GR 0 0\n")).contains("GR 0 0"));
}

#[test]
fn cr053_screen_pattern_without_a_screen() {
    // NEC-2 computes such a pattern without the screen, but prints
    // not-a-number where the image point of a segment lies on the axis.
    // Here every row is that of the plain ground.
    let deck = |mode: i32| {
        format!(
            "CM Probe\nCE\nGW 1 5 0 0 0 0 0 0.25 0.001\nGE 1\nGN 1\nEX 0 1 1 0 1 0\nFR 0 1 0 0 299.8 0\nRP {mode} 4 1 1000 0 0 30 0\nEN\n"
        )
    };
    let rows = |report: &str| -> Vec<Vec<String>> {
        tables(report, "RADIATION PATTERNS")[0].iter().map(|row| row.iter().map(|t| t.to_string()).collect()).collect()
    };
    let screen_mode = run_ok(&deck(4));
    assert!(warnings(&screen_mode).iter().any(|line| line.contains("COMPUTED WITHOUT A SCREEN")));
    assert_eq!(rows(&screen_mode).len(), 4);
    assert_eq!(rows(&screen_mode), rows(&run_ok(&deck(0))));
    assert!(!screen_mode.to_uppercase().contains("NAN"));
}

#[test]
fn cr053_a_screen_does_not_outlive_its_ground_card() {
    let monopole = "CM Probe\nCE\nGW 1 5 0 0 0 0 0 0.25 0.001\nGE 1\n";
    let screen = "GN 0 12 0 0 13 0.005 2 0.001\nEX 0 1 1 0 1 0\nFR 0 1 0 0 299.8 0\nRP 4 4 1 1000 0 0 30 0\n";
    let perfect = "GN 1\nEX 0 1 1 0 1 0\nFR 0 1 0 0 299.8 0\nRP 0 4 1 1000 0 0 30 0\nNE 0 1 1 3 0.2 0 0.1 0 0 0.1\n";
    // The last pattern and near-field tables of a report.
    let last = |report: &str| -> Vec<Vec<String>> {
        ["RADIATION PATTERNS", "NEAR ELECTRIC FIELDS"]
            .iter()
            .flat_map(|title| tables(report, title).pop().expect("table"))
            .map(|row| row.iter().map(|token| token.to_string()).collect())
            .collect()
    };

    let alone = last(&run_ok(&format!("{monopole}{perfect}EN\n")));
    assert_eq!(alone.len(), 4 + 3);
    // A perfect ground after a screen: in the next job and in the same job.
    assert_eq!(last(&run_ok(&format!("{monopole}{screen}NX\n{monopole}{perfect}EN\n"))), alone);
    assert_eq!(last(&run_ok(&format!("{monopole}{screen}{perfect}EN\n"))), alone);
}

// --- Surface patches ---

/// The rows of the patch-current table: number, centre, the two tangential
/// currents (magnitude, phase) and the x, y, z components.
fn patch_currents(report: &str) -> Vec<Vec<f64>> {
    tables(report, "SURFACE PATCH CURRENTS")[0].iter().map(|row| row.iter().map(|token| num(token)).collect()).collect()
}

#[test]
fn cr032_a_strip_of_patches_goes_on_with_sc_cards() {
    let report = run_ok(
        "CM Probe: strip of three quadrilateral patches entered with one SP card and three SC cards, plane wave \
         excitation\nCE\nSP 0 3 0 0 0 0.1 0 0\nSC 0 3 0.1 0.1 0 0 0.1 0\nSC 0 3 0.1 0.2 0 0 0.2 0\n\
         SC 0 3 0.1 0.3 0 0 0.3 0\nGE 0\nEX 1 1 1 0 45 0 0\nFR 0 1 0 0 299.8 0\nXQ 0\nEN\n",
    );
    let patches = &tables(&report, "SURFACE PATCH DATA")[0];
    assert_eq!(patches.len(), 3);
    for (row, y) in patches.iter().zip([0.05, 0.15, 0.25]) {
        assert_eq!((num(row[1]), num(row[2]), num(row[3])), (0.05, y, 0.0));
        assert_eq!(num(row[7]), 0.01);
    }
    let currents = patch_currents(&report);
    assert_eq!(currents.len(), 3);
    for row in &currents {
        assert_near(row[4], 5.3088e-3, 5.3088e-3 * REL);
        assert_near(row[5], 12.73, 0.01);
    }
    // An SC card that follows nothing it could add to stays an error.
    assert!(rejected("CM Probe\nCE\nSC 0 1 0.1 0.1 0\nGE 0\nEN\n").contains("PATCH DATA ERROR"));
}

#[test]
fn cr071_current_source_on_the_centre_of_a_patch() {
    let report = run_ok(
        "CM Probe\nCE\nSP 0 0 -0.2 0 0 90 0 0.01\nSP 0 0 0 0 0 90 0 0.01\nSP 0 0 0.2 0 0 90 0 0.01\nGE 0\n\
         EX 4 0 0 0 0 0 0 90 0 1\nFR 0 1 0 0 299.8 0\nXQ 0\nEN\n",
    );
    let currents = patch_currents(&report);
    // The patch under the source gets no field. The two others are mirror
    // images: the same current in opposite directions. (In nec2c the third
    // patch is computed with the geometry of the second and stays empty.)
    assert_eq!(currents[1][4], 0.0);
    assert!(currents[0][4] > 1.0);
    assert_near(currents[2][4], currents[0][4], currents[0][4] * 1.0e-4);
    assert_near((currents[2][5] - currents[0][5]).abs(), 180.0, 0.01);
}

#[test]
fn cr026_patch_data_must_be_finite() {
    // Corners so far apart that the sides of the patch are not finite.
    let error = rejected(
        "CM Probe\nCE\nSP 0 1 -1e308 0 0 1e308 0 0\nSC 0 0 1e308 1e308 0\nGE 0\nEX 1 1 1 0 45 0 0\nXQ 0\nEN\n",
    );
    assert!(error.contains("PATCH 1"), "{error}");
    // A number that the machine cannot hold is refused where it stands.
    let error = rejected("CM Probe\nCE\nSP 0 0 1e400 0 0 90 0 0.01\nGE 0\nEN\n");
    assert!(error.contains("VALUE 1e400 IN FIELD 3 IS TOO LARGE FOR THE MACHINE"), "{error}");
}

/// NEC-2 divides the patch under a wire end into four and looks at the four
/// again. A patch without area, or one far smaller than the segment, lies
/// under the wire end after every division, and NEC-2 does not end.
#[test]
fn a_patch_too_small_for_its_wire_is_an_error() {
    for area in ["0", "1e-300", "-0.04"] {
        let error = rejected(&format!(
            "CM Probe\nCE\nSP 0 0 0 0 0 90 0 {area}\nGW 1 3 0 0 0 0 0 0.2 0.002\nGE 0\nEX 0 1 1 0 1 0\nXQ 0\nEN\n"
        ));
        assert!(error.contains("TOO SMALL FOR THE END OF SEGMENT 1"), "{area}: {error}");
    }
    // Two wire ends on one patch: NEC-2 reads outside its arrays.
    let error = rejected(
        "CM Probe\nCE\nSP 0 0 0 0 0 90 0 0.04\nGW 1 3 0 0 0 0 0 0.25 0.002\nGW 2 3 0 0 0 0.2 0 0.15 0.002\nGE 0\n\
         EX 0 1 1 0 1 0\nXQ 0\nEN\n",
    );
    assert!(error.contains("WITHOUT FORMING A CLEAN JUNCTION"), "{error}");
}

/// A caller can limit the size of the structures it accepts.
#[test]
fn a_structure_above_the_size_limit_is_refused() {
    use nec2_deck::{Options, run_with};
    let deck = format!("{DIPOLE_9}XQ 0\nEN\n");
    let limit = |equations| Options { max_equations: Some(equations), ..Options::default() };
    assert_eq!(run_with(&deck, &limit(9)).error, None);
    let error = run_with(&deck, &limit(8)).error.expect("nine segments are too many");
    assert!(error.contains("NEEDS 9 EQUATIONS, MORE THAN THE LIMIT OF 8"), "{error}");
    // Two equations for every patch, and a repetition count is caught at its card.
    let patches = "CM Probe\nCE\nSM 2 2 0 0 0 0.2 0 0\nSC 0 0 0.2 0.2 0\nGM 0 360 0 0 30\nGE 0\nEN\n";
    let error = run_with(patches, &limit(100)).error.expect("too many patches");
    assert!(error.contains("GM 0 360") && error.contains("NEEDS 2888 EQUATIONS"), "{error}");
}

// --- The remaining robustness items ---

#[test]
fn cr026_segment_data_must_be_finite() {
    // A length that is not finite, no length, no radius.
    for wire in ["GW 1 3 0 0 -1e308 0 0 1e308 0.001", "GW 1 3 0 0 0 0 0 0 0.001", "GW 1 3 0 0 0 0 0 1 -0.001"] {
        let error = rejected(&format!("CM Probe\nCE\n{wire}\nGE 0\nEX 0 1 2 0 1 0\nXQ 0\nEN\n"));
        assert!(error.contains("SEGMENT 1") && error.contains("TAG 1"), "{wire}: {error}");
    }
}

const ACUTE_JUNCTION: &str = "\
CM Probe: acute junction, first segment centre of wire 2 inside the volume of fat wire 1
CE
GW 1 5 0 0 0 0 0 0.5 0.02
GW 2 10 0 0 0.5 0.03 0 0.2 0.001
GE 0
EX 0 1 3 0 1 0
FR 0 1 0 0 299.8 0
XQ 0
EN
";

#[test]
fn cr027_wires_in_the_same_space_are_reported() {
    use nec2_deck::{Options, OverlapCheck, run_with};
    let strict = Options { overlap: OverlapCheck::Error, ..Options::default() };
    let two_wires = |second: &str| {
        format!(
            "CM Probe\nCE\nGW 1 9 0 0 -0.25 0 0 0.25 0.01\n{second}\nGE 0\nEX 0 1 5 0 1 0\nFR 0 1 0 0 299.8 0\nXQ 0\nEN\n"
        )
    };
    // The same wire twice, two wires that cross without a junction, two
    // parallel wires closer than the radius of the fat one.
    for second in
        ["GW 2 9 0 0 -0.25 0 0 0.25 0.01", "GW 2 9 -0.25 0 0 0.25 0 0 0.001", "GW 2 9 0.005 0 -0.25 0.005 0 0.25 0.001"]
    {
        // NEC-2 solves the structure, and so does this program, with a warning.
        let report = run_ok(&two_wires(second));
        assert!(warnings(&report).iter().any(|line| line.contains("INSIDE THE VOLUME")), "{second}");
        // On request the structure is refused before anything is solved.
        let output = run_with(&two_wires(second), &strict);
        let error = output.error.expect("the structure must be refused");
        assert!(error.contains("INSIDE THE VOLUME"), "{second}: {error}");
        assert!(!output.report.contains("ANTENNA INPUT PARAMETERS"), "{second}");
    }

    let error = run_with(ACUTE_JUNCTION, &strict).error.expect("the structure must be refused");
    let pair =
        error.lines().find(|line| line.contains("SEGMENT 12 ") && line.contains("SEGMENT 4 ")).expect("pair 12 and 4");
    assert!(pair.contains("TAG 2") && pair.contains("TAG 1"), "{pair}");
    assert!(pair.contains("DISTANCE 1.9500E-2") && pair.contains("RADIUS 2.0000E-2"), "{pair}");
}

#[test]
fn cr027_the_overlap_test_can_warn_or_be_switched_off() {
    use nec2_deck::{Options, OverlapCheck, run_with};
    assert_eq!(Options::default().overlap, OverlapCheck::Warn);
    let warned = run_with(ACUTE_JUNCTION, &Options { overlap: OverlapCheck::Warn, ..Options::default() });
    assert_eq!(warned.error, None);
    assert!(warned.report.contains("WARNING") && warned.report.contains("INSIDE THE VOLUME"));
    assert_eq!(impedances(&warned.report).len(), 1);

    let unchecked = run_with(ACUTE_JUNCTION, &Options { overlap: OverlapCheck::Off, ..Options::default() });
    assert_eq!(unchecked.error, None);
    assert!(!unchecked.report.contains("INSIDE THE VOLUME"));
    assert_eq!(impedances(&unchecked.report), impedances(&warned.report));
}

#[test]
fn cr054_negative_conductivity_holds_for_every_frequency() {
    let deck = |frequency_card: &str| {
        format!(
            "CM Probe: imaginary part of the dielectric constant given as a negative conductivity\nCE\n\
             GW 1 9 0 -0.25 0.3 0 0.25 0.3 0.001\nGE 0\nGN 0 0 0 0 12 -3\nEX 0 1 5 0 1 0\n{frequency_card}\nXQ 0\nEN\n"
        )
    };
    let dielectric = |report: &str| -> Vec<String> {
        report
            .lines()
            .filter_map(|line| line.split_once("COMPLEX DIELECTRIC CONSTANT:"))
            .map(|(_, value)| value.trim().to_string())
            .collect()
    };
    let rows = |report: &str| -> Vec<Vec<String>> {
        tables(report, "ANTENNA INPUT PARAMETERS")
            .concat()
            .into_iter()
            .map(|row| row.into_iter().map(str::to_string).collect())
            .collect()
    };

    let sweep = run_ok(&deck("FR 0 3 0 0 250 50"));
    assert_eq!(dielectric(&sweep), vec!["1.2000E+01-3.0000E+00j"; 3]);
    let sweep_rows = rows(&sweep);
    assert_eq!(sweep_rows.len(), 3);
    for (k, frequency) in [250, 300, 350].into_iter().enumerate() {
        let single = run_ok(&deck(&format!("FR 0 1 0 0 {frequency} 0")));
        assert_eq!(rows(&single), vec![sweep_rows[k].clone()], "{frequency} MHz");
    }
}

/// The ratio of the smallest to the largest pivot printed with the matrix timing.
fn pivot_ratio(report: &str) -> f64 {
    values(report, "SMALLEST/LARGEST PIVOT")[0]
}

#[test]
fn cr080_ill_conditioned_matrices_are_reported() {
    use nec2_deck::{Options, OverlapCheck, run_with};
    // A healthy structure: no warning.
    let healthy = run_ok(&format!("{DIPOLE_9}XQ 0\nEN\n"));
    assert!(pivot_ratio(&healthy) > 0.05);
    assert!(!healthy.contains("ILL-CONDITIONED"));

    // The same wire entered twice, with the overlap test switched off.
    let twice = "CM Probe: coincident wires\nCE\nGW 1 9 0 0 -0.25 0 0 0.25 0.001\nGW 2 9 0 0 -0.25 0 0 0.25 0.001\n\
                 GE 0\nEX 0 1 5 0 1 0\nFR 0 1 0 0 299.8 0\nXQ 0\nEN\n";
    let output = run_with(twice, &Options { overlap: OverlapCheck::Off, ..Options::default() });
    assert_eq!(output.error, None);
    assert!(output.report.contains("ILL-CONDITIONED"));
    let ratio = pivot_ratio(&output.report);
    assert!((ratio / 1.7e-6 - 1.0).abs() < 0.2, "{ratio}");
}

/// Average power gain and the solid angle in units of pi, as printed.
fn average_gain(report: &str) -> (f64, f64) {
    let line = report.lines().find(|line| line.contains("AVERAGE POWER GAIN")).expect("average gain line");
    let gain = line.split("AVERAGE POWER GAIN:").nth(1).and_then(|rest| rest.split_whitespace().next()).expect("gain");
    let solid_angle = line.split('(').nth(1).and_then(|rest| rest.split(')').next()).expect("solid angle");
    (gain.parse().expect("gain is a number"), solid_angle.parse().expect("solid angle is a number"))
}

#[test]
fn cr090_average_gain_of_a_cut_through_the_zenith() {
    let (gain, solid_angle) = average_gain(&run_ok(&format!("{DIPOLE_9}RP 0 19 2 1002 -90 0 10 90\nEN\n")));
    assert_near(solid_angle, 1.0, 5.0e-5);
    assert_near(gain, 0.994, 0.003);

    // A full sphere gives what nec2c gives.
    let fine = run_ok(&format!("{DIPOLE_9}RP 0 37 73 1002 0 0 5 5\nEN\n"));
    let (gain, solid_angle) = average_gain(&fine);
    assert_near(solid_angle, 4.0, 5.0e-5);
    assert_near(gain, 0.99533, 1.0e-5);
    assert!(!fine.contains("WARNING"));

    // More power than the antenna takes: the three theta values miss the
    // nulls of the dipole.
    let coarse = run_ok(&format!("{DIPOLE_9}RP 0 3 3 1001 0 0 90 180\nEN\n"));
    assert!(average_gain(&coarse).0 > 1.1);
    assert!(coarse.contains("WARNING: THE AVERAGE POWER GAIN EXCEEDS ITS LIMIT"));

    // No solid angle at all: both phi cuts are the same.
    let flat = run_ok(&format!("{DIPOLE_9}RP 0 19 2 1002 -90 0 10 0\nEN\n"));
    assert!(flat.contains("AVERAGE POWER GAIN: UNDEFINED"));
}

const MONOPOLE_OVER_GROUND: &str = "\
CM Probe: pattern cut through both sides of zenith over perfect ground with normalized gain and average gain
CE
GW 1 9 0 0 0.1 0 0 0.6 0.001
GE 0
GN 1
EX 0 1 5 0 1 0
FR 0 1 0 0 299.8 0
";

#[test]
fn cr091_directions_below_the_ground_are_not_samples() {
    let report = run_ok(&format!("{MONOPOLE_OVER_GROUND}RP 0 7 2 1502 0 0 30 90\nEN\n"));
    let entries: Vec<[f64; 3]> = tables(&report, "NORMALIZED GAIN")[0]
        .iter()
        .flat_map(|row| row.chunks(3).map(|e| [num(e[0]), num(e[1]), num(e[2])]).collect::<Vec<_>>())
        .collect();
    let angles: Vec<[f64; 2]> = entries.iter().map(|e| [e[0], e[1]]).collect();
    let expected: Vec<[f64; 2]> =
        [0.0, 90.0].iter().flat_map(|&phi| [0.0, 30.0, 60.0, 90.0].map(|theta| [theta, phi])).collect();
    let mut sorted = angles.clone();
    sorted.sort_by(|a, b| (a[1], a[0]).partial_cmp(&(b[1], b[0])).expect("angles are numbers"));
    assert_eq!(sorted, expected);
    for entry in &entries {
        if entry[0] == 90.0 {
            assert_near(entry[2], 0.0, 0.005);
        }
        if entry[0] == 60.0 {
            assert_near(entry[2], -8.65, 0.01);
        }
    }
    let (_, solid_angle) = average_gain(&report);
    assert_near(solid_angle, 0.5, 5.0e-5);

    // The whole theta range in 1 degree steps: the average is taken over the upper hemisphere.
    let (gain, solid_angle) = average_gain(&run_ok(&format!("{MONOPOLE_OVER_GROUND}RP 0 181 5 1002 0 0 1 90\nEN\n")));
    assert_near(solid_angle, 2.0, 5.0e-5);
    assert_near(gain, 1.990, 0.005);
}
