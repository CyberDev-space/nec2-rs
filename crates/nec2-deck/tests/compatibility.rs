//! Errors that the compatibility rule asks for (`docs/nec2c-deviations.md`):
//! decks that NEC-2 runs to the end although the result means nothing. A
//! clear error is the better result. The decks that NEC-2 solves with a
//! result that means something are the `s7_*` decks of the fixture tests.

const MONOPOLE: &str = "\
CM Monopole over a ground
CE
GW 1 9 0 0 0.1 0 0 0.6 0.001
GE 0
";

/// The error of a deck and the report up to it.
fn stopped(deck: &str) -> (String, String) {
    let output = nec2_deck::run(deck);
    (output.error.expect("the deck must stop with an error"), output.report)
}

/// NEC-2 computes a pattern mode above 6 as mode 6. nec2c prints a null
/// pointer as the name of the cliff, and the Fortran what the name held last.
#[test]
fn a_pattern_mode_above_6_is_an_error() {
    for mode in [7, 90] {
        let deck = format!("{MONOPOLE}GN 1\nEX 0 1 5 0 1 0\nFR 0 1 0 0 299.8 0\nRP {mode} 4 2 1000 0 0 25 90\nEN\n");
        let (error, report) = stopped(&deck);
        assert!(error.contains("THE PATTERN MODE MUST BE 0 TO 6"), "{error}");
        assert!(!report.contains("RADIATION PATTERNS"));
    }
}

/// NEC-2 solves a structure without a source: the currents are zero, and the
/// efficiency and the gains are not numbers. The tables that do not depend
/// on a source are printed before the error.
#[test]
fn a_structure_without_a_source_is_an_error() {
    for request in ["XQ 0", "RP 0 2 1 1000 0 0 45 0", "NE 0 1 1 1 0.1 0 0.1"] {
        let (error, report) = stopped(&format!("{MONOPOLE}LD 5 0 0 0 5.8E7\nFR 0 1 0 0 299.8 0\n{request}\nEN\n"));
        assert!(error.contains("NO EX CARD DEFINES A SOURCE"), "{request}: {error}");
        assert!(report.contains("STRUCTURE IMPEDANCE LOADING") && report.contains("ANTENNA ENVIRONMENT"), "{request}");
        assert!(!report.contains("CURRENTS AND LOCATION") && !report.contains("POWER BUDGET"), "{request}");
    }
    // A source that a card before the request defines is enough.
    let output = nec2_deck::run(&format!("{MONOPOLE}FR 0 1 0 0 299.8 0\nEX 0 1 5 0 1 0\nXQ 0\nEN\n"));
    assert_eq!(output.error, None);
}

/// The cliff modes reflect from the medium beyond the cliff. Without its
/// parameters NEC-2 prints values that are not numbers.
#[test]
fn a_cliff_without_a_medium_beyond_it_is_an_error() {
    let deck = |ground: &str, mode: i32| {
        format!("{MONOPOLE}{ground}\nEX 0 1 5 0 1 0\nFR 0 1 0 0 299.8 0\nRP {mode} 4 1 1000 0 0 25 0\nEN\n")
    };
    for mode in [2, 3, 5, 6] {
        for ground in ["GN 1", "GN 0 0 0 0 13 0.005", "GN 0 16 0 0 13 0.005 0.5 0.002"] {
            let (error, _) = stopped(&deck(ground, mode));
            assert!(error.contains("THE PATTERN MODE NEEDS THE MEDIUM BEYOND THE CLIFF"), "{ground}, {mode}: {error}");
        }
        // The medium of the GN card, and no ground at all.
        for ground in ["GN 0 0 0 0 13 0.005 80 4 2 0.5", "GN -1"] {
            assert_eq!(nec2_deck::run(&deck(ground, mode)).error, None, "{ground}, {mode}");
        }
    }
}

/// At the origin the ground wave has the factor 1/R with R = 0.
#[test]
fn the_ground_wave_at_the_origin_is_an_error() {
    let deck = |points: &str| {
        format!("{MONOPOLE}GN 0 0 0 0 13 0.005\nEX 0 1 5 0 1 0\nFR 0 1 0 0 299.8 0\nRP 1 {points}\nEN\n")
    };
    let (error, _) = stopped(&deck("3 1 0 0 0 0.5 0"));
    assert!(error.contains("THE GROUND WAVE IS ASKED FOR AT THE ORIGIN"), "{error}");
    // A distance, or heights that leave out the origin.
    for points in ["3 1 0 0 0 0.5 0 100", "3 1 0 0.5 0 0.5 0"] {
        assert_eq!(nec2_deck::run(&deck(points)).error, None, "{points}");
    }
}

/// A count below zero. NEC-2 computes one frequency, prints tables without
/// rows, or damages what it has.
#[test]
fn a_negative_count_is_an_error() {
    let cards = [
        ("FR 0 -1 0 0 299.8 0\nEX 0 1 5 0 1 0\nXQ 0", "THE NUMBER OF FREQUENCIES IS NEGATIVE"),
        ("EX 0 1 5 0 1 0\nRP 0 -1 2 1000 0 0 30 90", "THE NUMBER OF ANGLES IS NEGATIVE"),
        ("EX 0 1 5 0 1 0\nRP 0 4 -2 1000 0 0 30 90", "THE NUMBER OF ANGLES IS NEGATIVE"),
        ("EX 0 1 5 0 1 0\nNE 0 1 -1 1 0.1 0 0.1", "THE NUMBER OF FIELD POINTS IS NEGATIVE"),
        ("EX 1 -2 1 0 30 0 0 30 0 0\nXQ 0", "THE NUMBER OF INCIDENCE ANGLES IS NEGATIVE"),
    ];
    for (cards, message) in cards {
        let (error, _) = stopped(&format!("{MONOPOLE}{cards}\nEN\n"));
        assert!(error.contains(message), "{cards}: {error}");
    }
}

/// The PT options above zero print the currents as a receiving pattern.
/// Without an incident wave nec2c reads beyond its table of names for the
/// heading and, for the normalized pattern, writes to a buffer it has not got.
#[test]
fn a_receiving_pattern_without_an_incident_wave_is_an_error() {
    for option in [1, 2, 3] {
        let (error, _) = stopped(&format!("{MONOPOLE}PT {option} 1 5 5\nEX 0 1 5 0 1 0\nXQ 0\nEN\n"));
        assert!(error.contains("IT NEEDS AN INCIDENT WAVE"), "{option}: {error}");
        let wave = format!("{MONOPOLE}PT {option} 1 5 5\nEX 1 2 1 0 30 0 0 30 0 0\nXQ 0\nEN\n");
        assert_eq!(nec2_deck::run(&wave).error, None, "{option}");
    }
}

/// The last line of defence: whatever the cause, a result that is not a
/// number ends the run. Here it is a ground wave over a ground whose
/// permittivity is beyond what the formulas can hold.
#[test]
fn a_result_that_is_not_a_number_ends_the_run() {
    let deck = format!(
        "{MONOPOLE}GN 0 0 0 0 -1e30 0.005\nEX 0 1 5 0 1 0\nFR 0 1 0 0 299.8 0\nRP 1 2 1 1000 0 0 0.5 0 100\nEN\n"
    );
    let (error, report) = stopped(&deck);
    assert!(error.contains("VALUES THAT ARE NOT NUMBERS"), "{error}");
    assert!(report.contains("RADIATED FIELDS NEAR GROUND"));
}

/// nec2c reads beyond its table of names to print such a type, crashes if
/// the type is large and else takes the card for a current source.
#[test]
fn an_excitation_type_that_nec2_does_not_have_is_an_error() {
    for kind in [6, 7, 100_000, -1] {
        let (error, _) = stopped(&format!("{MONOPOLE}EX {kind} 1 5 0 1 0\nXQ 0\nEN\n"));
        assert!(error.contains(&format!("EXCITATION TYPE {kind} IS NOT DEFINED")), "{error}");
    }
}
