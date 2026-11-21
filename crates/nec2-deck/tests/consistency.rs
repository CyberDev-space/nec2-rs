//! Checks that need no reference program: a result must not depend on what
//! was computed before it.

/// The line of the current table for segment 6 of tag 1, last occurrence.
fn current_on_segment_6(report: &str) -> String {
    report
        .lines()
        .rfind(|line| {
            let mut tokens = line.split_whitespace();
            tokens.next() == Some("6") && tokens.next() == Some("1") && tokens.count() == 8
        })
        .expect("current of segment 6")
        .to_string()
}

fn plane_wave_deck(angles: u32, first_theta: u32, line: bool) -> String {
    let line = if line { "TL 1 6 2 6 50 0 0 0 0 0\n" } else { "" };
    format!(
        "CM Plane wave on two dipoles joined by a transmission line\n\
         CE\n\
         GW 1 11 0 0 -0.25 0 0 0.25 0.001\n\
         GW 2 11 0.3 0 -0.25 0.3 0 0.25 0.001\n\
         GE 0\n\
         FR 0 1 0 0 299.8 0\n\
         {line}\
         EX 1 {angles} 1 0 {first_theta} 0 0 30 0 0\n\
         PT 0 1 6 6\n\
         XQ 0\n\
         EN\n"
    )
}

/// The network solution is kept from one incidence angle to the next, as in
/// NEC-2. (nec2c 1.3.1 drops it: from the second angle on it gives the
/// currents of the structure without the line.)
#[test]
fn networks_act_at_every_incidence_angle() {
    let sweep = nec2_deck::run(&plane_wave_deck(2, 60, true));
    let single = nec2_deck::run(&plane_wave_deck(1, 90, true));
    let without_line = nec2_deck::run(&plane_wave_deck(1, 90, false));
    assert_eq!(sweep.error, None);
    assert_eq!(single.error, None);
    assert_eq!(without_line.error, None);

    let second_angle = current_on_segment_6(&sweep.report);
    assert_eq!(second_angle, current_on_segment_6(&single.report));
    assert_ne!(second_angle, current_on_segment_6(&without_line.report));
}

/// Every request closes the group of EX cards before it, also an XQ card
/// that asks for a pattern. (NEC-2 leaves the group open there: the next EX
/// card adds to it, and the request after it is served with the old
/// currents, or with a loop over incidence angles that has already ended.)
#[test]
fn a_request_closes_the_group_of_excitation_cards() {
    let deck = |cards: &str| {
        format!(
            "CM Two dipoles\nCE\nGW 1 9 0 0 -0.25 0 0 0.25 0.001\nGW 2 9 0.2 0 -0.2 0.2 0 0.2 0.001\nGE 0\n\
             FR 0 1 0 0 299.8 0\n{cards}EN\n"
        )
    };
    // The report from the last table with this title on, without the echo of the cards.
    let last_table = |deck: &str, title: &str| {
        let output = nec2_deck::run(deck);
        assert_eq!(output.error, None);
        let report = &output.report[output.report.rfind(title).expect("the table")..];
        report.lines().filter(|line| !line.contains("DATA CARD No:")).collect::<Vec<_>>().join("\n")
    };

    // The second source alone, not the pattern of the first one again.
    let second = "EX 0 2 5 0 1 0\nXQ 1\n";
    let both = deck(&format!("EX 0 1 5 0 1 0\nXQ 1\n{second}"));
    let title = "ANTENNA INPUT PARAMETERS";
    assert_eq!(last_table(&both, title), last_table(&deck(second), title));

    // Both incidence angles again, not a third one.
    let wave = "EX 1 2 1 0 30 0 0 30 0 0\nPT 3 1 5 5\nXQ 0\n";
    let twice = deck(&format!("EX 1 2 1 0 30 0 0 30 0 0\nXQ 1\n{wave}"));
    let title = "NORMALIZED RECEIVING PATTERN";
    assert_eq!(last_table(&twice, title), last_table(&deck(wave), title));
}

/// Each frequency of a sweep over Sommerfeld ground gives what a run at
/// that frequency alone gives: nothing of the previous interpolation tables
/// is reused.
#[test]
fn sommerfeld_sweep_equals_single_frequencies() {
    let deck = |frequency_card: &str| {
        format!(
            "CM Dipole over average ground\n\
             CE\n\
             GW 1 15 0 -0.25 0.1 0 0.25 0.1 0.001\n\
             GE 0\n\
             GN 2 0 0 0 13 0.005\n\
             EX 0 1 8 0 1 0\n\
             {frequency_card}\n\
             PT -1\n\
             XQ 0\n\
             EN\n"
        )
    };
    let impedances = |report: &str| -> Vec<String> {
        let lines: Vec<&str> = report.lines().collect();
        lines
            .iter()
            .enumerate()
            .filter(|(_, line)| line.contains("ANTENNA INPUT PARAMETERS"))
            .map(|(i, _)| lines[i + 3].to_string())
            .collect()
    };

    let sweep = nec2_deck::run(&deck("FR 0 3 0 0 280 20"));
    assert_eq!(sweep.error, None);
    let sweep = impedances(&sweep.report);
    assert_eq!(sweep.len(), 3);
    for (k, frequency) in [280, 300, 320].into_iter().enumerate() {
        let single = nec2_deck::run(&deck(&format!("FR 0 1 0 0 {frequency} 0")));
        assert_eq!(impedances(&single.report), vec![sweep[k].clone()], "{frequency} MHz");
    }
}

/// Input impedance of the first source of a report.
fn impedance(deck: &str) -> (f64, f64) {
    let output = nec2_deck::run(deck);
    assert_eq!(output.error, None);
    let lines: Vec<&str> = output.report.lines().collect();
    let title = lines.iter().position(|line| line.contains("ANTENNA INPUT PARAMETERS")).expect("input parameters");
    let row: Vec<f64> = lines[title + 3].split_whitespace().map(|token| token.parse().expect("number")).collect();
    (row[6], row[7])
}

#[track_caller]
fn assert_same_impedance(a: (f64, f64), b: (f64, f64), relative: f64) {
    let size = a.0.hypot(a.1);
    assert!((a.0 - b.0).hypot(a.1 - b.1) <= relative * size, "{a:?} and {b:?}");
}

/// A wire on a patch gives the same result wherever its card stands in the
/// deck. (nec2c 1.3.1 takes the coefficient of the wire's own basis
/// function from the wire that was entered last.)
#[test]
fn a_wire_on_a_patch_does_not_depend_on_the_order_of_the_cards() {
    let plate = "SM 3 3 -0.3 -0.3 0 0.3 -0.3 0\nSC 0 0 0.3 0.3 0\n";
    let stub = "GW 1 1 0 0 0 0 0 0.05 0.002\n";
    let other = "GW 2 3 0.6 0 0.1 0.6 0 0.4 0.002\n";
    let run = "GE 0\nEX 0 1 1 0 1 0\nFR 0 1 0 0 299.8 0\nXQ 0\nEN\n";
    let first = impedance(&format!("CM Stub first\nCE\n{plate}{stub}{other}{run}"));
    let last = impedance(&format!("CM Stub last\nCE\n{plate}{other}{stub}{run}"));
    assert_same_impedance(first, last, 1.0e-9);
    // The reactance of the short stub, not a third of it.
    assert!((first.1 + 1091.0).abs() < 0.5, "{first:?}");
}

/// A structure with patches that is completed by rotation is solved through
/// its symmetry. The result is that of the same structure entered sector by
/// sector, which is solved as a whole.
#[test]
fn symmetry_with_patches_equals_the_structure_entered_in_full() {
    let run = "GE 0\nEX 0 1 1 0 1 0\nFR 0 1 0 0 299.8 0\nXQ 0\nEN\n";
    let rotated = format!(
        "CM Four sectors by rotation\nCE\nSP 0 0 0.3 0 0 90 0 0.06\nGW 1 3 0.3 0 0 0.3 0 0.2 0.002\nGR 0 4\n{run}"
    );
    let mut full = String::from("CM Four sectors entered one by one\nCE\n");
    let places = [(0.3, 0.0), (0.0, 0.3), (-0.3, 0.0), (0.0, -0.3)];
    for (x, y) in places {
        full += &format!("SP 0 0 {x} {y} 0 90 0 0.06\n");
    }
    for (x, y) in places {
        full += &format!("GW 1 3 {x} {y} 0 {x} {y} 0.2 0.002\n");
    }
    full += run;
    let symmetric = nec2_deck::run(&rotated);
    assert!(symmetric.report.contains("STRUCTURE HAS 4 FOLD ROTATIONAL SYMMETRY"));
    assert!(!nec2_deck::run(&full).report.contains("SYMMETRY\n"));
    assert_same_impedance(impedance(&rotated), impedance(&full), 1.0e-6);
}
