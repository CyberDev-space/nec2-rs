//! The results as numbers (`Output::jobs`) must be the numbers of the report.

mod support;

use std::path::PathBuf;

use nec2_deck::{Complex64, Output};
use support::{num, tables, values};

fn run(name: &str) -> Output {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(format!("../../tests/fixtures/decks/{name}.nec"));
    let deck = std::fs::read_to_string(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()));
    let output = nec2_deck::run(&deck);
    assert_eq!(output.error, None, "{name}");
    output
}

/// A number of the report has five significant digits.
#[track_caller]
fn assert_printed(value: f64, printed: &str) {
    let printed = num(printed);
    assert!((value - printed).abs() <= 1.0e-4 * printed.abs() + 1.0e-30, "{value} is printed as {printed}");
}

/// A complex number whose parts are printed in two columns: a part far
/// smaller than the magnitude has fewer digits that mean something.
#[track_caller]
fn assert_printed_complex(value: Complex64, re: &str, im: &str) {
    let printed = Complex64::new(num(re), num(im));
    assert!((value - printed).norm() <= 1.0e-4 * printed.norm() + 1.0e-30, "{value} is printed as {printed}");
}

#[test]
fn a_sweep_has_a_solution_for_every_frequency() {
    let output = run("s1_dipole_sweep");
    let report = &output.report;
    let [job] = &output.jobs[..] else { panic!("one job") };
    let frequencies: Vec<f64> = job.solutions.iter().map(|solution| solution.frequency).collect();
    assert_eq!(frequencies, [280.0, 290.0, 300.0, 310.0, 320.0]);

    let sources = tables(report, "ANTENNA INPUT PARAMETERS");
    let currents = tables(report, "CURRENTS AND LOCATION");
    let patterns = tables(report, "RADIATION PATTERNS");
    for (k, solution) in job.solutions.iter().enumerate() {
        // TAG SEG V V I I Z Z Y Y POWER
        let [source] = &solution.sources[..] else { panic!("one source") };
        let row = &sources[k][0];
        assert_eq!((source.tag, source.segment), (1, 11));
        assert_printed_complex(source.voltage, row[2], row[3]);
        assert_printed_complex(source.current, row[4], row[5]);
        assert_printed_complex(source.impedance, row[6], row[7]);
        assert_printed_complex(source.admittance, row[8], row[9]);
        assert_printed(source.power, row[10]);
        assert_eq!(solution.incident_wave, None);
        assert!(solution.network_ports.is_empty());

        // SEG TAG X Y Z LENGTH REAL IMAGINARY MAGNITUDE PHASE
        assert_eq!(solution.currents.len(), 21);
        for (current, row) in solution.currents.iter().zip(&currents[k]) {
            assert_printed_complex(*current, row[6], row[7]);
        }

        let power = solution.power.expect("a power budget");
        assert_printed(power.input, &format!("{:E}", values(report, "INPUT POWER")[k]));
        assert_printed(power.radiated, &format!("{:E}", values(report, "RADIATED POWER")[k]));
        assert_eq!((power.structure_loss, power.network_loss), (0.0, 0.0));

        // THETA PHI VERTC HORIZ TOTAL ...
        let [pattern] = &solution.patterns[..] else { panic!("one pattern") };
        assert_eq!(pattern.rows.len(), 7);
        for (row, printed) in pattern.rows.iter().zip(&patterns[k]) {
            assert_eq!((row.theta, row.phi), (num(printed[0]), num(printed[1])));
            if row.gain_total > -900.0 {
                assert!((row.gain_total - num(printed[4])).abs() < 0.00501, "{} dB", row.gain_total);
            }
        }
    }
}

#[test]
fn the_segments_of_the_structure() {
    let output = run("s1_dipole");
    let segments = &output.jobs[0].segments;
    assert_eq!(segments.len(), 21);
    assert!(segments.iter().all(|segment| segment.tag == 1 && segment.radius == 0.001));
    assert_eq!((segments[0].end1, segments[20].end2), ([0.0, 0.0, -0.25], [0.0, 0.0, 0.25]));
    for pair in segments.windows(2) {
        assert!((pair[0].end2[2] - pair[1].end1[2]).abs() < 1.0e-12);
    }
    // 19 angles of theta for two angles of phi.
    assert_eq!(output.jobs[0].solutions[0].patterns[0].rows.len(), 38);
    assert!(output.warnings.is_empty());
}

#[test]
fn an_incident_wave_has_a_solution_for_every_direction() {
    let output = run("s5_ex_plane_wave");
    let solutions = &output.jobs[0].solutions;
    let directions: Vec<[f64; 3]> = solutions.iter().map(|s| s.incident_wave.expect("a plane wave")).collect();
    let expected: Vec<[f64; 3]> = [30.0, 60.0, 90.0]
        .into_iter()
        .chain((0..7).map(|k| 15.0 * f64::from(k)))
        .map(|theta| [theta, 0.0, 0.0])
        .collect();
    assert_eq!(directions, expected);
    // The first three have the pattern of the scattered field, and no wave has a power budget.
    let patterns: Vec<usize> = solutions.iter().map(|s| s.patterns.len()).collect();
    assert_eq!(patterns, [1, 1, 1, 0, 0, 0, 0, 0, 0, 0]);
    assert!(solutions.iter().all(|s| s.power.is_none() && s.sources.is_empty() && s.currents.len() == 21));
}

#[test]
fn every_job_has_its_results() {
    let output = run("s5_two_jobs");
    let segments: Vec<usize> = output.jobs.iter().map(|job| job.segments.len()).collect();
    assert_eq!(segments, [11, 6]);
    let impedances = tables(&output.report, "ANTENNA INPUT PARAMETERS");
    for (job, table) in output.jobs.iter().zip(&impedances) {
        let [solution] = &job.solutions[..] else { panic!("one solution") };
        assert_printed_complex(solution.sources[0].impedance, table[0][6], table[0][7]);
        assert_eq!(solution.patterns.len(), 1);
    }
}

#[test]
fn warnings_and_what_was_computed_before_an_error() {
    let output = run("s7_ge1_free_space_touch");
    let [warning] = &output.warnings[..] else { panic!("one warning: {:?}", output.warnings) };
    assert!(warning.starts_with("SEGMENT ENDS JOINED TO THEIR IMAGE BY THE GE CARD: 1,"), "{warning}");

    // The second request is a pattern mode that NEC-2 does not have.
    let deck = "CM Dipole\nCE\nGW 1 9 0 0 -0.25 0 0 0.25 0.001\nGE 0\nEX 0 1 5 0 1 0\nFR 0 1 0 0 299.8 0\nXQ 0\n\
                RP 7 4 1 1000 0 0 30 0\nEN\n";
    let output = nec2_deck::run(deck);
    assert!(output.error.is_some());
    let [solution] = &output.jobs[0].solutions[..] else { panic!("one solution") };
    assert_eq!((solution.sources.len(), solution.patterns.len()), (1, 0));

    // An error in the geometry leaves no job.
    let output = nec2_deck::run("CM No structure\nCE\nGW 1 -3 0 0 0 0 0 1 0.001\nGE 0\nEN\n");
    assert!(output.error.is_some() && output.jobs.is_empty());
}
