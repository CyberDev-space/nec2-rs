//! Robustness check: runs the test decks after damaging them at random and
//! looks for panics, for runs that do not end and for reports with values
//! that are not numbers. A development aid, not part of the test suite:
//!
//! ```sh
//! cargo test --release -p nec2-deck --test fuzz -- --ignored --nocapture
//! ```
//!
//! `NEC2_FUZZ_RUNS` is the number of decks (3000 by default), `NEC2_FUZZ_SEED`
//! the seed of the random numbers (1) and `NEC2_FUZZ_OUT` the directory for
//! the decks that fail (a `nec2-fuzz` directory among the temporary files).
//! A damaged count can ask for a very large structure, so the decks run with
//! a limit on the number of equations.
//!
//! With `NEC2_FUZZ_REFERENCE` set to the path of a nec2c program, every deck
//! is also given to it, which checks the compatibility rule: where nec2c ends
//! with a report of finite numbers, this program prints the same numbers, or
//! it stops with one of the errors that the rule allows. The decks that do
//! neither are written to the directory, each with the reason in a file of
//! the same name that ends in `.txt`: `stricter_*` (another error here),
//! `differs_*` (other numbers) and `slight_*` (numbers that differ by less
//! than a hundred units of the last digit). They are for a person to look at:
//! a deck damaged at random can be so ill-conditioned that two correct
//! programs disagree. Besides such noise (a structure a millionth of a
//! wavelength long, currents of 1e-19), the decks that differ are those of
//! the differences that `docs/nec2c-deviations.md` describes, mostly these:
//! pattern angles beyond the horizon over a ground, a further request after
//! a loop over incidence angles, networks in such a loop, and print options
//! for which nec2c prints memory that it has not set.
//!
//! `NEC2_FUZZ_REPLAY` names a directory of decks to run instead of damaged
//! ones, for example the decks of an earlier run after a correction.

mod support;

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

/// Seconds after which a run counts as one that does not end.
const TIME_LIMIT: u64 = 60;

/// Largest structure that a damaged deck may ask for.
const MAX_EQUATIONS: usize = 600;

/// Lines of the reference report that are known to differ, see
/// `docs/nec2c-deviations.md`: the average gain has other weights at the poles.
const KNOWN_DIFFERENCES: [&str; 1] = ["AVERAGE POWER GAIN"];

/// Errors of this program for decks that nec2c runs to the end, with numbers
/// that mean nothing. `docs/nec2c-deviations.md` has every one of them.
const ALLOWED_ERRORS: [&str; 14] = [
    "IS NEGATIVE",
    "SEGMENT DATA ERROR",
    "EXCITATION TYPE",
    "IMPROPER LOAD TYPE",
    "PATCH DATA ERROR",
    "A FINITE GROUND NEEDS",
    "THE FREQUENCY",
    "CLOSER THAN THE CONNECTION TOLERANCE",
    "THE PATTERN MODE",
    "NO EX CARD DEFINES A SOURCE",
    "THE PRINT OPTION OF THE PT CARD",
    "THE GROUND WAVE IS ASKED FOR AT THE ORIGIN",
    "TOO LARGE FOR THE MACHINE",
    "ARE NOT NUMBERS",
];

/// Seconds that the reference program gets for one deck.
const REFERENCE_TIME_LIMIT: u64 = 20;

/// Largest report of the reference program, in bytes. An integration that
/// does not converge can print a message at every step, hundreds of
/// megabytes in all.
const REFERENCE_REPORT_LIMIT: u64 = 20 << 20;

/// xorshift64*: small, fast and the same on every machine.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }

    fn pick<'a, T>(&mut self, items: &'a [T]) -> &'a T {
        &items[self.below(items.len())]
    }
}

/// Values that take the place of a number on a card. The counts stay small:
/// the check is for faults of the program, not for the size of a structure.
const NUMBERS: [&str; 22] = [
    "0", "1", "2", "3", "7", "30", "-1", "-2", "0.5", "-0.5", "1e-3", "1e3", "1e-30", "1e30", "-1e30", "1e400", "90",
    "180", "360", "0.001", "299.8", "1.5",
];

/// Cards that are put into a deck at a random place.
const CARDS: [&str; 26] = [
    "GE 0",
    "GE 1",
    "GE -1",
    "GN 1",
    "GN 0 0 0 0 13 0.005",
    "GN 2 0 0 0 13 0.005",
    "GN -1",
    "GX 0 110",
    "GR 0 3",
    "GM 0 1 0 0 30 0 0 0.3",
    "GS 0 0 0.5",
    "SP 0 0 0.1 0.1 0.1 90 0 0.01",
    "SM 2 2 0 0 0 0.2 0 0",
    "SC 0 0 0.2 0.2 0",
    "EX 0 1 1 0 1 0",
    "EX 1 2 2 0 0 0 0 30 30 0",
    "EX 4 0 0 0 0.1 0.1 0.1 0 0 1",
    "EX 5 1 1 0 1 0",
    "LD 5 0 0 0 5.8E7",
    "TL 1 1 1 2 50 0",
    "NT 1 1 1 2 0 0.01 0 0.02 0 0.01",
    "RP 0 4 2 1502 0 0 30 90",
    "NE 0 2 1 2 0.1 0.1 0.1 0.1 0 0.1",
    "PT 2 1 1 1",
    "XQ 0",
    "FR 0 2 0 0 150 150",
];

fn is_number(token: &str) -> bool {
    token.parse::<f64>().is_ok()
}

/// Damages a deck in one to three places.
fn mutate(deck: &str, rng: &mut Rng) -> String {
    let mut lines: Vec<String> = deck.lines().map(str::to_string).collect();
    for _ in 0..1 + rng.below(3) {
        if lines.is_empty() {
            break;
        }
        let at = rng.below(lines.len());
        match rng.below(10) {
            // A number of a card gets another value.
            0..=5 => {
                let mut tokens: Vec<String> = lines[at].split_whitespace().map(str::to_string).collect();
                let numbers: Vec<usize> = (1..tokens.len()).filter(|&k| is_number(&tokens[k])).collect();
                if !numbers.is_empty() {
                    let k = *rng.pick(&numbers);
                    tokens[k] = rng.pick(&NUMBERS).to_string();
                    lines[at] = tokens.join(" ");
                }
            }
            6 => {
                lines.remove(at);
            }
            7 => {
                let copy = lines[at].clone();
                lines.insert(rng.below(lines.len() + 1), copy);
            }
            8 => {
                let other = rng.below(lines.len());
                lines.swap(at, other);
            }
            _ => lines.insert(at, rng.pick(&CARDS).to_string()),
        }
    }
    lines.join("\n") + "\n"
}

/// A value that is not a number, in any of the spellings of C and Rust.
fn not_finite(token: &str) -> bool {
    let token = token.to_ascii_uppercase();
    let bare = token.trim_start_matches(['-', '+']).trim_end_matches('J');
    bare == "NAN" || bare == "INF"
}

/// A deck with a number so large or so small that its results are rounding
/// errors, in which two correct programs differ.
fn extreme(deck: &str) -> bool {
    let numbers = deck.lines().filter(|line| !line.starts_with('C')).flat_map(str::split_whitespace);
    numbers.filter_map(|token| token.parse::<f64>().ok()).any(|x| x != 0.0 && !(1.0e-20..=1.0e20).contains(&x.abs()))
}

/// Runs the reference program on a deck. Returns its report, or `None` if it
/// refused the deck, crashed, did not end in time or printed without end.
fn run_reference(program: &Path, work: &Path, deck: &str) -> Option<String> {
    let report = work.join("a.out");
    let _ = std::fs::remove_file(&report);
    std::fs::write(work.join("a.nec"), deck).expect("deck for the reference program");
    let mut child = Command::new(program)
        .args(["-i", "a.nec", "-o", "a.out"])
        .current_dir(work)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .expect("reference program");
    let started = Instant::now();
    let too_large = || std::fs::metadata(&report).is_ok_and(|file| file.len() > REFERENCE_REPORT_LIMIT);
    let status = loop {
        match child.try_wait().expect("reference program") {
            Some(status) => break status,
            None if started.elapsed() > Duration::from_secs(REFERENCE_TIME_LIMIT) || too_large() => {
                let _ = child.kill();
                let _ = child.wait();
                let _ = std::fs::remove_file(&report);
                return None;
            }
            None => std::thread::sleep(Duration::from_millis(1)),
        }
    };
    if !status.success() || too_large() {
        let _ = std::fs::remove_file(&report);
        return None;
    }
    Some(String::from_utf8_lossy(&std::fs::read(report).ok()?).into_owned())
}

/// The decks of a directory with their names, in the order of the names.
fn read_decks(directory: &Path) -> Vec<(String, String)> {
    let mut decks: Vec<(String, String)> = std::fs::read_dir(directory)
        .expect("deck directory")
        .filter_map(|entry| {
            let path = entry.ok()?.path();
            let name = path.file_name()?.to_str()?.strip_suffix(".nec")?.to_string();
            Some((name, std::fs::read_to_string(&path).ok()?))
        })
        .collect();
    decks.sort();
    decks
}

#[test]
#[ignore = "development aid"]
fn damaged_decks_do_not_panic_or_hang() {
    let variable = |name: &str| std::env::var(name).ok();
    let runs: u64 = variable("NEC2_FUZZ_RUNS").and_then(|v| v.parse().ok()).unwrap_or(3000);
    let seed: u64 = variable("NEC2_FUZZ_SEED").and_then(|v| v.parse().ok()).unwrap_or(1);
    let out = variable("NEC2_FUZZ_OUT").map_or_else(|| std::env::temp_dir().join("nec2-fuzz"), PathBuf::from);
    std::fs::create_dir_all(&out).expect("directory for failing decks");
    let reference = variable("NEC2_FUZZ_REFERENCE").map(PathBuf::from);
    let work = out.join("work");
    std::fs::create_dir_all(&work).expect("directory for the reference program");

    let decks = read_decks(&PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/decks"));
    let replayed = variable("NEC2_FUZZ_REPLAY").map(|directory| read_decks(Path::new(&directory)));
    let runs = replayed.as_ref().map_or(runs, |decks| decks.len() as u64);

    // The deck that is running, for the watchdog.
    let current = Arc::new(Mutex::new((Instant::now(), String::new())));
    {
        let current = Arc::clone(&current);
        let out = out.clone();
        std::thread::spawn(move || {
            loop {
                std::thread::sleep(Duration::from_secs(1));
                let (started, deck) = &*current.lock().unwrap();
                if !deck.is_empty() && started.elapsed() > Duration::from_secs(TIME_LIMIT) {
                    let path = out.join(format!("hang_seed{seed}.nec"));
                    let _ = std::fs::write(&path, deck);
                    eprintln!("\nA RUN DOES NOT END AFTER {TIME_LIMIT} s: {}", path.display());
                    // Not `abort`: that would leave a core dump of the whole process.
                    std::process::exit(3);
                }
            }
        });
    }

    let options = nec2_deck::Options { max_equations: Some(MAX_EQUATIONS), ..Default::default() };
    let mut rng = Rng(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1);
    let (mut stopped, mut panics, mut nans) = (0u64, Vec::new(), Vec::new());
    // Against the reference: the decks it solved and what became of them here.
    let (mut solved, mut compared, mut allowed, mut ill_conditioned) = (0u64, 0u64, 0u64, 0u64);
    let (mut differs, mut slight, mut stricter) = (Vec::new(), Vec::new(), Vec::new());
    for run in 0..runs {
        let (name, damaged) = match &replayed {
            Some(decks) => decks[run as usize].clone(),
            None => {
                let (name, deck) = rng.pick(&decks);
                (name.clone(), mutate(deck, &mut rng))
            }
        };
        let name = name.as_str();
        *current.lock().unwrap() = (Instant::now(), damaged.clone());
        let result = std::panic::catch_unwind(|| nec2_deck::run_with(&damaged, &options));
        *current.lock().unwrap() = (Instant::now(), String::new());

        let save = |kind: &str, reason: &str| {
            let path = out.join(format!("{kind}_seed{seed}_{run}_{name}.nec"));
            let _ = std::fs::write(&path, &damaged);
            if !reason.is_empty() {
                let _ = std::fs::write(path.with_extension("txt"), reason);
            }
            path
        };
        let Ok(output) = result else {
            panics.push(save("panic", ""));
            continue;
        };
        if output.error.is_some() {
            stopped += 1;
        } else if output.report.split_whitespace().any(not_finite) {
            nans.push(save("nan", ""));
        }

        // The limit on the size is a matter of this check, not of the rule.
        let over_the_limit = output.error.as_deref().is_some_and(|error| error.contains("MORE THAN THE LIMIT"));
        let Some(program) = reference.as_deref().filter(|_| !over_the_limit) else { continue };
        let Some(expected) = run_reference(program, &work, &damaged) else { continue };
        if expected.split_whitespace().any(not_finite) {
            continue;
        }
        solved += 1;
        if let Some(error) = &output.error {
            if ALLOWED_ERRORS.iter().any(|allowed| error.contains(allowed)) {
                allowed += 1;
            } else {
                stricter.push(save("stricter", error));
            }
        } else if output.report.contains("ILL-CONDITIONED") || extreme(&damaged) {
            ill_conditioned += 1;
        } else {
            compared += 1;
            let actual = support::without_warnings(&output.report);
            let compare = |units: Option<f64>| {
                let mut mismatches = support::compare(&expected, &actual, units, false);
                mismatches.retain(|mismatch| !KNOWN_DIFFERENCES.iter().any(|known| mismatch.reference.contains(known)));
                mismatches
            };
            let mismatches = compare(None);
            if !mismatches.is_empty() {
                let reason = support::describe(name, &mismatches);
                if compare(Some(100.0)).is_empty() {
                    slight.push(save("slight", &reason));
                } else {
                    differs.push(save("differs", &reason));
                }
            }
        }
    }

    println!(
        "{runs} decks: {stopped} stopped with an error, {} printed values that are not numbers, {} panics",
        nans.len(),
        panics.len()
    );
    if reference.is_some() {
        println!(
            concat!(
                "the reference solved {} of them with finite numbers: {} compared ({} with other numbers, {} ",
                "slightly), {} stopped here with an error that the rule allows and {} with another, {} left out as ",
                "ill-conditioned or for a number beyond 1e20"
            ),
            solved,
            compared,
            differs.len(),
            slight.len(),
            allowed,
            stricter.len(),
            ill_conditioned
        );
    }
    for path in stricter.iter().chain(&differs).take(20) {
        println!("  to look at: {}", path.display());
    }
    for path in nans.iter().chain(&panics) {
        println!("  FAILED: {}", path.display());
    }
    assert!(panics.is_empty(), "{} decks made the program panic", panics.len());
    assert!(nans.is_empty(), "{} reports hold values that are not numbers", nans.len());
}
