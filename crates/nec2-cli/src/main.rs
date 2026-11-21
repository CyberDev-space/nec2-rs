//! Command-line NEC-2 solver.
//!
//! ```text
//! nec2 -i deck.nec [-o report.out] [--overlap error|warn|off]
//! ```
//!
//! Without `-o` the report goes next to the deck, with the extension `.out`,
//! as nec2c does. `-o -` writes it to the standard output. A deck with a PL
//! card also gets a plot file, the name of the deck with `.plt` added, as in
//! nec2c. `--overlap` says
//! what happens to wires that occupy the same space: they are listed as
//! warnings (the default), stop the run, or are not looked for.

use std::path::PathBuf;
use std::process::ExitCode;

const USAGE: &str = "usage: nec2 -i <input-file> [-o <output-file>] [--overlap error|warn|off]\n       nec2 -h | -v";

fn main() -> ExitCode {
    let mut input: Option<PathBuf> = None;
    let mut output: Option<PathBuf> = None;
    let mut options = nec2_deck::Options::default();
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-i" => input = args.next().map(PathBuf::from),
            "-o" => output = args.next().map(PathBuf::from),
            "--overlap" => {
                options.overlap = match args.next().as_deref() {
                    Some("error") => nec2_deck::OverlapCheck::Error,
                    Some("warn") => nec2_deck::OverlapCheck::Warn,
                    Some("off") => nec2_deck::OverlapCheck::Off,
                    _ => {
                        eprintln!("{USAGE}");
                        return ExitCode::FAILURE;
                    }
                };
            }
            "-h" => {
                println!("{USAGE}");
                return ExitCode::SUCCESS;
            }
            "-v" => {
                println!("nec2 {}", env!("CARGO_PKG_VERSION"));
                return ExitCode::SUCCESS;
            }
            _ => {
                eprintln!("{USAGE}");
                return ExitCode::FAILURE;
            }
        }
    }
    let Some(input) = input else {
        eprintln!("{USAGE}");
        return ExitCode::FAILURE;
    };

    let deck = match std::fs::read_to_string(&input) {
        Ok(deck) => deck,
        Err(error) => {
            eprintln!("nec2: {}: {error}", input.display());
            return ExitCode::FAILURE;
        }
    };
    let result = nec2_deck::run_with(&deck, &options);

    let output = output.unwrap_or_else(|| input.with_extension("out"));
    let written = if output.as_os_str() == "-" {
        println!("{}", result.report);
        Ok(())
    } else {
        std::fs::write(&output, &result.report)
    };
    if let Err(error) = written {
        eprintln!("nec2: {}: {error}", output.display());
        return ExitCode::FAILURE;
    }
    if let Some(plot) = &result.plot {
        let mut name = input.clone().into_os_string();
        name.push(".plt");
        if let Err(error) = std::fs::write(&name, plot) {
            eprintln!("nec2: {}: {error}", PathBuf::from(name).display());
            return ExitCode::FAILURE;
        }
    }

    match result.error {
        Some(message) => {
            eprintln!("nec2: {message}");
            ExitCode::FAILURE
        }
        None => ExitCode::SUCCESS,
    }
}
