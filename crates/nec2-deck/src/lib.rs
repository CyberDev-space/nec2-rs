//! Reader and interpreter for NEC-2 input decks.
//!
//! A deck is the classic input of NEC-2: comment cards, the geometry section
//! and the program cards that set parameters and request computations.
//! [`run`] interprets a deck with the `nec2` engine and returns the report in
//! the layout of the NEC-2 output file, as nec2c prints it.
//!
//! ```
//! let deck = "\
//! CM Half-wave dipole in free space
//! CE
//! GW 1 21 0 0 -0.25 0 0 0.25 0.001
//! GE 0
//! EX 0 1 11 0 1 0
//! FR 0 1 0 0 299.8 0
//! XQ 0
//! EN
//! ";
//! let output = nec2_deck::run(deck);
//! assert!(output.error.is_none());
//! assert!(output.report.contains("ANTENNA INPUT PARAMETERS"));
//! ```

mod cfmt;
mod parse;
mod run;

pub use nec2::OverlapCheck;
pub use run::{Options, Output, run, run_with};
