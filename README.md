# nec2-rs

[![crates.io](https://img.shields.io/crates/v/nec2.svg)](https://crates.io/crates/nec2)
[![docs.rs](https://img.shields.io/docsrs/nec2)](https://docs.rs/nec2)
[![license](https://img.shields.io/crates/l/nec2.svg)](LICENSE)

A Rust implementation of NEC-2, the Numerical Electromagnetics Code for wire antennas (method of moments).

It is a port of nec2c 1.3.1, the public-domain C translation of NEC-2, checked report by report against nec2c. The license is MIT.

## Status

Working and tested. The first release, 0.1.0, is on crates.io.

| Feature | State |
|---|---|
| Wires: GW (with GC taper), GA, GH, GM, GX, GR, GS, GE | done |
| Free space, perfect ground, finite ground by reflection coefficients | done |
| Finite ground by Sommerfeld-Norton, tables computed in the program | done |
| Radial-wire ground screen, second medium (cliff), ground wave | done |
| Sources: voltage (EX 0), current slope (EX 5), plane waves (EX 1 to 3), current element (EX 4) | done |
| Loads LD 0 to 5, transmission lines TL, networks NT | done |
| Far fields RP, near fields NE and NH, coupling CP, charges PQ, print control PT | done |
| Extended thin-wire kernel EK, interaction range KH, frequency sweeps | done |
| Structure symmetry in the solution (planes and rotation) | done |
| Surface patches SP, SM, SC, and wires that end on them | done |
| Plot files (PL) | done |
| Numerical Green's function GF, WG | not implemented (removed from nec2c too) |

The test suite runs 143 decks and three decks that must stop with an error, and compares every number of the reports with nec2c: 127 decks with the stock nec2c 1.3.1 and 16 with nec2c corrected in the places described in [docs/nec2c-deviations.md](docs/nec2c-deviations.md). Numbers must be printed in the same form and may differ by one unit of the last printed digit.

## Compatibility

nec2-rs is meant to replace NEC-2 engines without surprises. A deck that NEC-2 solves gives the same numbers, in the report layout of nec2c. Where NEC-2 is wrong the result is better, and where a model deserves no trust the report gains a line that begins with `WARNING:`. Where NEC-2 fails, or goes on to print values that are not numbers or numbers that mean nothing, the program stops with an error that names the cause. A run that ends without an error has only numbers in its report. Thirty-four of the test decks are there to hold this rule, and a test that damages decks at random compares nec2-rs with nec2c on the result. The details are in [docs/nec2c-deviations.md](docs/nec2c-deviations.md).

## Crates

| Crate | Contents |
|---|---|
| [`nec2`](https://crates.io/crates/nec2) | the engine: wires and patches, matrix fill and solution, fields. Computation only: no files, no threads. |
| [`nec2-deck`](https://crates.io/crates/nec2-deck) | reads a deck of NEC cards, runs it and returns the report in the nec2c layout, and the plot data |
| [`nec2-cli`](https://crates.io/crates/nec2-cli) | the command `nec2 -i deck.nec [-o report.out] [--overlap error\|warn\|off]` |

## Use

The command:

```sh
cargo install nec2-cli
nec2 -i deck.nec -o report.out
```

Or from this repository:

```sh
cargo build --release
./target/release/nec2 -i tests/fixtures/decks/s1_dipole.nec -o dipole.out
```

From Rust, a deck with `nec2-deck`:

```rust
let deck = std::fs::read_to_string("deck.nec")?;
let output = nec2_deck::run(&deck);
print!("{}", output.report);
```

Or the engine itself, with `nec2` and `num-complex`:

```rust
use nec2::{Engine, Excitation, Geometry};
use num_complex::Complex64;

let mut geometry = Geometry::new();
geometry.wire(1, 21, [0.0, 0.0, -0.25], [0.0, 0.0, 0.25], 0.001, 1.0, 1.0)?;
let structure = geometry.connect(0)?;

let mut engine = Engine::new(&structure);
engine.set_frequency(299.8);
engine.fill_and_factor(1.0, false)?;
engine.set_sources(&[(11, Complex64::new(1.0, 0.0))], &[]);
engine.excite(&Excitation::Sources)?;
let report = engine.solve(&[], false, false)?;
println!("{}", report.sources[0].impedance);
```

The engine API follows the routines of NEC-2 and will change before version 1.0.

## Tests

```sh
cargo test
```

The reference reports are stored in `tests/fixtures`. To regenerate them:

```sh
scripts/build-reference.sh   # downloads nec2c 1.3.1 and builds it, stock and fixed
scripts/gen-fixtures.sh
```

`cargo test -p nec2-deck --test fixtures survey -- --ignored --nocapture` prints how every deck compares.

`crates/nec2-deck/tests/fuzz.rs` runs the decks after damaging them at random, alone or against nec2c. Its first lines say how.

## Where the code comes from

See [PROVENANCE.md](PROVENANCE.md).
