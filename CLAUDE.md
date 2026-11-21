# nec2-rs

Rust port of NEC-2 (see README.md). The license is MIT, and it stays possible only under the rules of PROVENANCE.md.

## Clean-room rules

- Never open, search, fetch or summarize the source code of NEC2++ (necpp), xnec2c or any other copyleft NEC-2 derivative, on disk or on the web. Their documentation and change logs are fine.
- Never read the working directory of the dirty side. It holds copyleft sources and the notes of the dirty side, and it is kept outside this repository. The one document that crosses over is its `cleanroom.md`, and only after the project owner says in the conversation that they have reviewed it.
- Port only from `reference/nec2c-1.3.1` (public domain), the NEC-2 Fortran and the NEC-2 manuals. Where nec2c departs from the Fortran, follow the Fortran and record it in `docs/nec2c-deviations.md` and in `scripts/patch-nec2c.py`.

## Compatibility rule

The owner's rule: nec2-rs must never be worse than the programs it replaces. For the same deck it gives the same result as NEC-2 (the Fortran and nec2c) or a better one, and a clear error counts as better than garbage.

- NEC-2 solves the deck with a result that means something: print the same numbers, in the report layout of nec2c. Where the model or a card deserves no trust, add a line that begins with `WARNING:`. Add the deck to `tests/fixtures/decks` as `s7_*`.
- NEC-2 stops, crashes, does not end, reads or writes outside its arrays, or goes on to print values that are not numbers or numbers that mean nothing: stop with an error that names the cause. Add it to "Errors where NEC-2 goes on" in `docs/nec2c-deviations.md`, to `ALLOWED_ERRORS` in `crates/nec2-deck/tests/fuzz.rs` and to `crates/nec2-deck/tests/compatibility.rs`. Never print a value that is not a number.
- Before deciding, run the input through `reference/bin/nec2c` and look at what it really prints.

The first of these wins over an item of the clean-room specification that asks for an error.

After a change to the reading of cards, or to what is an error, run the damaged decks against nec2c, several seeds, and look at every deck it writes as `stricter_*` or `differs_*`:

```sh
NEC2_FUZZ_REFERENCE=$PWD/reference/bin/nec2c-fixed NEC2_FUZZ_SEED=1 cargo test --release -p nec2-deck --test fuzz -- --ignored --nocapture
```

Run it with `ulimit -c 0`, so that no crash leaves a core dump.

## Working rules

- `cargo test`, `cargo fmt --check` and `cargo clippy --all-targets` must be clean before a change is done.
- A change in results needs a deck in `tests/fixtures/decks` and reference reports made by `scripts/gen-fixtures.sh` (run `scripts/build-reference.sh` once before). Where no nec2c can give the reference (an item of the clean-room specification that the Fortran shares), it needs a test in `crates/nec2-deck/tests/cleanroom.rs` with the numbers of the specification instead.
- In `crates/nec2/src/engine` keep the routine names and the order of floating-point operations of nec2c: results are compared with nec2c to the last printed digit.
