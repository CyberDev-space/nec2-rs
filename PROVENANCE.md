# Provenance

nec2-rs is released under the MIT license. To keep that possible, the code is written only from public-domain sources, and code under copyleft licenses is kept away from the people and tools that write it. This file records where the code comes from and how that separation is kept.

## Sources the code is ported from

| Source | Status | Use |
|---|---|---|
| nec2c 1.3.1 by N. Kyriazis (5B4AZ) | public domain | the base of the port |
| NEC-2 Fortran, double precision with SOMNEC (`nec2dxs.f`) | public domain | cross-check where nec2c looks wrong |
| NEC-2 manuals, parts I to III (Lawrence Livermore Laboratory, 1981) | public domain | theory, program description, card formats |

- **nec2c 1.3.1.**
  - Archive: `nec2c_1.3.1.orig.tar.bz2` from the Debian pool, <https://deb.debian.org/debian/pool/main/n/nec2c/>.
  - SHA-256: `271bd2ea6037f896912139a66e3abb75981c72315bab5221fcf384a1736d0337`, the value of the Debian source control file.
  - License: section 7 of its README says "nec2c is Public Domain, same as the original FORTRAN source. Please keep any software you write incorporating nec2c in Public Domain or at least use an open license like GPL or BSD." The Debian copyright file records the same.
  - The archive also contains a `COPYING` file with the text of the GNU GPL version 3. It is the default file that the automake tool installs; no source file refers to it. The author's own statement in the README is taken as the license.
- **NEC-2 Fortran.** The file `nec2dxs.f` as distributed at <https://github.com/yeti01/nec2>. NEC-2 was written at Lawrence Livermore Laboratory for the United States government and is in the public domain.
- **Manuals.** <https://www.nec2.org/>.

`scripts/build-reference.sh` downloads the nec2c archive, checks the checksum and builds the reference programs. Nothing of nec2c is stored in this repository except the reports it prints for the test decks.

## Programs used only as black boxes

They are run on test decks and their output is compared; their source is not read.

- nec2c 1.3.1 built from the archive above, and the same with `scripts/patch-nec2c.py` applied. The patches are written for this project from the Fortran original (see `docs/nec2c-deviations.md`).
- NEC2++ through the npm package `@necpp-engine/wasm` (GPL-2.0-or-later).

## Sources that must not be read by whoever writes the code

Translating code into another language creates a derivative work, so the code of these projects is off limits to the "clean side" that writes nec2-rs:

- NEC2++ (`tmolteno/necpp` and its forks), GPL-2.0-or-later;
- xnec2c (`KJ7LNW/xnec2c`), GPL-3.0, including its changes to the nec2c files;
- any other NEC-2 derivative whose license is not public domain or permissive.

Their documentation, change logs, issue discussions and published papers may be read.

## Clean-room procedure for improvements

Improvements that exist in the copyleft projects reach nec2-rs only as a behavioural specification:

1. A separate "dirty side" reads the copyleft code and writes `cleanroom.md`: for each improvement, what nec2c does, what a correct program must do and how to verify it. The document contains no code, no identifiers and no description of how the other code is organized. A script checks this mechanically.
2. The project owner reviews the document before the clean side may read it.
3. The clean side implements each item from the specification and the public-domain sources, with a test that fails before the change and passes after it.

The dirty side works in a separate directory outside this repository, and its notes on where each finding comes from stay there.

Status: the specification was written on 2026-11-21 (62 items). The owner gave the go-ahead to implement it the same day and confirmed that the two lines flagged by `scripts/check-spec.py` contain no code. The clean side reads only `cleanroom.md` itself, not the other files of the dirty side.

## What the clean side has seen of copyleft projects

For the record, before the port started the following was looked at:

- NEC2++: the README and NEWS files, and the license lines at the top of one source file.
- xnec2c: the list of file names in its source directory, the line counts of some files, and a few lines of widget declarations of its optimizer window. A copy made to count lines was deleted without being read.
- `@necpp-engine/wasm`: its README and TypeScript type declarations, to run it as a black box.
- FieldLab: its README.

None of this is engine code.

## Tools

The port was written with an AI coding assistant (Claude Code). A language model may have seen public source code during its training. The safeguards are the ones above: the structure, names and order of operations of every routine are taken from nec2c, and each routine was checked against nec2c's printed results.
