# Where nec2-rs differs from nec2c 1.3.1

nec2-rs is ported from nec2c 1.3.1, the C translation of NEC-2. In a few places nec2c departs from the NEC-2 Fortran it was translated from. There nec2-rs follows the Fortran (`nec2dxs.f`, public domain).

## Compatibility rule

nec2-rs must never be worse than the programs it replaces: for the same deck it gives the same result or a better one. NEC-2 here means the Fortran original and nec2c 1.3.1.

- A deck that NEC-2 solves with a result that means something gives the same numbers in nec2-rs, in the report layout of nec2c.
- Where NEC-2 is wrong, nec2-rs gives the better result. Every such place is listed in this document.
- Where the model or a card deserves no trust but the result of NEC-2 still means something, nec2-rs gives that result and adds a line that begins with `WARNING:`. Without these lines the report has the layout of nec2c.
- Where NEC-2 stops, crashes, does not end, reads or writes outside its arrays, or goes on to print values that are not numbers or numbers that mean nothing, nec2-rs stops with an error that names the cause. A clear error is the better result. The section "Errors where NEC-2 goes on" lists these places. A run of nec2-rs that ends without an error has no value in its report that is not a number.

The decks `s7_*` in `tests/fixtures/decks` hold the rule: each is a deck out of the ordinary that NEC-2 solves, and its report is compared number by number with nec2c. `crates/nec2-deck/tests/compatibility.rs` has the errors.

The rule is also checked at random. `crates/nec2-deck/tests/fuzz.rs` damages the test decks in a few places (another number, a card removed, repeated, moved or added), runs them and looks for panics, for runs that do not end and for reports with values that are not numbers. Given a nec2c program, it runs every damaged deck there too: where nec2c ends with a report of finite numbers, nec2-rs must print the same numbers or stop with one of the errors listed below.

In a run of 20 000 damaged decks no deck made nec2-rs panic, run without end or print a value that is not a number. nec2c ended with finite numbers for 7110 of them. nec2-rs printed the same numbers for 6411 and stopped with one of those errors for 143; 480 were left out, because nec2-rs called the matrix ill-conditioned or the deck held a number above 1e20 or below 1e-20. The other 76 differed in the ways that this document describes (most of all pattern angles beyond the horizon over a ground, and requests that follow a loop over incidence angles) or in rounding errors: quantities that are zero in theory, and structures a millionth of a wavelength long.

## Results that differ where nec2c departs from the Fortran

Each difference below is also a patch in `scripts/patch-nec2c.py`. The script produces the "fixed" nec2c that the affected test decks are compared with (`tests/fixtures/nec2c-fixed`). All other decks are compared with the stock nec2c 1.3.1 (`tests/fixtures/nec2c`), and a test checks that the fixes do not change those reports.

### 1. Wire conductivity: `zint` (LD 5)

`zint` has three formulas for the internal impedance of a round wire, chosen by the ratio of the radius to the skin depth. In nec2c 1.3.1 the two `return` statements after the first two formulas are missing, so the function always ends with the third one, which is valid only when the skin depth is much smaller than the radius. The constant `cc5` also has the exponent `e4` instead of `e-4`.

The loss of thin or poorly conducting wires is too small in nec2c. Example (`s6_ld_wire_thin`, an 80 m band dipole of 0.2 mm steel wire):

| | Input impedance | Efficiency |
|---|---|---|
| nec2c 1.3.1 | 198 + j72 Ω | 36.4 % |
| nec2-rs, fixed nec2c | 688 − j129 Ω | 9.6 % |

The DC resistance of this wire is 1273 Ω, about half of which appears at the feed point of a half-wave dipole, so the second line is the plausible one. For ordinary antenna wire at HF the difference is small (`s6_ld_wire_hf`: 69.900 Ω against 69.907 Ω).

### 2. Sommerfeld integrals: integration path in `evlua`

For field points far from the source and close to the ground, the Fortran chooses between two integration paths. When the test `RMIS > 4*RHO/ZPH` holds it takes the path below the branch points. nec2c 1.3.1 takes neither path in that case and leaves the integral incomplete.

The affected grid points lie inside the interpolation tables when the ground is a good conductor. Example (`s6_som_sea_array`, two dipoles 0.2 and 0.35 wavelength above sea water):

| | Input impedance |
|---|---|
| nec2c 1.3.1, Sommerfeld | 97.35 + j89.22 Ω |
| nec2-rs and fixed nec2c, Sommerfeld | 97.84 + j83.40 Ω |
| NEC2++ (black-box run), Sommerfeld | 97.85 + j83.42 Ω |
| any of them, reflection coefficients | 97.69 + j84.09 Ω |

Over such a ground the Sommerfeld result should be close to the reflection-coefficient result, as the fixed value is. Over average ground (permittivity 13, conductivity 5 mS/m) the affected points are outside the tables and nothing changes.

### 3. Sommerfeld integrals: convergence test in `gshank`

The Fortran measures a complex difference as `|re| + |im|`. nec2c computes `|re + |im||`, which is smaller when the real part is negative, so the integration can stop too early. The effect is a few units of the fifth significant digit.

When the series does not converge in 20 steps, the Fortran prints a message and uses the last estimate; nec2c aborts the program. nec2-rs uses the last estimate.

### 4. Sommerfeld tables in a frequency sweep: `intrp`

`intrp` keeps the polynomial coefficients of the table region it used last. nec2c recomputes the tables at every frequency but does not discard the kept coefficients, so the first interactions at a new frequency can be interpolated from the tables of the previous one. nec2-rs discards them. The effect is one unit of the fifth significant digit in the decks tried (`s4_som_hf_dipole`). The test `sommerfeld_sweep_equals_single_frequencies` checks that a sweep gives what single-frequency runs give.

### 5. Networks: `netwk`

- **Loss accounting.** With several sources, one on a network port and another not, nec2c adds the power of the second source to the network loss. Example (`s6_tl_two_sources`, lossless line): nec2c reports a network loss of 5.95 mW, nec2-rs reports zero.
- **Asymmetry check without networks.** The EX option that prints the asymmetry of the admittance matrix makes nec2c 1.3.1 crash (segmentation fault) when the deck has no networks (`s6_asymmetry`). nec2-rs computes it.
- **Network solution between excitations.** NEC-2 keeps the factored network equations while the same structure is excited again (the incidence angles of a plane wave). nec2c frees them after every call, so from the second angle on it prints the currents of the structure without its networks. nec2-rs keeps them. The fixed nec2c does not include this change; the test `networks_act_at_every_incidence_angle` checks instead that the second angle of a sweep gives what a run at that angle alone gives.

### 6. Extended thin-wire kernel flags in `qdsrc` and `nefld`

The copies of the end-condition logic in `qdsrc` (EX 5) and `nefld` (NE) assign `ind1` where the Fortran assigns `ind2`, and `qdsrc` reads one connection entry with a wrong index. nec2-rs uses the logic of `cmww` everywhere. The decks `s6_ex5_ek` and `s6_near_ek` exercise these routines; their results happen to be the same with and without the fix.

### 7. Ground-field integration: `rom2`

The Fortran limits the halving of the integration step by the constant 65536. nec2c uses the number of segments instead, so small models cannot refine the integration of the Sommerfeld field. Example (`cr051`, a 2-segment wire 0.0005 wavelength above ground): nec2c gives 134.95 − j1258.6 Ω with "step size limited" messages, nec2-rs gives 127.96 − j1332.4 Ω.

### 8. Current source on the centre of a segment or of a patch: `etmns`

For EX 4 the entry of the right-hand side of a segment whose centre is the source point is not written and keeps the current of the previous solution. nec2-rs sets it to zero, and both entries of a patch. The Fortran has the same gap.

With patches nec2c has a second fault here: when a patch is skipped, its counter of the patches is not advanced, and every later patch is computed with the geometry of the patch before it. Example (`s8_patch_current_on_patch`, three patches in a row, the source on the middle one): nec2c gives no current on the third patch, nec2-rs the mirror image of the current on the first, 6.3899 A/m.

### 9. Wire on a patch: `cmsw`

Where a wire ends on a patch, a matrix element needs the coefficients of the basis function of the wire's end segment. The Fortran takes them after it has computed the basis functions on that segment (`AX(JSNO)` after `CALL TRIO`). nec2c fixes the index before, from the segment that was handled last, which is the last segment of the structure. If the two segments carry a different number of basis functions, the wrong coefficient is used, and the result depends on the order of the wire cards. Example (a one-segment stub on a plate, and a second wire of three segments elsewhere):

| | Input impedance of the stub |
|---|---|
| nec2c 1.3.1, stub entered first (`s8_patch_stub_first`) | 3.37 − j362.4 Ω |
| nec2c 1.3.1, stub entered last (`s8_patch_stub_last`) | 3.39 − j1091.0 Ω |
| nec2-rs and fixed nec2c, either order | 3.39 − j1091.0 Ω |

### 10. Strip of patches: SC cards in `datagn`

After a rectangular or quadrilateral patch, further SC cards add patches to a strip. nec2c tests for the wrong card when it decides whether an SC card is allowed, and stops at the first SC card of a strip with `PATCH DATA ERROR`. nec2-rs accepts the strip (`s8_patch_strip`).

## Behaviour that differs

- **Job without comment cards after NX.** nec2c rewinds the whole input file and runs the first job again. nec2-rs continues with the next job.
- **EX card after an XQ card with a pattern.** In NEC-2 a request closes the group of EX cards before it, so that the next EX card starts a new excitation. An XQ card that asks for a pattern (XQ 1 to 3) does not: the EX cards after it join the old group, and the next request prints the pattern of the old currents (voltage sources) or goes on with a loop over incidence angles that has ended, with writes beyond the buffer of a normalized receiving pattern. In nec2-rs every request closes the group.
- **Cone angle of a spiral that closes to a point.** For a GH card whose second radius is zero the formula of NEC-2 divides zero by zero, and the report shows a value that is not a number. nec2-rs prints the angle. The structure is the same.
- **Constant of the Hankel function.** nec2c has `-0.02457850915` for (1 − 2γ)/(2π); the Fortran value `-0.0245785095` is the correct one. The difference is beyond the printed precision.
- **Plot file.** `Output::plot` in `nec2-deck` holds what nec2c writes to its plot file for PL cards, and the command-line program writes it to the same file name (the deck name with `.plt` added).
- **Report.** The banner names nec2-rs, and run times are printed as zero. The line with the run times of the matrix also carries the ratio of the smallest to the largest pivot (see below).

## Errors where NEC-2 goes on

NEC-2 runs each of these decks to the end, or crashes on the way. What it prints is not a result: values that are not numbers, the contents of memory that does not belong to the table it reads, or numbers that belong to no physical problem. nec2-rs stops with an error that names the cause. The report holds everything up to the error.

- **A number too large for the machine**, such as `1e400`. nec2c reads it as infinity, prints `INF` for the card and goes on, in most fields to values that are not numbers. Where infinity would have a meaning the card is not needed: a wire without loss has no LD card, a perfect ground is `GN 1`.
- **No source.** A request (XQ, RP, NE, NH) without an EX card before it. NEC-2 prints currents of zero, and an efficiency and gains that are not numbers. nec2-rs prints the loading and the environment of the antenna, which do not depend on a source, and stops.
- **Excitation types.** An EX type other than 0 to 5. nec2c reads beyond its table of names to print the type, crashes if the type is large, and otherwise takes the card for a current source whose position is the angles of the card.
- **Pattern modes above 6.** NEC-2 computes them as mode 6. nec2c prints a null pointer as the name of the cliff, the Fortran whatever the name held before.
- **A cliff without a medium beyond it.** RP 2, 3, 5 and 6 over a ground need the second medium, from a GD card or from the GN card. Without it NEC-2 prints values that are not numbers for every direction whose reflection point lies beyond the cliff. Without any ground the modes give the pattern of free space, in NEC-2 and in nec2-rs.
- **The ground wave at the origin.** RP 1 with no distance and a field point at height zero: the field has the factor 1/R with R = 0.
- **A receiving pattern without an incident wave.** The PT options 1 to 3 print the currents as a receiving pattern. With a voltage source or a current element nec2c reads beyond its table of names for the heading and, for options 2 and 3, writes to a buffer that it has not allocated.
- **Counts below zero.** Segments on GW, GA and GH, repetitions on GM, a rotation count below one on GR, radial wires on GN, frequencies on FR, angles on RP and on EX 1 to 3, field points on NE and NH. NEC-2 damages the structure (a GW card with a count of -1 removes the last segment of the wire before it), prints tables without rows, or computes one frequency.
- **Load types.** A negative load type other than the leading `LD -1`. nec2c uses a stale impedance.
- **Patch shapes.** A shape on SP other than 0 to 3. nec2c reads beyond its table of names. For a shape above 3 it builds a patch that the manual does not define, for a shape below zero one from the corners of an earlier card.
- **Segments that do not exist.** A source, a transmission line, a network or a coupling port on a segment beyond the last one. nec2c reads outside its arrays and prints numbers that are not finite.
- **A frequency that is not positive**, and a sweep that reaches zero or goes below it. At zero nec2c prints values that are not numbers. Below zero it does not end if the structure has wires, and for patches alone it prints the numbers of a negative wavelength.
- **A finite ground without parameters.** A GN card of a type other than 1 with neither a dielectric constant nor a conductivity. NEC-2 divides by zero in the reflection coefficients; for a type above 2 it does so only in the fields.
- **A transmission line of zero impedance, a helix with zero turn spacing, a deck without segments.** NEC-2 prints values that are not numbers or fails.
- **A wire without length** that lies on a junction. NEC-2 refuses segments shorter than 1e-20, but it tests them after it has moved the segment ends to the junctions, where rounding gives such segments a length of some 1e-17. They pass and enter the matrix: the monopole with four radials of `s4_som_gp_radials` and a fifth wire of no length at its base has an input resistance of -5 ohms in nec2c. nec2-rs tests the segments as they were entered.
- **Segment ends within the connection tolerance that form no junction.** The tolerance is a thousandth of the length of the segment that is tested, so a very long segment takes ends for its neighbours that are far from each other. In nec2c the search around such a junction does not end, or it ends with ends joined that are far from each other.
- **Any other result that is not a number.** The lines that the structure data and every request print are searched for a value that is not a number, and the run ends if there is one. This is the last line of defence, for causes that have no error of their own: a matrix so ill-conditioned that the currents are zero, or ground parameters beyond what the formulas can hold.

## From the clean-room specification

These follow `cleanroom.md` (see PROVENANCE.md) and are tested in `crates/nec2-deck/tests/cleanroom.rs` and in the unit tests of the engine. nec2c 1.3.1 and the NEC-2 Fortran behave alike in most of these places, so there is no corrected nec2c to compare with. The tests use the numbers of the specification. Where the specification asks for an error but NEC-2 gives a finite result, the compatibility rule decides: nec2-rs gives the result of NEC-2 and a warning.

### Input handling

- A last line without a newline is a complete line, and the end of the file ends the job like EN. The end of the file before GE is an error.
- Indented cards are read; nec2c drops them silently.
- Lines have no length limit. A line that starts with an apostrophe is a comment, and text after `!` or an apostrophe on a card is ignored. A comment card after CE is ignored with a warning.
- An integer field may be written as a real number with an integral value; `D` is accepted as exponent letter. A fraction in an integer field stays an error.
- Errors name the line and the card. Geometry cards after GE, control cards before GE and unknown cards get separate messages.
- Blank counts mean one: incidence angles on EX 1 to 3, points on NE and NH, angles on RP.
- GR with one sector changes nothing. GH with zero length is a flat spiral whose first real field is the number of segments per turn.
- After a finished loop over more than one theta angle of an incident wave, a further XQ or RP starts the loop again; NEC-2 goes on to an angle beyond the last one and writes outside its buffers. With a single theta NEC-2 serves the request with the currents it has, and so does nec2-rs.

### Warnings where NEC-2 says nothing

NEC-2 solves all of these, and nec2-rs gives the same numbers.

- **Ground flag without a ground.** `GE 1` joins wire ends on the plane z = 0 to their image. Without a GN card the run is a free-space run, and the result belongs to no physical problem. A flag that joins nothing has no effect and gets no warning.
- **Wires below the ground or in it.** With `GE 1` NEC-2 refuses them, and so does nec2-rs. With `GE 0` or `GE -1` and a GN card NEC-2 does not look.
- **Wires in the same space.** A segment whose centre lies inside the volume of another wire. NEC-2 Part III, section II.1 says that such a model gives meaningless results. Two segments that meet at exactly one end form a junction and are not tested against each other. `Options::overlap` in `nec2-deck`, or `--overlap error` and `--overlap off` on the command line, turn the warnings into an error or switch the test off.
- **A wire without segments.** GW, GA and GH with a segment count of zero are ignored.
- **Ranges of segments.** On LD, PT and PQ a range that reaches beyond the last segment selects the segments that exist, without a warning. A range that selects no segment has no effect and gets a warning.
- **Ground types that are not defined.** For a type above 2 NEC-2 fills the matrix as for a perfect ground and uses the reflection coefficients of the ground in the fields. A type below -1 is taken as 0.
- **Pattern modes for a radial-wire screen without a screen.** RP 4, 5 and 6 are computed without the screen. nec2c prints not-a-number for a direction in which the image point of a segment lies on the axis; nec2-rs prints the value of the plain ground there.
- **A patch without area.** It carries a current but radiates nothing.
- **Coupling card for a segment that does not exist.** NEC-2 stops only when it computes the coupling, which needs a voltage source on one of the listed segments.
- **Negative reflection flags on GX.** NEC-2 ignores the sign.
- **Stepping types on FR that are not defined.** For a type other than 0 and 1 NEC-2 adds the step. nec2c then labels the table of input impedances with a number that an earlier computation left behind; nec2-rs prints the frequencies.
- **Pattern modes below zero.** With -1 NEC-2 computes no pattern; any other mode below zero is taken as mode 0.
- **Left-handed helix.** For a negative length NEC-2 builds the right-handed helix and exchanges x and y. The helix is left-handed, but it starts at (0, A, 0) instead of (A, 0, 0), and the radii of the card that are named for x lie along y. A wire that is meant to feed it at (A, 0, 0) stays unconnected. nec2-rs keeps this geometry, because decks written for NEC-2 rely on it, and names the start point in a warning.

Negative tags and negative tag increments on GX, GR and GM are plain labels in NEC-2 and get no warning.

### Checks of the numerical work

- **Segment and patch data.** After GE every segment must have finite ends, a finite length above 1e-20 and a finite positive radius, and every patch a finite centre, finite tangent vectors and a finite area. The error names the segment and its tag, or the patch. In nec2c a coordinate that is not a number passes the check.
- **Wires on patches.** Two wire ends on the centre of the same patch are an error; nec2c reads outside its arrays. A patch without area, or one so small that its quarters lie within the connection tolerance of the wire end on it, is an error; NEC-2 goes on dividing it without end.
- **A normalized receiving pattern without any received current** is printed as zeros; NEC-2 divides zero by zero.
- **Integrations.** `intx`, `hfk`, `rom2` and `rom1` have a bound on the number of steps, and they stop with an error when a value is not finite. The error names the routine and the segments. In nec2c these routines do not end when a limit is not a number.
- **Matrix.** A pivot that is zero or not finite is an error that names the equation and the frequency. The ratio of the smallest to the largest pivot is printed on the line of the matrix run times (`SMALLEST/LARGEST PIVOT`). Sound models give more than 0.01. Below 1e-5 a warning says that the matrix is ill-conditioned. A solution that is not finite is an error. nec2c prints `NAN` in these cases and ends with exit status 0.
- **Size.** A structure for whose matrix there is no memory is an error; nec2c ends with a message of the memory allocation, or it is ended by the system. `Options::max_equations` in `nec2-deck` sets a limit on the number of equations (one for a segment, two for a patch) for programs that run decks they did not write. Without it there is no limit, as in NEC-2.

### Results that differ

- **Negative ground conductivity in a frequency sweep.** A negative conductivity on GN stands for the imaginary part of the permittivity. nec2c and the Fortran turn it into a conductivity at the first frequency and keep that value, so in the rest of the sweep the permittivity is no longer the one on the card. nec2-rs applies the card at every frequency. The first frequency is not affected. Example: `GN 0 0 0 0 12 -3` at 250, 300 and 350 MHz gives 12 − j3 three times; nec2c gives 12 − j3, 12 − j2.5 and 12 − j2.14.
- **Radial-wire screen.** A GN card without a screen clears the screen of an earlier GN card, and the screen formulas are used only when the number of radials is above zero and the ground is not perfect.
- **Directions below the ground.** With a ground, a direction below it is not a sample of the pattern. nec2c leaves θ above 90.01° out of the pattern table but not out of the normalized-gain table, whose angles and values then no longer match, nor out of the solid angle of the average gain. nec2-rs lists only the directions it computed. This also holds for negative θ beyond −90°, for which nec2c prints values.
- **Average power gain.** The gain is averaged over the solid angle that the samples cover. Each sample has the weight Δφ·∫|sin θ|dθ over its cell, which is the weight of nec2c when the cell does not contain θ = 0° or 180°. nec2c divides by the nominal solid angle of the angle ranges, which is zero for a θ range that is symmetric about the zenith, and prints `INF`. When the samples cover no solid angle, nec2-rs prints `UNDEFINED`. A warning follows when the average of an antenna exceeds its limit, 4π divided by the solid angle, by more than 1 %: more power would leave the antenna than enters it. `AverageGain::efficiency` returns the average times the sampled fraction of the sphere. It is not printed, so that the line keeps the layout of nec2c.

### Not taken over

- The left-handed helix as a mirror image in the x-z plane: it would break decks written for NEC-2 (see the warning above). A helix with an explicit winding sense that starts at (A, 0, 0) belongs in the programming interface, not in the GH card.
- Exact physical constants. They would change every result of NEC-2 in the fourth or fifth digit and wait for a decision.
- The switch that brings back the frozen conductivity of nec2c in a sweep with a negative conductivity.

## Not implemented

- Numerical Green's function (GF and WG cards). The Fortran NEC-2 has it: WG saves the factored matrix of a structure to a file, and a later run that starts with GF adds wires and patches to that structure without solving it again. The result is that of one run with the whole structure; the gain is computing time. nec2c removed it on purpose and stops at these cards, and so does nec2-rs. A file written by one program cannot be read by another, so only the two-run way of working is missing.
- The plot flag of the GE card (its second integer). nec2c writes to a plot file that is not open yet and crashes. nec2-rs ignores the flag. The plot data that PL cards ask for is written as in nec2c.
