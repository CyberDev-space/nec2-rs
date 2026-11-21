//! Interpreter for NEC-2 decks.
//!
//! A deck is a small program: the cards after the geometry section set
//! parameters, and the cards XQ, RP, NE and NH start a computation that
//! reuses whatever earlier results are still valid. This module reproduces
//! the control flow and the report of nec2c 1.3.1 (`main` in `main.c` and
//! `datagn` in `geometry.c`).

use std::fmt::Write as _;

use nec2::consts::{PI, TA, TD, cang, db10, db20};
use nec2::{
    Engine, Excitation, Geometry, GroundParams, GroundType, HelixInfo, Load, Network, NetworkKind, OverlapCheck,
    PatternRequest, PortRecord, RadiationPattern, Structure,
};
use num_complex::Complex64;

use crate::cfmt::{e, e_plus, f, f_plus};
use crate::parse::{CONTROL_CARDS, GEOMETRY_CARDS, Line, Lines, mnemonic, parse_card};

/// Result of running a deck.
#[derive(Debug, Clone, PartialEq)]
pub struct Output {
    /// The report, in the layout of the NEC-2 output file.
    pub report: String,
    /// The message of the error that stopped the run, if any. The message is
    /// also the last thing in the report.
    pub error: Option<String>,
    /// The numbers that PL cards ask for, as nec2c writes them to its plot
    /// file: one line of columns for every current, field point or pattern
    /// direction. `None` if the deck has no PL card.
    pub plot: Option<String>,
}

/// Adds a line of numbers to the plot data.
fn plot_line(plot: &mut Option<String>, values: &[f64]) {
    if let Some(text) = plot {
        let columns: Vec<String> = values.iter().map(|&value| e(value, 12, 4)).collect();
        let _ = writeln!(text, "{}", columns.join(" "));
    }
}

/// Why the interpreter stops.
enum Stop {
    /// The deck ended normally.
    End,
    /// An error; the text is appended to the report.
    Error(String),
}

type Flow<T = ()> = Result<T, Stop>;

impl From<nec2::Error> for Stop {
    fn from(error: nec2::Error) -> Self {
        Stop::Error(format!("\n  {error}"))
    }
}

fn fail<T>(message: impl Into<String>) -> Flow<T> {
    Err(Stop::Error(message.into()))
}

/// Letters that mark the shape of a patch in the structure table: any shape,
/// rectangle, triangle, quadrilateral.
const PATCH_SHAPES: [char; 4] = ['P', 'R', 'T', 'Q'];

/// Warning for a wire card without segments, which NEC-2 ignores silently.
const NO_SEGMENTS: &str = "\n  WARNING: THE WIRE ABOVE HAS NO SEGMENTS AND IS IGNORED";

/// Choices for running a deck.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Options {
    /// What to do when wires occupy the same space. The default is to
    /// solve the structure, as NEC-2 does, and to warn.
    pub overlap: OverlapCheck,
    /// The largest number of equations a structure may have (one for every
    /// segment and two for every patch), or `None` for no limit. A program
    /// that runs decks it did not write sets a limit, so that a mistyped
    /// count cannot take all the memory and the time.
    pub max_equations: Option<usize>,
}

/// A value that is not a number, as the report prints it.
fn not_a_number(token: &str) -> bool {
    let bare = token.trim_start_matches(['-', '+']).trim_end_matches('j');
    bare.eq_ignore_ascii_case("nan") || bare.eq_ignore_ascii_case("inf")
}

/// Stops the run if the structure has more equations than the options allow.
fn check_size(segments: usize, patches: usize, options: &Options) -> Flow {
    let equations = segments + 2 * patches;
    match options.max_equations {
        Some(limit) if equations > limit => {
            fail(format!("\n  THE STRUCTURE NEEDS {equations} EQUATIONS, MORE THAN THE LIMIT OF {limit}"))
        }
        _ => Ok(()),
    }
}

/// Runs every job of a deck and returns the report.
pub fn run(input: &str) -> Output {
    run_with(input, &Options::default())
}

/// Runs every job of a deck with the given choices and returns the report.
pub fn run_with(input: &str, options: &Options) -> Output {
    let mut out = String::new();
    let mut plot = None;
    let mut lines = Lines::new(input);
    let mut first_job = true;
    let error = loop {
        match run_job(&mut lines, &mut out, &mut plot, first_job, options) {
            Ok(()) => first_job = false,
            Err(Stop::End) => break None,
            Err(Stop::Error(message)) => {
                // Name the card that was being processed.
                let message = match &lines.last {
                    Some(line) => format!("\n\n  ERROR AT LINE {}: {}{message}", line.number, line.text),
                    None => message,
                };
                out.push_str(&message);
                break Some(message.trim().to_string());
            }
        }
    };
    Output { report: out, error, plot }
}

/// Polarization names used in the excitation printouts.
const HPOL: [&str; 3] = ["LINEAR", "RIGHT", "LEFT"];

/// Runs one job: comments, geometry and program cards up to NX or EN.
fn run_job(lines: &mut Lines, out: &mut String, plot: &mut Option<String>, first_job: bool, options: &Options) -> Flow {
    let Some(first) = lines.next_line() else {
        // Nothing after NX: the deck is complete.
        return if first_job { fail("\n\n  ERROR: THE INPUT CONTAINS NO CARDS") } else { Err(Stop::End) };
    };
    out.push_str(concat!(
        "\n\n\n",
        "                               __________________________________________\n",
        "                              |                                          |\n",
        "                              |  NUMERICAL ELECTROMAGNETICS CODE (NEC-2) |\n",
        "                              |      nec2-rs, a Rust implementation      |\n",
        "                              |__________________________________________|\n",
    ));

    // Comment cards.
    let mut pending = None;
    let mut ain = mnemonic(&first.text).to_string();
    if ain == "CM" || ain == "CE" {
        out.push_str("\n\n\n                               ---------------- COMMENTS ----------------\n");
        let _ = writeln!(out, "                              {}", first.text.get(2..).unwrap_or(""));
        while ain == "CM" {
            let Some(line) = lines.next_line() else {
                return fail("\n\n  ERROR: UNEXPECTED END OF INPUT");
            };
            ain = mnemonic(&line.text).to_string();
            let _ = writeln!(out, "                              {}", line.text.get(2..).unwrap_or(""));
        }
        if ain != "CE" {
            return fail("\n\n  ERROR: INCORRECT LABEL FOR A COMMENT CARD");
        }
    } else {
        pending = Some(first);
    }

    let from = out.len();
    let structure = read_geometry(lines, pending, out, options)?;
    if out[from..].split_whitespace().any(not_a_number) {
        return fail("\n\n  THE STRUCTURE DATA ABOVE HOLD VALUES THAT ARE NOT NUMBERS");
    }
    out.push_str("\n\n\n");
    Job::new(&structure, out, plot).run(lines)
}

/// Reads the geometry section up to the GE card and prints the structure
/// specification and segmentation data.
fn read_geometry(lines: &mut Lines, mut pending: Option<Line>, out: &mut String, options: &Options) -> Flow<Structure> {
    let mut geometry = Geometry::new();
    let mut nwire = 0;
    let mut header_printed = false;
    // Corners 3 and 4 of the last patch and whether an SC card may add a
    // patch to it.
    let (mut c3, mut c4) = ([0.0; 3], [0.0; 3]);
    let mut isct = false;

    let ground_flag = loop {
        let Some(line) = pending.take().or_else(|| lines.next_line()) else {
            return fail("\n\n  ERROR: THE GEOMETRY SECTION IS NOT TERMINATED BY A GE CARD");
        };
        let line = line.text;
        let name = mnemonic(&line);
        if name == "CM" || name == "CE" {
            out.push_str("\n  WARNING: COMMENT CARD AFTER CE IS IGNORED");
            continue;
        }
        if CONTROL_CARDS.contains(&name) {
            return fail("\n  PROGRAM CONTROL CARD IN THE GEOMETRY SECTION (BEFORE THE GE CARD)");
        }
        if !GEOMETRY_CARDS.contains(&name) && name != "XT" {
            return fail("\n  UNKNOWN CARD IN THE GEOMETRY SECTION");
        }
        let card = parse_card(&line, 2, 7, "GEOMETRY").map_err(Stop::Error)?;
        if card.mnemonic == "XT" {
            out.push_str("\n\n  Exiting after an \"XT\" command");
            return Err(Stop::End);
        }
        if !header_printed {
            out.push_str(concat!(
                "\n\n\n",
                "                               -------- STRUCTURE SPECIFICATION --------\n",
                "                                     COORDINATES MUST BE INPUT IN\n",
                "                                     METERS OR BE SCALED TO METERS\n",
                "                                     BEFORE STRUCTURE INPUT IS ENDED\n",
                "\n",
                "  WIRE                                           ",
                "                                      SEG FIRST  LAST  TAG\n",
                "   No:        X1         Y1         Z1         X2      ",
                "   Y2         Z2       RADIUS   No:   SEG   SEG  No:",
            ));
            header_printed = true;
        }

        let (itg, ns) = (card.ints[0], card.ints[1]);
        let [xw1, yw1, zw1, xw2, yw2, zw2, rad] = card.floats[..] else { unreachable!("seven fields are parsed") };
        // A strip of patches goes on only with SC cards. (nec2c tests for
        // the wrong card here and refuses every SC card of a strip.)
        if card.mnemonic != "SC" {
            isct = false;
        }
        let i1 = geometry.segment_count() as i32 + 1;
        let i2 = geometry.segment_count() as i32 + ns;

        match card.mnemonic.as_str() {
            "GW" => {
                nwire += 1;
                let _ = write!(
                    out,
                    "\n {nwire:5}  {} {} {} {} {} {} {} {ns:5} {i1:5} {i2:5} {itg:4}",
                    f(xw1, 10, 5),
                    f(yw1, 10, 5),
                    f(zw1, 10, 5),
                    f(xw2, 10, 5),
                    f(yw2, 10, 5),
                    f(zw2, 10, 5),
                    f(rad, 10, 5)
                );
                let (mut rad, mut rdel, mut rrad) = (rad, 1.0, 1.0);
                if rad == 0.0 {
                    // A tapered wire: the next card gives the taper.
                    let line = lines.next_line().map(|line| line.text).unwrap_or_default();
                    let gc = parse_card(&line, 2, 7, "GEOMETRY").map_err(Stop::Error)?;
                    if gc.mnemonic != "GC" {
                        return fail("\n  GEOMETRY DATA CARD ERROR");
                    }
                    let (xs1, ys1, zs1) = (gc.floats[0], gc.floats[1], gc.floats[2]);
                    let _ = write!(
                        out,
                        "\n  ABOVE WIRE IS TAPERED.  SEGMENT LENGTH RATIO: {}\n                                 RADIUS FROM: {} TO: {}",
                        f(xs1, 9, 5),
                        f(ys1, 9, 5),
                        f(zs1, 9, 5)
                    );
                    if ys1 == 0.0 || zs1 == 0.0 {
                        return fail("\n  GEOMETRY DATA CARD ERROR");
                    }
                    rad = ys1;
                    rdel = xs1;
                    rrad = (zs1 / ys1).powf(1.0 / (f64::from(ns) - 1.0));
                }
                geometry.wire(itg, ns, [xw1, yw1, zw1], [xw2, yw2, zw2], rad, rdel, rrad)?;
                if ns == 0 {
                    out.push_str(NO_SEGMENTS);
                }
            }

            "GX" => {
                let mut iy = ns / 10;
                let iz = ns - iy * 10;
                let ix = iy / 10;
                iy -= ix * 10;
                let flag = |on: bool, axis: char| if on { axis } else { '*' };
                let _ = write!(
                    out,
                    "\n      STRUCTURE REFLECTED ALONG THE AXES {} {} {} - TAGS INCREMENTED BY {itg}",
                    flag(ix != 0, 'X'),
                    flag(iy != 0, 'Y'),
                    flag(iz != 0, 'Z')
                );
                if ns < 0 {
                    let _ = write!(
                        out,
                        "\n  WARNING: THE REFLECTION FLAGS {ns} ARE NEGATIVE - AS IN NEC-2, THE SIGN IS IGNORED"
                    );
                }
                if ns != 0 {
                    geometry.reflect(ix != 0, iy != 0, iz != 0, itg)?;
                } else {
                    geometry.reflect(false, false, false, itg)?;
                }
            }

            "GR" => {
                let _ = write!(out, "\n  STRUCTURE ROTATED ABOUT Z-AXIS {ns} TIMES - LABELS INCREMENTED BY {itg}");
                geometry.rotate(itg, ns)?;
            }

            "GS" => {
                geometry.scale(xw1);
                let _ = write!(out, "\n     STRUCTURE SCALED BY FACTOR: {}", f(xw1, 10, 5));
            }

            "GE" => break itg,

            "GM" => {
                let _ = write!(
                    out,
                    "\n     THE STRUCTURE HAS BEEN MOVED, MOVE DATA CARD IS:\n   {itg:3} {ns:5} {} {} {} {} {} {} {}",
                    f(xw1, 10, 5),
                    f(yw1, 10, 5),
                    f(zw1, 10, 5),
                    f(xw2, 10, 5),
                    f(yw2, 10, 5),
                    f(zw2, 10, 5),
                    f(rad, 10, 5)
                );
                geometry.move_structure(xw1 * TA, yw1 * TA, zw1 * TA, xw2, yw2, zw2, (rad + 0.5) as i32, ns, itg)?;
            }

            "GA" => {
                nwire += 1;
                let _ = write!(
                    out,
                    "\n {nwire:5} ARC RADIUS: {}  FROM: {} TO: {} DEGREES       {} {ns:5} {i1:5} {i2:5} {itg:4}",
                    f(xw1, 9, 5),
                    f(yw1, 8, 3),
                    f(zw1, 8, 3),
                    f(xw2, 11, 5)
                );
                geometry.arc(itg, ns, xw1, yw1, zw1, xw2)?;
                if ns == 0 {
                    out.push_str(NO_SEGMENTS);
                }
            }

            "GH" => {
                nwire += 1;
                let _ = write!(
                    out,
                    "\n {nwire:5} HELIX STRUCTURE - SPACING OF TURNS: {} AXIAL LENGTH: {}  {} {ns:5} {i1:5} {i2:5} {itg:4}\n       RADIUS X1:{} Y1:{} X2:{} Y2:{} ",
                    f(xw1, 8, 3),
                    f(yw1, 8, 3),
                    f(rad, 8, 3),
                    f(zw1, 8, 3),
                    f(xw2, 8, 3),
                    f(yw2, 8, 3),
                    f(zw2, 8, 3)
                );
                match geometry.helix(itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad)? {
                    Some(HelixInfo::Spiral { cone_angle }) => {
                        let _ = write!(out, "\n       THE CONE ANGLE OF THE SPIRAL IS {}", f(cone_angle, 10, 4));
                    }
                    Some(HelixInfo::Helix { pitch_angle, turn_length }) => {
                        let _ = write!(
                            out,
                            "\n       THE PITCH ANGLE IS: {pitch_angle:.4}    THE LENGTH OF WIRE/TURN IS: {turn_length:.4}"
                        );
                    }
                    None => {}
                }
                if ns == 0 {
                    out.push_str(NO_SEGMENTS);
                } else if yw1 < 0.0 {
                    // The start point and the radii of a left-handed helix.
                    let _ = write!(
                        out,
                        concat!(
                            "\n  WARNING: A LEFT-HANDED HELIX (NEGATIVE LENGTH) IS BUILT AS IN NEC-2, WITH X AND Y ",
                            "EXCHANGED: IT STARTS AT X = 0, Y = {}, AND THE X RADII OF THE CARD LIE ALONG Y"
                        ),
                        f(zw1, 0, 5).trim()
                    );
                }
            }

            // One patch: of any shape, a rectangle, a triangle or a quadrilateral.
            "SP" => {
                let shape = ns + 1;
                if itg != 0 || !(1..=4).contains(&shape) {
                    return fail("\n  PATCH DATA ERROR");
                }
                let _ = write!(
                    out,
                    "\n {:5}{} {} {} {} {} {} {}",
                    geometry.patch_count() + 1,
                    PATCH_SHAPES[(shape - 1) as usize],
                    f(xw1, 10, 5),
                    f(yw1, 10, 5),
                    f(zw1, 10, 5),
                    f(xw2, 10, 5),
                    f(yw2, 10, 5),
                    f(zw2, 10, 5)
                );
                if shape == 2 || shape == 4 {
                    isct = true;
                }
                if shape > 1 {
                    // The next card gives the other corners.
                    let line = lines.next_line().map(|line| line.text).unwrap_or_default();
                    let sc = parse_card(&line, 2, 7, "GEOMETRY").map_err(Stop::Error)?;
                    c3 = [sc.floats[0], sc.floats[1], sc.floats[2]];
                    c4 = [sc.floats[3], sc.floats[4], sc.floats[5]];
                    if shape == 2 {
                        c4 = [xw1 + c3[0] - xw2, yw1 + c3[1] - yw2, zw1 + c3[2] - zw2];
                    }
                    let _ = write!(
                        out,
                        "\n      {} {} {} {} {} {}",
                        f(c3[0], 11, 5),
                        f(c3[1], 11, 5),
                        f(c3[2], 11, 5),
                        f(c4[0], 11, 5),
                        f(c4[1], 11, 5),
                        f(c4[2], 11, 5)
                    );
                    if sc.mnemonic != "SC" {
                        return fail("\n  PATCH DATA ERROR");
                    }
                    geometry.patch(0, shape, [xw1, yw1, zw1], [xw2, yw2, zw2], c3, c4)?;
                } else {
                    geometry.patch(0, shape, [xw1, yw1, zw1], [xw2 * TA, yw2 * TA, zw2], c3, c4)?;
                }
            }

            // A rectangle covered by patches.
            "SM" => {
                let _ = write!(
                    out,
                    "\n {:5}R {} {} {} {} {} {}     SURFACE - {itg} BY {ns} PATCHES",
                    geometry.patch_count() + 1,
                    f(xw1, 10, 5),
                    f(yw1, 11, 5),
                    f(zw1, 11, 5),
                    f(xw2, 11, 5),
                    f(yw2, 11, 5),
                    f(zw2, 11, 5)
                );
                if itg < 1 || ns < 1 {
                    return fail("\n  PATCH DATA ERROR");
                }
                let line = lines.next_line().map(|line| line.text).unwrap_or_default();
                let sc = parse_card(&line, 2, 7, "GEOMETRY").map_err(Stop::Error)?;
                c3 = [sc.floats[0], sc.floats[1], sc.floats[2]];
                c4 = [xw1 + c3[0] - xw2, yw1 + c3[1] - yw2, zw1 + c3[2] - zw2];
                let _ = write!(
                    out,
                    "\n      {} {} {} {} {} {}",
                    f(c3[0], 11, 5),
                    f(c3[1], 11, 5),
                    f(c3[2], 11, 5),
                    f(c4[0], 11, 5),
                    f(c4[1], 11, 5),
                    f(c4[2], 11, 5)
                );
                if sc.mnemonic != "SC" {
                    return fail("\n  PATCH DATA ERROR");
                }
                geometry.patch(itg, ns, [xw1, yw1, zw1], [xw2, yw2, zw2], c3, c4)?;
            }

            // One more patch on the side of the last one that joins its
            // corners 4 and 3.
            "SC" => {
                let shape = ns + 1;
                if !isct || itg != 0 || (shape != 2 && shape != 4) {
                    return fail("\n  PATCH DATA ERROR");
                }
                let (c1, c2) = (c4, c3);
                c3 = [xw1, yw1, zw1];
                c4 = if shape == 4 {
                    [xw2, yw2, zw2]
                } else {
                    [c1[0] + c3[0] - c2[0], c1[1] + c3[1] - c2[1], c1[2] + c3[2] - c2[2]]
                };
                let _ = write!(
                    out,
                    "\n {:5}{} {} {} {} {} {} {}\n      {} {} {}  {} {} {}",
                    geometry.patch_count() + 1,
                    PATCH_SHAPES[(shape - 1) as usize],
                    f(c1[0], 10, 5),
                    f(c1[1], 11, 5),
                    f(c1[2], 11, 5),
                    f(c2[0], 11, 5),
                    f(c2[1], 11, 5),
                    f(c2[2], 11, 5),
                    f(c3[0], 11, 5),
                    f(c3[1], 11, 5),
                    f(c3[2], 11, 5),
                    f(c4[0], 11, 5),
                    f(c4[1], 11, 5),
                    f(c4[2], 11, 5)
                );
                geometry.patch(0, shape, c1, c2, c3, c4)?;
            }
            "GF" => return Err(nec2::Error::Unsupported("NUMERICAL GREEN'S FUNCTION (GF CARD)").into()),

            other => {
                return fail(format!(
                    "\n  GEOMETRY DATA CARD ERROR\n {other:>2} {itg:3} {ns:5} {} {} {} {} {} {} {}",
                    f(xw1, 10, 5),
                    f(yw1, 10, 5),
                    f(zw1, 10, 5),
                    f(xw2, 10, 5),
                    f(yw2, 10, 5),
                    f(zw2, 10, 5),
                    f(rad, 10, 5)
                ));
            }
        }
        check_size(geometry.segment_count(), geometry.patch_count(), options)?;
    };

    if ground_flag != 0 {
        out.push_str("\n\n     GROUND PLANE SPECIFIED.");
        if ground_flag > 0 {
            out.push_str(
                "\n     WHERE WIRE ENDS TOUCH GROUND, CURRENT WILL BE INTERPOLATED TO IMAGE IN GROUND PLANE.\n",
            );
        }
    }
    let structure = geometry.connect_with(ground_flag, options.overlap)?;
    check_size(structure.n(), structure.m(), options)?;
    for warning in &structure.warnings {
        let _ = write!(out, "\n  WARNING: {warning}");
    }
    let n = structure.n();
    let _ = write!(
        out,
        "\n\n     TOTAL SEGMENTS USED: {n}   SEGMENTS IN A SYMMETRIC CELL: {}   SYMMETRY FLAG: {}",
        structure.np, structure.ipsym
    );
    if structure.m() != 0 {
        let _ = write!(
            out,
            "\n       TOTAL PATCHES USED: {}   PATCHES IN A SYMMETRIC CELL: {}",
            structure.m(),
            structure.mp
        );
    }
    let sections = structure.symmetry_sections();
    if sections != 1 {
        if structure.ipsym < 0 {
            let _ = writeln!(out, "\n  STRUCTURE HAS {sections} FOLD ROTATIONAL SYMMETRY");
        } else {
            let planes = if sections == 8 { 3 } else { sections / 2 };
            let _ = writeln!(out, "\n  STRUCTURE HAS {planes} PLANES OF SYMMETRY");
        }
    }
    if !structure.junctions.is_empty() {
        out.push_str(concat!(
            "\n\n",
            "    ---------- MULTIPLE WIRE JUNCTIONS ----------\n",
            "    JUNCTION  SEGMENTS (- FOR END 1, + FOR END 2)",
        ));
        for (number, junction) in structure.junctions.iter().enumerate() {
            let _ = write!(out, "\n   {:5}      ", number + 1);
            for (k, segment) in junction.iter().enumerate() {
                let _ = write!(out, "{segment:5}");
                if (k + 1) % 20 == 0 {
                    out.push_str("\n              ");
                }
            }
        }
    }

    if n != 0 {
        out.push_str(concat!(
            "\n\n\n",
            "                               ---------- SEGMENTATION DATA ----------\n",
            "                                        COORDINATES IN METERS\n",
            "                            I+ AND I- INDICATE THE SEGMENTS BEFORE AND AFTER I\n",
            "\n",
            "   SEG    COORDINATES OF SEGM CENTER     SEGM    ORIENTATION",
            " ANGLES    WIRE    CONNECTION DATA   TAG\n",
            "   No:       X         Y         Z      LENGTH     ALPHA     ",
            " BETA    RADIUS    I-     I    I+   No:",
        ));
    }
    let s = &structure;
    for i in 0..n {
        let alpha = s.salp[i].asin() * TD;
        let beta = (s.y2[i] - s.y1[i]).atan2(s.x2[i] - s.x1[i]) * TD;
        let _ = write!(
            out,
            "\n {:5} {} {} {} {} {} {} {} {:5} {:5} {:5} {:5}",
            i + 1,
            f(s.x[i], 9, 4),
            f(s.y[i], 9, 4),
            f(s.z[i], 9, 4),
            f(s.si[i], 9, 4),
            f(alpha, 9, 4),
            f(beta, 9, 4),
            f(s.bi[i], 9, 4),
            s.icon1[i],
            i + 1,
            s.icon2[i],
            s.itag[i]
        );
    }

    if s.m() != 0 {
        out.push_str(concat!(
            "\n\n\n",
            "                                    --------- SURFACE PATCH DATA ---------\n",
            "                                             COORDINATES IN METERS\n\n",
            " PATCH      COORD. OF PATCH CENTER           UNIT NORMAL VECTOR      ",
            " PATCH           COMPONENTS OF UNIT TANGENT VECTORS\n",
            "  No:       X          Y          Z          X        Y        Z      ",
            " AREA         X1       Y1       Z1        X2       Y2      Z2",
        ));
        for (i, p) in s.patches.iter().enumerate() {
            let [xn, yn, zn] = p.normal();
            let _ = write!(
                out,
                "\n {:4} {} {} {}  {} {} {} {}  {} {} {}  {} {} {}",
                i + 1,
                f(p.x, 10, 5),
                f(p.y, 10, 5),
                f(p.z, 10, 5),
                f(xn, 8, 4),
                f(yn, 8, 4),
                f(zn, 8, 4),
                f(p.area, 10, 5),
                f(p.t1x, 8, 4),
                f(p.t1y, 8, 4),
                f(p.t1z, 8, 4),
                f(p.t2x, 8, 4),
                f(p.t2y, 8, 4),
                f(p.t2z, 8, 4)
            );
        }
    }
    Ok(structure)
}

/// Network card as entered: type 1 two-port, 2 line, 3 crossed line.
#[derive(Debug, Clone, Copy)]
struct NetCard {
    ntyp: i32,
    iseg1: usize,
    iseg2: usize,
    x: [f64; 6],
}

/// State of one job: the engine and the parameters set by the program cards.
struct Job<'a> {
    out: &'a mut String,
    /// Plot data of all jobs, and the four flags of the last PL card: what
    /// to plot (1 currents, 2 near fields, 3 patterns) and in which form.
    plot: &'a mut Option<String>,
    plot_flags: [i32; 4],
    structure: &'a Structure,
    engine: Engine,
    mpcnt: i32,
    /// How much of the previous solution is still valid: 1 nothing, 2 the
    /// frequency scaling, 3 the factored matrix, 4 and 5 the currents.
    igo: i32,
    /// The kind of card read last, which decides whether a card adds to a
    /// group or starts a new one.
    iflow: i32,
    ifrq: i32,
    nfrq: i32,
    mhz: i32,
    fmhz: f64,
    delfrq: f64,
    rkh: f64,
    iexk: i32,
    ixtyp: i32,
    masym: i32,
    iped: i32,
    zpnorm: f64,
    nthi: i32,
    nphi: i32,
    xpr: [f64; 6],
    thetis: f64,
    phiss: f64,
    loads: Vec<Load>,
    networks: Vec<NetCard>,
    ntsol: bool,
    sources: Vec<(usize, Complex64)>,
    slope_sources: Vec<(usize, Complex64)>,
    ground: Option<GroundParams>,
    /// The structure has been checked against the current ground.
    ground_checked: bool,
    epsr2: f64,
    sig2: f64,
    clt: f64,
    cht: f64,
    /// Near-field request: -1 none, 0 rectangular grid, 1 spherical grid.
    near: i32,
    /// Magnetic (NH) instead of electric (NE) near field.
    nfeh: bool,
    nrx: i32,
    nry: i32,
    nrz: i32,
    /// First point and increments of the near-field grid.
    xnr: [f64; 3],
    dxnr: [f64; 3],
    iptflg: i32,
    iptag: i32,
    iptagf: i32,
    iptagt: i32,
    iptflq: i32,
    iptaq: i32,
    iptaqf: i32,
    iptaqt: i32,
    ifar: i32,
    nth: i32,
    nph: i32,
    ipd: i32,
    iavp: i32,
    inor: i32,
    iax: i32,
    thets: f64,
    phis: f64,
    dth: f64,
    dph: f64,
    rfld: f64,
    gnor: f64,
    nctag: Vec<i32>,
    ncseg: Vec<i32>,
    icoup: usize,
    y11a: Vec<Complex64>,
    y12a: Vec<Complex64>,
    fnorm: Vec<f64>,
    isave: usize,
    /// Position in the sweep of plane-wave incidence angles; `inc` counts
    /// the excitations and `nprint` silences the repeated network printout.
    nthic: i32,
    nphic: i32,
    inc: i32,
    nprint: bool,
    /// The loop over the incidence angles of a plane wave has run to its end.
    incidence_done: bool,
}

/// Which card started a computation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Trigger {
    ExecuteOrPattern,
    NearField,
}

impl<'a> Job<'a> {
    fn new(structure: &'a Structure, out: &'a mut String, plot: &'a mut Option<String>) -> Self {
        Job {
            out,
            plot,
            plot_flags: [0; 4],
            structure,
            engine: Engine::new(structure),
            mpcnt: 0,
            igo: 1,
            iflow: 1,
            ifrq: 0,
            nfrq: 1,
            mhz: 1,
            fmhz: nec2::consts::CVEL,
            delfrq: 0.0,
            rkh: 1.0,
            iexk: 0,
            ixtyp: 0,
            masym: 0,
            iped: 0,
            zpnorm: 0.0,
            nthi: 0,
            nphi: 0,
            xpr: [0.0; 6],
            thetis: 0.0,
            phiss: 0.0,
            loads: Vec::new(),
            networks: Vec::new(),
            ntsol: false,
            sources: Vec::new(),
            slope_sources: Vec::new(),
            ground: None,
            ground_checked: false,
            epsr2: 0.0,
            sig2: 0.0,
            clt: 0.0,
            cht: 0.0,
            near: -1,
            nfeh: false,
            nrx: 0,
            nry: 0,
            nrz: 0,
            xnr: [0.0; 3],
            dxnr: [0.0; 3],
            iptflg: -2,
            iptag: 0,
            iptagf: 0,
            iptagt: 0,
            iptflq: -1,
            iptaq: 0,
            iptaqf: 0,
            iptaqt: 0,
            ifar: -1,
            nth: 0,
            nph: 0,
            ipd: 0,
            iavp: 0,
            inor: 0,
            iax: 0,
            thets: 0.0,
            phis: 0.0,
            dth: 0.0,
            dph: 0.0,
            rfld: 0.0,
            gnor: 0.0,
            nctag: Vec::new(),
            ncseg: Vec::new(),
            icoup: 0,
            y11a: Vec::new(),
            y12a: Vec::new(),
            fnorm: Vec::new(),
            isave: 0,
            nthic: 1,
            nphic: 1,
            inc: 1,
            nprint: false,
            incidence_done: false,
        }
    }

    /// Reads program cards until NX (returns `Ok`) or EN.
    fn run(mut self, lines: &mut Lines) -> Flow {
        let result = self.run_cards(lines);
        // A job that ends normally leaves no card to blame for later errors.
        if matches!(result, Ok(()) | Err(Stop::End)) {
            lines.last = None;
        }
        result
    }

    fn run_cards(&mut self, lines: &mut Lines) -> Flow {
        loop {
            let Some(line) = lines.next_line() else {
                // The end of the file ends the job like an EN card.
                self.out.push_str("\n\n  TOTAL RUN TIME: 0 msec");
                return Err(Stop::End);
            };
            let line = line.text;
            let name = mnemonic(&line);
            if name == "CM" || name == "CE" {
                self.out.push_str("\n  WARNING: COMMENT CARD AFTER CE IS IGNORED");
                continue;
            }
            if GEOMETRY_CARDS.contains(&name) {
                return fail("\n  GEOMETRY CARD AFTER THE GE CARD");
            }
            if !CONTROL_CARDS.contains(&name) && name != "XT" {
                return fail("\n  UNKNOWN CARD");
            }
            let card = parse_card(&line, 4, 6, "COMMAND").map_err(Stop::Error)?;
            if card.mnemonic == "XT" {
                self.out.push_str("\n\n  Exiting after an \"XT\" command");
                return Err(Stop::End);
            }
            let [itmp1, itmp2, itmp3, itmp4] = card.ints[..] else { unreachable!("four fields are parsed") };
            let tmp: [f64; 6] = card.floats[..].try_into().expect("six fields are parsed");

            self.mpcnt += 1;
            let _ = write!(
                self.out,
                "\n  DATA CARD No: {:3} {} {itmp1:3} {itmp2:5} {itmp3:5} {itmp4:5} {} {} {} {} {} {}",
                self.mpcnt,
                card.mnemonic,
                e(tmp[0], 12, 5),
                e(tmp[1], 12, 5),
                e(tmp[2], 12, 5),
                e(tmp[3], 12, 5),
                e(tmp[4], 12, 5),
                e(tmp[5], 12, 5)
            );

            let trigger = match card.mnemonic.as_str() {
                // Frequency.
                "FR" => {
                    if itmp2 < 0 {
                        return fail("\n  THE NUMBER OF FREQUENCIES IS NEGATIVE");
                    }
                    if itmp1 != 0 && itmp1 != 1 {
                        let _ = write!(
                            self.out,
                            "\n  WARNING: STEPPING TYPE {itmp1} IS NOT DEFINED - AS IN NEC-2, THE STEP IS ADDED"
                        );
                    }
                    let steps = itmp2.max(1);
                    let last = if itmp1 == 1 { tmp[0] } else { tmp[0] + f64::from(steps - 1) * tmp[1] };
                    if !(tmp[0].is_finite() && tmp[0] > 0.0) {
                        return fail("\n  THE FREQUENCY MUST BE GREATER THAN ZERO");
                    }
                    if steps > 1 && itmp1 == 1 && !(tmp[1].is_finite() && tmp[1] > 0.0) {
                        return fail("\n  THE FREQUENCY MULTIPLIER MUST BE GREATER THAN ZERO");
                    }
                    if !(last.is_finite() && last > 0.0) {
                        return fail("\n  THE FREQUENCY SWEEP REACHES ZERO OR A NEGATIVE FREQUENCY");
                    }
                    self.ifrq = itmp1;
                    self.nfrq = steps;
                    self.fmhz = tmp[0];
                    self.delfrq = tmp[1];
                    if self.iped == 1 {
                        self.zpnorm = 0.0;
                    }
                    self.igo = 1;
                    self.iflow = 1;
                    continue;
                }

                // Loading.
                "LD" => {
                    if self.iflow != 3 {
                        self.iflow = 3;
                        self.loads.clear();
                        if self.igo > 2 {
                            self.igo = 2;
                        }
                        if itmp1 == -1 {
                            continue;
                        }
                    }
                    let ldtagt = if itmp4 == 0 { itmp3 } else { itmp4 };
                    if ldtagt < itmp3 {
                        return fail(format!(
                            "\n\n  DATA FAULT ON LOADING CARD No: {}: ITAG STEP1: {itmp3} IS GREATER THAN ITAG STEP2: {ldtagt}",
                            self.loads.len() + 1
                        ));
                    }
                    // A tag that no segment carries is an error of the
                    // loading itself.
                    if itmp3 != 0 && (itmp2 == 0 || self.tagged(itmp2) != 0) {
                        self.warn_empty_range(itmp2, itmp3, ldtagt);
                    }
                    self.loads.push(Load {
                        ldtyp: itmp1,
                        ldtag: itmp2,
                        ldtagf: itmp3,
                        ldtagt,
                        zlr: tmp[0],
                        zli: tmp[1],
                        zlc: tmp[2],
                    });
                    continue;
                }

                // Ground parameters under the antenna.
                "GN" => {
                    self.iflow = 4;
                    if self.igo > 2 {
                        self.igo = 2;
                    }
                    if itmp1 == -1 {
                        self.ground = None;
                        self.ground_checked = false;
                        self.engine.set_ground(None)?;
                        continue;
                    }
                    // NEC-2 defines the types -1 to 2 and does not check.
                    if itmp1 > 2 {
                        let _ = write!(
                            self.out,
                            concat!(
                                "\n  WARNING: GROUND TYPE {} IS NOT DEFINED. AS IN NEC-2, THE MATRIX IS FILLED AS FOR A ",
                                "PERFECT GROUND AND THE FIELDS USE THE REFLECTION COEFFICIENTS OF THE GROUND"
                            ),
                            itmp1
                        );
                    } else if itmp1 < -1 {
                        let _ = write!(
                            self.out,
                            "\n  WARNING: GROUND TYPE {itmp1} IS NOT DEFINED. AS IN NEC-2, IT IS TAKEN AS TYPE 0"
                        );
                    }
                    self.ground_checked = false;
                    if itmp1 != 1 && tmp[0] == 0.0 && tmp[1] == 0.0 {
                        return fail("\n  A FINITE GROUND NEEDS A DIELECTRIC CONSTANT OR A CONDUCTIVITY");
                    }
                    if itmp2 < 0 {
                        return fail("\n  THE NUMBER OF RADIAL WIRES IS NEGATIVE");
                    }
                    // Without radials the two fields describe the second medium.
                    let (scrwlt, scrwrt) = if itmp2 != 0 { (tmp[2], tmp[3]) } else { (0.0, 0.0) };
                    let params = GroundParams { iperf: itmp1, nradl: itmp2, epsr: tmp[0], sig: tmp[1], scrwlt, scrwrt };
                    if let Err(error) = self.engine.set_ground(Some(params)) {
                        return fail(format!("\n\n  {error}"));
                    }
                    self.ground = Some(params);
                    if itmp2 == 0 {
                        self.epsr2 = tmp[2];
                        self.sig2 = tmp[3];
                        self.clt = tmp[4];
                        self.cht = tmp[5];
                    }
                    continue;
                }

                // Excitation.
                "EX" => {
                    // NEC-2 has no other excitation: it reads beyond its
                    // table of names to print the type, crashes if the type
                    // is large and else takes the card for a current source.
                    if !(0..=5).contains(&itmp1) {
                        return fail(format!("\n  EXCITATION TYPE {itmp1} IS NOT DEFINED"));
                    }
                    // A request closes the group of EX cards. In NEC-2 an XQ card
                    // that asks for a pattern does not: the EX cards after it
                    // join the old group, and the next request gets the old
                    // currents, or goes on with a loop over incidence angles
                    // that has ended.
                    if self.iflow != 5 || self.igo > 3 {
                        self.sources.clear();
                        self.slope_sources.clear();
                        self.iped = 0;
                        self.iflow = 5;
                        if self.igo > 3 {
                            self.igo = 3;
                        }
                    }
                    self.ixtyp = itmp1;
                    self.incidence_done = false;
                    self.masym = itmp4 / 10;
                    if itmp1 == 0 || itmp1 == 5 {
                        self.ntsol = false;
                        let segment = self.isegno(itmp2, itmp3)?;
                        let mut voltage = Complex64::new(tmp[0], tmp[1]);
                        if voltage.norm() < 1.0e-20 {
                            voltage = Complex64::new(1.0, 0.0);
                        }
                        if itmp1 == 5 {
                            self.slope_sources.push((segment, voltage));
                        } else {
                            self.sources.push((segment, voltage));
                        }
                        self.iped = itmp4 - self.masym * 10;
                        self.zpnorm = tmp[2];
                        if self.iped == 1 && self.zpnorm > 0.0 {
                            self.iped = 2;
                        }
                        continue;
                    }
                    if itmp1 != 4 && (itmp2 < 0 || itmp3 < 0) {
                        return fail("\n  THE NUMBER OF INCIDENCE ANGLES IS NEGATIVE");
                    }
                    // Blank angle counts mean one angle.
                    self.nthi = itmp2.max(1);
                    self.nphi = itmp3.max(1);
                    self.xpr = tmp;
                    self.sources.clear();
                    self.slope_sources.clear();
                    self.thetis = tmp[0];
                    self.phiss = tmp[1];
                    continue;
                }

                // Networks and transmission lines.
                "NT" | "TL" => {
                    if self.iflow != 6 {
                        self.networks.clear();
                        self.ntsol = false;
                        self.iflow = 6;
                        if self.igo > 3 {
                            self.igo = 3;
                        }
                        if itmp2 == -1 {
                            continue;
                        }
                    }
                    let mut net = NetCard {
                        ntyp: if card.mnemonic == "NT" { 1 } else { 2 },
                        iseg1: self.isegno(itmp1, itmp2)?,
                        iseg2: self.isegno(itmp3, itmp4)?,
                        x: tmp,
                    };
                    if net.ntyp == 2 && tmp[0] == 0.0 {
                        return fail("\n  THE CHARACTERISTIC IMPEDANCE OF A TRANSMISSION LINE MUST NOT BE ZERO");
                    }
                    // A line with a negative impedance is a crossed line.
                    if net.ntyp == 2 && tmp[0] <= 0.0 {
                        net.ntyp = 3;
                        net.x[0] = -tmp[0];
                    }
                    self.networks.push(net);
                    continue;
                }

                // Execute.
                "XQ" => {
                    if (self.iflow == 10 && itmp1 == 0) || (self.nfrq == 1 && itmp1 == 0 && self.iflow > 7) {
                        continue;
                    }
                    if itmp1 == 0 {
                        self.iflow = if self.iflow > 7 { 11 } else { 7 };
                    } else {
                        // A pattern cut is requested together with the execution.
                        self.ifar = 0;
                        self.rfld = 0.0;
                        self.ipd = 0;
                        self.iavp = 0;
                        self.inor = 0;
                        self.iax = 0;
                        self.nth = 91;
                        self.nph = 1;
                        self.thets = 0.0;
                        self.phis = if itmp1 == 2 { 90.0 } else { 0.0 };
                        self.dth = 1.0;
                        self.dph = 0.0;
                        if itmp1 == 3 {
                            self.nph = 2;
                            self.dph = 90.0;
                        }
                    }
                    Trigger::ExecuteOrPattern
                }

                // Second ground medium.
                "GD" => {
                    self.epsr2 = tmp[0];
                    self.sig2 = tmp[1];
                    self.clt = tmp[2];
                    self.cht = tmp[3];
                    self.iflow = 9;
                    continue;
                }

                // Radiation pattern.
                "RP" => {
                    // NEC-2 computes a mode above 6 as mode 6 and prints a name
                    // for its cliff that it does not have.
                    if itmp1 > 6 {
                        return fail("\n  THE PATTERN MODE MUST BE 0 TO 6");
                    }
                    if itmp2 < 0 || itmp3 < 0 {
                        return fail("\n  THE NUMBER OF ANGLES IS NEGATIVE");
                    }
                    if itmp1 < 0 {
                        let _ = write!(
                            self.out,
                            "\n  WARNING: PATTERN MODE {itmp1} IS NOT DEFINED - AS IN NEC-2, {}",
                            if itmp1 == -1 { "NO PATTERN IS COMPUTED" } else { "IT IS TAKEN AS MODE 0" }
                        );
                    }
                    self.ifar = itmp1;
                    if (4..=6).contains(&itmp1)
                        && self.ground.is_none_or(|ground| ground.iperf == 1 || ground.nradl <= 0)
                    {
                        self.out.push_str(concat!(
                            "\n  WARNING: THIS PATTERN MODE IS FOR A RADIAL WIRE GROUND SCREEN, BUT THE GN CARD ",
                            "DEFINES NONE - THE PATTERN IS COMPUTED WITHOUT A SCREEN"
                        ));
                    }
                    self.nth = itmp2.max(1);
                    self.nph = itmp3.max(1);
                    // XNDA digits: axis selection, normalization, gain kind, averaging.
                    self.ipd = itmp4 / 10;
                    self.iavp = itmp4 - self.ipd * 10;
                    self.inor = self.ipd / 10;
                    self.ipd -= self.inor * 10;
                    self.iax = self.inor / 10;
                    self.inor -= self.iax * 10;
                    if self.iax != 0 {
                        self.iax = 1;
                    }
                    if self.ipd != 0 {
                        self.ipd = 1;
                    }
                    if self.nth < 2 || self.nph < 2 || self.ifar == 1 {
                        self.iavp = 0;
                    }
                    self.thets = tmp[0];
                    self.phis = tmp[1];
                    self.dth = tmp[2];
                    self.dph = tmp[3];
                    self.rfld = tmp[4];
                    self.gnor = tmp[5];
                    // NEC-2 prints values that are not numbers for a cliff
                    // without a medium beyond it, and for the ground wave at
                    // no distance from the origin.
                    if matches!(itmp1, 2 | 3 | 5 | 6) && self.ground.is_some() && self.epsr2 == 0.0 && self.sig2 == 0.0
                    {
                        return fail(concat!(
                            "\n  THE PATTERN MODE NEEDS THE MEDIUM BEYOND THE CLIFF: A GD CARD, OR ITS VALUES ON THE ",
                            "GN CARD"
                        ));
                    }
                    if itmp1 == 1 && self.rfld == 0.0 {
                        // The heights of the field points, accumulated as the pattern does.
                        let mut height = self.thets - self.dth;
                        for _ in 0..self.nth {
                            height += self.dth;
                            if height == 0.0 {
                                return fail(
                                    "\n  THE GROUND WAVE IS ASKED FOR AT THE ORIGIN: GIVE THE FIELD POINTS A DISTANCE",
                                );
                            }
                        }
                    }
                    self.iflow = 10;
                    Trigger::ExecuteOrPattern
                }

                // Next job.
                "NX" => return Ok(()),

                // Print control for currents.
                "PT" => {
                    if itmp3 != 0 {
                        self.warn_empty_range(itmp2, itmp3, if itmp4 == 0 { itmp3 } else { itmp4 });
                    }
                    self.iptflg = itmp1;
                    self.iptag = itmp2;
                    self.iptagf = itmp3;
                    self.iptagt = itmp4;
                    if itmp3 == 0 && self.iptflg != -1 {
                        self.iptflg = -2;
                    }
                    if itmp4 == 0 {
                        self.iptagt = self.iptagf;
                    }
                    continue;
                }

                // Interaction approximation range.
                "KH" => {
                    self.rkh = tmp[0];
                    if self.igo > 2 {
                        self.igo = 2;
                    }
                    self.iflow = 1;
                    continue;
                }

                // Near fields.
                "NE" | "NH" => {
                    self.nfeh = card.mnemonic == "NH";
                    if self.iflow == 8 && self.nfrq != 1 {
                        self.out.push_str(concat!(
                            "\n\n  WHEN MULTIPLE FREQUENCIES ARE REQUESTED, ",
                            "ONLY ONE NEAR FIELD CARD CAN BE USED -",
                            "\n  LAST CARD READ WILL BE USED",
                        ));
                    }
                    if itmp2 < 0 || itmp3 < 0 || itmp4 < 0 {
                        return fail("\n  THE NUMBER OF FIELD POINTS IS NEGATIVE");
                    }
                    self.near = itmp1;
                    // Blank point counts mean one point.
                    self.nrx = itmp2.max(1);
                    self.nry = itmp3.max(1);
                    self.nrz = itmp4.max(1);
                    self.xnr = [tmp[0], tmp[1], tmp[2]];
                    self.dxnr = [tmp[3], tmp[4], tmp[5]];
                    self.iflow = 8;
                    if self.nfrq != 1 {
                        continue;
                    }
                    Trigger::NearField
                }

                // Print control for charges.
                "PQ" => {
                    if itmp3 != 0 {
                        self.warn_empty_range(itmp2, itmp3, if itmp4 == 0 { itmp3 } else { itmp4 });
                    }
                    self.iptflq = itmp1;
                    self.iptaq = itmp2;
                    self.iptaqf = itmp3;
                    self.iptaqt = itmp4;
                    if itmp3 == 0 && self.iptflq != -1 {
                        self.iptflq = -2;
                    }
                    if itmp4 == 0 {
                        self.iptaqt = self.iptaqf;
                    }
                    continue;
                }

                // Extended thin-wire kernel.
                "EK" => {
                    self.iexk = if itmp1 == -1 { 0 } else { 1 };
                    if self.igo > 2 {
                        self.igo = 2;
                    }
                    self.iflow = 1;
                    continue;
                }

                // Maximum coupling between segments.
                "CP" => {
                    if self.iflow != 2 {
                        self.nctag.clear();
                        self.ncseg.clear();
                        self.y11a.clear();
                        self.y12a.clear();
                    }
                    self.icoup = 0;
                    self.iflow = 2;
                    // A segment that does not exist stops the run only when
                    // the coupling is computed, as in NEC-2.
                    for (tag, position) in [(itmp1, itmp2), (itmp3, itmp4)] {
                        if position == 0 {
                            break;
                        }
                        if self.structure.isegno(tag, position).is_err() {
                            let _ = write!(
                                self.out,
                                "\n  WARNING: THE STRUCTURE HAS NO SEGMENT {position} OF TAG {tag} FOR THE COUPLING"
                            );
                        }
                        self.nctag.push(tag);
                        self.ncseg.push(position);
                    }
                    continue;
                }

                // Plot flags. The data is collected from the first PL card on.
                "PL" => {
                    self.plot_flags = [itmp1, itmp2, itmp3, itmp4];
                    self.plot.get_or_insert_with(String::new);
                    continue;
                }

                "WG" => return Err(nec2::Error::Unsupported("NUMERICAL GREEN'S FUNCTION (WG CARD)").into()),

                "EN" => {
                    self.out.push_str("\n\n  TOTAL RUN TIME: 0 msec");
                    return Err(Stop::End);
                }

                _ => return fail("\n\n  FAULTY DATA CARD LABEL AFTER GEOMETRY SECTION"),
            };

            let from = self.out.len();
            self.compute(trigger)?;
            // Whatever the cause, a result that is not a number ends the
            // run. NEC-2 prints `NAN` and goes on.
            if self.out[from..].split_whitespace().any(not_a_number) {
                return fail("\n\n  THE RESULTS ABOVE HOLD VALUES THAT ARE NOT NUMBERS");
            }
        }
    }

    /// Number of segments that a card can select by position: those that
    /// carry `tag`, or all of them for tag 0.
    fn tagged(&self, tag: i32) -> usize {
        if tag == 0 { self.structure.n() } else { self.structure.itag.iter().filter(|&&t| t == tag).count() }
    }

    /// Warns when the positions `first` to `last` among the segments of
    /// `tag` select no segment. A range that reaches beyond the last segment
    /// selects the segments that exist, as in NEC-2.
    fn warn_empty_range(&mut self, tag: i32, first: i32, last: i32) {
        let count = self.tagged(tag);
        if first <= last && last >= 1 && (first.max(1) as usize) <= count {
            return;
        }
        let of_tag = if tag == 0 { String::new() } else { format!(" WITH TAG {tag}") };
        let _ = write!(
            self.out,
            "\n  WARNING: SEGMENTS {first} TO {last} ARE NOT AMONG THE {count} SEGMENTS{of_tag} - THE CARD SELECTS NO SEGMENT"
        );
    }

    /// Warns, once for every ground, about a structure that does not fit it:
    /// wire ends joined to an image without a ground, and segments below a
    /// ground or in it. NEC-2 solves such structures without a message, and
    /// so does this program, but the result belongs to no physical problem.
    fn check_ground(&mut self) {
        if self.ground_checked {
            return;
        }
        self.ground_checked = true;
        let s = self.structure;
        if self.ground.is_none() {
            let joined = (0..s.n())
                .map(|i| usize::from(s.icon1[i] == (i + 1) as i32) + usize::from(s.icon2[i] == (i + 1) as i32))
                .sum::<usize>();
            if joined != 0 {
                let _ = write!(
                    self.out,
                    concat!(
                        "\n  WARNING: SEGMENT ENDS JOINED TO THEIR IMAGE BY THE GE CARD: {}, BUT NO GN CARD DEFINES A ",
                        "GROUND - THE RESULT BELONGS TO NO PHYSICAL PROBLEM"
                    ),
                    joined
                );
            }
            return;
        }
        let misplaced: Vec<usize> = (0..s.n())
            .filter(|&i| {
                let tolerance = s.si[i] * 1.0e-3;
                s.z1[i] < -tolerance || s.z2[i] < -tolerance || (s.z1[i] <= tolerance && s.z2[i] <= tolerance)
            })
            .collect();
        if let Some(first) = misplaced.first() {
            let _ = write!(
                self.out,
                concat!(
                    "\n  WARNING: SEGMENTS THAT REACH BELOW THE GROUND PLANE OR LIE IN IT: {}, THE FIRST IS SEGMENT {} - ",
                    "NEC-2 CANNOT MODEL SUCH WIRES AND THE RESULT BELONGS TO NO PHYSICAL PROBLEM"
                ),
                misplaced.len(),
                first + 1
            );
        }
    }

    fn isegno(&self, tag: i32, m: i32) -> Flow<usize> {
        self.structure.isegno(tag, m).map_err(|error| Stop::Error(format!("\n\n  {error}")))
    }

    /// Runs the frequency loop from the point that the validity flag `igo`
    /// allows: the whole solution, or only the parts that later cards changed.
    fn compute(&mut self, trigger: Trigger) -> Flow {
        if self.iped != 0 {
            let need = 4 * self.nfrq.max(1) as usize;
            if self.fnorm.len() < need {
                self.fnorm.resize(need, 0.0);
            }
        }
        if self.iptflg >= 2 {
            let need = (self.nthi.max(0) * self.nphi.max(0)) as usize;
            if self.fnorm.len() < need {
                self.fnorm.resize(need, 0.0);
            }
        }

        // With valid currents, XQ and RP go straight to the far field.
        let mut entry = if trigger == Trigger::ExecuteOrPattern && self.igo == 5 { 6 } else { self.igo };
        if entry >= 4 && (1..=3).contains(&self.ixtyp) && self.incidence_done && self.nthi > 1 {
            // The loop over incidence angles is finished: a further request
            // starts it again from the first angle. NEC-2 would go on to an
            // angle beyond the last one. With a single theta it does not,
            // and the currents of the last angle serve the request.
            entry = 3;
            self.xpr[0] = self.thetis;
            self.xpr[1] = self.phiss;
        }
        if entry == 1 {
            self.mhz = 1;
        }

        let mut back_to_cards = false;
        loop {
            // --- frequency loop ---
            let mut next_frequency = false;

            if entry == 1 {
                if self.mhz != 1 {
                    if self.ifrq == 1 {
                        self.fmhz *= self.delfrq;
                    } else {
                        self.fmhz += self.delfrq;
                    }
                }
                self.engine.set_frequency(self.fmhz);
                let _ = write!(
                    self.out,
                    concat!(
                        "\n\n\n",
                        "                               --------- FREQUENCY --------\n",
                        "                                FREQUENCY :{} MHz\n",
                        "                                WAVELENGTH:{} Mtr",
                        "\n\n",
                        "                        APPROXIMATE INTEGRATION EMPLOYED FOR SEGMENTS \n",
                        "                        THAT ARE MORE THAN {:.3} WAVELENGTHS APART"
                    ),
                    e(self.fmhz, 11, 4),
                    e(self.engine.wavelength(), 11, 4),
                    self.rkh
                );
                if self.iexk == 1 {
                    self.out.push_str("\n                        THE EXTENDED THIN WIRE KERNEL WILL BE USED");
                }
            }

            if entry <= 2 {
                self.check_ground();
                self.load_and_factor()?;
            }

            if entry <= 3 {
                self.incidence_done = false;
                self.nthic = 1;
                self.nphic = 1;
                self.inc = 1;
                self.nprint = false;
            }
            let mut stage = entry;

            loop {
                // --- excitation loop ---
                if stage <= 3 {
                    self.excite_and_solve()?;
                    self.igo = 4;
                    self.couple()?;

                    if self.iflow == 7 {
                        if self.ixtyp > 0 && self.ixtyp < 4 {
                            stage = 3;
                            if self.next_incidence_angle() {
                                continue;
                            }
                            self.incidence_done = true;
                            break;
                        }
                        if self.nfrq != 1 {
                            next_frequency = true;
                            break;
                        }
                        self.out.push_str("\n\n\n");
                        back_to_cards = true;
                        break;
                    }
                }

                if stage <= 4 {
                    self.igo = 5;
                }

                // Near fields.
                if stage <= 5 && self.near != -1 {
                    self.near_field_pattern()?;
                    if self.mhz == self.nfrq {
                        self.near = -1;
                    }
                    if self.nfrq == 1 {
                        self.out.push_str("\n\n\n");
                        back_to_cards = true;
                        break;
                    }
                }

                // Far fields.
                if self.ifar != -1 {
                    self.radiation_pattern();
                }
                if self.ixtyp == 0 || self.ixtyp >= 4 {
                    if self.mhz == self.nfrq {
                        self.ifar = -1;
                    }
                    if self.nfrq != 1 {
                        next_frequency = true;
                        break;
                    }
                    self.out.push_str("\n\n\n");
                    back_to_cards = true;
                    break;
                }

                // Next incidence angle of the plane wave.
                stage = 3;
                if !self.next_incidence_angle() {
                    self.incidence_done = true;
                    break;
                }
            }

            if back_to_cards {
                break;
            }
            if !next_frequency {
                self.nphic = 1;
                self.xpr[1] = self.phiss;
                if self.iptflg >= 2 {
                    self.normalized_receiving_pattern();
                }
                if self.mhz == self.nfrq {
                    self.ifar = -1;
                }
                if self.nfrq == 1 {
                    self.out.push_str("\n\n\n");
                    back_to_cards = true;
                    break;
                }
            }

            self.mhz += 1;
            if self.mhz > self.nfrq {
                break;
            }
            entry = 1;
        }

        if back_to_cards {
            return Ok(());
        }
        if self.iped != 0 {
            self.input_impedance_table();
        }
        self.nfrq = 1;
        self.mhz = 1;
        Ok(())
    }

    /// Steps the plane-wave excitation to the next incidence angle: theta
    /// first, then phi. Returns `false` when all angles are done.
    fn next_incidence_angle(&mut self) -> bool {
        self.nthic += 1;
        self.inc += 1;
        self.xpr[0] += self.xpr[3];
        if self.nthic <= self.nthi {
            return true;
        }
        self.nthic = 1;
        self.xpr[0] = self.thetis;
        self.xpr[1] += self.xpr[4];
        self.nphic += 1;
        self.nphic <= self.nphi
    }

    /// Loading, antenna environment, matrix fill and factorization.
    fn load_and_factor(&mut self) -> Flow {
        self.out.push_str("\n\n\n                          ------ STRUCTURE IMPEDANCE LOADING ------");
        if self.loads.is_empty() {
            self.engine.clear_loads();
            self.out.push_str("\n                                 THIS STRUCTURE IS NOT LOADED");
        } else {
            self.out.push_str(concat!(
                "\n",
                "  LOCATION        RESISTANCE  INDUCTANCE  CAPACITANCE   ",
                "  IMPEDANCE (OHMS)   CONDUCTIVITY  CIRCUIT\n",
                "  ITAG FROM THRU     OHMS       HENRYS      FARADS     ",
                "  REAL     IMAGINARY   MHOS/METER      TYPE",
            ));
            let report = self.engine.set_loads(&self.loads)?;
            for line in &report.lines {
                self.out.push_str("\n ");
                let ints = [line.tag, line.from, line.thru];
                let all = ints == [0, 0, 0];
                if all {
                    self.out.push_str(" ALL");
                }
                for &value in &ints[usize::from(all)..] {
                    if value == 0 {
                        self.out.push_str("     ");
                    } else {
                        let _ = write!(self.out, "{value:5}");
                    }
                }
                for &value in &line.values {
                    if value.abs() >= 1.0e-20 {
                        let _ = write!(self.out, " {}", e(value, 11, 4));
                    } else {
                        self.out.push_str("            ");
                    }
                }
                let _ = write!(self.out, "   {}", line.kind);
            }
            if report.loaded_twice {
                self.out.push_str("\n  NOTE, SOME OF THE ABOVE SEGMENTS HAVE BEEN LOADED TWICE - IMPEDANCES ADDED");
            }
        }

        self.out.push_str("\n\n\n                            -------- ANTENNA ENVIRONMENT --------");
        let info = match self.engine.prepare_ground() {
            Ok(info) => info,
            Err(error) => return fail(format!("\n {error}")),
        };
        match info.kind {
            GroundType::FreeSpace => self.out.push_str("\n                            FREE SPACE"),
            GroundType::Perfect => self.out.push_str("\n                            PERFECT GROUND"),
            GroundType::ReflectionCoefficient | GroundType::Sommerfeld => {
                if let Some((nradl, scrwlt, scrwrt)) = info.radials {
                    let _ = write!(
                        self.out,
                        concat!(
                            "\n",
                            "                            RADIAL WIRE GROUND SCREEN\n",
                            "                            {} WIRES\n",
                            "                            WIRE LENGTH: {} METERS\n",
                            "                            WIRE RADIUS: {} METERS",
                            "\n",
                            "                            MEDIUM UNDER SCREEN -"
                        ),
                        nradl,
                        f(scrwlt, 8, 2),
                        e(scrwrt, 10, 3)
                    );
                }
                self.out.push_str(if info.kind == GroundType::Sommerfeld {
                    "\n                            FINITE GROUND - SOMMERFELD SOLUTION"
                } else {
                    "\n                            FINITE GROUND - REFLECTION COEFFICIENT APPROXIMATION"
                });
                let _ = write!(
                    self.out,
                    concat!(
                        "\n",
                        "                            RELATIVE DIELECTRIC CONST: {:.3}\n",
                        "                            CONDUCTIVITY: {} MHOS/METER\n",
                        "                            COMPLEX DIELECTRIC CONSTANT: {}{}j"
                    ),
                    info.epsr,
                    e(info.sig, 10, 3),
                    e(info.epsc.re, 11, 4),
                    e_plus(info.epsc.im, 11, 4)
                );
            }
        }

        self.engine.fill_and_factor(self.rkh, self.iexk == 1)?;
        self.print_warnings();
        // The conditioning indicator shares the line of the run times, which
        // nec2c prints too.
        let _ = write!(
            self.out,
            concat!(
                "\n\n\n",
                "                             ---------- MATRIX TIMING ----------\n",
                "                               FILL: 0 msec  FACTOR: 0 msec  SMALLEST/LARGEST PIVOT: {}",
            ),
            e(self.engine.pivot_ratio(), 10, 4)
        );
        self.ntsol = false;
        Ok(())
    }

    fn print_warnings(&mut self) {
        for warning in self.engine.take_warnings() {
            let _ = write!(self.out, "\n  {warning}");
        }
    }

    /// Sets up the excitation, solves for the currents and prints them.
    fn excite_and_solve(&mut self) -> Flow {
        let inc = self.inc;
        let wlam = self.engine.wavelength();
        let [xpr1, xpr2, xpr3, xpr4, xpr5, xpr6] = self.xpr;
        let ixtyp = self.ixtyp;
        if self.iptflg > 0 && !(1..=3).contains(&ixtyp) {
            // NEC-2 reads beyond its table of names for the heading of the
            // receiving pattern, and for the normalized pattern it writes to
            // a buffer that it has not got.
            return fail(concat!(
                "\n\n  THE PRINT OPTION OF THE PT CARD IS FOR A RECEIVING PATTERN: ",
                "IT NEEDS AN INCIDENT WAVE (EX 1 TO 3)"
            ));
        }
        if matches!(ixtyp, 0 | 5) && self.sources.is_empty() && self.slope_sources.is_empty() {
            // NEC-2 goes on: all currents are zero, and the efficiency and
            // the gains it prints are not numbers.
            return fail("\n\n  NO EX CARD DEFINES A SOURCE: THERE IS NOTHING TO SOLVE");
        }

        let excitation = if ixtyp == 0 || ixtyp == 5 {
            Excitation::Sources
        } else {
            if self.iptflg <= 0 || ixtyp == 4 {
                self.out.push_str("\n\n\n                             ---------- EXCITATION ----------");
            }
            if ixtyp == 4 {
                let _ = write!(
                    self.out,
                    concat!(
                        "\n",
                        "                                      CURRENT SOURCE\n",
                        "                     -- POSITION (METERS) -- ",
                        "      ORIENTATION (DEG)\n",
                        "                     X          Y          Z ",
                        "      ALPHA        BETA   DIPOLE MOMENT\n",
                        "               {} {} {}  {}     {}    {}"
                    ),
                    f(xpr1, 10, 5),
                    f(xpr2, 10, 5),
                    f(xpr3, 10, 5),
                    f(xpr4, 7, 2),
                    f(xpr5, 7, 2),
                    f(xpr6, 8, 3)
                );
                Excitation::CurrentElement {
                    position: [xpr1 / wlam, xpr2 / wlam, xpr3 / wlam],
                    alpha: TA * xpr4,
                    beta: TA * xpr5,
                    moment: xpr6 / (wlam * wlam),
                }
            } else {
                if self.iptflg <= 0 {
                    let _ = write!(
                        self.out,
                        "\n  PLANE WAVE - THETA: {} deg, PHI: {} deg, ETA={} DEG, TYPE - {}  AXIAL RATIO: {}",
                        f(xpr1, 7, 2),
                        f(xpr2, 7, 2),
                        f(xpr3, 7, 2),
                        HPOL[(ixtyp - 1).clamp(0, 2) as usize],
                        f(xpr6, 6, 3)
                    );
                }
                Excitation::PlaneWave {
                    polarization: ixtyp,
                    theta: TA * xpr1,
                    phi: TA * xpr2,
                    eta: TA * xpr3,
                    axial_ratio: xpr6,
                }
            }
        };

        self.engine.set_sources(&self.sources, &self.slope_sources);
        self.engine.excite(&excitation)?;

        if !self.networks.is_empty() && inc <= 1 {
            self.print_network_data();
        }
        if inc > 1 && self.iptflg > 0 {
            self.nprint = true;
        }

        let networks: Vec<Network> = self
            .networks
            .iter()
            .map(|net| Network {
                kind: match net.ntyp {
                    1 => NetworkKind::TwoPort,
                    2 => NetworkKind::Line,
                    _ => NetworkKind::CrossedLine,
                },
                iseg1: net.iseg1,
                iseg2: net.iseg2,
                x11r: net.x[0],
                x11i: net.x[1],
                x12r: net.x[2],
                x12i: net.x[3],
                x22r: net.x[4],
                x22i: net.x[5],
            })
            .collect();
        let report = self.engine.solve(&networks, self.masym != 0, self.ntsol)?;
        self.ntsol = true;

        if let Some(asym) = report.asymmetry {
            let _ = write!(
                self.out,
                concat!(
                    "\n\n",
                    "   MAXIMUM RELATIVE ASYMMETRY OF THE DRIVING POINT ADMITTANCE\n",
                    "   MATRIX IS {} FOR SEGMENTS {} AND {}\n",
                    "   RMS RELATIVE ASYMMETRY IS {}"
                ),
                e(asym.max, 10, 3),
                asym.segment1,
                asym.segment2,
                e(asym.rms, 10, 3)
            );
        }
        // The messages of the network solution come after the asymmetry, as in nec2c.
        self.print_warnings();
        if !networks.is_empty() && !self.nprint {
            self.out.push_str(concat!(
                "\n\n\n",
                "                          ",
                "--------- STRUCTURE EXCITATION DATA AT NETWORK CONNECTION POINTS --------",
                "\n",
                "  TAG   SEG       VOLTAGE (VOLTS)          CURRENT (AMPS)        ",
                " IMPEDANCE (OHMS)       ADMITTANCE (MHOS)     POWER\n",
                "  No:   No:     REAL      IMAGINARY     REAL      IMAGINARY    ",
                " REAL      IMAGINARY     REAL      IMAGINARY   (WATTS)",
            ));
            for record in &report.network_ports {
                self.print_port(record);
            }
        }
        if !report.sources.is_empty() {
            self.out.push_str(concat!(
                "\n\n\n",
                "                        --------- ANTENNA INPUT PARAMETERS ---------",
                "\n",
                "  TAG   SEG       VOLTAGE (VOLTS)         ",
                "CURRENT (AMPS)         IMPEDANCE (OHMS)    ",
                "    ADMITTANCE (MHOS)     POWER\n",
                "  No:   No:     REAL      IMAGINARY",
                "     REAL      IMAGINARY     REAL      ",
                "IMAGINARY    REAL       IMAGINARY   (WATTS)",
            ));
            for record in &report.sources {
                self.print_port(record);
            }
        }

        // Impedance of the last source, for the table printed after a sweep.
        if self.iped != 0 {
            let zped = self.engine.source_impedance();
            let at = 4 * (self.mhz - 1) as usize;
            self.fnorm[at] = zped.re;
            self.fnorm[at + 1] = zped.im;
            self.fnorm[at + 2] = zped.norm();
            self.fnorm[at + 3] = cang(zped);
            if self.iped != 2 && self.fnorm[at + 2] > self.zpnorm {
                self.zpnorm = self.fnorm[at + 2];
            }
        }

        self.print_currents(inc);

        if ixtyp == 0 || ixtyp == 5 {
            let pin = self.engine.input_power();
            let pnls = self.engine.network_loss();
            let ploss = self.engine.structure_loss();
            let radiated = pin - pnls - ploss;
            let _ = write!(
                self.out,
                concat!(
                    "\n\n\n",
                    "                               ---------- POWER BUDGET ---------\n",
                    "                               INPUT POWER   = {} Watts\n",
                    "                               RADIATED POWER= {} Watts\n",
                    "                               STRUCTURE LOSS= {} Watts\n",
                    "                               NETWORK LOSS  = {} Watts\n",
                    "                               EFFICIENCY    = {} Percent"
                ),
                e(pin, 11, 4),
                e(radiated, 11, 4),
                e(ploss, 11, 4),
                e(pnls, 11, 4),
                f(100.0 * radiated / pin, 7, 2)
            );
        }
        Ok(())
    }

    fn print_port(&mut self, record: &PortRecord) {
        let _ = write!(
            self.out,
            "\n {:4} {:5} {} {} {} {} {} {} {} {} {}",
            record.tag,
            record.segment,
            e(record.voltage.re, 11, 4),
            e(record.voltage.im, 11, 4),
            e(record.current.re, 11, 4),
            e(record.current.im, 11, 4),
            e(record.impedance.re, 11, 4),
            e(record.impedance.im, 11, 4),
            e(record.admittance.re, 11, 4),
            e(record.admittance.im, 11, 4),
            e(record.power, 11, 4)
        );
    }

    /// Prints the network tables, two-port networks and transmission lines
    /// separately, and gives lines without a length the straight distance
    /// between their ends.
    fn print_network_data(&mut self) {
        const PNET: [&str; 3] = ["        ", "STRAIGHT", " CROSSED"];
        self.out.push_str("\n\n\n                                            ---------- NETWORK DATA ----------");
        let wlam = self.engine.wavelength();
        let mut other_family = 0;
        let mut family = self.networks[0].ntyp;
        for _ in 0..2 {
            if family == 3 {
                family = 2;
            }
            if family == 2 {
                self.out.push_str(concat!(
                    "\n",
                    "  -- FROM -  --- TO --      TRANSMISSION LINE       ",
                    " --------- SHUNT ADMITTANCES (MHOS) ---------   LINE\n",
                    "  TAG   SEG  TAG   SEG    IMPEDANCE      LENGTH    ",
                    " ----- END ONE -----      ----- END TWO -----   TYPE\n",
                    "  No:   No:  No:   No:         OHMS      METERS     ",
                    " REAL      IMAGINARY      REAL      IMAGINARY",
                ));
            } else if family == 1 {
                self.out.push_str(concat!(
                    "\n",
                    "  -- FROM -  --- TO --            --------",
                    " ADMITTANCE MATRIX ELEMENTS (MHOS) ---------\n",
                    "  TAG   SEG  TAG   SEG   ----- (ONE,ONE) ------  ",
                    " ----- (ONE,TWO) -----   ----- (TWO,TWO) -------\n",
                    "  No:   No:  No:   No:      REAL      IMAGINARY     ",
                    " REAL     IMAGINARY       REAL      IMAGINARY",
                ));
            }

            for j in 0..self.networks.len() {
                let net = self.networks[j];
                if net.ntyp / family != 1 {
                    other_family = net.ntyp;
                    continue;
                }
                let (seg1, seg2) = (net.iseg1 - 1, net.iseg2 - 1);
                if net.ntyp >= 2 && net.x[1] <= 0.0 {
                    let a = self.engine.segment_center(seg1);
                    let b = self.engine.segment_center(seg2);
                    let (xx, yy, zz) = (b[0] - a[0], b[1] - a[1], b[2] - a[2]);
                    self.networks[j].x[1] = wlam * (xx * xx + yy * yy + zz * zz).sqrt();
                }
                let x = self.networks[j].x;
                let _ = write!(
                    self.out,
                    "\n {:4} {:5} {:4} {:5}  {} {}  {} {}  {} {} {}",
                    self.engine.tag(seg1),
                    net.iseg1,
                    self.engine.tag(seg2),
                    net.iseg2,
                    e(x[0], 11, 4),
                    e(x[1], 11, 4),
                    e(x[2], 11, 4),
                    e(x[3], 11, 4),
                    e(x[4], 11, 4),
                    e(x[5], 11, 4),
                    PNET[(net.ntyp - 1) as usize]
                );
            }

            if other_family == 0 {
                break;
            }
            family = other_family;
        }
    }

    fn current_line(&mut self, i: usize, value: Complex64) {
        let [x, y, z] = self.engine.segment_center(i);
        let _ = write!(
            self.out,
            "\n {:5} {:4} {} {} {} {} {} {} {} {}",
            i + 1,
            self.engine.tag(i),
            f(x, 9, 4),
            f(y, 9, 4),
            f(z, 9, 4),
            f(self.engine.segment_length(i), 9, 5),
            e(value.re, 11, 4),
            e(value.im, 11, 4),
            e(value.norm(), 11, 4),
            f(cang(value), 8, 3)
        );
    }

    /// Prints the currents on the segments and on the patches.
    fn print_currents(&mut self, inc: i32) {
        if self.engine.n() != 0 {
            self.print_segment_currents(inc);
        }
        if self.engine.m() == 0 {
            return;
        }
        self.out.push_str(concat!(
            "\n\n\n",
            "                                       --------- SURFACE PATCH CURRENTS ---------\n",
            "                                                 DISTANCE IN WAVELENGTHS\n",
            "                                                 CURRENT IN AMPS/METER\n\n",
            "                                 --------- SURFACE COMPONENTS --------   ",
            " ---------------- RECTANGULAR COMPONENTS ----------------\n",
            "  PCH   --- PATCH CENTER ---     TANGENT VECTOR 1    ",
            " TANGENT VECTOR 2    ------- X ------    ------- Y ------   ",
            " ------- Z ------\n  No:    X       Y       Z       MAG.    ",
            "   PHASE     MAG.       PHASE    REAL   IMAGINARY    REAL  ",
            " IMAGINARY    REAL   IMAGINARY",
        ));
        let currents = self.engine.patch_currents();
        for (i, p) in self.engine.patches().iter().enumerate() {
            let (ex, ey, ez) = (currents[3 * i], currents[3 * i + 1], currents[3 * i + 2]);
            let eth = ex * p.t1x + ey * p.t1y + ez * p.t1z;
            let eph = ex * p.t2x + ey * p.t2y + ez * p.t2z;
            let _ = write!(
                self.out,
                "\n {:4} {} {} {} {} {} {} {} {} {} {} {} {} {}",
                i + 1,
                f(p.x, 7, 3),
                f(p.y, 7, 3),
                f(p.z, 7, 3),
                e(eth.norm(), 11, 4),
                f(cang(eth), 8, 2),
                e(eph.norm(), 11, 4),
                f(cang(eph), 8, 2),
                e(ex.re, 9, 2),
                e(ex.im, 9, 2),
                e(ey.re, 9, 2),
                e(ey.im, 9, 2),
                e(ez.re, 9, 2),
                e(ez.im, 9, 2)
            );
            if self.plot_flags[0] == 1 {
                match self.plot_flags[2] {
                    1 => plot_line(self.plot, &[ex.re, ex.im]),
                    2 => plot_line(self.plot, &[ey.re, ey.im]),
                    3 => plot_line(self.plot, &[ez.re, ez.im]),
                    4 => plot_line(self.plot, &[ex.re, ex.im, ey.re, ey.im, ez.re, ez.im]),
                    _ => {}
                }
            }
        }
    }

    /// Prints the segment currents (or the receiving pattern) and charges
    /// as the PT and PQ cards ask.
    fn print_segment_currents(&mut self, inc: i32) {
        let n = self.engine.n();
        let wlam = self.engine.wavelength();
        let [xpr1, xpr2, xpr3, _, _, xpr6] = self.xpr;

        if self.iptflg != -1 {
            if self.iptflg <= 0 {
                self.out.push_str(concat!(
                    "\n\n\n",
                    "                           -------- CURRENTS AND LOCATION --------\n",
                    "                                  DISTANCES IN WAVELENGTHS",
                    "\n\n",
                    "   SEG  TAG    COORDINATES OF SEGM CENTER     SEGM",
                    "    ------------- CURRENT (AMPS) -------------\n",
                    "   No:  No:       X         Y         Z      LENGTH",
                    "     REAL      IMAGINARY    MAGN        PHASE",
                ));
            } else if self.iptflg != 3 && inc <= 1 {
                let _ = write!(
                    self.out,
                    concat!(
                        "\n\n\n",
                        "                        -------- RECEIVING PATTERN PARAMETERS --------\n",
                        "                                 ETA: {} DEGREES\n",
                        "                                 TYPE: {}\n",
                        "                                 AXIAL RATIO: {}\n\n",
                        "                        THETA     PHI      ----- CURRENT ----    SEG\n",
                        "                        (DEG)    (DEG)     MAGNITUDE    PHASE    No:"
                    ),
                    f(xpr3, 7, 2),
                    HPOL[(self.ixtyp - 1).clamp(0, 2) as usize],
                    f(xpr6, 6, 3)
                );
            }
        }

        let mut count = 0;
        for i in 0..n {
            let curi = self.engine.currents()[i] * wlam;
            if self.iptflg == -1 {
                continue;
            }
            if self.iptflg < 0 {
                self.current_line(i, curi);
                if self.plot_flags[0] == 1 {
                    match self.plot_flags[1] {
                        1 => plot_line(self.plot, &[curi.re, curi.im]),
                        2 => plot_line(self.plot, &[curi.norm(), cang(curi)]),
                        _ => {}
                    }
                }
                continue;
            }
            if self.iptag != 0 && self.engine.tag(i) != self.iptag {
                continue;
            }
            count += 1;
            if count < self.iptagf || count > self.iptagt {
                continue;
            }
            if self.iptflg == 0 {
                self.current_line(i, curi);
                continue;
            }
            let cmag = curi.norm();
            if self.iptflg >= 2 {
                self.fnorm[(inc - 1) as usize] = cmag;
                self.isave = i + 1;
            }
            if self.iptflg != 3 {
                let _ = write!(
                    self.out,
                    "\n                      {}  {}   {}  {}  {:5}",
                    f(xpr1, 7, 2),
                    f(xpr2, 7, 2),
                    e(cmag, 11, 4),
                    f(cang(curi), 7, 2),
                    i + 1
                );
            }
        }

        if self.iptflq != -1 {
            self.out.push_str(concat!(
                "\n\n\n",
                "                                  ------ CHARGE DENSITIES ------\n",
                "                                     DISTANCES IN WAVELENGTHS\n\n",
                "   SEG   TAG    COORDINATES OF SEG CENTER     SEG        ",
                "  CHARGE DENSITY (COULOMBS/METER)\n",
                "   No:   No:     X         Y         Z       LENGTH   ",
                "  REAL      IMAGINARY     MAGN       PHASE",
            ));
            let fr = 1.0e-6 / self.fmhz;
            let mut count = 0;
            for i in 0..n {
                if self.iptflq != -2 {
                    if self.iptaq != 0 && self.engine.tag(i) != self.iptaq {
                        continue;
                    }
                    count += 1;
                    if count < self.iptaqf || count > self.iptaqt {
                        continue;
                    }
                }
                let (bir, bii) = self.engine.sine_term(i);
                self.current_line(i, fr * Complex64::new(-bii, bir));
            }
        }
    }

    /// Accumulates the admittances for the CP card and prints the maximum
    /// coupling once every listed segment has been excited (`couple`).
    fn couple(&mut self) -> Flow {
        let ncoup = self.nctag.len();
        if ncoup == 0 || self.icoup >= ncoup || self.sources.len() != 1 || !self.slope_sources.is_empty() {
            return Ok(());
        }
        let j = self.isegno(self.nctag[self.icoup], self.ncseg[self.icoup])?;
        if j != self.sources[0].0 {
            return Ok(());
        }
        let wlam = self.engine.wavelength();
        let zin = self.sources[0].1;
        self.icoup += 1;
        let current = |engine: &Engine, segment: usize| engine.currents()[segment - 1] * wlam / zin;
        self.y11a.push(current(&self.engine, j));
        for i in 0..ncoup {
            if i + 1 == self.icoup {
                continue;
            }
            let k = self.isegno(self.nctag[i], self.ncseg[i])?;
            self.y12a.push(current(&self.engine, k));
        }
        if self.icoup < ncoup {
            return Ok(());
        }

        self.out.push_str(concat!(
            "\n\n\n",
            "                        -----------",
            " ISOLATION DATA -----------\n\n",
            " ------- COUPLING BETWEEN ------     MAXIMUM    ",
            " ---------- FOR MAXIMUM COUPLING ----------\n",
            "            SEG              SEG    COUPLING  LOAD",
            " IMPEDANCE (2ND SEG)         INPUT IMPEDANCE \n",
            " TAG  SEG   No:   TAG  SEG   No:      (DB)       ",
            " REAL     IMAGINARY         REAL       IMAGINARY",
        ));
        let npm1 = ncoup - 1;
        for i in 0..npm1 {
            let (itt1, its1) = (self.nctag[i], self.ncseg[i]);
            let isg1 = self.isegno(itt1, its1)?;
            for j in i + 1..ncoup {
                let (itt2, its2) = (self.nctag[j], self.ncseg[j]);
                let isg2 = self.isegno(itt2, its2)?;
                let j1 = j + i * npm1 - 1;
                let j2 = i + j * npm1;
                let y11 = self.y11a[i];
                let y22 = self.y11a[j];
                let y12 = 0.5 * (self.y12a[j1] + self.y12a[j2]);
                let mut yin = y12 * y12;
                let dbc = yin.norm();
                let c = dbc / (2.0 * y11.re * y22.re - yin.re);
                if (0.0..=1.0).contains(&c) {
                    let gmax = if c >= 0.01 { (1.0 - (1.0 - c * c).sqrt()) / c } else { 0.5 * (c + 0.25 * c * c * c) };
                    let rho = gmax * yin.conj() / dbc;
                    let yl = ((1.0 - rho) / (1.0 + rho) + 1.0) * y22.re - y22;
                    let zl = 1.0 / yl;
                    yin = y11 - yin / (y22 + yl);
                    let zin = 1.0 / yin;
                    let _ = write!(
                        self.out,
                        "\n {itt1:4} {its1:4} {isg1:5}  {itt2:4} {its2:4} {isg2:5}  {}  {} {}  {} {}",
                        f(db10(gmax), 9, 3),
                        e(zl.re, 12, 5),
                        e(zl.im, 12, 5),
                        e(zin.re, 12, 5),
                        e(zin.im, 12, 5)
                    );
                } else {
                    let _ = write!(
                        self.out,
                        "\n {itt1:4} {its1:4} {isg1:5}   {itt2:4} {its2:4} {isg2:5}  **ERROR** COUPLING IS NOT BETWEEN 0 AND 1. (= {})",
                        e(c, 12, 5)
                    );
                }
            }
        }
        Ok(())
    }

    /// Computes and prints the near fields requested by an NE or NH card
    /// (`nfpat`): on a rectangular grid, or on a grid in spherical
    /// coordinates (radius, phi, theta).
    fn near_field_pattern(&mut self) -> Flow {
        if !self.nfeh {
            self.out.push_str(concat!(
                "\n\n\n",
                "                             -------- NEAR ELECTRIC FIELDS --------\n",
                "     ------- LOCATION -------     ------- EX ------    ------- EY ------    ------- EZ ------\n",
                "      X         Y         Z       MAGNITUDE   PHASE    MAGNITUDE   PHASE    MAGNITUDE   PHASE\n",
                "    METERS    METERS    METERS     VOLTS/M  DEGREES    VOLTS/M   DEGREES     VOLTS/M  DEGREES",
            ));
        } else {
            self.out.push_str(concat!(
                "\n\n\n",
                "                                   -------- NEAR MAGNETIC FIELDS ---------\n\n",
                "     ------- LOCATION -------     ------- HX ------    ------- HY ------    ------- HZ ------\n",
                "      X         Y         Z       MAGNITUDE   PHASE    MAGNITUDE   PHASE    MAGNITUDE   PHASE\n",
                "    METERS    METERS    METERS      AMPS/M  DEGREES      AMPS/M  DEGREES      AMPS/M  DEGREES",
            ));
        }

        let wlam = self.engine.wavelength();
        let spherical = self.near != 0;
        let (mut cth, mut sth, mut cph, mut sph) = (0.0, 0.0, 0.0, 0.0);
        let mut znrt = self.xnr[2] - self.dxnr[2];
        for _ in 0..self.nrz {
            znrt += self.dxnr[2];
            if spherical {
                cth = (TA * znrt).cos();
                sth = (TA * znrt).sin();
            }
            let mut ynrt = self.xnr[1] - self.dxnr[1];
            for _ in 0..self.nry {
                ynrt += self.dxnr[1];
                if spherical {
                    cph = (TA * ynrt).cos();
                    sph = (TA * ynrt).sin();
                }
                let mut xnrt = self.xnr[0] - self.dxnr[0];
                for _ in 0..self.nrx {
                    xnrt += self.dxnr[0];
                    let (xob, yob, zob) =
                        if spherical { (xnrt * sth * cph, xnrt * sth * sph, xnrt * cth) } else { (xnrt, ynrt, znrt) };
                    let field = if self.nfeh {
                        self.engine.near_magnetic_field(xob / wlam, yob / wlam, zob / wlam)?
                    } else {
                        self.engine.near_electric_field(xob / wlam, yob / wlam, zob / wlam)?
                    };
                    let _ = write!(
                        self.out,
                        "\n {} {} {}  {} {}  {} {}  {} {}",
                        f(xob, 9, 4),
                        f(yob, 9, 4),
                        f(zob, 9, 4),
                        e(field[0].norm(), 11, 4),
                        f(cang(field[0]), 7, 2),
                        e(field[1].norm(), 11, 4),
                        f(cang(field[1]), 7, 2),
                        e(field[2].norm(), 11, 4),
                        f(cang(field[2]), 7, 2)
                    );
                    if self.plot_flags[0] == 2 {
                        let [_, form, component, axis] = self.plot_flags;
                        let along = if axis < 0 {
                            xob
                        } else if axis == 0 {
                            yob
                        } else {
                            zob
                        };
                        // Magnitude and phase, or real and imaginary part.
                        let pair = |v: Complex64| if form == 2 { [v.norm(), cang(v)] } else { [v.re, v.im] };
                        let [x, y, z] = [pair(field[0]), pair(field[1]), pair(field[2])];
                        if form == 1 || form == 2 {
                            match component {
                                1 => plot_line(self.plot, &[along, x[0], x[1]]),
                                2 => plot_line(self.plot, &[along, y[0], y[1]]),
                                3 => plot_line(self.plot, &[along, z[0], z[1]]),
                                4 => plot_line(self.plot, &[along, x[0], x[1], y[0], y[1], z[0], z[1]]),
                                _ => {}
                            }
                        }
                    }
                }
            }
        }
        self.print_warnings();
        Ok(())
    }

    /// Computes and prints the radiation pattern requested by the RP card.
    fn radiation_pattern(&mut self) {
        let request = PatternRequest {
            ifar: self.ifar,
            nth: self.nth,
            nph: self.nph,
            directive: self.ipd != 0,
            iavp: self.iavp,
            inor: self.inor,
            thets: self.thets,
            phis: self.phis,
            dth: self.dth,
            dph: self.dph,
            rfld: self.rfld,
            gnor: self.gnor,
            ixtyp: self.ixtyp,
            xpr6: self.xpr[5],
            epsr2: self.epsr2,
            sig2: self.sig2,
            clt: self.clt,
            cht: self.cht,
        };
        let pattern = self.engine.radiation_pattern(&request);
        self.print_pattern(&pattern);
    }

    fn print_pattern(&mut self, pattern: &RadiationPattern) {
        const IGTP: [&str; 2] = ["----- POWER GAINS ----- ", "--- DIRECTIVE GAINS ---"];
        const IGAX: [&str; 4] = [" MAJOR", " MINOR", " VERTC", " HORIZ"];
        const IGNTP: [&str; 5] = [" MAJOR AXIS", "  MINOR AXIS", "    VERTICAL", "  HORIZONTAL", "       TOTAL "];
        let ifar = self.ifar;
        let out = &mut *self.out;

        if ifar > 1 {
            out.push_str("\n\n\n                               ------ FAR FIELD GROUND PARAMETERS ------\n\n");
            if ifar > 3 {
                let (nradl, scrwlt, scrwrt) =
                    self.ground.map_or((0, 0.0, 0.0), |ground| (ground.nradl, ground.scrwlt, ground.scrwrt));
                let _ = write!(
                    out,
                    concat!(
                        "\n",
                        "                               --- RADIAL WIRE GROUND SCREEN ---\n",
                        "                               NUM OF WIRES= {}\n",
                        "                               WIRE LENGTH= {} METERS\n",
                        "                               WIRE RADIUS= {} METERS"
                    ),
                    nradl,
                    f(scrwlt, 8, 2),
                    e(scrwrt, 10, 3)
                );
            }
            if ifar != 4 {
                let hclif = if ifar == 2 || ifar == 5 { "LINEAR" } else { "CIRCULAR" };
                let _ = write!(
                    out,
                    concat!(
                        "\n",
                        "                               --- {} CLIFF ---\n",
                        "                               EDGE DISTANCE= {} METERS\n",
                        "                                      HEIGHT= {} METERS\n",
                        "                               --- SECOND MEDIUM ---\n",
                        "                               RELATIVE DIELECTRIC CONST= {}\n",
                        "                                     GROUND CONDUCTIVITY= {} MHOS"
                    ),
                    hclif,
                    f(self.clt, 9, 2),
                    f(self.cht, 9, 2),
                    f(self.epsr2, 10, 3),
                    f(self.sig2, 10, 3)
                );
            }
        }

        if ifar == 1 {
            out.push_str(concat!(
                "\n\n\n",
                "                             ------- RADIATED FIELDS NEAR GROUND --------\n\n",
                "    ------- LOCATION -------     --- E(THETA) ---    ",
                " ---- E(PHI) ----    --- E(RADIAL) ---\n",
                "      RHO    PHI        Z           MAG    PHASE     ",
                "    MAG    PHASE        MAG     PHASE\n",
                "    METERS DEGREES    METERS      VOLTS/M DEGREES   ",
                "   VOLTS/M DEGREES     VOLTS/M  DEGREES",
            ));
            for row in &pattern.ground_wave {
                let _ = write!(
                    out,
                    "\n {} {} {}  {} {}  {} {}  {} {}",
                    f(row.rho, 9, 2),
                    f(row.phi, 7, 2),
                    f(row.z, 9, 2),
                    e(row.e_theta_mag, 11, 4),
                    f(row.e_theta_phase, 7, 2),
                    e(row.e_phi_mag, 11, 4),
                    f(row.e_phi_phase, 7, 2),
                    e(row.e_radial_mag, 11, 4),
                    f(row.e_radial_phase, 7, 2)
                );
            }
        } else {
            out.push_str("\n\n\n                             ---------- RADIATION PATTERNS -----------\n");
            if let Some((exrm, exra)) = pattern.range_factor {
                let _ = write!(
                    out,
                    concat!(
                        "\n",
                        "                             RANGE: {} METERS\n",
                        "                             EXP(-JKR)/R: {} AT PHASE: {} DEGREES\n"
                    ),
                    e(self.rfld, 13, 6),
                    e(exrm, 12, 5),
                    f(exra, 7, 2)
                );
            }
            let axes = 2 * self.iax as usize;
            let _ = write!(
                out,
                concat!(
                    "\n",
                    " ---- ANGLES -----     {:>23}      ---- POLARIZATION ----  ",
                    " ---- E(THETA) ----    ----- E(PHI) ------\n",
                    "  THETA      PHI      {:>6}   {:>6}    TOTAL       AXIAL    ",
                    "  TILT  SENSE   MAGNITUDE    PHASE    MAGNITUDE     PHASE\n",
                    " DEGREES   DEGREES        DB       DB       DB       RATIO  ",
                    " DEGREES            VOLTS/M   DEGREES     VOLTS/M   DEGREES"
                ),
                IGTP[self.ipd as usize],
                IGAX[axes],
                IGAX[axes + 1]
            );
            for row in &pattern.rows {
                let (gain1, gain2) = if self.iax != 1 {
                    (row.gain_major, row.gain_minor)
                } else {
                    (row.gain_vertical, row.gain_horizontal)
                };
                let _ = write!(
                    out,
                    "\n {} {}  {} {} {} {} {} {:>6} {} {} {} {}",
                    f(row.theta, 7, 2),
                    f(row.phi, 9, 2),
                    f(gain1, 8, 2),
                    f(gain2, 8, 2),
                    f(row.gain_total, 8, 2),
                    f(row.axial_ratio, 11, 4),
                    f(row.tilt, 9, 2),
                    row.sense,
                    e(row.e_theta_mag, 11, 4),
                    f(row.e_theta_phase, 9, 2),
                    e(row.e_phi_mag, 11, 4),
                    f(row.e_phi_phase, 9, 2)
                );
                if self.plot_flags[0] == 3 {
                    let [_, cut, component, gain] = self.plot_flags;
                    // The angle that runs along the cut.
                    let angle = match cut {
                        1 => row.theta,
                        2 => row.phi,
                        _ => continue,
                    };
                    match component {
                        1 => plot_line(self.plot, &[angle, row.e_theta_mag, row.e_theta_phase]),
                        2 => plot_line(self.plot, &[angle, row.e_phi_mag, row.e_phi_phase]),
                        _ => {}
                    }
                    match gain {
                        1 => plot_line(self.plot, &[angle, gain1]),
                        2 => plot_line(self.plot, &[angle, gain2]),
                        3 => plot_line(self.plot, &[angle, row.gain_total]),
                        _ => {}
                    }
                }
            }
        }

        if let Some(average) = pattern.average_gain {
            let _ = write!(
                out,
                "\n\n\n  AVERAGE POWER GAIN: {} - SOLID ANGLE USED IN AVERAGING: ({})*PI STERADIANS",
                average.gain.map_or("UNDEFINED".to_string(), |gain| e(gain, 11, 4)),
                f_plus(average.solid_angle / PI, 7, 4)
            );
            for warning in self.engine.take_warnings() {
                let _ = write!(out, "\n  {warning}");
            }
        }

        let Some(normalized) = &pattern.normalized else { return };
        let gmax = normalized.factor;
        let samples = &normalized.samples;
        let _ = write!(
            out,
            concat!(
                "\n\n\n",
                "                              ---------- NORMALIZED GAIN ----------\n",
                "                                      {:>6} GAIN\n",
                "                                   NORMALIZATION FACTOR: {:.2} db\n\n",
                "    ---- ANGLES ----                ---- ANGLES ----",
                "                ---- ANGLES ----\n",
                "    THETA      PHI        GAIN      THETA      PHI  ",
                "      GAIN      THETA      PHI       GAIN\n",
                "   DEGREES   DEGREES        DB     DEGREES   DEGREES ",
                "       DB     DEGREES   DEGREES       DB"
            ),
            IGNTP[(normalized.kind - 1).clamp(0, 4) as usize],
            gmax
        );

        // Three columns over the directions that were computed; the last
        // line may hold one or two entries.
        let total = samples.len();
        let rows = total.div_ceil(3);
        let missing = rows * 3 - total;
        let mut second = rows;
        let mut third = 2 * rows;
        if missing == 2 {
            third -= 1;
        }
        let entry = |index: usize| {
            let sample = samples[index - 1];
            format!("{} {} {}", f(sample.theta, 9, 2), f(sample.phi, 9, 2), f(sample.gain - gmax, 9, 2))
        };
        for i in 0..rows {
            second += 1;
            third += 1;
            if i + 1 == rows && missing != 0 {
                if missing != 2 {
                    let _ = write!(out, "\n {}   {}   ", entry(i + 1), entry(second));
                } else {
                    let _ = write!(out, "\n {}   ", entry(i + 1));
                }
                return;
            }
            let _ = write!(out, "\n {}   {}   {}", entry(i + 1), entry(second), entry(third));
        }
    }

    /// Prints the receiving pattern normalized to its maximum (PT 2 and 3).
    fn normalized_receiving_pattern(&mut self) {
        let total = (self.nthi * self.nphi).max(0) as usize;
        let maximum = self.fnorm[..total].iter().copied().fold(f64::MIN, f64::max);
        let [mut xpr1, mut xpr2, xpr3, xpr4, xpr5, xpr6] = self.xpr;
        let _ = write!(
            self.out,
            concat!(
                "\n\n\n",
                "                     ---- NORMALIZED RECEIVING PATTERN ----\n",
                "                      NORMALIZATION FACTOR: {}\n",
                "                      ETA: {} DEGREES\n",
                "                      TYPE: {}\n",
                "                      AXIAL RATIO: {}\n",
                "                      SEGMENT No: {}\n\n",
                "                      THETA     PHI       ---- PATTERN ----\n",
                "                      (DEG)    (DEG)       DB     MAGNITUDE"
            ),
            e(maximum, 11, 4),
            f(xpr3, 7, 2),
            HPOL[(self.ixtyp - 1).clamp(0, 2) as usize],
            f(xpr6, 6, 3),
            self.isave
        );
        for j in 0..self.nphi.max(0) as usize {
            for i in 0..self.nthi.max(0) as usize {
                // Without any received current there is nothing to normalize.
                let ratio = if maximum > 0.0 { self.fnorm[i + j * self.nthi as usize] / maximum } else { 0.0 };
                let _ = write!(
                    self.out,
                    "\n                    {}  {}   {}  {}",
                    f(xpr1, 7, 2),
                    f(xpr2, 7, 2),
                    f(db20(ratio), 7, 2),
                    e(ratio, 11, 4)
                );
                xpr1 += xpr4;
            }
            xpr1 = self.thetis;
            xpr2 += xpr5;
        }
        self.xpr[0] = xpr1;
        self.xpr[1] = self.phiss;
    }

    /// Prints the source impedance at every frequency of a sweep (EX card
    /// with the impedance-table option).
    fn input_impedance_table(&mut self) {
        let iss = self.slope_sources.last().or(self.sources.last()).map_or(0, |source| source.0);
        let _ = write!(
            self.out,
            concat!(
                "\n\n\n",
                "                             -------- INPUT IMPEDANCE DATA --------\n",
                "                                      SOURCE SEGMENT No: {}\n",
                "                                   NORMALIZATION FACTOR:{}\n\n",
                "              ----------- UNNORMALIZED IMPEDANCE ----------  ",
                "  ------------ NORMALIZED IMPEDANCE -----------\n",
                "      FREQ    RESISTANCE    REACTANCE    MAGNITUDE    PHASE  ",
                "  RESISTANCE    REACTANCE    MAGNITUDE    PHASE\n",
                "       MHz       OHMS         OHMS         OHMS     DEGREES  ",
                "     OHMS         OHMS         OHMS     DEGREES"
            ),
            iss,
            e(self.zpnorm, 12, 5)
        );

        let steps = f64::from(self.nfrq - 1);
        let mut frequency =
            if self.ifrq == 1 { self.fmhz / self.delfrq.powf(steps) } else { self.fmhz - steps * self.delfrq };
        for i in 0..self.nfrq.max(0) as usize {
            let z = &self.fnorm[4 * i..4 * i + 4];
            let _ = write!(
                self.out,
                "\n {}   {}  {}  {}  {}   {}  {}  {}  {}",
                f(frequency, 9, 3),
                e(z[0], 11, 4),
                e(z[1], 11, 4),
                e(z[2], 11, 4),
                f(z[3], 7, 2),
                e(z[0] / self.zpnorm, 11, 4),
                e(z[1] / self.zpnorm, 11, 4),
                e(z[2] / self.zpnorm, 11, 4),
                f(z[3], 7, 2)
            );
            if self.ifrq == 1 {
                frequency *= self.delfrq;
            } else {
                frequency += self.delfrq;
            }
        }
        self.out.push_str("\n\n\n");
    }
}
