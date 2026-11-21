//! Errors reported by the engine.

use std::fmt;

/// An error that stops a NEC-2 run.
///
/// The messages follow the wording of the NEC-2 output file so that they are
/// familiar to people who know the original program.
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    /// A wire, arc or helix was requested with a negative segment count.
    InvalidSegmentCount(i32),
    /// A count on a geometry card is negative.
    NegativeField(&'static str, i32),
    /// The turn spacing of a helix is zero.
    ZeroHelixSpacing,
    /// Segment ends lie within the connection tolerance of each other
    /// without forming a junction; the segment (1-based) where it was found.
    JunctionNotClosed(usize),
    /// The sweep of an arc is larger than 360 degrees.
    ArcAngleExceeds360,
    /// A rotation (GR) was requested with fewer than one sector.
    InvalidRotationCount(i32),
    /// A segment selected by position within a tag group used position zero.
    ZeroSegmentPosition,
    /// No segment carries the requested tag (or not enough of them).
    NoSegmentWithTag(i32),
    /// A segment number is larger than the number of segments.
    SegmentOutOfRange(i32),
    /// A segment (1-based) reaches below the ground plane.
    SegmentBelowGround(usize),
    /// A segment (1-based) lies in the ground plane.
    SegmentInGroundPlane(usize),
    /// A segment (1-based) lies in or crosses a plane of symmetry.
    SegmentInSymmetryPlane(usize),
    /// The centre of a patch (1-based) lies in a plane of symmetry.
    PatchInSymmetryPlane(usize),
    /// The four corners of a quadrilateral patch do not lie in a plane.
    PatchNotPlanar,
    /// A patch (1-based) has a centre, a tangent vector or an area that is
    /// not finite.
    PatchData(usize),
    /// The patch (1-based) under the end of a segment (1-based) has no area,
    /// or is so small that its quarters lie within the connection tolerance
    /// of that end.
    PatchTooSmall { patch: usize, segment: usize },
    /// A segment (1-based) has zero length or a non-positive radius.
    SegmentData { segment: usize, tag: i32 },
    /// Segment centres that lie inside the volume of another wire.
    Overlap(Vec<crate::geometry::Overlap>),
    /// The symmetric cell does not divide the structure.
    Symmetry { cell: usize, total: usize },
    /// The connection data of a segment (1-based) is inconsistent.
    Connection { routine: &'static str, segment: usize },
    /// The structure has no segments.
    EmptyStructure,
    /// A loading type outside 0..=5 was requested.
    ImproperLoadType(i32),
    /// A loading card names a tag that no segment carries.
    LoadTagNotFound(i32),
    /// A radial-wire ground screen cannot be combined with the Sommerfeld ground.
    RadialsWithSommerfeld,
    /// The Sommerfeld interpolation tables do not match the ground parameters.
    GroundParameters,
    /// An adaptive integration met a value that is not finite, or did not
    /// end within its bound on the number of steps. Segments are 1-based.
    /// The observation segment is absent for a field point, and both are
    /// absent for the Sommerfeld tables.
    Integration { routine: &'static str, source: Option<usize>, observation: Option<usize> },
    /// A pivot of the factorization is zero or not finite, so the matrix
    /// ("INTERACTION" or "NETWORK") is singular. The equation is 1-based.
    SingularMatrix { matrix: &'static str, equation: usize, frequency: f64 },
    /// The solution for the currents holds values that are not finite.
    NonFiniteSolution { frequency: f64 },
    /// There is no memory for the interaction matrix of this many equations.
    TooLarge { equations: usize },
    /// A feature of NEC-2 that this engine does not implement yet.
    Unsupported(&'static str),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::InvalidSegmentCount(ns) => write!(f, "GEOMETRY DATA CARD ERROR: SEGMENT COUNT {ns} IS NEGATIVE"),
            Error::NegativeField(what, value) => write!(f, "GEOMETRY DATA CARD ERROR: {what} {value} IS NEGATIVE"),
            Error::ZeroHelixSpacing => {
                write!(f, "GEOMETRY DATA CARD ERROR: THE TURN SPACING OF A HELIX MUST NOT BE ZERO")
            }
            Error::JunctionNotClosed(i) => write!(
                f,
                "GEOMETRY DATA ERROR -- SEGMENT ENDS NEAR SEGMENT {i} ARE CLOSER THAN THE CONNECTION TOLERANCE WITHOUT FORMING A CLEAN JUNCTION"
            ),
            Error::ArcAngleExceeds360 => write!(f, "ERROR -- ARC ANGLE EXCEEDS 360 DEGREES"),
            Error::InvalidRotationCount(n) => {
                write!(f, "GEOMETRY DATA CARD ERROR: {n} ROTATIONS REQUESTED, AT LEAST 1 IS NEEDED")
            }
            Error::ZeroSegmentPosition => {
                write!(f, "CHECK DATA, PARAMETER SPECIFYING SEGMENT POSITION IN A GROUP OF EQUAL TAGS MUST NOT BE ZERO")
            }
            Error::NoSegmentWithTag(tag) => write!(f, "NO SEGMENT HAS AN ITAG OF {tag}"),
            Error::SegmentOutOfRange(i) => write!(f, "SEGMENT NUMBER {i} IS LARGER THAN THE NUMBER OF SEGMENTS"),
            Error::SegmentBelowGround(i) => write!(f, "GEOMETRY DATA ERROR -- SEGMENT {i} EXTENDS BELOW GROUND"),
            Error::SegmentInGroundPlane(i) => write!(f, "GEOMETRY DATA ERROR -- SEGMENT {i} LIES IN GROUND PLANE"),
            Error::SegmentInSymmetryPlane(i) => write!(f, "GEOMETRY DATA ERROR--SEGMENT {i} LIES IN PLANE OF SYMMETRY"),
            Error::PatchInSymmetryPlane(i) => write!(f, "GEOMETRY DATA ERROR--PATCH {i} LIES IN PLANE OF SYMMETRY"),
            Error::PatchNotPlanar => write!(f, "ERROR -- CORNERS OF QUADRILATERAL PATCH DO NOT LIE IN A PLANE"),
            Error::PatchData(i) => {
                write!(f, "PATCH DATA ERROR: THE CENTRE, A TANGENT VECTOR OR THE AREA OF PATCH {i} IS NOT FINITE")
            }
            Error::PatchTooSmall { patch, segment } => write!(
                f,
                "PATCH DATA ERROR: PATCH {patch} IS TOO SMALL FOR THE END OF SEGMENT {segment} THAT LIES ON IT"
            ),
            Error::SegmentData { segment, tag } => {
                write!(f, "SEGMENT DATA ERROR: SEGMENT {segment} (TAG {tag}) HAS A ZERO OR NON-FINITE LENGTH OR RADIUS")
            }
            Error::Overlap(pairs) => {
                write!(f, "GEOMETRY DATA ERROR -- WIRES OCCUPY THE SAME SPACE")?;
                for line in crate::geometry::overlap_lines(pairs) {
                    write!(f, "\n  {line}")?;
                }
                Ok(())
            }
            Error::Symmetry { cell, total } => write!(f, "SYMMETRY ERROR - NROW: {cell} NCOL: {total}"),
            Error::Connection { routine, segment } => {
                write!(f, "{routine} - SEGMENT CONNECTION ERROR FOR SEGMENT {segment}")
            }
            Error::EmptyStructure => write!(f, "THE STRUCTURE HAS NO SEGMENTS"),
            Error::ImproperLoadType(t) => write!(f, "IMPROPER LOAD TYPE CHOSEN, REQUESTED TYPE IS {t}"),
            Error::LoadTagNotFound(tag) => write!(f, "LOADING DATA CARD ERROR, NO SEGMENT HAS AN ITAG = {tag}"),
            Error::RadialsWithSommerfeld => {
                write!(f, "RADIAL WIRE G.S. APPROXIMATION MAY NOT BE USED WITH SOMMERFELD GROUND OPTION")
            }
            Error::GroundParameters => write!(f, "ERROR IN GROUND PARAMETERS"),
            Error::Integration { routine, source, observation } => {
                write!(f, "{routine} -- THE INTEGRATION MET A VALUE THAT IS NOT FINITE OR DID NOT END")?;
                if let Some(j) = source {
                    write!(f, ", SOURCE SEGMENT {j}")?;
                }
                if let Some(i) = observation {
                    write!(f, ", OBSERVATION SEGMENT {i}")?;
                }
                Ok(())
            }
            Error::SingularMatrix { matrix, equation, frequency } => write!(
                f,
                "THE {matrix} MATRIX IS SINGULAR AT {frequency} MHZ: THE PIVOT OF EQUATION {equation} IS ZERO OR NOT FINITE"
            ),
            Error::TooLarge { equations } => {
                write!(f, "THE STRUCTURE IS TOO LARGE: NO MEMORY FOR THE MATRIX OF {equations} EQUATIONS")
            }
            Error::NonFiniteSolution { frequency } => {
                write!(f, "THE SOLUTION FOR THE CURRENTS AT {frequency} MHZ IS NOT FINITE")
            }
            Error::Unsupported(what) => write!(f, "NOT SUPPORTED: {what}"),
        }
    }
}

impl std::error::Error for Error {}

/// Result alias used across the crate.
pub type Result<T> = std::result::Result<T, Error>;
